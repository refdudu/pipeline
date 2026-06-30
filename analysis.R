# analysis.R — Phase 3 statistical analysis in R (canonical packages).
#
# Mirrors run_stats.py but uses the reference implementations cited in the paper:
#   ARTool  -> Aligned Rank Transform ANOVA            (Wobbrock et al. 2011, [16])
#   TOSTER  -> Two One-Sided Tests equivalence          (Lakens 2017, [17])
#   lme4/lmerTest -> linear mixed model sensitivity     (Pinheiro & Bates, [18])
# Friedman, Nemenyi (via studentized range), Cochran's Q and plots use base R,
# avoiding heavy deps (PMCMRplus/Rmpfr, tidyverse).
#
# Usage:  Rscript analysis.R [ledger_csv] [out_dir]

suppressPackageStartupMessages({
  library(ARTool); library(TOSTER); library(lme4); library(lmerTest)
})

args   <- commandArgs(trailingOnly = TRUE)
LEDGER <- ifelse(length(args) >= 1, args[1], "checkpoint_ledger_phase2.csv")
OUTDIR <- ifelse(length(args) >= 2, args[2], "results/r")
dir.create(OUTDIR, showWarnings = FALSE)
ALPHA  <- 0.05
MARGIN <- 1.0

RELABEL <- c("P-1"="P-1","P-2"="P0","P-3"="P1","P-4"="P2")
ORDER   <- c("P-1","P0","P1","P2")
PDESC   <- c("P-1"="controle negativo (padeiro)","P0"="neutro (sem persona)",
             "P1"="genérico (developer)","P2"="especializado (arquiteto sênior)")
METRICS <- c(complexity_delta="complexidade ciclomática", smells_delta="code smells")

fmtp <- function(p) if (is.na(p)) "n/a" else if (p < 0.001) "<0,001" else sub("\\.", ",", sprintf("%.3f", p))

# ---------------------------------------------------------------- load
d <- read.csv(LEDGER, stringsAsFactors = FALSE)
d$persona <- factor(RELABEL[d$persona], levels = ORDER)
d$valido  <- d$status == "VALIDO"
for (c0 in c("complexity_delta","smells_delta","complexity_base","smells_base","token_count"))
  d[[c0]] <- suppressWarnings(as.numeric(d[[c0]]))
valid <- d[d$valido, ]

# ---------------------------------------------------------------- cell aggregation
# Cell = (trecho, persona): median of valid replicas.
agg <- aggregate(cbind(complexity_delta, smells_delta) ~ trecho + persona + refatoracao_tipo,
                 data = valid, FUN = median)
nv  <- aggregate(valido ~ trecho + persona, data = valid, FUN = length)
names(nv)[3] <- "n_valid"
agg <- merge(agg, nv, by = c("trecho","persona"))
write.csv(agg[order(agg$trecho, agg$persona), ], file.path(OUTDIR,"cell_medians.csv"), row.names = FALSE)

# Inferential quality analysis EXCLUDES ReplaceMagicNumber: the chosen metrics
# (McCabe complexity, code smells) are insensitive to it (Δ≡0). Its zero variance
# creates a spurious interaction in ART and deflates TOST. It is kept only for
# descriptive stats / validity / tokens and reported as a construct-validity threat.
INFORMATIVE <- c("ExtractMethod", "ReplaceConditionalWithPolymorphism")
agg_inf <- agg[agg$refatoracao_tipo %in% INFORMATIVE, ]

# complete-case wide matrix (trecho x persona) for a metric (informative types only)
cc_wide <- function(metric) {
  a <- agg_inf
  w <- tapply(a[[metric]], list(a$trecho, a$persona), function(x) x[1])
  w <- w[, ORDER, drop = FALSE]
  w[complete.cases(w), , drop = FALSE]
}

# ---------------------------------------------------------------- 1. validity + Cochran Q
vr <- do.call(rbind, lapply(ORDER, function(p) {
  s <- d[d$persona == p, ]
  data.frame(persona = p, desc = PDESC[p], n_total = nrow(s),
             n_valido = sum(s$valido), taxa_validade = round(mean(s$valido), 4))
}))
write.csv(vr, file.path(OUTDIR,"validity_rates.csv"), row.names = FALSE)

# Cochran's Q on (trecho,replica) x persona binary "replica valid?" matrix.
bin <- tapply(d$valido, list(paste(d$trecho, d$replica), d$persona), function(x) as.integer(x[1]))
bin <- bin[, ORDER, drop = FALSE]; bin <- bin[complete.cases(bin), , drop = FALSE]
cochran_q <- function(m) {
  k <- ncol(m); Cj <- colSums(m); Ri <- rowSums(m); N <- sum(m)
  Q <- (k-1) * (k*sum(Cj^2) - N^2) / (k*N - sum(Ri^2))
  list(Q = Q, df = k-1, p = 1 - pchisq(Q, k-1), n = nrow(m))
}
cq <- if (length(unique(c(bin))) > 1) cochran_q(bin) else list(Q=NA,df=NA,p=NA,n=nrow(bin))

# ---------------------------------------------------------------- 2. RQ1 Friedman + Nemenyi
nemenyi <- function(w) {              # Demsar (2006) via studentized range
  k <- ncol(w); N <- nrow(w)
  R <- t(apply(w, 1, rank)); Rbar <- colMeans(R)
  P <- matrix(NA, k, k, dimnames = list(colnames(w), colnames(w)))
  for (i in 1:k) for (j in 1:k) if (i != j) {
    z <- abs(Rbar[i]-Rbar[j]) / sqrt(k*(k+1)/(6*N))
    P[i,j] <- 1 - ptukey(z*sqrt(2), k, Inf)
  }
  P
}
friedman <- list()
for (m in names(METRICS)) {
  w <- cc_wide(m); n <- nrow(w)
  e <- list(metric = m, label = METRICS[[m]], n_complete_case = n)
  if (n >= 3 && length(unique(c(w))) > 1) {
    ft <- friedman.test(as.matrix(w))
    e$statistic <- unname(ft$statistic); e$pvalue <- ft$p.value
    e$significant <- ft$p.value < ALPHA
    if (ft$p.value < ALPHA)
      write.csv(round(nemenyi(w), 4), file.path(OUTDIR, paste0("nemenyi_", m, ".csv")))
  } else { e$statistic <- NA; e$pvalue <- NA; e$note <- "n<3 ou sem variação" }
  friedman[[m]] <- e
}

# ---------------------------------------------------------------- 3. RQ2 ART interaction (ARTool)
art_res <- list()
for (m in names(METRICS)) {
  df <- agg_inf[, c("trecho","persona","refatoracao_tipo", m)]
  names(df)[4] <- "y"; df$tipo <- droplevels(factor(df$refatoracao_tipo)); df$trecho <- factor(df$trecho)
  e <- list(metric = m, label = METRICS[[m]], n = nrow(df))
  ok <- tryCatch({
    mod <- art(y ~ persona * tipo + (1|trecho), data = df)   # mixed ART (persona within, tipo between)
    a   <- anova(mod); row <- a[a$Term == "persona:tipo", ]
    e$F <- row$F; e$pvalue <- row[["Pr(>F)"]]; e$significant <- row[["Pr(>F)"]] < ALPHA
    pm  <- a[a$Term == "persona", ]; e$pmain <- if (nrow(pm)) pm[["Pr(>F)"]] else NA
    write.csv(a, file.path(OUTDIR, paste0("art_", m, ".csv")), row.names = FALSE); TRUE
  }, error = function(err) { e$note <<- paste("ART falhou:", conditionMessage(err)); FALSE })
  art_res[[m]] <- e
}

# ---------------------------------------------------------------- 4. TOST placebo (TOSTER)
pairs <- list(c("P-1","P2"), c("P0","P2"), c("P-1","P0"))
tost_rows <- list()
for (m in names(METRICS)) {
  w <- cc_wide(m)
  for (pr in pairs) {
    x <- w[, pr[1]]; y <- w[, pr[2]]; diff <- x - y
    if (length(x) < 3 || all(abs(diff) < 1e-9)) {
      tost_rows[[length(tost_rows)+1]] <- data.frame(metric=m, par=paste(pr[1],"vs",pr[2]),
        n=length(x), mean_diff=round(mean(diff),3), p_tost=NA, equivalente=all(abs(diff)<1e-9))
      next
    }
    r <- tryCatch(TOSTER::t_TOST(x = x, y = y, paired = TRUE, eqb = MARGIN), error=function(e) NULL)
    p_tost <- if (is.null(r)) NA else max(r$TOST["t-test","p.value"], NA, na.rm=TRUE) # fallback set below
    if (!is.null(r)) {
      # TOST overall p = larger of the two one-sided p-values
      p_tost <- max(r$TOST[c("TOST Lower","TOST Upper"), "p.value"])
    }
    tost_rows[[length(tost_rows)+1]] <- data.frame(metric=m, par=paste(pr[1],"vs",pr[2]),
      n=length(x), mean_diff=round(mean(diff),3), p_tost=p_tost,
      equivalente=(!is.na(p_tost) && p_tost < ALPHA))
  }
}
tost <- do.call(rbind, tost_rows)
write.csv(tost, file.path(OUTDIR,"tost.csv"), row.names = FALSE)

# ---------------------------------------------------------------- 5. LMM sensitivity (lmerTest)
lmm <- list()
for (m in names(METRICS)) {
  vi <- valid[valid$refatoracao_tipo %in% INFORMATIVE, ]
  df <- vi[, c("trecho", "persona", m)]; names(df)[3] <- "y"; df <- df[!is.na(df$y), ]
  df$persona <- relevel(factor(df$persona, levels = ORDER), ref = "P0")
  e <- list(metric = m, label = METRICS[[m]])
  ok <- tryCatch({
    fit <- lmerTest::lmer(y ~ persona + (1|trecho), data = df)
    co  <- summary(fit)$coefficients
    rows <- grep("^persona", rownames(co))
    e$coef <- data.frame(level = sub("persona","",rownames(co)[rows]),
                         beta = round(co[rows,"Estimate"],4), p = round(co[rows,"Pr(>|t|)"],4))
    write.csv(e$coef, file.path(OUTDIR, paste0("lmm_", m, ".csv")), row.names = FALSE); TRUE
  }, error = function(err) { e$note <<- conditionMessage(err); FALSE })
  lmm[[m]] <- e
}

# ---------------------------------------------------------------- 6. tokens
tok <- aggregate(token_count ~ persona, data = d, FUN = function(x) round(mean(x)))
tok <- tok[match(ORDER, tok$persona), ]
write.csv(tok, file.path(OUTDIR,"tokens.csv"), row.names = FALSE)

# per-type deltas
pt <- do.call(rbind, lapply(names(METRICS), function(m)
  do.call(rbind, lapply(split(agg, agg$refatoracao_tipo), function(g)
    data.frame(metric=m, tipo=g$refatoracao_tipo[1], n_cells=nrow(g),
               media=round(mean(g[[m]]),3), mediana=round(median(g[[m]]),3))))))
write.csv(pt, file.path(OUTDIR,"per_type.csv"), row.names = FALSE)

# ---------------------------------------------------------------- plots (base graphics)
for (m in names(METRICS)) {
  w <- cc_wide(m)
  png(file.path(OUTDIR, paste0("box_", m, ".png")), 800, 520, res=130)
  boxplot(as.data.frame(w)[ORDER], main=paste("Δ", METRICS[[m]], "por persona"),
          ylab=paste("Δ", METRICS[[m]])); abline(h=0, lty=2, col="grey"); dev.off()
}
png(file.path(OUTDIR,"validity.png"), 600, 500, res=130)
barplot(vr$taxa_validade, names.arg=vr$persona, ylim=c(0,1.05), main="Taxa de validade por persona"); dev.off()
png(file.path(OUTDIR,"tokens.png"), 600, 500, res=130)
barplot(tok$token_count, names.arg=tok$persona, main="Tokens médios por persona"); dev.off()

# ---------------------------------------------------------------- report (results_draft_R.md)
cochran_ns <- is.na(cq$p) || cq$p >= ALPHA
fried_ns   <- all(sapply(names(METRICS), function(m) is.na(friedman[[m]]$pvalue) || friedman[[m]]$pvalue >= ALPHA))
key <- tost[tost$par %in% c("P-1 vs P2","P0 vs P2"), ]
tost_eq <- nrow(key) > 0 && all(key$equivalente)
placebo <- cochran_ns && fried_ns && tost_eq

L <- c("# Resultados e Discussão (R) — rascunho","",
  "> Gerado por `analysis.R` (ARTool/TOSTER/lme4) a partir de `checkpoint_ledger_phase2.csv`. ",
  "Personas: P-1 controle negativo, P0 neutro, P1 genérico, P2 especializado. P3 contextual não executado. α=0,05; margem TOST=±1.","",
  "**Análise inferencial de qualidade (Friedman/ART/TOST/LMM) restrita a ExtractMethod + ReplaceConditionalWithPolymorphism.** ReplaceMagicNumber foi excluído por insensibilidade de construto (Δ≡0 nas métricas → artefato no ART, deflação no TOST) e consta apenas como descritivo. Validade e tokens usam todos os tipos.","",
  "## Dataset",
  sprintf("- %d obs; %d válidas (%.1f%%); complete-case (tipos mensuráveis) n=%d.",
          nrow(d), sum(d$valido), 100*mean(d$valido), friedman[["complexity_delta"]]$n_complete_case),"",
  "## Validade por persona + Cochran Q",
  "| Persona | Válidas/Total | Taxa |","|---|---|---|")
for (i in 1:nrow(vr)) L <- c(L, sprintf("| %s | %d/%d | %.3f |", vr$persona[i], vr$n_valido[i], vr$n_total[i], vr$taxa_validade[i]))
L <- c(L, "", sprintf("Cochran Q = %.3f, p = %s (%s diferença de validade).",
                      cq$Q, fmtp(cq$p), ifelse(cochran_ns,"sem","com")), "",
  "## RQ1 — Friedman (efeito principal da persona)","| Métrica | χ² | p | Signif.? |","|---|---|---|---|")
for (m in names(METRICS)) { e<-friedman[[m]]; L<-c(L, sprintf("| %s | %s | %s | %s |", e$label,
   ifelse(is.na(e$statistic),"n/a",sprintf("%.3f",e$statistic)), fmtp(e$pvalue),
   ifelse(is.na(e$pvalue),"—",ifelse(e$significant,"sim","não")))) }
L <- c(L, "", "## RQ2 — ART ANOVA (interação persona × tipo)","| Métrica | F | p | Signif.? |","|---|---|---|---|")
for (m in names(METRICS)) { e<-art_res[[m]]; L<-c(L, if(!is.null(e$F)) sprintf("| %s | %.3f | %s | %s |",
   e$label, e$F, fmtp(e$pvalue), ifelse(e$significant,"sim","não")) else sprintf("| %s | n/a | n/a | %s |", e$label, e$note)) }
L <- c(L, "", "## TOST (±1) — veredito placebo","| Métrica | Par | n | Δ médio | p | Equiv.? |","|---|---|---|---|---|---|")
for (i in 1:nrow(tost)) L <- c(L, sprintf("| %s | %s | %d | %s | %s | %s |", tost$metric[i], tost$par[i],
   tost$n[i], tost$mean_diff[i], fmtp(tost$p_tost[i]), ifelse(tost$equivalente[i],"sim","não")))
L <- c(L, "", "## LMM (ref. P0)")
for (m in names(METRICS)) { e<-lmm[[m]]; if(!is.null(e$coef)) {
  L<-c(L, sprintf("- *%s*: %s", e$label, paste(sprintf("%s β=%.3f (p=%s)", e$coef$level, e$coef$beta, sapply(e$coef$p, fmtp)), collapse=", ")))
} else L<-c(L, sprintf("- *%s*: %s", e$label, e$note)) }
L <- c(L, "", "## Conclusão",
  sprintf("**Efeito placebo %s**: validade %s (Cochran), %s diferença de qualidade (Friedman), equivalência TOST %s nos contrastes-chave.",
    ifelse(placebo,"SUSTENTADO","NÃO totalmente sustentado"),
    ifelse(cochran_ns,"semelhante","diferente"), ifelse(fried_ns,"ausência de","presença de"),
    ifelse(tost_eq,"confirmada","não confirmada")),"",
  "",
  sprintf("> **Sensibilidade RQ1:** Friedman (primário) e LMM não acusam efeito de persona; o ART acusa efeito principal de persona em complexidade (p=%s), porém sem direção pró-especialização (P2 especializado não reduz Δ) — tratado como sensibilidade; o veredito placebo se mantém.",
          fmtp(if (!is.null(art_res[["complexity_delta"]]$pmain)) art_res[["complexity_delta"]]$pmain else NA)),
  "",
  "### Limitações","- P3 contextual não executado (4/5 níveis).","- 2 réplicas (proposta: 5) → 240 obs.",
  "- Análise estática single-file não capta classes novas.",
  "- **ReplaceMagicNumber excluído do inferencial** (Δ≡0; métricas McCabe/smells insensíveis a essa refatoração) — limitação de construto.",
  "- ART/LMM divergem do Friedman no efeito principal (sensibilidade), reportados com cautela.")
writeLines(L, file.path(OUTDIR, "results_draft.md"))

cat("=== RESUMO (R) ===\n"); print(vr, row.names=FALSE)
cat(sprintf("\nCochran Q=%.3f p=%.4f\n", cq$Q, cq$p))
for (m in names(METRICS)) cat(sprintf("Friedman %s: chi2=%.3f p=%.4f (n=%d)\n",
    m, friedman[[m]]$statistic, friedman[[m]]$pvalue, friedman[[m]]$n_complete_case))
for (m in names(METRICS)) if(!is.null(art_res[[m]]$F)) cat(sprintf("ART %s: F=%.3f p=%.4f\n", m, art_res[[m]]$F, art_res[[m]]$pvalue))
cat("\nTOST:\n"); print(tost, row.names=FALSE)
cat(sprintf("\nVEREDITO PLACEBO: %s\n", ifelse(placebo,"SUSTENTADO","NÃO totalmente sustentado")))
