"""Phase 3: statistical analysis of the persona-prompting refactoring experiment.

Turns the raw per-round ledger (checkpoint_ledger_phase2.csv) into the answers to
RQ1/RQ2 and the placebo verdict, following Section V.E of the proposal:

  - cell aggregation (trecho x persona) by median of valid replicas;
  - validity rate per persona + Cochran's Q (attrition);
  - RQ1: complete-case Friedman + Nemenyi post-hoc (per quality metric);
  - RQ2: persona x refactoring-type interaction via Aligned Rank Transform ANOVA;
  - placebo: TOST equivalence (margin +/-1) on key persona pairs;
  - sensitivity: linear mixed model over all valid replicas;
  - cost: mean tokens per persona.

Outputs: results/*.csv, results/*.png, results_draft.md (Portuguese).

Usage:  .venv/bin/python run_stats.py [ledger_csv] [out_dir]
"""

import os
import sys
import json
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scipy.stats import friedmanchisquare, rankdata
import scikit_posthocs as sp
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm
from statsmodels.stats.contingency_tables import cochrans_q
from statsmodels.stats.weightstats import ttost_paired

# Directory persona label -> paper persona level (increasing specialization).
RELABEL = {"P-1": "P-1", "P-2": "P0", "P-3": "P1", "P-4": "P2"}
ORDER = ["P-1", "P0", "P1", "P2"]
PERSONA_DESC = {
    "P-1": "controle negativo (padeiro)",
    "P0": "neutro (sem persona)",
    "P1": "genérico (developer)",
    "P2": "especializado (arquiteto sênior)",
}
METRICS = [("complexity_delta", "complexidade ciclomática"), ("smells_delta", "code smells")]
ALPHA = 0.05
TOST_MARGIN = 1.0

# ----------------------------------------------------------------------------- helpers


def load_ledger(path):
    df = pd.read_csv(path)
    df["persona"] = df["persona"].map(RELABEL)
    if df["persona"].isna().any():
        raise ValueError("Unknown persona label in ledger.")
    df["persona"] = pd.Categorical(df["persona"], categories=ORDER, ordered=True)
    for col in ("complexity_delta", "smells_delta", "complexity_base", "smells_base",
                "token_count", "prompt_tokens", "candidates_tokens"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["valido"] = df["status"] == "VALIDO"
    return df


def aggregate_cells(df):
    """Cell = (trecho, persona). Median of valid-replica deltas + valid count."""
    valid = df[df["valido"]]
    agg = (valid.groupby(["trecho", "persona", "refatoracao_tipo"], observed=True)
           .agg(complexity_delta=("complexity_delta", "median"),
                smells_delta=("smells_delta", "median"),
                n_valid=("valido", "size"))
           .reset_index())
    return agg


def complete_case(cells, metric):
    """Wide matrix snippets x personas of the metric, keeping only snippets valid
    in ALL personas (preserves the within-subject pairing)."""
    wide = cells.pivot_table(index="trecho", columns="persona", values=metric, observed=True)
    wide = wide.reindex(columns=ORDER)
    return wide.dropna(axis=0, how="any")


# ----------------------------------------------------------------------------- analyses


def validity_analysis(df, out_dir):
    rows = []
    for p in ORDER:
        sub = df[df["persona"] == p]
        rows.append({"persona": p, "desc": PERSONA_DESC[p],
                     "n_total": len(sub), "n_valido": int(sub["valido"].sum()),
                     "taxa_validade": round(sub["valido"].mean(), 4)})
    tbl = pd.DataFrame(rows)
    tbl.to_csv(os.path.join(out_dir, "validity_rates.csv"), index=False)

    # Cochran's Q on the (trecho, replica) x persona binary "replica valid?" matrix
    # — the matched unit is one replica of one snippet, so attrition differences
    # between personas are visible (cell-level "any replica valid" is degenerate).
    binm = (df.pivot_table(index=["trecho", "replica"], columns="persona",
                           values="valido", observed=True)
            .reindex(columns=ORDER).dropna(axis=0, how="any").astype(int))
    if binm.shape[0] >= 1 and len(np.unique(binm.values)) > 1:
        q = cochrans_q(binm.values)
        cochran = {"statistic": float(q.statistic), "pvalue": float(q.pvalue),
                   "n_subjects": int(binm.shape[0]), "unidade": "(trecho, replica)"}
    else:
        cochran = {"statistic": float("nan"), "pvalue": float("nan"),
                   "n_subjects": int(binm.shape[0]),
                   "note": "sem variação -> Q indefinido"}
    with open(os.path.join(out_dir, "cochran_q.json"), "w") as f:
        json.dump(cochran, f, indent=2, ensure_ascii=False)
    return tbl, cochran


def friedman_nemenyi(cells, out_dir):
    out = {}
    for metric, label in METRICS:
        wide = complete_case(cells, metric)
        n = wide.shape[0]
        entry = {"metric": metric, "label": label, "n_complete_case": int(n)}
        # Zero-variance guard (e.g. all-zero deltas) -> Friedman undefined.
        if n < 3 or np.allclose(wide.values, wide.values[0, 0]):
            entry.update({"statistic": float("nan"), "pvalue": float("nan"),
                          "note": "n<3 ou sem variação -> Friedman não aplicável"})
            out[metric] = entry
            continue
        stat, p = friedmanchisquare(*[wide[c].values for c in ORDER])
        entry.update({"statistic": float(stat), "pvalue": float(p),
                      "significant": bool(p < ALPHA)})
        if p < ALPHA:
            nem = sp.posthoc_nemenyi_friedman(wide.values)
            nem.index = ORDER
            nem.columns = ORDER
            nem.to_csv(os.path.join(out_dir, f"nemenyi_{metric}.csv"))
            entry["nemenyi_file"] = f"nemenyi_{metric}.csv"
        out[metric] = entry
    with open(os.path.join(out_dir, "friedman.json"), "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    return out


def art_interaction(cells, out_dir):
    """Aligned Rank Transform ANOVA for the persona x type interaction (Wobbrock
    2011): align the response by removing the two main effects, rank, then run a
    factorial OLS ANOVA and read the interaction term. Exploratory."""
    out = {}
    for metric, label in METRICS:
        d = cells[["persona", "refatoracao_tipo", metric]].dropna().copy()
        d = d.rename(columns={metric: "y", "refatoracao_tipo": "tipo"})
        d["persona"] = d["persona"].astype(str)
        entry = {"metric": metric, "label": label, "n": int(len(d))}
        if d["y"].nunique() <= 1 or d["persona"].nunique() < 2 or d["tipo"].nunique() < 2:
            entry["note"] = "sem variação suficiente -> ART não aplicável"
            out[metric] = entry
            continue
        grand = d["y"].mean()
        a_eff = d.groupby("persona")["y"].transform("mean") - grand
        b_eff = d.groupby("tipo")["y"].transform("mean") - grand
        d["aligned"] = d["y"] - (grand + a_eff + b_eff)  # residual w.r.t. main effects
        d["rank"] = rankdata(d["aligned"].values)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model = smf.ols("rank ~ C(persona) * C(tipo)", data=d).fit()
            aov = anova_lm(model, typ=2)
        inter = "C(persona):C(tipo)"
        entry.update({"F": float(aov.loc[inter, "F"]), "pvalue": float(aov.loc[inter, "PR(>F)"]),
                      "significant": bool(aov.loc[inter, "PR(>F)"] < ALPHA)})
        aov.to_csv(os.path.join(out_dir, f"art_{metric}.csv"))
        out[metric] = entry
    with open(os.path.join(out_dir, "art.json"), "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    return out


def tost_placebo(cells, out_dir):
    """TOST equivalence (margin +/-1) on paired complete-case cell medians for the
    key persona contrasts that define the placebo reading."""
    pairs = [("P-1", "P2"), ("P0", "P2"), ("P-1", "P0")]
    rows = []
    for metric, label in METRICS:
        wide = complete_case(cells, metric)
        for a, b in pairs:
            x, y = wide[a].values, wide[b].values
            if len(x) < 3 or np.allclose(x - y, 0):
                rows.append({"metric": metric, "par": f"{a} vs {b}", "n": int(len(x)),
                             "mean_diff": float(np.mean(x - y)), "p_tost": float("nan"),
                             "equivalente": np.allclose(x - y, 0),
                             "note": "diferença nula ou n<3"})
                continue
            p, _, _ = ttost_paired(x, y, -TOST_MARGIN, TOST_MARGIN)
            rows.append({"metric": metric, "par": f"{a} vs {b}", "n": int(len(x)),
                         "mean_diff": round(float(np.mean(x - y)), 3),
                         "p_tost": float(p), "equivalente": bool(p < ALPHA)})
    tbl = pd.DataFrame(rows)
    tbl.to_csv(os.path.join(out_dir, "tost.csv"), index=False)
    return tbl


def lmm_sensitivity(df, out_dir):
    """Mixed model over ALL valid replicas: delta ~ persona, random intercept per
    snippet. Reference = P0 (neutral). Robustness check vs the complete-case tests."""
    out = {}
    valid = df[df["valido"]].copy()
    valid["persona"] = valid["persona"].astype(str)
    for metric, label in METRICS:
        d = valid[["trecho", "persona", metric]].dropna().rename(columns={metric: "y"})
        entry = {"metric": metric, "label": label, "n_obs": int(len(d))}
        if d["y"].nunique() <= 1:
            entry["note"] = "sem variação -> LMM não aplicável"
            out[metric] = entry
            continue
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                m = smf.mixedlm("y ~ C(persona, Treatment('P0'))", d, groups=d["trecho"]).fit()
            coefs = {}
            for name in m.params.index:
                if name.startswith("C(persona"):
                    coefs[name.split("[T.")[-1].rstrip("]")] = {
                        "coef": round(float(m.params[name]), 4),
                        "pvalue": round(float(m.pvalues[name]), 4)}
            entry["coeficientes_vs_P0"] = coefs
        except Exception as e:  # convergence / singular fits
            entry["note"] = f"falha no ajuste: {e}"
        out[metric] = entry
    with open(os.path.join(out_dir, "lmm.json"), "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    return out


def cost_analysis(df, out_dir):
    tbl = (df.groupby("persona", observed=True)["token_count"].mean().reindex(ORDER)
           .round(0).reset_index().rename(columns={"token_count": "tokens_medio"}))
    tbl.to_csv(os.path.join(out_dir, "tokens.csv"), index=False)
    return tbl


def per_type_table(cells, out_dir):
    rows = []
    for metric, _ in METRICS:
        for t, g in cells.groupby("refatoracao_tipo", observed=True):
            rows.append({"metric": metric, "tipo": t, "n_cells": int(len(g)),
                         "media": round(float(g[metric].mean()), 3),
                         "mediana": round(float(g[metric].median()), 3)})
    tbl = pd.DataFrame(rows)
    tbl.to_csv(os.path.join(out_dir, "per_type.csv"), index=False)
    cells.to_csv(os.path.join(out_dir, "cell_medians.csv"), index=False)
    return tbl


# ----------------------------------------------------------------------------- plots


def make_plots(cells, validity_tbl, tokens_tbl, out_dir):
    for metric, label in METRICS:
        wide = complete_case(cells, metric)
        if wide.shape[0]:
            fig, ax = plt.subplots(figsize=(6, 4))
            ax.boxplot([wide[c].values for c in ORDER], labels=ORDER)
            ax.axhline(0, color="grey", lw=0.8, ls="--")
            ax.set_title(f"Δ {label} por persona (complete-case, n={wide.shape[0]})")
            ax.set_ylabel(f"Δ {label} (negativo = melhora)")
            fig.tight_layout()
            fig.savefig(os.path.join(out_dir, f"box_{metric}.png"), dpi=130)
            plt.close(fig)

        # Interaction: delta by persona grouped by refactoring type.
        fig, ax = plt.subplots(figsize=(7, 4))
        types = sorted(cells["refatoracao_tipo"].unique())
        x = np.arange(len(ORDER))
        w = 0.8 / max(len(types), 1)
        for i, t in enumerate(types):
            means = [cells[(cells["persona"] == p) & (cells["refatoracao_tipo"] == t)][metric].mean()
                     for p in ORDER]
            ax.bar(x + i * w, means, w, label=t)
        ax.set_xticks(x + w * (len(types) - 1) / 2)
        ax.set_xticklabels(ORDER)
        ax.axhline(0, color="grey", lw=0.8, ls="--")
        ax.set_title(f"Δ {label}: persona × tipo")
        ax.set_ylabel(f"Δ {label} (média)")
        ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, f"interaction_{metric}.png"), dpi=130)
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(5, 4))
    ax.bar(validity_tbl["persona"], validity_tbl["taxa_validade"])
    ax.set_ylim(0, 1.05)
    ax.set_title("Taxa de validade por persona")
    ax.set_ylabel("proporção VALIDO")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "validity.png"), dpi=130)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5, 4))
    ax.bar(tokens_tbl["persona"], tokens_tbl["tokens_medio"])
    ax.set_title("Tokens médios por persona")
    ax.set_ylabel("tokens")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "tokens.png"), dpi=130)
    plt.close(fig)


# ----------------------------------------------------------------------------- report


def fmt_p(p):
    if p != p:  # NaN
        return "n/a"
    return "<0,001" if p < 0.001 else f"{p:.3f}".replace(".", ",")


def write_report(df, validity_tbl, cochran, friedman, art, tost, lmm, tokens_tbl,
                 per_type, cells, path):
    n_total = len(df)
    n_valid = int(df["valido"].sum())
    cc_n = friedman["complexity_delta"].get("n_complete_case", 0)

    fried_ns = all((friedman[m].get("pvalue") != friedman[m].get("pvalue")) or
                   (friedman[m].get("pvalue", 1) >= ALPHA) for m, _ in METRICS)
    cochran_ns = (cochran.get("pvalue") != cochran.get("pvalue")) or (cochran.get("pvalue", 1) >= ALPHA)
    tost_equiv = tost[tost["par"].isin(["P-1 vs P2", "P0 vs P2"])]
    tost_all_equiv = bool(tost_equiv["equivalente"].all()) if len(tost_equiv) else False
    placebo = fried_ns and cochran_ns and tost_all_equiv

    L = []
    L.append("# Resultados e Discussão (rascunho)\n")
    L.append("> Gerado por `run_stats.py` a partir de `checkpoint_ledger_phase2.csv`. "
             "Personas: **P-1** controle negativo (padeiro), **P0** neutro, **P1** genérico, "
             "**P2** especializado. O nível **P3 contextual** da proposta não foi executado "
             "(limitação). α = 0,05; margem de equivalência TOST = ±1.\n")

    L.append("## Visão geral do dataset\n")
    L.append(f"- Observações: **{n_total}** rodadas (30 trechos × 4 personas × 2 réplicas).")
    L.append(f"- Preservação de comportamento: **{n_valid}/{n_total}** válidas "
             f"({100*n_valid/n_total:.1f}%); {n_total-n_valid} falhas de compilação.")
    L.append(f"- Análise de qualidade em regime *complete-case*: **{cc_n}** trechos válidos em todas as personas.\n")

    L.append("## Taxa de validade por persona (RQ secundária)\n")
    L.append("| Persona | Descrição | Válidas/Total | Taxa |")
    L.append("|---|---|---|---|")
    for _, r in validity_tbl.iterrows():
        L.append(f"| {r['persona']} | {r['desc']} | {r['n_valido']}/{r['n_total']} | {r['taxa_validade']:.3f} |")
    L.append("")
    L.append(f"**Cochran Q** (atrito entre personas): Q = {cochran.get('statistic', float('nan')):.3f}, "
             f"p = {fmt_p(cochran.get('pvalue', float('nan')))} "
             f"({'sem' if cochran_ns else 'com'} diferença significativa de validade entre personas).\n")

    L.append("## RQ1 — Efeito principal da persona na qualidade\n")
    L.append("Teste de Friedman (pareado, complete-case) por métrica:\n")
    L.append("| Métrica | χ² | p | Significativo? |")
    L.append("|---|---|---|---|")
    for m, lab in METRICS:
        e = friedman[m]
        sig = "—" if e.get("pvalue") != e.get("pvalue") else ("sim" if e.get("significant") else "não")
        st = "n/a" if e.get("statistic") != e.get("statistic") else f"{e['statistic']:.3f}"
        L.append(f"| {lab} | {st} | {fmt_p(e.get('pvalue', float('nan')))} | {sig} |")
    L.append("")
    for m, lab in METRICS:
        if friedman[m].get("nemenyi_file"):
            L.append(f"- Post-hoc de Nemenyi para *{lab}*: ver `results/{friedman[m]['nemenyi_file']}`.")
    L.append("")

    L.append("## RQ2 — Interação persona × tipo de refatoração (exploratória)\n")
    L.append("ANOVA sobre postos alinhados (ART):\n")
    L.append("| Métrica | F (interação) | p | Significativo? |")
    L.append("|---|---|---|---|")
    for m, lab in METRICS:
        e = art[m]
        if "F" in e:
            sig = "sim" if e.get("significant") else "não"
            L.append(f"| {lab} | {e['F']:.3f} | {fmt_p(e['pvalue'])} | {sig} |")
        else:
            L.append(f"| {lab} | n/a | n/a | {e.get('note','—')} |")
    L.append("")

    L.append("## Deltas por tipo de refatoração\n")
    L.append("| Métrica | Tipo | n células | média Δ | mediana Δ |")
    L.append("|---|---|---|---|---|")
    for _, r in per_type.iterrows():
        L.append(f"| {r['metric']} | {r['tipo']} | {r['n_cells']} | {r['media']} | {r['mediana']} |")
    L.append("\n> **Ameaça à validade de construto:** ReplaceMagicNumber apresenta Δ≈0 em "
             "complexidade e smells — as métricas escolhidas (McCabe, contagem de smells) "
             "não capturam a melhoria semântica dessa refatoração.\n")

    L.append("## Veredito placebo — TOST (±1)\n")
    L.append("| Métrica | Par | n | Δ médio | p (TOST) | Equivalente? |")
    L.append("|---|---|---|---|---|---|")
    for _, r in tost.iterrows():
        L.append(f"| {r['metric']} | {r['par']} | {r['n']} | {r.get('mean_diff','')} | "
                 f"{fmt_p(r['p_tost'])} | {'sim' if r['equivalente'] else 'não'} |")
    L.append("")

    L.append("## Sensibilidade — Modelo linear de efeitos mistos\n")
    L.append("Efeito da persona (ref. = P0) sobre cada Δ, usando todas as réplicas válidas "
             "com intercepto aleatório por trecho:\n")
    for m, lab in METRICS:
        e = lmm[m]
        if "coeficientes_vs_P0" in e:
            parts = ", ".join(f"{k}: β={v['coef']} (p={fmt_p(v['pvalue'])})"
                              for k, v in e["coeficientes_vs_P0"].items())
            L.append(f"- *{lab}*: {parts}")
        else:
            L.append(f"- *{lab}*: {e.get('note','—')}")
    L.append("")

    L.append("## Custo (tokens)\n")
    L.append("| Persona | Tokens médios |")
    L.append("|---|---|")
    for _, r in tokens_tbl.iterrows():
        L.append(f"| {r['persona']} | {int(r['tokens_medio'])} |")
    L.append("")

    L.append("## Conclusão\n")
    verdict = ("**Efeito placebo SUSTENTADO**" if placebo else
               "**Efeito placebo NÃO totalmente sustentado**")
    L.append(f"{verdict}: leitura conjunta de (i) validade {'semelhante' if cochran_ns else 'diferente'} "
             f"entre personas (Cochran Q), (ii) {'ausência' if fried_ns else 'presença'} de diferença "
             f"significativa de qualidade (Friedman) e (iii) equivalência TOST "
             f"{'confirmada' if tost_all_equiv else 'não confirmada'} nos contrastes-chave "
             f"(controle/neutro vs especializada).")
    if placebo:
        L.append("\nOs dados não suportam H1: especializar a persona não produziu código "
                 "estruturalmente melhor, com custo de tokens semelhante entre níveis — "
                 "consistente com a hipótese de placebo da engenharia de prompts em refatoração.")
    L.append("\n### Limitações\n")
    L.append("- Nível **P3 contextual** não executado (4 de 5 níveis propostos).")
    L.append("- **2 réplicas** por condição (proposta previa 5) → 240 obs (não 750).")
    L.append("- Análise estática em **escopo single-file** não capta classes novas (ex.: "
             "ReplaceConditionalWithPolymorphism infla a complexidade do arquivo alterado).")
    L.append("- **ReplaceMagicNumber** com Δ≈0 nas métricas → baixa sensibilidade de construto.")
    L.append("- Um único modelo LLM (Gemini Flash) e linguagem (Java).")

    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")
    return placebo


# ----------------------------------------------------------------------------- main


def main(ledger_path="checkpoint_ledger_phase2.csv", out_dir="results"):
    os.makedirs(out_dir, exist_ok=True)
    df = load_ledger(ledger_path)
    cells = aggregate_cells(df)

    validity_tbl, cochran = validity_analysis(df, out_dir)
    per_type = per_type_table(cells, out_dir)
    friedman = friedman_nemenyi(cells, out_dir)
    art = art_interaction(cells, out_dir)
    tost = tost_placebo(cells, out_dir)
    lmm = lmm_sensitivity(df, out_dir)
    tokens_tbl = cost_analysis(df, out_dir)
    make_plots(cells, validity_tbl, tokens_tbl, out_dir)

    placebo = write_report(df, validity_tbl, cochran, friedman, art, tost, lmm,
                           tokens_tbl, per_type, cells, "results_draft.md")

    print("=== RESUMO ===")
    print(validity_tbl.to_string(index=False))
    print("\nCochran Q:", cochran)
    print("\nFriedman:", json.dumps(friedman, ensure_ascii=False))
    print("\nART:", json.dumps(art, ensure_ascii=False))
    print("\nTOST:\n", tost.to_string(index=False))
    print(f"\nVEREDITO PLACEBO: {'SUSTENTADO' if placebo else 'NÃO totalmente sustentado'}")
    print(f"\nSaídas em {out_dir}/ e results_draft.md")


if __name__ == "__main__":
    ledger = sys.argv[1] if len(sys.argv) > 1 else "checkpoint_ledger_phase2.csv"
    outd = sys.argv[2] if len(sys.argv) > 2 else "results"
    main(ledger, outd)
