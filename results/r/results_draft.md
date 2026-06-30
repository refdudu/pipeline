# Resultados e Discussão (R) — rascunho

> Gerado por `analysis.R` (ARTool/TOSTER/lme4) a partir de `checkpoint_ledger_phase2.csv`. 
Personas: P-1 controle negativo, P0 neutro, P1 genérico, P2 especializado. P3 contextual não executado. α=0,05; margem TOST=±1.

**Análise inferencial de qualidade (Friedman/ART/TOST/LMM) restrita a ExtractMethod + ReplaceConditionalWithPolymorphism.** ReplaceMagicNumber foi excluído por insensibilidade de construto (Δ≡0 nas métricas → artefato no ART, deflação no TOST) e consta apenas como descritivo. Validade e tokens usam todos os tipos.

## Dataset
- 240 obs; 236 válidas (98.3%); complete-case (tipos mensuráveis) n=20.

## Validade por persona + Cochran Q
| Persona | Válidas/Total | Taxa |
|---|---|---|
| P-1 | 59/60 | 0.983 |
| P0 | 60/60 | 1.000 |
| P1 | 59/60 | 0.983 |
| P2 | 58/60 | 0.967 |

Cochran Q = 2.400, p = 0,494 (sem diferença de validade).

## RQ1 — Friedman (efeito principal da persona)
| Métrica | χ² | p | Signif.? |
|---|---|---|---|
| complexidade ciclomática | 4.017 | 0,260 | não |
| code smells | 2.023 | 0,568 | não |

## RQ2 — ART ANOVA (interação persona × tipo)
| Métrica | F | p | Signif.? |
|---|---|---|---|
| complexidade ciclomática | 1.587 | 0,203 | não |
| code smells | 0.462 | 0,710 | não |

## TOST (±1) — veredito placebo
| Métrica | Par | n | Δ médio | p | Equiv.? |
|---|---|---|---|---|---|
| complexity_delta | P-1 vs P2 | 20 | -0.975 | 0,483 | não |
| complexity_delta | P0 vs P2 | 20 | -0.85 | 0,413 | não |
| complexity_delta | P-1 vs P0 | 20 | -0.125 | 0,041 | sim |
| smells_delta | P-1 vs P2 | 20 | 0.675 | 0,224 | não |
| smells_delta | P0 vs P2 | 20 | 0.625 | 0,212 | não |
| smells_delta | P-1 vs P0 | 20 | 0.05 | 0,005 | sim |

## LMM (ref. P0)
- *complexidade ciclomática*: P-1 β=-0.136 (p=0,811), P1 β=-0.649 (p=0,255), P2 β=0.843 (p=0,143)
- *code smells*: P-1 β=-0.027 (p=0,948), P1 β=-0.284 (p=0,494), P2 β=-0.684 (p=0,103)

## Conclusão
**Efeito placebo NÃO totalmente sustentado**: validade semelhante (Cochran), ausência de diferença de qualidade (Friedman), equivalência TOST não confirmada nos contrastes-chave.


> **Sensibilidade RQ1:** Friedman (primário) e LMM não acusam efeito de persona; o ART acusa efeito principal de persona em complexidade (p=0,003), porém sem direção pró-especialização (P2 especializado não reduz Δ) — tratado como sensibilidade; o veredito placebo se mantém.

### Limitações
- P3 contextual não executado (4/5 níveis).
- 2 réplicas (proposta: 5) → 240 obs.
- Análise estática single-file não capta classes novas.
- **ReplaceMagicNumber excluído do inferencial** (Δ≡0; métricas McCabe/smells insensíveis a essa refatoração) — limitação de construto.
- ART/LMM divergem do Friedman no efeito principal (sensibilidade), reportados com cautela.
