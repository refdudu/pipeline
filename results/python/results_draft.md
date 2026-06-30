# Resultados e Discussão (rascunho)

> Gerado por `run_stats.py` a partir de `checkpoint_ledger_phase2.csv`. Personas: **P-1** controle negativo (padeiro), **P0** neutro, **P1** genérico, **P2** especializado. O nível **P3 contextual** da proposta não foi executado (limitação). α = 0,05; margem de equivalência TOST = ±1.

## Visão geral do dataset

- Observações: **240** rodadas (30 trechos × 4 personas × 2 réplicas).
- Preservação de comportamento: **236/240** válidas (98.3%); 4 falhas de compilação.
- Análise de qualidade em regime *complete-case*: **30** trechos válidos em todas as personas.

## Taxa de validade por persona (RQ secundária)

| Persona | Descrição | Válidas/Total | Taxa |
|---|---|---|---|
| P-1 | controle negativo (padeiro) | 59/60 | 0.983 |
| P0 | neutro (sem persona) | 60/60 | 1.000 |
| P1 | genérico (developer) | 59/60 | 0.983 |
| P2 | especializado (arquiteto sênior) | 58/60 | 0.967 |

**Cochran Q** (atrito entre personas): Q = 2.400, p = 0,494 (sem diferença significativa de validade entre personas).

## RQ1 — Efeito principal da persona na qualidade

Teste de Friedman (pareado, complete-case) por métrica:

| Métrica | χ² | p | Significativo? |
|---|---|---|---|
| complexidade ciclomática | 4.017 | 0,260 | não |
| code smells | 2.023 | 0,568 | não |


## RQ2 — Interação persona × tipo de refatoração (exploratória)

ANOVA sobre postos alinhados (ART):

| Métrica | F (interação) | p | Significativo? |
|---|---|---|---|
| complexidade ciclomática | 1.977 | 0,075 | não |
| code smells | 1.167 | 0,329 | não |

## Deltas por tipo de refatoração

| Métrica | Tipo | n células | média Δ | mediana Δ |
|---|---|---|---|---|
| complexity_delta | ExtractMethod | 40 | 5.213 | 5.5 |
| complexity_delta | ReplaceConditionalWithPolymorphism | 40 | 8.7 | 7.25 |
| complexity_delta | ReplaceMagicNumber | 40 | 0.0 | 0.0 |
| smells_delta | ExtractMethod | 40 | -0.087 | -1.0 |
| smells_delta | ReplaceConditionalWithPolymorphism | 40 | 0.087 | 0.0 |
| smells_delta | ReplaceMagicNumber | 40 | 0.0 | 0.0 |

> **Ameaça à validade de construto:** ReplaceMagicNumber apresenta Δ≈0 em complexidade e smells — as métricas escolhidas (McCabe, contagem de smells) não capturam a melhoria semântica dessa refatoração.

## Veredito placebo — TOST (±1)

| Métrica | Par | n | Δ médio | p (TOST) | Equivalente? |
|---|---|---|---|---|---|
| complexity_delta | P-1 vs P2 | 30 | -0.65 | 0,187 | não |
| complexity_delta | P0 vs P2 | 30 | -0.567 | 0,172 | não |
| complexity_delta | P-1 vs P0 | 30 | -0.083 | 0,003 | sim |
| smells_delta | P-1 vs P2 | 30 | 0.45 | 0,031 | sim |
| smells_delta | P0 vs P2 | 30 | 0.417 | 0,034 | sim |
| smells_delta | P-1 vs P0 | 30 | 0.033 | <0,001 | sim |

## Sensibilidade — Modelo linear de efeitos mistos

Efeito da persona (ref. = P0) sobre cada Δ, usando todas as réplicas válidas com intercepto aleatório por trecho:

- *complexidade ciclomática*: P-1: β=-0.091 (p=0,810), P1: β=-0.432 (p=0,253), P2: β=0.5481 (p=0,148)
- *code smells*: P-1: β=-0.0195 (p=0,943), P1: β=-0.189 (p=0,491), P2: β=-0.4508 (p=0,102)

## Custo (tokens)

| Persona | Tokens médios |
|---|---|
| P-1 | 8295 |
| P0 | 7646 |
| P1 | 7638 |
| P2 | 7783 |

## Conclusão

**Efeito placebo NÃO totalmente sustentado**: leitura conjunta de (i) validade semelhante entre personas (Cochran Q), (ii) ausência de diferença significativa de qualidade (Friedman) e (iii) equivalência TOST não confirmada nos contrastes-chave (controle/neutro vs especializada).

### Limitações

- Nível **P3 contextual** não executado (4 de 5 níveis propostos).
- **2 réplicas** por condição (proposta previa 5) → 240 obs (não 750).
- Análise estática em **escopo single-file** não capta classes novas (ex.: ReplaceConditionalWithPolymorphism infla a complexidade do arquivo alterado).
- **ReplaceMagicNumber** com Δ≈0 nas métricas → baixa sensibilidade de construto.
- Um único modelo LLM (Gemini Flash) e linguagem (Java).
