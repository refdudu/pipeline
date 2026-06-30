# Seções de Resultados — Persona Prompting em Refatoração de Código

> Rascunho das seções a serem incorporadas ao artigo (converte a proposta em estudo
> executado). Texto em português acadêmico; valores extraídos de `results/r/`
> (análise canônica em R: Friedman/Nemenyi/Cochran Q em R base; ARTool, TOSTER,
> lme4/lmerTest). Figuras citadas estão em `results/r/*.png`.

---

## Nota sobre o protocolo executado (ajustes ao planejado)

A execução diferiu da proposta inicial em três pontos, que devem ser refletidos no
Abstract, na Seção V e na Tabela I:

1. **Quatro níveis de persona** foram executados — P-1 (controle negativo, padeiro),
   P0 (neutro), P1 (genérico) e P2 (especializado). O nível **P3 (contextual) não foi
   executado** e passa a constar como trabalho futuro.
2. **Duas réplicas** por condição (em vez de cinco), totalizando **240 observações**
   (30 trechos × 4 personas × 2 réplicas) em vez de 750.
3. A complexidade ciclomática de McCabe e a contagem de *code smells* mostraram-se
   **insensíveis à refatoração Substituir Número Mágico** (delta identicamente nulo).
   Para esse tipo adotou-se, como medida específica, o **delta de violações da regra
   `java:S109` ("Magic numbers should not be used")** do SonarQube, coletado em uma
   varredura estática adicional (regra de código-fonte, que dispensa recompilação).

O modelo foi um Gemini de *tier* Flash, e a preservação de comportamento foi verificada
pela suíte de testes do Defects4J (projeto Commons Math, versão corrigida).

---

## VI. RESULTADOS

### A. Visão geral e taxa de validade

Das 240 refatorações geradas, **236 (98,3%) preservaram o comportamento** (compilação
e suíte de testes íntegras); 4 falharam na compilação. A Tabela R-I apresenta a taxa de
validade por persona.

**Tabela R-I — Taxa de preservação de comportamento por persona.**

| Persona | Descrição | Válidas/Total | Taxa |
|---|---|---:|---:|
| P-1 | controle negativo (padeiro) | 59/60 | 0,983 |
| P0  | neutro (sem persona)        | 60/60 | 1,000 |
| P1  | genérico (*developer*)      | 59/60 | 0,983 |
| P2  | especializado (arquiteto)   | 58/60 | 0,967 |

O Teste Q de Cochran não detectou diferença significativa de validade entre as personas
(**Q = 2,40; p = 0,494**). Notavelmente, a persona de controle negativo (padeiro) preservou
o comportamento com a mesma frequência das demais, o que indica que a especialização da
identidade não reduziu o atrito operacional. (Figura: `validity.png`.)

A análise de qualidade foi conduzida em regime *complete-case* sobre os 20 trechos dos
tipos mensuráveis por McCabe/*code smells* (Extrair Método e Substituir Condicional por
Polimorfismo); o Substituir Número Mágico é analisado à parte na Subseção F.

### B. RQ1 — Efeito principal da persona na qualidade

O Teste de Friedman **não encontrou efeito significativo da persona** em nenhuma das duas
métricas (Tabela R-II). As distribuições dos deltas por persona são amplamente sobrepostas
(Figuras `box_complexity_delta.png` e `box_smells_delta.png`).

**Tabela R-II — Teste de Friedman (efeito principal da persona; *complete-case*, n = 20).**

| Métrica | χ² | p | Significativo (α=0,05)? |
|---|---:|---:|:--:|
| Complexidade ciclomática (Δ) | 4,017 | 0,260 | não |
| *Code smells* (Δ)            | 2,023 | 0,568 | não |

### C. Estatística descritiva por tipo de refatoração

A Tabela R-III mostra a mediana dos deltas por persona e tipo. Os valores são positivos em
complexidade (a refatoração tende a **aumentar** a complexidade ciclomática medida no escopo
do arquivo) e próximos de zero em *code smells*, sem padrão consistente de vantagem para as
personas especializadas — a persona P2 inclusive apresenta o maior aumento de complexidade
em ambos os tipos.

**Tabela R-III — Mediana dos deltas por persona × tipo de refatoração.**

| Tipo | Persona | Δ Complexidade (mediana) | Δ *Smells* (mediana) |
|---|:--:|---:|---:|
| Extrair Método | P-1 | 4,75 | −0,50 |
| Extrair Método | P0  | 5,25 | −1,00 |
| Extrair Método | P1  | 5,25 | −0,75 |
| Extrair Método | P2  | 6,00 | −1,00 |
| Subst. Cond. por Polimorfismo | P-1 | 8,25 | 0,25 |
| Subst. Cond. por Polimorfismo | P0  | 7,00 | 0,25 |
| Subst. Cond. por Polimorfismo | P1  | 6,50 | 0,00 |
| Subst. Cond. por Polimorfismo | P2  | 9,75 | 0,00 |

### D. RQ2 — Interação persona × tipo de refatoração

A ANOVA sobre postos alinhados (ART, via *ARTool*) **não indicou interação significativa**
entre persona e tipo de refatoração (Tabela R-IV). (Observação metodológica: a inclusão do
tipo Substituir Número Mágico — com delta identicamente nulo nessas métricas — produzia uma
interação espúria; sua exclusão dos testes de complexidade/*smells* eliminou o artefato.)

**Tabela R-IV — ART ANOVA, termo de interação persona × tipo (exploratória).**

| Métrica | F | p | Significativo? |
|---|---:|---:|:--:|
| Complexidade ciclomática (Δ) | 1,587 | 0,203 | não |
| *Code smells* (Δ)            | 0,462 | 0,710 | não |

### E. Teste de equivalência (TOST, margem ±1)

O TOST (via *TOSTER*) foi aplicado aos contrastes que definem a leitura de placebo
(Tabela R-V). A equivalência **não foi confirmada** entre o controle/neutro e a persona
especializada (P-1 vs P2 e P0 vs P2), embora as diferenças médias sejam pequenas: o
intervalo de confiança não cai integralmente dentro da margem ±1, refletindo poder
estatístico limitado (n = 20). Apenas o par P-1 vs P0 atingiu equivalência formal.

**Tabela R-V — TOST (equivalência prática, margem ±1).**

| Métrica | Contraste | n | Δ médio | p (TOST) | Equivalente? |
|---|---|:--:|---:|---:|:--:|
| Complexidade | P-1 vs P2 | 20 | −0,975 | 0,483 | não |
| Complexidade | P0 vs P2  | 20 | −0,850 | 0,413 | não |
| Complexidade | P-1 vs P0 | 20 | −0,125 | 0,041 | sim |
| *Smells*     | P-1 vs P2 | 20 | 0,675  | 0,224 | não |
| *Smells*     | P0 vs P2  | 20 | 0,625  | 0,212 | não |
| *Smells*     | P-1 vs P0 | 20 | 0,050  | 0,005 | sim |

### F. Substituir Número Mágico — remoção via `java:S109`

Como McCabe e *code smells* são cegos a essa refatoração (Δ ≡ 0), mediu-se o delta de
violações da regra `java:S109`. A linha de base somava **12 violações** nos 10 trechos.
**Todas as personas removeram os números mágicos de forma idêntica** (Tabela R-VI): o delta
médio é −1,2 para os quatro níveis, e em cada trecho as quatro personas produzem o mesmo
resultado (empate completo), de modo que o Teste de Friedman é inaplicável por ausência de
variação — isto é, **não há efeito de persona**, agora medido positivamente. A persona de
controle negativo (padeiro) removeu tantos números mágicos quanto a especializada.

**Tabela R-VI — Δ de violações `java:S109` por persona (média).**

| Persona | Δ S109 (média) |
|---|---:|
| P-1 | −1,2 |
| P0  | −1,2 |
| P1  | −1,2 |
| P2  | −1,2 |

### G. Análise de sensibilidade (modelo de efeitos mistos)

Um modelo linear de efeitos mistos (lme4/lmerTest; intercepto aleatório por trecho,
referência = P0, todas as réplicas válidas) **confirma a ausência de efeito**: nenhum
coeficiente de persona foi significativo (Tabela R-VII), corroborando o Friedman e
indicando que a conclusão não é artefato do regime *complete-case*.

**Tabela R-VII — LMM: efeito da persona (ref. P0).**

| Métrica | Persona | β | p |
|---|:--:|---:|---:|
| Complexidade | P-1 | −0,136 | 0,811 |
| Complexidade | P1  | −0,649 | 0,255 |
| Complexidade | P2  | 0,843  | 0,143 |
| *Smells*     | P-1 | −0,027 | 0,948 |
| *Smells*     | P1  | −0,284 | 0,494 |
| *Smells*     | P2  | −0,684 | 0,103 |

### H. Custo em tokens

O custo médio em tokens por execução foi semelhante entre as personas (Tabela R-VIII);
a persona de controle negativo (padeiro) foi inclusive a **mais cara**, sem qualquer
ganho de qualidade. (Figura: `tokens.png`.)

**Tabela R-VIII — Tokens médios por persona.**

| Persona | Tokens (média) |
|---|---:|
| P-1 | 8.295 |
| P0  | 7.646 |
| P1  | 7.638 |
| P2  | 7.783 |

---

## VII. DISCUSSÃO

**RQ1 — A persona não melhora a qualidade.** Em complexidade ciclomática e em *code smells*,
o Friedman não rejeitou *H*1₀ (p = 0,260 e p = 0,568), e o modelo de efeitos mistos
confirmou o resultado (todos os p > 0,10). A hipótese *H*1₁ — de que personas especializadas
gerariam código estruturalmente melhor — **não foi sustentada**; ao contrário, a persona
especializada (P2) exibiu, na descritiva, o maior aumento de complexidade. A medição
específica por `java:S109` reforça o achado no único tipo em que as métricas gerais eram
cegas: todas as personas, **inclusive o controle negativo**, removeram os números mágicos de
forma idêntica.

**RQ2 — Não há interação persona × tipo.** A ART não indicou interação significativa
(p = 0,203 e p = 0,710), de modo que *H*2₁ não se sustenta: a especialização da persona não
passou a "valer mais" nas tarefas estruturalmente mais complexas.

**Veredito de placebo: evidência consistente, equivalência formal pendente.** O protocolo
(Seção IV) condiciona o veredito de placebo a duas frentes: (i) taxas de validade semelhantes
entre as personas e (ii) equivalência nas métricas de qualidade. A primeira frente foi
plenamente atendida (Cochran Q n.s.; padeiro tão válido quanto arquiteto). A segunda, porém,
**não foi positivamente confirmada**: embora não haja diferença significativa (Friedman) e as
diferenças médias sejam pequenas, o TOST a ±1 não estabeleceu equivalência entre controle/neutro
e a persona especializada, em razão do poder estatístico limitado (n = 20 após a redução de
réplicas). A leitura honesta, portanto, é de **ausência de benefício detectável** — coerente
com a hipótese de placebo —, sem que se possa, ainda, afirmar equivalência estatística formal.
Esse ponto motiva diretamente a replicação com mais réplicas e o nível P3.

**Contribuição metodológica: a métrica importa.** Dois fenômenos merecem registro. Primeiro,
a complexidade ciclomática **no escopo do arquivo** tende a aumentar com Extrair Método
(cada método extraído adiciona um ponto de entrada) e com Substituir Condicional por
Polimorfismo (cujas classes novas escapam ao escopo de arquivo único analisado), de modo que
deltas positivos não significam, necessariamente, piora real. Segundo, a contagem total de
*code smells* é **cega** ao Substituir Número Mágico, exigindo a regra específica `java:S109`.
Em conjunto, esses pontos constituem um alerta prático: a avaliação de refatorações geradas
por LLM exige métricas alinhadas ao tipo de transformação, sob pena de medir artefato.

**Relevância prática.** Como o custo em tokens é semelhante entre os níveis — e a persona
incongruente foi a mais cara —, não há justificativa empírica, neste contexto, para investir
esforço na elaboração de preâmbulos de persona especializados em tarefas de refatoração:
o retorno em qualidade é nulo ou indetectável e o custo não diminui.

---

## VIII. CONCLUSÃO

Este estudo isolou e quantificou o efeito do nível de especialização da persona sobre a
qualidade estrutural de refatorações geradas por um LLM, em 240 refatorações de código Java
do Defects4J. Os resultados são consistentes em todas as frentes: o Teste de Friedman, o
modelo de efeitos mistos, a ART e o Teste Q de Cochran **não detectaram qualquer efeito da
persona** sobre a qualidade ou sobre a preservação de comportamento, e a remoção de números
mágicos (medida por `java:S109`) foi idêntica entre todas as personas — incluindo o controle
negativo. Não se confirmou, portanto, a hipótese de que personas especializadas produzem
código melhor; o comportamento observado é **compatível com um efeito placebo** da engenharia
de prompts, ressalvado que a equivalência estatística formal (TOST) não foi estabelecida sob o
poder amostral disponível. A implicação prática é direta: em refatoração de código, elaborar
prompts de persona especializados não trouxe ganho mensurável a um custo de tokens equivalente.
Como trabalhos futuros, destacam-se a execução do nível contextual (P3), o aumento do número
de réplicas para viabilizar a conclusão de equivalência, e a adoção de métricas por método e
multi-arquivo, alinhadas a cada tipo de refatoração.

---

## Atualizações necessárias nas seções existentes

- **Abstract / Seção V-B / Tabela I:** 5×3 → 4×3 executado; 750 → 240 observações; 5 → 2 réplicas.
- **Seção IV (Hipóteses):** manter; relatar *H*1₀ não rejeitada e *H*1₁/*H*2₁ não sustentadas.
- **Seção V-C (Variáveis dependentes):** acrescentar o Δ `java:S109` como medida específica do
  Substituir Número Mágico.
- **Seção VI (Ameaças à Validade):** acrescentar (a) cegueira das métricas gerais ao número
  mágico, contornada por `java:S109`; (b) inflação da complexidade de arquivo por Extrair
  Método / escopo de arquivo único; (c) poder limitado do TOST (n = 20) para conclusão de
  equivalência.
- **Seção VII (Plano/Cronograma):** substituída pela Seção de Resultados acima; remanescentes
  (P3, mais réplicas) movidos para Trabalhos Futuros.
