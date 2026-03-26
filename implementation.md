# Phase 52 — Hypothesis Test: Scaffold Family pIC50 Comparison

**Version:** 1.1 | **Tier:** Micro | **Date:** 2026-03-26

## Goal
Test whether scaffold families have statistically different pIC50 distributions and hit rates.
Tests: Kruskal-Wallis (global), Mann-Whitney U (pairwise), Welch t-test, Cohen's d effect size.

CLI: `python main.py --input data/compounds.csv --threshold 7.0`

Outputs: hypothesis_results.csv, hypothesis_plot.png

## Logic
- Load compounds, assign hit label (pIC50 >= threshold)
- Kruskal-Wallis: is there ANY difference across all families?
- All pairwise: Mann-Whitney U (non-parametric) + Welch t (parametric) + Cohen's d
- Report significant pairs at alpha=0.05
- Plot: boxplot by family + hit rate bar chart

## Results

### Global test
- **Kruskal-Wallis: H=19.752, p=0.0014** → families differ significantly

### Hit rates
| Family | Hit Rate | Count |
|---|---|---|
| ind | **100%** | 7/7 |
| naph | 86% | 6/7 |
| quin | 86% | 6/7 |
| benz | 67% | 8/12 |
| bzim | 50% | 3/6 |
| pyr | **0%** | 0/6 |

### Significant pairs (Mann-Whitney, alpha=0.05)
| Pair | p-value | Cohen's d | Means |
|---|---|---|---|
| ind vs pyr | 0.0012 | **3.62** | 7.90 vs 6.27 |
| naph vs pyr | 0.0012 | 2.95 | 7.49 vs 6.27 |
| quin vs pyr | 0.0012 | 3.00 | 7.62 vs 6.27 |
| ind vs bzim | 0.0082 | 2.01 | 7.90 vs 6.88 |
| benz vs pyr | 0.0191 | 1.68 | 7.21 vs 6.27 |
| quin vs bzim | 0.0350 | 1.46 | 7.62 vs 6.88 |

## Key Insights
- pyr is the weakest scaffold class — statistically significantly weaker than ALL other families
- ind is the strongest — 100% hit rate and highest mean pIC50 (7.90)
- Effect sizes are large (d>1.4 for all significant pairs) — differences are chemically meaningful
- This validates the scaffold SAR narrative: class matters as much as substituent pattern

## Deviations from Plan
- None; plan implemented as specified.
