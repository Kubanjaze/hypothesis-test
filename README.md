# hypothesis-test — Phase 52

Statistical comparison of pIC50 distributions across scaffold families. Uses Mann-Whitney U, Welch t-test, Cohen's d, and Kruskal-Wallis. Key finding: ind (100% hit rate) vs pyr (0% hit rate), Cohen's d=3.62.

## Usage

```bash
PYTHONUTF8=1 python main.py --input data/compounds.csv --threshold 7.0
```

## Outputs

- `output/hypothesis_results.csv` — all pairwise test results
- `output/hypothesis_plot.png` — boxplot by family + hit rate bar chart
