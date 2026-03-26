import sys
if sys.platform == "win32":
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")

import argparse, os, warnings
warnings.filterwarnings("ignore")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from scipy import stats
from itertools import combinations
RDLogger.DisableLog("rdApp.*")

FAMILY_COLORS = {"benz": "#4C72B0", "naph": "#DD8452", "ind": "#55A868",
                 "quin": "#C44E52", "pyr": "#8172B2", "bzim": "#937860", "other": "#808080"}


def load_compounds(path, threshold):
    df = pd.read_csv(path)
    records = []
    for _, row in df.iterrows():
        mol = Chem.MolFromSmiles(str(row["smiles"]))
        if mol is None:
            continue
        try:
            pic50 = float(row["pic50"])
        except (KeyError, ValueError):
            continue
        if np.isnan(pic50):
            continue
        fam = str(row["compound_name"]).split("_")[0]
        records.append({
            "compound_name": str(row["compound_name"]),
            "family": fam if fam in FAMILY_COLORS else "other",
            "pic50": pic50,
            "is_hit": int(pic50 >= threshold),
        })
    return pd.DataFrame(records)


def compare_families(df, fam_a, fam_b):
    """All pairwise statistical tests between two scaffold families."""
    a = df.loc[df["family"] == fam_a, "pic50"].values
    b = df.loc[df["family"] == fam_b, "pic50"].values

    results = {}
    # Mann-Whitney U (non-parametric, does not assume normality)
    if len(a) >= 3 and len(b) >= 3:
        stat_mw, p_mw = stats.mannwhitneyu(a, b, alternative="two-sided")
        results["mannwhitney_U"] = round(float(stat_mw), 3)
        results["mannwhitney_p"] = round(float(p_mw), 4)
    else:
        results["mannwhitney_U"] = None
        results["mannwhitney_p"] = None

    # Welch t-test (parametric, unequal variance)
    if len(a) >= 2 and len(b) >= 2:
        stat_t, p_t = stats.ttest_ind(a, b, equal_var=False)
        results["welch_t"] = round(float(stat_t), 3)
        results["welch_p"] = round(float(p_t), 4)
    else:
        results["welch_t"] = None
        results["welch_p"] = None

    # Effect size: Cohen's d
    pooled_std = np.sqrt((np.std(a, ddof=1) ** 2 + np.std(b, ddof=1) ** 2) / 2) if len(a) >= 2 and len(b) >= 2 else 0
    cohen_d = (np.mean(a) - np.mean(b)) / pooled_std if pooled_std > 0 else np.nan
    results["cohen_d"] = round(float(cohen_d), 3) if not np.isnan(cohen_d) else None

    results.update({
        "family_a": fam_a, "family_b": fam_b,
        "n_a": len(a), "n_b": len(b),
        "mean_a": round(float(np.mean(a)), 3) if len(a) > 0 else None,
        "mean_b": round(float(np.mean(b)), 3) if len(b) > 0 else None,
    })
    return results


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input", required=True)
    parser.add_argument("--threshold", type=float, default=7.0, help="pIC50 hit threshold")
    parser.add_argument("--alpha", type=float, default=0.05, help="Significance level")
    parser.add_argument("--output-dir", default="output")
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    print(f"\nLoading: {args.input}  (hit threshold: pIC50 >= {args.threshold})")
    df = load_compounds(args.input, args.threshold)

    families = [f for f in FAMILY_COLORS if f in df["family"].values and f != "other"]
    family_sizes = {f: len(df[df["family"] == f]) for f in families}
    print(f"  Families: {family_sizes}")

    # All pairwise comparisons
    pairs = list(combinations(families, 2))
    rows = []
    for fam_a, fam_b in pairs:
        r = compare_families(df, fam_a, fam_b)
        r["significant_mw"] = r["mannwhitney_p"] is not None and r["mannwhitney_p"] < args.alpha
        rows.append(r)

    res_df = pd.DataFrame(rows)
    res_df.to_csv(os.path.join(args.output_dir, "hypothesis_results.csv"), index=False)
    print(f"Saved: {args.output_dir}/hypothesis_results.csv")

    # Kruskal-Wallis: is there ANY difference across families?
    groups = [df.loc[df["family"] == f, "pic50"].values for f in families if len(df[df["family"] == f]) >= 2]
    h_stat, p_kw = stats.kruskal(*groups)
    print(f"\n  Kruskal-Wallis (all families): H={h_stat:.3f}, p={p_kw:.4f} ({'significant' if p_kw < args.alpha else 'not significant'})")

    # Hit rate per family
    hit_rates = {}
    for f in families:
        sub = df[df["family"] == f]
        hr = sub["is_hit"].mean()
        hit_rates[f] = hr

    # Plot: boxplot + hit rate
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Boxplot
    data_list = [df.loc[df["family"] == f, "pic50"].values for f in families]
    bp = ax1.boxplot(data_list, patch_artist=True, notch=False)
    for patch, fam in zip(bp["boxes"], families):
        patch.set_facecolor(FAMILY_COLORS[fam])
        patch.set_alpha(0.8)
    ax1.axhline(args.threshold, color="red", ls="--", lw=1.2, label=f"Hit threshold ({args.threshold})")
    ax1.set_xticks(range(1, len(families) + 1))
    ax1.set_xticklabels(families, fontsize=10)
    ax1.set_ylabel("pIC50", fontsize=11)
    ax1.set_title(f"pIC50 by Scaffold Family\n(KW H={h_stat:.2f}, p={p_kw:.4f})", fontsize=11, fontweight="bold")
    ax1.legend(fontsize=9)
    ax1.spines["top"].set_visible(False); ax1.spines["right"].set_visible(False)

    # Hit rate bar chart
    colors_bar = [FAMILY_COLORS[f] for f in families]
    rates = [hit_rates[f] for f in families]
    bars = ax2.bar(families, rates, color=colors_bar, edgecolor="white", alpha=0.85)
    for bar, rate in zip(bars, rates):
        ax2.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01, f"{rate:.0%}",
                 ha="center", fontsize=9)
    ax2.set_xlabel("Scaffold Family", fontsize=11)
    ax2.set_ylabel(f"Hit Rate (pIC50 >= {args.threshold})", fontsize=11)
    ax2.set_title("Hit Rate by Scaffold Family", fontsize=11, fontweight="bold")
    ax2.set_ylim(0, 1.1)
    ax2.spines["top"].set_visible(False); ax2.spines["right"].set_visible(False)

    plt.suptitle("Hypothesis Test: Scaffold Family pIC50 Comparison", fontsize=12, fontweight="bold")
    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, "hypothesis_plot.png"), dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {args.output_dir}/hypothesis_plot.png")

    # Summary: significant pairs
    sig = res_df[res_df["significant_mw"] == True]
    print(f"\n--- Significant Pairs (Mann-Whitney, alpha={args.alpha}) ---")
    if len(sig) == 0:
        print("  None")
    else:
        for _, row in sig.iterrows():
            print(f"  {row['family_a']:5s} vs {row['family_b']:5s}  p={row['mannwhitney_p']:.4f}  d={row['cohen_d']:.2f}  means={row['mean_a']:.2f} vs {row['mean_b']:.2f}")

    print(f"\n--- Hit Rates ---")
    for f in families:
        print(f"  {f:5s}: {hit_rates[f]:.0%} ({int(df[df['family']==f]['is_hit'].sum())}/{family_sizes[f]})")

    print("\nDone.")


if __name__ == "__main__":
    main()
