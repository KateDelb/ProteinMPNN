"""3D scatter of ProteinMPNN design quality: TM-score, pLDDT and MPNN score.

Reads designs_all.csv (columns pdb, chain, design, mpnn_score, plddt, tm_score,
sequence) and draws one point per designed variant. The pLDDT axis shows delta pLDDT (design minus native);
green points are the designs marked `included` in the CSV (TM > 0.5, delta pLDDT > -10,
pLDDT > 70; the cutoffs are set at the top of this script and written to the CSV on every run).
Saves a PNG to plots/ and optionally opens an interactive window (--show).

If native_baseline*.csv files (from fold_native_baseline.py) sit next to designs_all.csv, the
script also adds native_plddt, native_tm_score and delta_plddt (design pLDDT minus the pLDDT of
the native sequence folded the same way) to the CSV, and draws density plots of
TM-score and delta pLDDT, like panels B and C of Wittmann et al. Structures without a baseline
yet have empty delta values and are left out of those plots.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

DEFAULT_CSV = Path(__file__).resolve().parent / "protein_mpnn_results/designs_all.csv"
# selection cutoffs: a design is selected (green in the plots, `included` in the CSV) when all three hold
MIN_PLDDT, MIN_TM = 70, 0.5
MIN_DELTA_PLDDT = -10
PASS_COLOR, PASS_EDGE, FAIL_COLOR = "#2ca02c", "#1b6b1b", "#bdbdbd"
PLOTS_DIR = Path(__file__).resolve().parent / "protein_mpnn_results/plots"


def selected_mask(df: pd.DataFrame) -> pd.Series:
    """The designs marked green in the plots: TM > cutoff, delta pLDDT > cutoff and pLDDT > cutoff."""
    return (df["tm_score"] > MIN_TM) & (df["delta_plddt"] > MIN_DELTA_PLDDT) & (df["plddt"] > MIN_PLDDT)


SELECTION_LABEL = f"selected (TM>{MIN_TM}, \u0394pLDDT>{MIN_DELTA_PLDDT}, pLDDT>{MIN_PLDDT})"


def add_native(df: pd.DataFrame, csv_path: Path) -> pd.DataFrame:
    """Join native pLDDT/TM-score (native_baseline*.csv next to the CSV) and compute delta_plddt."""
    df = df.drop(columns=["native_plddt", "native_tm_score", "delta_plddt"], errors="ignore")
    files = sorted(csv_path.parent.glob("native_baseline*.csv"))
    if not files:
        return df
    native = pd.concat([pd.read_csv(f) for f in files]).drop_duplicates(["pdb", "chain"], keep="last")
    df = df.merge(native[["pdb", "chain", "native_plddt", "native_tm_score"]], on=["pdb", "chain"], how="left")
    df["delta_plddt"] = (df["plddt"] - df["native_plddt"]).round(2)
    return df


def annotate_csv(csv_path: Path) -> None:
    """Add (or refresh) an `included` column (= the selected designs, green in the plots) and the native-baseline columns."""
    df = add_native(pd.read_csv(csv_path), csv_path)
    df["included"] = selected_mask(df)
    df.to_csv(csv_path, index=False)
    print(f"Wrote 'included' to {csv_path} ({df['included'].sum()} of {len(df)} True)")
    if "delta_plddt" in df:
        print(f"  delta_plddt for {df['delta_plddt'].notna().sum()} of {len(df)} designs "
              f"({df.loc[df['delta_plddt'].isna(), 'pdb'].nunique()} structures still without a native baseline)")


def plot_density(df: pd.DataFrame, output: Path) -> None:
    """Density of TM-score and delta pLDDT over all designs (panels B and C of the paper)."""
    df = df.dropna(subset=["tm_score", "delta_plddt"])
    fig, axes = plt.subplots(2, 1, figsize=(7, 6.5))
    for ax, col, label, lo, hi in [
        (axes[0], "tm_score", "TM-score", 0, 1),
        (axes[1], "delta_plddt", "\u0394pLDDT (design \u2212 native)", df["delta_plddt"].min() - 10, df["delta_plddt"].max() + 10),
    ]:
        x = np.linspace(lo, hi, 400)
        ax.fill_between(x, gaussian_kde(df[col])(x), color=PASS_COLOR, alpha=0.3)
        ax.plot(x, gaussian_kde(df[col])(x), color=PASS_EDGE)
        ax.set_xlabel(label)
        ax.set_ylabel("Kernel density estimate")
        ax.set_yticks([])
    axes[1].axvline(0, color="#555555", lw=0.8, ls="--")
    axes[1].xaxis.set_major_locator(MultipleLocator(10))
    axes[0].set_title(f"ProteinMPNN designs (n={len(df)}, {df['pdb'].nunique()} structures)")
    fig.tight_layout()
    fig.savefig(output, dpi=200, bbox_inches="tight")
    print(f"Saved {output}")
    plt.close(fig)


def plot_tm_vs_delta(df: pd.DataFrame, output: Path) -> None:
    """2D scatter of TM-score against delta pLDDT (the 3D plot without the MPNN score)."""
    df = df.dropna(subset=["tm_score", "delta_plddt"])
    passed = df["included"]
    fig, ax = plt.subplots(figsize=(7, 5.5))
    for mask, color, edge, label in [
        (~passed, FAIL_COLOR, "#8c8c8c", "not selected"),
        (passed, PASS_COLOR, PASS_EDGE, SELECTION_LABEL),
    ]:
        sub = df[mask]
        ax.scatter(sub["tm_score"], sub["delta_plddt"], color=color, s=22, alpha=0.8, edgecolor=edge,
                   linewidth=0.3, label=f"{label} (n={len(sub)})")
    ax.axvline(MIN_TM, color="#555555", lw=0.8, ls="--")
    ax.axhline(0, color="#555555", lw=0.8, ls="--")
    ax.yaxis.set_major_locator(MultipleLocator(10))
    ax.set_xlabel("TM-score")
    ax.set_ylabel("\u0394pLDDT (design \u2212 native)")
    ax.set_title(f"TM-score vs \u0394pLDDT of ProteinMPNN designs (n={len(df)})")
    ax.legend(loc="lower left", frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(output, dpi=200, bbox_inches="tight")
    print(f"Saved {output}")
    plt.close(fig)


def plot_tm_delta_plddt(df: pd.DataFrame, output: Path) -> None:
    """3D scatter of TM-score, delta pLDDT and pLDDT, with the box TM > cutoff, delta pLDDT > cutoff, pLDDT > cutoff."""
    df = df.dropna(subset=["tm_score", "delta_plddt", "plddt"])
    passed = df["included"]
    label_pass = SELECTION_LABEL
    print(f"{passed.sum()} of {len(df)} designs pass the cutoffs")
    print(f"  TM>{MIN_TM}: {(df['tm_score'] > MIN_TM).sum()}, \u0394pLDDT>{MIN_DELTA_PLDDT}: "
          f"{(df['delta_plddt'] > MIN_DELTA_PLDDT).sum()}, pLDDT>{MIN_PLDDT}: {(df['plddt'] > MIN_PLDDT).sum()}")

    fig = plt.figure(figsize=(8, 6.5))
    ax = fig.add_subplot(projection="3d")
    for mask, color, edge, label in [
        (~passed, FAIL_COLOR, "#8c8c8c", "not selected"),
        (passed, PASS_COLOR, PASS_EDGE, label_pass),
    ]:
        sub = df[mask]
        ax.scatter(sub["tm_score"], sub["delta_plddt"], sub["plddt"], color=color, s=22, alpha=0.85,
                   edgecolor=edge, linewidth=0.3, label=f"{label} (n={len(sub)})")

    x0, x1 = MIN_TM, 1.0
    y0, y1 = MIN_DELTA_PLDDT, df["delta_plddt"].max()
    z0, z1 = MIN_PLDDT, 100
    c = [(x, y, z) for x in (x0, x1) for y in (y0, y1) for z in (z0, z1)]
    faces = [
        [c[0], c[1], c[3], c[2]], [c[4], c[5], c[7], c[6]],
        [c[0], c[1], c[5], c[4]], [c[2], c[3], c[7], c[6]],
        [c[0], c[2], c[6], c[4]], [c[1], c[3], c[7], c[5]],
    ]
    ax.add_collection3d(Poly3DCollection(faces, facecolors=PASS_COLOR, edgecolors=PASS_EDGE, linewidths=0.8, alpha=0.08))

    ax.set_xlabel("TM-score")
    ax.set_ylabel("\u0394pLDDT (design \u2212 native)", labelpad=14)
    ax.yaxis.set_major_locator(MultipleLocator(10))
    ax.tick_params(axis="y", labelsize=6.5)
    ax.set_zlabel("pLDDT")
    ax.set_title(f"TM-score, \u0394pLDDT and pLDDT of ProteinMPNN designs (n={len(df)})")
    ax.legend(loc="upper left", frameon=False, fontsize=8)
    ax.view_init(elev=22, azim=-55)
    fig.savefig(output, dpi=200, bbox_inches="tight")
    print(f"Saved {output}")
    plt.close(fig)


def plot_folding_scores(csv_path: Path, output: Path, show: bool = False) -> None:
    df = add_native(pd.read_csv(csv_path), csv_path)
    if "delta_plddt" not in df:
        raise SystemExit("No native_baseline*.csv found next to the CSV; run fold_native_baseline.py first")
    df = df.dropna(subset=["tm_score", "delta_plddt", "mpnn_score"])

    fig = plt.figure(figsize=(8, 6.5))
    ax = fig.add_subplot(projection="3d")

    passed = df["included"]
    print(f"{passed.sum()} of {len(df)} designs {SELECTION_LABEL}")
    for mask, color, label in [
        (~passed, FAIL_COLOR, "not selected"),
        (passed, PASS_COLOR, SELECTION_LABEL),
    ]:
        sub = df[mask]
        ax.scatter(
            sub["tm_score"], sub["delta_plddt"], sub["mpnn_score"],
            color=color, s=22, alpha=0.85, edgecolor=PASS_EDGE if mask is passed else "#8c8c8c",
            linewidth=0.3, label=f"{label} (n={len(sub)})",
        )

    ax.set_xlabel("TM-score")
    ax.set_ylabel("\u0394pLDDT (design \u2212 native)")
    ax.yaxis.set_major_locator(MultipleLocator(10))
    ax.set_zlabel("MPNN score")
    ax.set_title(f"Folding scores of ProteinMPNN designs (n={len(df)})")
    ax.legend(loc="upper left", frameon=False, fontsize=8)
    ax.view_init(elev=22, azim=-55)

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=200, bbox_inches="tight")
    print(f"Saved {output}")
    if show:
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV, help="designs_all.csv path")
    parser.add_argument("--output", type=Path, default=PLOTS_DIR / "folding_scores_3d.png")
    parser.add_argument("--no-annotate", action="store_true",
                        help="plot only: do not refresh the 'included'/delta columns in the CSV first (they must already be up to date)")
    parser.add_argument("--show", action="store_true", help="open an interactive window")
    args = parser.parse_args()
    if not args.no_annotate:
        annotate_csv(args.csv)  # the plots colour by the CSV's `included` column, so refresh it first
    if "included" not in pd.read_csv(args.csv, nrows=1):
        raise SystemExit("designs_all.csv has no 'included' column; run without --no-annotate first")
    plot_folding_scores(args.csv, args.output, args.show)
    with_native = add_native(pd.read_csv(args.csv), args.csv)
    if "delta_plddt" in with_native and with_native["delta_plddt"].notna().any():
        plot_density(with_native, args.output.with_name("delta_plddt_density.png"))
        plot_tm_vs_delta(with_native, args.output.with_name("tm_vs_delta_plddt.png"))
        plot_tm_delta_plddt(with_native, args.output.with_name("folding_scores_3d_tm_delta_plddt.png"))
