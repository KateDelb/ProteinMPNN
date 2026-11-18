"""Keep the best-scoring designs of every backbone.

ProteinMPNN scores depend on the backbone (a tightly packed helical bundle scores lower than
a loop-rich protein or an interface), so instead of one global MPNN threshold the usual
pipeline (e.g. RFdiffusion -> ProteinMPNN -> AlphaFold2/ESMFold) keeps the top-k lowest-scoring
sequences per backbone. This script does that for the table of designs:

    protein_mpnn_results/designs_all_full.csv   every design (kept as the archive, never modified here)
        -> protein_mpnn_results/designs_all.csv  the top-N designs per backbone (what plot_folding_scores.py reads)

A backbone is a pdb + chain. Rank 1 is the lowest (best) MPNN score; ties are broken by order in
the file. Designs without an MPNN score cannot be ranked and are left out of the top-N table. The
rank is written to an `mpnn_rank` column of the output.

    python select_top_designs.py            # top 3
    python select_top_designs.py --top-n 5
"""

import argparse
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "protein_mpnn_results"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", type=Path, default=RESULTS / "designs_all_full.csv")
    parser.add_argument("--output", type=Path, default=RESULTS / "designs_all.csv")
    parser.add_argument("--top-n", type=int, default=3, help="designs to keep per backbone (default 3)")
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    df["mpnn_rank"] = df.groupby(["pdb", "chain"])["mpnn_score"].rank(method="first", ascending=True)
    top = df[df["mpnn_rank"] <= args.top_n].copy()
    top["mpnn_rank"] = top["mpnn_rank"].astype(int)
    top.to_csv(args.output, index=False)
    print(f"{args.input.name}: {len(df)} designs, {df.groupby(['pdb', 'chain']).ngroups} backbones")
    print(f"Wrote {args.output} ({len(top)} designs, top {args.top_n} per backbone)")
    short = top.groupby(["pdb", "chain"]).size()
    if (short < args.top_n).any():
        print(f"  {(short < args.top_n).sum()} backbones have fewer than {args.top_n} rankable designs")
    if df["mpnn_score"].isna().any():
        print(f"  {df['mpnn_score'].isna().sum()} designs have no MPNN score and were not ranked")


if __name__ == "__main__":
    main()
