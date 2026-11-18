"""Fold the NATIVE sequences with ESMFold to get a baseline for the designs.

design_and_score.py reports absolute pLDDT and TM-score for designed sequences. Some
structures are hard for ESMFold even in their native form (peptides, chains of large
complexes, ...), so a low design pLDDT there says little about ProteinMPNN. This script
folds each native chain the same way the designs were folded (same chain, same TM-score
against the input structure) and writes one row per structure, so designs can be judged
relative to their native (delta pLDDT, as in Wittmann et al.).

Example (same machine/GPU environment as design_and_score.py):
    python fold_native_baseline.py pdb_results/structures

The chain per structure is taken from protein_mpnn_results/designs_all.csv when it has the
structure (so it matches the designs), else the PDB's first chain.

Output: protein_mpnn_results/native_baseline.csv
    pdb, chain, length, native_plddt, native_tm_score
"""

import argparse
import csv
import traceback
from pathlib import Path

from design_and_score import (
    ca_from_pdb_text,
    collect_pdbs,
    first_chain,
    fold_esmfold,
    load_esmfold,
    read_ca_chain,
    tm_score,
)

HERE = Path(__file__).resolve().parent


def design_chains(designs_csv: Path) -> dict[str, str]:
    """pdb -> chain that was designed, from designs_all.csv."""
    if not designs_csv.exists():
        return {}
    with designs_csv.open() as f:
        return {row["pdb"]: row["chain"] for row in csv.DictReader(f)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("pdbs", type=Path, nargs="+", help="PDB file(s) and/or folder(s) of .pdb files")
    parser.add_argument("--designs-csv", type=Path, default=HERE / "protein_mpnn_results/designs_all.csv")
    parser.add_argument("--output", type=Path, default=HERE / "protein_mpnn_results/native_baseline.csv")
    parser.add_argument("--chunk-size", type=int, default=None, help="ESMFold chunk size (e.g. 64) for long chains; less GPU memory, slower")
    parser.add_argument("--all", action="store_true", help="also fold structures that are not in designs_all.csv")
    args = parser.parse_args()

    chains = design_chains(args.designs_csv)
    pdbs = [p for p in collect_pdbs(args.pdbs) if args.all or p.stem in chains]
    if not pdbs:
        raise SystemExit("No matching .pdb files found")

    import torch

    model = load_esmfold(args.chunk_size)
    rows, failures = [], []
    for i, pdb in enumerate(pdbs, 1):
        chain = chains.get(pdb.stem) or first_chain(pdb)
        try:
            native_ca, native_seq = read_ca_chain(pdb, chain)
            with torch.no_grad():
                pdb_text, plddt = fold_esmfold(native_seq, model)
            tm = tm_score(ca_from_pdb_text(pdb_text), native_seq, native_ca, native_seq)
        except Exception as exc:  # keep the batch going
            traceback.print_exc()
            failures.append((pdb.name, f"{type(exc).__name__}: {exc}"))
            torch.cuda.empty_cache()
            continue
        rows.append(
            dict(pdb=pdb.stem, chain=chain, length=len(native_seq), native_plddt=round(plddt, 2), native_tm_score=round(tm, 4))
        )
        print(f"[{i}/{len(pdbs)}] {pdb.stem} chain {chain}: native pLDDT {plddt:.1f}, TM-score {tm:.3f}")
        # rewritten after every structure so an interrupted run keeps its results
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

    for name, reason in failures:
        print(f"FAILED {name}: {reason}")
    print(f"Wrote {args.output} ({len(rows)} structures)")


if __name__ == "__main__":
    main()
