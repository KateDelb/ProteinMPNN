"""Design sequences for one or many PDBs with ProteinMPNN, fold them, and score them.

For each designed sequence this reports:
- mpnn_score: ProteinMPNN's own score (average negative log-probability of the designed residues; lower is more confident)
- plddt: mean pLDDT (0-100) of the structure ESMFold predicts for the designed sequence
- x_removed: placeholder X residues (missing/non-standard residues in the input) dropped from the sequence before folding
- tm_score: TM-score (TM-align, normalised by the input structure's length) between that predicted structure and the input PDB

Pipeline: protein_mpnn_run.py  ->  ESMFold (transformers, needs a GPU in practice)  ->  tmtools.tm_align

Examples:
    python design_and_score.py inputs/PDB_monomers/pdbs/5L33.pdb --chains A --num-seqs 4 --temperature 0.1
    python design_and_score.py pdb_results/structures --num-seqs 4   # every .pdb in a folder

ESMFold folds one chain, so each PDB is designed and scored on one chain: --chains if given,
else its first chain. .cif files are skipped (ProteinMPNN reads PDB only). One ESMFold load
serves the whole batch; a failing PDB is logged and skipped.

Output, always next to this script (override with --out-dir):
    results/design_<pdb name>/      designs.csv, folded/<design>.pdb, seqs/<pdb name>.fa
    results/designs_all.csv         every design of every PDB (updated after each PDB)
    results.tar.gz                  all of results/, to download

Use --skip-fold to only run ProteinMPNN (no structure prediction, no scores).
"""

import argparse
import csv
import shutil
import subprocess
import sys
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODEL_CHOICES = ["v_48_002", "v_48_010", "v_48_020", "v_48_030"]
ESMFOLD_MODEL = "facebook/esmfold_v1"


def run_proteinmpnn(args: argparse.Namespace, pdb: Path, chains: str, out_dir: Path) -> list[tuple[str, float, str]]:
    """Run protein_mpnn_run.py and return (name, score, sequence) for each designed sequence."""
    cmd = [
        sys.executable,
        str(HERE / "protein_mpnn_run.py"),
        "--pdb_path", str(pdb),
        "--out_folder", str(out_dir),
        "--num_seq_per_target", str(args.num_seqs),
        "--sampling_temp", str(args.temperature),
        "--seed", str(args.seed),
        "--backbone_noise", str(args.backbone_noise),
        "--model_name", args.model_name,
        "--batch_size", "1",
    ]
    if chains:
        cmd += ["--pdb_path_chains", chains]
    if args.omit_aas:
        cmd += ["--omit_AAs", args.omit_aas]
    if args.ca_only:
        cmd += ["--ca_only"]
    if args.use_soluble_model:
        cmd += ["--use_soluble_model"]
    print("Running ProteinMPNN:", " ".join(cmd[1:]))
    subprocess.run(cmd, check=True, cwd=HERE)

    fasta = out_dir / "seqs" / f"{pdb.stem}.fa"
    designs = []
    header, sequence = None, []
    for line in fasta.read_text().splitlines() + [">end"]:
        if line.startswith(">"):
            # the first record is the native sequence; designs carry sample=... in their header
            if header is not None and "sample=" in header:
                score = float(header.split("score=")[1].split(",")[0])
                sample = header.split("sample=")[1].split(",")[0]
                designs.append((f"design_{int(sample):03d}", score, "".join(sequence)))
            header, sequence = line, []
        else:
            sequence.append(line.strip())
    return designs


def read_ca_chain(pdb: Path, chain_id: str) -> tuple[list[tuple[float, float, float]], str]:
    """CA coordinates and one-letter sequence of one chain (standard residues only)."""
    from Bio.PDB import PDBParser
    from Bio.PDB.Polypeptide import three_to_index, index_to_one

    structure = PDBParser(QUIET=True).get_structure("input", str(pdb))
    chain = structure[0][chain_id]
    coords, sequence = [], []
    for residue in chain:
        resname = residue.get_resname()
        if resname == "MSE":  # selenomethionine: ProteinMPNN reads it as methionine, so do we
            resname = "MET"
        elif residue.id[0] != " ":
            continue
        if "CA" not in residue:
            continue
        try:
            letter = index_to_one(three_to_index(resname))
        except KeyError:  # non-standard residue (e.g. UNK): ProteinMPNN turns it into an X, removed below
            continue
        coords.append(tuple(residue["CA"].coord))
        sequence.append(letter)
    return coords, "".join(sequence)


def fold_esmfold(sequence: str, model) -> tuple[str, float]:
    """Predict a structure with ESMFold; returns (PDB text, mean pLDDT)."""
    pdb_text = model.infer_pdb(sequence)
    plddt = [float(line[60:66]) for line in pdb_text.splitlines() if line.startswith("ATOM") and line[12:16].strip() == "CA"]
    mean = sum(plddt) / len(plddt)
    # some transformers versions write pLDDT as 0-1 in the B-factor column; report 0-100
    return pdb_text, mean * 100 if max(plddt) <= 1.0 else mean


def ca_from_pdb_text(pdb_text: str) -> list[tuple[float, float, float]]:
    return [
        (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        for line in pdb_text.splitlines()
        if line.startswith("ATOM") and line[12:16].strip() == "CA"
    ]


def tm_score(pred_ca, pred_seq, native_ca, native_seq) -> float:
    """TM-score normalised by the input (native) length, via TM-align (tmtools)."""
    import numpy as np
    from tmtools import tm_align

    result = tm_align(np.array(pred_ca), np.array(native_ca), pred_seq, native_seq)
    return float(result.tm_norm_chain2)


def first_chain(pdb: Path) -> str:
    """Id of the first chain with ATOM records."""
    for line in pdb.read_text().splitlines():
        if line.startswith("ATOM"):
            return line[21]
    raise ValueError(f"No ATOM records in {pdb}")


def collect_pdbs(inputs: list[Path]) -> list[Path]:
    """Expand directories to their .pdb files; warn about and skip anything else."""
    pdbs = []
    for path in inputs:
        files = sorted(path.glob("*.pdb")) if path.is_dir() else [path]
        for f in files:
            if f.suffix.lower() == ".pdb":
                pdbs.append(f.resolve())
            else:
                print(f"Skipping {f.name}: ProteinMPNN reads .pdb only")
        if path.is_dir():
            for f in sorted(path.glob("*.cif")):
                print(f"Skipping {f.name}: ProteinMPNN reads .pdb only")
    return pdbs


def load_esmfold(chunk_size: int | None = None):
    import torch
    from transformers import EsmForProteinFolding

    print(f"Loading {ESMFOLD_MODEL} (GPU recommended)")
    model = EsmForProteinFolding.from_pretrained(ESMFOLD_MODEL, low_cpu_mem_usage=True)
    model = model.cuda() if torch.cuda.is_available() else model.cpu()
    if chunk_size:  # fold long chains in chunks: slower, but needs far less GPU memory
        model.trunk.set_chunk_size(chunk_size)
    return model.eval()


def process_pdb(args: argparse.Namespace, pdb: Path, out_dir: Path, model) -> list[dict]:
    """Design, fold and score one PDB; returns one row per designed sequence."""
    chains = args.chains if args.skip_fold else (args.chains or first_chain(pdb))
    if not args.skip_fold and len(chains.split()) != 1:
        raise ValueError(f"ESMFold folds one chain at a time; pass a single chain with --chains (got {chains!r})")
    designs = run_proteinmpnn(args, pdb, chains, out_dir)
    if not designs:
        raise ValueError(f"ProteinMPNN wrote no designed sequences to {out_dir}")

    rows = []
    if args.skip_fold:
        return [dict(pdb=pdb.stem, chain=chains or "all", design=n, mpnn_score=sc, sequence=sq) for n, sc, sq in designs]

    import torch

    native_ca, native_seq = read_ca_chain(pdb, chains)
    # ProteinMPNN writes an X for every residue number missing from the structure (unresolved
    # residues, numbering gaps) and for non-standard residues. They are placeholders, not designed
    # positions, and read_ca_chain has no residue for them, so drop them before folding. (The fold
    # then sees the chain without the gap; x_removed records how many residues that was.)
    designs = [(name, score, sequence.replace("X", ""), sequence.count("X")) for name, score, sequence in designs]
    if len(designs[0][2]) != len(native_seq):
        raise ValueError(
            f"Designed length {len(designs[0][2])} (after removing {designs[0][3]} X placeholders) "
            f"!= input chain length {len(native_seq)}"
        )
    folded_dir = out_dir / "folded"
    folded_dir.mkdir(exist_ok=True)
    for name, score, sequence, x_removed in designs:
        with torch.no_grad():
            pdb_text, plddt = fold_esmfold(sequence, model)
        (folded_dir / f"{name}.pdb").write_text(pdb_text)
        tm = tm_score(ca_from_pdb_text(pdb_text), sequence, native_ca, native_seq)
        rows.append(
            dict(pdb=pdb.stem, chain=chains, design=name, mpnn_score=score, plddt=round(plddt, 2), tm_score=round(tm, 4), x_removed=x_removed, sequence=sequence)
        )
        print(f"{pdb.stem} {name}: pLDDT {plddt:.1f}, TM-score {tm:.3f}")
    return rows


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("pdbs", type=Path, nargs="+", help="PDB file(s) and/or folder(s) of .pdb files")
    parser.add_argument("--chains", default="", help="Chain to design, e.g. 'A' (default: each PDB's first chain; with --skip-fold, all chains)")
    parser.add_argument("--num-seqs", type=int, default=4, help="Sequences to design per PDB (default 4)")
    parser.add_argument("--temperature", type=float, default=0.1, help="Sampling temperature (default 0.1)")
    parser.add_argument("--seed", type=int, default=37, help="Random seed (default 37)")
    parser.add_argument("--model-name", default="v_48_020", choices=MODEL_CHOICES, help="ProteinMPNN weights (default v_48_020)")
    parser.add_argument("--use-soluble-model", action="store_true", help="Use weights trained on soluble proteins")
    parser.add_argument("--ca-only", action="store_true", help="CA-only model (input must be CA-only)")
    parser.add_argument("--backbone-noise", type=float, default=0.0, help="Gaussian noise on backbone atoms (default 0)")
    parser.add_argument("--omit-aas", default="X", help="Amino acids to leave out of designs (default X)")
    parser.add_argument("--out-dir", type=Path, default=HERE / "results", help="Output folder (default: results/ next to this script)")
    parser.add_argument("--chunk-size", type=int, default=None, help="ESMFold chunk size (e.g. 64) to fold long chains with less GPU memory; slower (default: off)")
    parser.add_argument("--skip-fold", action="store_true", help="Only run ProteinMPNN; skip ESMFold and TM-score")
    args = parser.parse_args()

    pdbs = collect_pdbs(args.pdbs)
    if not pdbs:
        raise SystemExit("No .pdb files found")
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    model = None if args.skip_fold else load_esmfold(args.chunk_size)

    all_rows, failures = [], []
    for i, pdb in enumerate(pdbs, 1):
        print(f"[{i}/{len(pdbs)}] {pdb.name}")
        try:
            rows = process_pdb(args, pdb, out_dir / f"design_{pdb.stem}", model)
        except Exception as exc:  # keep the batch going; the failure is listed at the end
            traceback.print_exc()
            failures.append((pdb.name, f"{type(exc).__name__}: {exc}"))
            if model is not None:
                import torch

                torch.cuda.empty_cache()
            continue
        write_csv(out_dir / f"design_{pdb.stem}" / "designs.csv", rows)
        all_rows += rows
        write_csv(out_dir / "designs_all.csv", all_rows)  # rewritten after every PDB so an interrupted run keeps its results

    if failures:
        print(f"\n{len(failures)} of {len(pdbs)} PDBs failed:")
        for name, reason in failures:
            print(f"  {name}: {reason}")
        (out_dir / "failures.txt").write_text("\n".join(f"{n}\t{r}" for n, r in failures) + "\n")
    if not all_rows:
        raise SystemExit("No PDB produced designs")
    archive = shutil.make_archive(str(out_dir), "gztar", out_dir.parent, out_dir.name)
    print(f"\nWrote {out_dir / 'designs_all.csv'} ({len(all_rows)} designs from {len(pdbs) - len(failures)} PDBs)")
    print(f"Download this file: {archive}")


if __name__ == "__main__":
    main()
