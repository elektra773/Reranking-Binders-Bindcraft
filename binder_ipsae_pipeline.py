#!/usr/bin/env python3
"""Generalized binder iPSAE workflow for arbitrary protein targets.

This script separates the Nipah-specific notebook into reusable steps:
1. Prepare Boltz YAML inputs for one target and many binders.
2. Optionally run Boltz predictions.
3. Score the resulting predictions with the bundled ipsae.py script.
4. Emit flat CSV summaries that are easy to sort and filter.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path


AMINO_ACID_3_TO_1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "MSE": "M",
}

SUMMARY_NUMERIC_COLUMNS = [
    "ipSAE",
    "ipSAE_d0chn",
    "ipSAE_d0dom",
    "ipTM_af",
    "ipTM_d0chn",
    "pDockQ",
    "pDockQ2",
    "LIS",
    "n0res",
    "n0chn",
    "n0dom",
    "d0res",
    "d0chn",
    "d0dom",
    "nres1",
    "nres2",
    "dist1",
    "dist2",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare, run, and score generalized binder iPSAE jobs."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare = subparsers.add_parser(
        "prepare",
        help="Create Boltz YAML inputs for one target and one or more binders.",
    )
    add_target_arguments(prepare)
    add_binder_arguments(prepare)
    add_prepare_output_arguments(prepare)

    predict = subparsers.add_parser(
        "predict-and-score",
        help="Prepare Boltz inputs, run predictions, and score them with iPSAE.",
    )
    add_target_arguments(predict)
    add_binder_arguments(predict)
    add_prepare_output_arguments(predict)
    add_predict_arguments(predict)
    add_score_common_arguments(predict)
    predict.add_argument(
        "--summary-csv",
        type=Path,
        default=Path("ipsae_summary.csv"),
        help="Where to write the flat summary CSV. Default: ipsae_summary.csv",
    )

    score = subparsers.add_parser(
        "score",
        help="Score one existing prediction pair (PAE + structure).",
    )
    score.add_argument("--pae", required=True, type=Path, help="Path to the PAE file.")
    score.add_argument(
        "--structure",
        required=True,
        type=Path,
        help="Path to the structure file (.cif or .pdb).",
    )
    score.add_argument(
        "--job-name",
        default="single_prediction",
        help="Name used in printed summaries. Default: single_prediction",
    )
    add_score_common_arguments(score)
    score.add_argument(
        "--output-json",
        type=Path,
        default=None,
        help="Optional path to save the nested JSON score summary.",
    )
    score.add_argument(
        "--output-csv",
        type=Path,
        default=None,
        help="Optional path to save the flat CSV summary.",
    )

    batch = subparsers.add_parser(
        "score-batch",
        help="Scan an existing Boltz predictions directory and score every model.",
    )
    batch.add_argument(
        "--predictions-root",
        required=True,
        type=Path,
        help="Boltz output directory containing the predictions/ folder.",
    )
    batch.add_argument(
        "--model-index",
        type=int,
        default=0,
        help="Model index to score. Default: 0",
    )
    batch.add_argument(
        "--summary-csv",
        type=Path,
        default=Path("ipsae_summary.csv"),
        help="Where to write the flat summary CSV. Default: ipsae_summary.csv",
    )
    add_score_common_arguments(batch)

    return parser.parse_args()


def add_target_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--target-sequence",
        default="",
        help="Target amino-acid sequence. Use this or --target-fasta or --target-pdb.",
    )
    parser.add_argument(
        "--target-fasta",
        type=Path,
        default=None,
        help="FASTA file containing the target sequence.",
    )
    parser.add_argument(
        "--target-pdb",
        type=Path,
        default=None,
        help="PDB file to extract the target sequence from.",
    )
    parser.add_argument(
        "--target-pdb-chain",
        default="A",
        help="Chain to extract from --target-pdb. Default: A",
    )
    parser.add_argument(
        "--target-chain-id",
        default="A",
        help="Chain ID to use for the target in the Boltz YAML. Default: A",
    )
    parser.add_argument(
        "--target-name",
        default="target",
        help="Human-readable target name stored in the manifest. Default: target",
    )
    parser.add_argument(
        "--target-msa",
        type=Path,
        default=None,
        help="Optional precomputed MSA for the target chain.",
    )


def add_binder_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--binder-sequence",
        default="",
        help="Single binder sequence.",
    )
    parser.add_argument(
        "--binder-name",
        default="binder",
        help="Name for --binder-sequence. Default: binder",
    )
    parser.add_argument(
        "--binder-fasta",
        type=Path,
        default=None,
        help="FASTA file with one or more binder sequences.",
    )
    parser.add_argument(
        "--binder-csv",
        type=Path,
        default=None,
        help="CSV file with one or more binder sequences, such as final_design_stats.csv.",
    )
    parser.add_argument(
        "--binder-name-column",
        default="Design",
        help="Column to use for binder names in --binder-csv. Default: Design",
    )
    parser.add_argument(
        "--binder-sequence-column",
        default="Sequence",
        help="Column to use for sequences in --binder-csv. Default: Sequence",
    )
    parser.add_argument(
        "--binder-rank-column",
        default="Rank",
        help="Optional rank/sort column in --binder-csv. Default: Rank",
    )
    parser.add_argument(
        "--binder-chain-id",
        default="B",
        help="Chain ID to use for the binder in the Boltz YAML. Default: B",
    )
    parser.add_argument(
        "--binder-msa",
        type=Path,
        default=None,
        help="Optional precomputed MSA for every binder chain.",
    )


def add_prepare_output_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--work-dir",
        type=Path,
        default=Path("binder_ipsae_runs"),
        help="Working directory for inputs, manifests, and outputs. Default: binder_ipsae_runs",
    )
    parser.add_argument(
        "--inputs-subdir",
        default="inputs",
        help="Subdirectory under --work-dir for YAML inputs. Default: inputs",
    )


def add_predict_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--boltz-executable",
        default="boltz",
        help="Boltz executable name or path. Default: boltz",
    )
    parser.add_argument(
        "--use-msa-server",
        action="store_true",
        help="Pass --use_msa_server to Boltz.",
    )
    parser.add_argument(
        "--use-potentials",
        action="store_true",
        help="Pass --use_potentials to Boltz.",
    )
    parser.add_argument(
        "--override",
        action="store_true",
        help="Pass --override to Boltz to rerun cached jobs.",
    )
    parser.add_argument(
        "--devices",
        type=int,
        default=1,
        help="Boltz --devices value. Default: 1",
    )
    parser.add_argument(
        "--recycling-steps",
        type=int,
        default=None,
        help="Optional Boltz --recycling_steps value.",
    )
    parser.add_argument(
        "--diffusion-samples",
        type=int,
        default=None,
        help="Optional Boltz --diffusion_samples value.",
    )
    parser.add_argument(
        "--boltz-out-dir",
        type=Path,
        default=None,
        help="Boltz output directory. Default: <work-dir>/boltz_output",
    )


def add_score_common_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--pae-cutoff",
        type=float,
        default=15.0,
        help="PAE cutoff passed into ipsae.py. Default: 15",
    )
    parser.add_argument(
        "--dist-cutoff",
        type=float,
        default=15.0,
        help="Distance cutoff passed into ipsae.py. Default: 15",
    )
    parser.add_argument(
        "--ipsae-script",
        type=Path,
        default=Path(__file__).with_name("ipsae.py"),
        help="Path to the bundled ipsae.py scorer.",
    )


def main() -> int:
    args = parse_args()

    if args.command == "prepare":
        manifest_path = prepare_inputs(args)
        print(f"Wrote manifest to {manifest_path}")
        return 0

    if args.command == "predict-and-score":
        manifest_path = prepare_inputs(args)
        run_boltz_predictions(args)
        records = score_from_manifest(
            manifest_path=manifest_path,
            boltz_out_dir=resolve_boltz_out_dir(args),
            model_index=0,
            pae_cutoff=args.pae_cutoff,
            dist_cutoff=args.dist_cutoff,
            ipsae_script=args.ipsae_script,
        )
        write_flat_summary(args.summary_csv, records)
        print(f"Wrote summary CSV to {args.summary_csv.resolve()}")
        return 0

    if args.command == "score":
        summary = score_prediction(
            pae_path=args.pae,
            structure_path=args.structure,
            pae_cutoff=args.pae_cutoff,
            dist_cutoff=args.dist_cutoff,
            ipsae_script=args.ipsae_script,
            job_name=args.job_name,
        )
        flat_rows = flatten_summary(args.job_name, summary)
        text = json.dumps(summary, indent=2)
        if args.output_json:
            args.output_json.parent.mkdir(parents=True, exist_ok=True)
            args.output_json.write_text(text + "\n", encoding="utf-8")
        if args.output_csv:
            write_flat_summary(args.output_csv, flat_rows)
        print(text)
        return 0

    if args.command == "score-batch":
        records = score_existing_predictions(
            predictions_root=args.predictions_root,
            model_index=args.model_index,
            pae_cutoff=args.pae_cutoff,
            dist_cutoff=args.dist_cutoff,
            ipsae_script=args.ipsae_script,
        )
        write_flat_summary(args.summary_csv, records)
        print(f"Wrote summary CSV to {args.summary_csv.resolve()}")
        return 0

    raise ValueError(f"Unsupported command: {args.command}")


def prepare_inputs(args: argparse.Namespace) -> Path:
    binders = load_binders(args)
    target_sequence = load_target_sequence(args)

    work_dir = args.work_dir.resolve()
    inputs_dir = work_dir / args.inputs_subdir
    inputs_dir.mkdir(parents=True, exist_ok=True)

    manifest_rows = []
    for index, binder in enumerate(binders, start=1):
        binder_name_sanitized = sanitize_name(binder["name"])
        if re.match(r"^\d+_", binder_name_sanitized):
            job_name = binder_name_sanitized
        else:
            job_name = f"{index:03d}_{binder_name_sanitized}"
        yaml_path = inputs_dir / f"{job_name}.yaml"
        yaml_path.write_text(
            render_boltz_yaml(
                target_sequence=target_sequence,
                binder_sequence=binder["sequence"],
                target_chain_id=args.target_chain_id,
                binder_chain_id=args.binder_chain_id,
                target_msa=args.target_msa,
                binder_msa=args.binder_msa,
            ),
            encoding="utf-8",
        )
        manifest_rows.append(
            {
                "job_name": job_name,
                "target_name": args.target_name,
                "target_chain_id": args.target_chain_id,
                "binder_name": binder["name"],
                "binder_chain_id": args.binder_chain_id,
                "binder_length": str(len(binder["sequence"])),
                "yaml_path": str(yaml_path.resolve()),
            }
        )

    manifest_path = work_dir / "binder_manifest.csv"
    write_csv(manifest_path, manifest_rows)
    return manifest_path


def run_boltz_predictions(args: argparse.Namespace) -> None:
    boltz_executable = shutil.which(args.boltz_executable) or args.boltz_executable
    if shutil.which(args.boltz_executable) is None and not Path(args.boltz_executable).exists():
        raise FileNotFoundError(
            "Boltz executable was not found. Install Boltz first or pass --boltz-executable."
        )

    input_path = (args.work_dir / args.inputs_subdir).resolve()
    out_dir = resolve_boltz_out_dir(args)
    out_dir.mkdir(parents=True, exist_ok=True)

    command = [
        boltz_executable,
        "predict",
        str(input_path),
        "--out_dir",
        str(out_dir),
        "--devices",
        str(args.devices),
    ]
    if args.use_msa_server:
        command.append("--use_msa_server")
    if args.use_potentials:
        command.append("--use_potentials")
    if args.override:
        command.append("--override")
    if args.recycling_steps is not None:
        command.extend(["--recycling_steps", str(args.recycling_steps)])
    if args.diffusion_samples is not None:
        command.extend(["--diffusion_samples", str(args.diffusion_samples)])

    print("Running:", " ".join(command))
    result = subprocess.run(command, text=True, capture_output=True)
    if result.stdout:
        print(result.stdout)
    if result.returncode != 0:
        if result.stderr:
            print(result.stderr, file=sys.stderr)
        raise RuntimeError(
            f"Boltz failed with exit code {result.returncode}. "
            "See the printed stdout/stderr above for details."
        )


def resolve_boltz_out_dir(args: argparse.Namespace) -> Path:
    if args.boltz_out_dir is not None:
        return args.boltz_out_dir.resolve()
    return (args.work_dir / "boltz_output").resolve()


def load_target_sequence(args: argparse.Namespace) -> str:
    if args.target_sequence:
        return normalize_sequence(args.target_sequence)

    if args.target_fasta is not None:
        records = read_fasta(args.target_fasta)
        if not records:
            raise ValueError(f"No sequences found in {args.target_fasta}")
        return records[0]["sequence"]

    if args.target_pdb is not None:
        sequence = extract_sequence_from_pdb(
            pdb_path=args.target_pdb,
            chain_id=args.target_pdb_chain,
        )
        if not sequence:
            raise ValueError(
                f"Could not extract a sequence from {args.target_pdb} chain {args.target_pdb_chain}"
            )
        return sequence

    raise ValueError(
        "Provide one of --target-sequence, --target-fasta, or --target-pdb."
    )


def load_binders(args: argparse.Namespace) -> list[dict[str, str]]:
    binders = []
    if args.binder_csv is not None:
        binders.extend(
            read_binder_csv(
                path=args.binder_csv,
                name_column=args.binder_name_column,
                sequence_column=args.binder_sequence_column,
                rank_column=args.binder_rank_column,
            )
        )
    if args.binder_fasta is not None:
        binders.extend(read_fasta(args.binder_fasta))

    if args.binder_sequence:
        binders.append(
            {
                "name": args.binder_name,
                "sequence": normalize_sequence(args.binder_sequence),
            }
        )

    if not binders:
        raise ValueError("Provide --binder-fasta or --binder-sequence.")

    return binders


def read_binder_csv(
    path: Path,
    name_column: str,
    sequence_column: str,
    rank_column: str,
) -> list[dict[str, str]]:
    rows = read_csv_rows(path)
    if not rows:
        raise ValueError(f"No binder rows found in {path}")

    available_columns = set(rows[0])
    missing = [
        column
        for column in (name_column, sequence_column)
        if column not in available_columns
    ]
    if missing:
        raise ValueError(
            f"Missing expected CSV columns {missing} in {path}. "
            f"Available columns start with: {sorted(list(available_columns))[:10]}"
        )

    records = []
    for row in rows:
        sequence = str(row.get(sequence_column, "")).strip()
        name = str(row.get(name_column, "")).strip()
        if not sequence or not name:
            continue

        rank_prefix = ""
        if rank_column and rank_column in row:
            rank_value = str(row.get(rank_column, "")).strip()
            if rank_value:
                rank_prefix = f"{rank_value}_"

        records.append(
            {
                "name": f"{rank_prefix}{name}",
                "sequence": normalize_sequence(sequence),
            }
        )

    if rank_column and records:
        def sort_key(record: dict[str, str]) -> tuple[int, str]:
            prefix, _, remainder = record["name"].partition("_")
            if prefix.isdigit():
                return (int(prefix), remainder)
            return (10**9, record["name"])

        records.sort(key=sort_key)

    return records


def read_fasta(path: Path) -> list[dict[str, str]]:
    records = []
    current_name = None
    current_lines: list[str] = []

    with path.expanduser().resolve().open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    records.append(
                        {
                            "name": current_name,
                            "sequence": normalize_sequence("".join(current_lines)),
                        }
                    )
                current_name = line[1:].strip() or f"binder_{len(records) + 1}"
                current_lines = []
                continue
            current_lines.append(line)

    if current_name is not None:
        records.append(
            {
                "name": current_name,
                "sequence": normalize_sequence("".join(current_lines)),
            }
        )

    return records


def extract_sequence_from_pdb(pdb_path: Path, chain_id: str) -> str:
    residues = []
    seen = set()

    with pdb_path.expanduser().resolve().open(encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            if line[21].strip() != chain_id:
                continue
            atom_name = line[12:16].strip()
            if atom_name != "CA":
                continue
            residue_name = line[17:20].strip().upper()
            residue_number = line[22:27].strip()
            if residue_number in seen:
                continue
            seen.add(residue_number)
            residues.append(AMINO_ACID_3_TO_1.get(residue_name, "X"))

    return "".join(residues)


def render_boltz_yaml(
    target_sequence: str,
    binder_sequence: str,
    target_chain_id: str,
    binder_chain_id: str,
    target_msa: Path | None,
    binder_msa: Path | None,
) -> str:
    def chain_block(chain_id: str, sequence: str, msa_path: Path | None) -> list[str]:
        lines = [
            "  - protein:",
            f"      id: {yaml_quote(chain_id)}",
            f"      sequence: {yaml_quote(sequence)}",
        ]
        if msa_path is not None:
            lines.append(f"      msa: {yaml_quote(str(msa_path.expanduser().resolve()))}")
        else:
            lines.append("      msa: empty")
        return lines

    lines = ["version: 1", "sequences:"]
    lines.extend(chain_block(target_chain_id, target_sequence, target_msa))
    lines.extend(chain_block(binder_chain_id, binder_sequence, binder_msa))
    lines.append("")
    return "\n".join(lines)


def yaml_quote(value: str) -> str:
    return "'" + value.replace("'", "''") + "'"


def normalize_sequence(sequence: str) -> str:
    cleaned = re.sub(r"\s+", "", sequence).upper()
    if not cleaned:
        raise ValueError("Encountered an empty sequence.")
    invalid = sorted({char for char in cleaned if not char.isalpha()})
    if invalid:
        raise ValueError(f"Sequence contains invalid characters: {''.join(invalid)}")
    return cleaned


def sanitize_name(name: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("_")
    return safe or "binder"


def score_from_manifest(
    manifest_path: Path,
    boltz_out_dir: Path,
    model_index: int,
    pae_cutoff: float,
    dist_cutoff: float,
    ipsae_script: Path,
) -> list[dict[str, str]]:
    records = []
    manifest = read_csv_rows(manifest_path)
    for row in manifest:
        prediction = find_prediction_artifacts(
            boltz_out_dir=boltz_out_dir,
            input_stem=Path(row["yaml_path"]).stem,
            model_index=model_index,
        )
        summary = score_prediction(
            pae_path=prediction["pae"],
            structure_path=prediction["structure"],
            pae_cutoff=pae_cutoff,
            dist_cutoff=dist_cutoff,
            ipsae_script=ipsae_script,
            job_name=row["job_name"],
        )
        for flat_row in flatten_summary(row["job_name"], summary):
            flat_row["target_name"] = row["target_name"]
            flat_row["target_chain_id"] = row["target_chain_id"]
            flat_row["binder_name"] = row["binder_name"]
            flat_row["binder_chain_id"] = row["binder_chain_id"]
            records.append(flat_row)
    return records


def score_existing_predictions(
    predictions_root: Path,
    model_index: int,
    pae_cutoff: float,
    dist_cutoff: float,
    ipsae_script: Path,
) -> list[dict[str, str]]:
    predictions_root = predictions_root.expanduser().resolve()
    records = []
    cif_paths = sorted(
        predictions_root.glob(f"predictions/*/*_model_{model_index}.cif")
    )
    for structure_path in cif_paths:
        job_name = structure_path.stem.replace(f"_model_{model_index}", "")
        pae_path = structure_path.with_name(f"pae_{structure_path.stem}.npz")
        if not pae_path.exists():
            raise FileNotFoundError(f"Missing PAE file for {structure_path}: {pae_path}")
        summary = score_prediction(
            pae_path=pae_path,
            structure_path=structure_path,
            pae_cutoff=pae_cutoff,
            dist_cutoff=dist_cutoff,
            ipsae_script=ipsae_script,
            job_name=job_name,
        )
        records.extend(flatten_summary(job_name, summary))
    return records


def find_prediction_artifacts(
    boltz_out_dir: Path,
    input_stem: str,
    model_index: int,
) -> dict[str, Path]:
    boltz_out_dir = boltz_out_dir.expanduser().resolve()
    structure_candidates = sorted(
        boltz_out_dir.glob(f"predictions/*/{input_stem}_model_{model_index}.cif")
    )
    if not structure_candidates:
        structure_candidates = sorted(
            boltz_out_dir.glob(f"predictions/{input_stem}/*_model_{model_index}.cif")
        )
    if not structure_candidates:
        raise FileNotFoundError(
            f"Could not find a structure for input '{input_stem}' under {boltz_out_dir}"
        )

    structure_path = structure_candidates[0]
    pae_path = structure_path.with_name(f"pae_{structure_path.stem}.npz")
    if not pae_path.exists():
        raise FileNotFoundError(f"Expected PAE file next to {structure_path}: {pae_path}")

    return {"structure": structure_path, "pae": pae_path}


def score_prediction(
    pae_path: Path,
    structure_path: Path,
    pae_cutoff: float,
    dist_cutoff: float,
    ipsae_script: Path,
    job_name: str,
) -> dict[str, object]:
    ipsae_script = ipsae_script.expanduser().resolve()
    if not ipsae_script.exists():
        raise FileNotFoundError(f"Could not find ipsae.py at {ipsae_script}")

    command = [
        sys.executable,
        str(ipsae_script),
        str(pae_path.expanduser().resolve()),
        str(structure_path.expanduser().resolve()),
        str(pae_cutoff),
        str(dist_cutoff),
    ]
    subprocess.run(command, check=True)

    table_path = predicted_table_path(structure_path, pae_cutoff, dist_cutoff)
    rows = read_score_table(table_path)
    return summarize_score_rows(job_name=job_name, rows=rows)


def predicted_table_path(structure_path: Path, pae_cutoff: float, dist_cutoff: float) -> Path:
    suffix = structure_path.suffix.lower()
    if suffix not in {".cif", ".pdb"}:
        raise ValueError(f"Unsupported structure suffix for scoring: {structure_path}")
    stem_path = structure_path.with_suffix("")
    return Path(
        f"{stem_path}_{format_cutoff(pae_cutoff)}_{format_cutoff(dist_cutoff)}.txt"
    )


def format_cutoff(value: float) -> str:
    integer = int(value)
    text = str(integer)
    if value < 10:
        return f"0{text}"
    return text


def read_score_table(path: Path) -> list[dict[str, object]]:
    lines = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if stripped:
                lines.append(stripped)

    reader = csv.DictReader(lines)
    rows = []
    for row in reader:
        parsed = dict(row)
        for column in SUMMARY_NUMERIC_COLUMNS:
            parsed[column] = float(parsed[column])
        rows.append(parsed)
    return rows


def summarize_score_rows(job_name: str, rows: list[dict[str, object]]) -> dict[str, object]:
    by_pair: dict[tuple[str, str], list[dict[str, object]]] = {}
    for row in rows:
        pair_key = tuple(sorted((str(row["Chn1"]), str(row["Chn2"]))))
        by_pair.setdefault(pair_key, []).append(row)

    summary: dict[str, object] = {"job_name": job_name, "pairs": {}}
    pairs = summary["pairs"]
    assert isinstance(pairs, dict)

    for pair_key, pair_rows in by_pair.items():
        pair_name = f"{pair_key[0]}-{pair_key[1]}"
        max_rows = [row for row in pair_rows if row["Type"] == "max"]
        asym_rows = [row for row in pair_rows if row["Type"] == "asym"]
        if len(max_rows) != 1 or not asym_rows:
            raise ValueError(f"Unexpected score table contents for pair {pair_name}")

        max_row = max_rows[0]
        min_summary = {
            column: min(float(row[column]) for row in asym_rows)
            for column in SUMMARY_NUMERIC_COLUMNS
        }
        max_summary = {
            column: float(max_row[column]) for column in SUMMARY_NUMERIC_COLUMNS
        }

        pairs[pair_name] = {
            "max": max_summary,
            "min": min_summary,
        }

    return summary


def flatten_summary(job_name: str, summary: dict[str, object]) -> list[dict[str, str]]:
    rows = []
    pairs = summary["pairs"]
    assert isinstance(pairs, dict)
    for pair_name, pair_summary in pairs.items():
        assert isinstance(pair_summary, dict)
        for summary_kind in ("max", "min"):
            metrics = pair_summary[summary_kind]
            assert isinstance(metrics, dict)
            row = {
                "job_name": job_name,
                "pair": pair_name,
                "summary_kind": summary_kind,
            }
            for column in SUMMARY_NUMERIC_COLUMNS:
                row[column] = str(metrics[column])
            rows.append(row)
    return rows


def write_flat_summary(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        raise ValueError("No summary rows were produced.")
    write_csv(path, rows)


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.expanduser().resolve().open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return [dict(row) for row in reader]


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path = path.expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    raise SystemExit(main())
