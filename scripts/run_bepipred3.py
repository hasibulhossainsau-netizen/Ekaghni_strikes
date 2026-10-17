#!/usr/bin/env python3
"""Run BepiPred 3 or produce schema-compatible mock prediction files."""

import argparse
import csv
import math
import subprocess
import sys
from pathlib import Path

from Bio import SeqIO

from pipeline_config import configured_tool_if_available, load_config, resolve_path


OUTPUT_FASTAS = (
    "Bcell_epitope_preds.fasta",
    "Bcell_epitope_top_20pct_preds.fasta",
    "Bcell_linepitope_top_20pct_preds.fasta",
)
RAW_COLUMNS = (
    "Accession",
    "BepiPred-3.0 score",
    "BepiPred-3.0 linear epitope score",
)


def write_mock_predictions(records, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    raw_path = output_dir / "raw_output.csv"
    with raw_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RAW_COLUMNS)
        writer.writeheader()
        for record in records:
            sequence = str(record.seq).upper()
            positive_count = min(len(sequence), max(1, math.ceil(len(sequence) * 0.2)))
            for index in range(len(sequence)):
                score = 0.9 if index < positive_count else 0.1
                writer.writerow({
                    "Accession": record.id,
                    "BepiPred-3.0 score": score,
                    "BepiPred-3.0 linear epitope score": score,
                })

    all_predictions = []
    selected_predictions = []
    for record in records:
        sequence = str(record.seq).upper()
        positive_count = min(len(sequence), max(1, math.ceil(len(sequence) * 0.2)))
        selected = sequence[:positive_count].upper() + sequence[positive_count:].lower()
        all_predictions.append(f">{record.id}\n{sequence}\n")
        selected_predictions.append(f">{record.id}\n{selected}\n")

    for filename in OUTPUT_FASTAS[:1]:
        (output_dir / filename).write_text("".join(all_predictions))
    for filename in OUTPUT_FASTAS[1:]:
        (output_dir / filename).write_text("".join(selected_predictions))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Protein FASTA input")
    parser.add_argument("--output-dir", required=True, help="BepiPred output directory")
    parser.add_argument("--prediction-mode", required=True)
    parser.add_argument("--threshold", type=float, required=True)
    parser.add_argument("--window-size", type=int, required=True)
    parser.add_argument("--config", default=None)
    parser.add_argument("--mock", action="store_true", help="Generate deterministic mock output files")
    args = parser.parse_args()

    config = load_config(args.config)
    input_path = resolve_path(args.input, config)
    output_dir = resolve_path(args.output_dir, config)
    if not input_path.is_file():
        parser.error(f"Input FASTA does not exist: {input_path}")

    records = list(SeqIO.parse(input_path, "fasta"))
    executable = configured_tool_if_available(config, "bepipred3")
    use_mock = args.mock or executable is None
    if use_mock:
        reason = "mock mode requested" if args.mock else "tools.bepipred3 is unavailable in config.yaml"
        print(f"[WARNING] {reason}; writing schema-compatible mock BepiPred 3 outputs.")
        write_mock_predictions(records, output_dir)
        return

    output_dir.mkdir(parents=True, exist_ok=True)
    command = [sys.executable, executable] if executable.endswith(".py") else [executable]
    command.extend([
        "-i", str(input_path),
        "-o", str(output_dir),
        "-pred", args.prediction_mode,
        "-t", str(args.threshold),
        "-rolling_window_size", str(args.window_size),
    ])
    try:
        subprocess.run(command, check=True)
    except FileNotFoundError:
        print("[WARNING] BepiPred 3 disappeared before execution; writing mock outputs.")
        write_mock_predictions(records, output_dir)


if __name__ == "__main__":
    main()
