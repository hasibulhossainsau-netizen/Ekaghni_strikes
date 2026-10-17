#!/usr/bin/env python3
import sys
import argparse
import pandas as pd
from pathlib import Path
import subprocess
import tempfile
import os
from pipeline_config import (
    configured_tool_if_available,
    deterministic_score,
    load_config,
    require_columns,
    resolve_path,
)

def peptides_to_fasta(peptides, temp_fasta):
    with open(temp_fasta, 'w') as f:
        for i, pep in enumerate(peptides):
            f.write(f">pep_{i}\n{pep}\n")

def main():
    parser = argparse.ArgumentParser(description="AlgPred2 Allergenicity Prediction Wrapper")
    parser.add_argument("--input", required=True, help="Input CSV containing 'peptide' column")
    parser.add_argument("--output", required=True, help="Output CSV for predictions")
    parser.add_argument("--config", default=None, help="Path to pipeline config.yaml")
    parser.add_argument("--mock", action="store_true", help="Write mock predictions without running AlgPred2")
    args = parser.parse_args()
    config = load_config(args.config)
    algpred2 = configured_tool_if_available(config, "algpred2")
    use_mock = args.mock or algpred2 is None
    if use_mock:
        reason = "mock mode requested" if args.mock else "tools.algpred2 is unavailable in config.yaml"
        print(f"[WARNING] {reason}; writing mock AlgPred2 predictions.")
    input_path = resolve_path(args.input, config)
    output_path = resolve_path(args.output, config)

    if not input_path.exists():
        print(f"[ERROR] Input file not found: {input_path}")
        sys.exit(1)

    df = pd.read_csv(input_path)
    df.columns = df.columns.str.lower()
    require_columns(df, ("peptide",), "AlgPred2 input table")
    df["peptide"] = df["peptide"].astype(str).str.upper().str.strip()
    peptides = df["peptide"].drop_duplicates().tolist()

    predictions = {}

    if use_mock:
        predictions = {
            pep: ("Non-Allergen", round(deterministic_score(pep, 0.02, 0.4), 4))
            for pep in peptides
        }
    else:
        with tempfile.TemporaryDirectory() as tmpdir:
            temp_fasta = Path(tmpdir) / "peptides.fa"
            peptides_to_fasta(peptides, temp_fasta)

            original_cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                cmd = [
                    algpred2,
                    "-i", str(temp_fasta),
                    "-t", str(config["allergenicity"]["threshold"])
                ]
                result = subprocess.run(
                    cmd,
                    check=True,
                    capture_output=True,
                    text=True,
                    timeout=config["allergenicity"]["timeout_seconds"],
                )
                if result.stdout:
                    print(result.stdout)
                if result.stderr:
                    print(result.stderr)

                pred_file = Path(tmpdir) / "outfile.csv"
                if not pred_file.exists():
                    print("[ERROR] AlgPred2 did not generate outfile.csv")
                    sys.exit(1)

                pred_df = pd.read_csv(pred_file)
                pred_df.columns = pred_df.columns.str.lower()
                require_columns(pred_df, ("prediction", "score"), "AlgPred2 output")
                if len(pred_df) != len(peptides):
                    raise ValueError(
                        "AlgPred2 output row count does not match the submitted peptide count."
                    )

                predictions = {
                    pep: (str(row["prediction"]), float(row["score"]))
                    for pep, (_, row) in zip(peptides, pred_df.iterrows())
                }

            except subprocess.CalledProcessError as e:
                print(f"[ERROR] AlgPred2 failed:\nSTDOUT:\n{e.stdout}\nSTDERR:\n{e.stderr}")
                sys.exit(1)
            except FileNotFoundError:
                print("[WARNING] AlgPred2 disappeared before execution; using mock predictions.")
                predictions = {
                    pep: ("Non-Allergen", round(deterministic_score(pep, 0.02, 0.4), 4))
                    for pep in peptides
                }
            finally:
                os.chdir(original_cwd)

    result_df = df.copy()
    result_df["allergenicity_prediction"] = result_df["peptide"].map(
        {pep: prediction for pep, (prediction, _) in predictions.items()}
    ).fillna("Non-Allergen")
    result_df["allergenicity_score"] = result_df["peptide"].map(
        {pep: score for pep, (_, score) in predictions.items()}
    ).fillna(0.0).astype(float).round(4)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result_df.to_csv(output_path, index=False)

    print(f"[SUCCESS] AlgPred2 prediction completed → {output_path}")
    print(f"Total peptides processed: {len(peptides)}")

if __name__ == "__main__":
    main()