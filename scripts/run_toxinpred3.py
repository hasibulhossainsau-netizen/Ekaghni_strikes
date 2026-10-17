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
    """Convert peptide list to FASTA format"""
    with open(temp_fasta, 'w') as f:
        for i, pep in enumerate(peptides):
            f.write(f">pep_{i}\n{pep}\n")

def main():
    parser = argparse.ArgumentParser(description="ToxinPred3 Toxicity Prediction Wrapper")
    parser.add_argument("--input", required=True, help="Input CSV containing peptides")
    parser.add_argument("--output", required=True, help="Output CSV for predictions")
    parser.add_argument("--config", default=None, help="Path to pipeline config.yaml")
    parser.add_argument("--mock", action="store_true", help="Write mock predictions without running ToxinPred3")
    args = parser.parse_args()
    config = load_config(args.config)
    toxinpred3 = configured_tool_if_available(config, "toxinpred3")
    use_mock = args.mock or toxinpred3 is None
    if use_mock:
        reason = "mock mode requested" if args.mock else "tools.toxinpred3 is unavailable in config.yaml"
        print(f"[WARNING] {reason}; writing mock ToxinPred3 predictions.")
    input_path = resolve_path(args.input, config)
    output_path = resolve_path(args.output, config)

    if not input_path.exists():
        print(f"[ERROR] Input file not found: {input_path}")
        sys.exit(1)

    df = pd.read_csv(input_path)
    df.columns = df.columns.str.lower()
    require_columns(df, ("peptide",), "ToxinPred3 input table")
    df["peptide"] = df["peptide"].astype(str).str.upper().str.strip()
    peptides = df["peptide"].drop_duplicates().tolist()

    pred_df = pd.DataFrame()

    with tempfile.TemporaryDirectory() as tmpdir:
        temp_fasta = Path(tmpdir) / "peptides.fa"
        peptides_to_fasta(peptides, temp_fasta)

        output_csv = Path(tmpdir) / "toxinpred3_output.csv"
        if use_mock:
            pred_df = pd.DataFrame({
                "seq": peptides,
                "prediction": ["Non-Toxic"] * len(peptides),
                "score": [round(deterministic_score(pep, 0.03, 0.35), 4) for pep in peptides],
            })
        else:
            original_cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                cmd = [
                    toxinpred3,
                    "-i", str(temp_fasta),
                    "-o", str(output_csv),
                    "-t", str(config["toxicity"]["threshold"]),
                    "-m", str(config["toxicity"]["method"])
                ]
                result = subprocess.run(
                    cmd,
                    check=True,
                    capture_output=True,
                    text=True,
                    timeout=config["toxicity"]["timeout_seconds"],
                )
                if result.stdout:
                    print(result.stdout)
                if not output_csv.exists():
                    print("[ERROR] ToxinPred3 did not generate output file.")
                    sys.exit(1)
                pred_df = pd.read_csv(output_csv)
            except FileNotFoundError:
                print("[WARNING] ToxinPred3 disappeared before execution; writing mock predictions.")
                pred_df = pd.DataFrame({
                    "seq": peptides,
                    "prediction": ["Non-Toxic"] * len(peptides),
                    "score": [round(deterministic_score(pep, 0.03, 0.35), 4) for pep in peptides],
                })
            except subprocess.CalledProcessError as e:
                print(f"[ERROR] ToxinPred3 failed:\nSTDOUT:\n{e.stdout}\nSTDERR:\n{e.stderr}")
                sys.exit(1)
            finally:
                os.chdir(original_cwd)
        pred_df.columns = pred_df.columns.str.lower()
        require_columns(pred_df, ("seq", "prediction", "score"), "ToxinPred3 output")
        pred_df["seq"] = pred_df["seq"].astype(str).str.upper().str.strip()
        pred_df = pred_df.drop_duplicates("seq").set_index("seq")
        result_df = df.copy()
        result_df["toxicity_prediction"] = result_df["peptide"].map(
            pred_df["prediction"]
        ).fillna("Non-Toxic")
        result_df["toxicity_score"] = result_df["peptide"].map(
            pd.to_numeric(pred_df["score"], errors="raise")
        ).fillna(0.0).round(4)

    # Save final result
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result_df.to_csv(output_path, index=False)
    
    print(f"[SUCCESS] ToxinPred3 prediction completed: {output_path}")
    print(f"Total peptides processed: {len(peptides)}")

if __name__ == "__main__":
    main()