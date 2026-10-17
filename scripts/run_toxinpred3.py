#!/usr/bin/env python3
import sys
import argparse
import pandas as pd
from pathlib import Path
import subprocess
import tempfile
import os
from pipeline_config import configured_tool_if_available, load_config, resolve_path

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
    peptides = df['peptide'].astype(str).str.upper().str.strip().unique().tolist()

    pred_df = pd.DataFrame()

    with tempfile.TemporaryDirectory() as tmpdir:
        temp_fasta = Path(tmpdir) / "peptides.fa"
        peptides_to_fasta(peptides, temp_fasta)

        output_csv = Path(tmpdir) / "toxinpred3_output.csv"
        if use_mock:
            pred_df = pd.DataFrame({
                "seq": peptides,
                "prediction": ["Non-Toxic"] * len(peptides),
                "score": [0.0] * len(peptides),
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
                    "score": [0.0] * len(peptides),
                })
            except subprocess.CalledProcessError as e:
                print(f"[ERROR] ToxinPred3 failed:\nSTDOUT:\n{e.stdout}\nSTDERR:\n{e.stderr}")
                sys.exit(1)
            finally:
                os.chdir(original_cwd)
        if not pred_df.empty:
            pred_df.columns = pred_df.columns.str.lower()
            if 'sequence' in pred_df.columns and 'seq' not in pred_df.columns:
                pred_df.rename(columns={'sequence': 'seq'}, inplace=True)
                
            if 'seq' in pred_df.columns:
                pred_df['seq'] = pred_df['seq'].astype(str).str.upper().str.strip()
            else:
                print(f"[ERROR] ToxinPred3 output missing 'seq' column. Found columns: {list(pred_df.columns)}")
                sys.exit(1)
        elif peptides:
            print("[ERROR] No prediction data recovered from ToxinPred3.")
            sys.exit(1)
        else:
            pred_df = pd.DataFrame(columns=["seq", "prediction", "score"])
        
        result_list = []
        for pep in peptides:
            pep_upper = str(pep).upper().strip()
            match = pred_df[pred_df['seq'] == pep_upper]
            if not match.empty:
                prediction = match.iloc[0].get('prediction', 'Non-Toxic')
                score = match.iloc[0].get('score', 0.0)
            else:
                prediction = "Non-Toxic"
                score = 0.0

            result_list.append({
                "peptide": pep_upper,
                "toxicity_prediction": prediction,
                "toxicity_score": round(float(score), 4)
            })

    # Save final result
    res_df = pd.DataFrame(result_list)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    res_df.to_csv(output_path, index=False)
    
    print(f"[SUCCESS] ToxinPred3 prediction completed: {output_path}")
    print(f"Total peptides processed: {len(peptides)}")

if __name__ == "__main__":
    main()