#!/usr/bin/env python3
import sys
import argparse
import pandas as pd
from pathlib import Path
import subprocess
import tempfile
import os

def peptides_to_fasta(peptides, temp_fasta):
    with open(temp_fasta, 'w') as f:
        for i, pep in enumerate(peptides):
            f.write(f">pep_{i}\n{pep}\n")

def main():
    parser = argparse.ArgumentParser(description="AlgPred2 Allergenicity Prediction Wrapper")
    parser.add_argument("--input", required=True, help="Input CSV containing 'peptide' column")
    parser.add_argument("--output", required=True, help="Output CSV for predictions")
    args = parser.parse_args()

    if not Path(args.input).exists():
        print(f"[ERROR] Input file not found: {args.input}")
        sys.exit(1)

    df = pd.read_csv(args.input)
    df.columns = df.columns.str.lower()
    if 'peptide' not in df.columns:
        print(f"[ERROR] Input CSV must have 'peptide' column. Found: {list(df.columns)}")
        sys.exit(1)

    peptides = df['peptide'].astype(str).str.upper().str.strip().unique().tolist()

    predictions = {}

    with tempfile.TemporaryDirectory() as tmpdir:
        temp_fasta = Path(tmpdir) / "peptides.fa"
        peptides_to_fasta(peptides, temp_fasta)

        original_cwd = os.getcwd()
        os.chdir(tmpdir)

        try:
            cmd = [
                "algpred2",
                "-i", str(temp_fasta),
                "-t", "0.5"
            ]
            result = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True
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

            for i, row in pred_df.iterrows():
                seq_id = f"pep_{i}"
                prediction = row.get("prediction", row.iloc[1] if len(row) > 1 else "Non-Allergen")
                score = row.get("score", row.iloc[2] if len(row) > 2 else 0.0)
                predictions[seq_id] = (str(prediction), float(score))

        except subprocess.CalledProcessError as e:
            print(f"[ERROR] AlgPred2 failed:\nSTDOUT:\n{e.stdout}\nSTDERR:\n{e.stderr}")
            sys.exit(1)
        except FileNotFoundError:
            print("[ERROR] 'algpred2' command not found. Please install with: pip install algpred2")
            sys.exit(1)
        finally:
            os.chdir(original_cwd)

        # Build final output
        result_list = []
        for i, pep in enumerate(peptides):
            seq_id = f"pep_{i}"
            pred, score = predictions.get(seq_id, ("Non-Allergen", 0.0))
            result_list.append({
                "peptide": pep,
                "allergenicity_prediction": pred,
                "allergenicity_score": round(float(score), 4)
            })

    res_df = pd.DataFrame(result_list)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    res_df.to_csv(args.output, index=False)

    print(f"[SUCCESS] AlgPred2 prediction completed → {args.output}")
    print(f"Total peptides processed: {len(peptides)}")

if __name__ == "__main__":
    main()