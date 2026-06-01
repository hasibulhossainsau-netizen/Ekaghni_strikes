#!/usr/bin/env python3
import sys
import argparse
import pandas as pd
from pathlib import Path
import subprocess
import tempfile
import os

def peptides_to_fasta(peptides, temp_fasta):
    """Convert peptide list to FASTA format"""
    with open(temp_fasta, 'w') as f:
        for i, pep in enumerate(peptides):
            f.write(f">pep_{i}\n{pep}\n")

def main():
    parser = argparse.ArgumentParser(description="ToxinPred3 Toxicity Prediction Wrapper")
    parser.add_argument("--input", required=True, help="Input CSV containing peptides")
    parser.add_argument("--output", required=True, help="Output CSV for predictions")
    args = parser.parse_args()

    if not Path(args.input).exists():
        print(f"[ERROR] Input file not found: {args.input}")
        sys.exit(1)

    df = pd.read_csv(args.input)
    df.columns = df.columns.str.lower()
    peptides = df['peptide'].astype(str).str.upper().str.strip().unique().tolist()

    pred_df = pd.DataFrame()

    with tempfile.TemporaryDirectory() as tmpdir:
        temp_fasta = Path(tmpdir) / "peptides.fa"
        peptides_to_fasta(peptides, temp_fasta)

        output_csv = Path(tmpdir) / "toxinpred3_output.csv"
        original_cwd = os.getcwd()
        os.chdir(tmpdir)
        
        try:
            cmd = [
                "toxinpred3",
                "-i", str(temp_fasta),
                "-o", str(output_csv),
                "-t", "0.5",
                "-m", "1"
            ]
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            if result.stdout:
                print(result.stdout)

            if not output_csv.exists():
                print("[ERROR] ToxinPred3 did not generate output file.")
                sys.exit(1)

            pred_df = pd.read_csv(output_csv)

        except subprocess.CalledProcessError as e:
            print(f"[ERROR] ToxinPred3 failed:\nSTDOUT:\n{e.stdout}\nSTDERR:\n{e.stderr}")
            sys.exit(1)
        except FileNotFoundError:
            print("[ERROR] 'toxinpred3' command not found. Please install with: pip install toxinpred3")
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
        else:
            print("[ERROR] No prediction data recovered from ToxinPred3.")
            sys.exit(1)
        
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
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    res_df.to_csv(args.output, index=False)
    
    print(f"[SUCCESS] ToxinPred3 prediction completed: {args.output}")
    print(f"Total peptides processed: {len(peptides)}")

if __name__ == "__main__":
    main()