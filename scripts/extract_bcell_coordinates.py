#!/usr/bin/env python3
import sys
import re
import pandas as pd
from Bio import SeqIO
import argparse
from pipeline_config import load_config, resolve_path

def extract_epitopes(fasta_path, csv_path, output_path):
    print(f"[INFO] Reading residue-level scores from: {csv_path}")
    try:
        df_raw = pd.read_csv(csv_path)
        # Clean columns and clear tracking whitespaces
        df_raw.columns = [c.strip() for c in df_raw.columns]
        if 'Accession' in df_raw.columns:
            df_raw['Accession'] = df_raw['Accession'].astype(str).str.strip()
            
        print("[INFO] Generating per-protein sequence positional index tracks...")
        # Since indices reset per protein, calculate a 1-based residue counter for each group
        df_raw['Residue_idx'] = df_raw.groupby('Accession').cumcount() + 1
        
    except Exception as e:
        print(f"[ERROR] Failed to read raw CSV file: {e}")
        sys.exit(1)

    records = []
    print(f"[INFO] Parsing linear epitope topology from: {fasta_path}")
    
    with open(fasta_path, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            protein_id = record.id.strip()
            seq_str = str(record.seq)
            
            # Identify continuous ranges of uppercase characters (Predicted Epitopes)
            for match in re.finditer(r'[A-Z]+', seq_str):
                peptide = match.group()
                start_pos = match.start() + 1  # 1-based biological coordinate match
                end_pos = match.end()
                length = len(peptide)
                
                # Filter cleanly on both Accession AND our group-specific 1-based position column
                sample_raw = df_raw[
                    (df_raw['Accession'] == protein_id) & 
                    (df_raw['Residue_idx'].between(start_pos, end_pos))
                ]
                
                # Fallback slice mechanism using native offsets if the custom counter yields empty matches
                if sample_raw.empty:
                    protein_block = df_raw[df_raw['Accession'] == protein_id]
                    if not protein_block.empty:
                        sample_raw = protein_block.iloc[match.start():match.end()]

                mean_bp3 = sample_raw['BepiPred-3.0 score'].mean() if not sample_raw.empty else 0.0
                mean_linear = sample_raw['BepiPred-3.0 linear epitope score'].mean() if not sample_raw.empty else 0.0
                
                records.append({
                    "Protein_ID": protein_id,
                    "Peptide": peptide,
                    "Start": start_pos,
                    "End": end_pos,
                    "Length": length,
                    "Mean_BepiPred_Score": round(mean_bp3, 4),
                    "Mean_Linear_Epitope_Score": round(mean_linear, 4)
                })

    out_df = pd.DataFrame(records)
    if out_df.empty:
        out_df = pd.DataFrame(columns=["Protein_ID", "Peptide", "Start", "End", "Length", "Mean_BepiPred_Score", "Mean_Linear_Epitope_Score"])
        
    out_df.to_csv(output_path, index=False)
    print(f"[SUCCESS] Written {len(out_df)} linear B-cell epitope coordinates to: {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_fasta")
    parser.add_argument("raw_output_csv")
    parser.add_argument("output_summary_csv")
    parser.add_argument("--config", default=None)
    args = parser.parse_args()
    config = load_config(args.config)
    extract_epitopes(
        resolve_path(args.input_fasta, config),
        resolve_path(args.raw_output_csv, config),
        resolve_path(args.output_summary_csv, config),
    )