#!/usr/bin/env python3
import argparse
import json
import pandas as pd
from pathlib import Path
import sys

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="+", required=True)
    parser.add_argument("--output_csv", required=True)
    parser.add_argument("--output_json", required=True)
    args = parser.parse_args()

    blacklisted_positions = []
    scanned_files_count = 0
    blast_cols = ['qseqid', 'sseqid', 'pident', 'length', 'qstart', 'qend', 'evalue']

    for file_path in args.inputs:
        path = Path(file_path)
        if not path.exists() or path.stat().st_size == 0:
            continue
            
        scanned_files_count += 1
        try:
            df = pd.read_csv(path, sep="\t", header=None, names=blast_cols)
            if df.empty:
                continue

            high_hits = df[df['pident'] > 95.0]
            
            for _, row in high_hits.iterrows():
                blacklisted_positions.append({
                    "protein_id": str(row['qseqid']).strip(),
                    "host_match": str(row['sseqid']).strip(),
                    "identity": round(float(row['pident']), 2),
                    "qstart": int(row['qstart']),
                    "qend": int(row['qend'])
                })
        except Exception as e:
            print(f"[WARNING] Crash in {path.name}: {str(e)}", file=sys.stderr)

    if blacklisted_positions:
        df_out = pd.DataFrame(blacklisted_positions)
        df_out = df_out.drop_duplicates(subset=["protein_id", "qstart", "qend"])
    else:
        df_out = pd.DataFrame([{"protein_id": "None", "host_match": "None", "identity": 0.0, "qstart": 0, "qend": 0}])
        
    Path(args.output_csv).parent.mkdir(parents=True, exist_ok=True)
    df_out.to_csv(args.output_csv, index=False)

    summary_report = {
        "scanned_files": scanned_files_count,
        "total_host_hits": len(blacklisted_positions)
    }
    with open(args.output_json, "w") as f:
        json.dump(summary_report, f, indent=4)

    print(f"[SUCCESS] Host filtering finished. Positions blacklisted.")

if __name__ == "__main__":
    main()