#!/usr/bin/env python3
import argparse
import json
import pandas as pd
from pathlib import Path
import sys
from pipeline_config import load_config, resolve_path

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="+", required=True)
    parser.add_argument("--output_csv", required=True)
    parser.add_argument("--output_json", required=True)
    parser.add_argument("--config", default=None, help="Path to pipeline config.yaml")
    args = parser.parse_args()
    config = load_config(args.config)
    thresholds = config["host_homology"]["thresholds"]
    args.inputs = [resolve_path(path, config) for path in args.inputs]
    args.output_csv = resolve_path(args.output_csv, config)
    args.output_json = resolve_path(args.output_json, config)

    blacklisted_positions = []
    scanned_files_count = 0
    blast_cols = [
        'qseqid', 'sseqid', 'pident', 'length', 'qstart', 'qend', 'evalue', 'qcovs'
    ]

    for file_path in args.inputs:
        path = Path(file_path)
        if not path.exists() or path.stat().st_size == 0:
            continue
            
        scanned_files_count += 1
        try:
            df = pd.read_csv(path, sep="\t", header=None, names=blast_cols)
            if df.empty:
                continue

            high_hits = df[
                (df['pident'] >= thresholds["identity"]) &
                (df['qcovs'] >= thresholds["coverage"]) &
                (df['evalue'] <= thresholds["evalue"])
            ]
            
            for _, row in high_hits.iterrows():
                blacklisted_positions.append({
                    "protein_id": str(row['qseqid']).strip(),
                    "host_match": str(row['sseqid']).strip(),
                    "identity": round(float(row['pident']), 2),
                    "query_coverage": round(float(row['qcovs']), 2),
                    "evalue": float(row['evalue']),
                    "qstart": int(row['qstart']),
                    "qend": int(row['qend'])
                })
        except Exception as e:
            print(f"[WARNING] Crash in {path.name}: {str(e)}", file=sys.stderr)

    if blacklisted_positions:
        df_out = pd.DataFrame(blacklisted_positions)
        df_out = df_out.drop_duplicates(subset=["protein_id", "qstart", "qend"])
    else:
        df_out = pd.DataFrame([{
            "protein_id": "None", "host_match": "None", "identity": 0.0,
            "query_coverage": 0.0, "evalue": None, "qstart": 0, "qend": 0
        }])
        
    Path(args.output_csv).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output_json).parent.mkdir(parents=True, exist_ok=True)
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