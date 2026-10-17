#!/usr/bin/env python3
import argparse
import json
import pandas as pd
from pathlib import Path
import sys
from pipeline_config import load_config, require_columns, resolve_path


BLAST_COLUMNS = (
    "qseqid", "sseqid", "pident", "length", "qstart", "qend", "evalue", "qcovs"
)


def filter_blast_hits(blast_df, thresholds):
    """Return unique query regions meeting all configured host-hit thresholds."""
    require_columns(blast_df, BLAST_COLUMNS, "BLAST tabular results")
    high_hits = blast_df[
        (blast_df["pident"] >= thresholds["identity"])
        & (blast_df["qcovs"] >= thresholds["coverage"])
        & (blast_df["evalue"] <= thresholds["evalue"])
    ]
    results = pd.DataFrame({
        "protein_id": high_hits["qseqid"].astype(str).str.strip(),
        "host_match": high_hits["sseqid"].astype(str).str.strip(),
        "identity": high_hits["pident"].astype(float).round(2),
        "query_coverage": high_hits["qcovs"].astype(float).round(2),
        "evalue": high_hits["evalue"].astype(float),
        "qstart": high_hits["qstart"].astype(int),
        "qend": high_hits["qend"].astype(int),
    })
    return results.drop_duplicates(subset=["protein_id", "qstart", "qend"])

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
    for file_path in args.inputs:
        path = Path(file_path)
        if not path.exists() or path.stat().st_size == 0:
            continue
            
        scanned_files_count += 1
        df = pd.read_csv(path, sep="\t", header=None, names=BLAST_COLUMNS)
        if df.empty:
            continue
        blacklisted_positions.extend(
            filter_blast_hits(df, thresholds).to_dict(orient="records")
        )

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