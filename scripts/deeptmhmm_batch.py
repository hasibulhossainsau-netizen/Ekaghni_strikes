#!/usr/bin/env python3
"""
DeepTMHMM Optimized Batch Pipeline for Snakemake Integration.
Reads master protein sequences from FASTA, runs biolib bulk prediction,
parses the results, and projects topology down to individual epitope windows.
"""

import sys
import os
import subprocess
from pathlib import Path
import pandas as pd
import argparse
from pipeline_config import configured_tool, load_config, require_columns, resolve_path

def run_deeptmhmm_bulk(
    fasta_path: Path, output_dir: Path, executable: str, model: str, timeout: int
):
    """Runs DeepTMHMM bulk prediction over the full protein FASTA via BioLib."""
    print(f"\n[PROCESSING] Running DeepTMHMM bulk prediction via BioLib...")
    
    absolute_fasta_path = fasta_path.resolve()
    
    command = [
        executable,
        "run",
        model,
        "--fasta",
        str(absolute_fasta_path),
    ]

    try:
        subprocess.run(
            command,
            cwd=output_dir,
            check=True,
            text=True,
            timeout=timeout,
        )
        print("[SUCCESS] DeepTMHMM bulk run completed successfully!")
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] DeepTMHMM run failed: {e}")
        sys.exit(1)

def parse_tmhmm_output(output_dir: Path, final_csv_path: Path, strong_binders_path: Path = None):
    """
    Parses the 3line topology format, fragments the localization string,
    and maps scores directly to individual target peptides to satisfy downstream rules.
    """
    biolib_file = output_dir / "biolib_results" / "predicted_topologies.3line"
    if not biolib_file.is_file():
        raise FileNotFoundError(
            f"DeepTMHMM output is missing: {biolib_file}. "
            "The prediction command must produce predicted_topologies.3line."
        )

    protein_records = {}
    lines = biolib_file.read_text().splitlines()
    if len(lines) % 3:
        raise ValueError(f"DeepTMHMM output has an incomplete 3-line record: {biolib_file}")
    for index in range(0, len(lines), 3):
        header, sequence, topology = (line.strip() for line in lines[index:index + 3])
        if not header.startswith(">") or not sequence or not topology:
            raise ValueError(
                f"DeepTMHMM output record at line {index + 1} is missing a header, "
                "sequence, or topology string."
            )
        if len(sequence) != len(topology):
            raise ValueError(
                f"DeepTMHMM sequence/topology lengths differ for {header}: "
                f"{len(sequence)} != {len(topology)}."
            )
        protein_records[header[1:].split()[0]] = (sequence.upper(), topology.upper())

    # --- Downstream Sub-fragment Projection Engine ---
    # We read your existing binders list to map full proteins down to explicit sub-peptides
    data = []
    
    # Check if we can intercept the sliding windows from your existing strong binders list
    if strong_binders_path is None or not Path(strong_binders_path).is_file():
        raise FileNotFoundError(
            f"DeepTMHMM binder input is missing: {strong_binders_path}"
        )
    binders_df = pd.read_csv(strong_binders_path)
    require_columns(binders_df, ("peptide",), "DeepTMHMM binder input")
    unique_peptides = binders_df["peptide"].dropna().astype(str).str.upper().str.strip().unique()
    for peptide in unique_peptides:
        matches = [
            (sequence, topology)
            for sequence, topology in protein_records.values()
            if peptide in sequence
        ]
        if not matches:
            raise ValueError(
                f"DeepTMHMM protein sequences do not contain binder peptide {peptide!r}."
            )
        sequence, topology = matches[0]
        start = sequence.index(peptide)
        sub_topology = topology[start:start + len(peptide)]
        data.append({
            "peptide": peptide,
            "topology_score": sub_topology.count("O") / len(peptide),
            "topology_string": sub_topology,
        })

    df = pd.DataFrame(
        data, columns=("peptide", "topology_score", "topology_string")
    )
    df.to_csv(final_csv_path, index=False)
    print(f"[SUCCESS] Projected topology CSV successfully written with explicit 'peptide' column at: {final_csv_path}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_fasta")
    parser.add_argument("output_topology_csv")
    parser.add_argument("strong_binders_csv", nargs="?")
    parser.add_argument("--config", default=None)
    args = parser.parse_args()
    config = load_config(args.config)
    fasta_path = resolve_path(args.input_fasta, config)
    output_csv_path = resolve_path(args.output_topology_csv, config)
    strong_binders_path = (
        resolve_path(args.strong_binders_csv, config)
        if args.strong_binders_csv
        else resolve_path(config["paths"]["mhc1_predictions_csv"], config)
    )
    result_dir = output_csv_path.parent
    result_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("Production-Grade DeepTMHMM Protein Pipeline (Unified Coordinates)")
    print("=" * 60)

    if not fasta_path.exists():
        print(f"[ERROR] Input sequence FASTA file not found at: {fasta_path}")
        sys.exit(1)

    # Execute BioLib calculation layer
    tool_config = config["tools"]["deeptmhmm"]
    run_deeptmhmm_bulk(
        fasta_path,
        result_dir,
        configured_tool(config, "deeptmhmm"),
        tool_config["model"],
        tool_config["timeout_seconds"],
    )
    
    # Process and restructure the mapping arrays
    parse_tmhmm_output(result_dir, output_csv_path, strong_binders_path)

if __name__ == "__main__":
    main()