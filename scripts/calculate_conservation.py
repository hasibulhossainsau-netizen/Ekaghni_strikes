#!/usr/bin/env python3
import sys
import pandas as pd
from collections import Counter
from pathlib import Path
import argparse
from pipeline_config import load_config, resolve_path


def parse_fasta(fasta_path):
    sequences = []
    current_seq = []
    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_seq:
                    sequences.append("".join(current_seq))
                    current_seq = []
            else:
                current_seq.append(line)
        if current_seq:
            sequences.append("".join(current_seq))
    return sequences


def calculate_conservation(sequences, window_size):
    """Return per-column consensus and conservation scores for aligned sequences."""
    if not sequences:
        raise ValueError("At least one aligned sequence is required.")
    lengths = {len(sequence) for sequence in sequences}
    if len(lengths) != 1:
        raise ValueError("Aligned FASTA sequences must all have the same length.")
    if window_size < 1:
        raise ValueError("Conservation window size must be a positive integer.")

    conservation_data = []
    for position in range(next(iter(lengths))):
        residues = [sequence[position].upper() for sequence in sequences]
        non_gap_residues = [residue for residue in residues if residue != "-"]
        if not non_gap_residues:
            consensus_aa, conservation_score = "X", 0.0
        else:
            consensus_aa, count = Counter(non_gap_residues).most_common(1)[0]
            conservation_score = count / len(residues)
        conservation_data.append({
            "position": position + 1,
            "consensus_aa": consensus_aa,
            "raw_score": round(conservation_score, 4),
        })

    df = pd.DataFrame(conservation_data)
    df[f"window_{window_size}mer_avg"] = df["raw_score"].rolling(
        window=window_size, min_periods=1
    ).mean().round(4)
    return df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_aligned_fasta")
    parser.add_argument("output_csv")
    parser.add_argument("--config", default=None)
    args = parser.parse_args()
    config = load_config(args.config)
    input_fasta = resolve_path(args.input_aligned_fasta, config)
    output_csv = resolve_path(args.output_csv, config)
    
    if not Path(input_fasta).exists():
        parser.error(f"Input file {input_fasta} not found.")
        
    sequences = parse_fasta(input_fasta)
    if not sequences:
        parser.error(f"No sequences found in FASTA: {input_fasta}")
        
    num_seqs = len(sequences)
    
    print(f"[INFO] Calculating FIX-CONSERVATION for {num_seqs} sequences...")
    window_size = config["conservation"]["window_size"]
    try:
        df = calculate_conservation(sequences, window_size)
    except ValueError as error:
        parser.error(str(error))
    
    Path(output_csv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_csv, index=False)
    print(f"[SUCCESS] Ultra-Safe Conservation report generated at: {output_csv}")

if __name__ == "__main__":
    main()