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
        print(f"Error: Input file {input_fasta} not found.")
        sys.exit(1)
        
    sequences = parse_fasta(input_fasta)
    if not sequences:
        print("Error: No sequences found in FASTA.")
        sys.exit(1)
        
    num_seqs = len(sequences)
    seq_len = len(sequences[0])
    
    print(f"[INFO] Calculating FIX-CONSERVATION for {num_seqs} sequences...")
    
    conservation_data = []
    
    for pos in range(seq_len):
        residues_at_pos = [seq[pos] for seq in sequences if pos < len(seq)]
        counts = Counter(residues_at_pos)
        most_common_aa, count = counts.most_common(1)[0]
        if most_common_aa == '-':
            remaining = [r for r in residues_at_pos if r != '-']
            if remaining:
                most_common_aa, count = Counter(remaining).most_common(1)[0]
                conservation_score = (count / num_seqs) * (len(remaining) / num_seqs) 
            else:
                most_common_aa = 'X'
                conservation_score = 0.0
        else:
            conservation_score = count / num_seqs
            
        conservation_data.append({
            "position": pos + 1,
            "consensus_aa": most_common_aa,
            "raw_score": round(conservation_score, 4)
        })
        
    df = pd.DataFrame(conservation_data)
    
    window_size = config["conservation"]["window_size"]
    df[f'window_{window_size}mer_avg'] = df['raw_score'].rolling(
        window=window_size, min_periods=1
    ).mean().round(4)
    
    Path(output_csv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_csv, index=False)
    print(f"[SUCCESS] Ultra-Safe Conservation report generated at: {output_csv}")

if __name__ == "__main__":
    main()