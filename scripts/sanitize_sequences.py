#!/usr/bin/env python3
import argparse
from pipeline_config import load_config, resolve_path

def sanitize_fasta(input_path, output_path, minimum_length):
    with open(input_path, "r") as infile, open(output_path, "w") as outfile:
        current_header = ""
        current_seq = []
        
        for line in infile:
            line = line.strip()
            if line.startswith(">"):
                if current_header:
                    seq_str = "".join(current_seq).upper().replace(" ", "")
                    if len(seq_str) > minimum_length:
                        outfile.write(f"{current_header}\n{seq_str}\n")
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)
                
        if current_header:
            seq_str = "".join(current_seq).upper().replace(" ", "")
            if len(seq_str) > minimum_length:
                outfile.write(f"{current_header}\n{seq_str}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("output")
    parser.add_argument("--config", default=None)
    args = parser.parse_args()
    config = load_config(args.config)
    min_length = config["qc"]["minimum_sequence_length"]
    sanitize_fasta(
        resolve_path(args.input, config),
        resolve_path(args.output, config),
        min_length,
    )