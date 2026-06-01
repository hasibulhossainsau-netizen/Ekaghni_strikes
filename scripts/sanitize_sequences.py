#!/usr/bin/env python3
import sys

def sanitize_fasta(input_path, output_path):
    with open(input_path, "r") as infile, open(output_path, "w") as outfile:
        current_header = ""
        current_seq = []
        
        for line in infile:
            line = line.strip()
            if line.startswith(">"):
                if current_header:
                    seq_str = "".join(current_seq).upper().replace(" ", "")
                    if len(seq_str) > 100: 
                        outfile.write(f"{current_header}\n{seq_str}\n")
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)
                
        if current_header:
            seq_str = "".join(current_seq).upper().replace(" ", "")
            if len(seq_str) > 100:
                outfile.write(f"{current_header}\n{seq_str}\n")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python sanitize_sequences.py <input> <output>")
        sys.exit(1)
    sanitize_fasta(sys.argv[1], sys.argv[2])