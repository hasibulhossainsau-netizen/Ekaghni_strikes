#!/usr/bin/env python3
import sys
import time
import subprocess
from pathlib import Path
from Bio import SeqIO

def main():
    if len(sys.argv) < 4:
        print("Usage: python vp2_blastp.py <input_fasta> <db_base_path> <output_txt>")
        sys.exit(1)

    fasta_file = Path(sys.argv[1])
    db_base_path = sys.argv[2]
    outfile = Path(sys.argv[3])

    if not fasta_file.exists() or fasta_file.stat().st_size == 0:
        print(f"❌ Error: Input file {fasta_file.name} is empty.")
        sys.exit(1)

    outfile.parent.mkdir(parents=True, exist_ok=True)
    outfmt_str = "6 qseqid sseqid pident length qstart qend evalue"

    try:
        start_time = time.time()
        cmd = [
            "blastp",
            "-query", str(fasta_file),
            "-db", str(db_base_path),
            "-out", str(outfile),
            "-max_target_seqs", "10",
            "-evalue", "1e-5",
            "-outfmt", outfmt_str
        ]
        
        print(f"   Running BLASTp for {fasta_file.name} against Host DB...")
        subprocess.run(cmd, check=True, timeout=120)
        print(f"   ✅ Success (Time: {time.time() - start_time:.2f}s) → {outfile.name}")

    except Exception as e:
        print(f"   ❌ BLASTp Runtime Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()