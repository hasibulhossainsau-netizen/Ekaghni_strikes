#!/usr/bin/env python3
import argparse
import time
import subprocess
from pathlib import Path
from pipeline_config import configured_tool, load_config, resolve_path

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_fasta")
    parser.add_argument("db_base_path")
    parser.add_argument("output_txt")
    parser.add_argument("--config", default=None)
    args = parser.parse_args()
    config = load_config(args.config)
    fasta_file = resolve_path(args.input_fasta, config)
    db_base_path = resolve_path(args.db_base_path, config)
    outfile = resolve_path(args.output_txt, config)
    blast_config = config["host_homology"]["blastp"]

    if not fasta_file.exists() or fasta_file.stat().st_size == 0:
        print(f"❌ Error: Input file {fasta_file.name} is empty.")
        sys.exit(1)

    outfile.parent.mkdir(parents=True, exist_ok=True)
    outfmt_str = "6 qseqid sseqid pident length qstart qend evalue qcovs"
    homology_thresholds = config["host_homology"]["thresholds"]

    try:
        start_time = time.time()
        cmd = [
            configured_tool(config, "blastp"),
            "-query", str(fasta_file),
            "-db", str(db_base_path),
            "-out", str(outfile),
            "-max_target_seqs", str(blast_config["max_target_seqs"]),
            "-evalue", str(homology_thresholds["evalue"]),
            "-outfmt", outfmt_str
        ]
        
        print(f"   Running BLASTp for {fasta_file.name} against Host DB...")
        subprocess.run(cmd, check=True, timeout=blast_config["timeout_seconds"])
        print(f"   ✅ Success (Time: {time.time() - start_time:.2f}s) → {outfile.name}")

    except Exception as e:
        print(f"   ❌ BLASTp Runtime Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()