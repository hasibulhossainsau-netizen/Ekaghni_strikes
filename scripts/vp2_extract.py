#!/usr/bin/env python3
import sys
import subprocess
import yaml
from pathlib import Path
from Bio import SeqIO

def run_blastn(query_fasta, ref_db):
    outfmt = "6 qseqid sseqid pident length qstart qend sstart send evalue sstrand"
    cmd = ["blastn", "-query", str(query_fasta), "-db", str(ref_db), "-outfmt", outfmt, "-max_target_seqs", "1", "-evalue", "1e-5"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"BLASTN error: {result.stderr}")
    lines = result.stdout.strip().split("\n")
    for line in lines:
        if line and not line.startswith("#"):
            return line
    return None

def main():
    PROJECT_ROOT = Path(__file__).parent.parent
    DATA_DIR = PROJECT_ROOT / "data"
    FASTA_DIR = DATA_DIR / "fasta"
    
    with open(PROJECT_ROOT / "config.yaml", "r") as c_file:
        config = yaml.safe_load(c_file)
    
    samples = config["samples"]
    ref_genome_name = Path(config["reference"]["genome"]).stem

    ref_fasta = DATA_DIR / "reference" / f"{ref_genome_name}.fasta"
    if not ref_fasta.exists():
        ref_fasta = FASTA_DIR / f"{ref_genome_name}.fasta"

    ref_db_prefix = FASTA_DIR / f"{ref_genome_name}_db"
    input_files = [FASTA_DIR / f"{s}.fasta" for s in samples if s != ref_genome_name]

    print("[*] Creating BLAST database...")
    make_cmd = ["makeblastdb", "-in", str(ref_fasta), "-dbtype", "nucl", "-out", str(ref_db_prefix)]
    subprocess.run(make_cmd, check=True, stdout=subprocess.DEVNULL)

    for query_file in input_files:
        if not query_file.exists():
            print(f"[-] file wasn't found: {query_file.name}")
            continue

        print(f"[*] processing: {query_file.name}")
        try:
            blast_result = run_blastn(query_file, ref_db_prefix)
            if blast_result is None:
                print(f"[-] {query_file.name}: No hit was found against VP2")
                continue

            parts = blast_result.strip().split("\t")
            qstart, qend = int(parts[4]), int(parts[5])
            strand = parts[9]

            record = SeqIO.read(str(query_file), "fasta")
            start_idx = min(qstart, qend) - 1
            end_idx = max(qstart, qend)
            
            if strand == "plus" or strand == "plus":
                vp2_nt = record.seq[start_idx:end_idx]
            else:
                vp2_nt = record.seq[start_idx:end_idx].reverse_complement()

            vp2_protein = vp2_nt.translate()
            protein_str = str(vp2_protein).rstrip('*')

            out_file = FASTA_DIR / f"{query_file.stem}_VP2.fasta"
            header = f"{query_file.stem}_VP2"
            
            with open(out_file, "w") as f:
                f.write(f">{header}\n{protein_str}\n")
            print(f"[✓] saved (Protein FASTA): {out_file}\n")

        except Exception as e:
            print(f"[x] Error {query_file.name}: {e}")
    for ext in [".nhr", ".nin", ".nsq"]:
        try:
            (ref_db_prefix.with_suffix(ext)).unlink()
        except OSError:
            pass

if __name__ == "__main__":
    main()