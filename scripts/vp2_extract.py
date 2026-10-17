#!/usr/bin/env python3
import argparse
import subprocess
from pathlib import Path
from Bio import SeqIO
from pipeline_config import configured_tool, load_config, resolve_path


def run_blastn(
    query_fasta, ref_db, blastn, evalue, max_target_seqs, min_identity, timeout
):
    outfmt = "6 qseqid sseqid pident length qstart qend sstart send evalue sstrand"
    cmd = [
        blastn, "-query", str(query_fasta), "-db", str(ref_db),
        "-outfmt", outfmt, "-max_target_seqs", str(max_target_seqs),
        "-evalue", str(evalue), "-perc_identity", str(min_identity)
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    if result.returncode != 0:
        raise RuntimeError(f"BLASTN error: {result.stderr}")
    lines = result.stdout.strip().split("\n")
    for line in lines:
        if line and not line.startswith("#"):
            return line
    return None

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default=None, help="Path to pipeline config.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    samples = config["samples"]
    ref_genome_name = Path(config["reference"]["genome"]).stem
    fasta_dir = resolve_path(config["paths"]["raw_genomes"])
    ref_fasta = resolve_path(config["reference"]["genome"])
    if not ref_fasta.exists():
        ref_fasta = fasta_dir / f"{ref_genome_name}.fasta"

    ref_db_prefix = resolve_path(config["paths"]["blast"]) / f"{ref_genome_name}_nucleotide_db"
    input_files = [fasta_dir / f"{s}.fasta" for s in samples if s != ref_genome_name]
    blastn = configured_tool(config, "blastn")
    makeblastdb = configured_tool(config, "makeblastdb")
    extraction = config["sequence_extraction"]
    ref_db_prefix.parent.mkdir(parents=True, exist_ok=True)

    print("[*] Creating BLAST database...")
    make_cmd = [makeblastdb, "-in", str(ref_fasta), "-dbtype", "nucl", "-out", str(ref_db_prefix)]
    subprocess.run(
        make_cmd,
        check=True,
        stdout=subprocess.DEVNULL,
        timeout=extraction["timeout_seconds"],
    )

    for query_file in input_files:
        if not query_file.exists():
            print(f"[-] file wasn't found: {query_file.name}")
            continue

        print(f"[*] processing: {query_file.name}")
        try:
            blast_result = run_blastn(
                query_file,
                ref_db_prefix,
                blastn,
                extraction["evalue"],
                extraction["max_target_seqs"],
                extraction["min_identity"],
                extraction["timeout_seconds"],
            )
            if blast_result is None:
                print(f"[-] {query_file.name}: No hit was found against VP2")
                continue

            parts = blast_result.strip().split("\t")
            qstart, qend = int(parts[4]), int(parts[5])
            strand = parts[9]

            record = SeqIO.read(str(query_file), "fasta")
            start_idx = min(qstart, qend) - 1
            end_idx = max(qstart, qend)
            
            if strand == "plus":
                vp2_nt = record.seq[start_idx:end_idx]
            else:
                if not extraction["keep_reverse_complement"]:
                    print(f"[-] {query_file.name}: Reverse-strand hit skipped by configuration")
                    continue
                vp2_nt = record.seq[start_idx:end_idx].reverse_complement()

            vp2_protein = vp2_nt.translate()
            protein_str = str(vp2_protein).rstrip('*')

            out_file = fasta_dir / f"{query_file.stem}_VP2.fasta"
            header = f"{query_file.stem}_VP2"
            
            with open(out_file, "w") as f:
                f.write(f">{header}\n{protein_str}\n")
            print(f"[✓] saved (Protein FASTA): {out_file}\n")

        except Exception as e:
            print(f"[x] Error {query_file.name}: {e}")
    for ext in [".nhr", ".nin", ".nsq"]:
        try:
            Path(f"{ref_db_prefix}{ext}").unlink()
        except OSError:
            pass

if __name__ == "__main__":
    main()