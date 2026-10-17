#!/usr/bin/env python3
import argparse
from Bio import SeqIO
from pipeline_config import load_config, resolve_path

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("genome_fasta", nargs="?")
    parser.add_argument("start", nargs="?", type=int, help="1-based inclusive VP2 start coordinate")
    parser.add_argument("end", nargs="?", type=int, help="1-based inclusive VP2 end coordinate")
    parser.add_argument("output_fasta", nargs="?")
    parser.add_argument("--config", default=None)
    args = parser.parse_args()
    config = load_config(args.config)
    genome_fasta = resolve_path(args.genome_fasta or config["reference"]["genome"], config)
    start = args.start if args.start is not None else config["reference"]["vp2_gene"]["start"]
    end = args.end if args.end is not None else config["reference"]["vp2_gene"]["end"]
    output_fasta = resolve_path(
        args.output_fasta or config["paths"]["reference_vp2_fasta"], config
    )

    ref = SeqIO.read(genome_fasta, "fasta")
    if start < 1 or end < start or end > len(ref.seq):
        parser.error("VP2 coordinates must be a valid 1-based inclusive range in the input genome")

    vp2_nt_seq = ref.seq[start - 1:end]
    # Translating the nucleotide
    vp2_protein_seq = vp2_nt_seq.translate()
    protein_str = str(vp2_protein_seq).rstrip('*')

    # writing out in fasta format
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with open(output_fasta, "w") as f:
        f.write(f">{ref.id}_VP2\n{protein_str}\n")
        
    print(f"Success: Extracted and translated VP2 coordinates → {output_fasta}")

if __name__ == "__main__":
    main()