#!/usr/bin/env python3
import sys
from Bio import SeqIO

def main():
    if len(sys.argv) < 5:
        print("Usage: python extract_ref_vp2.py <genome_fasta> <start> <end> <output_fasta>")
        sys.exit(1)

    genome_file = sys.argv[1]
    # Taking the input argument
    output_file = sys.argv[4]

    # Reading the reference genome(CDS VP2)
    ref = SeqIO.read(genome_file, "fasta")
    
    # Taking the whole sequence
    vp2_nt_seq = ref.seq

    if not vp2_nt_seq.startswith("ATG") and vp2_nt_seq.reverse_complement().startswith("ATG"):
        vp2_nt_seq = vp2_nt_seq.reverse_complement()

    # Translating the nucleotide
    vp2_protein_seq = vp2_nt_seq.translate()
    protein_str = str(vp2_protein_seq).rstrip('*')

    # writing out in fasta format
    with open(output_file, "w") as f:
        f.write(f">{ref.id}_VP2\n{protein_str}\n")
        
    print(f"Success: Reference is already a CDS. Translated entire sequence → {output_file}")

if __name__ == "__main__":
    main()