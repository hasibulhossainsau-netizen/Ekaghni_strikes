#!/usr/bin/env python3
"""
Robust T-cell epitope prediction using NetMHCpan 4.1 for Snakemake.
Dynamically accepts alleles, lengths, and thresholds from Snakemake params.
"""
import csv
import os
import sys
import subprocess
import tempfile
from pathlib import Path
from Bio import SeqIO

NETMHCPAN_PATH = Path.home() / "netMHCpan-4.1" / "netMHCpan"

def run_netmhcpan(allele, sequence, length):
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp:
        tmp.write(f">seq\n{sequence}\n")
        tmp_path = tmp.name

    cmd = [
        str(NETMHCPAN_PATH),
        "-f", tmp_path,
        "-a", allele,
        "-l", str(length),
        "-BA",
        "-t", "-0.5"  
    ]

    try:
        print(f"   → {allele} len={length}...", end=" ", flush=True)
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600, check=True)
        
        predictions = []
        header_found = False
        header = []

        for line in result.stdout.splitlines():
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            if "Pos" in line and "Peptide" in line and "Core" in line:
                header = line.split()
                header_found = True
                continue

            if header_found and line:
                values = line.split()
                if len(values) >= 8 and values[0].isdigit():
                    pred = {}
                    for i, col in enumerate(header):
                        if i < len(values):
                            key = col.lower().replace('%', '').replace('(', '').replace(')', '').replace('.', '')
                            if key:
                                pred[key] = values[i]
                    predictions.append(pred)

        print(f"Done ({len(predictions)} peptides)")
        return predictions

    except subprocess.CalledProcessError as e:
        print(f"CLI Error: {e.stderr}")
        return []
    finally:
        if os.path.exists(tmp_path):
            Path(tmp_path).unlink(missing_ok=True)

def main():
    if not NETMHCPAN_PATH.exists():
        print(f"❌ NetMHCpan binary missing at: {NETMHCPAN_PATH}")
        sys.exit(1)

    if len(sys.argv) < 6:
        print("Usage: python netmhcpan_runner.py <input_fasta> <alleles_str> <lengths_str> <rank_threshold> <output_csv>")
        sys.exit(1)

    input_fasta = Path(sys.argv[1])
    dla_alleles = sys.argv[2].split(",")          
    peptide_lengths = [int(l) for l in sys.argv[3].split(",")] 
    rank_threshold = float(sys.argv[4])           
    output_csv = Path(sys.argv[5])

    output_csv.parent.mkdir(parents=True, exist_ok=True)

    if not input_fasta.exists():
        print(f"❌ Input file missing: {input_fasta}")
        sys.exit(1)

    fieldnames = ["protein_id", "allele", "peptide", "start", "end", "length", "ic50", "percentile_rank"]
    strong_binders_pool = []


    for record in SeqIO.parse(input_fasta, "fasta"):
        protein_id = record.id.replace("_VP2_protein", "").replace("_VP2", "")
        seq = str(record.seq).strip()

        print(f"\n📊 Processing {protein_id} ({len(seq)} aa)")

        all_preds = []
        for allele in dla_alleles:
            for length in peptide_lengths:
                preds = run_netmhcpan(allele, seq, length)
                for p in preds:
                    p["target_allele"] = allele
                    p["target_length"] = length
                all_preds.extend(preds)


        for p in all_preds:
            try:
                rank_str = p.get("rank_el", p.get("rank_ba", p.get("rank", "100")))
                rank = float(rank_str)
                
                if rank < rank_threshold:
                    start = int(p.get("pos", 0))
                    length_val = int(p.get("target_length", 9))
                    
                    strong_binders_pool.append({
                        "protein_id": protein_id,
                        "allele": p.get("target_allele"),
                        "peptide": p.get("peptide"),
                        "start": start,
                        "end": start + length_val - 1,
                        "length": length_val,
                        "ic50": p.get("affinitynm", p.get("ic50", "N/A")),
                        "percentile_rank": rank
                    })
            except (ValueError, TypeError):
                continue

    print(f"\n🎯 Total strong binders found: {len(strong_binders_pool)}")


    if strong_binders_pool:
        df_out = sorted(strong_binders_pool, key=lambda x: x["percentile_rank"])
    else:

        df_out = [{"protein_id": "None", "allele": "None", "peptide": "None", "start": 0, "end": 0, "length": 0, "ic50": "N/A", "percentile_rank": 99.0}]

    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(df_out)

    print(f"✅ [SUCCESS] Results safely saved to: {output_csv}")

if __name__ == "__main__":
    main()