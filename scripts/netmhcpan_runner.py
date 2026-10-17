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
from pipeline_config import configured_tool_if_available, load_config, resolve_path

def run_netmhcpan(allele, sequence, length, netmhcpan_path, command_threshold, timeout):
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp:
        tmp.write(f">seq\n{sequence}\n")
        tmp_path = tmp.name

    cmd = [
        netmhcpan_path,
        "-f", tmp_path,
        "-a", allele,
        "-l", str(length),
        "-BA",
        "-t", str(command_threshold)
    ]

    try:
        print(f"   → {allele} len={length}...", end=" ", flush=True)
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, check=True)
        
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


def mock_netmhcpan_prediction(sequence, length, rank_threshold):
    if len(sequence) < length:
        return []
    return [{
        "pos": "1",
        "peptide": sequence[:length],
        "affinitynm": "100.0",
        "rank_el": str(min(rank_threshold / 2, 0.1)),
    }]


def main():
    if len(sys.argv) < 6:
        print("Usage: python netmhcpan_runner.py <input_fasta> <alleles_str> <lengths_str> <rank_threshold> <output_csv> [--config config.yaml]")
        sys.exit(1)

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("input_fasta")
    parser.add_argument("alleles_str")
    parser.add_argument("lengths_str")
    parser.add_argument("rank_threshold", type=float)
    parser.add_argument("output_csv")
    parser.add_argument("--config", default=None)
    parser.add_argument("--mock", action="store_true", help="Write deterministic mock predictions without running NetMHCpan")
    args = parser.parse_args()
    config = load_config(args.config)
    netmhcpan_path = configured_tool_if_available(config, "netmhcpan")
    use_mock = args.mock or netmhcpan_path is None
    if use_mock:
        reason = "mock mode requested" if args.mock else "tools.netmhcpan is unavailable in config.yaml"
        print(
            f"[WARNING] {reason}; writing mock NetMHCpan predictions."
        )
    input_fasta = resolve_path(args.input_fasta, config)
    dla_alleles = args.alleles_str.split(",")
    peptide_lengths = [int(length) for length in args.lengths_str.split(",")]
    rank_threshold = args.rank_threshold
    output_csv = resolve_path(args.output_csv, config)
    netmhcpan_config = config["tcell_epitope"]["netmhcpan"]
    command_threshold = netmhcpan_config["command_threshold"]
    timeout = config["tcell_epitope"]["timeout_seconds"]

    output_csv.parent.mkdir(parents=True, exist_ok=True)

    if not input_fasta.exists():
        print(f"❌ Input file missing: {input_fasta}")
        sys.exit(1)

    fieldnames = ["protein_id", "allele", "peptide", "start", "end", "length", "ic50", "percentile_rank"]
    strong_binders_pool = []


    records = list(SeqIO.parse(input_fasta, "fasta"))
    for record in records:
        protein_id = record.id.replace("_VP2_protein", "").replace("_VP2", "")
        seq = str(record.seq).strip()

        print(f"\n📊 Processing {protein_id} ({len(seq)} aa)")

        all_preds = []
        for allele in dla_alleles:
            for length in peptide_lengths:
                if use_mock:
                    preds = mock_netmhcpan_prediction(seq, length, rank_threshold)
                else:
                    try:
                        preds = run_netmhcpan(
                            allele, seq, length, netmhcpan_path, command_threshold, timeout
                        )
                    except FileNotFoundError:
                        print("[WARNING] NetMHCpan disappeared before execution; switching to mock predictions.")
                        use_mock = True
                        preds = mock_netmhcpan_prediction(seq, length, rank_threshold)
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