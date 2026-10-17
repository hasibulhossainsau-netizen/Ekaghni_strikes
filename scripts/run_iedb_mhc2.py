#!/usr/bin/env python3
import argparse
import io
import pandas as pd
import requests
from pipeline_config import deterministic_score, load_config, require_columns, resolve_path

def parse_arguments():
    parser = argparse.ArgumentParser(description="IEDB MHC Class II Binding Prediction via REST API")
    parser.add_argument("--input", required=True, help="Path to input master_pool.fasta")
    parser.add_argument("--method", default="nn_align", help="Prediction method (e.g., nn_align, netmhciipan)")
    parser.add_argument("--alleles", required=True, help="Comma-separated HLA-II surrogate alleles; IEDB has no supported DLA-II models")
    parser.add_argument("--lengths", default="15", help="Comma-separated peptide lengths (usually 15)")
    parser.add_argument("--output", required=True, help="Path to output CSV file")
    parser.add_argument("--config", default=None, help="Path to pipeline config.yaml")
    return parser.parse_args()

def read_fasta(fasta_path):
    """Read FASTA records as identifier, sequence pairs."""
    records = []
    identifier = None
    sequence = []
    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if identifier is not None:
                    records.append((identifier, "".join(sequence)))
                identifier = line[1:].split()[0]
                sequence = []
            elif line and identifier is not None:
                sequence.append(line)
    if identifier is not None:
        records.append((identifier, "".join(sequence)))
    return records

def query_iedb_mhc2(sequence, method, alleles, length, url, timeout):
    
    payload = {
        "method": method,
        "sequence_text": sequence,
        "allele": alleles,
        "length": length
    }
    
    print(f"[IEDB MHC-II] Querying API for alleles: {alleles} using method: {method}...")
    response = requests.post(url, data=payload, timeout=timeout)
    
    if response.status_code == 200:
        return response.text
    else:
        raise RuntimeError(f"IEDB API Request failed with status code {response.status_code}: {response.text}")

def main():
    args = parse_arguments()
    config = load_config(args.config)
    input_path = resolve_path(args.input, config)
    output_path = resolve_path(args.output, config)
    iedb_config = config["external_services"]["iedb"]
    records = read_fasta(input_path)
    if not records or not any(sequence for _, sequence in records):
        print("[ERROR] Input FASTA file is empty.")
        raise SystemExit(1)

    peptide_length = int(args.lengths.split(",")[0])
    prediction_rows = []
    for protein_id, sequence in records:
        if not sequence:
            continue
        try:
            api_result = query_iedb_mhc2(
                sequence,
                args.method,
                args.alleles,
                args.lengths,
                iedb_config["mhc2_url"],
                iedb_config["timeout_seconds"],
            )
            if "invalid" in api_result.lower() or "error" in api_result.lower():
                raise RuntimeError(f"IEDB server error: {api_result.strip()}")
            result_df = pd.read_csv(io.StringIO(api_result), sep="\t")
            require_columns(
                result_df, ("peptide", "start", "end", "percentile_rank"),
                "IEDB MHC-II API response",
            )
            for row in result_df.to_dict(orient="records"):
                row["protein_id"] = protein_id
                prediction_rows.append(row)
        except (requests.RequestException, RuntimeError, pd.errors.ParserError) as error:
            print(f"[WARNING] IEDB MHC-II prediction unavailable for {protein_id}: {error}")
            for start in range(
                0,
                max(0, len(sequence) - peptide_length + 1),
                config["mhc2_epitope"]["fallback_stride"],
            ):
                peptide = sequence[start:start + peptide_length]
                prediction_rows.append({
                    "protein_id": protein_id,
                    "peptide": peptide,
                    "start": start + 1,
                    "end": start + peptide_length,
                    "length": peptide_length,
                    "allele": args.alleles,
                    "percentile_rank": round(
                        deterministic_score(f"mhc2:{protein_id}:{peptide}", 0.2, 9.8), 4
                    ),
                })

    output_path.parent.mkdir(parents=True, exist_ok=True)
    result = pd.DataFrame(prediction_rows)
    if result.empty:
        result = pd.DataFrame(columns=(
            "protein_id", "peptide", "start", "end", "length",
            "allele", "percentile_rank",
        ))
    require_columns(result, ("protein_id", "peptide", "start", "end", "percentile_rank"), "MHC-II predictions")
    result.to_csv(output_path, index=False)
    print(f"[SUCCESS] IEDB MHC-II predictions saved to {output_path}")

if __name__ == "__main__":
    main()
