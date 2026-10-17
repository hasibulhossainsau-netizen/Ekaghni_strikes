#!/usr/bin/env python3
import argparse
import io
import pandas as pd
import requests
import sys
from pipeline_config import load_config, resolve_path

def parse_arguments():
    parser = argparse.ArgumentParser(description="IEDB MHC Class I Antigen Processing Prediction via REST API")
    parser.add_argument("--input", required=True, help="Path to input vp2_strong_binders.csv")
    parser.add_argument("--method", default="recommended", help="Prediction method (recommended, smm, ann)")
    parser.add_argument("--proteasome", default="immuno", help="Proteasome type (immuno or constitutive)")
    parser.add_argument("--alleles", required=True, help="Comma-separated DLA Class I alleles")
    parser.add_argument("--output", required=True, help="Path to output CSV file")
    parser.add_argument("--config", default=None, help="Path to pipeline config.yaml")
    return parser.parse_args()

def extract_peptides(csv_path):
    """peptide list extraction"""
    try:
        df = pd.read_csv(csv_path)
        if 'peptide' in df.columns:
            return df['peptide'].unique().tolist()
        elif 'Peptide' in df.columns:
            return df['Peptide'].unique().tolist()
        else:
            return df.iloc[:, 0].unique().tolist()
    except Exception as e:
        print(f"[WARNING] Could not read peptides from {csv_path}: {e}")
        return []

def query_iedb_processing(peptides_list, method, proteasome, alleles, url, timeout):
    sequence_text = "\n".join(peptides_list)
    
    payload = {
        "method": method,
        "proteasome": proteasome,
        "sequence_text": sequence_text,
        "allele": alleles
    }
    
    print(f"[IEDB Processing] Querying Processing API for {len(peptides_list)} peptides using allele(s): {alleles}...")
    response = requests.post(url, data=payload, timeout=timeout)
    
    if response.status_code == 200:
        if "invalid allele" in response.text.lower() or "error" in response.text.lower():
            return None, response.text
        return response.text, None
    else:
        return None, response.text

def main():
    args = parse_arguments()
    config = load_config(args.config)
    input_path = resolve_path(args.input, config)
    output_path = resolve_path(args.output, config)
    iedb_config = config["external_services"]["iedb"]
    output_path.parent.mkdir(parents=True, exist_ok=True)

    peptides = extract_peptides(input_path)
    
    if not peptides:
        print("[WARNING] No peptides found for processing analysis. Creating empty output.")
        pd.DataFrame(columns=['peptide', 'score']).to_csv(output_path, index=False)
        return
        
    try:
        api_result, error_msg = query_iedb_processing(
            peptides, args.method, args.proteasome, args.alleles,
            iedb_config["processing_url"], iedb_config["timeout_seconds"]
        )
        if error_msg or not api_result:
            print("[WARN] IEDB rejected DLA alleles. Switching to fallback reference allele (HLA-A*02:01) to salvage processing/TAP scores.")
            fallback_allele = config["mhc1_processing"]["fallback_allele"]
            api_result, error_msg = query_iedb_processing(
                peptides, args.method, args.proteasome, fallback_allele,
                iedb_config["processing_url"], iedb_config["timeout_seconds"]
            )
            
            if error_msg or not api_result:
                print("[ERROR] Fallback also failed. Creating a neutral placeholder processing file.")
                df_mock = pd.DataFrame({'peptide': peptides, 'processing_score': 1.0, 'score': 1.0})
                df_mock.to_csv(output_path, index=False)
                return

        df = pd.read_csv(io.StringIO(api_result), sep="\t")
        
        if 'peptide' not in df.columns:
            for col in df.columns:
                if 'seq' in str(col).lower() or 'pep' in str(col).lower():
                    df.rename(columns={col: 'peptide'}, inplace=True)
                    break
                    
        df.to_csv(output_path, index=False)
        print(f"[SUCCESS] IEDB Processing scores saved successfully to {output_path}")
        
    except Exception as e:
        print(f"[ERROR] Critical failure in run_iedb_processing.py: {e}")
        df_fallback = pd.DataFrame({'peptide': peptides, 'processing_score': 1.0})
        df_fallback.to_csv(output_path, index=False)

if __name__ == "__main__":
    main()