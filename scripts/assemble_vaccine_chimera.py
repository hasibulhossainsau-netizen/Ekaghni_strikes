#!/usr/bin/env python3
import argparse
import pandas as pd
import os
import sys
from pipeline_config import load_config, resolve_path

def parse_args():
    parser = argparse.ArgumentParser(
        description="Multi-Epitope Vaccine Chimera Structural Assembly Engine - Production Grade"
    )
    parser.add_argument("--intersected", required=True, help="Overlapping consensus regions CSV")
    parser.add_argument("--mhc1", required=True, help="Filtered MHC I epitopes CSV")
    parser.add_argument("--mhc2", required=True, help="Filtered MHC II epitopes CSV")
    parser.add_argument("--bcell", required=True, help="Filtered B-cell linear epitopes CSV")
    parser.add_argument("--output", required=True, help="Target path for final structural FASTA export")
    parser.add_argument("--config", default=None, help="Path to pipeline config.yaml")
    return parser.parse_args()

def get_peptide_list(df, max_count):
    """
    Extracts valid amino acid sequence strings based on exact case-insensitive matches.
    Defensively strips out structural IDs, metadata rows, and header artifacts.
    """
    if df.empty:
        return []
    
    # Map headers to lowercase to handle casing variants (e.g., 'Peptide' vs 'peptide')
    col_mapping = {col.lower(): col for col in df.columns}
    chosen_col = None

    # Explicitly check for preferred sequence tracking targets
    for target in ["peptide", "sequence", "seq", "epitope"]:
        if target in col_mapping:
            chosen_col = col_mapping[target]
            break

    # Robust fallback: find the first column that resembles valid uppercase peptide letters
    if not chosen_col:
        for col in df.columns:
            sample = df[col].dropna().astype(str).str.strip().iloc[:3]
            # Skip columns that contain metadata markers like underscores, dots, or numerals
            if any("_" in s or "." in s or any(char.isdigit() for char in s) for s in sample):
                continue
            chosen_col = col
            break

    if not chosen_col:
        print("[WARN] Could not strictly isolate a pure sequence column. Defaulting to first column.")
        chosen_col = df.columns[0]

    # Extract clean, uppercase, unique sequence strings
    raw_list = df[chosen_col].dropna().astype(str).str.upper().str.strip().unique()
    clean_peptides = []
    
    for pep in raw_list:
        # Strict biological string validation: alphabet only, no IDs (like FJ197845.1)
        if pep.isalpha() and "_" not in pep and "." not in pep and len(pep) >= 4:
            clean_peptides.append(pep)
            
    return clean_peptides[:max_count]

def main():
    args = parse_args()
    config = load_config(args.config)
    args.intersected = resolve_path(args.intersected, config)
    args.mhc1 = resolve_path(args.mhc1, config)
    args.mhc2 = resolve_path(args.mhc2, config)
    args.bcell = resolve_path(args.bcell, config)
    args.output = resolve_path(args.output, config)
    construct = config["vaccine_construct"]
    max_epitopes = construct["max_epitopes_per_class"]

    print("[PROCESSING] Loading filtered consensus epitope datasets...")
    mhc1_df = pd.read_csv(args.mhc1)
    mhc2_df = pd.read_csv(args.mhc2)
    bcell_df = pd.read_csv(args.bcell)

    # Extract sequence collections cleanly
    mhc1_peps = get_peptide_list(mhc1_df, max_count=max_epitopes)
    mhc2_peps = get_peptide_list(mhc2_df, max_count=max_epitopes)
    bcell_peps = get_peptide_list(bcell_df, max_count=max_epitopes)

    print(f"[INFO] Parsed MHC-I Epitopes:  {mhc1_peps}")
    print(f"[INFO] Parsed MHC-II Epitopes: {mhc2_peps}")
    print(f"[INFO] Parsed B-Cell Epitopes: {bcell_peps}")

    # 1. Adjuvant: Flagellin FliC from Salmonella enterica (Full immunogenic segment)
    adjuvant = construct["adjuvant"]
    
    # 2. Universal T-helper Element: PADRE Sequence 
    padre = construct["helper_epitope"]

    # 3. Spacers & Linkers
    linkers = construct["linkers"]
    adjuvant_linker = linkers["adjuvant_helper"]
    padre_linker = linkers["helper_mhc1"]

    # 4. Construct Domain Sub-chains using explicit joint mechanics
    mhc1_chain = linkers["mhc1"].join(mhc1_peps) if mhc1_peps else ""
    mhc2_chain = linkers["mhc2"].join(mhc2_peps) if mhc2_peps else ""
    bcell_chain = linkers["bcell"].join(bcell_peps) if bcell_peps else ""

    # 5. Cascaded Vaccine String Assembly
    # Joins the segments strictly in the requested engineering order
    final_vaccine_sequence = (
        adjuvant + 
        adjuvant_linker + 
        padre + 
        padre_linker + 
        mhc1_chain + 
        linkers["mhc1_mhc2"] +
        mhc2_chain + 
        linkers["mhc2_bcell"] +
        bcell_chain
    )

    # 6. Export to production FASTA format
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    
    with open(args.output, "w") as f:
        f.write(">CPV2_FliC_PADRE_MultiEpitope_Vaccine_Construct\n")
        # Break sequence strings into clean, standardized 60-character FASTA lines
        for i in range(0, len(final_vaccine_sequence), 60):
            f.write(final_vaccine_sequence[i:i+60] + "\n")

    print("[SUCCESS] Vaccine construct assembled with absolute architectural precision!")
    print(f"[SUCCESS] Total Chimera Sequence Length: {len(final_vaccine_sequence)} AA")
    print(f"[SUCCESS] Clean FASTA sequence written to: {args.output}")

if __name__ == "__main__":
    main()