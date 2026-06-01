#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import os
import sys

def parse_args():
    parser = argparse.ArgumentParser(
        description="Consensus Epitope Intersector and Multi-Antigen Mapping Engine"
    )
    parser.add_argument("--mhc1", required=True, help="Ranked MHC Class I epitopes CSV")
    parser.add_argument("--mhc2", required=True, help="IEDB MHC Class II predictions CSV")
    parser.add_argument("--bcell", required=True, help="B-cell linear epitopes summary CSV")
    parser.add_argument("--tox", required=True, help="ToxinPred3 output CSV")
    parser.add_argument("--alg", required=True, help="AlgPred2 output CSV")
    parser.add_argument("--cons", required=True, help="Positional conservation CSV")
    parser.add_argument("--blast", required=True, help="Host blast identity exclusion CSV")
    parser.add_argument("--out_dir", required=True, help="Target directory for consensus outputs")
    return parser.parse_args()

def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    print("[PROCESSING] Loading safety and exclusion profiles...")
    
    # 1. Load Toxicity Profiles
    tox_df = pd.read_csv(args.tox)
    tox_col = "toxicity_prediction" if "toxicity_prediction" in tox_df.columns else "prediction"
    tox_set = set(
        tox_df[tox_df[tox_col].astype(str).str.lower().str.contains("toxic") & 
               ~tox_df[tox_col].astype(str).str.lower().str.contains("non")]["peptide"]
        .astype(str).str.upper().str.strip()
    )

    # 2. Load Allergenicity Profiles
    alg_df = pd.read_csv(args.alg)
    alg_col = "allergenicity_prediction" if "allergenicity_prediction" in alg_df.columns else "prediction"
    alg_set = set(
        alg_df[alg_df[alg_col].astype(str).str.lower().str.contains("allergen") & 
               ~alg_df[alg_col].astype(str).str.lower().str.contains("non")]["peptide"]
        .astype(str).str.upper().str.strip()
    )

    # 3. Load Host Homology Protein Exclusions
    blast_df = pd.read_csv(args.blast)
    has_homology = blast_df[(blast_df["host_match"].notna()) & (blast_df["host_match"] != "None")]
    excluded_proteins = set(has_homology["protein_id"].unique())

    # Helper function to screen structural safety
    def is_safe(peptide, protein_id=None):
        pep_str = str(peptide).upper().strip()
        if pep_str in tox_set or pep_str in alg_set:
            return False
        if protein_id and protein_id in excluded_proteins:
            return False
        return True

    print("[PROCESSING] Filtering and intersecting vaccine lineages...")

    # 4. Filter MHC I Candidates
    mhc1 = pd.read_csv(args.mhc1)
    mhc1["passed_consensus_safety"] = mhc1.apply(
        lambda row: is_safe(row["peptide"], row.get("protein_id")), axis=1
    )
    mhc1_filtered = mhc1[mhc1["passed_consensus_safety"] == True].copy()
    mhc1_filtered.to_csv(os.path.join(args.out_dir, "filtered_mhc1_epitopes.csv"), index=False)

    # 5. Filter MHC II Predictions
    mhc2 = pd.read_csv(args.mhc2)
    mhc2_pep_col = "peptide" if "peptide" in mhc2.columns else mhc2.columns[0]
    mhc2["passed_consensus_safety"] = mhc2.apply(
        lambda row: is_safe(row[mhc2_pep_col], row.get("protein_id")), axis=1
    )
    mhc2_filtered = mhc2[mhc2["passed_consensus_safety"] == True].copy()
    mhc2_filtered.to_csv(os.path.join(args.out_dir, "filtered_mhc2_epitopes.csv"), index=False)

    # 6. Filter B-Cell Linear Tracks
    bcell = pd.read_csv(args.bcell)
    bcell_pep_col = "peptide" if "peptide" in bcell.columns else ("sequence" if "sequence" in bcell.columns else bcell.columns[0])
    bcell["passed_consensus_safety"] = bcell.apply(
        lambda row: is_safe(row[bcell_pep_col], row.get("protein_id")), axis=1
    )
    bcell_filtered = bcell[bcell["passed_consensus_safety"] == True].copy()
    bcell_filtered.to_csv(os.path.join(args.out_dir, "filtered_bcell_epitopes.csv"), index=False)

    # 7. Generate Overlapping Coordinate Maps (Consensus Blocks)
    print("[PROCESSING] Constructing regional coordinate intersections...")
    
    # Track hotspots along the sequence lengths
    cons_df = pd.read_csv(args.cons)
    max_position = int(cons_df["position"].max())
    
    mhc1_coverage = np.zeros(max_position + 2)
    mhc2_coverage = np.zeros(max_position + 2)
    bcell_coverage = np.zeros(max_position + 2)

    for _, row in mhc1_filtered.iterrows():
        if "start" in row and "end" in row:
            mhc1_coverage[int(row["start"]):int(row["end"]) + 1] = 1

    for _, row in mhc2_filtered.iterrows():
        if "start" in row and "end" in row:
            mhc2_coverage[int(row["start"]):int(row["end"]) + 1] = 1

    for _, row in bcell_filtered.iterrows():
        if "start" in row and "end" in row:
            bcell_coverage[int(row["start"]):int(row["end"]) + 1] = 1

    # Overlaps where multiple epitope frameworks hit the same residue window
    overlapping_regions = []
    in_region = False
    start_idx = None

    for pos in range(1, max_position + 1):
        # A position scores if hit by at least two distinct structural lineages
        score = int(mhc1_coverage[pos]) + int(mhc2_coverage[pos]) + int(bcell_coverage[pos])
        if score >= 2:
            if not in_region:
                in_region = True
                start_idx = pos
        else:
            if in_region:
                overlapping_regions.append({"start": start_idx, "end": pos - 1, "length": (pos - start_idx)})
                in_region = False

    if in_region:
        overlapping_regions.append({"start": start_idx, "end": max_position, "length": (max_position - start_idx) + 1})

    overlap_df = pd.DataFrame(overlapping_regions)
    if overlap_df.empty:
        overlap_df = pd.DataFrame(columns=["start", "end", "length"])
        
    overlap_df.to_csv(os.path.join(args.out_dir, "overlapping_consensus_regions.csv"), index=False)
    
    print(f"[SUCCESS] Consensus intersections written to {args.out_dir}/")

if __name__ == "__main__":
    main()