#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import sys
import os
from pipeline_config import load_config, resolve_path

def parse_args():
    parser = argparse.ArgumentParser(
        description="Canine Vaccine Prioritization Engine: CPV2 Production Edition"
    )
    parser.add_argument("--mhc", required=True, help="Path to MHC Class I binders / Core base dataframe")
    parser.add_argument("--topology", required=True, help="DeepTMHMM membrane topology file")
    parser.add_argument("--conservation", required=True, help="Sequence conservation tracking file")
    parser.add_argument("--blast", required=True, help="Host blast identity exclusions file")
    parser.add_argument("--immunogenicity", required=True, help="Integrated Ensemble immunogenicity file")
    parser.add_argument("--toxicity", required=True, help="ToxinPred3 output predictions")
    parser.add_argument("--allergenicity", required=True, help="AlgPred2 allergen profiles")
    parser.add_argument("--config", required=True, help="Path to pipeline config.yaml containing weights")
    parser.add_argument("--output", required=True, help="Target final path for ranked T-Cell CSV export")
    return parser.parse_args()

def main():
    args = parse_args()

    # 1. Load config weights safely
    config = load_config(args.config)
    for argument in (
        "mhc", "topology", "conservation", "blast", "immunogenicity",
        "toxicity", "allergenicity", "output",
    ):
        setattr(args, argument, resolve_path(getattr(args, argument), config))
    
    weights = config["ranking"]["weights"]

    # 2. Load Core MHC Dataframe
    print("[PROCESSING] Loading core MHC binders...")
    mhc_df = pd.read_csv(args.mhc)
    mhc_df["peptide_upper"] = mhc_df["peptide"].astype(str).str.upper().str.strip()

    # 3. Load Membrane Topologies
    print("[PROCESSING] Mapping DeepTMHMM membrane topologies...")
    topo_map = {}
    if os.path.exists(args.topology):
        try:
            topo_df = pd.read_csv(args.topology)
            if "peptide" in topo_df.columns:
                topo_map = dict(zip(topo_df["peptide"].astype(str).str.upper().str.strip(), topo_df.iloc[:, -1]))
        except Exception:
            pass

    # 4. Load Host BLAST Identity Exclusions
    print("[PROCESSING] Parsing host homology screens...")
    blast_df = pd.read_csv(args.blast)
    
    blast_id_col = next((c for c in ["protein_id", "id", "ID", "sample", "sample_id", "accession"] if c in blast_df.columns), None)
    host_match_col = next((c for c in ["host_match", "identity", "alignment_length"] if c in blast_df.columns), None)
    
    if host_match_col and blast_id_col:
        has_homology = blast_df[(blast_df[host_match_col].notna()) & (blast_df[host_match_col] != "None")]
        host_homologue_proteins = set(has_homology[blast_id_col].unique())
    else:
        host_homologue_proteins = set()

    # 5. Load Ensemble Immunogenicity Scores
    print("[PROCESSING] Parsing ensemble immunogenicity markers...")
    ens_df = pd.read_csv(args.immunogenicity)
    ens_score_col = next((c for c in ["ensemble_immunogenicity_score", "score", "immunogenicity_score"] if c in ens_df.columns), ens_df.columns[-1])
    ens_map = dict(zip(ens_df["peptide"].astype(str).str.upper().str.strip(), ens_df[ens_score_col]))

    # 6. Load Toxicity Profiles (Checking verified column: toxicity_prediction)
    print("[PROCESSING] Checking toxin profiles...")
    tox_df = pd.read_csv(args.toxicity)
    pred_tox_col = "toxicity_prediction" if "toxicity_prediction" in tox_df.columns else tox_df.columns[1]
    toxic_peptides = set(tox_df[tox_df[pred_tox_col].astype(str).str.lower().str.contains("toxin") & 
                                ~tox_df[pred_tox_col].astype(str).str.lower().str.contains("non")]["peptide"].astype(str).str.upper().str.strip())

    # 7. Load Allergenicity Profiles (Checking verified column: allergenicity_prediction)
    print("[PROCESSING] Checking allergen profiles...")
    alg_df = pd.read_csv(args.allergenicity)
    pred_alg_col = "allergenicity_prediction" if "allergenicity_prediction" in alg_df.columns else alg_df.columns[1]
    allergen_peptides = set(alg_df[alg_df[pred_alg_col].astype(str).str.lower().str.contains("allergen") & 
                                   ~alg_df[pred_alg_col].astype(str).str.lower().str.contains("non")]["peptide"].astype(str).str.upper().str.strip())

    # 8. Load Sequence Conservation Tracking Matrix
    print("[PROCESSING] Processing alignment positional conservation matrix...")
    cons_df = pd.read_csv(args.conservation)
    pos_col = next((c for c in ["position", "pos", "Pos"] if c in cons_df.columns), cons_df.columns[0])
    window_column = f"window_{config['conservation']['window_size']}mer_avg"
    val_col = next((c for c in [window_column, "score", "conservation"] if c in cons_df.columns), cons_df.columns[-1])
    pos_scores = dict(zip(cons_df[pos_col].astype(int), cons_df[val_col].astype(float)))

    # 9. Assemble and Score Final Combined Dataframe
    final_df = mhc_df.copy()
    
    cons_scores = []
    start_col = next((c for c in ["start", "pos", "position"] if c in final_df.columns), None)
    end_col = next((c for c in ["end", "stop"] if c in final_df.columns), None)

    for _, row in final_df.iterrows():
        try:
            if start_col and end_col:
                start_pos, end_pos = int(row[start_col]), int(row[end_col])
                window_vals = [pos_scores[p] for p in range(start_pos, end_pos + 1) if p in pos_scores]
                cons_scores.append(np.mean(window_vals) if window_vals else 1.0)
            else:
                cons_scores.append(1.0)
        except Exception:
            cons_scores.append(1.0)
            
    final_df["conservation_score"] = cons_scores
    final_df["topology_score"] = final_df["peptide_upper"].map(topo_map).fillna(1.0)
    final_df["ensemble_immunogenicity"] = final_df["peptide_upper"].map(ens_map).fillna(0.0)
    
    # Check host homologue safely
    mhc_id_col = next((c for c in ["protein_id", "id", "ID", "sample", "sample_id", "accession"] if c in final_df.columns), None)
    if mhc_id_col:
        final_df["is_host_homologue"] = final_df[mhc_id_col].isin(host_homologue_proteins)
    else:
        final_df["is_host_homologue"] = False
        
    final_df["is_toxic"] = final_df["peptide_upper"].isin(toxic_peptides)
    final_df["is_allergen"] = final_df["peptide_upper"].isin(allergen_peptides)
    final_df["passed_safety_filter"] = ~(final_df["is_host_homologue"] | final_df["is_toxic"] | final_df["is_allergen"])
    final_df["safety_score"] = final_df["passed_safety_filter"].astype(float)

    # Dynamic column resolver for MHC Percentile Rank (Prevents KeyError)
    rank_col = next((c for c in ["percentile_rank", "rank", "Rank", "%Rank", "score_rank", "ic50"] if c in final_df.columns), None)

    if "mhc1_score" not in final_df.columns:
        if rank_col:
            # If rank is based on IC50 or binding score directly where higher is better, logic alters. For standard %Rank, lower is better:
            if "score" in rank_col.lower():
                final_df["mhc1_score"] = final_df[rank_col].fillna(0.0)
            else:
                final_df["mhc1_score"] = 1 - (final_df[rank_col].fillna(100.0) / 100.0)
        else:
            final_df["mhc1_score"] = 0.5  # Safe default fallback boundary

    # Map secondary features if columns exist in immunogenicity ensemble file
    if "mhc2_score" not in final_df.columns and "peptide" in ens_df.columns:
        mhc2_c = next((c for c in ["mhc2_score", "mhc2"] if c in ens_df.columns), None)
        if mhc2_c:
            final_df["mhc2_score"] = final_df["peptide_upper"].map(dict(zip(ens_df["peptide"].astype(str).str.upper(), ens_df[mhc2_c]))).fillna(0.0)
        else:
            final_df["mhc2_score"] = 0.0

    if "processing_efficiency" not in final_df.columns and "peptide" in ens_df.columns:
        proc_c = next((c for c in ["proc_score", "processing_efficiency", "cleavage_score"] if c in ens_df.columns), None)
        if proc_c:
            final_df["processing_efficiency"] = final_df["peptide_upper"].map(dict(zip(ens_df["peptide"].astype(str).str.upper(), ens_df[proc_c]))).fillna(0.0)
        else:
            final_df["processing_efficiency"] = 0.0

    # 10. Multi-Criteria Optimization Weight Scoring
    mhc2_key = "weights" if "weights" in weights else "mhc2_binding"
    
    final_df["final_vaccine_score"] = (
        (final_df["mhc1_score"] * weights.get("mhc_binding", 0.25)) +
        (final_df.get("mhc2_score", 0.0) * weights.get(mhc2_key, 0.20)) +
        (final_df.get("processing_efficiency", 0.0) * weights.get("processing_efficiency", 0.15)) +
        (final_df["conservation_score"] * weights.get("conservation", 0.15)) +
        (final_df["ensemble_immunogenicity"] * weights.get("immunogenicity", 0.15)) +
        (final_df["safety_score"] * weights.get("safety", 0.10))
    )

    final_df.drop(columns=["peptide_upper"], inplace=True, errors='ignore')

    # Sort candidates by safety filter first, then highest vaccine score
    final_df = final_df.sort_values(by=["passed_safety_filter", "final_vaccine_score"], ascending=[False, False]).reset_index(drop=True)
    final_df.insert(0, "prioritization_rank", final_df.index + 1)
    
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    final_df.to_csv(args.output, index=False)
    print(f"[SUCCESS] Exported {len(final_df)} vaccine target candidates safely to: {args.output}")

if __name__ == "__main__":
    main()