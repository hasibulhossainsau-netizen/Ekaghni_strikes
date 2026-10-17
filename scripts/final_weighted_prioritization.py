#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import sys
import os
from pipeline_config import load_config, require_columns, resolve_path


def rank_candidates(
    mhc_df, topology_df, conservation_df, blast_df, ens_df, tox_df, alg_df,
    weights, conservation_window,
):
    """Validate upstream schemas, calculate features, and rank MHC-I candidates."""
    required_schemas = (
        (mhc_df, ("protein_id", "peptide", "start", "end", "percentile_rank"), "MHC-I predictions"),
        (topology_df, ("peptide", "topology_score"), "topology predictions"),
        (conservation_df, ("position", conservation_window), "conservation scores"),
        (blast_df, ("protein_id", "host_match"), "host homology results"),
        (
            ens_df,
            ("peptide", "mhc2_score", "proc_score", "ensemble_immunogenicity_score"),
            "ensemble immunogenicity scores",
        ),
        (tox_df, ("peptide", "toxicity_prediction"), "toxicity predictions"),
        (alg_df, ("peptide", "allergenicity_prediction"), "allergenicity predictions"),
    )
    for dataframe, columns, name in required_schemas:
        require_columns(dataframe, columns, name)
    required_weights = {
        "mhc_binding", "mhc2_binding", "processing_efficiency",
        "conservation", "immunogenicity", "safety",
    }
    missing_weights = required_weights - set(weights)
    if missing_weights:
        raise ValueError(f"Ranking weights are missing: {', '.join(sorted(missing_weights))}.")

    final_df = mhc_df.copy()
    final_df["peptide_upper"] = final_df["peptide"].astype(str).str.upper().str.strip()
    final_df["mhc1_score"] = 1.0 - (
        pd.to_numeric(final_df["percentile_rank"], errors="raise") / 100.0
    )

    topology_map = dict(zip(
        topology_df["peptide"].astype(str).str.upper().str.strip(),
        pd.to_numeric(topology_df["topology_score"], errors="raise"),
    ))
    final_df["topology_score"] = final_df["peptide_upper"].map(topology_map)

    position_scores = dict(zip(
        pd.to_numeric(conservation_df["position"], errors="raise").astype(int),
        pd.to_numeric(conservation_df[conservation_window], errors="raise"),
    ))
    conservation_scores = []
    for start, end in zip(final_df["start"], final_df["end"]):
        positions = range(int(start), int(end) + 1)
        scores = [position_scores[position] for position in positions if position in position_scores]
        if not scores:
            raise ValueError(f"No conservation scores found for candidate range {start}-{end}.")
        conservation_scores.append(float(np.mean(scores)))
    final_df["conservation_score"] = conservation_scores

    homologues = set(
        blast_df.loc[
            blast_df["host_match"].notna()
            & (blast_df["host_match"].astype(str) != "None"),
            "protein_id",
        ].astype(str)
    )
    final_df["is_host_homologue"] = final_df["protein_id"].astype(str).isin(homologues)

    toxic = set(
        tox_df.loc[
            tox_df["toxicity_prediction"].astype(str).str.lower().str.contains("toxic")
            & ~tox_df["toxicity_prediction"].astype(str).str.lower().str.contains("non"),
            "peptide",
        ].astype(str).str.upper().str.strip()
    )
    allergens = set(
        alg_df.loc[
            alg_df["allergenicity_prediction"].astype(str).str.lower().str.contains("allergen")
            & ~alg_df["allergenicity_prediction"].astype(str).str.lower().str.contains("non"),
            "peptide",
        ].astype(str).str.upper().str.strip()
    )
    final_df["is_toxic"] = final_df["peptide_upper"].isin(toxic)
    final_df["is_allergen"] = final_df["peptide_upper"].isin(allergens)
    final_df["passed_safety_filter"] = ~(
        final_df["is_host_homologue"] | final_df["is_toxic"] | final_df["is_allergen"]
    )
    final_df["safety_score"] = final_df["passed_safety_filter"].astype(float)

    ensemble = ens_df.copy()
    ensemble["peptide_upper"] = ensemble["peptide"].astype(str).str.upper().str.strip()
    ensemble = ensemble.drop_duplicates("peptide_upper").set_index("peptide_upper")
    final_df["mhc2_score"] = final_df["peptide_upper"].map(ensemble["mhc2_score"]).fillna(0.0)
    final_df["processing_efficiency"] = final_df["peptide_upper"].map(ensemble["proc_score"]).fillna(0.0)
    final_df["ensemble_immunogenicity"] = final_df["peptide_upper"].map(
        ensemble["ensemble_immunogenicity_score"]
    ).fillna(0.0)

    final_df["final_vaccine_score"] = sum(
        final_df[column] * weights[weight]
        for column, weight in (
            ("mhc1_score", "mhc_binding"),
            ("mhc2_score", "mhc2_binding"),
            ("processing_efficiency", "processing_efficiency"),
            ("conservation_score", "conservation"),
            ("ensemble_immunogenicity", "immunogenicity"),
            ("safety_score", "safety"),
        )
    )
    final_df.drop(columns=["peptide_upper"], inplace=True)
    final_df = final_df.sort_values(
        by=["passed_safety_filter", "final_vaccine_score"],
        ascending=[False, False],
    ).reset_index(drop=True)
    final_df.insert(0, "prioritization_rank", final_df.index + 1)
    return final_df

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

    window_column = f"window_{config['conservation']['window_size']}mer_avg"
    final_df = rank_candidates(
        pd.read_csv(args.mhc),
        pd.read_csv(args.topology),
        pd.read_csv(args.conservation),
        pd.read_csv(args.blast),
        pd.read_csv(args.immunogenicity),
        pd.read_csv(args.toxicity),
        pd.read_csv(args.allergenicity),
        weights,
        window_column,
    )
    
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    final_df.to_csv(args.output, index=False)
    print(f"[SUCCESS] Exported {len(final_df)} vaccine target candidates safely to: {args.output}")

if __name__ == "__main__":
    main()