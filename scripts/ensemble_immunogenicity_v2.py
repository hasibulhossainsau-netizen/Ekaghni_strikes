#!/usr/bin/env python3
import argparse
import pandas as pd
import sys
from pipeline_config import load_config, require_columns, resolve_path

def parse_args():
    parser = argparse.ArgumentParser(description="Ensemble Immunogenicity Scoring Framework v2 (IEDB Integrated)")
    parser.add_argument("--mhc2-scores", required=True, help="Path to IEDB MHC-II predictions CSV")
    parser.add_argument("--processing", required=True, help="Path to IEDB MHC-I processing CSV")
    parser.add_argument("--breed-weights", required=True, help="Path to breed weights config CSV")
    parser.add_argument("--output", required=True, help="Path to save ensemble immunogenicity output CSV")
    parser.add_argument("--weight-mhc2", type=float, default=None, help="Weight for MHC-II predictions")
    parser.add_argument("--weight-processing", type=float, default=None, help="Weight for MHC-I processing score")
    parser.add_argument("--weight-population", type=float, default=None, help="Weight for Breed Population/Allele coverage")
    parser.add_argument("--config", default=None, help="Path to pipeline config.yaml")
    return parser.parse_args()


def main():
    args = parse_args()
    config = load_config(args.config)
    args.mhc2_scores = resolve_path(args.mhc2_scores, config)
    args.processing = resolve_path(args.processing, config)
    args.breed_weights = resolve_path(args.breed_weights, config)
    args.output = resolve_path(args.output, config)
    args.weight_mhc2 = args.weight_mhc2 if args.weight_mhc2 is not None else config["ensemble"]["weight_mhc2"]
    args.weight_processing = args.weight_processing if args.weight_processing is not None else config["ensemble"]["weight_processing"]
    args.weight_population = args.weight_population if args.weight_population is not None else config["ensemble"]["weight_population_coverage"]
    
    try:
        mhc2_df = pd.read_csv(args.mhc2_scores)
        proc_df = pd.read_csv(args.processing)
        breed_df = pd.read_csv(args.breed_weights)
    except Exception as e:
        print(f"Error reading input files: {e}", file=sys.stderr)
        sys.exit(1)

    require_columns(mhc2_df, ("peptide", "percentile_rank"), "IEDB MHC-II predictions")
    require_columns(proc_df, ("peptide", "processing_score"), "IEDB processing predictions")
    require_columns(breed_df, ("weight",), "breed weights")

    mhc2_df['peptide'] = mhc2_df['peptide'].astype(str).str.strip()
    proc_df['peptide'] = proc_df['peptide'].astype(str).str.strip()

    mhc2_df["percentile_rank"] = pd.to_numeric(mhc2_df["percentile_rank"], errors="raise")
    proc_df["processing_score"] = pd.to_numeric(proc_df["processing_score"], errors="raise")

    mhc2_df["mhc2_score"] = 1.0 - (mhc2_df["percentile_rank"] / 100.0)
    mhc2_df = mhc2_df.drop_duplicates(subset=["peptide"])
    
    mhc2_grouped = mhc2_df.groupby('peptide', as_index=False)['mhc2_score'].mean()

    proc_min = proc_df["processing_score"].min()
    proc_max = proc_df["processing_score"].max()
    if proc_max != proc_min:
        proc_df["proc_score"] = (
            (proc_df["processing_score"] - proc_min) / (proc_max - proc_min)
        )
    else:
        proc_df["proc_score"] = 1.0
    proc_df = proc_df.drop_duplicates(subset=["peptide"])
        
    proc_grouped = proc_df.groupby('peptide', as_index=False)['proc_score'].mean()

    # 3. merge
    master_peptides = pd.merge(mhc2_grouped, proc_grouped, on='peptide', how='outer').fillna(0)

    # ৪. Population Score
    if breed_df.empty:
        raise ValueError("Breed weights table must contain at least one weight.")
    pop_factor = pd.to_numeric(breed_df["weight"], errors="raise").mean()
    master_peptides['population_score'] = pop_factor

    # ৫. Ensemble Score
    master_peptides['ensemble_immunogenicity_score'] = (
        (master_peptides['mhc2_score'] * args.weight_mhc2) +
        (master_peptides['proc_score'] * args.weight_processing) +
        (master_peptides['population_score'] * args.weight_population)
    )

    # save
    master_peptides.to_csv(args.output, index=False)
    print(f"✅ Ensemble immunogenicity scoring completed successfully!")
    print(f"   Output: {args.output}")
    print(f"   Total peptides: {len(master_peptides)}")


if __name__ == "__main__":
    main()