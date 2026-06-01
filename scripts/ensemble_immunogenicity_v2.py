#!/usr/bin/env python3
import argparse
import pandas as pd
import sys

def parse_args():
    parser = argparse.ArgumentParser(description="Ensemble Immunogenicity Scoring Framework v2 (IEDB Integrated)")
    parser.add_argument("--mhc2-scores", required=True, help="Path to IEDB MHC-II predictions CSV")
    parser.add_argument("--processing", required=True, help="Path to IEDB MHC-I processing CSV")
    parser.add_argument("--breed-weights", required=True, help="Path to breed weights config CSV")
    parser.add_argument("--output", required=True, help="Path to save ensemble immunogenicity output CSV")
    parser.add_argument("--weight-mhc2", type=float, default=0.40, help="Weight for MHC-II predictions")
    parser.add_argument("--weight-processing", type=float, default=0.35, help="Weight for MHC-I processing score")
    parser.add_argument("--weight-population", type=float, default=0.25, help="Weight for Breed Population/Allele coverage")
    return parser.parse_args()


def standardize_peptide_column(df):
    """Standardizing peptide column"""
    possible_names = ['peptide', 'sequence', 'description', 'seq', 'peptide_seq', 'epitope', 'Peptide']
    for col in df.columns:
        if str(col).lower() in [x.lower() for x in possible_names] or any(p in str(col).lower() for p in ['seq', 'pep']):
            df.rename(columns={col: 'peptide'}, inplace=True)
            return True
    return False


def main():
    args = parse_args()
    
    try:
        mhc2_df = pd.read_csv(args.mhc2_scores)
        proc_df = pd.read_csv(args.processing)
        breed_df = pd.read_csv(args.breed_weights)
    except Exception as e:
        print(f"Error reading input files: {e}", file=sys.stderr)
        sys.exit(1)

    # column standardization
    mhc2_ok = standardize_peptide_column(mhc2_df)
    proc_ok = standardize_peptide_column(proc_df)
            
    if not mhc2_ok or not proc_ok:
        print(f"Error: Could not find peptide column.", file=sys.stderr)
        print(f"MHC2 columns: {list(mhc2_df.columns)}", file=sys.stderr)
        print(f"Proc columns: {list(proc_df.columns)}", file=sys.stderr)
        sys.exit(1)

    # removal of duplicate columns
    mhc2_df = mhc2_df.loc[:, ~mhc2_df.columns.duplicated(keep='first')]
    proc_df = proc_df.loc[:, ~proc_df.columns.duplicated(keep='first')]

    # converting into str
    mhc2_df['peptide'] = mhc2_df['peptide'].astype(str).str.strip()
    proc_df['peptide'] = proc_df['peptide'].astype(str).str.strip()

    # removal of duplicate columns
    mhc2_df = mhc2_df.drop_duplicates(subset=['peptide'])
    proc_df = proc_df.drop_duplicates(subset=['peptide'])

    # 1. MHC-II score
    if 'percentile_rank' in mhc2_df.columns:
        mhc2_df['mhc2_score'] = 1 - (mhc2_df['percentile_rank'] / 100.0)
    else:
        mhc2_score_col = [c for c in mhc2_df.columns if any(x in c.lower() for x in ['score', 'rank', 'ic50', 'affinity'])]
        if mhc2_score_col:
            col = mhc2_score_col[0]
            mhc2_df['mhc2_score'] = (mhc2_df[col] - mhc2_df[col].min()) / (mhc2_df[col].max() - mhc2_df[col].min() + 1e-9)
        else:
            mhc2_df['mhc2_score'] = 1.0
    
    mhc2_grouped = mhc2_df.groupby('peptide', as_index=False)['mhc2_score'].mean()

    # 2. MHC-I Processing score
    proc_score_col = [c for c in proc_df.columns if any(x in c.lower() for x in ['score', 'processing', 'total', 'tap'])]
    if proc_score_col:
        proc_df['proc_score'] = proc_df[proc_score_col[0]]
    else:
        proc_df['proc_score'] = 1.0
        
    if proc_df['proc_score'].max() != proc_df['proc_score'].min():
        proc_df['proc_score'] = (proc_df['proc_score'] - proc_df['proc_score'].min()) / (proc_df['proc_score'].max() - proc_df['proc_score'].min() + 1e-9)
    else:
        proc_df['proc_score'] = 1.0
        
    proc_grouped = proc_df.groupby('peptide', as_index=False)['proc_score'].mean()

    # 3. merge
    master_peptides = pd.merge(mhc2_grouped, proc_grouped, on='peptide', how='outer').fillna(0)

    # ৪. Population Score
    if not breed_df.empty and 'weight' in breed_df.columns and len(breed_df) > 0:
        pop_factor = breed_df['weight'].mean()
    else:
        pop_factor = 1.0
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