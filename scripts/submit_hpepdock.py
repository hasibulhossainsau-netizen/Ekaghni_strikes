#!/usr/bin/env python3
import os
import time
import argparse
import pandas as pd
import requests

def parse_args():
    parser = argparse.ArgumentParser(description="Submit validated pipeline epitopes to HPEPDOCK server with resume support.")
    parser.add_argument("--csv", required=True, help="Path to clean consensus MHC-I target file")
    parser.add_argument("--pdb-dir", required=True, help="Relative path to local directory containing receptor crystal structures")
    parser.add_argument("--output", required=True, help="Output submission log tracking path")
    parser.add_argument("--url", required=True, help="Target HPEPDOCK operational endpoint URL")
    parser.add_argument("--email", required=True, help="User context tracking notifications address")
    parser.add_argument("--retries", type=int, default=3, help="Maximum structural dispatch timeout fallback steps")
    parser.add_argument("--receptors", nargs="+", required=True, help="Space-separated list of target receptor codes (e.g. 5F1I 7CJQ)")
    return parser.parse_args()

def verify_and_get_pdb(pdb_id, pdb_dir):
    """Verifies existence of PDB inside structural data context, fetches via RCSB fallback if missing."""
    pdb_path = os.path.join(pdb_dir, f"{pdb_id}.pdb")
    if not os.path.exists(pdb_path):
        print(f"Structure {pdb_id}.pdb not found in local path '{pdb_dir}'. Attempting remote fetching from RCSB...")
        os.makedirs(pdb_dir, exist_ok=True)
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        response = requests.get(url, timeout=30)
        if response.status_code == 200:
            with open(pdb_path, "w") as f:
                f.write(response.text)
            print(f"Successfully downloaded structural file {pdb_id}.pdb to {pdb_path}")
        else:
            raise FileNotFoundError(f"Could not locate or download structural file context for reference sequence: {pdb_id}")
    return pdb_path

def submit_to_hpepdock(receptor_pdb_path, peptide_seq, job_name, pdb_id, server_url, email, max_retries=3):
    """Submits the structural screening job to the target HPEPDOCK server instance using robust connection context handling."""
    print(f"Submitting computational docking run: '{job_name}' (Target Core Sequence: {peptide_seq})...")

    data = {
        "jobname": job_name,
        "email": email,
        "pdbid1": pdb_id,
        "fastaseq2": peptide_seq,
        "docktyp": "dock",
        "nmodpep": "1000",
        "upload": "Submit"
    }

    with open(receptor_pdb_path, "rb") as f:
        files = {
            "pdbfile1": (os.path.basename(receptor_pdb_path), f, "application/octet-stream")
        }

        for attempt in range(1, max_retries + 1):
            try:
                response = requests.post(server_url, data=data, files=files, timeout=180)
                if response.status_code == 200:
                    print(f"Job context '{job_name}' accepted on target framework loop {attempt}.")
                    return response.text
                else:
                    print(f"Endpoint warning: received server status code {response.status_code} during retry cycle {attempt}.")
            except requests.exceptions.Timeout:
                print(f"Connection timeout occurred on attempt {attempt}. Retrying execution block in 15 seconds...")
                time.sleep(15)
            except Exception as e:
                print(f"Runtime execution error caught on attempt {attempt}: {e}")
                time.sleep(15)

        print(f"Critical Error: Execution failed across all {max_retries} connection loops for run context reference: '{job_name}'.")
        return None

def main():
    args = parse_args()
    
    # Ensure directory tracking structures exist safely
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        
    # --- RESUME MECHANISM: Load state tracking context ---
    completed_jobs = set()
    existing_records = []
    
    if os.path.exists(args.output) and os.path.getsize(args.output) > 0:
        try:
            old_df = pd.read_csv(args.output)
            # Filter to find rows that genuinely succeeded
            success_df = old_df[old_df['Status'].str.strip().str.lower() == 'submitted']
            for _, row in success_df.iterrows():
                # Store structural key tuple as state tracker
                completed_jobs.add((str(row['Receptor']).strip(), str(row['Peptide']).strip()))
            
            # Keep track of everything else to append safely or reconstruct
            existing_records = old_df.to_dict(orient='records')
            print(f"State recovery found: Loaded {len(completed_jobs)} jobs previously 'Submitted'. These will be skipped.")
        except Exception as e:
            print(f"Could not parse existing tracking file due to formatting, starting crisp execution log. Detail: {e}")

    df = pd.read_csv(args.csv)
    pep_col = 'peptide' if 'peptide' in df.columns else df.columns[0]
    
    # Slice the unique array context here to restrict analysis to the top 5 values
    top_peptides = df[pep_col].unique()[:5]
    
    print(f"Restricted evaluation to the top {len(top_peptides)} unique candidate epitopes for structural alignment checks.")

    # Locate/Download the receptors using the relative directory specified
    receptor_files = {}
    for pdb_id in args.receptors:
        receptor_files[pdb_id] = verify_and_get_pdb(pdb_id, args.pdb_dir)

    # Initialize tracking data structures across execution loops
    tracked_runs = { (str(r['Receptor']).strip(), str(r['Peptide']).strip()): r for r in existing_records }

    for pdb_id, pdb_path in receptor_files.items():
        for i, pep in enumerate(top_peptides, start=1):
            
            # Evaluate against state memory
            if (str(pdb_id), str(pep)) in completed_jobs:
                print(f"Skipping: {pdb_id} + Peptide {pep} (Already successfully submitted).")
                continue
                
            job_name = f"DLA_{pdb_id}_Pep_{i}"
            response_text = submit_to_hpepdock(
                receptor_pdb_path=pdb_path, 
                peptide_seq=pep, 
                job_name=job_name, 
                pdb_id=pdb_id, 
                server_url=args.url, 
                email=args.email, 
                max_retries=args.retries
            )

            if response_text and "error" not in response_text.lower():
                status_state = "Submitted"
            else:
                status_state = "Failed"
            
            # Save or update state dict
            tracked_runs[(str(pdb_id), str(pep))] = {
                "Job_Name": job_name,
                "Receptor": pdb_id,
                "Peptide": pep,
                "Status": status_state
            }
            
            # Real-time state flushing so data isn't lost if interrupted again
            updated_df = pd.DataFrame(tracked_runs.values())
            updated_df[["Job_Name", "Receptor", "Peptide", "Status"]].to_csv(args.output, index=False)

            time.sleep(5) # Prevent webserver rate-limiting triggers

    print(f"\nPipeline sync completed successfully. Execution matrix state preserved in log: '{args.output}'")

if __name__ == "__main__":
    main()