#!/usr/bin/env python3
import os
import time
import argparse
import pandas as pd
import requests
from pipeline_config import load_config, resolve_path

def parse_args():
    parser = argparse.ArgumentParser(description="Submit validated pipeline epitopes to HPEPDOCK server with resume support.")
    parser.add_argument("--csv", required=True, help="Path to clean consensus MHC-I target file")
    parser.add_argument("--pdb-dir", help="Path to local directory containing receptor crystal structures")
    parser.add_argument("--output", required=True, help="Output submission log tracking path")
    parser.add_argument("--url", help="HPEPDOCK endpoint URL (defaults to config.yaml)")
    parser.add_argument("--email", help="Notification email (defaults to config.yaml)")
    parser.add_argument("--retries", type=int, help="Maximum submission attempts (defaults to config.yaml)")
    parser.add_argument("--receptors", nargs="+", help="Receptor IDs (defaults to config.yaml)")
    parser.add_argument("--config", default=None, help="Path to pipeline config.yaml")
    return parser.parse_args()

def verify_and_get_pdb(pdb_id, pdb_dir, download_url, timeout):
    """Verifies existence of PDB inside structural data context, fetches via RCSB fallback if missing."""
    pdb_path = os.path.join(pdb_dir, f"{pdb_id}.pdb")
    if not os.path.exists(pdb_path):
        print(f"Structure {pdb_id}.pdb not found in local path '{pdb_dir}'. Attempting remote fetching from RCSB...")
        os.makedirs(pdb_dir, exist_ok=True)
        url = download_url.format(pdb_id=pdb_id)
        response = requests.get(url, timeout=timeout)
        if response.status_code == 200:
            with open(pdb_path, "w") as f:
                f.write(response.text)
            print(f"Successfully downloaded structural file {pdb_id}.pdb to {pdb_path}")
        else:
            raise FileNotFoundError(f"Could not locate or download structural file context for reference sequence: {pdb_id}")
    return pdb_path

def submit_to_hpepdock(receptor_pdb_path, peptide_seq, job_name, pdb_id, server_url, email,
                       max_retries, request_timeout, retry_delay, peptide_models):
    """Submits the structural screening job to the target HPEPDOCK server instance using robust connection context handling."""
    print(f"Submitting computational docking run: '{job_name}' (Target Core Sequence: {peptide_seq})...")

    data = {
        "jobname": job_name,
        "email": email,
        "pdbid1": pdb_id,
        "fastaseq2": peptide_seq,
        "docktyp": "dock",
        "nmodpep": str(peptide_models),
        "upload": "Submit"
    }

    with open(receptor_pdb_path, "rb") as f:
        files = {
            "pdbfile1": (os.path.basename(receptor_pdb_path), f, "application/octet-stream")
        }

        for attempt in range(1, max_retries + 1):
            try:
                response = requests.post(server_url, data=data, files=files, timeout=request_timeout)
                if response.status_code == 200:
                    print(f"Job context '{job_name}' accepted on target framework loop {attempt}.")
                    return response.text
                else:
                    print(f"Endpoint warning: received server status code {response.status_code} during retry cycle {attempt}.")
            except requests.exceptions.Timeout:
                print(f"Connection timeout occurred on attempt {attempt}. Retrying execution block in 15 seconds...")
                time.sleep(retry_delay)
            except Exception as e:
                print(f"Runtime execution error caught on attempt {attempt}: {e}")
                time.sleep(retry_delay)

        print(f"Critical Error: Execution failed across all {max_retries} connection loops for run context reference: '{job_name}'.")
        return None

def main():
    args = parse_args()
    config = load_config(args.config)
    args.csv = resolve_path(args.csv, config)
    args.output = resolve_path(args.output, config)
    args.pdb_dir = resolve_path(args.pdb_dir or config["paths"]["mhc_structures"], config)
    service = config["external_services"]["hpepdock"]
    rcsb_url = config["external_services"]["rcsb"]["download_url"]
    rcsb_timeout = config["external_services"]["rcsb"]["timeout_seconds"]
    args.url = args.url or config["docking"]["server_url"]
    args.email = args.email or config["docking"]["email"]
    args.retries = args.retries if args.retries is not None else config["docking"]["max_retries"]
    args.receptors = args.receptors or config["docking"]["receptors"]
    
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
    top_peptides = df[pep_col].unique()[:config["docking"]["max_epitopes"]]
    
    print(f"Restricted evaluation to the top {len(top_peptides)} unique candidate epitopes for structural alignment checks.")

    # Locate/Download the receptors using the relative directory specified
    receptor_files = {}
    for pdb_id in args.receptors:
        receptor_files[pdb_id] = verify_and_get_pdb(
            pdb_id, args.pdb_dir, rcsb_url, rcsb_timeout
        )

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
                max_retries=args.retries,
                request_timeout=config["docking"]["timeout_seconds"],
                retry_delay=service["retry_delay_seconds"],
                peptide_models=service["peptide_models"],
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

            time.sleep(service["rate_limit_delay_seconds"])

    print(f"\nPipeline sync completed successfully. Execution matrix state preserved in log: '{args.output}'")

if __name__ == "__main__":
    main()