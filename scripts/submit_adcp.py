#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess
import pandas as pd
import shutil

def parse_arguments():
    parser = argparse.ArgumentParser(description="Automated High-Throughput Batch Screening via ADCP")
    parser.add_argument("--csv", required=True, help="Path to filtered_mhc1_epitopes.csv")
    parser.add_argument("--output", required=True, help="Path to write the final adcp_screening_summary.txt")
    parser.add_argument("--bin-path", required=True, help="Path to ADFRsuite bin directory")
    parser.add_argument("--replicas", type=int, default=40, help="Number of independent GA runs (-N)")
    parser.add_argument("--steps", type=int, default=1000000, help="Number of search steps per run (-n)")
    parser.add_argument("--threads", type=int, default=4, help="Number of CPU cores to assign per run (-c)")
    return parser.parse_args()

def extract_scores_from_output(output_prefix):
    """
    Parses ADCP output format to extract the best binding energy (totalE)
    by checking the generated .pdb and .out replica files.
    """
    possible_files = [
        f"{output_prefix}_1.out",
        f"{output_prefix}_1.pdb",
        f"{output_prefix}.pdb"
    ]
    
    target_file = None
    for f_path in possible_files:
        if os.path.exists(f_path):
            target_file = f_path
            break
            
    if not target_file:
        dir_name = os.path.dirname(output_prefix)
        base_name = os.path.basename(output_prefix)
        if os.path.exists(dir_name):
            matched = [os.path.join(dir_name, f) for f in os.listdir(dir_name) if f.startswith(base_name) and (f.endswith(".out") or f.endswith(".pdb"))]
            if matched:
                target_file = matched[0]

    if not target_file:
        return "Missing Output File"
        
    try:
        best_score = None
        with open(target_file, 'r') as f:
            for line in f:
                if "totalE" in line:
                    parts = line.split()
                    for i, p in enumerate(parts):
                        if "totalE" in p:
                            try:
                                score_val = float(parts[i+1].strip())
                                if best_score is None or score_val < best_score:
                                    best_score = score_val
                            except (ValueError, IndexError):
                                continue
        
        if best_score is not None:
            return best_score
            
        with open(target_file, 'r') as f:
            for line in f:
                if "ENERGY" in line.upper():
                    parts = line.split()
                    for p in parts:
                        try:
                            val = float(p.strip())
                            if val < 50.0: 
                                if best_score is None or val < best_score:
                                    best_score = val
                        except ValueError:
                            continue
                            
        return best_score if best_score is not None else "Score Missing"
    except Exception as e:
        return f"Parse Error: {str(e)}"

def main():
    args = parse_arguments()
    
    # 1. Validate inputs
    if not os.path.exists(args.csv):
        print(f"[-] Error: Input epitope CSV not found at {args.csv}")
        sys.exit(1) 
        
    adcp_executable = os.path.join(args.bin_path, "adcp")
    if not os.path.exists(adcp_executable):
        print(f"[-] Error: adcp binary not found at path: {adcp_executable}")
        sys.exit(1)

    # 2. Load the pipeline-generated epitopes
    try:
        df = pd.read_csv(args.csv)
        seq_col = [col for col in df.columns if col.lower() in ['sequence', 'peptide', 'epitope']][0]
        epitopes = df[seq_col].dropna().unique().tolist()[:5]
        print(f"[+] Loaded {len(epitopes)} unique target epitopes for screening.")
    except Exception as e:
        print(f"[-] Error parsing CSV file columns: {e}")
        sys.exit(1)

    # 3. Discover generated grid targets in the data directory
    grid_dir = "data/MHC_I"
    if not os.path.exists(grid_dir):
        print(f"[-] Error: Grid data storage path does not exist: {grid_dir}")
        sys.exit(1)
        
    grid_files = [f for f in os.listdir(grid_dir) if f.endswith("_receptor.trg")]
    if not grid_files:
        print("[-] Error: No compiled receptor grids (*_receptor.trg) discovered.")
        sys.exit(1)
        
    print(f"[+] Discovered {len(grid_files)} target receptor structural grids.")

    # Prepare work space root directory
    screening_root = "results/docking/virtual_screening"
    os.makedirs(screening_root, exist_ok=True)
    
    summary_results = []

    # 4. Run the high-throughput loop matrix (All Epitopes x All Receptors)
    for grid_file in grid_files:
        receptor_id = grid_file.replace("_receptor.trg", "")
        src_trg_path = os.path.join(grid_dir, grid_file)
        
        # Localize target grid to run path to prevent structural parsing crashes
        local_trg_path = os.path.join(screening_root, grid_file)
        shutil.copy2(src_trg_path, local_trg_path)
        
        for pep in epitopes:
            output_name = f"screening_{receptor_id}_{pep}"
            output_prefix = os.path.join(screening_root, output_name)
            
            print(f"\n[Running] Target: {receptor_id} | Peptide: {pep}")
            
            cmd = [
                adcp_executable,
                "-t", grid_file,            # Pass filename directly due to working dir switch below
                "-s", pep,
                "-N", str(args.replicas),
                "-n", str(args.steps),
                "-c", str(args.threads),
                "-o", output_name
            ]
            
            try:
                # Execute from inside the localized directory to accommodate the tool's path expectations
                process = subprocess.run(
                    cmd, 
                    cwd=screening_root, 
                    stdout=subprocess.PIPE, 
                    stderr=subprocess.PIPE, 
                    text=True, 
                    check=True
                )
                
                # Extract structural score
                top_score = extract_scores_from_output(output_prefix)
                print(f"[Success] {receptor_id} + {pep} -> Top Binding Energy: {top_score} kcal/mol")
                summary_results.append({
                    "Receptor": receptor_id,
                    "Peptide": pep,
                    "Affinity_Score": top_score,
                    "Status": "Completed"
                })
                
            except subprocess.CalledProcessError as e:
                print(f"[-] Execution Failed for {receptor_id} + {pep}")
                print(f"Error Details: {e.stderr}")
                summary_results.append({
                    "Receptor": receptor_id,
                    "Peptide": pep,
                    "Affinity_Score": "N/A",
                    "Status": f"Failed: {str(e)}"
                })
                
        # Cleanup localized target layout files when target receptor batch finishes
        if os.path.exists(local_trg_path):
            os.remove(local_trg_path)

    # 5. Write the final master execution log requested by Snakemake's rule all
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    summary_df = pd.DataFrame(summary_results)
    
    with open(args.output, "w") as f:
        f.write("=========================================================\n")
        f.write("      AUTODOCK CRANKPEP VIRTUAL SCREENING SUMMARY        \n")
        f.write("=========================================================\n\n")
        f.write(summary_df.to_string(index=False))
        f.write("\n\n=== End of Pipeline Report ===")
        
    print(f"\n[+] Master screening operations complete. Summary logs stored at: {args.output}")

if __name__ == "__main__":
    main()