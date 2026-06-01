# =============================================================================
# Snakefile — CPV2 Canine Vaccine Design Pipeline v0.1.0 
# =============================================================================

configfile: "config.yaml"

# 1. Extract the chosen tool dynamically from config
CHOSEN_TOOL = config.get("docking", {}).get("tool", "hpepdock").lower()

# 2. Define the input router function
def get_final_docking_target(wildcards):
    if CHOSEN_TOOL == "adcp":
        return "results/docking/adcp_screening_summary.txt"
    else:
        return "results/docking/hpepdock_submissions.txt"

rule all:
    input:
        "results/filtering/host_homology_exclusion.csv",
        "data/align/VP2_conservation.csv",
        "results/safety/mhc1_algpred2_checked.csv",
        "results/safety/mhc2_algpred2_checked.csv",
        "results/safety/bcell_algpred2_checked.csv",
        "results/final_ranked_tcell_epitopes.csv",
        "results/consensus/overlapping_consensus_regions.csv",
        
        # CHANGED: Replaced the hardcoded path with dynamic function
        get_final_docking_target,  
        
        "results/vaccine_construct/final_vaccine_sequence.fasta"


# =============================================================================
# MODULE 1: DATA PREPARATION
# =============================================================================

rule extract_ref_vp2:
    input:
        ref = config["reference"]["genome"]
    output:
        "data/fasta/FJ197845_VP2.fasta"
    params:
        ref_start = config["reference"]["vp2_gene"]["start"],
        ref_end   = config["reference"]["vp2_gene"]["end"]
    conda:
        "envs/blast.yaml"
    shell:
        "python scripts/extract_ref_vp2.py {input.ref} {params.ref_start} {params.ref_end} {output}"

rule extract_query_vp2:
    input:
        genomes  = expand("data/fasta/{sample}.fasta", sample=config["samples"]),
        ref_done = "data/fasta/FJ197845_VP2.fasta"
    output:
        status = expand(
            "data/fasta/{sample}_VP2.fasta",
            sample=[s for s in config["samples"] if s != "FJ197845"]
        )
    conda:
        "envs/blast.yaml"
    shell:
        "python scripts/vp2_extract.py"

rule sanitize_sequences:
    input:
        "data/fasta/{sample}_VP2.fasta"
    output:
        "data/fasta_cleaned/{sample}_VP2_qc.fasta"
    conda:
        "envs/alignment.yaml"
    shell:
        "python scripts/sanitize_sequences.py {input} {output}"

rule build_alignment_pool:
    input:
        expand("data/fasta_cleaned/{sample}_VP2_qc.fasta", sample=config["samples"])
    output:
        "data/protein/master_pool.fasta"
    conda:
        "envs/alignment.yaml"
    shell:
        "cat {input} > {output}"


# =============================================================================
# MODULE 2: HOST FILTERING (Clean Peptide Universe Generation)
# =============================================================================

rule make_canine_blastdb:
    input:
        fasta = "data/reference/dog_proteome.fa"
    output:
        phr = "data/reference/dog_db.phr",
        pin = "data/reference/dog_db.pin",
        psq = "data/reference/dog_db.psq"
    conda:
        "envs/blast.yaml"
    shell:
        "makeblastdb -in {input.fasta} -dbtype prot -out data/reference/dog_db"

rule blastp:
    input:
        query  = "data/fasta_cleaned/{sample}_VP2_qc.fasta",
        db_phr = "data/reference/dog_db.phr",
        db_pin = "data/reference/dog_db.pin",
        db_psq = "data/reference/dog_db.psq"
    output:
        txt = "data/blast/{sample}_blastp.txt"
    params:
        db_base = "data/reference/dog_db"
    conda:
        "envs/blast.yaml"
    shell:
        "python scripts/vp2_blastp.py {input.query} {params.db_base} {output.txt}"

rule filter_host_homology:
    input:
        blast_txt = expand("data/blast/{sample}_blastp.txt", sample=config["samples"])
    output:
        csv         = "results/filtering/host_homology_exclusion.csv",
        json_report = "results/host_exclusion_report.json"
    conda:
        "envs/blast.yaml"
    shell:
        """
        python scripts/filter_host_homology.py \
            --inputs {input.blast_txt} \
            --output_csv {output.csv} \
            --output_json {output.json_report}
        """


# =============================================================================
# MODULE 3: CONSERVATION ANALYSIS (Orthogonal Property Channel)
# =============================================================================

rule mafft:
    input:
        "data/protein/master_pool.fasta"
    output:
        "data/align/All_VP2_aligned.fasta"
    params:
        threads = config["alignment"].get("threads", 4)
    conda:
        "envs/alignment.yaml"
    shell:
        "mafft --maxiterate 1000 --globalpair --thread {params.threads} {input} > {output}"

rule conservation:
    input:
        aligned = "data/align/All_VP2_aligned.fasta"
    output:
        csv = "data/align/VP2_conservation.csv"
    conda:
        "envs/alignment.yaml"
    shell:
        "python scripts/calculate_conservation.py {input.aligned} {output.csv}"


# =============================================================================
# MODULE 4: EPITOPE DISCOVERY (Parallel Predictors From Clean Universe)
# =============================================================================

rule netmhcpan:
    input:
        fasta = "data/protein/master_pool.fasta",
        universe = "results/filtering/host_homology_exclusion.csv"
    output:
        csv = "data/align/vp2_strong_binders.csv"
    params:
        alleles = ",".join(config["tcell_epitope"]["dla_alleles"]),
        lengths = ",".join(str(l) for l in config["tcell_epitope"]["peptide_lengths"]),
        rank    = config["tcell_epitope"]["thresholds"]["strong_binder_rank"]
    conda:
        "envs/netmhcpan.yaml"
    shell:
        """
        python scripts/netmhcpan_runner.py \
            {input.fasta} \
            '{params.alleles}' \
            '{params.lengths}' \
            {params.rank} \
            {output.csv}
        """

rule iedb_mhc2:
    input:
        fasta = "data/protein/master_pool.fasta",
        universe = "results/filtering/host_homology_exclusion.csv"
    output:
        csv = "results/immunogenicity/iedb_mhc2_predictions.csv"
    params:
        method  = config["mhc2_epitope"].get("method", "nn_align"),
        alleles = ",".join(config["mhc2_epitope"]["dla_mhc2_alleles"]),
        lengths = ",".join(str(l) for l in config["mhc2_epitope"].get("lengths", [15]))
    conda:
        "envs/iedb_tools.yaml"
    shell:
        """
        python scripts/run_iedb_mhc2.py \
            --input {input.fasta} \
            --method {params.method} \
            --alleles '{params.alleles}' \
            --lengths '{params.lengths}' \
            --output {output.csv}
        """

rule bepipred3:
    input:
        fasta = "data/protein/master_pool.fasta",
        universe = "results/filtering/host_homology_exclusion.csv"
    output:
        preds        = "results/bcell/Bcell_epitope_preds.fasta",
        top_pct_pred = "results/bcell/Bcell_epitope_top_20pct_preds.fasta",
        lin_pred     = "results/bcell/Bcell_linepitope_top_20pct_preds.fasta",
        raw_csv      = "results/bcell/raw_output.csv"
    params:
        script      = config["bcell_epitope"]["script_path"],
        out_dir     = "results/bcell",
        pred_mode   = config["bcell_epitope"]["prediction_mode"],
        threshold   = config["bcell_epitope"]["threshold"],
        window_size = config["bcell_epitope"]["rolling_window_size"]
    conda:
        "envs/bepipred3.yaml"
    shell:
        """
        python {params.script} \
            -i {input.fasta} \
            -o {params.out_dir} \
            -pred {params.pred_mode} \
            -t {params.threshold} \
            -rolling_window_size {params.window_size}
        """

rule parse_bcell_coords:
    input:
        fasta = "results/bcell/Bcell_linepitope_top_20pct_preds.fasta",
        csv   = "results/bcell/raw_output.csv"
    output:
        summary = "results/bcell/bcell_linear_epitopes_summary.csv"
    conda:
        "envs/canine_immuno.yaml"
    shell:
        "python scripts/extract_bcell_coordinates.py {input.fasta} {input.csv} {output.summary}"


# =============================================================================
# MODULE 5: SAFETY SCREENING (Parallel Chains for Independent Gatekeeping)
# =============================================================================

rule toxinpred3_mhc1:
    input:  csv = "data/align/vp2_strong_binders.csv"
    output: csv = "results/safety/mhc1_toxinpred3_checked.csv"
    conda:  "envs/toxinpred3.yaml"
    shell:  "python scripts/run_toxinpred3.py --input {input.csv} --output {output.csv}"

rule algpred2_mhc1:
    input:  csv = "results/safety/mhc1_toxinpred3_checked.csv"
    output: csv = "results/safety/mhc1_algpred2_checked.csv"
    conda:  "envs/algpred2.yaml"
    shell:  "python scripts/run_algpred2.py --input {input.csv} --output {output.csv}"

rule toxinpred3_mhc2:
    input:  csv = "results/immunogenicity/iedb_mhc2_predictions.csv"
    output: csv = "results/safety/mhc2_toxinpred3_checked.csv"
    conda:  "envs/toxinpred3.yaml"
    shell:  "python scripts/run_toxinpred3.py --input {input.csv} --output {output.csv}"

rule algpred2_mhc2:
    input:  csv = "results/safety/mhc2_toxinpred3_checked.csv"
    output: csv = "results/safety/mhc2_algpred2_checked.csv"
    conda:  "envs/algpred2.yaml"
    shell:  "python scripts/run_algpred2.py --input {input.csv} --output {output.csv}"

rule toxinpred3_bcell:
    input:  csv = "results/bcell/bcell_linear_epitopes_summary.csv"
    output: csv = "results/safety/bcell_toxinpred3_checked.csv"
    conda:  "envs/toxinpred3.yaml"
    shell:  "python scripts/run_toxinpred3.py --input {input.csv} --output {output.csv}"

rule algpred2_bcell:
    input:  csv = "results/safety/bcell_toxinpred3_checked.csv"
    output: csv = "results/safety/bcell_algpred2_checked.csv"
    conda:  "envs/algpred2.yaml"
    shell:  "python scripts/run_algpred2.py --input {input.csv} --output {output.csv}"


# =============================================================================
# MODULE 6: FUNCTIONAL SCORING (No Filtering - Parallel Stream Aggregation)
# =============================================================================

rule iedb_mhc1_processing:
    input:   peptides = "data/align/vp2_strong_binders.csv"
    output:  csv = "results/immunogenicity/iedb_mhc1_processing.csv"
    params:
        method     = config["mhc1_processing"].get("method", "recommended"),
        proteasome = config["mhc1_processing"].get("proteasome", "immuno"),
        alleles    = ",".join(config["tcell_epitope"]["dla_alleles"])
    conda:
        "envs/iedb_tools.yaml"
    shell:
        """
        python scripts/run_iedb_processing.py \
            --input {input.peptides} \
            --method {params.method} \
            --proteasome {params.proteasome} \
            --alleles '{params.alleles}' \
            --output {output.csv}
        """

rule deeptmhmm:
    input:   fasta = "data/protein/master_pool.fasta"
    output:  csv = "results/topology/vp2_tm_topology.csv"
    conda:   "envs/deeptmhmm.yaml"
    shell:   "python scripts/deeptmhmm_batch.py {input.fasta} {output.csv}"

rule ensemble_immunogenicity:
    input:
        mhc2  = "results/immunogenicity/iedb_mhc2_predictions.csv",
        proc  = "results/immunogenicity/iedb_mhc1_processing.csv",
        breed = "config/breed_weights.csv"
    output:
        ensemble = "results/immunogenicity/ensemble_immunogenicity.csv"
    params:
        w_mhc2 = config["ensemble"]["weight_mhc2"],
        w_proc = config["ensemble"]["weight_processing"],
        w_pop  = config["ensemble"]["weight_population_coverage"]
    conda:
        "envs/iedb_tools.yaml"
    shell:
        """
        python scripts/ensemble_immunogenicity_v2.py \
            --mhc2-scores {input.mhc2} \
            --processing {input.proc} \
            --breed-weights {input.breed} \
            --output {output.ensemble} \
            --weight-mhc2 {params.w_mhc2} \
            --weight-processing {params.w_proc} \
            --weight-population {params.w_pop}
        """


# =============================================================================
# MODULE 7: FINAL RANKING & PRIORITIZATION
# =============================================================================

rule canine_epitope_rank:
    input:
        mhc_csv      = "results/safety/mhc1_algpred2_checked.csv",
        tm_csv       = "results/topology/vp2_tm_topology.csv",
        cons_csv     = "data/align/VP2_conservation.csv",
        blast_filter = "results/filtering/host_homology_exclusion.csv",
        ensemble_csv = "results/immunogenicity/ensemble_immunogenicity.csv",
        tox_csv      = "results/safety/mhc1_toxinpred3_checked.csv",  
        allergen_csv = "results/safety/mhc1_algpred2_checked.csv"    
    output:
        tcell_out    = "results/final_ranked_tcell_epitopes.csv"
    conda:
        "envs/canine_immuno.yaml"
    shell:
        """
        python scripts/final_weighted_prioritization.py \
            --mhc           {input.mhc_csv} \
            --topology      {input.tm_csv} \
            --conservation  {input.cons_csv} \
            --blast         {input.blast_filter} \
            --immunogenicity {input.ensemble_csv} \
            --toxicity      {input.tox_csv} \
            --allergenicity {input.allergen_csv} \
            --config        config.yaml \
            --output        {output.tcell_out}
        """


# =============================================================================
# MODULE 8: INTERSECTION & DOWNSTREAM VALIDATION
# =============================================================================

rule clean_and_intersect_epitopes:
    input:
        mhc1_raw    = "results/final_ranked_tcell_epitopes.csv",
        mhc2_raw    = "results/safety/mhc2_algpred2_checked.csv",
        bcell_raw   = "results/safety/bcell_algpred2_checked.csv",
        cons_csv    = "data/align/VP2_conservation.csv",
        blast_csv   = "results/filtering/host_homology_exclusion.csv",
        mhc1_tox    = "results/safety/mhc1_toxinpred3_checked.csv",   
        mhc1_alg    = "results/safety/mhc1_algpred2_checked.csv"      
    output:
        mhc1_clean  = "results/consensus/filtered_mhc1_epitopes.csv",
        mhc2_clean  = "results/consensus/filtered_mhc2_epitopes.csv",
        bcell_clean = "results/consensus/filtered_bcell_epitopes.csv",
        overlapping = "results/consensus/overlapping_consensus_regions.csv"
    conda:
        "envs/canine_immuno.yaml"
    shell:
        """
        python scripts/find_consensus_epitopes.py \
            --mhc1 {input.mhc1_raw} \
            --mhc2 {input.mhc2_raw} \
            --bcell {input.bcell_raw} \
            --cons {input.cons_csv} \
            --blast {input.blast_csv} \
            --tox {input.mhc1_tox} \
            --alg {input.mhc1_alg} \
            --out_dir results/consensus
        """

# -----------------------------------------------------------------------------
# SWITCH OPTION A: REMOTE HPEPDOCK STRUCTURAL DOCKING ENGINE
# -----------------------------------------------------------------------------
rule hpepdock_screening:
    input:
        csv = "results/consensus/filtered_mhc1_epitopes.csv",
        pdb_dir = config["paths"]["mhc_structures"]
    output:
        log = "results/docking/hpepdock_submissions.txt"
    params:
        receptors = config.get("docking", {}).get("receptors", []),
        url       = config.get("docking", {}).get("server_url", ""),
        email     = config.get("docking", {}).get("email", ""),
        retries   = config.get("docking", {}).get("max_retries", 3)
    conda:
        "envs/iedb_tools.yaml"
    shell:
        """
        python scripts/submit_hpepdock.py \
            --csv {input.csv} \
            --pdb-dir {input.pdb_dir} \
            --output {output.log} \
            --url '{params.url}' \
            --email '{params.email}' \
            --retries {params.retries} \
            --receptors "{params.receptors}"
        """

# =============================================================================
# MODULE 8.1: LOCAL ADCP AUTODOCK CRANKPEP STRUCTURAL WORKFLOW
# =============================================================================

rule adcp_extract_native_peptide:
    input:
        pdb = "data/MHC_I/{receptor}.pdb"
    output:
        peptide_pdbqt = "data/MHC_I/{receptor}_native_peptide.pdbqt"
    params:
        bin_path = config.get("docking", {}).get("adfr_suite_path", "")
    conda:
        "envs/iedb_tools.yaml"
    shell:
        """
        cd data/MHC_I
        
        awk 'BEGIN {{FS=""}} /^ATOM/ && substr($0, 22, 1)=="C" {{print $0}}' {wildcards.receptor}.pdb > {wildcards.receptor}_temp_pep.pdb
        
        
        {params.bin_path}/prepare_ligand \
            -l {wildcards.receptor}_temp_pep.pdb \
            -o {wildcards.receptor}_native_peptide.pdbqt \
            -A hydrogens

        rm {wildcards.receptor}_temp_pep.pdb
        """

rule adcp_prepare_receptor:
    input:
        pdb = "data/MHC_I/{receptor}.pdb"
    output:
        pdbqt = "data/MHC_I/{receptor}_receptor.pdbqt"
    params:
        bin_path = config.get("docking", {}).get("adfr_suite_path", "")
    conda:
        "envs/iedb_tools.yaml"
    shell:
        """
        awk 'BEGIN {{FS=""}} /^ATOM/ && (substr($0, 22, 1)=="A" || substr($0, 22, 1)=="B") {{print $0}}' {input.pdb} > data/MHC_I/{wildcards.receptor}_clean_target.pdb
        {params.bin_path}/prepare_receptor \
            -r data/MHC_I/{wildcards.receptor}_clean_target.pdb \
            -o {output.pdbqt} \
            -A hydrogens \
            -U nsw
        rm data/MHC_I/{wildcards.receptor}_clean_target.pdb
        """

rule adcp_generate_affinity_grid:
    input:
        receptor_pdbqt = "data/MHC_I/{receptor}_receptor.pdbqt",
        peptide_pdbqt  = "data/MHC_I/{receptor}_native_peptide.pdbqt"
    output:
        trg = "data/MHC_I/{receptor}_receptor.trg"
    params:
        bin_path = config.get("docking", {}).get("adfr_suite_path", ""),
        padding  = config.get("docking", {}).get("grid_padding", 4.0)
    conda:
        "envs/iedb_tools.yaml"
    shell:
        """
        {params.bin_path}/agfr \
            -r {input.receptor_pdbqt} \
            -l {input.peptide_pdbqt} \
            -b ligand \
            -P {params.padding} \
            -o data/MHC_I/{wildcards.receptor}_receptor
        """

rule adcp_screening:
    input:
        csv = "results/consensus/filtered_mhc1_epitopes.csv",
        # Ensures all receptors in your config are fully mapped and ready before docking
        trgs = expand("data/MHC_I/{receptor}_receptor.trg", receptor=config.get("docking", {}).get("receptors", []))
    output:
        summary = "results/docking/adcp_screening_summary.txt"
    params:
        bin_path = config.get("docking", {}).get("adfr_suite_path", ""),
        replicas = config.get("docking", {}).get("replicas", 40),
        steps    = config.get("docking", {}).get("steps", 1000000),
        threads  = config.get("docking", {}).get("threads", 4)
    conda:
        "envs/mhc_vax.yaml"
    shell:
        """
        python scripts/submit_adcp.py \
            --csv {input.csv} \
            --output {output.summary} \
            --bin-path {params.bin_path} \
            --replicas {params.replicas} \
            --steps {params.steps} \
            --threads {params.threads}
        """
# =============================================================================
# MODULE 9: VACCINE ASSEMBLY CHIMERA
# =============================================================================

rule construct_vaccine:
    input:
        overlapping = "results/consensus/overlapping_consensus_regions.csv",
        mhc1_clean  = "results/consensus/filtered_mhc1_epitopes.csv",
        mhc2_clean  = "results/consensus/filtered_mhc2_epitopes.csv",
        bcell_clean = "results/consensus/filtered_bcell_epitopes.csv"
    output:
        fasta = "results/vaccine_construct/final_vaccine_sequence.fasta"
    conda:
        "envs/canine_immuno.yaml"
    shell:
        """
        python scripts/assemble_vaccine_chimera.py \
            --intersected {input.overlapping} \
            --mhc1 {input.mhc1_clean} \
            --mhc2 {input.mhc2_clean} \
            --bcell {input.bcell_clean} \
            --output {output.fasta}
        """