# =============================================================================
# Snakefile — CPV2 Canine Vaccine Design Pipeline v0.3.0
# =============================================================================

configfile: "config.yaml"
from pathlib import Path

PATHS = config["paths"]
TOOLS = config["tools"]
RAW_FASTA = PATHS["raw_genomes"]
CLEAN_FASTA = PATHS["cleaned_genomes"]
PROTEIN_DIR = PATHS["proteins"]
ALIGN_DIR = PATHS["alignment"]
BLAST_DIR = PATHS["blast"]
FILTER_DIR = PATHS["filtering"]
IMMUNO_DIR = PATHS["immunogenicity"]
BCELL_DIR = PATHS["bcell"]
SAFETY_DIR = PATHS["safety"]
TOPOLOGY_DIR = PATHS["topology"]
CONSENSUS_DIR = PATHS["consensus"]
VACCINE_DIR = PATHS["vaccine_construct"]
MHC_DIR = PATHS["mhc_structures"]
DOG_DB = config["databases"]["dog_blast_db"]
REFERENCE_ID = Path(config["reference"]["genome"]).stem

# 1. Extract the chosen tool dynamically from config
CHOSEN_TOOL = config.get("docking", {}).get("tool", "hpepdock").lower()

# 2. Define the input router function
def get_final_docking_target(wildcards):
    if CHOSEN_TOOL == "adcp":
        return f"{PATHS['docking']}/adcp_screening_summary.txt"
    else:
        return f"{PATHS['docking']}/hpepdock_submissions.txt"

rule all:
    input:
        PATHS["host_exclusion_csv"],
        PATHS["conservation_csv"],
        f"{SAFETY_DIR}/mhc1_algpred2_checked.csv",
        f"{SAFETY_DIR}/mhc2_algpred2_checked.csv",
        f"{SAFETY_DIR}/bcell_algpred2_checked.csv",
        PATHS["final_ranked_tcell_csv"],
        f"{CONSENSUS_DIR}/overlapping_consensus_regions.csv",
        
        # CHANGED: Replaced the hardcoded path with dynamic function
        get_final_docking_target,  
        
        f"{VACCINE_DIR}/final_vaccine_sequence.fasta"


# =============================================================================
# MODULE 1: DATA PREPARATION
# =============================================================================

rule extract_ref_vp2:
    input:
        ref = config["reference"]["genome"]
    output:
        PATHS["reference_vp2_fasta"]
    params:
        ref_start = config["reference"]["vp2_gene"]["start"],
        ref_end   = config["reference"]["vp2_gene"]["end"]
    conda:
        "envs/blast.yaml"
    shell:
        "python scripts/extract_ref_vp2.py {input.ref} {params.ref_start} {params.ref_end} {output} --config config.yaml"

rule extract_query_vp2:
    input:
        genomes  = expand(f"{RAW_FASTA}/{{sample}}.fasta", sample=config["samples"]),
        ref_done =         PATHS["reference_vp2_fasta"]
    output:
        status = expand(
            f"{RAW_FASTA}/{{sample}}_VP2.fasta",
            sample=[s for s in config["samples"] if s != REFERENCE_ID]
        )
    conda:
        "envs/blast.yaml"
    shell:
        "python scripts/vp2_extract.py --config config.yaml"

rule sanitize_sequences:
    input:
        f"{RAW_FASTA}/{{sample}}_VP2.fasta"
    output:
        f"{CLEAN_FASTA}/{{sample}}_VP2_qc.fasta"
    conda:
        "envs/alignment.yaml"
    shell:
        "python scripts/sanitize_sequences.py {input} {output} --config config.yaml"

rule build_alignment_pool:
    input:
        expand(f"{CLEAN_FASTA}/{{sample}}_VP2_qc.fasta", sample=config["samples"])
    output:
        PATHS["master_protein_fasta"]
    conda:
        "envs/alignment.yaml"
    shell:
        "cat {input} > {output}"


# =============================================================================
# MODULE 2: HOST FILTERING (Clean Peptide Universe Generation)
# =============================================================================

rule make_canine_blastdb:
    input:
        fasta = config["databases"]["dog_proteome_fasta"]
    output:
        phr = f"{DOG_DB}.phr",
        pin = f"{DOG_DB}.pin",
        psq = f"{DOG_DB}.psq"
    conda:
        "envs/blast.yaml"
    params:
        db_base = DOG_DB,
        makeblastdb = TOOLS["makeblastdb"].get("path") or TOOLS["makeblastdb"]["command"]
    shell:
        "{params.makeblastdb} -in {input.fasta} -dbtype prot -out {params.db_base}"

rule blastp:
    input:
        query  = f"{CLEAN_FASTA}/{{sample}}_VP2_qc.fasta",
        db_phr = f"{DOG_DB}.phr",
        db_pin = f"{DOG_DB}.pin",
        db_psq = f"{DOG_DB}.psq"
    output:
        txt = f"{BLAST_DIR}/{{sample}}_blastp.txt"
    params:
        db_base = DOG_DB,
        config_path = "config.yaml"
    conda:
        "envs/blast.yaml"
    shell:
        "python scripts/vp2_blastp.py {input.query} {params.db_base} {output.txt} --config {params.config_path}"

rule filter_host_homology:
    input:
        blast_txt = expand(f"{BLAST_DIR}/{{sample}}_blastp.txt", sample=config["samples"])
    output:
        csv         = PATHS["host_exclusion_csv"],
        json_report = PATHS["host_exclusion_report"]
    conda:
        "envs/blast.yaml"
    shell:
        """
        python scripts/filter_host_homology.py \
            --inputs {input.blast_txt} \
            --output_csv {output.csv} \
            --output_json {output.json_report} \
            --config config.yaml
        """


# =============================================================================
# MODULE 3: CONSERVATION ANALYSIS (Orthogonal Property Channel)
# =============================================================================

rule mafft:
    input:
        PATHS["master_protein_fasta"]
    output:
        PATHS["aligned_proteins"]
    params:
        threads = config["alignment"]["threads"],
        mode = config["alignment"]["mode"],
        maxiterate = config["alignment"]["maxiterate"],
        gap_penalty = config["alignment"]["gap_penalty"],
        mafft = TOOLS["mafft"].get("path") or TOOLS["mafft"]["command"]
    conda:
        "envs/alignment.yaml"
    shell:
        "{params.mafft} {params.mode} --maxiterate {params.maxiterate} --op {params.gap_penalty} --thread {params.threads} {input} > {output}"

rule conservation:
    input:
        aligned = PATHS["aligned_proteins"]
    output:
        csv = PATHS["conservation_csv"]
    conda:
        "envs/alignment.yaml"
    shell:
        "python scripts/calculate_conservation.py {input.aligned} {output.csv} --config config.yaml"


# =============================================================================
# MODULE 4: EPITOPE DISCOVERY (Parallel Predictors From Clean Universe)
# =============================================================================

rule netmhcpan:
    input:
        fasta = PATHS["master_protein_fasta"],
        universe = PATHS["host_exclusion_csv"]
    output:
        csv = PATHS["mhc1_predictions_csv"]
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
            {output.csv} \
            --config config.yaml
        """

rule iedb_mhc2:
    input:
        fasta = PATHS["master_protein_fasta"],
        universe = PATHS["host_exclusion_csv"]
    output:
        csv = f"{IMMUNO_DIR}/iedb_mhc2_predictions.csv"
    params:
        method  = config["mhc2_epitope"].get("method", "nn_align"),
        alleles = ",".join(config["mhc2_epitope"]["hla_surrogate_alleles"]),
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
            --output {output.csv} \
            --config config.yaml
        """

rule bepipred3:
    input:
        fasta = PATHS["master_protein_fasta"],
        universe = PATHS["host_exclusion_csv"]
    output:
        preds        = f"{BCELL_DIR}/Bcell_epitope_preds.fasta",
        top_pct_pred = f"{BCELL_DIR}/Bcell_epitope_top_20pct_preds.fasta",
        lin_pred     = f"{BCELL_DIR}/Bcell_linepitope_top_20pct_preds.fasta",
        raw_csv      = f"{BCELL_DIR}/raw_output.csv"
    params:
        out_dir     = BCELL_DIR,
        pred_mode   = config["bcell_epitope"]["prediction_mode"],
        threshold   = config["bcell_epitope"]["threshold"],
        window_size = config["bcell_epitope"]["rolling_window_size"]
    conda:
        "envs/bepipred3.yaml"
    shell:
        """
        python scripts/run_bepipred3.py \
            --input {input.fasta} \
            --output-dir {params.out_dir} \
            --prediction-mode {params.pred_mode} \
            --threshold {params.threshold} \
            --window-size {params.window_size} \
            --config config.yaml
        """

rule parse_bcell_coords:
    input:
        fasta = f"{BCELL_DIR}/Bcell_linepitope_top_20pct_preds.fasta",
        csv   = f"{BCELL_DIR}/raw_output.csv"
    output:
        summary = f"{BCELL_DIR}/bcell_linear_epitopes_summary.csv"
    conda:
        "envs/canine_immuno.yaml"
    shell:
        "python scripts/extract_bcell_coordinates.py {input.fasta} {input.csv} {output.summary} --config config.yaml"


# =============================================================================
# MODULE 5: SAFETY SCREENING (Parallel Chains for Independent Gatekeeping)
# =============================================================================

rule toxinpred3_mhc1:
    input:  csv = PATHS["mhc1_predictions_csv"]
    output: csv = f"{SAFETY_DIR}/mhc1_toxinpred3_checked.csv"
    conda:  "envs/toxinpred3.yaml"
    shell:  "python scripts/run_toxinpred3.py --input {input.csv} --output {output.csv} --config config.yaml"

rule algpred2_mhc1:
    input:  csv = f"{SAFETY_DIR}/mhc1_toxinpred3_checked.csv"
    output: csv = f"{SAFETY_DIR}/mhc1_algpred2_checked.csv"
    conda:  "envs/algpred2.yaml"
    shell:  "python scripts/run_algpred2.py --input {input.csv} --output {output.csv} --config config.yaml"

rule toxinpred3_mhc2:
    input:  csv = f"{IMMUNO_DIR}/iedb_mhc2_predictions.csv"
    output: csv = f"{SAFETY_DIR}/mhc2_toxinpred3_checked.csv"
    conda:  "envs/toxinpred3.yaml"
    shell:  "python scripts/run_toxinpred3.py --input {input.csv} --output {output.csv} --config config.yaml"

rule algpred2_mhc2:
    input:  csv = f"{SAFETY_DIR}/mhc2_toxinpred3_checked.csv"
    output: csv = f"{SAFETY_DIR}/mhc2_algpred2_checked.csv"
    conda:  "envs/algpred2.yaml"
    shell:  "python scripts/run_algpred2.py --input {input.csv} --output {output.csv} --config config.yaml"

rule toxinpred3_bcell:
    input:  csv = f"{BCELL_DIR}/bcell_linear_epitopes_summary.csv"
    output: csv = f"{SAFETY_DIR}/bcell_toxinpred3_checked.csv"
    conda:  "envs/toxinpred3.yaml"
    shell:  "python scripts/run_toxinpred3.py --input {input.csv} --output {output.csv} --config config.yaml"

rule algpred2_bcell:
    input:  csv = f"{SAFETY_DIR}/bcell_toxinpred3_checked.csv"
    output: csv = f"{SAFETY_DIR}/bcell_algpred2_checked.csv"
    conda:  "envs/algpred2.yaml"
    shell:  "python scripts/run_algpred2.py --input {input.csv} --output {output.csv} --config config.yaml"


# =============================================================================
# MODULE 6: FUNCTIONAL SCORING (No Filtering - Parallel Stream Aggregation)
# =============================================================================

rule iedb_mhc1_processing:
    input:   peptides = PATHS["mhc1_predictions_csv"]
    output:  csv = f"{IMMUNO_DIR}/iedb_mhc1_processing.csv"
    params:
        method     = config["mhc1_processing"].get("method", "recommended"),
        proteasome = config["mhc1_processing"].get("proteasome", "immuno"),
        alleles    = ",".join(config["mhc1_processing"]["hla_surrogate_alleles"])
    conda:
        "envs/iedb_tools.yaml"
    shell:
        """
        python scripts/run_iedb_processing.py \
            --input {input.peptides} \
            --method {params.method} \
            --proteasome {params.proteasome} \
            --alleles '{params.alleles}' \
            --output {output.csv} \
            --config config.yaml
        """

rule deeptmhmm:
    input:   fasta = PATHS["master_protein_fasta"]
    output:  csv = PATHS["topology_csv"]
    conda:   "envs/deeptmhmm.yaml"
    shell:   "python scripts/deeptmhmm_batch.py {input.fasta} {output.csv} --config config.yaml"

rule ensemble_immunogenicity:
    input:
        mhc2  = f"{IMMUNO_DIR}/iedb_mhc2_predictions.csv",
        proc  = f"{IMMUNO_DIR}/iedb_mhc1_processing.csv",
        breed = config["ensemble"]["breed_weights_file"]
    output:
        ensemble = f"{IMMUNO_DIR}/ensemble_immunogenicity.csv"
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
            --weight-population {params.w_pop} \
            --config config.yaml
        """


# =============================================================================
# MODULE 7: FINAL RANKING & PRIORITIZATION
# =============================================================================

rule canine_epitope_rank:
    input:
        mhc_csv      = f"{SAFETY_DIR}/mhc1_algpred2_checked.csv",
        tm_csv       = PATHS["topology_csv"],
        cons_csv     = PATHS["conservation_csv"],
        blast_filter = PATHS["host_exclusion_csv"],
        ensemble_csv = f"{IMMUNO_DIR}/ensemble_immunogenicity.csv",
        tox_csv      = f"{SAFETY_DIR}/mhc1_toxinpred3_checked.csv",
        allergen_csv = f"{SAFETY_DIR}/mhc1_algpred2_checked.csv"
    output:
        tcell_out    = PATHS["final_ranked_tcell_csv"]
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
        mhc1_raw    = PATHS["final_ranked_tcell_csv"],
        mhc2_raw    = f"{SAFETY_DIR}/mhc2_algpred2_checked.csv",
        bcell_raw   = f"{SAFETY_DIR}/bcell_algpred2_checked.csv",
        cons_csv    = PATHS["conservation_csv"],
        blast_csv   = PATHS["host_exclusion_csv"],
        mhc1_tox    = f"{SAFETY_DIR}/mhc1_toxinpred3_checked.csv",
        mhc1_alg    = f"{SAFETY_DIR}/mhc1_algpred2_checked.csv"
    output:
        mhc1_clean  = f"{CONSENSUS_DIR}/filtered_mhc1_epitopes.csv",
        mhc2_clean  = f"{CONSENSUS_DIR}/filtered_mhc2_epitopes.csv",
        bcell_clean = f"{CONSENSUS_DIR}/filtered_bcell_epitopes.csv",
        overlapping = f"{CONSENSUS_DIR}/overlapping_consensus_regions.csv"
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
            --out_dir {CONSENSUS_DIR} \
            --config config.yaml
        """

# -----------------------------------------------------------------------------
# SWITCH OPTION A: REMOTE HPEPDOCK STRUCTURAL DOCKING ENGINE
# -----------------------------------------------------------------------------
rule hpepdock_screening:
    input:
        csv = f"{CONSENSUS_DIR}/filtered_mhc1_epitopes.csv",
        pdb_dir = config["paths"]["mhc_structures"]
    output:
        log = f"{PATHS['docking']}/hpepdock_submissions.txt"
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
            --receptors {params.receptors} \
            --config config.yaml
        """

# =============================================================================
# MODULE 8.1: LOCAL ADCP AUTODOCK CRANKPEP STRUCTURAL WORKFLOW
# =============================================================================

rule adcp_extract_native_peptide:
    input:
        pdb = f"{PATHS['mhc_structures']}/{{receptor}}.pdb"
    output:
        peptide_pdbqt = f"{PATHS['mhc_structures']}/{{receptor}}_native_peptide.pdbqt"
    conda:
        "envs/iedb_tools.yaml"
    shell:
        """
        awk 'BEGIN {{FS=""}} /^ATOM/ && substr($0, 22, 1)=="C" {{print $0}}' {input.pdb} > {MHC_DIR}/{wildcards.receptor}_temp_pep.pdb
        python scripts/run_adfrsuite.py \
            --tool prepare_ligand \
            --output {output.peptide_pdbqt} \
            --config config.yaml \
            -- \
            -l {MHC_DIR}/{wildcards.receptor}_temp_pep.pdb \
            -o {output.peptide_pdbqt} \
            -A hydrogens

        rm {MHC_DIR}/{wildcards.receptor}_temp_pep.pdb
        """

rule adcp_prepare_receptor:
    input:
        pdb = f"{PATHS['mhc_structures']}/{{receptor}}.pdb"
    output:
        pdbqt = f"{PATHS['mhc_structures']}/{{receptor}}_receptor.pdbqt"
    conda:
        "envs/iedb_tools.yaml"
    shell:
        """
        awk 'BEGIN {{FS=""}} /^ATOM/ && (substr($0, 22, 1)=="A" || substr($0, 22, 1)=="B") {{print $0}}' {input.pdb} > {MHC_DIR}/{wildcards.receptor}_clean_target.pdb
        python scripts/run_adfrsuite.py \
            --tool prepare_receptor \
            --output {output.pdbqt} \
            --config config.yaml \
            -- \
            -r {MHC_DIR}/{wildcards.receptor}_clean_target.pdb \
            -o {output.pdbqt} \
            -A hydrogens \
            -U nsw
        rm {MHC_DIR}/{wildcards.receptor}_clean_target.pdb
        """

rule adcp_generate_affinity_grid:
    input:
        receptor_pdbqt = f"{PATHS['mhc_structures']}/{{receptor}}_receptor.pdbqt",
        peptide_pdbqt  = f"{PATHS['mhc_structures']}/{{receptor}}_native_peptide.pdbqt"
    output:
        trg = f"{PATHS['mhc_structures']}/{{receptor}}_receptor.trg"
    conda:
        "envs/iedb_tools.yaml"
    shell:
        """
        python scripts/run_adfrsuite.py \
            --tool agfr \
            --output {output.trg} \
            --config config.yaml \
            -- \
            -r {input.receptor_pdbqt} \
            -l {input.peptide_pdbqt} \
            -b ligand \
            -P {params.padding} \
            -o {MHC_DIR}/{wildcards.receptor}_receptor
        """

rule adcp_screening:
    input:
        csv = f"{CONSENSUS_DIR}/filtered_mhc1_epitopes.csv",
        # Ensures all receptors in your config are fully mapped and ready before docking
        trgs = expand(f"{PATHS['mhc_structures']}/{{receptor}}_receptor.trg", receptor=config.get("docking", {}).get("receptors", []))
    output:
        summary = f"{PATHS['docking']}/adcp_screening_summary.txt"
    params:
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
            --replicas {params.replicas} \
            --steps {params.steps} \
            --threads {params.threads} \
            --config config.yaml
        """
# =============================================================================
# MODULE 9: VACCINE ASSEMBLY CHIMERA
# =============================================================================

rule construct_vaccine:
    input:
        overlapping = f"{CONSENSUS_DIR}/overlapping_consensus_regions.csv",
        mhc1_clean  = f"{CONSENSUS_DIR}/filtered_mhc1_epitopes.csv",
        mhc2_clean  = f"{CONSENSUS_DIR}/filtered_mhc2_epitopes.csv",
        bcell_clean = f"{CONSENSUS_DIR}/filtered_bcell_epitopes.csv"
    output:
        fasta = f"{VACCINE_DIR}/final_vaccine_sequence.fasta"
    conda:
        "envs/canine_immuno.yaml"
    shell:
        """
        python scripts/assemble_vaccine_chimera.py \
            --intersected {input.overlapping} \
            --mhc1 {input.mhc1_clean} \
            --mhc2 {input.mhc2_clean} \
            --bcell {input.bcell_clean} \
            --output {output.fasta} \
            --config config.yaml
        """