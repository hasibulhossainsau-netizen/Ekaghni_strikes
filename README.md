# Ekaghni_strikes: CPV2 Canine Vaccine Design Pipeline

[![CI](https://github.com/hasibulhossainsau-netizen/Ekaghni_strikes/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/hasibulhossainsau-netizen/Ekaghni_strikes/actions/workflows/ci.yml)
[![Snakemake](https://img.shields.io/badge/snakemake-≥7.0-brightgreen)](https://snakemake.readthedocs.io)
[![Python](https://img.shields.io/badge/python-3.10-blue)](https://www.python.org/)
[![Pipeline Version](https://img.shields.io/badge/version-0.3.0-orange)](https://github.com/hasibulhossainsau-netizen/Ekaghni_strikes)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)


> A **Snakemake-based reverse vaccinology pipeline** for multi-epitope vaccine candidate discovery against **Canine Parvovirus Type 2 (CPV-2)** VP2 capsid protein.

---

## Table of Contents

- [Project Overview](#project-overview)
- [Pipeline Architecture](#pipeline-architecture)
- [Prerequisites](#prerequisites)
- [Installation & Setup](#installation--setup)
- [Configuration Guide](#configuration-guide)
- [Usage](#usage)
- [Input & Output Description](#input--output-description)
- [Workflow Summary](#workflow-summary)
- [Estimated Runtime](#estimated-runtime)
- [Troubleshooting & Common Errors](#troubleshooting--common-errors)
- [Contribution Guidelines](#contribution-guidelines)
- [License & Citation](#license--citation)

---

## Project Overview

**Biological Target:** VP2 capsid protein of CPV-2 (genomic coordinates 2771–4535 on reference `FJ197845`)  
**Host Species:** *Canis lupus familiaris* (domestic dog)  
**Pipeline Version:** v0.3.0

This pipeline takes raw CPV-2 genome FASTA files as input and outputs a ranked list of epitope candidates scored across immunogenicity, safety, conservation, and host-specificity criteria. It uses canine DLA alleles where NetMHCpan supports them. IEDB does not currently expose canine DLA models for MHC-II binding or antigen processing, so those two analyses use explicitly documented HLA surrogate alleles; these are proxies, not validated canine equivalents.

**Key capabilities:**
- Automatic VP2 gene extraction via BLASTn from multi-sample genome inputs
- Host homology filtering using a dog proteome BLAST database to remove self-similar peptides
- Parallel MHC Class I (NetMHCpan 4.1), MHC Class II (IEDB), and B-cell (BepiPred 3.0) epitope prediction
- Dual safety screening via ToxinPred3 (toxicity) and AlgPred2 (allergenicity)
- Breed-aware ensemble immunogenicity scoring weighted across Labrador, German Shepherd, and mixed populations
- Switchable docking engine: **local ADCP** (AutoDock CrankPep) or **remote HPepDock** via a single config toggle
- Automated vaccine chimera assembly from consensus epitope regions

---

## Pipeline Architecture

```
Genome FASTA (multi-sample)
        │
        ▼
[1] VP2 Extraction (BLASTn)
        │
        ▼
[2] Sequence QC & Sanitation
        │
        ▼
[3] MAFFT Multiple Sequence Alignment (--linsi, maxiterate 1000)
        │
        ├─────────────────────────────────────────────────┐
        ▼                                                 ▼
[4a] Conservation Analysis                     [4b] Build Canine BLAST DB
(Shannon entropy, residue freq)                (dog_proteome.fa → makeblastdb)
        │                                                 │
        └─────────────┬───────────────────────────────────┘
                      ▼
        [5] Host Homology Filter (BLASTp)
              (exclude self-similar peptides)
                      │
        ┌─────────────┼──────────────────┐
        ▼             ▼                  ▼
[6a] NetMHCpan   [6b] IEDB MHC II    [6c] BepiPred 3.0
  MHC Class I      Class II binding    B-cell epitope
  (DLA alleles)    (15-mer, nn_align)  prediction
        │             │                  │
        ▼             ▼                  ▼
[7] Safety Screening (ToxinPred3 → AlgPred2 — parallel per epitope class)
        │
        ▼
[8] Ensemble Immunogenicity Scoring
   (MHC-II weight: 40% | Processing: 35% | Breed coverage: 25%)
        │
        ▼
[9] Final Weighted Ranking
   results/final_ranked_tcell_epitopes.csv
        │
        ▼
[10] Consensus Epitope Intersection & Docking
   (ADCP local  OR  HPepDock remote — config toggle)
        │
        ▼
[11] Vaccine Chimera Assembly
   results/vaccine_construct/final_vaccine_sequence.fasta
```

---

## Prerequisites

### System Requirements

| Requirement | Minimum | Recommended |
|-------------|---------|-------------|
| OS | Ubuntu 20.04 / macOS 12 | Ubuntu 22.04 LTS |
| RAM | 8 GB | 16 GB+ |
| CPU cores | 4 | 8+ |
| Disk space | 20 GB | 50 GB |
| Conda/Mamba | Any recent | Mamba (faster) |
| Snakemake | ≥ 7.0 | ≥ 8.0 |

### External Tools and Optional Manual Installation

NetMHCpan, BepiPred 3.0, and ADFRsuite require separate downloads and are not
included in the Docker image. In the default test/mock configuration they may be
absent. For scientific runs, install them and configure their executable paths.

**1. NetMHCpan 4.1** (requires free academic registration)
- Download from [DTU Health Tech](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/)
- After extracting, edit the launch script to set `NMHOME`:
```bash
# Inside the NetMHCpan launcher, ensure NMHOME points to the extracted
# NetMHCpan installation directory.
setenv NMHOME /path/to/netMHCpan-4.1
```

**2. ADFRsuite v1.1dev** *(required only when `docking.tool: "adcp"` is set)*
- Download from [ADFR Suite](https://ccsb.scripps.edu/adfr/downloads/)
- Verify that these three binaries are executable:
```bash
/your/ADFRsuite/bin/prepare_receptor
/your/ADFRsuite/bin/prepare_ligand
/your/ADFRsuite/bin/agfr
/your/ADFRsuite/bin/adcp
```

**3. BepiPred 3.0** (for B-cell epitope prediction)
- Follow instructions at [BepiPred 3.0 GitHub](https://github.com/mnielLab/bepipred3)
- Note the full path to `bepipred3_CLI.py` — you will need it for `config.yaml`

---

## Installation & Setup

### Step 1: Clone the repository

```bash
git clone https://github.com/hasibulhossainsau-netizen/Ekaghni_strikes.git
cd Ekaghni_strikes
```

### Step 2: Place input data

Download the configured CPV-2 genomes and the canine reference proteome:

```bash
python scripts/download_data.py
```

The helper downloads each accession in `config.yaml` from NCBI Nucleotide and
the canine reference proteome UP000002254 from UniProt. Existing valid FASTA
files are kept; pass `--force` to replace them. You can customize sources and
destination paths in `config.yaml`. Alternatively, place matching genome files
in `data/fasta/` and the canine reference at `data/reference/dog_proteome.fa`.

The configured genome files have this layout (one accession per FASTA):

```
data/fasta/
├── FJ197845.fasta      ← Reference strain (required)
├── PV611245.fasta
├── PV611246.fasta
└── PV611247.fasta
```

BLAST indexes are generated by the workflow from the downloaded proteome.

The repository includes a tiny synthetic genome at `data/fasta/test_sample.fasta`
for checking workflow setup. It is only a smoke-test fixture, not biological input.

### Step 3: Place MHC structural data *(ADCP docking only)*

When running with `docking.tool: "adcp"`, place receptor PDB files here:

```
data/MHC_I/
├── 5F1I.pdb            ← Receptor crystal structure
└── 7CJQ.pdb            ← Receptor crystal structure
```

> **Note:**  
> PDB files can be downloaded from [RCSB](https://www.rcsb.org). Each file must contain Chain A & B (receptor) and Chain C (native peptide).

### Step 4: Configure paths in `config.yaml`

Configure executable paths and database locations in `config.yaml`. Relative paths are resolved from the repository root, and bare executable names are looked up on `PATH` (including the Conda environment):

```yaml
tools:
  netmhcpan:
    path: "/path/to/netMHCpan-4.1/netMHCpan"
  bepipred3:
    path: "/path/to/bepipred3_CLI.py"
  prepare_receptor:
    path: "/path/to/ADFRsuite/bin/prepare_receptor"
  prepare_ligand:
    path: "/path/to/ADFRsuite/bin/prepare_ligand"
  agfr:
    path: "/path/to/ADFRsuite/bin/agfr"
  adcp:
    path: "/path/to/ADFRsuite/bin/adcp"

databases:
  dog_proteome_fasta: "data/reference/dog_proteome.fa"
  dog_blast_db: "data/reference/dog_db"
```

### Step 5: Install Conda environments

Run the fixture-based DAG check, then verify the workflow environment declarations:

```bash
make smoke-test
make check-env
make test
make mock-test
```

After placing real input genomes and configuring external tools, run the pipeline
with `make run`. Snakemake creates the per-rule Conda environments automatically.

Or install individual environments manually:

```bash
for env in blast alignment netmhcpan iedb_tools toxinpred3 algpred2 \
           bepipred3 deeptmhmm canine_immuno mhc_vax; do
  conda env create -f envs/${env}.yaml
done
```

### Quick workflow checks

```bash
make smoke-test  # Uses the bundled fixture and an isolated temporary config/output tree.
make check-env   # Checks the workflow driver and required packages declared per rule env.
make test        # Runs pytest unit tests for conservation, host filtering, and ranking.
make mock-test   # Checks mock predictors, output schemas, and downstream parsing.
make run         # Checks dependencies, then executes Snakemake with --use-conda.
```

The smoke test performs a complete Snakemake DAG dry run only; it does not launch
predictors, contact external services, or create pipeline outputs. The dependency
check verifies driver prerequisites and confirms that workflow environments declare
the Python packages used by their scripts (including Biopython); `make run` lets
Snakemake provision those per-rule environments from `envs/`.

Run `make mock-test` to exercise external-predictor fallback output files and
their immediate downstream parsers. The IEDB API wrappers also write synthetic
scores if remote requests fail or reject the requested surrogate alleles.
Docker Compose's default test mode runs the DAG dry run, unit tests, and
mock-predictor checks.

NetMHCpan, BepiPred 3, ToxinPred3, AlgPred2, and ADFRsuite wrappers accept
`--mock`. They also switch to deterministic, sequence-dependent synthetic outputs
with a warning when their configured executable is absent. ToxinPred3 and
AlgPred2 are installed in the Docker runtime but use the same fallback if their
commands are unavailable. ADFRsuite fallback files are structural placeholders.
Mock scores exist only to exercise DAG execution and output schemas; they are not
scientific predictions and must never be interpreted or reported as biological
results. If a tool is present and no `--mock` flag is supplied, its real
executable is used.

### Run with Docker Compose

Docker Compose provides a prebuilt runtime containing Snakemake, the workflow's
Python packages, BLAST+, and MAFFT. No local Python, Conda, or bioinformatics
packages are needed. By default, Compose builds the image and runs the safe smoke
test, unit tests, and external-tool mock checks:

```bash
docker-compose up --build
```

To execute the full pipeline instead:

```bash
PIPELINE_MODE=full docker-compose up --build
```

Full mode runs against the repository's `config.yaml`; acquire inputs first with
`python scripts/download_data.py`. Manually licensed/downloaded tools
(NetMHCpan, BepiPred 3.0, and ADFRsuite for local ADCP docking) are not bundled.
The Compose file mounts `${EXTERNAL_TOOLS_DIR:-./tools}` read-only at
`/opt/external-tools`. Put the installations under `tools/` or set the variable
to another host directory:

```bash
EXTERNAL_TOOLS_DIR=/srv/bio-tools PIPELINE_MODE=full docker-compose up --build
```

Use container-visible paths in `config.yaml`, for example
`/opt/external-tools/netMHCpan-4.1/netMHCpan`,
`/opt/external-tools/bepipred3/bepipred3_CLI.py`, and
`/opt/external-tools/ADFRsuite-1.1/bin/adcp` (also configure its other three
binaries). These mounts are optional for smoke/mock runs. Generated files are
written to the mounted repository. For non-default host account IDs, set
`LOCAL_UID` and `LOCAL_GID` so generated files have the expected ownership.

---

## Configuration Guide

All executable paths, database locations, and pipeline parameters are controlled via `config.yaml`. The full, annotated [config.example.yaml](config.example.yaml) is provided as a starting template. The in-repository Python scripts accept `--config` (defaulting to the repository's `config.yaml`), and the workflow passes it explicitly. Relative paths are rooted at the repository. The `tools` section supports either an explicit `path` or a bare `command` resolved through `PATH`; configure NetMHCpan and the ADFRsuite executables explicitly when they are not available on `PATH`.

### Docking Engine Toggle

```yaml
# The most important switch — choose local or remote docking
docking:
  tool: "adcp"        # Options: "adcp" (local) | "hpepdock" (remote server)
  receptors:
    - "5F1I"
    - "7CJQ"
```

| Setting | `"adcp"` | `"hpepdock"` |
|---------|-----------|--------------|
| Execution | Local (your machine) | Remote (Huang Lab server) |
| Requires | ADFRsuite binaries | Internet + valid email |
| Output | `results/docking/adcp_screening_summary.txt` | `results/docking/hpepdock_submissions.txt` |
| Speed | Slower (CPU-bound) | Faster (offloaded) |

### Sample List

```yaml
samples:
  - FJ197845    # Must match filename in data/fasta/ (without .fasta extension)
  - PV611245
```

Add any additional CPV-2 accession IDs here. Corresponding FASTA files must exist in `data/fasta/`.

### Allele coverage and HLA surrogate limitations

```yaml
# NetMHCpan uses the configured canine DLA-I alleles.
tcell_epitope:
  dla_alleles:
    - "DLA-8803401"
mhc2_epitope:
  # IEDB currently lacks supported canine DLA-II models.
  hla_surrogate_alleles: ["HLA-DRB1*01:01"]
mhc1_processing:
  # IEDB processing also uses human HLA-I surrogates; not canine allele calls.
  hla_surrogate_alleles: ["HLA-A*02:01"]
```

The HLA-I/II values above are computational surrogates because the corresponding
IEDB DLA models are unavailable; they are not biologically validated
equivalencies. Interpret those results accordingly.

### Ensemble Scoring Weights

```yaml
# Total weight must be 1
ensemble:
  weight_mhc2: 0.40                # MHC Class II antigenicity
  weight_processing: 0.35          # MHC Class I processing efficiency
  weight_population_coverage: 0.25 # Breed-aware DLA frequency coverage
```

### Host Homology Filter Thresholds

```yaml
host_homology:
  thresholds:
    identity: 80       # BLASTp identity % — peptides above this are excluded
    coverage: 70       # Query coverage % threshold
    evalue: 1e-3
```

These three thresholds are applied to BLASTp results and control host-hit exclusion.

### Ranking Weights

```yaml
# Final ranking 
ranking:
  weights:
    mhc_binding: 0.25
    mhc2_binding: 0.20
    processing_efficiency: 0.10
    conservation: 0.15
    immunogenicity: 0.20
    safety: 0.10
```

---

## Usage

### Dry run (verify DAG without executing)

```bash
# Dry run first - No file will be created
snakemake -n --use-conda
```

### Full pipeline run

```bash
snakemake --use-conda --cores 4
```

### Run a specific rule or target file

```bash
# Run up to and including the ensemble scoring step
snakemake results/immunogenicity/ensemble_immunogenicity.csv --use-conda --cores 4

# Run only safety screening for MHC1 epitopes
snakemake results/safety/mhc1_algpred2_checked.csv --use-conda --cores 2
```

### Visualize the workflow DAG

```bash
snakemake --dag | dot -Tsvg > dag.svg
open dag.svg    # macOS
xdg-open dag.svg  # Linux
```

### Force re-run of specific rules

```bash
# Re-run conservation analysis even if output exists
snakemake data/align/VP2_conservation.csv --use-conda --cores 4 --forcerun conservation
```

### Clean all outputs

```bash
snakemake --delete-all-output
```

---

## Input & Output Description

### Required Input Files

| File | Location | Description |
|------|----------|-------------|
| `{sample}.fasta` | `data/fasta/` | Raw CPV-2 genome FASTA per accession |
| `FJ197845.fasta` | `data/fasta/` | Reference genome (mandatory) |
| `dog_proteome.fa` | `data/reference/` | Canine proteome for host homology filtering |
| `{receptor}.pdb` | `data/MHC_I/` | MHC-I crystal structures (ADCP mode only) |
| `breed_weights.csv` | `config/` | DLA allele frequencies by breed |

### Key Output Files

| File | Description |
|------|-------------|
| `results/final_ranked_tcell_epitopes.csv` | **Primary output** — T-cell epitopes ranked by composite score |
| `results/vaccine_construct/final_vaccine_sequence.fasta` | **Final output** — assembled multi-epitope vaccine chimera |
| `results/consensus/overlapping_consensus_regions.csv` | Epitopes shared across MHC-I, MHC-II, and B-cell predictions |
| `results/immunogenicity/ensemble_immunogenicity.csv` | Weighted ensemble immunogenicity scores |
| `results/filtering/host_homology_exclusion.csv` | Peptides excluded due to dog proteome similarity |
| `data/align/VP2_conservation.csv` | Per-position Shannon entropy conservation scores |
| `data/align/vp2_strong_binders.csv` | NetMHCpan strong binders (rank ≤ 0.5%) |
| `results/safety/mhc1_algpred2_checked.csv` | MHC-I epitopes post toxicity + allergenicity filter |
| `results/docking/adcp_screening_summary.txt` | ADCP docking results (when `tool: "adcp"`) |
| `results/docking/hpepdock_submissions.txt` | HPepDock submission log (when `tool: "hpepdock"`), you will receive the result via email from HPEPDOCK |
| `results/host_exclusion_report.json` | Summary statistics of host filtering |

### Final Ranking Columns

| Column | Description |
|--------|-------------|
| `peptide` | Amino acid sequence |
| `best_rank` | NetMHCpan best rank (%) across all DLA alleles |
| `n_strong_alleles` | Number of DLA alleles with strong binding |
| `mhc_score` | Normalized MHC binding score |
| `conservation_score` | Shannon entropy-derived conservation (higher = more conserved) |
| `surface_exposure` | DeepTMHMM surface exposure probability |
| `toxicity_prediction` | ToxinPred3 class: `Non-Toxin` / `Toxin` |
| `allergenicity_prediction` | AlgPred2 class: `Non-Allergen` / `Allergen` |
| `host_penalty` | BLASTp similarity penalty score against dog proteome |
| `composite_score` | **Final weighted score — higher is better** |

### Project Directory Structure

```
Ekaghni_strikes/
├── Snakefile                        # Main workflow definition
├── config.yaml                      # All parameters and paths
├── README.md
├── data/
│   ├── fasta/                       # Input genomes (one per accession)
│   ├── fasta_cleaned/               # QC-passed sequences
│   ├── protein/                     # Pooled VP2 protein pool
│   ├── align/                       # MSA and conservation outputs
│   ├── blast/                       # BLASTp result files
│   ├── MHC_I/                       # Receptor PDB + PDBQT files (ADCP)
│   └── reference/                   # Reference genome + dog proteome
├── results/
│   ├── filtering/                   # Host homology exclusion results
│   ├── immunogenicity/              # MHC-I, MHC-II, ensemble scores
│   ├── safety/                      # Toxicity + allergenicity outputs
│   ├── topology/                    # DeepTMHMM transmembrane topology
│   ├── bcell/                       # BepiPred 3.0 B-cell predictions
│   ├── consensus/                   # Filtered and intersected epitopes
│   ├── docking/                     # ADCP or HPepDock outputs
│   └── vaccine_construct/           # Final FASTA vaccine sequence
├── scripts/                         # All Python analysis scripts
├── config/
│   ├── dla_alleles.txt              # DLA alleles + breed frequencies
│   └── breed_weights.csv
└── envs/                            # Conda environment YAML files
```

---

## Workflow Summary

The Snakefile is organized into 9 logical modules:

**Module 1 — Data Preparation**  
Extracts VP2 coding sequences from raw genomes using BLASTn against the FJ197845 reference, runs quality control (invalid amino acid removal, header validation), and builds a pooled protein FASTA for downstream analysis.

**Module 2 — Host Filtering**  
Builds a canine proteome BLAST database and runs BLASTp for all VP2 sequences. Peptides with ≥80% identity and ≥70% query coverage to dog proteins are flagged and excluded from the clean peptide universe.

**Module 3 — Conservation Analysis**  
Performs MAFFT multiple sequence alignment (--linsi, 1000 iterations) and calculates per-position Shannon entropy, residue frequency, and variability index. Only positions with ≥95% conservation across a minimum 8-residue window are retained.

**Module 4 — Epitope Discovery**  
Three parallel prediction streams run independently from the filtered sequence universe: NetMHCpan 4.1 for DLA-aware MHC Class I binding (8–11-mers), IEDB nn_align for MHC Class II (15-mers), and BepiPred 3.0 for linear B-cell epitopes.

**Module 5 — Safety Screening**  
ToxinPred3 and AlgPred2 are applied sequentially to each epitope class (MHC-I, MHC-II, B-cell) in independent parallel chains. Peptides predicted as toxic or allergenic are excluded before scoring.

**Module 6 — Functional Scoring**  
IEDB antigen processing scores and DeepTMHMM topological data are computed. These feed into a breed-aware ensemble immunogenicity score combining MHC-II binding (40%), processing efficiency (35%), and DLA population coverage (25%).

**Module 7 — Final Ranking**  
All evidence streams converge into a weighted composite score. The seven criteria (MHC binding, conservation, immunogenicity, safety, topology, processing, population coverage) are min-max normalized and combined per the weights in `config.yaml`.

**Module 8 — Intersection & Docking**  
Epitopes appearing across MHC-I, MHC-II, and B-cell predictions are identified as consensus regions. The docking tool (ADCP or HPepDock) is selected dynamically based on `config.docking.tool`. ADCP runs a full local receptor preparation → affinity grid generation → peptide screening workflow. HPepDock submits jobs to the Huang Lab remote server.

**Module 9 — Vaccine Assembly**  
Consensus overlapping regions and top-ranked epitopes from all three classes are assembled into a multi-epitope chimera FASTA using linker logic in `assemble_vaccine_chimera.py`.

## Estimated Runtime

Approximate wall-clock ranges for a small four-genome run on 4 CPU cores, with
inputs already downloaded; actual time depends on network/API availability,
system resources, predictor installation, and remote queue times. Independent
branches run concurrently. These estimates exclude Conda environment creation
and any large reference-data downloads.

| Step | Work | Typical time |
|---|---|---:|
| 1 | Input validation and genome quality control | 1–5 min |
| 2 | VP2 reference/query extraction with BLASTn | 1–10 min |
| 3 | Sequence sanitization and pooled protein preparation | <1–3 min |
| 4 | Dog BLAST database creation and host-homology filtering | 1–10 min |
| 5 | MAFFT alignment and conservation scoring | 1–10 min |
| 6 | MHC-I binding prediction (NetMHCpan) | 1–20 min |
| 7 | MHC-II binding prediction (IEDB; remote service) | 1–15 min |
| 8 | B-cell prediction (BepiPred 3) | 1–20 min |
| 9 | Safety, processing, topology, and ensemble scoring | 2–30 min |
| 10 | Candidate ranking, consensus, and reporting | <1–5 min |
| 11 | Docking and vaccine construct assembly | 1–60+ min |

Smoke/mock checks usually finish in a few minutes and skip scientific prediction
services. HPepDock server wait times and real-tool setup/runtime are especially
variable; a full run should be budgeted at tens of minutes to several hours.

---

## Troubleshooting & Common Errors

### `FileNotFoundError` for `netMHCpan`

**Cause:** `tools.netmhcpan.path` in `config.yaml` is empty or points to an invalid executable, or `NMHOME` inside the NetMHCpan installation is not set to its installation directory.

```bash
# Verify that the executable configured as tools.netmhcpan.path is executable
test -x /path/to/netMHCpan-4.1/netMHCpan
```

---

### ADCP docking: `No such file or directory` for `.trg` file

**Cause:** The `agfr` grid generation step failed silently, often because the receptor PDB contains multiple assemblies or non-standard chain labels.

```bash
# Verify chain layout of your PDB
grep "^ATOM" data/MHC_I/5F1I.pdb | awk '{print substr($0,22,1)}' | sort -u
# Expected output: A B C
# Chain A+B = receptor, Chain C = native peptide
```

---

### BLASTp DB error: `BLAST Database error: No alias or index file found`

**Cause:** The `makeblastdb` rule did not complete or the `.phr/.pin/.psq` files are missing.

```bash
# Rebuild the BLAST database manually
makeblastdb -in data/reference/dog_proteome.fa -dbtype prot -out data/reference/dog_db
ls data/reference/dog_db.*   # Should show .phr, .pin, .psq
```

---

### MHC-II alleles mismatch warning in IEDB output

**Cause:** The `config.yaml` uses `HLA-DRB1` alleles as proxies (IEDB does not yet have a full DLA-II database). This is expected behavior. The pipeline uses these as phylogenetically informative surrogates.

---

### HPepDock: `ConnectionError` or no response after submission

**Cause:** The Huang Lab server may be under load or temporarily unavailable.

```bash
# Increase the timeout and retry count in config.yaml
docking:
  max_retries: 5
  timeout_seconds: 300
```

---

### Snakemake locks the working directory

```bash
# If a previous run crashed and left a lock:
snakemake --unlock
```

---

### conda env creation fails with `ResolvePackageNotFound`

```bash
# Use mamba for faster and more reliable dependency resolution
conda install -n base -c conda-forge mamba
mamba env create -f envs/canine_immuno.yaml
```

---

## Contribution Guidelines

Contributions are welcome from the computational biology and veterinary vaccinology community.

**How to contribute:**

1. Fork the repository on GitHub.
2. Create a feature branch: `git checkout -b feature/your-feature-name`
3. Keep configuration in `config.yaml`; update `config.example.yaml` and related
   documentation when adding/changing settings. Do not commit local tool paths,
   credentials, or downloaded/licensed tools.
4. Follow existing Python conventions and use explicit input/output schemas;
   update the relevant Conda or Docker dependencies when introducing packages.
5. Add/update focused pytest coverage for behavior and schema changes.
6. Run `make test`, `make check-env`, `make smoke-test`, and `make mock-test`.
7. Open a pull request describing the change, validation performed, and
   biological interpretation or limitations where applicable.

**Areas where contributions are especially valued:**

- Adding support for additional DLA allele databases as they become publicly available
- Integrating structural quality metrics (e.g., Rosetta binding energy) into the ensemble score
- Adding support for CPV-2 variants beyond VP2 (e.g., VP1)
- Improving Windows/WSL compatibility for the ADCP preparation steps
- Writing automated tests (`pytest`) for the scoring scripts

The built-in mock outputs are deterministic synthetic fixtures for software
testing only. Contributions must not present mock outputs as real predictions.

**Reporting Issues:**  
Open a GitHub issue with your OS, Snakemake version, the full error traceback, and the output of `snakemake --version`.

---

## License & Citation

This pipeline is released under the **MIT License**. See [LICENSE](LICENSE) for full terms.
Citation metadata is provided in [CITATION.cff](CITATION.cff); GitHub can use
this file to generate a formatted software citation.

### Please cite the following tools when publishing results from this pipeline:

| Tool | Reference |
|------|-----------|
| **Snakemake** | Mölder et al. (2021). *F1000Research.* https://doi.org/10.12688/f1000research.29032.2 |
| **NetMHCpan 4.1** | Reynisson et al. (2020). *Nucleic Acids Research.* https://doi.org/10.1093/nar/gkaa379 |
| **MAFFT** | Katoh & Standley (2013). *Mol. Biol. Evol.* https://doi.org/10.1093/molbev/mst010 |
| **BLAST+** | Camacho et al. (2009). *BMC Bioinformatics.* https://doi.org/10.1186/1471-2105-10-421 |
| **BepiPred 3.0** | Clifford et al. (2022). *Protein Science.* https://doi.org/10.1093/nar/gkac476 |
| **DeepTMHMM** | Hallgren et al. (2022). *bioRxiv.* https://doi.org/10.1101/2022.04.08.487609 |
| **ToxinPred3** | Rathore et al. (2023). *Computers in Biology and Medicine.* https://doi.org/10.1016/j.compbiomed.2024.108926 |
| **AlgPred2** | Sharma et al. (2021). *Briefings in Bioinformatics.* https://doi.org/10.1093/bib/bbaa294 |
| **ADCP / ADFRsuite** | Zhang et al. (2019). *Bioinformatics.* https://doi.org/10.1093/bioinformatics/btz459 |
| **HPepDock** | Zhou et al. (2018). *Nucleic Acids Research.* https://doi.org/10.1093/nar/gky357 |

---

*Developed as part of independent computational veterinary vaccinology research . For questions or collaborations, open an issue on GitHub.*  