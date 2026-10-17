#!/usr/bin/env python3
"""Run an isolated Snakemake DAG dry run using the bundled minimal genome."""

from copy import deepcopy
from pathlib import Path
import shutil
import subprocess
import tempfile

import yaml


ROOT = Path(__file__).resolve().parent.parent


def configure_smoke_run(config, temporary_root, genome_path, dog_proteome):
    config = deepcopy(config)
    config["samples"] = ["test_sample", "test_query"]
    config["reference"]["genome"] = str(genome_path)
    config["reference"]["vp2_gene"] = {"start": 1, "end": 180}
    config["qc"]["minimum_sequence_length"] = 0

    paths = config["paths"]
    paths.update({
        "raw_genomes": str(genome_path.parent),
        "cleaned_genomes": str(temporary_root / "fasta_cleaned"),
        "proteins": str(temporary_root / "protein"),
        "alignment": str(temporary_root / "align"),
        "blast": str(temporary_root / "blast"),
        "filtering": str(temporary_root / "filtering"),
        "immunogenicity": str(temporary_root / "immunogenicity"),
        "bcell": str(temporary_root / "bcell"),
        "safety": str(temporary_root / "safety"),
        "topology": str(temporary_root / "topology"),
        "consensus": str(temporary_root / "consensus"),
        "vaccine_construct": str(temporary_root / "vaccine_construct"),
        "docking": str(temporary_root / "docking"),
        "mhc_structures": str(temporary_root / "mhc_structures"),
        "master_protein_fasta": str(temporary_root / "protein/master_pool.fasta"),
        "reference_vp2_fasta": str(genome_path.parent / "test_sample_VP2.fasta"),
        "aligned_proteins": str(temporary_root / "align/aligned.fasta"),
        "conservation_csv": str(temporary_root / "align/conservation.csv"),
        "mhc1_predictions_csv": str(temporary_root / "align/mhc1.csv"),
        "topology_csv": str(temporary_root / "topology/topology.csv"),
        "final_ranked_tcell_csv": str(temporary_root / "final_ranked.csv"),
        "docking_screening": str(temporary_root / "docking/screening"),
        "host_exclusion_csv": str(temporary_root / "filtering/host_exclusion.csv"),
        "host_exclusion_report": str(temporary_root / "filtering/report.json"),
    })
    config["databases"]["dog_proteome_fasta"] = str(dog_proteome)
    config["databases"]["dog_blast_db"] = str(temporary_root / "blast/dog_db")
    config["ensemble"]["breed_weights_file"] = str(ROOT / "config/breed_weights.csv")
    return config


def main():
    genome = ROOT / "data/fasta/test_sample.fasta"
    if not genome.is_file():
        raise FileNotFoundError(f"Bundled smoke-test genome is missing: {genome}")

    with tempfile.TemporaryDirectory(prefix="vacdes-smoke-") as temp_dir:
        temporary_root = Path(temp_dir)
        genome_dir = temporary_root / "input"
        genome_dir.mkdir()
        genome_path = genome_dir / "test_sample.fasta"
        shutil.copy2(genome, genome_path)
        query_genome = genome_dir / "test_query.fasta"
        shutil.copy2(genome, query_genome)
        dog_proteome = temporary_root / "dog_proteome.fa"
        dog_proteome.write_text(">test_dog_protein\nMKWVTFISLLFLFSSAYS\n")
        (temporary_root / "outputs/mhc_structures").mkdir(parents=True)

        with (ROOT / "config.yaml").open() as stream:
            config = yaml.safe_load(stream)
        smoke_config = configure_smoke_run(
            config, temporary_root / "outputs", genome_path, dog_proteome
        )

        config_path = temporary_root / "config.smoke.yaml"
        config_path.write_text(yaml.safe_dump(smoke_config, sort_keys=False))

        command = [
            "snakemake",
            "--snakefile", str(ROOT / "Snakefile"),
            "--directory", str(ROOT),
            "--configfile", str(config_path),
            "--cores", "1",
            "--dry-run",
        ]
        print("Running isolated workflow DAG dry run with test_sample.fasta...")
        result = subprocess.run(
            command, cwd=ROOT, capture_output=True, text=True, check=False
        )
        output = result.stdout + result.stderr
        if result.returncode:
            print(output, end="")
            return result.returncode
        if "This was a dry-run" not in output:
            print(output, end="")
            raise RuntimeError("Snakemake did not report completing a DAG dry run.")
        print("Smoke-test DAG dry run passed. No jobs were executed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
