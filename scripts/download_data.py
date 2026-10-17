#!/usr/bin/env python3
"""Download configured CPV-2 accessions and the UniProt canine reference proteome."""

import argparse
import os
import tempfile
import time
from pathlib import Path
from urllib.error import URLError
from urllib.parse import quote
from urllib.request import Request, urlopen

from pipeline_config import load_config, resolve_path


def validate_fasta(path):
    saw_header = False
    saw_residue = False
    with path.open(encoding="utf-8") as fasta:
        for line in fasta:
            stripped = line.strip()
            if stripped.startswith(">"):
                saw_header = True
            elif stripped and saw_header:
                saw_residue = True
    if not saw_header or not saw_residue:
        raise ValueError(f"Downloaded file is not a non-empty FASTA: {path}")


def download_fasta(url, destination, force=False, attempts=3):
    destination = Path(destination)
    if destination.exists() and not force:
        validate_fasta(destination)
        print(f"[SKIP] Existing FASTA: {destination}")
        return

    destination.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(1, attempts + 1):
        temp_path = None
        try:
            request = Request(url, headers={"User-Agent": "Ekaghni-strikes/0.3.0"})
            with urlopen(request, timeout=90) as response:
                with tempfile.NamedTemporaryFile(
                    mode="wb", dir=destination.parent, delete=False
                ) as temp_file:
                    temp_path = Path(temp_file.name)
                    while chunk := response.read(1024 * 1024):
                        temp_file.write(chunk)
            validate_fasta(temp_path)
            os.replace(temp_path, destination)
            print(f"[DOWNLOADED] {destination}")
            return
        except (OSError, URLError, ValueError) as error:
            if temp_path and temp_path.exists():
                temp_path.unlink()
            if attempt == attempts:
                raise RuntimeError(
                    f"Failed to download a valid FASTA from {url} after {attempts} attempts: {error}"
                ) from error
            time.sleep(attempt)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=None, help="Pipeline YAML config (defaults to config.yaml)")
    parser.add_argument("--force", action="store_true", help="Replace existing FASTA downloads")
    args = parser.parse_args()

    config = load_config(args.config)
    sources = config["data_sources"]
    for accession in config["samples"]:
        accession = str(accession)
        url = sources["genome_url_template"].format(accession=quote(accession, safe=""))
        destination = resolve_path(config["paths"]["raw_genomes"], config) / f"{accession}.fasta"
        download_fasta(url, destination, force=args.force)

    dog_proteome = resolve_path(config["databases"]["dog_proteome_fasta"], config)
    download_fasta(sources["dog_proteome_url"], dog_proteome, force=args.force)
    print("Reference genome and dog proteome downloads are ready.")


if __name__ == "__main__":
    main()
