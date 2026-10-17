#!/usr/bin/env python3
"""Check workflow-driver imports and required dependencies in rule environments."""

import importlib.util
import re
import os
import sys
from pathlib import Path
import shutil

import yaml


ROOT = Path(__file__).resolve().parent.parent

# The smoke-test and environment-check scripts use PyYAML in the driver
# interpreter. Analysis packages are checked against their Snakemake env files.
DRIVER_MODULES = {"yaml": "pyyaml"}

# Dependencies needed by scripts in each Snakemake rule environment.
RULE_ENV_MODULES = {
    "blast": ("Bio", "pandas", "yaml"),
    "alignment": ("Bio", "pandas", "yaml"),
    "netmhcpan": ("Bio", "numpy", "pandas", "yaml"),
    "iedb_tools": ("pandas", "requests", "yaml"),
    "toxinpred3": ("Bio", "pandas", "sklearn", "yaml"),
    "algpred2": ("Bio", "pandas", "sklearn", "yaml"),
    "bepipred3": ("Bio", "numpy", "pandas", "yaml", "torch", "esm"),
    "deeptmhmm": ("Bio", "numpy", "pandas", "biolib", "yaml"),
    "canine_immuno": ("Bio", "numpy", "pandas", "yaml"),
    "mhc_vax": ("numpy", "pandas", "yaml"),
}

MODULE_PACKAGES = {
    "Bio": "biopython",
    "numpy": "numpy",
    "pandas": "pandas",
    "requests": "requests",
    "yaml": "pyyaml",
    "sklearn": "scikit-learn",
    "torch": "pytorch",
    "esm": "fair-esm",
    "biolib": "pybiolib",
}

MOCK_TEST_MODULES = ("Bio", "pandas", "yaml")


def declared_package_names(environment_path):
    environment = yaml.safe_load(environment_path.read_text())
    package_names = set()
    for dependency in environment.get("dependencies", []):
        if isinstance(dependency, str):
            package_names.add(re.split(r"[<>=!~]", dependency, maxsplit=1)[0].lower())
        elif isinstance(dependency, dict):
            for pip_dependency in dependency.get("pip", []):
                package_names.add(
                    re.split(r"[<>=!~]", pip_dependency, maxsplit=1)[0].lower()
                )
    return package_names


def main():
    errors = []
    use_conda = os.environ.get("VACDES_USE_CONDA", "1") == "1"
    check_scope = os.environ.get("VACDES_CHECK_SCOPE")
    runtime_modules = (
        MOCK_TEST_MODULES
        if check_scope == "mock-test"
        else tuple(MODULE_PACKAGES)
    )
    if sys.version_info < (3, 10):
        errors.append("Python 3.10 or newer is required to run workflow scripts.")

    for module, package in DRIVER_MODULES.items():
        if importlib.util.find_spec(module) is None:
            errors.append(f"Missing driver dependency: install {package} for {sys.executable}.")

    if shutil.which("snakemake") is None:
        errors.append("Snakemake is not on PATH.")

    if not use_conda:
        for module in runtime_modules:
            package = MODULE_PACKAGES[module]
            if importlib.util.find_spec(module) is None:
                errors.append(
                    f"Missing shared-runtime dependency: install {package} "
                    f"for {sys.executable} or enable Snakemake Conda environments."
                )

    for environment, modules in RULE_ENV_MODULES.items():
        environment_path = ROOT / "envs" / f"{environment}.yaml"
        if not environment_path.is_file():
            errors.append(f"Missing rule environment file: {environment_path.relative_to(ROOT)}.")
            continue
        declared = declared_package_names(environment_path)
        for module in modules:
            package = MODULE_PACKAGES[module]
            if package not in declared:
                errors.append(
                    f"envs/{environment}.yaml must declare {package} "
                    f"(required by the workflow scripts)."
                )
            if (
                not use_conda
                and module in runtime_modules
                and importlib.util.find_spec(module) is None
            ):
                errors.append(
                    f"Missing shared-runtime dependency for {environment}: "
                    f"install {package} for {sys.executable}."
                )

    if errors:
        print("Workflow environment check failed:", file=sys.stderr)
        for error in errors:
            print(f"  - {error}", file=sys.stderr)
        return 1

    print(f"Workflow environment verified with {sys.executable}.")
    print("Required imports and per-rule environment dependencies are present.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
