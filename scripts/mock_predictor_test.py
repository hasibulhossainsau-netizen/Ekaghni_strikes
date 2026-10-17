#!/usr/bin/env python3
"""Exercise external-tool mock outputs and their immediate downstream parsers."""

import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd
import yaml

from pipeline_config import PROJECT_ROOT, load_config


def run_script(script, *args):
    subprocess.run(
        [sys.executable, str(PROJECT_ROOT / "scripts" / script), *map(str, args)],
        cwd=PROJECT_ROOT,
        check=True,
    )


def main():
    with tempfile.TemporaryDirectory(prefix="ekaghni-mock-tools-") as temp:
        temp_dir = Path(temp)
        config = load_config()
        for tool_name in (
            "netmhcpan", "bepipred3", "toxinpred3", "algpred2",
            "prepare_ligand", "prepare_receptor", "agfr", "adcp",
        ):
            config["tools"][tool_name] = {"path": str(temp_dir / f"missing-{tool_name}")}
        config["paths"]["mhc_structures"] = str(temp_dir / "missing-structures")
        config["paths"]["docking_screening"] = str(temp_dir / "docking-screening")
        config["docking"]["receptors"] = ["TEST"]
        config_path = temp_dir / "config.yaml"
        config_path.write_text(yaml.safe_dump(config))

        protein_fasta = temp_dir / "proteins.fasta"
        protein_sequence = "MKWVTFISLLFLFSSAYS" * 4
        protein_fasta.write_text(f">test_protein\n{protein_sequence}\n")
        netmhc_csv = temp_dir / "netmhcpan.csv"
        run_script(
            "netmhcpan_runner.py",
            protein_fasta, "DLA-TEST", "8,9", 0.5, netmhc_csv,
            "--config", config_path,
        )
        netmhc = pd.read_csv(netmhc_csv)
        assert {"protein_id", "allele", "peptide", "start", "end", "length", "ic50", "percentile_rank"} <= set(netmhc.columns)
        assert not netmhc.empty
        assert netmhc["percentile_rank"].nunique() > 1

        toxin_csv = temp_dir / "toxinpred3.csv"
        run_script("run_toxinpred3.py", "--input", netmhc_csv, "--output", toxin_csv, "--config", config_path)
        toxin = pd.read_csv(toxin_csv)
        assert {"peptide", "toxicity_prediction", "toxicity_score"} <= set(toxin.columns)
        assert {"protein_id", "percentile_rank"} <= set(toxin.columns)
        assert len(toxin) == len(netmhc)
        assert toxin["toxicity_score"].nunique() > 1

        allergen_csv = temp_dir / "algpred2.csv"
        run_script("run_algpred2.py", "--input", toxin_csv, "--output", allergen_csv, "--config", config_path)
        allergen = pd.read_csv(allergen_csv)
        assert {"peptide", "allergenicity_prediction", "allergenicity_score"} <= set(allergen.columns)
        assert len(allergen) == len(toxin)
        assert {"protein_id", "percentile_rank"} <= set(allergen.columns)
        assert allergen["allergenicity_score"].nunique() > 1

        bcell_dir = temp_dir / "bcell"
        run_script(
            "run_bepipred3.py",
            "--input", protein_fasta,
            "--output-dir", bcell_dir,
            "--prediction-mode", "vt_pred",
            "--threshold", "0.15",
            "--window-size", "9",
            "--config", config_path,
        )
        raw_bcell = pd.read_csv(bcell_dir / "raw_output.csv")
        assert {"Accession", "BepiPred-3.0 score", "BepiPred-3.0 linear epitope score"} <= set(raw_bcell.columns)
        assert raw_bcell["BepiPred-3.0 score"].nunique() > 1
        bcell_summary = temp_dir / "bcell-summary.csv"
        run_script(
            "extract_bcell_coordinates.py",
            bcell_dir / "Bcell_linepitope_top_20pct_preds.fasta",
            bcell_dir / "raw_output.csv",
            bcell_summary,
            "--config", config_path,
        )
        assert not pd.read_csv(bcell_summary).empty

        for tool_name, output_name in (
            ("prepare_ligand", "peptide.pdbqt"),
            ("prepare_receptor", "receptor.pdbqt"),
            ("agfr", "receptor_receptor.trg"),
        ):
            output = temp_dir / output_name
            run_script(
                "run_adfrsuite.py",
                "--tool", tool_name,
                "--output", output,
                "--config", config_path,
                "--", "-o", output,
            )
            assert output.is_file() and output.stat().st_size > 0

        fake_adfr = temp_dir / "fake-adfr"
        fake_adfr.write_text(
            "#!/usr/bin/env python3\n"
            "import sys\n"
            "from pathlib import Path\n"
            "Path(sys.argv[sys.argv.index('-o') + 1]).write_text('fake ADFRsuite result\\n')\n"
        )
        fake_adfr.chmod(0o755)
        config["tools"]["prepare_ligand"] = {"path": str(fake_adfr)}
        config_path.write_text(yaml.safe_dump(config))
        forwarded_output = temp_dir / "forwarded.pdbqt"
        run_script(
            "run_adfrsuite.py",
            "--tool", "prepare_ligand",
            "--output", forwarded_output,
            "--config", config_path,
            "--", "-l", "input.pdb", "-o", forwarded_output,
        )
        assert forwarded_output.read_text() == "fake ADFRsuite result\n"

        adcp_csv = temp_dir / "adcp-input.csv"
        pd.DataFrame({"peptide": ["PEPTIDE"]}).to_csv(adcp_csv, index=False)
        adcp_report = temp_dir / "adcp-report.txt"
        run_script(
            "submit_adcp.py",
            "--csv", adcp_csv,
            "--output", adcp_report,
            "--config", config_path,
        )
        assert all(column in adcp_report.read_text() for column in ("Receptor", "Peptide", "Affinity_Score", "Status"))

        print("Mock predictor tests passed.")


if __name__ == "__main__":
    main()
