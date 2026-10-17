#!/usr/bin/env python3
"""Run one ADFRsuite executable, or create a downstream-compatible mock file."""

import argparse
import subprocess
from pathlib import Path

from pipeline_config import configured_tool_if_available, load_config, resolve_path


def write_mock_output(tool, output_path):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if tool in {"prepare_ligand", "prepare_receptor"}:
        output_path.write_text(
            "ATOM      1  C   ALA A   1       0.000   0.000   0.000"
            "  1.00  0.00     0.000 C\nEND\n"
        )
    else:
        output_path.write_text("MOCK ADFRsuite affinity grid\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tool", required=True, choices=(
        "prepare_ligand", "prepare_receptor", "agfr",
    ))
    parser.add_argument("--output", required=True, help="Expected tool output file")
    parser.add_argument("--config", default=None)
    parser.add_argument("--mock", action="store_true", help="Generate a mock output file")
    parser.add_argument("tool_args", nargs=argparse.REMAINDER)
    args = parser.parse_args()

    config = load_config(args.config)
    output_path = resolve_path(args.output, config)
    executable = configured_tool_if_available(config, args.tool)
    use_mock = args.mock or executable is None
    if use_mock:
        reason = "mock mode requested" if args.mock else f"tools.{args.tool} is unavailable in config.yaml"
        print(f"[WARNING] {reason}; writing a mock ADFRsuite output.")
        write_mock_output(args.tool, output_path)
        return

    command = [executable, *args.tool_args]
    if command[1:2] == ["--"]:
        command.pop(1)
    try:
        subprocess.run(command, check=True)
    except FileNotFoundError:
        print(f"[WARNING] {args.tool} disappeared before execution; writing a mock output.")
        write_mock_output(args.tool, output_path)


if __name__ == "__main__":
    main()
