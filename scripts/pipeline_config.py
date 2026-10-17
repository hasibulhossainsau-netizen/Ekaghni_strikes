"""Shared configuration loading and path resolution for workflow scripts."""

import os
import shutil
from pathlib import Path

import yaml


PROJECT_ROOT = Path(__file__).resolve().parent.parent


def load_config(config_path=None):
    path = Path(config_path) if config_path else PROJECT_ROOT / "config.yaml"
    if not path.is_absolute():
        path = PROJECT_ROOT / path
    with path.open() as config_file:
        config = yaml.safe_load(config_file)
    if not isinstance(config, dict):
        raise ValueError(f"Configuration must contain a YAML mapping: {path}")
    return config


def resolve_path(path, config=None):
    """Resolve configured relative paths against the repository, not the caller's cwd."""
    resolved = Path(path).expanduser()
    if not resolved.is_absolute():
        resolved = PROJECT_ROOT / resolved
    return resolved


def configured_tool(config, name):
    try:
        value = config["tools"][name]
    except KeyError as exc:
        raise KeyError(f"Missing tools.{name} in config.yaml") from exc
    if isinstance(value, dict):
        value = value.get("path") or value.get("command")
    if not value:
        raise ValueError(f"Set tools.{name}.path or tools.{name}.command in config.yaml")
    return str(resolve_path(value, config)) if "/" in str(value) or str(value).startswith("~") else str(value)


def configured_tool_if_available(config, name):
    """Return a configured executable when it resolves to an executable file."""
    value = config.get("tools", {}).get(name)
    if isinstance(value, dict):
        value = value.get("path") or value.get("command")
    if not value:
        return None

    executable = configured_tool(config, name)
    if "/" in executable or executable.startswith("~"):
        path = resolve_path(executable, config)
        if path.is_file() and (path.suffix == ".py" or os.access(path, os.X_OK)):
            return str(path)
        return None
    return shutil.which(executable)
