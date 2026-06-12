"""Configuration loading for the LBD pipeline.

All tunable thresholds live in ``config/default.json``. A user can override any
subset with ``--config path/to/overrides.json``; the override file is deep-merged
on top of the defaults so partial overrides are allowed.
"""
from __future__ import annotations

import json
import os
from copy import deepcopy
from typing import Any, Dict

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(_HERE)
DEFAULT_CONFIG_PATH = os.path.join(_REPO_ROOT, "config", "default.json")


def _deep_merge(base: Dict[str, Any], override: Dict[str, Any]) -> Dict[str, Any]:
    out = deepcopy(base)
    for key, value in override.items():
        if key in out and isinstance(out[key], dict) and isinstance(value, dict):
            out[key] = _deep_merge(out[key], value)
        else:
            out[key] = deepcopy(value)
    return out


def load_config(override_path: str | None = None) -> Dict[str, Any]:
    with open(DEFAULT_CONFIG_PATH, "r", encoding="utf-8") as fh:
        cfg = json.load(fh)
    if override_path:
        with open(override_path, "r", encoding="utf-8") as fh:
            override = json.load(fh)
        cfg = _deep_merge(cfg, override)
    # Allow an API key from the environment without editing config files.
    env_key = os.environ.get("NCBI_API_KEY")
    if env_key:
        cfg["eutils"]["api_key"] = env_key
        # With a key NCBI permits 10 req/s.
        cfg["eutils"]["requests_per_second"] = max(cfg["eutils"]["requests_per_second"], 10)
    env_email = os.environ.get("NCBI_EMAIL")
    if env_email:
        cfg["eutils"]["email"] = env_email
    return cfg


def repo_root() -> str:
    return _REPO_ROOT
