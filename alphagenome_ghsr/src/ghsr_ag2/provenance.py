"""Provenance and reproducibility recording.

Captures software versions, model/API identifiers, database access dates and
per-step execution metadata. Never records the API key.
"""

from __future__ import annotations

import datetime
import importlib.metadata
import json
import platform
import sys
from typing import Any

from . import config

_PACKAGES = (
    "alphagenome",
    "anndata",
    "grpcio",
    "matplotlib",
    "numpy",
    "pandas",
    "pyarrow",
    "requests",
    "scipy",
    "seaborn",
)

ALPHAGENOME_API_HOST = "gdmscience.googleapis.com:443"
ALPHAGENOME_MODEL_FAMILY = "AlphaGenome 2 (alphagenome.models.v2)"


def utc_now() -> str:
    return datetime.datetime.now(datetime.timezone.utc).isoformat(
        timespec="seconds"
    )


def package_versions() -> dict[str, str]:
    out: dict[str, str] = {}
    for name in _PACKAGES:
        try:
            out[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            out[name] = "not installed"
    return out


def environment_record() -> dict[str, Any]:
    return {
        "recorded_utc": utc_now(),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "packages": package_versions(),
        "alphagenome_model_family": ALPHAGENOME_MODEL_FAMILY,
        "alphagenome_api_host": ALPHAGENOME_API_HOST,
        "genome_build": config.GENOME_BUILD,
        "api_key_source": f"environment variable {config.API_KEY_ENV} "
        "(value never recorded)",
    }


def write(name: str, payload: dict[str, Any]) -> None:
    """Write a provenance record to data/raw/provenance_<name>.json."""
    record = {"step": name, **environment_record(), **payload}
    path = config.DATA_RAW / f"provenance_{name}.json"
    path.write_text(json.dumps(record, indent=2, default=str))
    print(f"[provenance] wrote {path}")
