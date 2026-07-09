"""Shared catchment/source paths for SUMMA post-processing NetCDF products."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

Product = Literal["mixed_precip", "snow_metrics"]

DEFAULT_GPEP_ROOT = Path("/anvil/projects/x-ees240082/users/dcasson/gpep")

# Input/output layout differs slightly between products (see notebooks).
_SOURCE_LAYOUT: dict[Product, dict[str, dict[str, object]]] = {
    "mixed_precip": {
        "rf": {
            "label": "RF Ensemble",
            "input_subdir": "rf_best_regression_static_dynamic/summa/output/complete",
            "output_suffix": "rf_mixed_precip",
            "deterministic": False,
        },
        "casr": {
            "label": "CASR",
            "input_subdir": "casr_remapped/summa_output/complete",
            "output_suffix": "casr_mixed_precip",
            "deterministic": True,
        },
        "era5": {
            "label": "ERA5",
            "input_subdir": "era5_remapped/summa_output/complete",
            "output_suffix": "era5_mixed_precip",
            "deterministic": True,
        },
    },
    "snow_metrics": {
        "rf": {
            "label": "RF Ensemble",
            "input_subdir": "rf_best_regression_static_dynamic/summa/output/complete",
            "output_suffix": "rf_snow_metrics",
            "deterministic": False,
        },
        "casr": {
            "label": "CASR",
            "input_subdir": "casr_remapped/summa_output",
            "output_suffix": "casr_snow_metrics",
            "deterministic": True,
        },
        "era5": {
            "label": "ERA5",
            "input_subdir": "era5_remapped/summa_output",
            "output_suffix": "era5_snow_metrics",
            "deterministic": True,
        },
    },
}

_OUTPUT_SUBDIR = {
    "mixed_precip": "mixed_precip",
    "snow_metrics": "snow_metrics",
}


def ensemble_sources_for_catchment(
    catchment: str,
    *,
    product: Product,
    gpep_root: Path | str = DEFAULT_GPEP_ROOT,
    source_keys: list[str] | None = None,
) -> list[dict]:
    """Return source specs for one catchment and product."""
    gpep_root = Path(gpep_root)
    base = gpep_root / catchment / "ensemble_generation"
    out_subdir = _OUTPUT_SUBDIR[product]
    layout = _SOURCE_LAYOUT[product]

    keys = source_keys or list(layout.keys())
    unknown = [k for k in keys if k not in layout]
    if unknown:
        raise ValueError(f"Unknown source keys {unknown}; choose from {list(layout.keys())}")

    specs: list[dict] = []
    for key in keys:
        meta = layout[key]
        specs.append(
            {
                "key": key,
                "label": meta["label"],
                "input_dir": base / str(meta["input_subdir"]),
                "output_nc": base / out_subdir / f"{catchment}_{meta['output_suffix']}.nc",
                "deterministic": bool(meta["deterministic"]),
            }
        )
    return specs


def output_paths_for_catchment(
    catchment: str,
    *,
    product: Product,
    gpep_root: Path | str = DEFAULT_GPEP_ROOT,
    source_keys: list[str] | None = None,
) -> dict[str, Path]:
    """Map source label → saved NetCDF path for one catchment."""
    return {
        spec["label"]: Path(spec["output_nc"])
        for spec in ensemble_sources_for_catchment(
            catchment,
            product=product,
            gpep_root=gpep_root,
            source_keys=source_keys,
        )
    }


def build_output_nc_manifest(
    catchments: list[str],
    *,
    product: Product,
    gpep_root: Path | str = DEFAULT_GPEP_ROOT,
    source_keys: list[str] | None = None,
) -> dict[str, dict[str, Path]]:
    """Map catchment → source label → NetCDF path (no I/O)."""
    return {
        catchment: output_paths_for_catchment(
            catchment,
            product=product,
            gpep_root=gpep_root,
            source_keys=source_keys,
        )
        for catchment in catchments
    }
