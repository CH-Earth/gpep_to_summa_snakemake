"""Load per-GRU snow metrics CSV exports for comparison with model NetCDF outputs."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

import compute_water_year_snow_metrics as wysm

_CSV_METRIC_PREFIX = {
    wysm.VAR_SCF_DAYS: "SCF",
    wysm.VAR_SNOW_DISAPPEARANCE: "SDD",
}

_COL_PATTERN = re.compile(r"^(SCF|SDD)_(\d{4})_mean$")
_COL_PATTERN_STAT = re.compile(r"^(SCF|SDD)_(\d{4})_(mean|median|max|min)$")

_CSV_ID_FIELDS = ("HRU_ID", "GRU_ID", "gruId", "gru_id", "hruId", "hru_id", "fid")


def obs_csv_metric_prefix(var: str) -> str:
    if var not in _CSV_METRIC_PREFIX:
        raise KeyError(
            f"Unsupported var {var!r}; choose from {list(_CSV_METRIC_PREFIX)}"
        )
    return _CSV_METRIC_PREFIX[var]


def list_obs_csv_water_years(
    csv_path: Path | str,
    var: str,
    *,
    stat: str = "mean",
) -> list[int]:
    """Return sorted water years available for ``var`` and ``stat`` in the CSV."""
    prefix = obs_csv_metric_prefix(var)
    csv_path = Path(csv_path)
    header = pd.read_csv(csv_path, nrows=0).columns
    years: list[int] = []
    for col in header:
        m = _COL_PATTERN_STAT.match(col)
        if m and m.group(1) == prefix and m.group(3) == stat:
            years.append(int(m.group(2)))
    years.sort()
    if not years:
        raise ValueError(
            f"No {prefix}_YYYY_{stat} columns found in {csv_path}; "
            f"columns sample: {list(header[:8])}"
        )
    return years


def infer_csv_id_field(df: pd.DataFrame) -> str:
    for field in _CSV_ID_FIELDS:
        if field in df.columns:
            return field
    raise ValueError(
        "CSV has no GRU/HRU id column; "
        f"expected one of {_CSV_ID_FIELDS}; got {list(df.columns[:12])}"
    )


def load_obs_csv_spatial_values(
    csv_path: Path | str,
    var: str,
    *,
    water_year: int | None = None,
    stat: str = "mean",
    id_field: str | None = None,
) -> pd.DataFrame:
    """
    Per-polygon values from a wide GEE snow-metrics CSV export.

    Returns a DataFrame with columns ``[id_field, "value"]``.

    When ``water_year`` is None, values are averaged across all available years
    for the requested ``stat`` (e.g. ``SCF_2009_mean`` … ``SCF_2019_mean``).
    """
    prefix = obs_csv_metric_prefix(var)
    csv_path = Path(csv_path)
    if not csv_path.is_file():
        raise FileNotFoundError(f"Obs CSV not found: {csv_path}")

    df = pd.read_csv(csv_path)
    join_field = id_field or infer_csv_id_field(df)
    if join_field not in df.columns:
        raise KeyError(f"{csv_path} missing id column {join_field!r}")

    if water_year is not None:
        col = f"{prefix}_{int(water_year)}_{stat}"
        if col not in df.columns:
            available = list_obs_csv_water_years(csv_path, var, stat=stat)
            raise KeyError(
                f"{csv_path} missing column {col!r}; available years: {available}"
            )
        values = pd.to_numeric(df[col], errors="coerce")
    else:
        years = list_obs_csv_water_years(csv_path, var, stat=stat)
        cols = [f"{prefix}_{year}_{stat}" for year in years]
        values = df[cols].apply(pd.to_numeric, errors="coerce").mean(axis=1)

    out = pd.DataFrame(
        {
            join_field: pd.to_numeric(df[join_field], errors="coerce").astype("Int64"),
            "value": values.astype(float),
        }
    )
    return out.dropna(subset=[join_field])


def load_obs_csv_gru_mean_series(
    csv_path: Path | str,
    var: str,
    *,
    gru_ids: np.ndarray | list[int] | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Return (water_years, values) averaged across GRU rows in the CSV.

    Pass ``gru_ids`` to restrict the mean to a subset of polygons (e.g. high-elevation GRUs).

    Expects columns like ``SCF_2009_mean`` and ``SDD_2009_mean``.
    """
    if var not in _CSV_METRIC_PREFIX:
        raise KeyError(
            f"Unsupported var {var!r}; choose from {list(_CSV_METRIC_PREFIX)}"
        )
    prefix = _CSV_METRIC_PREFIX[var]
    csv_path = Path(csv_path)
    if not csv_path.is_file():
        raise FileNotFoundError(f"Obs CSV not found: {csv_path}")

    df = pd.read_csv(csv_path)
    if gru_ids is not None:
        join_field = infer_csv_id_field(df)
        gru_ids = np.asarray(gru_ids, dtype=np.int64)
        df = df[df[join_field].astype(np.int64).isin(gru_ids)]
        if df.empty:
            raise ValueError(
                f"No CSV rows match gru_ids filter ({len(gru_ids)} ids) in {csv_path}"
            )

    pairs: list[tuple[int, str]] = []
    for col in df.columns:
        m = _COL_PATTERN.match(col)
        if m and m.group(1) == prefix:
            pairs.append((int(m.group(2)), col))
    pairs.sort(key=lambda x: x[0])
    if not pairs:
        raise ValueError(
            f"No {prefix}_YYYY_mean columns found in {csv_path}; "
            f"columns sample: {list(df.columns[:8])}"
        )

    years = np.array([year for year, _ in pairs], dtype=int)
    cols = [col for _, col in pairs]
    vals = np.asarray(df[cols].values, dtype=float)
    gru_mean = np.nanmean(vals, axis=0)
    return years, gru_mean
