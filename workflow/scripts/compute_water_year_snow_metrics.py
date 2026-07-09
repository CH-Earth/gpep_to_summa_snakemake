"""
Water-year snow metrics from SUMMA timestep outputs.

Input layouts
-------------

1. One NetCDF per GRU
   Filenames match ``*_G<gru>-<gru>_timestep.nc``. Each file holds a single
   HRU/GRU time series.

2. All GRUs/HRUs in each file
   No ``_G<id>-<id>_timestep`` pattern. ``scalarSWE`` uses a ``gru`` or ``hru``
   dimension. Example ensemble files:

       chena_001_timestep.nc
       chena_002_timestep.nc
       chena_003_timestep.nc

   With ``deterministic=False``, the member id is inferred from the filename
   suffix before ``_timestep``. For example, ``chena_002_timestep.nc`` becomes
   member ``002``.

Water-year definition
---------------------

US/NOAA water year:

    Oct 1 through Sep 30

The water year is labeled by the calendar year in which it ends. For example:

    WY 2024 = 2023-10-01 through 2024-09-30

Metrics
-------

1. Snow cover frequency, SCF

   Fraction of valid daily values in the water year where daily maximum
   ``scalarSWE`` is greater than ``swe_threshold``.

2. Snow disappearance day, SDD

   Numeric day of water year, not datetime.

   By default:

       Oct 1 = day 1

   The SDD algorithm searches backward within each water year and finds the
   latest transition where:

       run_days consecutive valid snow days
       are immediately followed by
       run_days consecutive valid snow-free days

   The reported value is the day-of-water-year of the first snow-free day.

Important hourly-data behavior
------------------------------

Inputs are assumed hourly by default. Hourly/subdaily SWE is first aggregated to
daily maximum SWE. A daily value is only considered valid when it has at least
``min_daily_timesteps`` non-missing timestep values.

Missing or insufficiently sampled days are not treated as snow-free.

This avoids the common bug where:

    NaN > threshold

evaluates to False and missing data accidentally becomes no-snow.

Output SDD encoding
-------------------

SDD is stored as float32 day-of-water-year values. Missing/undetected SDD is
stored as NaN.

This avoids the datetime NaT integer artifact:

    -9223372036854775808
"""

from __future__ import annotations

import logging
import re
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple, Union

import numpy as np
import pandas as pd
import xarray as xr

# Allow ``import compute_mixed_precip_fractions`` when cwd is not workflow/scripts.
_SCRIPTS_DIR = Path(__file__).resolve().parent
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

from compute_mixed_precip_fractions import (
    ParsedTimestepName,
    _expand_gru,
    _expand_member,
    _squeeze_hru,
    _time_encoding_key,
    discover_nc_files,
    parse_timestep_nc_path,
)

logger = logging.getLogger(__name__)

VAR_SWE = "scalarSWE"
"""Name of SWE in SUMMA timestep NetCDF files."""

VAR_SCF_DAYS = "snow_cover_frequency_days"
"""
Snow cover frequency.

The historical variable name is retained for compatibility, but the value is now
a fraction, not a day count.
"""

VAR_SNOW_DISAPPEARANCE = "snow_disappearance_day_of_water_year"
"""Snow disappearance day as numeric day of water year. Oct 1 = 1 by default."""

VAR_VALID_DAYS = "valid_snow_metric_days"
"""Number of valid daily SWE values used per water year."""

_GRU_FILENAME = re.compile(r"_G(\d+)-(\d+)_timestep\.nc$", re.IGNORECASE)


def classify_timestep_path_layout(paths: Sequence[Path]) -> str:
    """
    Classify timestep file layout.

    Returns
    -------
    str
        ``"per_gru_filenames"`` if every path matches
        ``*_G<gru>-<gru>_timestep.nc``.

        ``"multi_gru_files"`` if no path matches that pattern.

    Raises
    ------
    ValueError
        If file naming conventions are mixed.
    """
    matched = [_GRU_FILENAME.search(p.name) is not None for p in paths]

    if all(matched):
        return "per_gru_filenames"

    if not any(matched):
        return "multi_gru_files"

    bad = [p.name for p, m in zip(paths, matched) if m != matched[0]]

    raise ValueError(
        "Mixed timestep filename conventions in one batch. All files must either "
        "match *_G<gru>-<gru>_timestep.nc or none must. "
        f"Examples: {bad[:3]}"
    )


def infer_member_label_from_path(path: Path, *, deterministic: bool) -> Optional[str]:
    """
    Infer ensemble member id from filename.

    Examples
    --------
    ``chena_002_timestep.nc`` -> ``"002"``

    ``chena_010_timestep.nc`` -> ``"010"``

    If ``deterministic=True``, returns None.
    """
    if deterministic:
        return None

    stem = Path(path).stem

    if stem.lower().endswith("_timestep"):
        stem = stem[: -len("_timestep")]

    if "_" in stem:
        return stem.rsplit("_", 1)[1]

    return stem if stem else "unknown"


def us_water_year_label(time: xr.DataArray) -> xr.DataArray:
    """
    NOAA-style US water-year index.

    Oct 1 through Sep 30, labeled by the calendar year in which the period ends.

    Examples
    --------
    2023-10-01 -> 2024

    2024-09-30 -> 2024

    2024-10-01 -> 2025
    """
    year = time.dt.year
    month = time.dt.month
    return xr.where(month >= 10, year + 1, year).astype(np.int32)


def _water_year_start_date(water_year: int) -> np.datetime64:
    """
    Return Oct 1 start date for a NOAA-style water year.

    Example
    -------
    WY 2024 starts on 2023-10-01.
    """
    return np.datetime64(f"{int(water_year) - 1:04d}-10-01", "D")


def _water_year_end_date(water_year: int) -> np.datetime64:
    """Return Sep 30 end date for a NOAA-style water year."""
    return np.datetime64(f"{int(water_year):04d}-09-30", "D")


def _day_of_water_year(
    date: Union[np.datetime64, pd.Timestamp],
    *,
    water_year: int,
    one_based: bool = True,
) -> float:
    """
    Convert a calendar date to numeric day of water year.

    With one_based=True:

        Oct 1 = 1

    With one_based=False:

        Oct 1 = 0
    """
    d = np.datetime64(date, "D")
    wy_start = _water_year_start_date(water_year)
    offset = int((d - wy_start) / np.timedelta64(1, "D"))

    if one_based:
        return float(offset + 1)

    return float(offset)


def _source_1d_values(
    source: xr.Dataset,
    name: str,
    dim: str,
) -> Optional[np.ndarray]:
    """
    Return 1D values for ``name`` from a source dataset when it exists on ``dim``.

    Handles both coordinates and data variables.
    """
    if name not in source:
        return None

    da = source[name]

    if da.dims != (dim,):
        return None

    return np.asarray(da.values)


def _open_timestep_dataset(path: Path) -> xr.Dataset:
    """Open a SUMMA timestep NetCDF with a clearer error for corrupt HDF5 files."""
    try:
        return xr.open_dataset(path)
    except RuntimeError as exc:
        if "HDF error" not in str(exc):
            raise
        hint = (
            f"Could not read NetCDF file (corrupt or truncated HDF5): {path}. "
            "If a sibling ``complete/`` directory exists, its merged files may be "
            "incomplete; retry with ``prefer_multi_gru_files=False`` and per-GRU "
            "filenames such as ``*_G<id>-<id>_timestep.nc``."
        )
        raise RuntimeError(hint) from exc


def load_one_timestep_swe(path: Union[str, Path]) -> xr.Dataset:
    """
    Load ``scalarSWE`` from one timestep file in the one-file-per-GRU layout.

    A single HRU is squeezed by the helper imported from
    ``compute_mixed_precip_fractions``.
    """
    path = Path(path)

    with _open_timestep_dataset(path) as ds:
        if VAR_SWE not in ds:
            raise KeyError(
                f"{path} must contain {VAR_SWE!r}; found {list(ds.data_vars)}"
            )

        swe = _squeeze_hru(ds[VAR_SWE].load(), path)

    out = xr.Dataset({VAR_SWE: swe.astype(np.float32)})
    out[VAR_SWE].attrs.update({k: v for k, v in swe.attrs.items() if k != "_FillValue"})

    return out


def normalize_swe_to_gru_dimension(
    swe: xr.DataArray,
    path: Path,
    *,
    source_ds: Optional[xr.Dataset] = None,
) -> xr.DataArray:
    """
    Ensure SWE has a ``gru`` dimension.

    SUMMA files often store ``scalarSWE(time, hru)`` even when ``gru`` and
    ``gruId`` variables are also present. For this module, the spatial dimension
    is normalized to ``gru`` for downstream stacking.

    If ``hruId`` or ``gruId`` are available as 1D variables, they are copied onto
    the output ``gru`` coordinate where possible.
    """
    path = Path(path)
    source = source_ds

    if "gru" in swe.dims and "hru" in swe.dims:
        raise ValueError(
            f"{path}: {VAR_SWE} has both 'gru' and 'hru' dimensions; expected one spatial dim."
        )

    if "gru" in swe.dims:
        out = swe

        if source is not None:
            gid = _source_1d_values(source, "gruId", "gru")
            if gid is not None:
                out = out.assign_coords(gru=("gru", gid.astype(np.int64)))
                return out

        if "gruId" in out.coords and out.coords["gruId"].dims == ("gru",):
            gid = np.asarray(out.coords["gruId"].values)
            out = out.assign_coords(gru=("gru", gid.astype(np.int64)))

        return out

    if "hru" in swe.dims:
        out = swe.rename({"hru": "gru"})

        if source is not None:
            hid = _source_1d_values(source, "hruId", "hru")
            if hid is not None:
                out = out.assign_coords(gru=("gru", hid.astype(np.int64)))
                return out

            # Some SUMMA files have one HRU per GRU. Use gruId only if it has
            # the same length as the hru dimension.
            gid = _source_1d_values(source, "gruId", "gru")
            if gid is not None and gid.size == out.sizes["gru"]:
                out = out.assign_coords(gru=("gru", gid.astype(np.int64)))
                return out

        if "hruId" in swe.coords and swe.coords["hruId"].dims == ("hru",):
            hid = np.asarray(swe.coords["hruId"].values)
            out = out.assign_coords(gru=("gru", hid.astype(np.int64)))
            return out

        out = out.assign_coords(
            gru=("gru", np.arange(out.sizes["gru"], dtype=np.int64))
        )
        return out

    # Scalar site: keep explicit spatial dimension.
    return swe.expand_dims(gru=[np.int64(0)])


def load_timestep_swe_multi_gru_file(path: Union[str, Path]) -> xr.Dataset:
    """
    Load ``scalarSWE`` from a file that may contain many GRUs/HRUs.

    Example expected input:

        scalarSWE(time, hru)

    Output:

        scalarSWE(time, gru)
    """
    path = Path(path)

    with _open_timestep_dataset(path) as ds:
        if VAR_SWE not in ds:
            raise KeyError(
                f"{path} must contain {VAR_SWE!r}; found {list(ds.data_vars)}"
            )

        swe = normalize_swe_to_gru_dimension(
            ds[VAR_SWE].load(),
            path,
            source_ds=ds,
        )

    out = xr.Dataset({VAR_SWE: swe.astype(np.float32)})
    out[VAR_SWE].attrs.update({k: v for k, v in swe.attrs.items() if k != "_FillValue"})

    return out


def _assert_same_gru_coord(a: xr.Dataset, b: xr.Dataset, pa: Path, pb: Path) -> None:
    """Require identical ``gru`` coordinate when concatenating along time."""
    ga = a[VAR_SWE].coords["gru"].values
    gb = b[VAR_SWE].coords["gru"].values

    if ga.shape != gb.shape or not np.array_equal(ga, gb):
        raise ValueError(
            "When joining multi-GRU files along time, the gru coordinate must match. "
            f"Mismatch between {pa!s} and {pb!s}."
        )


def _sort_by_time_if_present(ds: xr.Dataset) -> xr.Dataset:
    """Sort a dataset by time when a time coordinate is present."""
    if "time" not in ds.dims and "time" not in ds.coords:
        return ds

    try:
        return ds.sortby("time")
    except Exception:
        logger.warning("Could not sort dataset by time.", exc_info=True)
        return ds


def merge_timestep_swe_files_multi_gru_in_file(
    paths: Sequence[Path],
    *,
    deterministic: bool = False,
    time_join: str = "outer",
    force_member_dim: bool = True,
    reference_time_encoding: Optional[Tuple[Any, Any, Any]] = None,
) -> xr.Dataset:
    """
    Merge files where each file contains ``scalarSWE(time, gru)``.

    Files are grouped by inferred ensemble member. Within each member, multiple
    files are concatenated along time. Across members, datasets are concatenated
    along ``member``.

    For files like ``chena_002_timestep.nc``, use ``deterministic=False`` so
    member ``002`` is preserved.
    """
    if time_join not in ("inner", "outer"):
        raise ValueError("time_join must be 'inner' or 'outer'")

    ref_enc = reference_time_encoding

    by_member: Dict[Optional[str], List[Path]] = defaultdict(list)

    for p in paths:
        member_id = infer_member_label_from_path(p, deterministic=deterministic)
        by_member[member_id].append(Path(p))

    member_blocks: List[xr.Dataset] = []

    for member_id in sorted(by_member.keys(), key=lambda m: (m is None, str(m))):
        plist = sorted(by_member[member_id])
        pieces: List[xr.Dataset] = []

        for filepath in plist:
            piece = load_timestep_swe_multi_gru_file(filepath)

            if ref_enc is None:
                ref_enc = _time_encoding_key(piece)
            else:
                cur = _time_encoding_key(piece)
                if cur != ref_enc:
                    logger.warning(
                        "Time encoding differs from first file: %s vs ref %s (%s)",
                        cur,
                        ref_enc,
                        filepath,
                    )

            if pieces:
                prev_path = plist[len(pieces) - 1]
                _assert_same_gru_coord(pieces[-1], piece, prev_path, filepath)

            pieces.append(piece)

        if len(pieces) == 1:
            ds_member = pieces[0]
        else:
            ds_member = xr.concat(
                pieces,
                dim="time",
                join=time_join,
                combine_attrs="drop",
            )

        ds_member = _sort_by_time_if_present(ds_member)

        if member_id is not None:
            ds_member = _expand_member(ds_member, member_id)

        member_blocks.append(ds_member)

    if not member_blocks:
        raise ValueError("No input files to merge.")

    if len(member_blocks) == 1:
        out = member_blocks[0]
    else:
        out = xr.concat(
            member_blocks,
            dim="member",
            join=time_join,
            combine_attrs="drop",
        )

    out = _sort_by_time_if_present(out)

    n_member = out.sizes.get("member", 1)
    if n_member == 1 and not force_member_dim and "member" in out.dims:
        out = out.squeeze("member", drop=True)

    if "hru" in out.dims:
        raise AssertionError("Merged SWE output must not contain an 'hru' dimension.")

    logger.info(
        "Merged multi-GRU-in-file SWE: dims=%s, n_files=%s",
        dict(out.sizes),
        len(paths),
    )

    return out


def merge_timestep_swe_files(
    paths: Sequence[Path],
    *,
    deterministic: bool = False,
    time_join: str = "outer",
    force_member_dim: bool = True,
    reference_time_encoding: Optional[Tuple[Any, Any, Any]] = None,
    multi_gru_per_file: Optional[bool] = None,
) -> xr.Dataset:
    """
    Merge timestep SWE into ``scalarSWE(time, gru[, member])``.

    If ``multi_gru_per_file`` is None, layout is inferred from filenames.
    """
    if time_join not in ("inner", "outer"):
        raise ValueError("time_join must be 'inner' or 'outer'")

    paths = [Path(p) for p in paths]

    if multi_gru_per_file is None:
        layout = classify_timestep_path_layout(paths)
        multi_gru_per_file = layout == "multi_gru_files"

    if multi_gru_per_file:
        return merge_timestep_swe_files_multi_gru_in_file(
            paths,
            deterministic=deterministic,
            time_join=time_join,
            force_member_dim=force_member_dim,
            reference_time_encoding=reference_time_encoding,
        )

    parsed: List[Tuple[Path, ParsedTimestepName]] = []

    for p in paths:
        parsed.append((p, parse_timestep_nc_path(p, deterministic=deterministic)))

    seen: Set[Tuple[Optional[str], int]] = set()

    for p, info in parsed:
        key = (info.member_id, info.gru_id)
        if key in seen:
            raise ValueError(f"Duplicate GRU/member combination for files including {p.name}")
        seen.add(key)

    ref_enc = reference_time_encoding
    groups: Dict[Optional[str], List[Tuple[Path, int]]] = defaultdict(list)

    for p, info in parsed:
        groups[info.member_id].append((p, info.gru_id))

    member_blocks: List[xr.Dataset] = []

    for member_id in sorted(groups.keys(), key=lambda m: (m is None, str(m))):
        items = sorted(groups[member_id], key=lambda t: t[1])
        gru_pieces: List[xr.Dataset] = []

        for filepath, gru_id in items:
            piece = load_one_timestep_swe(filepath)

            if ref_enc is None:
                ref_enc = _time_encoding_key(piece)
            else:
                cur = _time_encoding_key(piece)
                if cur != ref_enc:
                    logger.warning(
                        "Time encoding differs from first file: %s vs ref %s (%s)",
                        cur,
                        ref_enc,
                        filepath,
                    )

            gru_pieces.append(_expand_gru(piece, gru_id))

        ds_member = xr.concat(
            gru_pieces,
            dim="gru",
            join=time_join,
            combine_attrs="drop",
        )

        ds_member = _sort_by_time_if_present(ds_member)

        if member_id is not None:
            ds_member = _expand_member(ds_member, member_id)

        member_blocks.append(ds_member)

    if not member_blocks:
        raise ValueError("No input files to merge.")

    if len(member_blocks) == 1:
        out = member_blocks[0]
    else:
        out = xr.concat(
            member_blocks,
            dim="member",
            join=time_join,
            combine_attrs="drop",
        )

    out = _sort_by_time_if_present(out)

    n_member = out.sizes.get("member", 1)
    if n_member == 1 and not force_member_dim and "member" in out.dims:
        out = out.squeeze("member", drop=True)

    if "hru" in out.dims:
        raise AssertionError("Merged SWE output must not contain an 'hru' dimension.")

    logger.info("Merged SWE dataset: dims=%s, n_files=%s", dict(out.sizes), len(paths))

    return out


def _daily_swe_summary(
    swe: xr.DataArray,
    *,
    expected_timesteps_per_day: int = 24,
    min_daily_timesteps: Optional[int] = None,
) -> Tuple[xr.DataArray, xr.DataArray, xr.DataArray]:
    """
    Collapse hourly/subdaily SWE to daily values.

    Returns
    -------
    daily_max_swe
        Daily maximum SWE.

    daily_valid_count
        Number of non-missing timestep values contributing to each daily value.

    daily_valid
        True where the day has enough non-missing timestep values.

    Notes
    -----
    This is critical because missing data must not be interpreted as no-snow.
    """
    if "time" not in swe.dims:
        raise ValueError("SWE must have a time dimension.")

    if expected_timesteps_per_day < 1:
        raise ValueError("expected_timesteps_per_day must be >= 1")

    if min_daily_timesteps is None:
        min_daily_timesteps = expected_timesteps_per_day

    if min_daily_timesteps < 1:
        raise ValueError("min_daily_timesteps must be >= 1")

    daily_max_swe = swe.resample(time="1D").max(keep_attrs=True)
    daily_valid_count = swe.notnull().resample(time="1D").sum()
    daily_valid = daily_valid_count >= int(min_daily_timesteps)

    daily_max_swe.attrs.update({k: v for k, v in swe.attrs.items() if k != "_FillValue"})
    daily_valid_count.name = "daily_valid_timestep_count"
    daily_valid.name = "daily_valid"

    return daily_max_swe, daily_valid_count, daily_valid


def _snow_disappearance_day_for_year(
    has_snow_daily: np.ndarray,
    valid_daily: np.ndarray,
    day_times: np.ndarray,
    *,
    water_year: int,
    run_days: int = 5,
    require_complete_wy: bool = True,
    min_valid_wy_days: int = 360,
    require_terminal_snow_free: bool = False,
    allow_late_transient_snow: bool = True,
    one_based_day: bool = True,
) -> float:
    """
    Snow disappearance day for one water year.

    Returns
    -------
    float
        Numeric day of water year. Missing/undetected = np.nan.

    Logic
    -----
    Search backward through the water year. Find the latest transition where:

        run_days consecutive valid snow days
        are immediately followed by
        run_days consecutive valid no-snow days

    The reported value is the day-of-water-year of the first no-snow day.

    Missing or invalid days cannot satisfy either the snow or no-snow window.
    """
    if run_days < 1:
        raise ValueError("run_days must be >= 1")

    snow = np.asarray(has_snow_daily, dtype=bool)
    valid = np.asarray(valid_daily, dtype=bool)
    day_times = np.asarray(day_times)

    n = snow.size

    if valid.size != n:
        raise ValueError("valid_daily must have the same length as has_snow_daily")

    if day_times.size != n:
        raise ValueError("day_times must have the same length as has_snow_daily")

    if n < 2 * run_days:
        return np.nan

    if require_complete_wy and int(valid.sum()) < int(min_valid_wy_days):
        return np.nan

    wy_start = _water_year_start_date(water_year)
    wy_end = _water_year_end_date(water_year)

    # Latest possible start of the snow-free run is n - run_days.
    # Earliest possible start is run_days because the prior snow run must fit.
    for no_snow_start in range(n - run_days, run_days - 1, -1):
        prior = slice(no_snow_start - run_days, no_snow_start)
        after = slice(no_snow_start, no_snow_start + run_days)

        if not bool(valid[prior].all() and valid[after].all()):
            continue

        prior_snow = snow[prior].all()
        following_no_snow = (~snow[after]).all()

        if not bool(prior_snow and following_no_snow):
            continue

        sdd_date = np.datetime64(day_times[no_snow_start], "D")

        # Safety check: detected date must belong to the current WY.
        if sdd_date < wy_start or sdd_date > wy_end:
            continue

        if require_terminal_snow_free:
            later = slice(no_snow_start + run_days, n)

            if not bool(valid[later].all()):
                continue

            later_snow = snow[later]

            if allow_late_transient_snow:
                later_has_sustained_snow = any(
                    bool(later_snow[j : j + run_days].all())
                    for j in range(0, max(0, later_snow.size - run_days + 1))
                )
                if later_has_sustained_snow:
                    continue
            else:
                if bool(later_snow.any()):
                    continue

        return _day_of_water_year(
            sdd_date,
            water_year=water_year,
            one_based=one_based_day,
        )

    return np.nan


def _assert_member_dimension_if_requested(
    ds: xr.Dataset,
    *,
    require_member_dim: bool,
) -> None:
    """
    Ensure ensemble members are represented explicitly when requested.

    Metrics are computed separately for each member only when ``scalarSWE`` has
    a ``member`` dimension.
    """
    if not require_member_dim:
        return

    if VAR_SWE not in ds:
        raise KeyError(f"Merged dataset does not contain {VAR_SWE!r}")

    if "member" not in ds[VAR_SWE].dims:
        raise ValueError(
            "Expected ensemble members, but merged scalarSWE has no 'member' dimension. "
            "For files like chena_002_timestep.nc, use deterministic=False, "
            "force_member_dim=True, multi_gru_per_file=True, and require_member_dim=True."
        )

    if ds.sizes.get("member", 0) < 1:
        raise ValueError("Merged dataset has a member dimension but no members.")


def compute_water_year_metrics(
    swe_merged: xr.Dataset,
    *,
    swe_threshold: float = 1.0,
    run_days: int = 5,
    expected_timesteps_per_day: int = 24,
    min_daily_timesteps: Optional[int] = None,
    require_complete_wy: bool = True,
    min_valid_wy_days: int = 360,
    require_terminal_snow_free: bool = False,
    allow_late_transient_snow: bool = True,
    one_based_sdd: bool = True,
) -> xr.Dataset:
    """
    Compute water-year SCF and SDD from merged timestep SWE.

    SDD is calculated separately for every non-time coordinate combination,
    including every ensemble member when a ``member`` dimension is present.

    Parameters
    ----------
    swe_merged
        Dataset containing ``scalarSWE``.

    swe_threshold
        A day is snow-covered when daily maximum SWE is strictly greater than
        this value.

    run_days
        Required length of consecutive snow and no-snow runs.

    expected_timesteps_per_day
        Expected timestep count per day. Use 24 for hourly SUMMA output.

    min_daily_timesteps
        Minimum number of valid timestep values required for a day to be usable.
        If None, defaults to ``expected_timesteps_per_day``.

    require_complete_wy
        If True, SDD is only computed when the water year has at least
        ``min_valid_wy_days`` valid daily values.

    min_valid_wy_days
        Minimum valid daily values required in a water year for SDD.

    require_terminal_snow_free
        If True, require the selected transition to lead into the terminal
        snow-free part of the water year.

    allow_late_transient_snow
        Only used when ``require_terminal_snow_free=True``. If True, isolated
        later snow days are tolerated, but later sustained snow runs reject the
        candidate.

    one_based_sdd
        If True, Oct 1 = day 1. If False, Oct 1 = day 0.
    """
    if run_days < 1:
        raise ValueError("run_days must be >= 1")

    if expected_timesteps_per_day < 1:
        raise ValueError("expected_timesteps_per_day must be >= 1")

    if min_daily_timesteps is None:
        min_daily_timesteps = expected_timesteps_per_day

    if min_daily_timesteps < 1:
        raise ValueError("min_daily_timesteps must be >= 1")

    if min_valid_wy_days < 1:
        raise ValueError("min_valid_wy_days must be >= 1")

    if VAR_SWE not in swe_merged:
        raise KeyError(f"Input dataset must contain {VAR_SWE!r}")

    swe = swe_merged[VAR_SWE]

    swe_daily, daily_valid_count, valid_daily = _daily_swe_summary(
        swe,
        expected_timesteps_per_day=expected_timesteps_per_day,
        min_daily_timesteps=min_daily_timesteps,
    )

    # Do not use has_snow alone to infer no-snow, because NaN > threshold is False.
    # valid_daily must always be used to distinguish observed no-snow from missing data.
    has_snow = swe_daily > float(swe_threshold)

    wy = us_water_year_label(swe_daily["time"])

    has_snow_wy = has_snow.assign_coords(water_year=("time", wy.data))
    valid_wy = valid_daily.assign_coords(water_year=("time", wy.data))

    # -------------------------------------------------------------------------
    # SCF: fraction of valid days in each water year with daily max SWE > threshold
    # -------------------------------------------------------------------------
    observed_snow_days = (
        (has_snow_wy & valid_wy)
        .groupby("water_year")
        .sum(dim="time")
        .astype(float)
    )

    observed_valid_days = (
        valid_wy
        .groupby("water_year")
        .sum(dim="time")
        .astype(float)
    )

    scf_fraction = xr.where(
        observed_valid_days > 0,
        observed_snow_days / observed_valid_days,
        np.nan,
    ).rename(VAR_SCF_DAYS)

    scf_fraction.attrs.clear()
    scf_fraction.attrs["long_name"] = (
        f"Snow cover frequency: fraction of valid days with daily max {VAR_SWE} > {swe_threshold}"
    )
    scf_fraction.attrs["units"] = "1"
    scf_fraction.attrs["description"] = (
        "Fraction of valid days per US water year, Oct 1 to Sep 30, with snow present. "
        f"Hourly/subdaily {VAR_SWE} is collapsed to daily maximum before thresholding. "
        f"A day is valid only if it has at least {min_daily_timesteps} non-missing timestep values. "
        "Missing or insufficiently sampled days are excluded from the denominator."
    )

    valid_days_per_wy = observed_valid_days.rename(VAR_VALID_DAYS)
    valid_days_per_wy.attrs.clear()
    valid_days_per_wy.attrs["long_name"] = "Number of valid daily SWE values per water year"
    valid_days_per_wy.attrs["units"] = "days"
    valid_days_per_wy.attrs["description"] = (
        f"Number of daily values with at least {min_daily_timesteps} valid timestep samples."
    )

    # -------------------------------------------------------------------------
    # SDD: numeric day of water year, separately for each member/GRU combination
    # -------------------------------------------------------------------------
    non_time_dims = [d for d in has_snow.dims if d != "time"]

    wy_values = np.asarray(wy.values)
    times_np = pd.to_datetime(swe_daily["time"].values)
    uniq_wy = np.unique(wy_values).astype(np.int32)
    n_wy = int(uniq_wy.size)

    if non_time_dims:
        stacked_snow = has_snow.stack(stacked=tuple(non_time_dims))
        stacked_valid = valid_daily.stack(stacked=tuple(non_time_dims))

        n_loc = int(stacked_snow.sizes["stacked"])

        sdd_flat = np.full(
            (n_loc, n_wy),
            np.nan,
            dtype=np.float32,
        )

        snow_values = stacked_snow.values
        valid_values = stacked_valid.values

        for idx in range(n_loc):
            isnow = snow_values[:, idx]
            ivalid = valid_values[:, idx]

            for j, wyi in enumerate(uniq_wy):
                sel = wy_values == wyi
                if not np.any(sel):
                    continue

                order = np.argsort(times_np[sel])

                t_sub = times_np[sel][order]
                snow_sub = isnow[sel][order]
                valid_sub = ivalid[sel][order]

                sdd_flat[idx, j] = _snow_disappearance_day_for_year(
                    snow_sub,
                    valid_sub,
                    t_sub,
                    water_year=int(wyi),
                    run_days=run_days,
                    require_complete_wy=require_complete_wy,
                    min_valid_wy_days=min_valid_wy_days,
                    require_terminal_snow_free=require_terminal_snow_free,
                    allow_late_transient_snow=allow_late_transient_snow,
                    one_based_day=one_based_sdd,
                )

        sdd = (
            xr.DataArray(
                sdd_flat,
                dims=("stacked", "water_year"),
                coords={
                    "stacked": stacked_snow["stacked"],
                    "water_year": uniq_wy,
                },
                name=VAR_SNOW_DISAPPEARANCE,
            )
            .unstack("stacked")
        )

    else:
        sdd_values = np.full(n_wy, np.nan, dtype=np.float32)

        snow_values = has_snow.values
        valid_values = valid_daily.values

        for j, wyi in enumerate(uniq_wy):
            sel = wy_values == wyi
            if not np.any(sel):
                continue

            order = np.argsort(times_np[sel])

            t_sub = times_np[sel][order]
            snow_sub = snow_values[sel][order]
            valid_sub = valid_values[sel][order]

            sdd_values[j] = _snow_disappearance_day_for_year(
                snow_sub,
                valid_sub,
                t_sub,
                water_year=int(wyi),
                run_days=run_days,
                require_complete_wy=require_complete_wy,
                min_valid_wy_days=min_valid_wy_days,
                require_terminal_snow_free=require_terminal_snow_free,
                allow_late_transient_snow=allow_late_transient_snow,
                one_based_day=one_based_sdd,
            )

        sdd = xr.DataArray(
            sdd_values,
            dims=("water_year",),
            coords={"water_year": uniq_wy},
            name=VAR_SNOW_DISAPPEARANCE,
        )

    sdd.attrs.clear()
    sdd.attrs["long_name"] = "Snow disappearance day of water year"
    sdd.attrs["units"] = "day"
    sdd.attrs["description"] = (
        "Snow disappearance date expressed as numeric day of water year, not as datetime. "
        f"Oct 1 is {'day 1' if one_based_sdd else 'day 0'}. "
        f"Hourly/subdaily {VAR_SWE} is first aggregated to daily maximum. "
        f"A valid snow day requires daily max {VAR_SWE} > {swe_threshold}; "
        f"a valid snow-free day requires daily max {VAR_SWE} <= {swe_threshold}. "
        f"Each daily value must have at least {min_daily_timesteps} valid timestep samples. "
        f"The SDD rule searches backward within each water year for {run_days} consecutive "
        f"valid snow days immediately followed by {run_days} consecutive valid snow-free days. "
        "The reported value is the first day of the qualifying snow-free run. "
        "Missing or undetected SDD is stored as NaN."
    )

    spatial_dims = tuple(d for d in sdd.dims if d != "water_year")
    dim_order = ("water_year",) + spatial_dims

    out = xr.Dataset(
        {
            VAR_SCF_DAYS: scf_fraction.transpose(*dim_order),
            VAR_SNOW_DISAPPEARANCE: sdd.transpose(*dim_order),
            VAR_VALID_DAYS: valid_days_per_wy.transpose(*dim_order),
        }
    )

    return out


def attach_cf_metadata(
    ds: xr.Dataset,
    *,
    title: str = "Water-year snow metrics from SUMMA",
    source: str = "SUMMA timestep output",
    history_note: str = "",
    institution: str = "",
    swe_threshold: float = 1.0,
    run_days: int = 5,
    expected_timesteps_per_day: int = 24,
    min_daily_timesteps: Optional[int] = None,
    require_complete_wy: bool = True,
    min_valid_wy_days: int = 360,
    require_terminal_snow_free: bool = False,
    allow_late_transient_snow: bool = True,
    one_based_sdd: bool = True,
) -> xr.Dataset:
    """Attach CF-style metadata to output dataset."""
    ds = ds.copy()

    if min_daily_timesteps is None:
        min_daily_timesteps = expected_timesteps_per_day

    hist = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    if history_note:
        hist = f"{hist} {history_note}"

    ds.attrs["Conventions"] = "CF-1.8"
    ds.attrs["title"] = title
    ds.attrs["source"] = source
    ds.attrs["history"] = hist

    if institution:
        ds.attrs["institution"] = institution

    ds.attrs["water_year_definition"] = (
        "US water year: October 1 00:00 through September 30, labeled by the calendar "
        "year in which the period ends."
    )

    ds.attrs["snow_metrics_parameters"] = (
        f"swe_threshold={swe_threshold} using strict daily max {VAR_SWE} > threshold for snow; "
        f"run_days={run_days}; "
        f"expected_timesteps_per_day={expected_timesteps_per_day}; "
        f"min_daily_timesteps={min_daily_timesteps}; "
        f"require_complete_wy={require_complete_wy}; "
        f"min_valid_wy_days={min_valid_wy_days}; "
        f"require_terminal_snow_free={require_terminal_snow_free}; "
        f"allow_late_transient_snow={allow_late_transient_snow}; "
        f"one_based_sdd={one_based_sdd}."
    )

    ds.attrs["daily_aggregation"] = (
        f"Input timestep SWE is aggregated to daily maximum before computing snow metrics. "
        f"A day is valid only when at least {min_daily_timesteps} non-missing timestep values "
        f"are available. Missing or insufficiently sampled days are not treated as snow-free."
    )

    ds.attrs["sdd_encoding"] = (
        f"{VAR_SNOW_DISAPPEARANCE} is stored as a numeric day-of-water-year value. "
        f"Oct 1 is {'1' if one_based_sdd else '0'}. Missing values are NaN."
    )

    if VAR_SCF_DAYS in ds:
        ds[VAR_SCF_DAYS].attrs.setdefault("units", "1")

    if VAR_SNOW_DISAPPEARANCE in ds:
        ds[VAR_SNOW_DISAPPEARANCE].attrs.setdefault("units", "day")

    if VAR_VALID_DAYS in ds:
        ds[VAR_VALID_DAYS].attrs.setdefault("units", "days")

    return ds


def build_netcdf_encoding(
    ds: xr.Dataset,
    *,
    compress: bool = False,
) -> Dict[str, Dict[str, Any]]:
    """
    Build NetCDF encoding.

    SDD is float32 with NaN missing values, avoiding datetime NaT integer artifacts.
    """
    encoding: Dict[str, Dict[str, Any]] = {}

    for name, da in ds.data_vars.items():
        enc: Dict[str, Any] = {}

        if compress:
            enc["zlib"] = True
            enc["complevel"] = 4

        if name in (VAR_SCF_DAYS, VAR_SNOW_DISAPPEARANCE, VAR_VALID_DAYS):
            enc["dtype"] = "float32"
            enc["_FillValue"] = np.float32(np.nan)
        elif np.issubdtype(da.dtype, np.floating):
            enc["dtype"] = "float32"
            enc["_FillValue"] = np.float32(np.nan)

        encoding[name] = enc

    return encoding


def write_water_year_snow_metrics_netcdf(
    output_path: Union[str, Path],
    *,
    input_dir: Union[str, Path],
    glob_pattern: str = "**/*.nc",
    case_name_filter: Optional[str] = None,
    deterministic: bool = False,
    time_join: str = "outer",
    force_member_dim: bool = True,
    require_member_dim: bool = False,
    prefer_multi_gru_files: bool = True,
    multi_gru_per_file: Optional[bool] = None,
    swe_threshold: float = 1.0,
    run_days: int = 5,
    expected_timesteps_per_day: int = 24,
    min_daily_timesteps: Optional[int] = None,
    require_complete_wy: bool = True,
    min_valid_wy_days: int = 360,
    require_terminal_snow_free: bool = False,
    allow_late_transient_snow: bool = True,
    one_based_sdd: bool = True,
    compress: bool = False,
    title: str = "Water-year snow metrics from SUMMA",
    source: str = "SUMMA timestep output",
    history_note: str = "",
    institution: str = "",
) -> xr.Dataset:
    """
    Discover timestep files, merge SWE, compute water-year metrics, and write NetCDF.

    For ensemble files such as:

        chena_001_timestep.nc
        chena_002_timestep.nc
        chena_003_timestep.nc

    recommended settings are:

        deterministic=False
        force_member_dim=True
        require_member_dim=True
        multi_gru_per_file=True

    SDD is written as numeric day of water year, not datetime.
    """
    if min_daily_timesteps is None:
        min_daily_timesteps = expected_timesteps_per_day

    paths = discover_nc_files(
        input_dir,
        glob_pattern=glob_pattern,
        case_name_filter=case_name_filter,
        prefer_multi_gru_files=prefer_multi_gru_files,
    )

    if not paths:
        raise FileNotFoundError(
            f"No NetCDF files matched under {input_dir!r} with glob {glob_pattern!r}"
        )

    if deterministic and require_member_dim:
        raise ValueError(
            "deterministic=True removes ensemble member inference. "
            "For ensemble files like chena_002_timestep.nc, use deterministic=False."
        )

    swe_ds = merge_timestep_swe_files(
        paths,
        deterministic=deterministic,
        time_join=time_join,
        force_member_dim=force_member_dim,
        multi_gru_per_file=multi_gru_per_file,
    )

    logger.info("Merged SWE dims: %s", dict(swe_ds[VAR_SWE].sizes))
    logger.info("Merged SWE dimension order: %s", swe_ds[VAR_SWE].dims)

    if "member" in swe_ds[VAR_SWE].dims:
        logger.info("Members found: %s", swe_ds["member"].values)
    else:
        logger.warning(
            "Merged SWE has no member dimension. Metrics will not be separated by ensemble member."
        )

    _assert_member_dimension_if_requested(
        swe_ds,
        require_member_dim=require_member_dim,
    )

    ds = compute_water_year_metrics(
        swe_ds,
        swe_threshold=swe_threshold,
        run_days=run_days,
        expected_timesteps_per_day=expected_timesteps_per_day,
        min_daily_timesteps=min_daily_timesteps,
        require_complete_wy=require_complete_wy,
        min_valid_wy_days=min_valid_wy_days,
        require_terminal_snow_free=require_terminal_snow_free,
        allow_late_transient_snow=allow_late_transient_snow,
        one_based_sdd=one_based_sdd,
    )

    ds = attach_cf_metadata(
        ds,
        title=title,
        source=source,
        history_note=history_note,
        institution=institution,
        swe_threshold=swe_threshold,
        run_days=run_days,
        expected_timesteps_per_day=expected_timesteps_per_day,
        min_daily_timesteps=min_daily_timesteps,
        require_complete_wy=require_complete_wy,
        min_valid_wy_days=min_valid_wy_days,
        require_terminal_snow_free=require_terminal_snow_free,
        allow_late_transient_snow=allow_late_transient_snow,
        one_based_sdd=one_based_sdd,
    )

    out_path = Path(output_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    encoding = build_netcdf_encoding(ds, compress=compress)
    ds.to_netcdf(out_path, encoding=encoding)

    logger.info("Wrote %s", out_path)

    return ds


def debug_water_year_snow_series(
    swe_merged: xr.Dataset,
    *,
    water_year: int,
    gru: Optional[Union[int, str]] = None,
    member: Optional[Union[int, str]] = None,
    swe_threshold: float = 1.0,
    expected_timesteps_per_day: int = 24,
    min_daily_timesteps: Optional[int] = None,
) -> pd.DataFrame:
    """
    Return a daily diagnostic table for one member/GRU/water year.

    This is useful for checking why SDD is early, missing, or unexpected.

    Returns columns:

        time
        water_year
        daily_max_swe
        valid_timestep_count
        valid_day
        has_snow
        day_of_water_year
    """
    if min_daily_timesteps is None:
        min_daily_timesteps = expected_timesteps_per_day

    if VAR_SWE not in swe_merged:
        raise KeyError(f"Input dataset must contain {VAR_SWE!r}")

    swe = swe_merged[VAR_SWE]

    if member is not None and "member" in swe.dims:
        swe = swe.sel(member=member)

    if gru is not None and "gru" in swe.dims:
        swe = swe.sel(gru=gru)

    swe_daily, valid_count, valid_day = _daily_swe_summary(
        swe,
        expected_timesteps_per_day=expected_timesteps_per_day,
        min_daily_timesteps=min_daily_timesteps,
    )

    has_snow = swe_daily > float(swe_threshold)
    wy = us_water_year_label(swe_daily["time"])

    sel = np.asarray(wy.values) == int(water_year)

    times = pd.to_datetime(swe_daily["time"].values[sel])
    wy_start = _water_year_start_date(water_year)

    day_of_wy = np.array(
        [
            int((np.datetime64(t, "D") - wy_start) / np.timedelta64(1, "D")) + 1
            for t in times
        ],
        dtype=int,
    )

    return pd.DataFrame(
        {
            "time": times,
            "water_year": np.asarray(wy.values)[sel],
            "daily_max_swe": np.asarray(swe_daily.values)[sel],
            "valid_timestep_count": np.asarray(valid_count.values)[sel],
            "valid_day": np.asarray(valid_day.values)[sel],
            "has_snow": np.asarray(has_snow.values)[sel],
            "day_of_water_year": day_of_wy,
        }
    )