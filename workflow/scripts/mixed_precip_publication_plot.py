"""Publication-quality snow-fraction climatology panels by water-year day."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.patches import Patch

import compute_mixed_precip_fractions as mpf

DEFAULT_SOURCE_COLORS = ("#2166ac", "#b2182b", "#4daf4a", "#984ea3")
ENSEMBLE_BAND_ALPHA = 0.30

_WYDOY_TICKS = (1, 92, 183, 274, 365)
_WYDOY_LABELS = ("Oct 1", "Jan 1", "Apr 1", "Jul 1", "Sep 30")

WyDoyClimatology = dict[str, dict[str, dict[str, np.ndarray | None]]]
SpatialGruCache = dict[str, dict[str, dict[str, np.ndarray]]]


def _apply_publication_style() -> None:
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.labelsize": 10,
            "axes.titlesize": 11,
            "axes.titleweight": "normal",
            "axes.titlepad": 8,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
            "axes.linewidth": 0.8,
            "grid.linewidth": 0.5,
            "grid.alpha": 0.35,
        }
    )


def _format_catchment_label(catchment: str) -> str:
    return catchment.strip().replace("_", " ").title()


def _labels_consistent_across_catchments(
    output_nc_by_catchment: dict[str, dict[str, Path]],
) -> list[str]:
    labels: list[str] | None = None
    for catchment, path_dict in output_nc_by_catchment.items():
        keys = list(path_dict.keys())
        if labels is None:
            labels = keys
        elif keys != labels:
            raise ValueError(
                f"Catchment {catchment!r} has sources {keys!r}; expected {labels!r} "
                "for all catchments."
            )
    if labels is None:
        raise ValueError("output_nc_by_catchment is empty.")
    return labels


def _water_year_doy_vectorized(time_vals: np.ndarray) -> np.ndarray:
    """Water-year day (1–365) for each timestamp; vectorized over ``time_vals``."""
    ts = pd.DatetimeIndex(time_vals)
    years = ts.year.to_numpy()
    oct1_year = np.where(ts.month.to_numpy() >= 10, years, years - 1)
    oct1 = pd.DatetimeIndex(pd.to_datetime({"year": oct1_year, "month": 10, "day": 1}))
    return (ts - oct1).days + 1


def _leap_day_mask(time_vals: np.ndarray) -> np.ndarray:
    ts = pd.DatetimeIndex(time_vals)
    return (ts.month.to_numpy() == 2) & (ts.day.to_numpy() == 29)


def _mean_by_wydoy(values: np.ndarray, wydoy: np.ndarray) -> np.ndarray:
    """
    Mean of ``values`` grouped by water-year day.

    ``values`` shape ``(n_times,)`` or ``(n_times, n_members)`` → ``(365,)`` or ``(365, n_members)``.
    """
    wydoy = np.asarray(wydoy, dtype=np.int32)
    valid = (wydoy >= 1) & (wydoy <= 365)
    wydoy = wydoy[valid]
    values = values[valid] if values.ndim == 1 else values[valid, :]

    if values.ndim == 1:
        finite = np.isfinite(values)
        if not finite.any():
            return np.full(365, np.nan)
        w = wydoy[finite]
        v = values[finite]
        sums = np.bincount(w, weights=v, minlength=366)[1:]
        counts = np.bincount(w, minlength=366)[1:]
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(counts > 0, sums / counts, np.nan)

    _, n_members = values.shape
    out = np.full((365, n_members), np.nan)
    w_idx = wydoy - 1
    for j in range(n_members):
        col = values[:, j]
        finite = np.isfinite(col)
        if not finite.any():
            continue
        w = w_idx[finite]
        v = col[finite]
        sums = np.bincount(w, weights=v, minlength=365)
        counts = np.bincount(w, minlength=365)
        with np.errstate(invalid="ignore", divide="ignore"):
            out[:, j] = np.where(counts > 0, sums / counts, np.nan)
    return out


def _doy_composite_by_member(
    da: xr.DataArray,
    *,
    resample_for_plot: str | None = None,
    handle_leap_days: str = "drop",
) -> np.ndarray | None:
    """Per-member snow-fraction composite by water-year day; shape (365, n_members) or None."""
    if resample_for_plot:
        da = da.resample(time=resample_for_plot).mean(skipna=True)
    if "member" not in da.dims:
        return None
    da = da.transpose("time", "member")

    time_vals = da["time"].values
    all_wydoy = _water_year_doy_vectorized(time_vals)
    mask = np.ones(len(time_vals), dtype=bool)
    if handle_leap_days == "drop":
        mask = ~_leap_day_mask(time_vals)

    values = da.values
    if values.ndim == 1:
        values = values[:, None]
    return _mean_by_wydoy(values[mask, :], all_wydoy[mask])


def _mean_per_wydoy(
    da: xr.DataArray,
    *,
    resample_for_plot: str | None = None,
    handle_leap_days: str = "drop",
) -> tuple[np.ndarray, np.ndarray | None, np.ndarray | None, np.ndarray | None]:
    """
    Composite snow fraction by water-year day (1–365), pooled over calendar years.

    Returns ``(mean, p10, p90, member_std)``; band/std entries are ``None`` for deterministic runs.
    """
    per_member = _doy_composite_by_member(
        da,
        resample_for_plot=resample_for_plot,
        handle_leap_days=handle_leap_days,
    )
    if per_member is not None:
        ens_mean = np.nanmean(per_member, axis=1)
        ens_p10 = np.nanpercentile(per_member, 10, axis=1)
        ens_p90 = np.nanpercentile(per_member, 90, axis=1)
        ens_std = np.nanstd(per_member, axis=1)
        return ens_mean, ens_p10, ens_p90, ens_std

    if resample_for_plot:
        da = da.resample(time=resample_for_plot).mean(skipna=True)

    time_vals = da["time"].values
    all_wydoy = _water_year_doy_vectorized(time_vals)
    mask = np.ones(len(time_vals), dtype=bool)
    if handle_leap_days == "drop":
        mask = ~_leap_day_mask(time_vals)

    values = da.values[mask]
    mean_per_doy = _mean_by_wydoy(values, all_wydoy[mask])
    return mean_per_doy, None, None, None


def _manifest_fingerprint(
    output_nc_by_catchment: dict[str, dict[str, Path]],
    catchments_list: list[str],
    labels_ordered: list[str],
) -> str:
    entries: list[tuple[str, str, int, int]] = []
    for catchment in catchments_list:
        for label in labels_ordered:
            path = Path(output_nc_by_catchment[catchment][label])
            stat = path.stat()
            entries.append((catchment, label, stat.st_mtime_ns, stat.st_size))
    payload = json.dumps(entries, sort_keys=True)
    return hashlib.sha256(payload.encode()).hexdigest()[:16]


def _spatial_cache_max_gru(
    cache: SpatialGruCache,
    catchments: list[str],
    labels: list[str],
) -> int:
    max_gru = 0
    for catchment in catchments:
        for label in labels:
            entry = cache.get(catchment, {}).get(label)
            if entry is not None:
                max_gru = max(max_gru, int(entry["gru"].size))
    return max_gru


def _spatial_reference_gru(
    cache: SpatialGruCache,
    catchment: str,
    labels: list[str],
) -> np.ndarray:
    for label in labels:
        entry = cache.get(catchment, {}).get(label)
        if entry is not None:
            return entry["gru"]
    raise ValueError(f"No spatial ensemble data cached for catchment {catchment!r}.")


def _spatial_cache_to_arrays(
    cache: SpatialGruCache,
    catchments: list[str],
    labels: list[str],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    max_gru = _spatial_cache_max_gru(cache, catchments, labels)
    gru_id = np.full((len(catchments), max_gru), -1, dtype=np.int64)
    values = np.full((len(catchments), len(labels), max_gru), np.nan, dtype=np.float32)
    n_gru = np.zeros(len(catchments), dtype=np.int32)
    for i, catchment in enumerate(catchments):
        gru = _spatial_reference_gru(cache, catchment, labels)
        n = int(gru.size)
        n_gru[i] = n
        gru_id[i, :n] = gru
        for j, label in enumerate(labels):
            entry = cache.get(catchment, {}).get(label)
            if entry is not None:
                values[i, j, :n] = entry["value"]
    return gru_id, values, n_gru


def _spatial_arrays_to_cache(
    gru_id: np.ndarray,
    values: np.ndarray,
    n_gru: np.ndarray,
    catchments: list[str],
    labels: list[str],
) -> SpatialGruCache:
    cache: SpatialGruCache = {}
    for i, catchment in enumerate(catchments):
        cache[catchment] = {}
        n = int(n_gru[i])
        gru = np.asarray(gru_id[i, :n], dtype=np.int64)
        for j, label in enumerate(labels):
            vals = np.asarray(values[i, j, :n], dtype=float)
            if np.all(np.isnan(vals)):
                continue
            cache[catchment][label] = {
                "gru": gru,
                "value": vals,
            }
    return cache


def _year_range_attr(start_year: int | None, end_year: int | None) -> tuple[str, str]:
    return (
        "" if start_year is None else str(start_year),
        "" if end_year is None else str(end_year),
    )


def compute_snow_fraction_wydoy_climatology(
    output_nc_by_catchment: dict[str, dict[str, Path]],
    *,
    catchments_list: list[str] | None = None,
    resample_for_plot: str | None = None,
    handle_leap_days: str = "drop",
    wy_snowfall_start_year: int | None = None,
    wy_snowfall_end_year: int | None = None,
) -> WyDoyClimatology:
    """Compute snow-fraction water-year-day climatology for all catchments and sources."""
    climatology, _, _ = compute_mixed_precip_publication_cache(
        output_nc_by_catchment,
        catchments_list=catchments_list,
        resample_for_plot=resample_for_plot,
        handle_leap_days=handle_leap_days,
        wy_snowfall_start_year=wy_snowfall_start_year,
        wy_snowfall_end_year=wy_snowfall_end_year,
    )
    return climatology


def compute_mixed_precip_publication_cache(
    output_nc_by_catchment: dict[str, dict[str, Path]],
    *,
    catchments_list: list[str] | None = None,
    resample_for_plot: str | None = None,
    handle_leap_days: str = "drop",
    wy_snowfall_start_year: int | None = None,
    wy_snowfall_end_year: int | None = None,
) -> tuple[WyDoyClimatology, SpatialGruCache, SpatialGruCache]:
    """
    Compute wydoy climatology and spatial member-std caches in one pass over source NetCDFs.

    Returns ``(wydoy_climatology, spatial_member_std, spatial_wy_snowfall_member_std)``.
    """
    from mixed_precip_spatial_plot import (
        spatial_mean_wy_snowfall_member_std_per_gru,
        spatial_member_std_per_gru,
    )

    catchments_list = catchments_list or list(output_nc_by_catchment.keys())
    labels_ordered = _labels_consistent_across_catchments(output_nc_by_catchment)
    climatology: WyDoyClimatology = {}
    spatial_member_std: SpatialGruCache = {}
    spatial_wy_snowfall: SpatialGruCache = {}

    for catchment in catchments_list:
        path_dict = output_nc_by_catchment[catchment]
        climatology[catchment] = {}
        spatial_member_std[catchment] = {}
        spatial_wy_snowfall[catchment] = {}
        opened: dict[str, xr.Dataset] = {}
        try:
            for label in labels_ordered:
                opened[label] = xr.open_dataset(path_dict[label], decode_times=True)
            for label in labels_ordered:
                ds = opened[label]
                da = ds[mpf.VAR_FRAC].mean(dim="gru", skipna=True)
                mean, p10, p90, member_std = _mean_per_wydoy(
                    da,
                    resample_for_plot=resample_for_plot,
                    handle_leap_days=handle_leap_days,
                )
                climatology[catchment][label] = {
                    "mean": mean,
                    "p10": p10,
                    "p90": p90,
                    "member_std": member_std,
                }

                if "member" in ds[mpf.VAR_FRAC].dims:
                    frac_std_da = spatial_member_std_per_gru(ds, mpf.VAR_FRAC)
                    spatial_member_std[catchment][label] = {
                        "gru": np.asarray(frac_std_da["gru"].values, dtype=np.int64),
                        "value": np.asarray(frac_std_da.values, dtype=float),
                    }

                if "member" in ds[mpf.VAR_SNOW].dims:
                    wy_std_da = spatial_mean_wy_snowfall_member_std_per_gru(
                        ds,
                        start_year=wy_snowfall_start_year,
                        end_year=wy_snowfall_end_year,
                        var=mpf.VAR_SNOW,
                    )
                    spatial_wy_snowfall[catchment][label] = {
                        "gru": np.asarray(wy_std_da["gru"].values, dtype=np.int64),
                        "value": np.asarray(wy_std_da.values, dtype=float),
                    }
        finally:
            for ds in opened.values():
                ds.close()

    return climatology, spatial_member_std, spatial_wy_snowfall


def save_snow_fraction_wydoy_climatology(
    climatology: WyDoyClimatology,
    cache_path: Path | str,
    *,
    output_nc_by_catchment: dict[str, dict[str, Path]] | None = None,
    catchments_list: list[str] | None = None,
    resample_for_plot: str | None = None,
    handle_leap_days: str = "drop",
    spatial_member_std: SpatialGruCache | None = None,
    spatial_wy_snowfall_member_std: SpatialGruCache | None = None,
    wy_snowfall_start_year: int | None = None,
    wy_snowfall_end_year: int | None = None,
) -> Path:
    """Write climatology and optional spatial caches to NetCDF for fast reload during plotting."""
    cache_path = Path(cache_path)
    cache_path.parent.mkdir(parents=True, exist_ok=True)

    catchments = catchments_list or list(climatology.keys())
    labels = _labels_consistent_across_catchments(
        {c: {lbl: Path(".") for lbl in climatology[c]} for c in catchments}
    )
    wydoy = np.arange(1, 366, dtype=np.int16)

    mean_arr = np.full((len(catchments), len(labels), 365), np.nan, dtype=np.float32)
    p10_arr = mean_arr.copy()
    p90_arr = mean_arr.copy()
    std_arr = mean_arr.copy()

    for i, catchment in enumerate(catchments):
        for j, label in enumerate(labels):
            stats = climatology[catchment][label]
            mean_arr[i, j, :] = stats["mean"]
            if stats["p10"] is not None:
                p10_arr[i, j, :] = stats["p10"]
            if stats["p90"] is not None:
                p90_arr[i, j, :] = stats["p90"]
            if stats["member_std"] is not None:
                std_arr[i, j, :] = stats["member_std"]

    attrs: dict[str, str] = {
        "handle_leap_days": handle_leap_days,
        "resample_for_plot": resample_for_plot or "",
        "title": "Mixed precip publication cache (wydoy climatology and spatial member std)",
    }
    wy_start_attr, wy_end_attr = _year_range_attr(wy_snowfall_start_year, wy_snowfall_end_year)
    attrs["wy_snowfall_start_year"] = wy_start_attr
    attrs["wy_snowfall_end_year"] = wy_end_attr
    if output_nc_by_catchment is not None:
        attrs["source_fingerprint"] = _manifest_fingerprint(
            output_nc_by_catchment, catchments, labels
        )

    data_vars = {
        "mean": (["catchment", "source", "wydoy"], mean_arr),
        "p10": (["catchment", "source", "wydoy"], p10_arr),
        "p90": (["catchment", "source", "wydoy"], p90_arr),
        "member_std": (["catchment", "source", "wydoy"], std_arr),
    }
    if spatial_member_std is not None:
        spatial_gru_id, spatial_values, spatial_n_gru = _spatial_cache_to_arrays(
            spatial_member_std, catchments, labels
        )
        data_vars["spatial_gru_id"] = (["catchment", "gru"], spatial_gru_id)
        data_vars["spatial_n_gru"] = (["catchment"], spatial_n_gru)
        data_vars["spatial_member_std"] = (
            ["catchment", "source", "gru"],
            spatial_values,
        )
    if spatial_wy_snowfall_member_std is not None:
        wy_gru_id, wy_values, wy_n_gru = _spatial_cache_to_arrays(
            spatial_wy_snowfall_member_std, catchments, labels
        )
        data_vars["wy_snowfall_spatial_gru_id"] = (["catchment", "gru"], wy_gru_id)
        data_vars["wy_snowfall_spatial_n_gru"] = (["catchment"], wy_n_gru)
        data_vars["wy_snowfall_spatial_member_std"] = (
            ["catchment", "source", "gru"],
            wy_values,
        )

    ds = xr.Dataset(
        data_vars=data_vars,
        coords={
            "catchment": catchments,
            "source": labels,
            "wydoy": wydoy,
        },
        attrs=attrs,
    )
    ds.to_netcdf(cache_path)
    return cache_path


def load_snow_fraction_wydoy_climatology(cache_path: Path | str) -> WyDoyClimatology:
    """Load climatology written by ``save_snow_fraction_wydoy_climatology``."""
    with xr.open_dataset(cache_path) as ds:
        catchments = [str(c) for c in ds["catchment"].values]
        labels = [str(s) for s in ds["source"].values]
        climatology: WyDoyClimatology = {}
        for i, catchment in enumerate(catchments):
            climatology[catchment] = {}
            for j, label in enumerate(labels):
                p10 = ds["p10"].isel(catchment=i, source=j).values
                p90 = ds["p90"].isel(catchment=i, source=j).values
                std = ds["member_std"].isel(catchment=i, source=j).values
                climatology[catchment][label] = {
                    "mean": ds["mean"].isel(catchment=i, source=j).values,
                    "p10": None if np.all(np.isnan(p10)) else p10,
                    "p90": None if np.all(np.isnan(p90)) else p90,
                    "member_std": None if np.all(np.isnan(std)) else std,
                }
    return climatology


def load_spatial_member_std_cache(cache_path: Path | str) -> SpatialGruCache:
    """Load spatial snow-fraction member std from the publication cache NetCDF."""
    with xr.open_dataset(cache_path) as ds:
        if "spatial_member_std" not in ds:
            raise KeyError(
                f"{cache_path} has no spatial_member_std variable; "
                "recompute with ensure_snow_fraction_wydoy_climatology(..., overwrite=True)."
            )
        catchments = [str(c) for c in ds["catchment"].values]
        labels = [str(s) for s in ds["source"].values]
        return _spatial_arrays_to_cache(
            ds["spatial_gru_id"].values,
            ds["spatial_member_std"].values,
            ds["spatial_n_gru"].values,
            catchments,
            labels,
        )


def load_spatial_wy_snowfall_member_std_cache(
    cache_path: Path | str,
    *,
    start_year: int | None = None,
    end_year: int | None = None,
) -> SpatialGruCache:
    """Load spatial mean WY snowfall member std from the publication cache NetCDF."""
    with xr.open_dataset(cache_path) as ds:
        if "wy_snowfall_spatial_member_std" not in ds:
            raise KeyError(
                f"{cache_path} has no wy_snowfall_spatial_member_std variable; "
                "recompute with ensure_snow_fraction_wydoy_climatology(..., overwrite=True)."
            )
        expected_start, expected_end = _year_range_attr(start_year, end_year)
        if ds.attrs.get("wy_snowfall_start_year", "") != expected_start:
            raise ValueError(
                "Cache wy_snowfall_start_year "
                f"{ds.attrs.get('wy_snowfall_start_year', '')!r} != requested {expected_start!r}."
            )
        if ds.attrs.get("wy_snowfall_end_year", "") != expected_end:
            raise ValueError(
                "Cache wy_snowfall_end_year "
                f"{ds.attrs.get('wy_snowfall_end_year', '')!r} != requested {expected_end!r}."
            )
        catchments = [str(c) for c in ds["catchment"].values]
        labels = [str(s) for s in ds["source"].values]
        return _spatial_arrays_to_cache(
            ds["wy_snowfall_spatial_gru_id"].values,
            ds["wy_snowfall_spatial_member_std"].values,
            ds["wy_snowfall_spatial_n_gru"].values,
            catchments,
            labels,
        )


def wydoy_climatology_cache_is_current(
    cache_path: Path | str,
    output_nc_by_catchment: dict[str, dict[str, Path]],
    *,
    catchments_list: list[str] | None = None,
    resample_for_plot: str | None = None,
    handle_leap_days: str = "drop",
    wy_snowfall_start_year: int | None = None,
    wy_snowfall_end_year: int | None = None,
) -> bool:
    """Return True if ``cache_path`` matches inputs and is newer than source NetCDFs."""
    cache_path = Path(cache_path)
    if not cache_path.is_file():
        return False

    catchments_list = catchments_list or list(output_nc_by_catchment.keys())
    labels_ordered = _labels_consistent_across_catchments(output_nc_by_catchment)
    expected_start, expected_end = _year_range_attr(
        wy_snowfall_start_year, wy_snowfall_end_year
    )

    with xr.open_dataset(cache_path) as ds:
        if ds.attrs.get("handle_leap_days", "drop") != handle_leap_days:
            return False
        if ds.attrs.get("resample_for_plot", "") != (resample_for_plot or ""):
            return False
        if list(ds["catchment"].values) != catchments_list:
            return False
        if list(ds["source"].values) != labels_ordered:
            return False
        if ds.attrs.get("wy_snowfall_start_year", "") != expected_start:
            return False
        if ds.attrs.get("wy_snowfall_end_year", "") != expected_end:
            return False
        if "spatial_member_std" not in ds:
            return False
        if "wy_snowfall_spatial_member_std" not in ds:
            return False
        expected_fp = _manifest_fingerprint(
            output_nc_by_catchment, catchments_list, labels_ordered
        )
        if ds.attrs.get("source_fingerprint") != expected_fp:
            return False

    cache_mtime = cache_path.stat().st_mtime
    for catchment in catchments_list:
        for label in labels_ordered:
            src_mtime = Path(output_nc_by_catchment[catchment][label]).stat().st_mtime
            if src_mtime > cache_mtime:
                return False
    return True


def ensure_snow_fraction_wydoy_climatology(
    output_nc_by_catchment: dict[str, dict[str, Path]],
    cache_path: Path | str,
    *,
    catchments_list: list[str] | None = None,
    resample_for_plot: str | None = None,
    handle_leap_days: str = "drop",
    wy_snowfall_start_year: int | None = None,
    wy_snowfall_end_year: int | None = None,
    overwrite: bool = False,
) -> WyDoyClimatology:
    """
    Load cached climatology when valid; otherwise compute from NetCDFs and save.

    The cache file also stores spatial member-std maps used by the spatial plot helpers.
    Subsequent plot calls can pass the returned dict via ``wydoy_climatology=`` for instant plotting.
    """
    cache_path = Path(cache_path)
    if not overwrite and wydoy_climatology_cache_is_current(
        cache_path,
        output_nc_by_catchment,
        catchments_list=catchments_list,
        resample_for_plot=resample_for_plot,
        handle_leap_days=handle_leap_days,
        wy_snowfall_start_year=wy_snowfall_start_year,
        wy_snowfall_end_year=wy_snowfall_end_year,
    ):
        return load_snow_fraction_wydoy_climatology(cache_path)

    climatology, spatial_member_std, spatial_wy_snowfall = compute_mixed_precip_publication_cache(
        output_nc_by_catchment,
        catchments_list=catchments_list,
        resample_for_plot=resample_for_plot,
        handle_leap_days=handle_leap_days,
        wy_snowfall_start_year=wy_snowfall_start_year,
        wy_snowfall_end_year=wy_snowfall_end_year,
    )
    save_snow_fraction_wydoy_climatology(
        climatology,
        cache_path,
        output_nc_by_catchment=output_nc_by_catchment,
        catchments_list=catchments_list,
        resample_for_plot=resample_for_plot,
        handle_leap_days=handle_leap_days,
        spatial_member_std=spatial_member_std,
        spatial_wy_snowfall_member_std=spatial_wy_snowfall,
        wy_snowfall_start_year=wy_snowfall_start_year,
        wy_snowfall_end_year=wy_snowfall_end_year,
    )
    return climatology


def _style_axis(ax: Axes, *, show_xlabel: bool, show_ylabel: bool) -> None:
    if show_ylabel:
        ax.set_ylabel("Snowfall fraction")
    else:
        ax.set_ylabel("")
    if show_xlabel:
        ax.set_xlabel("Day of water year")
    else:
        ax.set_xlabel("")
        ax.tick_params(labelbottom=False)

    ax.set_ylim(-0.02, 1.02)
    ax.set_xlim(1, 365)
    ax.set_xticks(_WYDOY_TICKS)
    ax.set_xticklabels(_WYDOY_LABELS)
    ax.grid(True, axis="y")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _style_std_axis(
    ax: Axes,
    *,
    show_ylabel: bool,
    ymax: float,
) -> None:
    if show_ylabel:
        ax.set_ylabel("Ensemble Std Dev")
    else:
        ax.set_ylabel("")
    ax.set_ylim(0.0, ymax)
    ax.set_xlim(1, 365)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(True)
    ax.tick_params(axis="y", colors="#333333")


def _max_std_from_climatology(
    climatology: WyDoyClimatology,
    catchments_list: list[str],
    labels_ordered: list[str],
) -> float:
    peak = 0.0
    for catchment in catchments_list:
        for label in labels_ordered:
            member_std = climatology[catchment][label]["member_std"]
            if member_std is not None:
                peak = max(peak, float(np.nanmax(member_std)))
    return peak if peak > 0 else 0.1


def _plot_catchment_from_climatology(
    ax: Axes,
    ax_std: Axes | None,
    catchment_stats: dict[str, dict[str, np.ndarray | None]],
    *,
    labels_ordered: list[str],
    colors: tuple[str, ...],
) -> list:
    handles: list = []
    t_ax = np.arange(1, 366)

    for i, label in enumerate(labels_ordered):
        stats = catchment_stats[label]
        mean = stats["mean"]
        p10 = stats["p10"]
        p90 = stats["p90"]
        member_std = stats["member_std"]
        color = colors[i % len(colors)]

        if p10 is not None and p90 is not None:
            ax.fill_between(
                t_ax,
                p10,
                p90,
                alpha=ENSEMBLE_BAND_ALPHA,
                color=color,
                linewidth=0,
                zorder=2,
            )
            handles.append(
                Patch(
                    facecolor=color,
                    edgecolor="none",
                    alpha=ENSEMBLE_BAND_ALPHA,
                    label=f"{label} 10–90%",
                )
            )
        (line,) = ax.plot(
            t_ax,
            mean,
            color=color,
            lw=1.8,
            label=label,
            zorder=3,
        )
        handles.append(line)

        if ax_std is not None and member_std is not None:
            (std_line,) = ax_std.plot(
                t_ax,
                member_std,
                color=color,
                ls="--",
                lw=1.5,
                label=f"{label} Std Dev",
                zorder=2,
            )
            handles.append(std_line)

    return handles


def plot_snow_fraction_water_year_doy_publication(
    output_nc_by_catchment: dict[str, dict[str, Path]] | None = None,
    *,
    wydoy_climatology: WyDoyClimatology | None = None,
    cache_path: Path | str | None = None,
    recompute_cache: bool = False,
    catchments_list: list[str] | None = None,
    resample_for_plot: str | None = None,
    handle_leap_days: str = "drop",
    source_colors: tuple[str, ...] = DEFAULT_SOURCE_COLORS,
    figsize: tuple[float, float] | None = None,
    legend_ncol: int | None = None,
) -> tuple[Figure, np.ndarray]:
    """
    Publication panel: snow fraction by water-year day, one row per catchment.

    Each panel shows RF Ensemble (with 10–90% member band when present), CASR, and ERA5.
    Ensemble member std dev per water-year day is overlaid on a right y-axis (dashed).

    Pass precomputed ``wydoy_climatology`` or ``cache_path`` to skip re-reading source NetCDFs.
    When ``cache_path`` is set, climatology is loaded if current; otherwise computed and saved.
    """
    if wydoy_climatology is None:
        if output_nc_by_catchment is None:
            if cache_path is None:
                raise ValueError(
                    "Pass output_nc_by_catchment, wydoy_climatology, or cache_path."
                )
            wydoy_climatology = load_snow_fraction_wydoy_climatology(cache_path)
        elif cache_path is not None:
            wydoy_climatology = ensure_snow_fraction_wydoy_climatology(
                output_nc_by_catchment,
                cache_path,
                catchments_list=catchments_list,
                resample_for_plot=resample_for_plot,
                handle_leap_days=handle_leap_days,
                overwrite=recompute_cache,
            )
        else:
            wydoy_climatology = compute_snow_fraction_wydoy_climatology(
                output_nc_by_catchment,
                catchments_list=catchments_list,
                resample_for_plot=resample_for_plot,
                handle_leap_days=handle_leap_days,
            )

    _apply_publication_style()

    catchments_list = catchments_list or list(wydoy_climatology.keys())
    labels_ordered = list(wydoy_climatology[catchments_list[0]].keys())
    n_rows = len(catchments_list)

    if figsize is None:
        figsize = (10.0, 2.1 * n_rows)

    fig, axes = plt.subplots(n_rows, 1, figsize=figsize, squeeze=False, sharex=True, sharey=True)
    axes = axes.ravel()

    legend_handles: list = []
    legend_labels: list[str] = []
    seen_legend: set[str] = set()

    std_ymax = _max_std_from_climatology(wydoy_climatology, catchments_list, labels_ordered) * 1.05

    for row, catchment in enumerate(catchments_list):
        ax = axes[row]
        ax_std = ax.twinx()
        handles = _plot_catchment_from_climatology(
            ax,
            ax_std,
            wydoy_climatology[catchment],
            labels_ordered=labels_ordered,
            colors=source_colors,
        )
        _style_axis(
            ax,
            show_xlabel=(row == n_rows - 1),
            show_ylabel=True,
        )
        ax.tick_params(axis="y", labelleft=True)
        _style_std_axis(ax_std, show_ylabel=True, ymax=std_ymax)
        ax_std.tick_params(axis="y", labelright=True)
        if row != n_rows - 1:
            ax_std.tick_params(labelbottom=False)
        ax.text(
            0.0,
            1.03,
            _format_catchment_label(catchment),
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=11,
            clip_on=False,
        )

        for h in handles:
            label = h.get_label()
            if label not in seen_legend:
                seen_legend.add(label)
                legend_handles.append(h)
                legend_labels.append(label)

    if legend_handles:
        if legend_ncol is None:
            legend_ncol = min(len(legend_handles), 5)
        fig.legend(
            legend_handles,
            legend_labels,
            loc="lower center",
            bbox_to_anchor=(0.5, 0.01),
            ncol=legend_ncol,
            frameon=False,
            handlelength=2.4,
            columnspacing=1.4,
        )

    fig.subplots_adjust(left=0.12, right=0.88, top=0.97, bottom=0.15, hspace=0.28)
    return fig, axes
