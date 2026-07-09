"""Spatial choropleth panel for mixed-precip ensemble spread."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

import compute_mixed_precip_fractions as mpf
from compute_water_year_snow_metrics import us_water_year_label
from summa_postprocess_specs import DEFAULT_GPEP_ROOT

if TYPE_CHECKING:
    from matplotlib.axes import Axes


def catchment_gpkg_path(catchment: str, gpep_root: Path | str) -> Path:
    return Path(gpep_root) / catchment / "gis" / f"{catchment}_tdx.gpkg"


def infer_gpkg_join_field(gru_values: np.ndarray, gdf: gpd.GeoDataFrame) -> str:
    gru_values = np.asarray(gru_values, dtype=np.int64)
    candidates: list[tuple[int, str]] = []

    for field in ("HRU_ID", "GRU_ID", "gruId", "gru_id", "hruId", "hru_id"):
        if field not in gdf.columns:
            continue
        gpkg_ids = gdf[field].astype(np.int64).values
        n_match = int(np.isin(gru_values, gpkg_ids).sum())
        if n_match:
            candidates.append((n_match, field))

    if not candidates:
        raise ValueError(
            "Could not match id values to the GeoPackage. "
            f"id sample={gru_values[:8].tolist()} …; gpkg columns={list(gdf.columns)}"
        )

    candidates.sort(reverse=True)
    return candidates[0][1]


def _format_catchment_label(catchment: str) -> str:
    return catchment.strip().replace("_", " ").title()


def _apply_spatial_publication_style() -> None:
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.titleweight": "normal",
        }
    )


def _resolve_year_range(
    start_year: int | None,
    end_year: int | None,
    default_start: int | None,
    default_end: int | None,
) -> tuple[int | None, int | None]:
    if start_year is None:
        start_year = default_start
    if end_year is None:
        end_year = default_end
    if start_year is not None and end_year is not None and start_year > end_year:
        raise ValueError(f"start_year ({start_year}) must be <= end_year ({end_year})")
    return start_year, end_year


def spatial_member_std_per_gru(
    ds: xr.Dataset,
    var: str = mpf.VAR_FRAC,
) -> xr.DataArray:
    """Std dev across ensemble members of the time-mean snow fraction, per GRU."""
    if var not in ds:
        raise KeyError(f"{var!r} not in dataset; found {list(ds.data_vars)}")
    da = ds[var]
    if "member" not in da.dims:
        raise ValueError(f"{var} has no member dimension; std map requires an ensemble NetCDF.")
    if "time" not in da.dims:
        raise ValueError(f"{var} has no time dimension.")
    out = da.mean(dim="time", skipna=True).std(dim="member", skipna=True)
    if "gru" not in out.dims:
        raise ValueError(f"Expected gru dimension after reduction; got {out.dims}")
    return out


def _wy_total_snowfall_by_water_year(
    ds: xr.Dataset,
    *,
    start_year: int | None = None,
    end_year: int | None = None,
    var: str = mpf.VAR_SNOW,
) -> xr.DataArray:
    """Sum ``var`` over time per water year, GRU, and member."""
    if var not in ds:
        raise KeyError(f"{var!r} not in dataset; found {list(ds.data_vars)}")
    da = ds[var]
    if "member" not in da.dims:
        raise ValueError(f"{var} has no member dimension; requires an ensemble NetCDF.")
    if "time" not in da.dims:
        raise ValueError(f"{var} has no time dimension.")

    wy = us_water_year_label(da["time"])
    da = da.assign_coords(water_year=("time", wy.data))
    wy_totals = da.groupby("water_year").sum(skipna=True)

    if start_year is not None or end_year is not None:
        lo = start_year if start_year is not None else int(wy_totals["water_year"].min())
        hi = end_year if end_year is not None else int(wy_totals["water_year"].max())
        wy_totals = wy_totals.sel(water_year=slice(lo, hi))
    return wy_totals


def spatial_mean_wy_snowfall_member_std_per_gru(
    ds: xr.Dataset,
    *,
    start_year: int | None = None,
    end_year: int | None = None,
    var: str = mpf.VAR_SNOW,
) -> xr.DataArray:
    """
    Mean across water years of normalized member std dev of total snowfall.

    For each water year: sum ``scalarSnowfall_total`` over time per GRU and member,
    compute std across members, divide by the mean member total for that water year,
    then average over the selected water years.
    """
    wy_totals = _wy_total_snowfall_by_water_year(
        ds,
        start_year=start_year,
        end_year=end_year,
        var=var,
    )
    wy_member_std = wy_totals.std(dim="member", skipna=True)
    wy_member_mean = wy_totals.mean(dim="member", skipna=True)
    wy_normalized = (wy_member_std / wy_member_mean).where(wy_member_mean > 0)
    out = wy_normalized.mean(dim="water_year", skipna=True)
    if "gru" not in out.dims:
        raise ValueError(f"Expected gru dimension after reduction; got {out.dims}")
    return out


def spatial_mean_wy_snowfall_per_gru(
    ds: xr.Dataset,
    *,
    start_year: int | None = None,
    end_year: int | None = None,
    var: str = mpf.VAR_SNOW,
) -> xr.DataArray:
    """Mean across water years of ensemble-mean total snowfall."""
    wy_totals = _wy_total_snowfall_by_water_year(
        ds,
        start_year=start_year,
        end_year=end_year,
        var=var,
    )
    out = wy_totals.mean(dim="member", skipna=True).mean(dim="water_year", skipna=True)
    if "gru" not in out.dims:
        raise ValueError(f"Expected gru dimension after reduction; got {out.dims}")
    return out


def build_spatial_member_std_gdf(
    nc_path: Path | str | None = None,
    *,
    catchment: str,
    gpep_root: Path | str = DEFAULT_GPEP_ROOT,
    var: str = mpf.VAR_FRAC,
    gpkg_path: Path | str | None = None,
    gpkg_id_field: str | None = None,
    gru_values: np.ndarray | None = None,
    std_values: np.ndarray | None = None,
) -> gpd.GeoDataFrame:
    if gru_values is None or std_values is None:
        if nc_path is None:
            raise ValueError("Pass nc_path or precomputed gru_values and std_values.")
        nc_path = Path(nc_path)
        ds = xr.open_dataset(nc_path)
        try:
            std_da = spatial_member_std_per_gru(ds, var).load()
            gru_values = np.asarray(std_da["gru"].values, dtype=np.int64)
            std_values = np.asarray(std_da.values, dtype=float)
        finally:
            ds.close()

    return _merge_gru_values_gdf(
        catchment,
        gru_values,
        std_values,
        gpep_root=gpep_root,
        gpkg_path=gpkg_path,
        gpkg_id_field=gpkg_id_field,
    )


def build_spatial_mean_wy_snowfall_member_std_gdf(
    nc_path: Path | str | None = None,
    *,
    catchment: str,
    gpep_root: Path | str = DEFAULT_GPEP_ROOT,
    var: str = mpf.VAR_SNOW,
    start_year: int | None = None,
    end_year: int | None = None,
    gpkg_path: Path | str | None = None,
    gpkg_id_field: str | None = None,
    gru_values: np.ndarray | None = None,
    std_values: np.ndarray | None = None,
) -> gpd.GeoDataFrame:
    if gru_values is None or std_values is None:
        if nc_path is None:
            raise ValueError("Pass nc_path or precomputed gru_values and std_values.")
        nc_path = Path(nc_path)
        ds = xr.open_dataset(nc_path, decode_times=True)
        try:
            std_da = spatial_mean_wy_snowfall_member_std_per_gru(
                ds,
                start_year=start_year,
                end_year=end_year,
                var=var,
            ).load()
            gru_values = np.asarray(std_da["gru"].values, dtype=np.int64)
            std_values = np.asarray(std_da.values, dtype=float)
        finally:
            ds.close()

    return _merge_gru_values_gdf(
        catchment,
        gru_values,
        std_values,
        gpep_root=gpep_root,
        gpkg_path=gpkg_path,
        gpkg_id_field=gpkg_id_field,
    )


def build_spatial_mean_wy_snowfall_gdf(
    nc_path: Path | str | None = None,
    *,
    catchment: str,
    gpep_root: Path | str = DEFAULT_GPEP_ROOT,
    var: str = mpf.VAR_SNOW,
    start_year: int | None = None,
    end_year: int | None = None,
    gpkg_path: Path | str | None = None,
    gpkg_id_field: str | None = None,
    gru_values: np.ndarray | None = None,
    mean_values: np.ndarray | None = None,
) -> gpd.GeoDataFrame:
    if gru_values is None or mean_values is None:
        if nc_path is None:
            raise ValueError("Pass nc_path or precomputed gru_values and mean_values.")
        nc_path = Path(nc_path)
        ds = xr.open_dataset(nc_path, decode_times=True)
        try:
            mean_da = spatial_mean_wy_snowfall_per_gru(
                ds,
                start_year=start_year,
                end_year=end_year,
                var=var,
            ).load()
            gru_values = np.asarray(mean_da["gru"].values, dtype=np.int64)
            mean_values = np.asarray(mean_da.values, dtype=float)
        finally:
            ds.close()

    return _merge_gru_values_gdf(
        catchment,
        gru_values,
        mean_values,
        gpep_root=gpep_root,
        gpkg_path=gpkg_path,
        gpkg_id_field=gpkg_id_field,
    )


def _merge_gru_values_gdf(
    catchment: str,
    gru_values: np.ndarray,
    std_values: np.ndarray,
    *,
    gpep_root: Path | str = DEFAULT_GPEP_ROOT,
    gpkg_path: Path | str | None = None,
    gpkg_id_field: str | None = None,
) -> gpd.GeoDataFrame:
    gpkg_path = (
        Path(gpkg_path)
        if gpkg_path is not None
        else catchment_gpkg_path(catchment, gpep_root)
    )
    if not gpkg_path.is_file():
        raise FileNotFoundError(f"GeoPackage not found: {gpkg_path}")

    gdf = gpd.read_file(gpkg_path)
    join_field = gpkg_id_field or infer_gpkg_join_field(gru_values, gdf)
    merge_df = pd.DataFrame(
        {
            join_field: np.asarray(gru_values, dtype=np.int64),
            "value": np.asarray(std_values, dtype=float),
        }
    )

    gdf = gdf.copy()
    gdf[join_field] = gdf[join_field].astype(np.int64)
    return gdf.merge(merge_df, on=join_field, how="left")


def _snowfall_std_legend_label() -> str:
    return "(-)"


def _snowfall_mm_legend_label() -> str:
    return "mm"


def _shared_vmax(gdf_by_catchment: dict[str, gpd.GeoDataFrame], panel_catchments: list[str]) -> float:
    vals = np.concatenate(
        [gdf_by_catchment[c]["value"].to_numpy(dtype=float) for c in panel_catchments]
    )
    finite = vals[np.isfinite(vals)]
    return float(np.nanpercentile(finite, 98)) if finite.size else 1.0


def _scale_gdf_values_for_plot(
    gdf_by_catchment: dict[str, gpd.GeoDataFrame],
    panel_catchments: list[str],
    factor: float,
) -> dict[str, gpd.GeoDataFrame]:
    scaled: dict[str, gpd.GeoDataFrame] = {}
    for catchment in panel_catchments:
        plot_gdf = gdf_by_catchment[catchment].copy()
        plot_gdf["value"] = plot_gdf["value"] * factor
        scaled[catchment] = plot_gdf
    return scaled


def _scale_limit(limit: float | None, factor: float) -> float | None:
    if limit is None:
        return None
    return limit * factor


def _plot_spatial_catchment_panel(
    gdf_by_catchment: dict[str, gpd.GeoDataFrame],
    panel_catchments: list[str],
    *,
    figsize: tuple[float, float] | None,
    cmap: str,
    vmin: float | None,
    vmax: float | None,
    legend_label: str,
    edgecolor: str,
    linewidth: float,
) -> tuple:
    n_cols = len(panel_catchments)
    if figsize is None:
        figsize = (3.6 * n_cols, 3.6)

    fig, axes = plt.subplots(1, n_cols, figsize=figsize, squeeze=False, constrained_layout=True)
    axes = axes.ravel()

    for col, catchment in enumerate(panel_catchments):
        _plot_gdf_column(
            gdf_by_catchment[catchment],
            axes[col],
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            legend_label=legend_label,
            edgecolor=edgecolor,
            linewidth=linewidth,
            show_legend=(col == n_cols - 1),
        )

    _add_spatial_catchment_titles(fig, axes, panel_catchments, n_cols=n_cols)
    return fig, axes


def _add_spatial_catchment_titles(
    fig,
    axes,
    panel_catchments: list[str],
    *,
    n_cols: int,
) -> None:
    fig.set_constrained_layout_pads(rect=(0, 0, 1, 0.94))
    fig.canvas.draw()
    top_y = max(axes[col].get_position().y1 for col in range(n_cols))
    title_pad = 0.012
    for col, catchment in enumerate(panel_catchments):
        pos = axes[col].get_position()
        fig.text(
            (pos.x0 + pos.x1) / 2,
            top_y + title_pad,
            _format_catchment_label(catchment),
            ha="center",
            va="bottom",
            fontsize=11,
            clip_on=False,
        )


def _plot_spatial_multi_row_catchment_panel(
    gdf_by_row: tuple[dict[str, gpd.GeoDataFrame], ...],
    panel_catchments: list[str],
    *,
    row_labels: tuple[str, ...],
    figsize: tuple[float, float] | None,
    cmap_by_row: tuple[str, ...],
    vmin_by_row: tuple[float | None, ...],
    vmax_by_row: tuple[float | None, ...],
    legend_labels: tuple[str, ...],
    edgecolor: str,
    linewidth: float,
    row_height: float = 2.2,
    row_label_x: float = -0.10,
) -> tuple:
    n_rows = len(gdf_by_row)
    n_cols = len(panel_catchments)
    if len(row_labels) != n_rows:
        raise ValueError(f"row_labels length {len(row_labels)} != n_rows {n_rows}")
    if figsize is None:
        figsize = (3.6 * n_cols, row_height * n_rows)

    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=figsize,
        squeeze=False,
        constrained_layout=True,
    )

    for col, catchment in enumerate(panel_catchments):
        for row, (gdf_by_catchment, row_label, cmap, vmin, vmax, legend_label) in enumerate(
            zip(
                gdf_by_row,
                row_labels,
                cmap_by_row,
                vmin_by_row,
                vmax_by_row,
                legend_labels,
            )
        ):
            ax = axes[row, col]
            _plot_gdf_column(
                gdf_by_catchment[catchment],
                ax,
                cmap=cmap,
                vmin=vmin,
                vmax=vmax,
                legend_label=legend_label,
                edgecolor=edgecolor,
                linewidth=linewidth,
                show_legend=(col == n_cols - 1),
            )
            if col == 0:
                ax.text(
                    row_label_x,
                    0.5,
                    row_label,
                    transform=ax.transAxes,
                    va="center",
                    ha="right",
                    rotation=0,
                    fontsize=10,
                    clip_on=False,
                )

    _add_spatial_catchment_titles(fig, axes[0], panel_catchments, n_cols=n_cols)
    return fig, axes


def _plot_gdf_column(
    plot_gdf: gpd.GeoDataFrame,
    ax: Axes,
    *,
    cmap: str,
    vmin: float | None,
    vmax: float | None,
    legend_label: str,
    edgecolor: str,
    linewidth: float,
    show_legend: bool,
) -> None:
    plot_kw = dict(
        column="value",
        cmap=cmap,
        legend=show_legend,
        ax=ax,
        edgecolor=edgecolor,
        linewidth=linewidth,
        missing_kwds={"color": "lightgrey", "label": "no data"},
    )
    if show_legend:
        plot_kw["legend_kwds"] = {"label": legend_label, "shrink": 0.85}
    if vmin is not None:
        plot_kw["vmin"] = vmin
    if vmax is not None:
        plot_kw["vmax"] = vmax
    plot_gdf.plot(**plot_kw)
    ax.set_axis_off()


def _resolve_catchments_list(
    catchments_list: list[str] | None,
    *,
    spatial_cache: dict | None,
    output_nc_by_catchment: dict | None,
) -> list[str]:
    if catchments_list is not None:
        return catchments_list
    if spatial_cache is not None:
        return list(spatial_cache.keys())
    if output_nc_by_catchment is not None:
        return list(output_nc_by_catchment.keys())
    raise ValueError("catchments_list is required when using cache_path alone.")


def _collect_spatial_member_std_gdfs(
    catchments_list: list[str],
    *,
    output_nc_by_catchment: dict[str, dict[str, Path | str]] | None,
    spatial_member_std: dict | None,
    source_label: str,
    gpep_root: Path | str,
    var: str,
    gpkg_id_field: str | None,
) -> tuple[list[str], dict[str, gpd.GeoDataFrame]]:
    panel_catchments: list[str] = []
    gdf_by_catchment: dict[str, gpd.GeoDataFrame] = {}
    for catchment in catchments_list:
        if spatial_member_std is not None:
            catchment_cache = spatial_member_std.get(catchment, {})
            if source_label not in catchment_cache:
                continue
            entry = catchment_cache[source_label]
            plot_gdf = build_spatial_member_std_gdf(
                catchment=catchment,
                gpep_root=gpep_root,
                var=var,
                gpkg_id_field=gpkg_id_field,
                gru_values=entry["gru"],
                std_values=entry["value"],
            )
        else:
            path_dict = output_nc_by_catchment.get(catchment, {})
            if source_label not in path_dict:
                continue
            plot_gdf = build_spatial_member_std_gdf(
                path_dict[source_label],
                catchment=catchment,
                gpep_root=gpep_root,
                var=var,
                gpkg_id_field=gpkg_id_field,
            )
        panel_catchments.append(catchment)
        gdf_by_catchment[catchment] = plot_gdf
    return panel_catchments, gdf_by_catchment


def _collect_spatial_wy_snowfall_gdfs(
    catchments_list: list[str],
    *,
    output_nc_by_catchment: dict[str, dict[str, Path | str]] | None,
    spatial_wy_snowfall_member_std: dict | None,
    source_label: str,
    gpep_root: Path | str,
    start_year: int | None,
    end_year: int | None,
    var: str,
    gpkg_id_field: str | None,
) -> tuple[list[str], dict[str, gpd.GeoDataFrame]]:
    panel_catchments: list[str] = []
    gdf_by_catchment: dict[str, gpd.GeoDataFrame] = {}
    for catchment in catchments_list:
        if spatial_wy_snowfall_member_std is not None:
            catchment_cache = spatial_wy_snowfall_member_std.get(catchment, {})
            if source_label not in catchment_cache:
                continue
            entry = catchment_cache[source_label]
            plot_gdf = build_spatial_mean_wy_snowfall_member_std_gdf(
                catchment=catchment,
                gpep_root=gpep_root,
                var=var,
                start_year=start_year,
                end_year=end_year,
                gpkg_id_field=gpkg_id_field,
                gru_values=entry["gru"],
                std_values=entry["value"],
            )
        else:
            path_dict = output_nc_by_catchment.get(catchment, {})
            if source_label not in path_dict:
                continue
            plot_gdf = build_spatial_mean_wy_snowfall_member_std_gdf(
                Path(path_dict[source_label]),
                catchment=catchment,
                gpep_root=gpep_root,
                var=var,
                start_year=start_year,
                end_year=end_year,
                gpkg_id_field=gpkg_id_field,
            )
        panel_catchments.append(catchment)
        gdf_by_catchment[catchment] = plot_gdf
    return panel_catchments, gdf_by_catchment


def _collect_spatial_wy_snowfall_mean_gdfs(
    catchments_list: list[str],
    *,
    output_nc_by_catchment: dict[str, dict[str, Path | str]] | None,
    source_label: str,
    gpep_root: Path | str,
    start_year: int | None,
    end_year: int | None,
    var: str,
    gpkg_id_field: str | None,
) -> tuple[list[str], dict[str, gpd.GeoDataFrame]]:
    if output_nc_by_catchment is None:
        raise ValueError(
            "output_nc_by_catchment is required to compute mean WY snowfall spatial maps."
        )

    panel_catchments: list[str] = []
    gdf_by_catchment: dict[str, gpd.GeoDataFrame] = {}
    for catchment in catchments_list:
        path_dict = output_nc_by_catchment.get(catchment, {})
        if source_label not in path_dict:
            continue
        plot_gdf = build_spatial_mean_wy_snowfall_gdf(
            Path(path_dict[source_label]),
            catchment=catchment,
            gpep_root=gpep_root,
            var=var,
            start_year=start_year,
            end_year=end_year,
            gpkg_id_field=gpkg_id_field,
        )
        panel_catchments.append(catchment)
        gdf_by_catchment[catchment] = plot_gdf
    return panel_catchments, gdf_by_catchment


def plot_spatial_member_std_panel(
    output_nc_by_catchment: dict[str, dict[str, Path | str]] | None = None,
    *,
    cache_path: Path | str | None = None,
    spatial_member_std: dict | None = None,
    source_label: str = "RF Ensemble",
    catchments_list: list[str] | None = None,
    gpep_root: Path | str = DEFAULT_GPEP_ROOT,
    var: str = mpf.VAR_FRAC,
    gpkg_id_field: str | None = None,
    figsize: tuple[float, float] | None = None,
    cmap: str = "YlOrRd",
    vmin: float | None = 0.0,
    vmax: float | None = None,
    edgecolor: str = "black",
    linewidth: float = 0.15,
    legend_label: str = "Member std dev",
) -> tuple:
    """
    One row × catchments: std dev of time-mean snow fraction across ensemble members.

    Pass ``cache_path`` (same file as ``wydoy_cache``) or preloaded ``spatial_member_std``
    to skip re-reading source NetCDFs.
    """
    _apply_spatial_publication_style()

    if spatial_member_std is None:
        if cache_path is None:
            if output_nc_by_catchment is None:
                raise ValueError(
                    "Pass output_nc_by_catchment, cache_path, or spatial_member_std."
                )
        else:
            from mixed_precip_publication_plot import load_spatial_member_std_cache

            spatial_member_std = load_spatial_member_std_cache(cache_path)

    catchments_list = _resolve_catchments_list(
        catchments_list,
        spatial_cache=spatial_member_std,
        output_nc_by_catchment=output_nc_by_catchment,
    )

    panel_catchments, gdf_by_catchment = _collect_spatial_member_std_gdfs(
        catchments_list,
        output_nc_by_catchment=output_nc_by_catchment,
        spatial_member_std=spatial_member_std,
        source_label=source_label,
        gpep_root=gpep_root,
        var=var,
        gpkg_id_field=gpkg_id_field,
    )

    if not panel_catchments:
        raise ValueError(
            f"No catchments with source {source_label!r} were found in output_nc_by_catchment."
        )

    if vmax is None:
        vmax = _shared_vmax(gdf_by_catchment, panel_catchments)

    fig, axes = _plot_spatial_catchment_panel(
        gdf_by_catchment,
        panel_catchments,
        figsize=figsize,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        legend_label=legend_label,
        edgecolor=edgecolor,
        linewidth=linewidth,
    )

    return fig, axes, gdf_by_catchment


def plot_spatial_wy_snowfall_member_std_panel(
    output_nc_by_catchment: dict[str, dict[str, Path | str]] | None = None,
    *,
    cache_path: Path | str | None = None,
    spatial_wy_snowfall_member_std: dict | None = None,
    source_label: str = "RF Ensemble",
    catchments_list: list[str] | None = None,
    gpep_root: Path | str = DEFAULT_GPEP_ROOT,
    start_year: int | None = None,
    end_year: int | None = None,
    default_start_year: int | None = None,
    default_end_year: int | None = None,
    var: str = mpf.VAR_SNOW,
    gpkg_id_field: str | None = None,
    figsize: tuple[float, float] | None = None,
    cmap: str = "YlOrRd",
    vmin: float | None = 0.0,
    vmax: float | None = None,
    edgecolor: str = "black",
    linewidth: float = 0.15,
    legend_label: str | None = None,
) -> tuple:
    """
    One row × catchments: mean per-WY normalized member std of total snowfall (RF).

    For each water year, total ``scalarSnowfall_total`` is summed per GRU and member;
    member std is divided by the mean member total for that water year, then averaged
    over the selected water years.

    Pass ``cache_path`` (same file as ``wydoy_cache``) or preloaded
    ``spatial_wy_snowfall_member_std`` to skip re-reading source NetCDFs.
    """
    _apply_spatial_publication_style()

    start_year, end_year = _resolve_year_range(
        start_year,
        end_year,
        default_start_year,
        default_end_year,
    )

    if spatial_wy_snowfall_member_std is None:
        if cache_path is None:
            if output_nc_by_catchment is None:
                raise ValueError(
                    "Pass output_nc_by_catchment, cache_path, or spatial_wy_snowfall_member_std."
                )
        else:
            from mixed_precip_publication_plot import load_spatial_wy_snowfall_member_std_cache

            spatial_wy_snowfall_member_std = load_spatial_wy_snowfall_member_std_cache(
                cache_path,
                start_year=start_year,
                end_year=end_year,
            )

    catchments_list = _resolve_catchments_list(
        catchments_list,
        spatial_cache=spatial_wy_snowfall_member_std,
        output_nc_by_catchment=output_nc_by_catchment,
    )

    panel_catchments, gdf_by_catchment = _collect_spatial_wy_snowfall_gdfs(
        catchments_list,
        output_nc_by_catchment=output_nc_by_catchment,
        spatial_wy_snowfall_member_std=spatial_wy_snowfall_member_std,
        source_label=source_label,
        gpep_root=gpep_root,
        start_year=start_year,
        end_year=end_year,
        var=var,
        gpkg_id_field=gpkg_id_field,
    )

    if not panel_catchments:
        raise ValueError(
            f"No catchments with source {source_label!r} were found in output_nc_by_catchment."
        )

    if vmax is None:
        vmax = _shared_vmax(gdf_by_catchment, panel_catchments)

    if legend_label is None:
        legend_label = _snowfall_std_legend_label()

    fig, axes = _plot_spatial_catchment_panel(
        gdf_by_catchment,
        panel_catchments,
        figsize=figsize,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        legend_label=legend_label,
        edgecolor=edgecolor,
        linewidth=linewidth,
    )

    return fig, axes, gdf_by_catchment


def plot_spatial_ensemble_spread_panel(
    output_nc_by_catchment: dict[str, dict[str, Path | str]] | None = None,
    *,
    cache_path: Path | str | None = None,
    spatial_member_std: dict | None = None,
    spatial_wy_snowfall_member_std: dict | None = None,
    source_label: str = "RF Ensemble",
    catchments_list: list[str] | None = None,
    gpep_root: Path | str = DEFAULT_GPEP_ROOT,
    start_year: int | None = None,
    end_year: int | None = None,
    default_start_year: int | None = None,
    default_end_year: int | None = None,
    frac_var: str = mpf.VAR_FRAC,
    snow_var: str = mpf.VAR_SNOW,
    gpkg_id_field: str | None = None,
    figsize: tuple[float, float] | None = None,
    frac_cmap: str = "YlOrRd",
    snow_cmap: str = "BuPu",
    ratio_cmap: str = "YlGn",
    cmap: str | None = None,
    frac_vmin: float | None = 0.0,
    frac_vmax: float | None = None,
    snow_vmin: float | None = 0.0,
    snow_vmax: float | None = None,
    edgecolor: str = "black",
    linewidth: float = 0.15,
    frac_legend_label: str = "mm",
    snow_legend_label: str | None = None,
    snow_mean_legend_label: str | None = None,
    include_wy_mean: bool = True,
) -> tuple:
    """
    Multi-row × catchments: snow-fraction member std, mean WY snowfall, and
    normalized member std / WY mean (-).

    Pass ``cache_path`` (same file as ``wydoy_cache``) or preloaded spatial caches
    to skip re-reading source NetCDFs for cached metrics. ``output_nc_by_catchment``
    is required when ``include_wy_mean`` is True.
    """
    _apply_spatial_publication_style()

    start_year, end_year = _resolve_year_range(
        start_year,
        end_year,
        default_start_year,
        default_end_year,
    )

    if include_wy_mean and output_nc_by_catchment is None:
        raise ValueError(
            "output_nc_by_catchment is required when include_wy_mean is True."
        )

    if spatial_member_std is None and cache_path is not None:
        from mixed_precip_publication_plot import load_spatial_member_std_cache

        spatial_member_std = load_spatial_member_std_cache(cache_path)

    if spatial_wy_snowfall_member_std is None and cache_path is not None:
        from mixed_precip_publication_plot import load_spatial_wy_snowfall_member_std_cache

        spatial_wy_snowfall_member_std = load_spatial_wy_snowfall_member_std_cache(
            cache_path,
            start_year=start_year,
            end_year=end_year,
        )

    if spatial_member_std is None and spatial_wy_snowfall_member_std is None:
        if output_nc_by_catchment is None:
            raise ValueError(
                "Pass output_nc_by_catchment, cache_path, or preloaded spatial caches."
            )

    catchments_list = _resolve_catchments_list(
        catchments_list,
        spatial_cache=spatial_member_std or spatial_wy_snowfall_member_std,
        output_nc_by_catchment=output_nc_by_catchment,
    )

    _, gdf_member_std = _collect_spatial_member_std_gdfs(
        catchments_list,
        output_nc_by_catchment=output_nc_by_catchment,
        spatial_member_std=spatial_member_std,
        source_label=source_label,
        gpep_root=gpep_root,
        var=frac_var,
        gpkg_id_field=gpkg_id_field,
    )
    _, gdf_wy_snowfall = _collect_spatial_wy_snowfall_gdfs(
        catchments_list,
        output_nc_by_catchment=output_nc_by_catchment,
        spatial_wy_snowfall_member_std=spatial_wy_snowfall_member_std,
        source_label=source_label,
        gpep_root=gpep_root,
        start_year=start_year,
        end_year=end_year,
        var=snow_var,
        gpkg_id_field=gpkg_id_field,
    )

    gdf_wy_snowfall_mean: dict[str, gpd.GeoDataFrame] = {}
    if include_wy_mean:
        _, gdf_wy_snowfall_mean = _collect_spatial_wy_snowfall_mean_gdfs(
            catchments_list,
            output_nc_by_catchment=output_nc_by_catchment,
            source_label=source_label,
            gpep_root=gpep_root,
            start_year=start_year,
            end_year=end_year,
            var=snow_var,
            gpkg_id_field=gpkg_id_field,
        )

    required_gdfs = [gdf_member_std, gdf_wy_snowfall]
    if include_wy_mean:
        required_gdfs.append(gdf_wy_snowfall_mean)

    panel_catchments = [
        catchment
        for catchment in catchments_list
        if all(catchment in gdf for gdf in required_gdfs)
    ]
    if not panel_catchments:
        raise ValueError(
            f"No catchments with source {source_label!r} were found for all spatial metrics."
        )

    if frac_vmax is None:
        frac_vmax = _shared_vmax(gdf_member_std, panel_catchments)
    if snow_vmax is None:
        snow_vmax = _shared_vmax(gdf_wy_snowfall, panel_catchments)
    if snow_legend_label is None:
        snow_legend_label = _snowfall_std_legend_label()
    if snow_mean_legend_label is None:
        snow_mean_legend_label = _snowfall_mm_legend_label()

    if cmap is not None:
        frac_cmap = snow_cmap = cmap

    gdf_by_row: list[dict[str, gpd.GeoDataFrame]] = [gdf_member_std]
    row_labels: list[str] = ["Snowfall Fraction Std Dev"]
    cmap_by_row: list[str] = [frac_cmap]
    vmin_by_row: list[float | None] = [frac_vmin]
    vmax_by_row: list[float | None] = [frac_vmax]
    legend_labels: list[str] = [frac_legend_label]

    if include_wy_mean:
        snow_mean_scale = 1.0
        gdf_wy_snowfall_mean_plot = _scale_gdf_values_for_plot(
            gdf_wy_snowfall_mean,
            panel_catchments,
            snow_mean_scale,
        )
        snow_mean_vmax = _shared_vmax(gdf_wy_snowfall_mean, panel_catchments)
        gdf_by_row.append(gdf_wy_snowfall_mean_plot)
        row_labels.append("Snowfall WY Mean")
        cmap_by_row.append(snow_cmap)
        vmin_by_row.append(_scale_limit(snow_vmin, snow_mean_scale))
        vmax_by_row.append(_scale_limit(snow_mean_vmax, snow_mean_scale))
        legend_labels.append(snow_mean_legend_label)

    gdf_by_row.append(gdf_wy_snowfall)
    row_labels.append("Std Dev / WY Mean")
    cmap_by_row.append(ratio_cmap)
    vmin_by_row.append(snow_vmin)
    vmax_by_row.append(snow_vmax)
    legend_labels.append(snow_legend_label)

    fig, axes = _plot_spatial_multi_row_catchment_panel(
        tuple(gdf_by_row),
        panel_catchments,
        row_labels=tuple(row_labels),
        figsize=figsize,
        cmap_by_row=tuple(cmap_by_row),
        vmin_by_row=tuple(vmin_by_row),
        vmax_by_row=tuple(vmax_by_row),
        legend_labels=tuple(legend_labels),
        edgecolor=edgecolor,
        linewidth=linewidth,
        row_label_x=-0.14,
    )

    gdf_by_metric = {
        "member_std": {c: gdf_member_std[c] for c in panel_catchments},
        "wy_snowfall_std_ratio": {c: gdf_wy_snowfall[c] for c in panel_catchments},
    }
    if include_wy_mean:
        gdf_by_metric["wy_snowfall_mean"] = {
            c: gdf_wy_snowfall_mean[c] for c in panel_catchments
        }
    return fig, axes, gdf_by_metric
