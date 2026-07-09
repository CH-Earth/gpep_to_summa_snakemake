"""Spatial choropleth panel for snow-metrics model vs MOD10A1 comparison."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

import compute_water_year_snow_metrics as wysm
from snow_metrics_obs_csv import load_obs_csv_spatial_values
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


def _default_cmap_for_var(var: str, cmap: str) -> str:
    if cmap != "viridis":
        return cmap
    if var == wysm.VAR_SCF_DAYS:
        return "YlGnBu"
    if var == wysm.VAR_SNOW_DISAPPEARANCE:
        return "plasma"
    return cmap


def _display_scale(var: str) -> float:
    if var == wysm.VAR_SCF_DAYS:
        return 100.0
    return 1.0


def _scale_gdf_for_display(plot_gdf: gpd.GeoDataFrame, var: str) -> gpd.GeoDataFrame:
    scale = _display_scale(var)
    if scale == 1.0:
        return plot_gdf
    out = plot_gdf.copy()
    for col in ("model_value", "obs_value", "diff_value"):
        if col in out.columns:
            out[col] = out[col] * scale
    return out


def _spatial_var_legend_label(var: str) -> str:
    if var == wysm.VAR_SCF_DAYS:
        return "SCF (%)"
    if var == wysm.VAR_SNOW_DISAPPEARANCE:
        return "SDD (doy)"
    return var


def _default_diff_label(var: str) -> str:
    if var == wysm.VAR_SCF_DAYS:
        return "Difference SCF (%)"
    if var == wysm.VAR_SNOW_DISAPPEARANCE:
        return "Difference SDD (doy)"
    return "Difference"


def _default_value_limits(var: str, values: np.ndarray) -> tuple[float | None, float | None]:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return None, None
    if var == wysm.VAR_SCF_DAYS:
        return 0.0, 100.0
    if var == wysm.VAR_SNOW_DISAPPEARANCE:
        return float(np.nanmin(finite)), float(np.nanmax(finite))
    return float(np.nanmin(finite)), float(np.nanmax(finite))


def _default_diff_limits(diff: np.ndarray) -> tuple[float, float]:
    finite = np.asarray(diff, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return -1.0, 1.0
    cap = float(np.nanpercentile(np.abs(finite), 98))
    if cap <= 0:
        cap = 1.0
    return -cap, cap


def load_model_spatial_values(
    nc_path: Path | str,
    var: str,
    *,
    water_year: int | None = None,
) -> pd.DataFrame:
    nc_path = Path(nc_path)
    if not nc_path.is_file():
        raise FileNotFoundError(f"Snow-metrics NetCDF not found: {nc_path}")

    ds = xr.open_dataset(nc_path)
    try:
        if var not in ds:
            raise KeyError(f"{var!r} not in {nc_path}; found {list(ds.data_vars)}")
        da = ds[var]
        if water_year is not None:
            da = da.sel(water_year=int(water_year))
        reduce_dims = tuple(d for d in ("water_year", "member") if d in da.dims)
        if reduce_dims:
            da = da.mean(dim=reduce_dims, skipna=True)
        if "gru" not in da.dims:
            raise ValueError(f"Expected gru dimension in {nc_path}; got {da.dims}")
        return pd.DataFrame(
            {
                "gru": np.asarray(da["gru"].values, dtype=np.int64),
                "model_value": np.asarray(da.values, dtype=float),
            }
        )
    finally:
        ds.close()


def build_spatial_model_obs_gdf(
    nc_path: Path | str,
    csv_path: Path | str,
    *,
    catchment: str,
    gpep_root: Path | str = DEFAULT_GPEP_ROOT,
    var: str = wysm.VAR_SCF_DAYS,
    water_year: int | None = None,
    stat: str = "mean",
    gpkg_path: Path | str | None = None,
    gpkg_id_field: str | None = None,
    csv_id_field: str | None = None,
) -> tuple[gpd.GeoDataFrame, str]:
    gpkg_path = (
        Path(gpkg_path)
        if gpkg_path is not None
        else catchment_gpkg_path(catchment, gpep_root)
    )
    if not gpkg_path.is_file():
        raise FileNotFoundError(f"GeoPackage not found: {gpkg_path}")

    model_df = load_model_spatial_values(nc_path, var, water_year=water_year)
    obs_df = load_obs_csv_spatial_values(
        csv_path,
        var,
        water_year=water_year,
        stat=stat,
        id_field=csv_id_field,
    )
    obs_id_field = obs_df.columns[0]

    gdf = gpd.read_file(gpkg_path)
    model_join = gpkg_id_field or infer_gpkg_join_field(model_df["gru"].values, gdf)
    obs_join = gpkg_id_field or infer_gpkg_join_field(
        obs_df[obs_id_field].astype(np.int64).values,
        gdf,
    )
    join_field = model_join
    if obs_join != model_join:
        print(
            f"Note: model joins on {model_join!r}, obs on {obs_join!r}; using {join_field!r}."
        )

    model_merge = model_df.rename(columns={"gru": join_field, "model_value": "model_value"})
    obs_merge = obs_df.rename(columns={obs_id_field: join_field, "value": "obs_value"})

    gdf = gdf.copy()
    gdf[join_field] = gdf[join_field].astype(np.int64)
    plot_gdf = gdf.merge(model_merge[[join_field, "model_value"]], on=join_field, how="left")
    plot_gdf = plot_gdf.merge(obs_merge[[join_field, "obs_value"]], on=join_field, how="left")
    plot_gdf["diff_value"] = plot_gdf["model_value"] - plot_gdf["obs_value"]
    return _scale_gdf_for_display(plot_gdf, var), join_field


def _plot_gdf_column(
    plot_gdf: gpd.GeoDataFrame,
    column: str,
    ax: Axes,
    *,
    cmap: str,
    vmin: float | None = None,
    vmax: float | None = None,
    legend_label: str | None = None,
    edgecolor: str = "black",
    linewidth: float = 0.15,
    show_legend: bool = True,
) -> None:
    plot_kw = dict(
        column=column,
        cmap=cmap,
        legend=show_legend,
        ax=ax,
        edgecolor=edgecolor,
        linewidth=linewidth,
        missing_kwds={"color": "lightgrey", "label": "no data"},
    )
    if show_legend:
        plot_kw["legend_kwds"] = {"label": legend_label or column, "shrink": 0.75}
    if vmin is not None:
        plot_kw["vmin"] = vmin
    if vmax is not None:
        plot_kw["vmax"] = vmax
    plot_gdf.plot(**plot_kw)
    ax.set_axis_off()


def _resolve_spatial_panel_catchments(
    output_nc_by_catchment: dict[str, dict[str, Path | str]],
    obs_csv_by_catchment: dict[str, Path | str | None],
    catchments_list: list[str],
    model_label: str,
) -> list[str]:
    panel_catchments: list[str] = []
    for catchment in catchments_list:
        csv_path = obs_csv_by_catchment.get(catchment)
        if csv_path is None or not Path(csv_path).is_file():
            continue
        path_dict = output_nc_by_catchment.get(catchment, {})
        if model_label not in path_dict:
            continue
        panel_catchments.append(catchment)
    return panel_catchments


def _render_spatial_model_obs_axes(
    axes,
    gdf_by_catchment: dict[str, gpd.GeoDataFrame],
    panel_catchments: list[str],
    *,
    var: str,
    model_row_label: str,
    obs_label: str,
    diff_label: str,
    cmap: str = "viridis",
    diff_cmap: str = "RdBu_r",
    vmin: float | None = None,
    vmax: float | None = None,
    diff_vmin: float | None = None,
    diff_vmax: float | None = None,
    edgecolor: str = "black",
    linewidth: float = 0.15,
    include_obs_row: bool = True,
    row_label_x: float = -0.12,
) -> tuple[float, float, float, float]:
    if vmin is None or vmax is None:
        if include_obs_row:
            shared_vals = np.concatenate(
                [
                    gdf_by_catchment[c][col].to_numpy(dtype=float)
                    for c in panel_catchments
                    for col in ("model_value", "obs_value")
                ]
            )
        else:
            shared_vals = np.concatenate(
                [gdf_by_catchment[c]["model_value"].to_numpy(dtype=float) for c in panel_catchments]
            )
        auto_vmin, auto_vmax = _default_value_limits(var, shared_vals)
        vmin = auto_vmin if vmin is None else vmin
        vmax = auto_vmax if vmax is None else vmax

    diff_vals = np.concatenate(
        [gdf_by_catchment[c]["diff_value"].to_numpy(dtype=float) for c in panel_catchments]
    )
    if diff_vmin is None or diff_vmax is None:
        auto_dmin, auto_dmax = _default_diff_limits(diff_vals)
        diff_vmin = auto_dmin if diff_vmin is None else diff_vmin
        diff_vmax = auto_dmax if diff_vmax is None else diff_vmax

    value_cmap = _default_cmap_for_var(var, cmap)
    value_legend = _spatial_var_legend_label(var)
    row_specs: list[tuple[str, str, str, float | None, float | None, str]] = [
        ("model_value", model_row_label, value_cmap, vmin, vmax, value_legend),
    ]
    if include_obs_row:
        row_specs.append(("obs_value", obs_label, value_cmap, vmin, vmax, value_legend))
    row_specs.append(
        ("diff_value", diff_label, diff_cmap, diff_vmin, diff_vmax, diff_label)
    )
    n_rows = len(row_specs)
    n_cols = len(panel_catchments)

    for col, catchment in enumerate(panel_catchments):
        plot_gdf = gdf_by_catchment[catchment]
        for row, (column, row_label, row_cmap, row_vmin, row_vmax, legend_label) in enumerate(
            row_specs
        ):
            ax = axes[row, col]
            _plot_gdf_column(
                plot_gdf,
                column,
                ax,
                cmap=row_cmap,
                vmin=row_vmin,
                vmax=row_vmax,
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

    return vmin, vmax, diff_vmin, diff_vmax


def plot_spatial_model_obs_panel(
    output_nc_by_catchment: dict[str, dict[str, Path | str]],
    obs_csv_by_catchment: dict[str, Path | str | None],
    *,
    model_label: str = "RF Ensemble",
    model_row_label: str = "RF ensemble mean",
    catchments_list: list[str] | None = None,
    gpep_root: Path | str = DEFAULT_GPEP_ROOT,
    var: str = wysm.VAR_SCF_DAYS,
    water_year: int | None = None,
    stat: str = "mean",
    obs_label: str = "MOD10A1",
    diff_label: str | None = None,
    gpkg_id_field: str | None = None,
    csv_id_field: str | None = None,
    figsize: tuple[float, float] | None = None,
    cmap: str = "viridis",
    diff_cmap: str = "RdBu_r",
    vmin: float | None = None,
    vmax: float | None = None,
    diff_vmin: float | None = None,
    diff_vmax: float | None = None,
    edgecolor: str = "black",
    linewidth: float = 0.15,
    include_obs_row: bool = True,
    fig: "Figure | None" = None,
    axes: np.ndarray | None = None,
    show_column_titles: bool = True,
    row_label_x: float = -0.12,
) -> tuple:
    """
    Multi-catchment choropleth panel: rows = model | obs | difference, columns = catchments.

    Set ``include_obs_row=False`` for a two-row panel (model and difference only).

    Catchments without a configured observation CSV are skipped.
    """
    _apply_spatial_publication_style()

    if diff_label is None:
        diff_label = _default_diff_label(var)

    if catchments_list is None:
        catchments_list = list(output_nc_by_catchment.keys())

    panel_catchments = _resolve_spatial_panel_catchments(
        output_nc_by_catchment,
        obs_csv_by_catchment,
        catchments_list,
        model_label,
    )
    gdf_by_catchment: dict[str, gpd.GeoDataFrame] = {}
    for catchment in panel_catchments:
        csv_path = obs_csv_by_catchment[catchment]
        path_dict = output_nc_by_catchment[catchment]
        plot_gdf, _ = build_spatial_model_obs_gdf(
            path_dict[model_label],
            csv_path,
            catchment=catchment,
            gpep_root=gpep_root,
            var=var,
            water_year=water_year,
            stat=stat,
            gpkg_id_field=gpkg_id_field,
            csv_id_field=csv_id_field,
        )
        gdf_by_catchment[catchment] = plot_gdf

    if not panel_catchments:
        raise ValueError(
            "No catchments with both model NetCDF and observation CSV were found."
        )

    n_rows = 2 if not include_obs_row else 3
    n_cols = len(panel_catchments)
    embed = axes is not None
    if embed:
        if axes.shape != (n_rows, n_cols):
            raise ValueError(
                f"Expected axes shape ({n_rows}, {n_cols}); got {axes.shape}"
            )
        if fig is None:
            fig = axes[0, 0].figure
    else:
        if figsize is None:
            figsize = (3.6 * n_cols, 2.2 * n_rows)
        fig, axes = plt.subplots(
            n_rows,
            n_cols,
            figsize=figsize,
            squeeze=False,
            constrained_layout=True,
        )

    _render_spatial_model_obs_axes(
        axes,
        gdf_by_catchment,
        panel_catchments,
        var=var,
        model_row_label=model_row_label,
        obs_label=obs_label,
        diff_label=diff_label,
        cmap=cmap,
        diff_cmap=diff_cmap,
        vmin=vmin,
        vmax=vmax,
        diff_vmin=diff_vmin,
        diff_vmax=diff_vmax,
        edgecolor=edgecolor,
        linewidth=linewidth,
        include_obs_row=include_obs_row,
        row_label_x=row_label_x,
    )

    if show_column_titles and not embed:
        fig.set_constrained_layout_pads(rect=(0, 0, 1, 0.96), h_pad=0.04, w_pad=0.04)
        fig.canvas.draw()
        top_y = max(axes[0, col].get_position().y1 for col in range(n_cols))
        title_pad = 0.012
        for col, catchment in enumerate(panel_catchments):
            pos = axes[0, col].get_position()
            fig.text(
                (pos.x0 + pos.x1) / 2,
                top_y + title_pad,
                _format_catchment_label(catchment),
                ha="center",
                va="bottom",
                fontsize=11,
                clip_on=False,
            )

    return fig, axes, gdf_by_catchment
