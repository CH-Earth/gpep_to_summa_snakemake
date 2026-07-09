"""Publication-quality multi-panel snow metrics time-series figures."""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib.axes import Axes
from matplotlib.figure import Figure

import compute_water_year_snow_metrics as wysm
from snow_metrics_obs_csv import load_obs_csv_gru_mean_series
from snow_metrics_spatial_plot import catchment_gpkg_path, infer_gpkg_join_field
from summa_postprocess_specs import DEFAULT_GPEP_ROOT

_ELEV_FIELDS = ("elev", "elevation", "ELEV", "Elev", "avg_elev", "mean_elev", "z")

DEFAULT_VAR_ORDER = (
    wysm.VAR_SNOW_DISAPPEARANCE,
    wysm.VAR_SCF_DAYS,
)

DEFAULT_VAR_YLABELS = {
    wysm.VAR_SNOW_DISAPPEARANCE: "SDD (doy)",
    wysm.VAR_SCF_DAYS: "SCF (%)",
}

DEFAULT_SOURCE_COLORS = ("#2166ac", "#b2182b", "#4daf4a", "#984ea3")


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


def infer_gpkg_elev_field(gdf: gpd.GeoDataFrame) -> str:
    for field in _ELEV_FIELDS:
        if field in gdf.columns:
            return field
    for col in gdf.columns:
        if "elev" in col.lower():
            return col
    raise ValueError(
        f"Could not find an elevation column in GeoPackage; columns={list(gdf.columns)}"
    )


def top_elev_gru_ids(
    catchment: str,
    gru_values: np.ndarray,
    *,
    gpep_root: Path | str,
    top_fraction: float,
    elev_field: str | None = None,
    gpkg_id_field: str | None = None,
) -> np.ndarray:
    """GRU ids in the highest-elevation ``top_fraction`` of polygons (per catchment)."""
    if not 0.0 < top_fraction <= 1.0:
        raise ValueError(f"top_fraction must be in (0, 1]; got {top_fraction}")

    gpkg_path = catchment_gpkg_path(catchment, gpep_root)
    if not gpkg_path.is_file():
        raise FileNotFoundError(f"GeoPackage not found: {gpkg_path}")

    gru_values = np.asarray(gru_values, dtype=np.int64)
    gdf = gpd.read_file(gpkg_path)
    join_field = gpkg_id_field or infer_gpkg_join_field(gru_values, gdf)
    elev_col = elev_field or infer_gpkg_elev_field(gdf)

    sub = gdf[[join_field, elev_col]].copy()
    sub[join_field] = sub[join_field].astype(np.int64)
    sub[elev_col] = pd.to_numeric(sub[elev_col], errors="coerce")
    sub = sub[sub[join_field].isin(gru_values)].dropna(subset=[elev_col])
    if sub.empty:
        raise ValueError(f"No overlapping GRUs with elevation in {gpkg_path}")

    n_keep = max(1, int(np.ceil(len(sub) * top_fraction)))
    cutoff = float(sub[elev_col].nlargest(n_keep).min())
    return sub.loc[sub[elev_col] >= cutoff, join_field].astype(np.int64).values


def _resolve_top_gru_ids_by_catchment(
    output_nc_by_catchment: dict[str, dict[str, Path]],
    catchments_list: list[str],
    labels_ordered: list[str],
    *,
    top_elev_fraction: float | None,
    gpep_root: Path | str,
    elev_field: str | None,
    gpkg_id_field: str | None,
) -> dict[str, np.ndarray | None]:
    if top_elev_fraction is None:
        return {catchment: None for catchment in catchments_list}

    gru_by_catchment: dict[str, np.ndarray | None] = {}
    for catchment in catchments_list:
        ref_path = output_nc_by_catchment[catchment][labels_ordered[0]]
        with xr.open_dataset(ref_path, decode_times=False) as ds:
            if "gru" in ds.coords:
                gru_values = np.asarray(ds["gru"].values, dtype=np.int64)
            else:
                for var in ds.data_vars:
                    if "gru" in ds[var].dims:
                        gru_values = np.asarray(ds[var]["gru"].values, dtype=np.int64)
                        break
                else:
                    raise ValueError(f"{ref_path} has no gru coordinate")

        gru_by_catchment[catchment] = top_elev_gru_ids(
            catchment,
            gru_values,
            gpep_root=gpep_root,
            top_fraction=top_elev_fraction,
            elev_field=elev_field,
            gpkg_id_field=gpkg_id_field,
        )
    return gru_by_catchment


def _sel_grus(da: xr.DataArray, gru_ids: np.ndarray | None) -> xr.DataArray:
    if gru_ids is None:
        return da
    return da.sel(gru=np.asarray(gru_ids, dtype=np.int64))


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


def _sel_water_year_range(
    da: xr.DataArray,
    start_year: int | None,
    end_year: int | None,
) -> xr.DataArray:
    if start_year is None and end_year is None:
        return da
    if "water_year" not in da.dims:
        return da
    lo = start_year if start_year is not None else int(da["water_year"].min())
    hi = end_year if end_year is not None else int(da["water_year"].max())
    return da.sel(water_year=slice(lo, hi))


def _filter_obs_series(
    wy: np.ndarray,
    mean: np.ndarray,
    start_year: int | None,
    end_year: int | None,
) -> tuple[np.ndarray, np.ndarray]:
    wy = np.asarray(wy, dtype=int)
    mean = np.asarray(mean, dtype=float)
    mask = np.ones(wy.shape, dtype=bool)
    if start_year is not None:
        mask &= wy >= int(start_year)
    if end_year is not None:
        mask &= wy <= int(end_year)
    return wy[mask], mean[mask]


def _gru_mean_series(da: xr.DataArray, gru_ids: np.ndarray | None = None) -> xr.DataArray:
    if "gru" not in da.dims:
        raise ValueError(f"Expected gru dimension; got {da.dims}")
    da = _sel_grus(da, gru_ids)
    return da.mean(dim="gru", skipna=True)


def _ensemble_stats(
    da_wy_member: xr.DataArray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    wy = np.asarray(da_wy_member["water_year"].values)
    arr = np.asarray(da_wy_member.values, dtype=float)
    if "member" in da_wy_member.dims:
        member_axis = da_wy_member.dims.index("member")
        mean = np.nanmean(arr, axis=member_axis)
        p10 = np.nanpercentile(arr, 10, axis=member_axis)
        p90 = np.nanpercentile(arr, 90, axis=member_axis)
    else:
        mean = arr
        p10 = arr
        p90 = arr
    return wy, mean, p10, p90


def _display_scale(var: str) -> float:
    if var == wysm.VAR_SCF_DAYS:
        return 100.0
    return 1.0


def _default_ylim(var: str) -> tuple[float | None, float | None]:
    if var == wysm.VAR_SCF_DAYS:
        return 0.0, 100.0
    return None, None


def _style_axis(
    ax: Axes,
    *,
    start_year: int | None,
    end_year: int | None,
    show_xlabel: bool,
    show_ylabel: bool,
    ylabel: str | None,
) -> None:
    if show_ylabel and ylabel:
        ax.set_ylabel(ylabel)
    else:
        ax.set_ylabel("")
    if show_xlabel:
        ax.set_xlabel("Water year")
    else:
        ax.set_xlabel("")
        ax.tick_params(labelbottom=False)

    if start_year is not None or end_year is not None:
        lo = start_year if start_year is not None else None
        hi = end_year if end_year is not None else None
        if lo is not None and hi is not None:
            ax.set_xlim(lo - 0.5, hi + 0.5)
        elif lo is not None:
            ax.set_xlim(left=lo - 0.5)
        elif hi is not None:
            ax.set_xlim(right=hi + 0.5)

    ax.grid(True, axis="y")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _plot_panel_series(
    ax: Axes,
    opened: dict[str, xr.Dataset],
    *,
    var: str,
    labels_ordered: list[str],
    colors: tuple[str, ...],
    start_year: int | None,
    end_year: int | None,
    obs_series: tuple[np.ndarray, np.ndarray] | None,
    obs_label: str,
    gru_ids: np.ndarray | None = None,
) -> list:
    handles: list = []
    for i, label in enumerate(labels_ordered):
        da = opened[label][var]
        da = _sel_water_year_range(da, start_year, end_year)
        da_m = _gru_mean_series(da, gru_ids=gru_ids)
        wy, mean, p10, p90 = _ensemble_stats(da_m)
        scale = _display_scale(var)
        mean = mean * scale
        p10 = p10 * scale
        p90 = p90 * scale
        color = colors[i % len(colors)]
        if "member" in da_m.dims:
            ax.fill_between(wy, p10, p90, alpha=0.18, color=color, linewidth=0)
        (line,) = ax.plot(
            wy,
            mean,
            color=color,
            lw=1.8,
            marker="o",
            ms=4.5,
            markerfacecolor=color,
            markeredgecolor="white",
            markeredgewidth=0.6,
            label=label,
            zorder=3,
        )
        handles.append(line)

    if obs_series is not None:
        obs_wy, obs_mean = obs_series
        obs_mean = obs_mean * _display_scale(var)
        (obs_line,) = ax.plot(
            obs_wy,
            obs_mean,
            color="#222222",
            ls="--",
            lw=1.8,
            marker="o",
            ms=4.5,
            markerfacecolor="#222222",
            markeredgecolor="white",
            markeredgewidth=0.6,
            label=obs_label,
            zorder=4,
        )
        handles.append(obs_line)

    ymin, ymax = _default_ylim(var)
    if ymin is not None and ymax is not None:
        ax.set_ylim(ymin, ymax)

    return handles


def _labels_for_plotting(
    output_nc_by_catchment: dict[str, dict[str, Path]],
    catchments_list: list[str],
    source_labels: tuple[str, ...] | list[str] | None,
) -> list[str]:
    """Ordered source labels present in every catchment."""
    if source_labels is not None:
        labels = list(source_labels)
        for catchment in catchments_list:
            path_dict = output_nc_by_catchment[catchment]
            missing = [label for label in labels if label not in path_dict]
            if missing:
                raise KeyError(
                    f"{catchment}: missing source(s) {missing!r}; "
                    f"available {list(path_dict.keys())!r}"
                )
        return labels

    resolved: list[str] | None = None
    for catchment in catchments_list:
        keys = list(output_nc_by_catchment[catchment].keys())
        if resolved is None:
            resolved = keys
        elif keys != resolved:
            raise ValueError(
                f"Catchment {catchment!r} has sources {keys!r}; expected {resolved!r} "
                "for all catchments."
            )
    if resolved is None:
        raise ValueError("No catchments to plot.")
    return resolved


def plot_snow_metrics_publication(
    output_nc_by_catchment: dict[str, dict[str, Path]],
    *,
    vars: tuple[str, ...] | list[str] | None = None,
    var_ylabels: dict[str, str] | None = None,
    catchments_list: list[str] | None = None,
    start_year: int | None = None,
    end_year: int | None = None,
    default_start_year: int | None = None,
    default_end_year: int | None = None,
    include_obs: bool = True,
    obs_csv_by_catchment: dict[str, Path | None] | None = None,
    obs_label: str = "MOD10A1",
    source_labels: tuple[str, ...] | list[str] | None = None,
    top_elev_fraction: float | None = None,
    gpep_root: Path | str = DEFAULT_GPEP_ROOT,
    elev_field: str | None = None,
    gpkg_id_field: str | None = None,
    source_colors: tuple[str, ...] = DEFAULT_SOURCE_COLORS,
    figsize: tuple[float, float] | None = None,
    legend_ncol: int | None = None,
    fig: Figure | None = None,
    axes: np.ndarray | None = None,
    show_legend: bool = True,
    show_column_titles: bool = True,
    legend_bbox_anchor: tuple[float, float] = (0.5, -0.02),
    subplots_adjust: dict[str, float] | None = None,
) -> tuple[Figure, np.ndarray]:
    """
    Publication panel: rows = snow metrics, columns = catchments.

    MOD10A1 observations from ``obs_csv_by_catchment`` are overlaid when a CSV
    path exists for that catchment.

    Pass ``source_labels`` (e.g. ``["RF Ensemble"]``) to plot a subset of sources.
    Pass ``top_elev_fraction=0.5`` to average only the highest-elevation half of GRUs
    (elevation from ``{catchment}_tdx.gpkg``).

    Returns ``(fig, axes)`` with shape ``(n_vars, n_catchments)``.
    """
    _apply_publication_style()

    vars = tuple(vars or DEFAULT_VAR_ORDER)
    var_ylabels = {**DEFAULT_VAR_YLABELS, **(var_ylabels or {})}
    catchments_list = catchments_list or list(output_nc_by_catchment.keys())
    labels_ordered = _labels_for_plotting(
        output_nc_by_catchment,
        catchments_list,
        source_labels,
    )
    gru_ids_by_catchment = _resolve_top_gru_ids_by_catchment(
        output_nc_by_catchment,
        catchments_list,
        labels_ordered,
        top_elev_fraction=top_elev_fraction,
        gpep_root=gpep_root,
        elev_field=elev_field,
        gpkg_id_field=gpkg_id_field,
    )

    start_year, end_year = _resolve_year_range(
        start_year,
        end_year,
        default_start_year,
        default_end_year,
    )

    n_rows = len(vars)
    n_cols = len(catchments_list)
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
            figsize = (3.8 * n_cols, 2.0 * n_rows)
        fig, axes = plt.subplots(
            n_rows,
            n_cols,
            figsize=figsize,
            squeeze=False,
            sharex=True,
            sharey="row",
        )

    legend_handles: list = []
    legend_labels: list[str] = []
    seen_legend: set[str] = set()
    opened_by_var: dict[str, dict[str, dict[str, xr.Dataset]]] = {}

    try:
        for var in vars:
            opened_by_var[var] = {}
            for catchment in catchments_list:
                path_dict = output_nc_by_catchment[catchment]
                opened: dict[str, xr.Dataset] = {}
                for label in labels_ordered:
                    opened[label] = xr.open_dataset(path_dict[label], decode_times=False)
                opened_by_var[var][catchment] = opened

        for row, var in enumerate(vars):
            for col, catchment in enumerate(catchments_list):
                ax = axes[row, col]

                obs_series = None
                gru_ids = gru_ids_by_catchment[catchment]
                if include_obs and obs_csv_by_catchment:
                    csv_path = obs_csv_by_catchment.get(catchment)
                    if csv_path is not None and Path(csv_path).is_file():
                        wy, mean = load_obs_csv_gru_mean_series(
                            csv_path,
                            var,
                            gru_ids=gru_ids,
                        )
                        obs_series = _filter_obs_series(wy, mean, start_year, end_year)

                handles = _plot_panel_series(
                    ax,
                    opened_by_var[var][catchment],
                    var=var,
                    labels_ordered=labels_ordered,
                    colors=source_colors,
                    start_year=start_year,
                    end_year=end_year,
                    obs_series=obs_series,
                    obs_label=obs_label,
                    gru_ids=gru_ids,
                )

                if row == 0 and show_column_titles:
                    ax.set_title(_format_catchment_label(catchment))

                _style_axis(
                    ax,
                    start_year=start_year,
                    end_year=end_year,
                    show_xlabel=(row == n_rows - 1),
                    show_ylabel=(col == 0),
                    ylabel=var_ylabels.get(var),
                )

                for h in handles:
                    label = h.get_label()
                    if label not in seen_legend:
                        seen_legend.add(label)
                        legend_handles.append(h)
                        legend_labels.append(label)

        if legend_handles and show_legend:
            if legend_ncol is None:
                legend_ncol = min(len(legend_handles), 4)
            fig.legend(
                legend_handles,
                legend_labels,
                loc="lower center",
                bbox_to_anchor=legend_bbox_anchor,
                ncol=legend_ncol,
                frameon=False,
                handlelength=2.4,
                columnspacing=1.4,
            )

        if not embed:
            layout = subplots_adjust or {
                "left": 0.12,
                "right": 0.98,
                "top": 0.92,
                "bottom": 0.18,
                "wspace": 0.22,
                "hspace": 0.06,
            }
            fig.subplots_adjust(**layout)
        return fig, axes

    finally:
        for by_catchment in opened_by_var.values():
            for opened in by_catchment.values():
                for ds in opened.values():
                    ds.close()
