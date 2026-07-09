"""
Build a CF-oriented NetCDF of rainfall, snowfall, snow fraction of total precip,
and a mixed-phase flag from SUMMA per-GRU timestep outputs.

Typical inputs are one NetCDF per GRU with ``*_G<id>-<id>_timestep.nc`` filenames.
Files without that pattern (e.g. ``chena_001_timestep.nc``) are treated as one or more
GRUs along a ``gru`` or ``hru`` dimension in the file; a single-site file without a
``gru`` dimension uses ``gruId`` / ``hruId`` to label ``gru``.

Designed to be called from a Jupyter notebook (no CLI).
"""


import gc
import logging
import re
import shutil
import tempfile
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, NamedTuple, Optional, Sequence, Set, Tuple, Union

import numpy as np
import xarray as xr

logger = logging.getLogger(__name__)

VAR_RAIN = "scalarRainfall_total"
VAR_SNOW = "scalarSnowfall_total"
VAR_FRAC = "snow_fraction_of_total_precip"
VAR_FLAG = "mixed_phase"

_DEFAULT_GRU_BATCH_SIZE = 128
_DEFAULT_CONCAT_BATCH_SIZE = 32

_GRU_SUFFIX = re.compile(r"_G(\d+)-(\d+)_timestep\.nc$", re.IGNORECASE)


def _concat_in_batches(
    pieces: Sequence[xr.Dataset],
    *,
    dim: str,
    batch_size: int,
    join: str = "outer",
) -> xr.Dataset:
    """Concatenate datasets in batches to limit peak memory from long piece lists."""
    if not pieces:
        raise ValueError("No datasets to concatenate.")
    if len(pieces) == 1:
        return pieces[0]
    if batch_size <= 0 or len(pieces) <= batch_size:
        return xr.concat(list(pieces), dim=dim, join=join, combine_attrs="drop")

    work: List[xr.Dataset] = list(pieces)
    while len(work) > 1:
        nxt: List[xr.Dataset] = []
        for i in range(0, len(work), batch_size):
            batch = work[i : i + batch_size]
            if len(batch) == 1:
                nxt.append(batch[0])
            else:
                nxt.append(xr.concat(batch, dim=dim, join=join, combine_attrs="drop"))
        work = nxt
    return work[0]


def _merge_along_dim_incrementally(
    blocks: Sequence[xr.Dataset],
    *,
    dim: str,
    join: str,
    concat_batch_size: int,
) -> xr.Dataset:
    """Fold datasets along ``dim`` without keeping every block in memory."""
    out: Optional[xr.Dataset] = None
    pending: List[xr.Dataset] = []
    for block in blocks:
        pending.append(block)
        if len(pending) < concat_batch_size:
            continue
        chunk = _concat_in_batches(pending, dim=dim, batch_size=concat_batch_size, join=join)
        pending.clear()
        if out is None:
            out = chunk
        else:
            merged = xr.concat([out, chunk], dim=dim, join=join, combine_attrs="drop")
            del out, chunk
            out = merged
        gc.collect()

    if pending:
        chunk = _concat_in_batches(pending, dim=dim, batch_size=concat_batch_size, join=join)
        pending.clear()
        if out is None:
            out = chunk
        else:
            merged = xr.concat([out, chunk], dim=dim, join=join, combine_attrs="drop")
            del out, chunk
            out = merged
        gc.collect()

    if out is None:
        raise ValueError("No datasets to merge.")
    return out


def classify_timestep_path_layout(paths: Sequence[Path]) -> str:
    """Return ``per_gru_filenames`` or ``multi_gru_files`` from filename patterns."""
    matched = [_GRU_SUFFIX.search(p.name) is not None for p in paths]
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
    """Infer ensemble member id from filenames like ``chena_002_timestep.nc``."""
    if deterministic:
        return None
    stem = Path(path).stem
    if stem.lower().endswith("_timestep"):
        stem = stem[: -len("_timestep")]
    if "_" in stem:
        return stem.rsplit("_", 1)[1]
    return stem if stem else "unknown"


class ParsedTimestepName(NamedTuple):
    gru_id: int
    member_id: Optional[str]


def _is_per_gru_timestep_filename(name: str) -> bool:
    return _GRU_SUFFIX.search(name) is not None


def resolve_timestep_paths(
    paths: Sequence[Path],
    *,
    prefer_multi_gru_files: bool = True,
) -> List[Path]:
    """
    Choose which timestep files to merge when more than one layout is present.

    When ``prefer_multi_gru_files`` is True and both layouts are discovered
    (e.g. ``chena_chena_timestep.nc`` plus ``*_G016-016_timestep.nc``), keep
    only the multi-GRU-in-file paths.
    """
    paths = [Path(p) for p in paths]
    if not paths or not prefer_multi_gru_files:
        return paths

    multi_gru = [p for p in paths if not _is_per_gru_timestep_filename(p.name)]
    per_gru = [p for p in paths if _is_per_gru_timestep_filename(p.name)]
    if multi_gru and per_gru:
        logger.info(
            "Found %d multi-GRU and %d per-GRU timestep files; using multi-GRU only.",
            len(multi_gru),
            len(per_gru),
        )
        return sorted(multi_gru)
    return paths


def discover_nc_files(
    input_dir: Union[str, Path],
    glob_pattern: str = "**/*.nc",
    case_name_filter: Optional[str] = None,
    prefer_multi_gru_files: bool = True,
) -> List[Path]:
    root = Path(input_dir)
    if not root.is_dir():
        raise NotADirectoryError(f"input_dir is not a directory: {root}")

    complete_dir = root / "complete"
    if prefer_multi_gru_files and complete_dir.is_dir():
        complete_paths = sorted(complete_dir.glob("*.nc"))
        if complete_paths:
            logger.info(
                "Using %d timestep file(s) from %s",
                len(complete_paths),
                complete_dir,
            )
            paths = complete_paths
        else:
            paths = sorted(root.glob(glob_pattern))
    else:
        paths = sorted(root.glob(glob_pattern))

    if case_name_filter:
        paths = [p for p in paths if case_name_filter in p.name]
    return resolve_timestep_paths(paths, prefer_multi_gru_files=prefer_multi_gru_files)


def parse_timestep_nc_path(path: Path, *, deterministic: bool) -> ParsedTimestepName:
    m = _GRU_SUFFIX.search(path.name)
    if not m:
        raise ValueError(
            f"Filename does not match expected SUMMA timestep pattern "
            f"*_G<gru>-<gru>_timestep.nc: {path.name}"
        )
    g1, g2 = int(m.group(1)), int(m.group(2))
    if g1 != g2:
        raise ValueError(f"GRU indices differ in filename (expected equal): {path.name}")
    prefix = path.name[: m.start()]
    if deterministic:
        member_id: Optional[str] = None
    else:
        if "_" in prefix:
            member_id = prefix.rsplit("_", 1)[1]
        else:
            member_id = prefix if prefix else "unknown"
    return ParsedTimestepName(gru_id=g1, member_id=member_id)


def _source_1d_values(source: xr.Dataset, name: str, dim: str) -> Optional[np.ndarray]:
    if name not in source:
        return None
    da = source[name]
    if da.dims != (dim,):
        return None
    return np.asarray(da.values)


def normalize_to_gru_dimension(
    da: xr.DataArray,
    path: Path,
    *,
    source_ds: Optional[xr.Dataset] = None,
) -> xr.DataArray:
    """Ensure a spatial dimension is named ``gru`` with integer ids when available."""
    path = Path(path)
    source = source_ds

    if "gru" in da.dims and "hru" in da.dims:
        raise ValueError(
            f"{path}: variable has both 'gru' and 'hru' dimensions; expected one spatial dim."
        )

    if "gru" in da.dims:
        out = da
        if source is not None:
            gid = _source_1d_values(source, "gruId", "gru")
            if gid is not None:
                return out.assign_coords(gru=("gru", gid.astype(np.int64)))
        if "gruId" in out.coords and out.coords["gruId"].dims == ("gru",):
            gid = np.asarray(out.coords["gruId"].values)
            return out.assign_coords(gru=("gru", gid.astype(np.int64)))
        return out

    if "hru" in da.dims:
        out = da.rename({"hru": "gru"})
        if source is not None:
            hid = _source_1d_values(source, "hruId", "hru")
            if hid is not None:
                return out.assign_coords(gru=("gru", hid.astype(np.int64)))
            gid = _source_1d_values(source, "gruId", "gru")
            if gid is not None and gid.size == out.sizes["gru"]:
                return out.assign_coords(gru=("gru", gid.astype(np.int64)))
        if "hruId" in da.coords and da.coords["hruId"].dims == ("hru",):
            hid = np.asarray(da.coords["hruId"].values)
            return out.assign_coords(gru=("gru", hid.astype(np.int64)))
        return out.assign_coords(
            gru=("gru", np.arange(out.sizes["gru"], dtype=np.int64))
        )

    return da.expand_dims(gru=[np.int64(0)])


def _squeeze_hru(da: xr.DataArray, path: Path) -> xr.DataArray:
    if "hru" not in da.dims:
        return da
    n = da.sizes.get("hru", 0)
    if n != 1:
        raise ValueError(
            f"Expected a single HRU (hru dimension size 1) in {path}, got hru={n}. "
            "HRU–GRU mapping is assumed 1:1."
        )
    return da.squeeze("hru", drop=True)


def _time_encoding_key(ds: xr.Dataset) -> Tuple[Any, Any, Any]:
    t = ds["time"]
    return (t.attrs.get("units"), t.attrs.get("calendar"), str(t.dtype))


def load_one_timestep_file(path: Union[str, Path]) -> xr.Dataset:
    path = Path(path)
    with xr.open_dataset(path) as ds:
        if VAR_RAIN not in ds or VAR_SNOW not in ds:
            raise KeyError(
                f"{path} must contain {VAR_RAIN!r} and {VAR_SNOW!r}; "
                f"found {list(ds.data_vars)}"
            )
        R = _squeeze_hru(ds[VAR_RAIN].load(), path)
        S = _squeeze_hru(ds[VAR_SNOW].load(), path)
    return _build_precip_phase_dataset(R, S)


def _build_precip_phase_dataset(R: xr.DataArray, S: xr.DataArray) -> xr.Dataset:
    if R.dims != S.dims or not np.array_equal(R["time"].values, S["time"].values):
        raise ValueError("Rain and snow time series differ in shape or time")

    total = R + S
    snow_frac = xr.where(total > 0.0, S / total, np.nan).astype(np.float32)
    mixed = ((R > 0.0) & (S > 0.0)).astype(np.int8)

    out = xr.Dataset(
        {
            VAR_RAIN: R.astype(np.float32),
            VAR_SNOW: S.astype(np.float32),
            VAR_FRAC: snow_frac,
            VAR_FLAG: mixed,
        }
    )
    out[VAR_RAIN].attrs.update({k: v for k, v in R.attrs.items() if k != "_FillValue"})
    out[VAR_SNOW].attrs.update({k: v for k, v in S.attrs.items() if k != "_FillValue"})
    return out


def load_one_timestep_multi_gru_file(
    path: Union[str, Path],
    *,
    gru_batch_size: Optional[int] = None,
) -> xr.Dataset:
    """Load rain/snow from a file that may contain many GRUs/HRUs along ``gru``."""
    path = Path(path)
    batch_size = gru_batch_size if gru_batch_size is not None else _DEFAULT_GRU_BATCH_SIZE
    with xr.open_dataset(path) as ds:
        if VAR_RAIN not in ds or VAR_SNOW not in ds:
            raise KeyError(
                f"{path} must contain {VAR_RAIN!r} and {VAR_SNOW!r}; "
                f"found {list(ds.data_vars)}"
            )
        rain_da = normalize_to_gru_dimension(ds[VAR_RAIN], path, source_ds=ds)
        snow_da = normalize_to_gru_dimension(ds[VAR_SNOW], path, source_ds=ds)
        n_gru = int(rain_da.sizes.get("gru", 1))

        if batch_size <= 0 or n_gru <= batch_size:
            R = rain_da.load()
            S = snow_da.load()
            return _build_precip_phase_dataset(R, S)

        chunks: List[xr.Dataset] = []
        for start in range(0, n_gru, batch_size):
            end = min(start + batch_size, n_gru)
            R = rain_da.isel(gru=slice(start, end)).load()
            S = snow_da.isel(gru=slice(start, end)).load()
            chunks.append(_build_precip_phase_dataset(R, S))
            del R, S
        return _concat_in_batches(
            chunks,
            dim="gru",
            batch_size=max(1, batch_size // 4),
            join="outer",
        )


def _assert_same_gru_coord(a: xr.Dataset, b: xr.Dataset, pa: Path, pb: Path) -> None:
    ga = a[VAR_RAIN].coords["gru"].values
    gb = b[VAR_RAIN].coords["gru"].values
    if ga.shape != gb.shape or not np.array_equal(ga, gb):
        raise ValueError(
            "When joining multi-GRU files along time, the gru coordinate must match. "
            f"Mismatch between {pa!s} and {pb!s}."
        )


def _sort_by_time_if_present(ds: xr.Dataset) -> xr.Dataset:
    if "time" not in ds.dims and "time" not in ds.coords:
        return ds
    try:
        return ds.sortby("time")
    except Exception:
        logger.warning("Could not sort dataset by time.", exc_info=True)
        return ds


def _build_multi_gru_member_dataset(
    plist: Sequence[Path],
    *,
    member_id: Optional[str],
    time_join: str,
    gru_batch_size: Optional[int],
    concat_batch_size: int,
    reference_time_encoding: Optional[Tuple[Any, Any, Any]],
) -> Tuple[xr.Dataset, Tuple[Any, Any, Any]]:
    ref_enc = reference_time_encoding
    pieces: List[xr.Dataset] = []
    for filepath in plist:
        piece = load_one_timestep_multi_gru_file(filepath, gru_batch_size=gru_batch_size)
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
        ds_member = _concat_in_batches(
            pieces,
            dim="time",
            batch_size=concat_batch_size,
            join=time_join,
        )

    ds_member = _sort_by_time_if_present(ds_member)
    if member_id is not None:
        ds_member = _expand_member(ds_member, member_id)
    return ds_member, ref_enc if ref_enc is not None else _time_encoding_key(ds_member)


def merge_timestep_files_multi_gru_in_file(
    paths: Sequence[Path],
    *,
    deterministic: bool = False,
    time_join: str = "outer",
    force_member_dim: bool = False,
    reference_time_encoding: Optional[Tuple[Any, Any, Any]] = None,
    gru_batch_size: Optional[int] = None,
    concat_batch_size: int = _DEFAULT_CONCAT_BATCH_SIZE,
    cleanup_paths: Optional[List[Path]] = None,
) -> xr.Dataset:
    """Merge files like ``chena_001_timestep.nc`` where each file holds many GRUs."""
    if time_join not in ("inner", "outer"):
        raise ValueError("time_join must be 'inner' or 'outer'")

    ref_enc = reference_time_encoding
    by_member: Dict[Optional[str], List[Path]] = defaultdict(list)
    for p in paths:
        by_member[infer_member_label_from_path(p, deterministic=deterministic)].append(Path(p))

    member_ids = sorted(by_member.keys(), key=lambda m: (m is None, str(m)))
    if not member_ids:
        raise ValueError("No input files to merge.")

    if len(member_ids) > 1:
        temp_dir = Path(tempfile.mkdtemp(prefix="mixed_precip_members_"))
        if cleanup_paths is not None:
            cleanup_paths.append(temp_dir)
        temp_paths: List[Path] = []
        for member_id in member_ids:
            plist = sorted(by_member[member_id])
            ds_member, ref_enc = _build_multi_gru_member_dataset(
                plist,
                member_id=member_id,
                time_join=time_join,
                gru_batch_size=gru_batch_size,
                concat_batch_size=concat_batch_size,
                reference_time_encoding=ref_enc,
            )
            label = str(member_id) if member_id is not None else "det"
            temp_path = temp_dir / f"member_{label}.nc"
            ds_member.to_netcdf(temp_path, encoding=build_netcdf_encoding(ds_member, compress=False))
            temp_paths.append(temp_path)
            del ds_member
            gc.collect()

        out = xr.open_mfdataset(
            [str(p) for p in temp_paths],
            combine="nested",
            concat_dim="member",
            join=time_join,
        )
    else:
        member_id = member_ids[0]
        plist = sorted(by_member[member_id])
        out, _ = _build_multi_gru_member_dataset(
            plist,
            member_id=member_id,
            time_join=time_join,
            gru_batch_size=gru_batch_size,
            concat_batch_size=concat_batch_size,
            reference_time_encoding=ref_enc,
        )

    out = _sort_by_time_if_present(out)
    n_member = out.sizes.get("member", 1)
    if n_member == 1 and not force_member_dim and "member" in out.dims:
        out = out.squeeze("member", drop=True)

    _sanity_check_output(out, check_values=not _dataset_is_dask_backed(out))
    logger.info("Merged multi-GRU-in-file dataset: dims=%s, n_files=%s", dict(out.sizes), len(paths))
    return out


def _expand_gru(ds: xr.Dataset, gru_id: int) -> xr.Dataset:
    return ds.expand_dims(gru=[int(gru_id)])


def _expand_member(ds: xr.Dataset, member_id: str) -> xr.Dataset:
    return ds.expand_dims(member=[str(member_id)])


def merge_timestep_files_to_dataset(
    paths: Sequence[Path],
    *,
    deterministic: bool = False,
    time_join: str = "outer",
    force_member_dim: bool = False,
    reference_time_encoding: Optional[Tuple[Any, Any, Any]] = None,
    multi_gru_per_file: Optional[bool] = None,
    gru_batch_size: Optional[int] = None,
    concat_batch_size: int = _DEFAULT_CONCAT_BATCH_SIZE,
    cleanup_paths: Optional[List[Path]] = None,
) -> xr.Dataset:
    if time_join not in ("inner", "outer"):
        raise ValueError("time_join must be 'inner' or 'outer'")

    paths = [Path(p) for p in paths]
    if multi_gru_per_file is None:
        multi_gru_per_file = classify_timestep_path_layout(paths) == "multi_gru_files"
    if multi_gru_per_file:
        return merge_timestep_files_multi_gru_in_file(
            paths,
            deterministic=deterministic,
            time_join=time_join,
            force_member_dim=force_member_dim,
            reference_time_encoding=reference_time_encoding,
            gru_batch_size=gru_batch_size,
            concat_batch_size=concat_batch_size,
            cleanup_paths=cleanup_paths,
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
        for path, gru_id in items:
            piece = load_one_timestep_file(path)
            if ref_enc is None:
                ref_enc = _time_encoding_key(piece)
            else:
                cur = _time_encoding_key(piece)
                if cur != ref_enc:
                    logger.warning(
                        "Time encoding differs from first file: %s vs ref %s (%s)",
                        cur,
                        ref_enc,
                        path,
                    )
            gru_pieces.append(_expand_gru(piece, gru_id))
            if len(gru_pieces) >= concat_batch_size:
                merged = _concat_in_batches(
                    gru_pieces,
                    dim="gru",
                    batch_size=concat_batch_size,
                    join=time_join,
                )
                gru_pieces = [merged]
                gc.collect()

        ds_member = _concat_in_batches(
            gru_pieces,
            dim="gru",
            batch_size=concat_batch_size,
            join=time_join,
        )
        if member_id is not None:
            ds_member = _expand_member(ds_member, member_id)
        member_blocks.append(ds_member)
        gc.collect()

    if not member_blocks:
        raise ValueError("No input files to merge.")

    if len(member_blocks) == 1:
        out = member_blocks[0]
    else:
        out = _merge_along_dim_incrementally(
            member_blocks,
            dim="member",
            join=time_join,
            concat_batch_size=max(2, concat_batch_size // 4),
        )

    n_member = out.sizes.get("member", 1)
    if n_member == 1 and not force_member_dim and "member" in out.dims:
        out = out.squeeze("member", drop=True)

    _sanity_check_output(out)
    logger.info("Merged dataset: dims=%s, n_files=%s", dict(out.sizes), len(paths))
    return out


def _dataset_is_dask_backed(ds: xr.Dataset) -> bool:
    return any(hasattr(ds[name].data, "compute") for name in ds.data_vars)


def _sanity_check_output(ds: xr.Dataset, *, check_values: bool = True) -> None:
    if "hru" in ds.dims:
        raise AssertionError("Output must not contain an 'hru' dimension.")
    for v in (VAR_RAIN, VAR_SNOW, VAR_FRAC, VAR_FLAG):
        if v not in ds:
            raise AssertionError(f"Missing variable {v!r}")
    dims0 = ds[VAR_RAIN].dims
    for v in (VAR_SNOW, VAR_FRAC, VAR_FLAG):
        if ds[v].dims != dims0:
            raise AssertionError(
                f"Dimension mismatch: {VAR_RAIN} {dims0} vs {v} {ds[v].dims}"
            )
    sf = ds[VAR_FRAC]
    if not check_values:
        return
    if hasattr(sf.data, "compute"):
        lo = float(sf.min().compute())
        hi = float(sf.max().compute())
    else:
        arr = sf.values
        mask = np.isfinite(arr)
        if not mask.any():
            return
        lo, hi = float(np.nanmin(arr[mask])), float(np.nanmax(arr[mask]))
    if lo < -1e-4 or hi > 1.0 + 1e-4:
        logger.warning(
            "snow_fraction_of_total_precip outside [0, 1] in places (min=%s max=%s)",
            lo,
            hi,
        )


def attach_cf_metadata(
    ds: xr.Dataset,
    *,
    title: str = "Mixed-phase precipitation from SUMMA",
    source: str = "SUMMA timestep output",
    history_note: str = "",
    institution: str = "",
    copy: bool = False,
) -> xr.Dataset:
    if copy:
        ds = ds.copy()
    hist = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    if history_note:
        hist = f"{hist} {history_note}"
    ds.attrs["Conventions"] = "CF-1.8"
    ds.attrs["title"] = title
    ds.attrs["source"] = source
    ds.attrs["history"] = hist
    if institution:
        ds.attrs["institution"] = institution
    ds.attrs["mixed_precip_note"] = (
        f"{VAR_FRAC} is snowfall / (rainfall + snowfall); missing where total is zero. "
        f"{VAR_FLAG} is 1 where rainfall and snowfall are both strictly positive."
    )

    t = ds["time"]
    if "standard_name" not in t.attrs:
        t.attrs["standard_name"] = "time"

    ds[VAR_FRAC].attrs.setdefault(
        "long_name",
        "Snow fraction of total modeled precipitation (snowfall / (rainfall + snowfall))",
    )
    ds[VAR_FRAC].attrs.setdefault("units", "1")

    ds[VAR_FLAG].attrs["long_name"] = (
        "Mixed-phase precipitation (rain and snow simultaneously, strictly positive)"
    )
    ds[VAR_FLAG].attrs["flag_values"] = np.array([0, 1], dtype=np.int8)
    ds[VAR_FLAG].attrs["flag_meanings"] = "not_mixed mixed"

    return ds


def build_netcdf_encoding(ds: xr.Dataset, *, compress: bool = False) -> Dict[str, Dict[str, Any]]:
    """Encoding for to_netcdf. NaNs in floats use IEEE representation under NetCDF4."""
    encoding: Dict[str, Dict[str, Any]] = {}
    for name, da in ds.data_vars.items():
        enc: Dict[str, Any] = {}
        if compress:
            enc["zlib"] = True
            enc["complevel"] = 4
        if name == VAR_FLAG:
            enc["dtype"] = "int8"
        elif np.issubdtype(da.dtype, np.floating):
            enc["dtype"] = "float32"
        encoding[name] = enc
    return encoding


def write_precip_phase_netcdf(
    output_path: Union[str, Path],
    *,
    input_dir: Union[str, Path],
    glob_pattern: str = "**/*.nc",
    case_name_filter: Optional[str] = None,
    deterministic: bool = False,
    time_join: str = "outer",
    force_member_dim: bool = False,
    prefer_multi_gru_files: bool = True,
    multi_gru_per_file: Optional[bool] = None,
    gru_batch_size: Optional[int] = None,
    concat_batch_size: int = _DEFAULT_CONCAT_BATCH_SIZE,
    compress: bool = False,
    overwrite: bool = False,
    return_dataset: bool = False,
    title: str = "Mixed-phase precipitation from SUMMA",
    source: str = "SUMMA timestep output",
    history_note: str = "",
    institution: str = "",
) -> Optional[xr.Dataset]:
    out_path = Path(output_path)
    if out_path.exists() and not overwrite:
        logger.info("Output already exists, skipping: %s", out_path)
        if return_dataset:
            with xr.open_dataset(out_path) as ds:
                return ds.load()
        return None

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

    cleanup_paths: List[Path] = []
    ds = merge_timestep_files_to_dataset(
        paths,
        deterministic=deterministic,
        time_join=time_join,
        force_member_dim=force_member_dim,
        multi_gru_per_file=multi_gru_per_file,
        gru_batch_size=gru_batch_size,
        concat_batch_size=concat_batch_size,
        cleanup_paths=cleanup_paths,
    )
    try:
        ds = attach_cf_metadata(
            ds,
            title=title,
            source=source,
            history_note=history_note,
            institution=institution,
        )
        out_path.parent.mkdir(parents=True, exist_ok=True)
        encoding = build_netcdf_encoding(ds, compress=compress)
        ds.to_netcdf(out_path, encoding=encoding)
        logger.info("Wrote %s", out_path)
        if return_dataset:
            return ds.load()
        return None
    finally:
        if hasattr(ds, "close"):
            try:
                ds.close()
            except Exception:
                pass
        for path in cleanup_paths:
            if path.is_dir():
                shutil.rmtree(path, ignore_errors=True)
            elif path.is_file():
                path.unlink(missing_ok=True)
        del ds
        gc.collect()