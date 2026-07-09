#!/usr/bin/env python3
"""Aggregate MOD10A1.061 NDSI_Snow_Cover to HRU polygon means and write NetCDF.

This version avoids rasterio opening HDF4-EOS subdatasets directly.
Instead:
    HDF4-EOS subdataset -> temporary GeoTIFF via GDAL -> rasterio
"""

import argparse
import os
import re
import sys
import tempfile
import xml.etree.ElementTree as ET
from collections import defaultdict
from datetime import datetime, timedelta, timezone

import geopandas as gpd
import numpy as np
import rasterio
import xarray as xr
from osgeo import gdal
from rasterio.mask import mask

# -----------------------------------------------------------------------------
# Constants
# -----------------------------------------------------------------------------

MOD10_GRID = "MOD_Grid_Snow_500m"
MOD10_SNOW_FIELD = "NDSI_Snow_Cover"
MOD10_BASIC_QA_FIELD = "NDSI_Snow_Cover_Basic_QA"


# -----------------------------------------------------------------------------
# Metadata helpers
# -----------------------------------------------------------------------------

def _local_tag(tag: str) -> str:
    return tag.split("}")[-1] if "}" in tag else tag


def read_modis_hdf_files(directory: str, pattern: str = "MOD10A1") -> list[str]:
    """Return sorted MODIS HDF file paths."""
    paths = []

    for f in os.listdir(directory):
        if not f.endswith(".hdf"):
            continue
        if pattern not in f:
            continue

        paths.append(os.path.join(directory, f))

    return sorted(paths)


def granule_xml_path(hdf_path: str) -> str | None:
    for candidate in (
        hdf_path + ".xml",
        hdf_path.replace(".hdf", ".hdf.xml"),
    ):
        if os.path.isfile(candidate):
            return candidate

    return None


def parse_granule_datetime_from_xml(xml_path: str) -> datetime | None:
    try:
        root = ET.parse(xml_path).getroot()

        want = {
            "RangeBeginningDateTime",
            "SingleDateTime",
            "BeginningDateTime",
            "RangeEndingDateTime",
        }

        for el in root.iter():
            if _local_tag(el.tag) in want and el.text:
                txt = el.text.strip()

                if txt.endswith("Z"):
                    txt = txt[:-1] + "+00:00"

                return datetime.fromisoformat(txt)

        return None

    except Exception:
        return None


def extract_date_from_modis_filename(hdf_path: str) -> datetime:
    base = os.path.basename(hdf_path)

    m = re.search(r"\.A(\d{4})(\d{3})\.", base)

    if not m:
        raise ValueError(f"Cannot parse date from {base}")

    year = int(m.group(1))
    doy = int(m.group(2))

    return datetime(year, 1, 1, tzinfo=timezone.utc) + timedelta(days=doy - 1)


def granule_date_for_hdf(hdf_path: str) -> datetime:
    xml_p = granule_xml_path(hdf_path)

    d_xml = parse_granule_datetime_from_xml(xml_p) if xml_p else None
    d_fn = extract_date_from_modis_filename(hdf_path)

    if d_xml:
        if d_xml.tzinfo is None:
            d_xml = d_xml.replace(tzinfo=timezone.utc)

        return d_xml

    return d_fn


def granule_calendar_date_utc(hdf_path: str) -> datetime:
    dt = granule_date_for_hdf(hdf_path)

    d = dt.astimezone(timezone.utc).date()

    return datetime(
        d.year,
        d.month,
        d.day,
        tzinfo=timezone.utc,
    )


# -----------------------------------------------------------------------------
# GDAL helpers
# -----------------------------------------------------------------------------

def hdf4_subdataset_uri(hdf_path: str, field: str) -> str:
    return (
        f'HDF4_EOS:EOS_GRID:"{hdf_path}":'
        f"{MOD10_GRID}:{field}"
    )


def hdf4_subdataset_to_temp_tif(hdf_path: str, field: str) -> str:
    """Translate HDF4-EOS subdataset to temporary GeoTIFF."""

    sd = hdf4_subdataset_uri(hdf_path, field)

    src = gdal.Open(sd, gdal.GA_ReadOnly)

    if src is None:
        raise RuntimeError(f"GDAL could not open subdataset:\n{sd}")

    tmp = tempfile.NamedTemporaryFile(
        suffix=f"_{field}.tif",
        delete=False,
    )

    tmp.close()

    out_path = tmp.name

    out_ds = gdal.Translate(
        out_path,
        src,
        format="GTiff",
        creationOptions=[
            "COMPRESS=LZW",
            "TILED=YES",
        ],
    )

    if out_ds is None:
        raise RuntimeError(f"GDAL Translate failed:\n{sd}")

    out_ds = None
    src = None

    return out_path


# -----------------------------------------------------------------------------
# QA helpers
# -----------------------------------------------------------------------------

def _valid_mod10_snow_mask(snow: np.ndarray) -> np.ndarray:
    return (snow >= 0) & (snow <= 100)


def _valid_mod10_basic_qa_mask(
    qa: np.ndarray,
    allowed: tuple[int, ...],
) -> np.ndarray:
    return np.isin(qa, np.asarray(allowed, dtype=qa.dtype))


# -----------------------------------------------------------------------------
# HRU aggregation
# -----------------------------------------------------------------------------

def calculate_hru_mean_ndsi_snow_cover(
    hdf_path: str,
    hrus: gpd.GeoDataFrame,
    *,
    apply_basic_qa: bool = True,
    allowed_basic_qa: tuple[int, ...] = (0, 1, 2),
) -> list[tuple]:
    """Mean NDSI snow cover per HRU."""

    snow_tif = hdf4_subdataset_to_temp_tif(
        hdf_path,
        MOD10_SNOW_FIELD,
    )

    qa_tif = None

    if apply_basic_qa:
        qa_tif = hdf4_subdataset_to_temp_tif(
            hdf_path,
            MOD10_BASIC_QA_FIELD,
        )

    results = []

    qa_ds = None

    try:
        with rasterio.open(snow_tif) as snow_src:

            hrus_r = hrus.to_crs(snow_src.crs)

            if apply_basic_qa and qa_tif is not None:
                qa_ds = rasterio.open(qa_tif)

            for hru in hrus_r.itertuples():

                geom = [hru.geometry]

                snow_arr, _ = mask(
                    snow_src,
                    geom,
                    crop=True,
                    filled=False,
                )

                snow2d = snow_arr[0]

                if np.ma.isMaskedArray(snow2d):
                    snow2d = snow2d.filled(255)

                snow2d = np.asarray(snow2d)

                if apply_basic_qa and qa_ds is not None:

                    qa_arr, _ = mask(
                        qa_ds,
                        geom,
                        crop=True,
                        filled=False,
                    )

                    qa2d = qa_arr[0]

                    if np.ma.isMaskedArray(qa2d):
                        qa2d = qa2d.filled(255)

                    qa2d = np.asarray(qa2d)

                    qa_ok = _valid_mod10_basic_qa_mask(
                        qa2d,
                        allowed_basic_qa,
                    )

                else:
                    qa_ok = np.ones(snow2d.shape, dtype=bool)

                valid = _valid_mod10_snow_mask(snow2d) & qa_ok

                vals = snow2d[valid]

                if vals.size == 0:
                    results.append((hru.HRU_ID, np.nan))
                else:
                    results.append(
                        (hru.HRU_ID, float(np.mean(vals)))
                    )

    finally:

        if qa_ds is not None:
            qa_ds.close()

        for p in (snow_tif, qa_tif):

            if p is None:
                continue

            try:
                os.remove(p)
            except OSError:
                pass

    return results


def aggregate_modis_by_hru(
    hdf_paths: list[str],
    hrus: gpd.GeoDataFrame,
    **kwargs,
) -> dict:

    by_hru_day = defaultdict(lambda: defaultdict(list))

    for hdf_path in hdf_paths:

        day = granule_calendar_date_utc(hdf_path)

        means = calculate_hru_mean_ndsi_snow_cover(
            hdf_path,
            hrus,
            **kwargs,
        )

        for hru_id, val in means:
            by_hru_day[hru_id][day].append(val)

    data = {}

    for hru_id, day_map in by_hru_day.items():

        series = []

        for day in sorted(day_map):

            vals = [
                v for v in day_map[day]
                if not np.isnan(v)
            ]

            if vals:
                series.append((day, float(np.mean(vals))))
            else:
                series.append((day, np.nan))

        data[hru_id] = series

    return data


# -----------------------------------------------------------------------------
# NetCDF writing
# -----------------------------------------------------------------------------

def write_ndsi_snow_cover_netcdf(
    hru_time_data: dict,
    output_file: str,
    *,
    apply_basic_qa: bool,
    allowed_basic_qa: tuple[int, ...],
) -> None:

    hru_keys_sorted = sorted(
        hru_time_data.keys(),
        key=lambda k: str(k),
    )

    hru_ids_str = [str(k) for k in hru_keys_sorted]

    times = sorted({
        t
        for series in hru_time_data.values()
        for t, _ in series
    })

    np_time = np.array([
        np.datetime64(
            t.astimezone(timezone.utc)
             .replace(tzinfo=None)
             .isoformat()
        )
        for t in times
    ])

    arr = np.full(
        (len(hru_ids_str), len(times)),
        np.nan,
        dtype=np.float32,
    )

    for i, k in enumerate(hru_keys_sorted):

        series = hru_time_data[k]

        for t, v in series:
            j = times.index(t)
            arr[i, j] = v

    ds = xr.Dataset(
        {
            "ndsi_snow_cover": (
                ["hru", "time"],
                arr,
                {
                    "long_name": "NDSI snow cover",
                    "units": "percent",
                },
            )
        },
        coords={
            "HRU_ID": ("hru", hru_ids_str),
            "time": ("time", np_time),
        },
        attrs={
            "title": "MOD10A1 HRU snow cover",
            "source": "MOD10A1.061",
            "history": (
                f"Created {datetime.now(timezone.utc).isoformat()}"
            ),
            "mod10_basic_qa_applied": str(apply_basic_qa),
            "mod10_basic_qa_allowed_values": ",".join(
                str(x) for x in allowed_basic_qa
            ),
        },
    )

    encoding = {
        "ndsi_snow_cover": {
            "zlib": True,
            "complevel": 4,
            "_FillValue": np.nan,
        }
    }

    ds.to_netcdf(
        output_file,
        format="NETCDF4",
        encoding=encoding,
    )


# -----------------------------------------------------------------------------
# Filtering
# -----------------------------------------------------------------------------

def filter_readable_hdf(
    hdf_paths: list[str],
    *,
    apply_basic_qa: bool,
) -> list[str]:

    valid = []

    for hdf_path in hdf_paths:

        snow_sd = hdf4_subdataset_uri(
            hdf_path,
            MOD10_SNOW_FIELD,
        )

        snow_ds = gdal.Open(snow_sd)

        if snow_ds is None:
            print(
                f"Skipping unreadable snow field: {hdf_path}",
                file=sys.stderr,
            )
            continue

        if apply_basic_qa:

            qa_sd = hdf4_subdataset_uri(
                hdf_path,
                MOD10_BASIC_QA_FIELD,
            )

            qa_ds = gdal.Open(qa_sd)

            if qa_ds is None:
                print(
                    f"Skipping unreadable QA field: {hdf_path}",
                    file=sys.stderr,
                )
                continue

        valid.append(hdf_path)

    return valid


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def parse_args(argv=None):

    p = argparse.ArgumentParser(
        description=(
            "Aggregate MOD10A1 snow cover over HRUs"
        )
    )

    p.add_argument(
        "--modis-dir",
        required=True,
    )

    p.add_argument(
        "--hrus",
        required=True,
    )

    p.add_argument(
        "--output",
        required=True,
    )

    p.add_argument(
        "--filename-pattern",
        default="MOD10A1",
    )

    p.add_argument(
        "--no-basic-qa",
        action="store_true",
    )

    p.add_argument(
        "--allowed-basic-qa",
        default="0,1,2",
    )

    return p.parse_args(argv)


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main(argv=None):

    gdal.UseExceptions()

    args = parse_args(argv)

    apply_basic_qa = not args.no_basic_qa

    allowed_basic_qa = tuple(
        int(x.strip())
        for x in args.allowed_basic_qa.split(",")
    )

    hrus = gpd.read_file(args.hrus)

    hdf_paths = read_modis_hdf_files(
        args.modis_dir,
        pattern=args.filename_pattern,
    )

    if not hdf_paths:
        print(
            f"No HDF files found in {args.modis_dir}",
            file=sys.stderr,
        )
        return 1

    hdf_paths = filter_readable_hdf(
        hdf_paths,
        apply_basic_qa=apply_basic_qa,
    )

    if not hdf_paths:
        print(
            "No readable MODIS granules found.",
            file=sys.stderr,
        )
        return 1

    print(
        f"Processing {len(hdf_paths)} granules...",
        flush=True,
    )

    hru_data = aggregate_modis_by_hru(
        hdf_paths,
        hrus,
        apply_basic_qa=apply_basic_qa,
        allowed_basic_qa=allowed_basic_qa,
    )

    os.makedirs(
        os.path.dirname(os.path.abspath(args.output)) or ".",
        exist_ok=True,
    )

    write_ndsi_snow_cover_netcdf(
        hru_data,
        args.output,
        apply_basic_qa=apply_basic_qa,
        allowed_basic_qa=allowed_basic_qa,
    )

    print(f"Wrote {args.output}", flush=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())