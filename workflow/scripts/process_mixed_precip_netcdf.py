#!/usr/bin/env python3
"""Build mixed-phase precipitation NetCDF files from SUMMA timestep outputs."""

from __future__ import annotations

import argparse
import gc
import logging
from pathlib import Path

import compute_mixed_precip_fractions as mpf
from summa_postprocess_specs import DEFAULT_GPEP_ROOT, ensemble_sources_for_catchment

logger = logging.getLogger(__name__)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Write mixed-phase precipitation NetCDF from SUMMA timestep output.",
    )
    p.add_argument(
        "--catchment",
        required=True,
        help="Catchment name (e.g. chena, bow, tuolumne).",
    )
    p.add_argument(
        "--gpep-root",
        type=Path,
        default=DEFAULT_GPEP_ROOT,
        help=f"Root data directory (default: {DEFAULT_GPEP_ROOT}).",
    )
    p.add_argument(
        "--sources",
        default="rf,casr,era5",
        help="Comma-separated source keys: rf, casr, era5 (default: all).",
    )
    p.add_argument("--glob-pattern", default="**/*.nc")
    p.add_argument("--time-join", choices=("inner", "outer"), default="inner")
    p.add_argument("--compress", action="store_true", help="Enable zlib compression.")
    p.add_argument("--no-compress", action="store_true", help="Disable zlib compression.")
    p.add_argument("--case-name-filter", default=None)
    p.add_argument("--overwrite", action="store_true")
    p.add_argument(
        "--gru-batch-size",
        type=int,
        default=None,
        help=(
            "Load/compute this many GRUs at a time from multi-GRU timestep files "
            "(default: 128). Use smaller values if memory is tight."
        ),
    )
    p.add_argument(
        "--concat-batch-size",
        type=int,
        default=None,
        help="Concatenate this many GRUs/files at a time while merging (default: 32).",
    )
    return p.parse_args(argv)


def _resolve_compress(args: argparse.Namespace) -> bool:
    if args.compress and args.no_compress:
        raise ValueError("Pass at most one of --compress and --no-compress.")
    if args.compress:
        return True
    if args.no_compress:
        return False
    return True


def process_catchment(args: argparse.Namespace) -> int:
    source_keys = [s.strip() for s in args.sources.split(",") if s.strip()]
    compress = _resolve_compress(args)

    specs = ensemble_sources_for_catchment(
        args.catchment,
        product="mixed_precip",
        gpep_root=args.gpep_root,
        source_keys=source_keys,
    )

    for spec in specs:
        out_path = Path(spec["output_nc"])
        out_path.parent.mkdir(parents=True, exist_ok=True)

        if out_path.exists() and not args.overwrite:
            logger.info("Skip (exists): %s", out_path)
            continue

        input_dir = Path(spec["input_dir"])
        if not input_dir.is_dir():
            logger.error("Input directory not found: %s", input_dir)
            return 1

        logger.info(
            "Processing %s / %s -> %s",
            args.catchment,
            spec["label"],
            out_path,
        )
        write_kwargs: dict = {
            "glob_pattern": args.glob_pattern,
            "case_name_filter": args.case_name_filter,
            "deterministic": spec["deterministic"],
            "time_join": args.time_join,
            "compress": compress,
            "overwrite": args.overwrite,
            "title": f"Mixed-phase precip — {spec['label']} ({args.catchment})",
            "source": f"SUMMA output: {input_dir}",
            "history_note": (
                f"process_mixed_precip_netcdf.py {args.catchment} {spec['label']}"
            ),
        }
        if args.gru_batch_size is not None:
            write_kwargs["gru_batch_size"] = args.gru_batch_size
        if args.concat_batch_size is not None:
            write_kwargs["concat_batch_size"] = args.concat_batch_size

        mpf.write_precip_phase_netcdf(
            out_path,
            input_dir=input_dir,
            **write_kwargs,
        )
        gc.collect()

    return 0


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    args = parse_args(argv)
    return process_catchment(args)


if __name__ == "__main__":
    raise SystemExit(main())
