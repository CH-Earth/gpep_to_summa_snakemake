#!/usr/bin/env python3
"""Build water-year snow metrics NetCDF files from SUMMA timestep outputs."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import compute_water_year_snow_metrics as wysm
from summa_postprocess_specs import DEFAULT_GPEP_ROOT, ensemble_sources_for_catchment

logger = logging.getLogger(__name__)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Write water-year snow metrics NetCDF from SUMMA timestep output.",
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

    p.add_argument("--swe-threshold", type=float, default=1.0)
    p.add_argument("--run-days", type=int, default=5)
    p.add_argument("--expected-timesteps-per-day", type=int, default=24)
    p.add_argument("--min-daily-timesteps", type=int, default=None)
    p.add_argument("--require-complete-wy", action=argparse.BooleanOptionalAction, default=True)
    p.add_argument("--min-valid-wy-days", type=int, default=360)
    p.add_argument(
        "--require-terminal-snow-free",
        action=argparse.BooleanOptionalAction,
        default=False,
    )
    p.add_argument(
        "--allow-late-transient-snow",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    p.add_argument("--one-based-sdd", action=argparse.BooleanOptionalAction, default=True)
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
        product="snow_metrics",
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
        wysm.write_water_year_snow_metrics_netcdf(
            out_path,
            input_dir=input_dir,
            glob_pattern=args.glob_pattern,
            case_name_filter=args.case_name_filter,
            deterministic=spec["deterministic"],
            time_join=args.time_join,
            force_member_dim=True,
            require_member_dim=not spec["deterministic"],
            prefer_multi_gru_files=True,
            multi_gru_per_file=None,
            swe_threshold=args.swe_threshold,
            run_days=args.run_days,
            expected_timesteps_per_day=args.expected_timesteps_per_day,
            min_daily_timesteps=args.min_daily_timesteps,
            require_complete_wy=args.require_complete_wy,
            min_valid_wy_days=args.min_valid_wy_days,
            require_terminal_snow_free=args.require_terminal_snow_free,
            allow_late_transient_snow=args.allow_late_transient_snow,
            one_based_sdd=args.one_based_sdd,
            compress=compress,
            title=f"Water-year snow metrics — {spec['label']} ({args.catchment})",
            source=f"SUMMA output: {input_dir}",
            history_note=(
                f"process_water_year_snow_metrics_netcdf.py {args.catchment} {spec['label']}"
            ),
        )

    return 0


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    args = parse_args(argv)
    return process_catchment(args)


if __name__ == "__main__":
    raise SystemExit(main())
