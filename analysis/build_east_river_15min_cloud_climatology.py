#!/usr/bin/env python3
"""Calculate pixelwise 15-minute cloud-frequency climatologies."""

from __future__ import annotations

import argparse
import calendar
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


EXPECTED_HOURS = {0, 1, *range(14, 24)}


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument(
        "--mask-dir", type=Path,
        default=Path("/glade/derecho/scratch/cdalden/east_river_goes/goes16/cloud_mask_tempbin_tree"),
    )
    p.add_argument(
        "--output-dir", type=Path,
        default=Path("/glade/derecho/scratch/cdalden/east_river_goes/climatology_15min"),
    )
    p.add_argument("--start", default="2017-04")
    p.add_argument("--end", default="2024-12")
    p.add_argument("--domain-token", default="east_river")
    p.add_argument("--output-prefix", default="east_river")
    p.add_argument("--title-domain", default="East River")
    p.add_argument(
        "--include-partial-months", action="store_true",
        help="Use all available accepted days instead of requiring every calendar day.",
    )
    p.add_argument("--skip-existing", action="store_true")
    args = p.parse_args()
    pattern = re.compile(
        rf"goes16_C02_C05_C13_rgb_{re.escape(args.domain_token)}_(\d{{8}})_cloud_binary_tempbin_tree\.nc$"
    )

    files_by_month: dict[tuple[int, int], dict[int, Path]] = defaultdict(dict)
    for path in args.mask_dir.glob("*.nc"):
        match = pattern.match(path.name)
        if match:
            date = pd.Timestamp(match.group(1))
            files_by_month[(date.year, date.month)][date.day] = path

    complete: dict[int, list[tuple[int, int, list[Path]]]] = defaultdict(list)
    for period in pd.period_range(args.start, args.end, freq="M"):
        ndays = calendar.monthrange(period.year, period.month)[1]
        daily = files_by_month.get((period.year, period.month), {})
        if daily and (args.include_partial_months or set(daily) == set(range(1, ndays + 1))):
            complete[period.month].append(
                (period.year, len(daily), [daily[day] for day in sorted(daily)])
            )

    n_complete = sum(len(items) for items in complete.values())
    print(f"Found {n_complete} usable year-months", flush=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    intervals = np.array(
        [hour * 60 + quarter * 15 for hour in sorted(EXPECTED_HOURS) for quarter in range(4)],
        dtype=np.int16,
    )

    for calendar_month in range(1, 13):
        output = args.output_dir / f"{args.output_prefix}_cloud_frequency_15min_climatology_month_{calendar_month:02d}.nc"
        if args.skip_existing and output.is_file() and output.stat().st_size > 0:
            print(f"REUSE {output}", flush=True)
            continue
        groups = complete.get(calendar_month, [])
        if not groups:
            print(f"SKIP calendar month {calendar_month:02d}: no complete inputs", flush=True)
            continue
        first_path = groups[0][2][0]
        with xr.open_dataset(first_path) as first:
            lat = np.asarray(first.latitude, dtype=np.float64)
            lon = np.asarray(first.longitude, dtype=np.float64)
        shape = (len(intervals), len(lat), len(lon))
        cloudy = np.zeros(shape, dtype=np.uint32)
        valid = np.zeros(shape, dtype=np.uint32)
        interval_lookup = {int(value): index for index, value in enumerate(intervals)}

        for year, _, paths in groups:
            print(f"Month {calendar_month:02d}: reading {year} ({len(paths)} days)", flush=True)
            for path in paths:
                with xr.open_dataset(path) as ds:
                    if ds.attrs.get("tree_logic") != "AND within each leaf; OR across cloudy leaves":
                        raise ValueError(f"Wrong mask provenance: {path}")
                    times = pd.DatetimeIndex(pd.to_datetime(ds.t.values))
                    values = np.asarray(ds.cloud_binary.values)
                    if values.shape[1:] != shape[1:]:
                        raise ValueError(f"Grid mismatch: {path}")
                    for interval_minute in np.unique(times.hour * 60 + (times.minute // 15) * 15):
                        if int(interval_minute) not in interval_lookup:
                            continue
                        index = interval_lookup[int(interval_minute)]
                        selected = (times.hour * 60 + (times.minute // 15) * 15) == interval_minute
                        subset = values[selected]
                        finite = np.isfinite(subset)
                        valid[index] += finite.sum(axis=0, dtype=np.uint32)
                        cloudy[index] += ((subset == 1) & finite).sum(axis=0, dtype=np.uint32)

        frequency = np.divide(
            cloudy, valid, out=np.full(shape, np.nan, dtype=np.float32), where=valid > 0
        )
        source_months = [f"{year}-{calendar_month:02d}" for year, _, _ in groups]
        source_days = sum(n_days for _, n_days, _ in groups)
        out = xr.Dataset(
            data_vars={
                "cloud_frequency": (("interval", "latitude", "longitude"), frequency),
                "cloudy_observation_count": (("interval", "latitude", "longitude"), cloudy),
                "valid_observation_count": (("interval", "latitude", "longitude"), valid),
            },
            coords={
                "interval": np.arange(len(intervals), dtype=np.int16),
                "interval_start_minute_utc": ("interval", intervals),
                "latitude": lat,
                "longitude": lon,
            },
            attrs={
                "title": f"{args.title_domain} pixelwise 15-minute cloud-frequency climatology",
                "calendar_month": calendar_month,
                "interval_definition": "left-inclusive UTC bins, e.g. 1080=[18:00,18:15)",
                "source_complete_year_months": ",".join(source_months),
                "source_year_months": ",".join(source_months),
                "source_accepted_days": source_days,
                "partial_months_included": str(args.include_partial_months),
                "n_complete_year_months": len(groups),
                "mask_logic": "AND within each decision-tree leaf; OR across cloudy leaves",
            },
        )
        out.cloud_frequency.attrs.update(long_name="cloud occurrence frequency", units="1")
        temporary = output.with_suffix(".nc.part")
        out.to_netcdf(
            temporary,
            encoding={
                "cloud_frequency": {"zlib": True, "complevel": 4, "dtype": "float32"},
                "cloudy_observation_count": {"zlib": True, "complevel": 4, "dtype": "uint32"},
                "valid_observation_count": {"zlib": True, "complevel": 4, "dtype": "uint32"},
            },
        )
        temporary.replace(output)
        print(f"WROTE {output} from {len(groups)} complete year-months", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
