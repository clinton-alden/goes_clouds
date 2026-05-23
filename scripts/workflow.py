"""Small public API for the GOES VINTAGE workflow."""

from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

import pandas as pd
from tqdm.auto import tqdm
import xarray as xr


@dataclass
class WorkflowConfig:
    domain: str
    goes: str
    start_date: str
    end_date: str
    lon_min: float
    lat_min: float
    lon_max: float
    lat_max: float
    base_dir: Path
    goes_hours: str = "20"
    goes_timesteps_per_hour: int | None = None
    start_hour_utc: float | None = None
    end_hour_utc: float | None = None
    threshold_csv: Path | None = None
    overwrite_mask: bool = True
    keep_mask_diagnostics: bool = False
    run: bool = True

    def __post_init__(self) -> None:
        self.base_dir = Path(self.base_dir)
        if self.threshold_csv is not None:
            self.threshold_csv = Path(self.threshold_csv)
        if self.start_hour_utc is None or self.end_hour_utc is None:
            start_hour, end_hour = _mask_window_from_goes_hours(self.goes_hours)
            if self.start_hour_utc is None:
                self.start_hour_utc = start_hour
            if self.end_hour_utc is None:
                self.end_hour_utc = end_hour

    @property
    def repo_dir(self) -> Path:
        return Path(__file__).resolve().parents[1]

    @property
    def scripts_dir(self) -> Path:
        return self.repo_dir / "scripts"

    @property
    def rgb_dir(self) -> Path:
        return self.base_dir / self.goes / "rgb_composite"

    @property
    def era5_dir(self) -> Path:
        return self.base_dir / "era5_land" / "t2m_hourly"

    @property
    def mask_dir(self) -> Path:
        return self.base_dir / self.goes / "vintage_mask"

    @property
    def gif_dir(self) -> Path:
        return self.base_dir / self.goes / "vintage_gif_loops"

    @property
    def resolved_threshold_csv(self) -> Path:
        if self.threshold_csv is not None:
            return self.threshold_csv
        return self.repo_dir / "thresholds" / "gothic_vintage_rgb_tree_rules_5c_sw_kt050_090.csv"

    @property
    def dates(self) -> pd.DatetimeIndex:
        return pd.date_range(self.start_date, self.end_date, freq="D")

    def date_ymd(self, timestamp: pd.Timestamp) -> str:
        return timestamp.strftime("%Y%m%d")

    def rgb_path(self, timestamp: pd.Timestamp) -> Path:
        date = self.date_ymd(timestamp)
        return self.rgb_dir / f"{self.goes}_C02_C05_C13_rgb_{self.domain}_{date}.nc"

    def mask_path(self, timestamp: pd.Timestamp) -> Path:
        rgb_path = self.rgb_path(timestamp)
        return self.mask_dir / f"{rgb_path.stem}_vintage_mask.nc"

    def env(self) -> dict[str, str]:
        env = os.environ.copy()
        env.update(
            {
                "GOES_HOURS": self.goes_hours,
                "LON_MIN": str(self.lon_min),
                "LAT_MIN": str(self.lat_min),
                "LON_MAX": str(self.lon_max),
                "LAT_MAX": str(self.lat_max),
                "PYTHONPATH": f"{self.scripts_dir}:{env.get('PYTHONPATH', '')}",
            }
        )
        if self.goes_timesteps_per_hour is not None:
            env["GOES_TIMESTEPS_PER_HOUR"] = str(self.goes_timesteps_per_hour)
        return env


def validate_credentials(config: WorkflowConfig) -> None:
    """Check the credentials that are required once downloads begin."""
    if not config.run:
        _status("Dry run: credential checks skipped")
        return
    if not os.environ.get("OPENTOPOGRAPHY_API_KEY"):
        raise RuntimeError(
            "Set your own OPENTOPOGRAPHY_API_KEY before orthorectifying GOES data."
        )
    has_cds_env = os.environ.get("CDSAPI_URL") and os.environ.get("CDSAPI_KEY")
    has_cds_rc = os.environ.get("CDSAPI_RC") or (Path.home() / ".cdsapirc").exists()
    if not (has_cds_env or has_cds_rc):
        raise RuntimeError(
            "Missing CDS credentials. Copy example_private_env.sh to private_env.sh, "
            "set CDSAPI_URL/CDSAPI_KEY, and source it before running."
        )


def download_goes(config: WorkflowConfig) -> list[Path]:
    """Download GOES C02/C05/C13 files for the configured date range."""
    _ensure_dirs(config)
    download_script = config.scripts_dir / "download-goes.py"
    outputs: list[Path] = []
    tasks = [(date, channel) for date in config.dates for channel in ("C02", "C05", "C13")]
    hours = _parse_goes_hours(config.goes_hours)
    timesteps = config.goes_timesteps_per_hour or 12
    expected_files = len(tasks) * len(hours) * timesteps
    _status(
        f"Download request: {len(config.dates)} day(s), 3 channels, "
        f"{len(hours)} inclusive UTC hour(s), roughly {expected_files} files"
    )
    for date, channel in tqdm(tasks, desc="download", unit="channel"):
        year, month, day = date.year, date.month, date.day
        _status(f"Downloading {config.goes} {channel} for {date:%Y-%m-%d}")
        _run(
            [
                sys.executable,
                download_script,
                "-B",
                f"noaa-{config.goes}",
                "-Y",
                year,
                "-M",
                month,
                "-D",
                day,
                day,
                "-p",
                "ABI-L1b-RadC",
                "-c",
                channel,
                "-b",
                config.lon_min,
                config.lat_min,
                config.lon_max,
                config.lat_max,
                "-d",
                config.base_dir,
            ],
            config,
        )
        outputs.append(config.base_dir / config.goes / str(year) / str(month) / str(day))
    _status("Download step complete")
    return outputs


def orthorectify(config: WorkflowConfig) -> None:
    """Orthorectify downloaded GOES NetCDF files."""
    ortho_script = config.scripts_dir / "batch_ortho.py"
    for date in tqdm(config.dates, desc="ortho", unit="day"):
        month_dir = config.base_dir / config.goes / str(date.year) / str(date.month)
        pending = _count_pending_ortho_inputs(month_dir)
        _status(f"Orthorectifying {date:%Y-%m-%d}: {pending} file(s) pending")
        _run_with_ortho_progress(
            [sys.executable, ortho_script, month_dir, config.domain],
            config,
            month_dir,
            pending,
        )
    _status("Orthorectification step complete")


def build_zarr(config: WorkflowConfig) -> None:
    """Convert daily orthorectified NetCDF files to channel Zarr stores."""
    zarr_script = config.scripts_dir / "zarr_v2_days.py"
    for date in tqdm(config.dates, desc="zarr", unit="day"):
        _status(f"Building Zarr stores for {date:%Y-%m-%d}")
        _run(
            [
                sys.executable,
                zarr_script,
                str(config.base_dir) + "/",
                date.year,
                date.month,
                date.day,
                date.day,
                config.goes,
                config.domain,
            ],
            config,
        )
        if config.run:
            _dedupe_zarr_timestamps(config, date)
    _status("Zarr step complete")


def build_rgb(config: WorkflowConfig) -> list[Path]:
    """Build daily RGB NetCDF files from C02/C05/C13 Zarr stores."""
    if str(config.scripts_dir) not in sys.path:
        sys.path.insert(0, str(config.scripts_dir))
    import utils

    outputs = []
    for date in tqdm(config.dates, desc="rgb", unit="day"):
        _status(f"Building RGB for {date:%Y-%m-%d}")
        if config.run:
            utils.goes_rad_to_rgb(
                str(config.base_dir / config.goes) + "/",
                config.date_ymd(date),
                config.goes,
                location=config.domain,
                complete_day=config.goes_timesteps_per_hour is None,
            )
        outputs.append(config.rgb_path(date))
    _status("RGB step complete")
    return outputs


def apply_mask(config: WorkflowConfig) -> list[Path]:
    """Download/reuse ERA5-Land temperature and apply the GOES VINTAGE mask."""
    mask_script = config.scripts_dir / "apply_tempbin_thresholds.py"
    outputs = []
    for date in tqdm(config.dates, desc="mask", unit="day"):
        rgb_path = config.rgb_path(date)
        _status(f"Applying VINTAGE mask for {date:%Y-%m-%d}")
        cmd = [
            sys.executable,
            mask_script,
            "--rgb-file",
            rgb_path,
            "--threshold-csv",
            config.resolved_threshold_csv,
            "--era5-dir",
            config.era5_dir,
            "--mask-dir",
            config.mask_dir,
            "--gif-dir",
            config.gif_dir,
            "--domain",
            config.domain,
            "--lon-min",
            config.lon_min,
            "--lat-min",
            config.lat_min,
            "--lon-max",
            config.lon_max,
            "--lat-max",
            config.lat_max,
            "--start-hour-utc",
            config.start_hour_utc,
            "--end-hour-utc",
            config.end_hour_utc,
        ]
        if config.overwrite_mask:
            cmd.append("--overwrite")
        if config.keep_mask_diagnostics:
            cmd.append("--keep-diagnostics")
        _run(cmd, config)
        outputs.append(config.mask_path(date))
    _status("Mask step complete")
    return outputs


def _ensure_dirs(config: WorkflowConfig) -> None:
    for path in (config.base_dir, config.era5_dir, config.mask_dir, config.gif_dir):
        path.mkdir(parents=True, exist_ok=True)


def _run(cmd: list[object], config: WorkflowConfig) -> None:
    cmd = [str(part) for part in cmd]
    if not config.run:
        _status("Dry run: " + " ".join(cmd))
        return
    result = subprocess.run(
        cmd,
        cwd=config.repo_dir,
        env=config.env(),
        text=True,
        capture_output=True,
    )
    if result.returncode != 0:
        tail = "\n".join((result.stdout + "\n" + result.stderr).splitlines()[-20:])
        raise RuntimeError(f"Command failed: {' '.join(cmd)}\n\nLast output lines:\n{tail}")


def _run_with_ortho_progress(
    cmd: list[object], config: WorkflowConfig, month_dir: Path, total_files: int
) -> None:
    cmd = [str(part) for part in cmd]
    if not config.run:
        _status("Dry run: " + " ".join(cmd))
        return
    if total_files == 0:
        _status("Orthorectification already complete")
        return

    with tempfile.TemporaryFile(mode="w+t") as log:
        process = subprocess.Popen(
            cmd,
            cwd=config.repo_dir,
            env=config.env(),
            text=True,
            stdout=log,
            stderr=subprocess.STDOUT,
        )
        with tqdm(total=total_files, desc="ortho files", unit="file", leave=False) as progress:
            while process.poll() is None:
                completed = total_files - _count_pending_ortho_inputs(month_dir)
                progress.update(max(0, completed - progress.n))
                time.sleep(2)
            completed = total_files - _count_pending_ortho_inputs(month_dir)
            progress.update(max(0, completed - progress.n))

        if process.returncode != 0:
            log.seek(0)
            tail = "\n".join(log.read().splitlines()[-20:])
            raise RuntimeError(f"Command failed: {' '.join(cmd)}\n\nLast output lines:\n{tail}")


def _parse_goes_hours(goes_hours: str) -> list[int]:
    """Parse the download hour syntax. Ranges are inclusive, e.g. 18-22 is 5 hours."""
    value = str(goes_hours).strip()
    if not value:
        return list(range(24))
    if "-" in value and "," not in value:
        start, end = (int(part) for part in value.split("-", 1))
        if end < start:
            raise ValueError(f"GOES_HOURS range must increase, got {goes_hours!r}")
        hours = list(range(start, end + 1))
    else:
        hours = [int(part) for part in value.split(",") if part.strip()]
    if any(hour < 0 or hour > 23 for hour in hours):
        raise ValueError(f"GOES_HOURS values must be between 0 and 23, got {goes_hours!r}")
    return hours


def _mask_window_from_goes_hours(goes_hours: str) -> tuple[float, float]:
    """Use the downloaded GOES hour range as the mask time window."""
    hours = _parse_goes_hours(goes_hours)
    if not hours:
        return 0.0, 24.0
    return float(min(hours)), float(min(max(hours) + 1, 24))


def _count_pending_ortho_inputs(month_dir: Path) -> int:
    if not month_dir.exists():
        return 0
    count = 0
    for path in month_dir.rglob("*.nc"):
        if path.name.endswith("_ortho.nc"):
            continue
        if path.with_name(path.stem + "_ortho.nc").exists():
            continue
        count += 1
    return count


def _dedupe_zarr_timestamps(config: WorkflowConfig, date: pd.Timestamp) -> None:
    date_ymd = config.date_ymd(date)
    for channel in ("C02", "C05", "C13"):
        zarr_path = config.base_dir / config.goes / channel / f"{config.goes}_{channel}_{config.domain}_{date_ymd}.zarr"
        if not zarr_path.is_dir():
            continue
        ds = xr.open_dataset(zarr_path)
        try:
            if "t" not in ds.coords or ds.indexes["t"].is_unique:
                continue
            keep = ~ds.indexes["t"].duplicated(keep="first")
            tmp = Path(str(zarr_path) + ".tmp")
            if tmp.exists():
                shutil.rmtree(tmp)
            ds.isel(t=keep).to_zarr(tmp, mode="w")
        finally:
            ds.close()
        shutil.rmtree(zarr_path)
        tmp.rename(zarr_path)


def _status(message: str) -> None:
    try:
        from IPython.display import clear_output

        clear_output(wait=True)
    except Exception:
        pass
    print(message, flush=True)
