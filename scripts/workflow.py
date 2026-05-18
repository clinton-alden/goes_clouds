"""Small public API for the GOES RGB + ERA5 cloud-mask workflow."""

from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path
import shutil
import subprocess
import sys

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
    goes_hours: str = "14-23"
    start_hour_utc: float = 14
    end_hour_utc: float = 24
    threshold_csv: Path | None = None
    overwrite_mask: bool = True
    run: bool = True

    def __post_init__(self) -> None:
        self.base_dir = Path(self.base_dir)
        if self.threshold_csv is not None:
            self.threshold_csv = Path(self.threshold_csv)

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
        return self.base_dir / self.goes / "cloud_mask_tempbin_10c"

    @property
    def gif_dir(self) -> Path:
        return self.base_dir / self.goes / "gif_loops_tempbin_10c"

    @property
    def resolved_threshold_csv(self) -> Path:
        if self.threshold_csv is not None:
            return self.threshold_csv
        return self.repo_dir / "thresholds" / "gothic_temp_bin_rgb_thresholds_10c.csv"

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
        return self.mask_dir / f"{rgb_path.stem}_cloud_binary_tempbin10c.nc"

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
    if not (Path.home() / ".cdsapirc").exists():
        raise RuntimeError(
            "Missing ~/.cdsapirc. Copy example_cdsapirc to ~/.cdsapirc and add your Copernicus/CDS key."
        )


def download_goes(config: WorkflowConfig) -> list[Path]:
    """Download GOES C02/C05/C13 files for the configured date range."""
    _ensure_dirs(config)
    download_script = config.scripts_dir / "download-goes.py"
    outputs: list[Path] = []
    tasks = [(date, channel) for date in config.dates for channel in ("C02", "C05", "C13")]
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
        _status(f"Orthorectifying {date:%Y-%m-%d}")
        _run([sys.executable, ortho_script, month_dir, config.domain], config)
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
                config.base_dir,
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
            )
        outputs.append(config.rgb_path(date))
    _status("RGB step complete")
    return outputs


def apply_mask(config: WorkflowConfig) -> list[Path]:
    """Download/reuse ERA5-Land temperature and apply the RGB cloud mask."""
    mask_script = config.scripts_dir / "apply_tempbin_thresholds.py"
    outputs = []
    for date in tqdm(config.dates, desc="mask", unit="day"):
        rgb_path = config.rgb_path(date)
        _status(f"Applying cloud mask for {date:%Y-%m-%d}")
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
            "--start-hour-utc",
            config.start_hour_utc,
            "--end-hour-utc",
            config.end_hour_utc,
        ]
        if config.overwrite_mask:
            cmd.append("--overwrite")
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
