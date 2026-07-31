"""Direct normal irradiance by date and hour.

DNI never affects the ray trace. The Quadoa source is fixed at 38484.5 W over a
38.4845 m^2 aperture -- exactly 1000 W/m^2 -- so real DNI enters as a scale
factor at analysis time:

    flux = count * W_per_ray * throughput * eta_shade * eta_block
           * (dni / 1000) / bin_area

That means the DNI model can be swapped, refined, or replaced with measured data
at any point **without re-tracing anything**.

Providers
---------
``ConstantDNI``        fixed value; what reproduces prior work.
``TableDNI``           measured or downloaded hourly series, matched by exact
                        calendar day. Carries real day-to-day weather, which
                        means it also carries whatever weather the TMY splice
                        happened to donate for that one day (see
                        ``MonthlyProfileDNI``).
``MonthlyProfileDNI``   mean diurnal curve per calendar month, averaged over
                        every day in that month -- the default. See
                        config.toml [dni] for the reasoning.

Two free online sources are supported by :func:`fetch`, neither needing an API
key:

``pvgis``  EU JRC PVGIS typical meteorological year. Purpose-built for annual
           energy estimates, global coverage, returns Gb(n) = DNI.
``nasa``   NASA POWER hourly ALLSKY_SFC_SW_DNI for a specific historical year.

Both return UTC timestamps, which are converted to the site's local clock time
using the configured timezone offset.
"""

from __future__ import annotations

import datetime as _dt
from abc import ABC, abstractmethod
from pathlib import Path

import numpy as np
import pandas as pd

STANDARD_DNI = 1000.0  # W/m^2, the DNI the ray trace is normalised to


class DNIProvider(ABC):
    """Returns DNI in W/m^2 for a given date and local clock hour."""

    @abstractmethod
    def dni(self, date: _dt.date, hour: float) -> float: ...

    def scale(self, date: _dt.date, hour: float) -> float:
        """Multiplier to apply to trace-derived flux."""
        return self.dni(date, hour) / STANDARD_DNI

    def describe(self) -> str:
        return self.__class__.__name__


class ConstantDNI(DNIProvider):
    def __init__(self, value: float = STANDARD_DNI):
        self.value = float(value)

    def dni(self, date: _dt.date, hour: float) -> float:
        return self.value

    def describe(self) -> str:
        return f"ConstantDNI({self.value:g} W/m2)"


class TableDNI(DNIProvider):
    """Hourly DNI from a table, matched by day-of-year and hour.

    The table is matched on (month, day, hour) rather than absolute date, so a
    TMY or any single historical year can drive a sweep configured for any year.
    Hours are linearly interpolated; missing days fall back to ``default``.
    """

    def __init__(self, frame: pd.DataFrame, default: float = STANDARD_DNI, source: str = ""):
        required = {"month", "day", "hour", "dni_w_m2"}
        missing = required - set(frame.columns)
        if missing:
            raise ValueError(f"DNI table missing columns: {sorted(missing)}")
        self.default = float(default)
        self.source = source
        self._by_day: dict[tuple[int, int], tuple[np.ndarray, np.ndarray]] = {}
        for (m, d), grp in frame.groupby(["month", "day"], sort=False):
            grp = grp.sort_values("hour")
            self._by_day[(int(m), int(d))] = (
                grp["hour"].to_numpy(float),
                grp["dni_w_m2"].to_numpy(float),
            )

    def dni(self, date: _dt.date, hour: float) -> float:
        entry = self._by_day.get((date.month, date.day))
        if entry is None:
            return self.default
        hours, values = entry
        if hours.size == 1:
            return float(values[0])
        return float(np.interp(hour, hours, values, left=values[0], right=values[-1]))

    def describe(self) -> str:
        return f"TableDNI({len(self._by_day)} days from {self.source or 'table'})"


class MonthlyProfileDNI(DNIProvider):
    """Mean diurnal DNI curve per calendar month, averaged across that month's days.

    Built by grouping the source table on ``(month, hour)`` and averaging over
    every row sharing that key -- deliberately *not* built by first averaging
    each day down to a scalar and then averaging days together, and not keyed
    on individual calendar days at query time at all. Hours are interpolated
    within the resulting profile exactly as :class:`TableDNI` interpolates
    within a day.

    The ``(month, hour)`` grouping also sidesteps the PVGIS TMY month-splice
    stub for free: (month=2, day=29) is a leftover 3-row bucket (local
    21:00-23:59) produced when the Feb/Mar splice boundary walks past the end
    of February in a leap donor year (see :class:`TableDNI` and
    ``scripts/check_dni.py``), and the matching (2, 28) bucket is short by
    those same 3 hours. Grouping by day first and then averaging days would
    have to special-case this (a near-empty "day" skews a per-day mean, and
    naively averaging 29 day-means for a 28-day month is wrong). Grouping by
    (month, hour) instead needs no such logic: at hour 21:00-23:59 the 3-row
    (2, 29) stub simply supplies the count that (2, 28) is missing, so every
    hour of February is averaged over exactly 28 rows -- the true day count --
    regardless of which day-bucket those rows happen to live in.
    """

    def __init__(self, frame: pd.DataFrame, default: float = STANDARD_DNI, source: str = ""):
        required = {"month", "day", "hour", "dni_w_m2"}
        missing = required - set(frame.columns)
        if missing:
            raise ValueError(f"DNI table missing columns: {sorted(missing)}")
        self.default = float(default)
        self.source = source
        self._by_month: dict[int, tuple[np.ndarray, np.ndarray]] = {}
        for m, grp in frame.groupby("month", sort=False):
            profile = grp.groupby("hour")["dni_w_m2"].mean().sort_index()
            self._by_month[int(m)] = (
                profile.index.to_numpy(float),
                profile.to_numpy(float),
            )

    def dni(self, date: _dt.date, hour: float) -> float:
        entry = self._by_month.get(date.month)
        if entry is None:
            return self.default
        hours, values = entry
        if hours.size == 1:
            return float(values[0])
        return float(np.interp(hour, hours, values, left=values[0], right=values[-1]))

    def day_kwh_m2(self, month: int) -> float:
        """Day-integrated energy under this month's average diurnal curve."""
        hours, values = self._by_month[month]
        return float(np.trapz(values, hours) / 1000.0)

    def describe(self) -> str:
        return f"MonthlyProfileDNI({len(self._by_month)} months from {self.source or 'table'})"


def _tidy(frame: pd.DataFrame, tz_offset_hours: float) -> pd.DataFrame:
    """Convert a UTC timestamp + DNI frame into month/day/hour local form."""
    local = frame["timestamp_utc"] + pd.to_timedelta(tz_offset_hours, unit="h")
    return pd.DataFrame(
        {
            "month": local.dt.month,
            "day": local.dt.day,
            "hour": local.dt.hour + local.dt.minute / 60.0,
            "dni_w_m2": frame["dni_w_m2"].to_numpy(float),
        }
    )


_TABLE_MODES = {"table": TableDNI, "monthly": MonthlyProfileDNI}


def provider_for(cfg, mode: str | None = None) -> DNIProvider:
    """Build a DNI provider by name.

    ``mode`` overrides ``cfg.dni.mode`` for this call only; pass ``None`` (the
    default) to build whatever config.toml's ``[dni] mode`` currently says.
    This is the one-line switch a live GUI control needs -- it takes the
    already-loaded ``cfg`` and a mode string, and never touches config.toml or
    re-reads it, so flipping modes at runtime doesn't require a config reload.
    """
    spec = getattr(cfg, "dni", None)
    mode = mode or (spec.mode if spec is not None else "constant")

    if mode == "constant":
        value = STANDARD_DNI if spec is None else spec.constant_w_m2
        return ConstantDNI(value)

    if mode in _TABLE_MODES:
        if spec is None:
            raise ValueError(f"DNI mode {mode!r} requires a [dni] section in config.toml")
        path = cfg.path(spec.table_file)
        if not path.exists():
            raise FileNotFoundError(
                f"DNI table {path} not found. Run `python -m beamdown fetch-dni` "
                f"to download it, or set [dni] mode = \"constant\"."
            )
        frame = pd.read_excel(path) if path.suffix.lower() in (".xlsx", ".xls") else pd.read_csv(path)
        return _TABLE_MODES[mode](frame, default=spec.constant_w_m2, source=path.name)

    raise ValueError(f"unknown DNI mode {mode!r}; use one of "
                     f"{sorted({'constant'} | set(_TABLE_MODES))}")


def load_dni_provider(cfg) -> DNIProvider:
    """Build the provider described by the ``[dni]`` section of config.toml.

    Thin wrapper over :func:`provider_for` kept for existing call sites.
    """
    return provider_for(cfg)


# --------------------------------------------------------------------------
# Online fetchers. Network access only happens when these are called directly.
# --------------------------------------------------------------------------

def fetch(source: str, cfg, out_path: Path | None = None, year: int | None = None) -> Path:
    """Download an hourly DNI series and cache it as a tidy CSV.

    Neither source requires an API key. Returns the path written.
    """
    site = cfg.site
    out_path = Path(out_path) if out_path else cfg.path(f"dni_{source}.csv")

    if source == "pvgis":
        frame = _fetch_pvgis_tmy(site.latitude, site.longitude)
    elif source == "nasa":
        frame = _fetch_nasa_power(site.latitude, site.longitude, year)
    else:
        raise ValueError(f"unknown DNI source {source!r}; use 'pvgis' or 'nasa'")

    tidy = _tidy(frame, site.timezone)
    tidy.to_csv(out_path, index=False)
    return out_path


def _fetch_pvgis_tmy(lat: float, lon: float) -> pd.DataFrame:
    """PVGIS typical meteorological year. Column ``Gb(n)`` is DNI in W/m^2."""
    import io
    import urllib.request

    url = (
        "https://re.jrc.ec.europa.eu/api/v5_2/tmy"
        f"?lat={lat}&lon={lon}&outputformat=csv"
    )
    with urllib.request.urlopen(url, timeout=120) as resp:
        text = resp.read().decode("utf-8", errors="replace")

    lines = text.splitlines()
    start = next(i for i, ln in enumerate(lines) if ln.startswith("time(UTC)"))
    end = start + 1
    while end < len(lines) and lines[end] and lines[end][0].isdigit():
        end += 1
    table = pd.read_csv(io.StringIO("\n".join(lines[start:end])))

    return pd.DataFrame(
        {
            "timestamp_utc": pd.to_datetime(table["time(UTC)"], format="%Y%m%d:%H%M"),
            "dni_w_m2": table["Gb(n)"].to_numpy(float),
        }
    )


def _fetch_nasa_power(lat: float, lon: float, year: int | None) -> pd.DataFrame:
    """NASA POWER hourly all-sky DNI for one historical year."""
    import json
    import urllib.request

    if year is None:
        year = _dt.date.today().year - 1
    url = (
        "https://power.larc.nasa.gov/api/temporal/hourly/point"
        "?parameters=ALLSKY_SFC_SW_DNI&community=RE"
        f"&latitude={lat}&longitude={lon}"
        f"&start={year}0101&end={year}1231&format=JSON"
    )
    with urllib.request.urlopen(url, timeout=180) as resp:
        payload = json.loads(resp.read().decode("utf-8"))

    series = payload["properties"]["parameter"]["ALLSKY_SFC_SW_DNI"]
    stamps = pd.to_datetime(list(series.keys()), format="%Y%m%d%H")
    values = np.array(list(series.values()), dtype=float)
    values[values < -900] = 0.0  # POWER fill value
    return pd.DataFrame({"timestamp_utc": stamps, "dni_w_m2": values})
