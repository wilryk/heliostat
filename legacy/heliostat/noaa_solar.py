import numpy as np
from datetime import datetime

def solar_position_calculator(Lat=None, Lon=None, TZone=None,
                                 Year=None, Month=None, Day=None, Time=None):
    """
    Python translation of SolarPositionCalculatorV3 (MATLAB).
    Returns:
      - If Time is None (6 args): sunrise and sunset (fraction of day)
      - If Time is provided (7 args): solar azimuth and elevation (deg)
      - If no args: default Tucson example
    """

    # -----------------------------
    # Handle inputs (MATLAB nargin)
    # -----------------------------
    if Lat is None:
        Lat = 32.0944
        Lon = -110.8147
        TZone = -7
        # Year = 2020
        # Month = 6
        # Day = 21
        # Time = 12 / 24  # Noon, fraction of day
        nargs = 0
    elif Time is None:
        Time = 0.5
        nargs = 6
    else:
        nargs = 7

    # -----------------------------
    # Helper functions
    # -----------------------------
    radians = np.deg2rad
    degrees = np.rad2deg

    # -----------------------------
    # Date calculations
    # -----------------------------
    dt = datetime(Year, Month, Day)
    matlab_datenum = dt.toordinal() + 366  # MATLAB serial date
    ex_date = matlab_datenum - 693960      # Excel serial date

    jul_day = ex_date + 2415018.5 + Time - TZone / 24
    jul_cen = (jul_day - 2451545) / 36525

    # -----------------------------
    # Solar geometry
    # -----------------------------
    GMLS = np.mod(
        280.46646 + jul_cen * (36000.76983 + jul_cen * 0.0003032), 360
    )

    GMAS = 357.52911 + jul_cen * (35999.05029 - 0.0001537 * jul_cen)

    EEO = 0.016708634 - jul_cen * (0.000042037 + 0.0000001267 * jul_cen)

    SEOC = (
        np.sin(radians(GMAS)) *
        (1.914602 - jul_cen * (0.004817 + 0.000014 * jul_cen))
        + np.sin(radians(2 * GMAS)) * (0.019993 - 0.000101 * jul_cen)
        + np.sin(radians(3 * GMAS)) * 0.000289
    )

    STL = GMLS + SEOC
    STA = GMAS + SEOC

    SRV = (1.000001018 * (1 - EEO**2)) / (1 + EEO * np.cos(radians(STA)))

    SAL = STL - 0.00569 - 0.00478 * np.sin(radians(125.04 - 1934.136 * jul_cen))

    MOE = (
        23 +
        (26 + (21.448 - jul_cen * (46.815 + jul_cen *
         (0.00059 - jul_cen * 0.001813))) / 60) / 60
    )

    OC = MOE + 0.00256 * np.cos(radians(125.04 - 1934.136 * jul_cen))

    SRA = degrees(np.arctan2(
        np.cos(radians(OC)) * np.sin(radians(SAL)),
        np.cos(radians(SAL))
    ))

    SDE = degrees(np.arcsin(
        np.sin(radians(OC)) * np.sin(radians(SAL))
    ))

    vy = np.tan(radians(OC / 2))**2

    EOT = 4 * degrees(
        vy * np.sin(2 * radians(GMLS))
        - 2 * EEO * np.sin(radians(GMAS))
        + 4 * EEO * vy * np.sin(radians(GMAS)) * np.cos(2 * radians(GMLS))
        - 0.5 * vy**2 * np.sin(4 * radians(GMLS))
        - 1.25 * EEO**2 * np.sin(2 * radians(GMAS))
    )

    HASR = degrees(np.arccos(
        np.cos(radians(90.833)) /
        (np.cos(radians(Lat)) * np.cos(radians(SDE)))
        - np.tan(radians(Lat)) * np.tan(radians(SDE))
    ))

    SN = (720 - 4 * Lon - EOT + TZone * 60) / 1440
    SRT = (SN * 1440 - HASR * 4) / 1440
    SST = (SN * 1440 + HASR * 4) / 1440

    TST = np.mod(Time * 1440 + EOT + 4 * Lon - 60 * TZone, 1440)

    HA = np.where(TST / 4 < 0, TST / 4 + 180, TST / 4 - 180)

    SZA = degrees(np.arccos(
        np.sin(radians(Lat)) * np.sin(radians(SDE))
        + np.cos(radians(Lat)) * np.cos(radians(SDE)) * np.cos(radians(HA))
    ))

    solelu = 90 - SZA

    # -----------------------------
    # Atmospheric refraction
    # -----------------------------
    AAR = np.zeros_like(solelu)

    C1 = solelu >= 85
    C2 = (solelu >= 5) & (solelu < 85)
    C3 = (solelu < 5) & (solelu >= -0.575)
    C4 = solelu < -0.575

    AAR[C1] = 0

    AAR[C2] = (
        (58.1 / np.tan(radians(solelu[C2]))
         - 0.07 / np.tan(radians(solelu[C2]))**3
         + 0.000086 / np.tan(radians(solelu[C2]))**5) / 3600
    )

    AAR[C3] = (
        (1735 + solelu[C3] *
         (-518.2 + solelu[C3] *
          (103.4 + solelu[C3] *
           (-12.79 + solelu[C3] * 0.711)))) / 3600
    )

    AAR[C4] = -20.772 / np.tan(radians(solelu[C4])) / 3600

    # -----------------------------
    # Azimuth & corrected elevation
    # -----------------------------
    theta = degrees(np.arccos(
        ((np.sin(radians(Lat)) * np.cos(radians(SZA))) - np.sin(radians(SDE)))
        / (np.cos(radians(Lat)) * np.sin(radians(SZA)))
    ))

    solaz = np.mod(180 + (1 - np.sign(HA)) * 180 + np.sign(HA) * theta, 360)
    solel = solelu + AAR

    # -----------------------------
    # Outputs
    # -----------------------------
    if nargs == 6:
        return SRT, SST
    else:
        return solaz, solel
