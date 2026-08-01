"""Licence-free analytic evaluation of beam-down secondary designs.

One source of truth for the design math that both the GUI's **Design** tab and
the ``scripts/`` design tools use.  Nothing here imports ``quadoa`` or opens a
:class:`beamdown.session.Session`, and nothing here writes a file, so every
function is safe to call while a multi-hour sweep holds the licence seat.

What lives here
---------------
*Geometry*, moved verbatim from the scripts that grew it, so the scripts and the
GUI can never drift apart:

* :class:`Cassegrain`, :func:`close_design`, :func:`rays_to_f1` --- the
  hyperboloid closure, from ``scripts/design_cassegrain.py`` (which now imports
  them back; its printed output is unchanged).
* :func:`geometry_terms` --- the vectorised, sun-independent half of the axicon
  solve, from ``scripts/quick_axicon_iterate.py`` (same arrangement).

*Evaluators*, one per layout, each returning a plain ``dict`` of numbers with no
formatting opinions:

* :func:`eval_axicon` (tip height, half-angle)
* :func:`eval_cassegrain` (dish rim height, prime-focus height F1)
* :func:`eval_prime_focus` (focus height, no secondary at all)

How the energy index is built, and what it is not
-------------------------------------------------
``energy_index`` is a RELATIVE number: 1.0000 is the built axicon (tip 27,000
mm, half-angle 20 deg).  It is **not** MWh and must never be printed as such.
It is the product of two analytic factors, both evaluated on the sweeps' own
94-step time grid and DNI provider:

1. ``sum over (heliostat, timestep) of cos(AOI) * DNI`` --- the cosine loss of
   aiming each mirror at that layout's aim point, integrated the way the sweep
   integrates it;
2. the occlusion transfer curve below, read at the design's mean aim-ray
   axis-crossing height;

times the reflectivity ratio when the layout has a different number of
reflections than the axicon (prime focus has one bounce, not two).

The exact (slow, still licence-free) cross-checks are
``scripts/scan_axicon_annual.py`` and ``scripts/scan_cassegrain_annual.py``.
Against traced runs the analytic route has been vetted to within ~0.6%.

Spot sizes here are IDEAL-OPTICS solar images: the 0.5 deg sun projected
through perfect geometry, with the receiver-plane obliquity stretch applied and
nothing else.  Real traces run larger --- slope error, the astigmatic residual
the axicon needs correcting for, and finite mirror size are all absent.  Use
them to compare designs, not to predict a trace.
"""

from __future__ import annotations

import datetime as _dt
from dataclasses import dataclass
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent

# --------------------------------------------------------------------------
# Fixed plant geometry.  These are properties of the existing installation,
# not free parameters: the receiver aperture and the 15 m rim are what any
# candidate secondary has to work with.
# --------------------------------------------------------------------------
F2_MM = 7000.0                 # receiver on the axis (secondary_height - offset)
RIM_RADIUS_MM = 15000.0        # aperture radius, the axicon's, the 30 m cap
TIP_RAD_MM = 500.0             # the axicon tip is truncated: beams inside this miss
SUN_FULL_DEG = 0.5             # solar full angle used for ideal image sizes

# The built axicon, the reference every index is quoted against.
BUILT_TIP_MM = 27000.0
BUILT_ANGLE_DEG = 20.0
BUILT = (BUILT_TIP_MM, BUILT_ANGLE_DEG)

# The seven declinations the sweeps trace.  The annual cosine integral is taken
# over the same grid so an index computed here is comparable with a sweep.
TRACED_DATES = ["2026-12-21", "2026-01-21", "2026-02-20", "2026-03-20",
                "2026-04-21", "2026-05-21", "2026-06-21"]

# The rim height of the settled cassegrain design and the trade sweep around it.
DEFAULT_RIM_HEIGHT_MM = 32460.0
TRADE_RIM_HEIGHTS_MM = np.unique(
    np.concatenate([np.arange(24000.0, 34001.0, 1000.0), [DEFAULT_RIM_HEIGHT_MM]])
)

# --------------------------------------------------------------------------
# The occlusion transfer curve.
#
# PROVENANCE: these are EXACT union occlusion efficiencies (lit AND unblocked,
# the form the sweep's scalar branch uses) computed by
# ``scripts/scan_prime_focus_height.py`` over the full 645-mirror field at the
# representative instant 2026-02-20 09:27 -- the average-AOI instant -- with the
# cassegrain disc (15 m radius at z 32,460) as the shading body.  The x axis is
# the height at which the aim rays cross the tower axis: lower focus, flatter
# beams, more mirror-on-mirror blocking.
#
# It is used as a RELATIVE transfer curve over mean aim height.  What licenses
# that shortcut is that the axicon's own exact occlusion at the same instant --
# computed independently, with its cone rather than a disc -- lands on this
# curve at its own mean aim height to about 0.1%.  Monotone by construction;
# clamped at both ends rather than extrapolated.
# --------------------------------------------------------------------------
OCC_H = np.array([30000., 32000., 33000., 34000., 35000., 36000., 38000.,
                  40000., 43000., 47000.])
OCC_V = np.array([0.9226, 0.9307, 0.9345, 0.9382, 0.9412, 0.9445, 0.9494,
                  0.9534, 0.9582, 0.9626])

# F1 at which cassegrain blocking matches the built axicon's, from the same
# scan: the axicon reference row sits between the 35,000 and 36,000 rows, so a
# cassegrain at F1 = 36,000 is the "same blocking" comparison point.
BLOCKING_COMPARABLE_F1_MM = 36000.0


def occlusion_at(mean_aim_height_mm: float) -> float:
    """Union occlusion efficiency for a design whose aim rays cross at this height."""
    return float(np.interp(np.clip(float(mean_aim_height_mm), OCC_H[0], OCC_H[-1]),
                           OCC_H, OCC_V))


# ==========================================================================
# GEOMETRY -- the axicon
# ==========================================================================

def geometry_terms(R, tip, ang_deg):
    """Vectorised copy of the axicon solve's sun-independent geometry.

    Returns ``(x_r, y_r, x_a, y_a, s, s_prime)``: the aim-point offset in the
    radial/height plane, the point where the ray meets the cone, and the two
    sagittal conjugate distances whose reciprocal difference is the extra
    curvature the heliostat has to carry.
    """
    from beamdown.secondary.axicon import receiver_correction

    drop = tip - F2_MM
    x_r, y_r, x_a, y_a = receiver_correction(R, tip, drop, ang_deg)
    alpha = np.deg2rad(ang_deg)
    cone_dist = np.hypot(x_a, y_a)
    s_prime = -np.hypot(drop + y_a, x_a)
    radius_axicon = cone_dist / np.tan(alpha)
    axicon_aoi = np.arctan2(x_a, drop + y_a) + alpha
    s = 1.0 / (2.0 * np.cos(axicon_aoi) / radius_axicon - 1.0 / s_prime)
    return x_r, y_r, x_a, y_a, s, s_prime


# ==========================================================================
# GEOMETRY -- the cassegrain hyperboloid
# ==========================================================================

@dataclass(frozen=True)
class Cassegrain:
    """A closed hyperboloid secondary design. All lengths mm."""

    z1: float            # upper (prime) focus height
    z2: float            # lower focus = receiver height
    rim_r: float         # rim radius
    rim_z: float         # rim height
    a: float             # semi-transverse axis
    c: float             # half focus separation
    z_c: float           # centre (midpoint of foci)
    field_radius_mm: float
    aperture_filled: bool

    # -- derived conic quantities -------------------------------------------
    @property
    def b2(self) -> float:
        return self.c ** 2 - self.a ** 2

    @property
    def b(self) -> float:
        return float(np.sqrt(self.b2))

    @property
    def e(self) -> float:
        return self.c / self.a

    @property
    def K(self) -> float:
        return -self.e ** 2

    @property
    def R_v(self) -> float:
        """Vertex radius of curvature, positive for the sag convention above."""
        return self.b2 / self.a

    @property
    def vertex_z(self) -> float:
        """Vertex of the USED sheet: the one enclosing F1, so centre + a."""
        return self.z_c + self.a

    @property
    def sag_mm(self) -> float:
        """Depth of the dish, vertex to rim."""
        return self.rim_z - self.vertex_z

    # -- surface, exactly -----------------------------------------------------
    def sag_at(self, r):
        """Sag above the vertex, algebraic hyperbola form."""
        r = np.asarray(r, float)
        return self.a * (np.sqrt(1.0 + r ** 2 / self.b2) - 1.0)

    def sag_at_conic(self, r):
        """Same sag via the standard conic formula -- a cross-check on K/R_v."""
        r = np.asarray(r, float)
        R, K = self.R_v, self.K
        return r ** 2 / (R * (1.0 + np.sqrt(1.0 - (1.0 + K) * r ** 2 / R ** 2)))

    def surface_z(self, r):
        return self.vertex_z + self.sag_at(r)

    def implicit(self, p):
        """f(P) = (z-z_c)^2/a^2 - (x^2+y^2)/b^2 - 1 ; zero on the quadric."""
        p = np.atleast_2d(np.asarray(p, float))
        return ((p[:, 2] - self.z_c) ** 2 / self.a ** 2
                - (p[:, 0] ** 2 + p[:, 1] ** 2) / self.b2 - 1.0)

    def normal(self, p):
        """Unit normal, oriented to point DOWN (out of the reflective underside).

        Direction is irrelevant to the reflection formula ``d - 2(d.n)n`` but
        fixing it makes the printed incidence angles unambiguous.
        """
        p = np.atleast_2d(np.asarray(p, float))
        g = np.column_stack((
            -2.0 * p[:, 0] / self.b2,
            -2.0 * p[:, 1] / self.b2,
            2.0 * (p[:, 2] - self.z_c) / self.a ** 2,
        ))
        g /= np.linalg.norm(g, axis=1, keepdims=True)
        flip = g[:, 2] > 0.0
        g[flip] *= -1.0
        return g

    @property
    def F1(self) -> np.ndarray:
        return np.array([0.0, 0.0, self.z1])

    @property
    def F2(self) -> np.ndarray:
        return np.array([0.0, 0.0, self.z2])

    # -- ray work -------------------------------------------------------------
    def intersect(self, origins, dirs):
        """Exact ray/quadric intersection on the USED (upper, F1) sheet.

        Substituting P = O + t*d into ``implicit`` gives a quadratic in t. Both
        roots are returned so the caller can see the other sheet, plus the root
        selected as the physical hit: smallest t > 0 lying on the upper sheet
        (``z > z_c``).

        Returns ``(t_hit, hit_point, ok, t_roots)``.
        """
        o = np.atleast_2d(np.asarray(origins, float))
        d = np.atleast_2d(np.asarray(dirs, float))
        d = d / np.linalg.norm(d, axis=1, keepdims=True)

        ia2, ib2 = 1.0 / self.a ** 2, 1.0 / self.b2
        oz = o[:, 2] - self.z_c

        A = d[:, 2] ** 2 * ia2 - (d[:, 0] ** 2 + d[:, 1] ** 2) * ib2
        B = 2.0 * (oz * d[:, 2] * ia2 - (o[:, 0] * d[:, 0] + o[:, 1] * d[:, 1]) * ib2)
        C = oz ** 2 * ia2 - (o[:, 0] ** 2 + o[:, 1] ** 2) * ib2 - 1.0

        n = o.shape[0]
        roots = np.full((n, 2), np.nan)
        with np.errstate(invalid="ignore", divide="ignore"):
            disc = B ** 2 - 4.0 * A * C
            quad = (np.abs(A) > 1e-18) & (disc >= 0.0)
            sq = np.sqrt(np.where(disc >= 0.0, disc, 0.0))
            # numerically stable pair, then sorted
            q = -0.5 * (B + np.sign(np.where(B == 0.0, 1.0, B)) * sq)
            r1 = np.where(quad, q / A, np.nan)
            r2 = np.where(quad & (np.abs(q) > 1e-30), C / q, np.nan)
            lin = (np.abs(A) <= 1e-18) & (np.abs(B) > 1e-18)
            r1 = np.where(lin, -C / B, r1)
            roots[:, 0] = np.minimum(r1, r2)
            roots[:, 1] = np.maximum(r1, r2)
            roots[lin, 0] = r1[lin]

        t_hit = np.full(n, np.nan)
        for col in (0, 1):                    # smallest positive first
            t = roots[:, col]
            p = o + t[:, None] * d
            good = np.isfinite(t) & (t > 1e-9) & (p[:, 2] > self.z_c)
            take = good & ~np.isfinite(t_hit)
            t_hit[take] = t[take]

        hit = o + np.where(np.isfinite(t_hit), t_hit, 0.0)[:, None] * d
        return t_hit, hit, np.isfinite(t_hit), roots

    def reflect(self, dirs, points):
        d = np.atleast_2d(np.asarray(dirs, float))
        d = d / np.linalg.norm(d, axis=1, keepdims=True)
        n = self.normal(points)
        return d - 2.0 * np.sum(d * n, axis=1)[:, None] * n


def close_design(rim_z: float,
                 field_radius_mm: float,
                 z2: float,
                 rim_r: float = RIM_RADIUS_MM,
                 z1: float | None = None) -> Cassegrain:
    """Close the design. ``z1=None`` -> aperture fill; otherwise z1 is taken as given.

    Step 2 of the closure (see ``scripts/design_cassegrain.py``'s module
    docstring): z1 = z_r * R_f / (R_f - r_rim).
    Step 3: 2a = | |P-F1| - |P-F2| | at the rim point P.
    """
    filled = z1 is None
    if filled:
        if field_radius_mm <= rim_r:
            raise ValueError(
                f"field radius {field_radius_mm:.1f} mm must exceed the rim radius "
                f"{rim_r:.1f} mm for aperture fill to have a solution"
            )
        z1 = rim_z * field_radius_mm / (field_radius_mm - rim_r)
    z1 = float(z1)

    if z1 <= rim_z:
        raise ValueError(
            f"F1 (z={z1:.1f}) must be above the rim (z={rim_z:.1f}): the reflector "
            "has to intercept the bundle BELOW the prime focus"
        )

    P = np.array([rim_r, 0.0, rim_z])
    d1 = float(np.linalg.norm(P - np.array([0.0, 0.0, z1])))
    d2 = float(np.linalg.norm(P - np.array([0.0, 0.0, z2])))
    a = abs(d1 - d2) / 2.0
    c = (z1 - z2) / 2.0
    if not (0.0 < a < c):
        raise ValueError(f"degenerate conic: a={a:.3f}, c={c:.3f} (need 0 < a < c)")

    return Cassegrain(
        z1=z1, z2=float(z2), rim_r=float(rim_r), rim_z=float(rim_z),
        a=a, c=c, z_c=(z1 + z2) / 2.0,
        field_radius_mm=float(field_radius_mm), aperture_filled=filled,
    )


def rays_to_f1(x_mm, y_mm, des: Cassegrain):
    """Origins and unit directions of every heliostat's ray aimed at F1."""
    o = np.column_stack((np.asarray(x_mm, float),
                         np.asarray(y_mm, float),
                         np.zeros(np.size(x_mm))))
    d = des.F1[None, :] - o
    d /= np.linalg.norm(d, axis=1, keepdims=True)
    return o, d


# ==========================================================================
# The cached inputs: field, time grid, DNI.
# ==========================================================================

@dataclass(frozen=True)
class FieldData:
    """Everything the evaluators need about the plant, loaded once.

    Dragging a slider re-evaluates a design a few times a second, and loading
    the 645-row field file plus building the 94-step time grid takes far longer
    than the evaluation itself -- so this is built lazily and kept.
    """

    x: np.ndarray            # heliostat centres, mm
    y: np.ndarray
    R: np.ndarray            # field radius of each, mm
    sun: np.ndarray          # (steps, 3) unit sun vectors
    dni: np.ndarray          # (steps,) W/m^2
    reflectivity: float      # per-bounce mirror reflectivity, from config
    n_steps: int
    n_heliostats: int

    @property
    def R_min(self) -> float:
        return float(self.R.min())

    @property
    def R_max(self) -> float:
        return float(self.R.max())


_FIELD_CACHE: dict[str, FieldData] = {}
# id(field) -> (the FieldData itself, its evaluation).  The FieldData is kept in
# the value so the id can never be recycled onto a different object.
_BUILT_CACHE: dict[int, tuple[FieldData, dict]] = {}


def load_field_data(config_path=None) -> FieldData:
    """Load (and cache) the field, the sweeps' time grid and its DNI provider."""
    key = str(config_path) if config_path is not None else "<default>"
    hit = _FIELD_CACHE.get(key)
    if hit is not None:
        return hit

    from beamdown import field as F
    from beamdown.config import load_config
    from beamdown.dni import provider_for
    from beamdown.solar import build_time_grid

    cfg = load_config(config_path)
    full = F.load_field(cfg)
    x = np.asarray(full.x_mm, float)
    y = np.asarray(full.y_mm, float)

    dates = [_dt.date.fromisoformat(d) for d in TRACED_DATES]
    steps = build_time_grid(cfg, dates)
    prov = provider_for(cfg)
    azr = np.deg2rad([s.solar_az_deg for s in steps])
    elr = np.deg2rad([s.solar_el_deg for s in steps])
    sun = np.stack([np.cos(elr) * np.cos(np.pi / 2 - azr),
                    np.cos(elr) * np.sin(np.pi / 2 - azr), np.sin(elr)], axis=1)
    dni = np.array([prov.dni(s.date, s.hour) for s in steps])

    data = FieldData(
        x=x, y=y, R=np.hypot(x, y), sun=sun, dni=dni,
        reflectivity=float(cfg.optics.mirror_reflectivity),
        n_steps=len(steps), n_heliostats=len(x),
    )
    _FIELD_CACHE[key] = data
    return data


def _cos_annual(field: FieldData, aim: np.ndarray) -> tuple[float, np.ndarray]:
    """DNI-weighted sum of cos(AOI) over the grid, and the per-heliostat weights.

    ``aim`` is (n, 3): each heliostat's aim point.  The mirror normal bisects
    the sun vector and the vector to the aim point, so the incidence angle is
    half the angle between them, exactly as ``mirror.heliostat_orientation``
    computes it.
    """
    mirror = np.stack([field.x, field.y, np.zeros(field.n_heliostats)], axis=1)
    to_aim = aim - mirror
    f = np.linalg.norm(to_aim, axis=1)
    t = to_aim / f[:, None]
    cosaoi = np.cos(0.5 * np.arccos(np.clip(t @ field.sun.T, -1.0, 1.0)))
    weights = (cosaoi * field.dni[None, :]).sum(axis=1)
    return float(weights.sum()), weights


def built_axicon(field: FieldData | None = None) -> dict:
    """The built design's own evaluation, computed (never hardcoded), cached.

    Every ``*_vs_built`` ratio in this module is against this dict.  The
    sagittal-correction cap the owner's design rule uses is
    ``built_axicon()["max_sagittal_correction"]`` = 7.115e-06 /mm.
    """
    field = field or load_field_data()
    hit = _BUILT_CACHE.get(id(field))
    if hit is None:
        hit = (field, _eval_axicon_raw(BUILT_TIP_MM, BUILT_ANGLE_DEG, field))
        _BUILT_CACHE[id(field)] = hit
    return hit[1]


# ==========================================================================
# EVALUATORS
# ==========================================================================

def _eval_axicon_raw(tip_mm: float, angle_deg: float, field: FieldData) -> dict:
    """The axicon evaluation without the ``*_vs_built`` ratios.

    Split out so :func:`built_axicon` can compute the reference without
    recursing into itself.
    """
    theta_half = np.deg2rad(0.5 * SUN_FULL_DEG)
    n = field.n_heliostats
    R = field.R

    x_r, y_r, x_a, y_a, s, s_prime = geometry_terms(R, tip_mm, angle_deg)
    aim = np.stack([field.x / R * x_r, field.y / R * x_r,
                    np.full(n, tip_mm + y_r)], axis=1)
    mirror = np.stack([field.x, field.y, np.zeros(n)], axis=1)
    to_aim = aim - mirror
    f = np.linalg.norm(to_aim, axis=1)
    f_s = f + s + s_prime
    dpow = np.abs(1.0 / f_s - 1.0 / f)          # extra sagittal power, 1/mm

    # Coverage: every beam must land on the cone flank -- outside the truncated
    # tip, inside the rim -- and the sagittal conjugate must stay real.
    cov_ok = bool(((x_a > TIP_RAD_MM) & (x_a < RIM_RADIUS_MM)).all()
                  and (f_s > 0).all())

    cos_annual, w = _cos_annual(field, aim)

    # The aim rays are pushed off-axis by x_r (< 0, i.e. inward), so they cross
    # the tower axis above the aim point; that crossing height is what the
    # occlusion curve is indexed on.
    z_axis = aim[:, 2] * R / (R - x_r)
    mean_aim_height = float(z_axis.mean())
    occ = occlusion_at(mean_aim_height)

    # Ideal solar image at the receiver: the sun's angular radius through the
    # slant distance, stretched by the obliquity of arrival at the receiver
    # plane, and taken at the 90%-enclosed radius of a uniform disc.
    stretch = 1.0 / np.cos(np.arctan2(x_a, (tip_mm + y_a) - F2_MM))
    r90 = np.sqrt(0.9) * theta_half * f * np.sqrt(stretch)

    rim_z = tip_mm + RIM_RADIUS_MM * np.tan(np.deg2rad(angle_deg))
    return {
        "layout": "axicon",
        "tip_mm": float(tip_mm),
        "angle_deg": float(angle_deg),
        "feasible": cov_ok,
        "max_sagittal_correction": float(dpow.max()),
        "inner_hit_mm": float(x_a.min()),
        "outer_hit_mm": float(x_a.max()),
        "cos_annual": cos_annual,
        "occlusion": occ,
        "mean_aim_height_mm": mean_aim_height,
        "energy_raw": cos_annual * occ,
        "r90_mm": float(np.average(r90, weights=w)),
        "r90_inner_mm": float(r90[np.argmin(R)]),
        "r90_outer_mm": float(r90[np.argmax(R)]),
        "rim_z_mm": float(rim_z),
        "cone_depth_mm": float(rim_z - tip_mm),
    }


def eval_axicon(tip_mm: float, angle_deg: float,
                field: FieldData | None = None) -> dict:
    """Evaluate a cone secondary: tip height and half-angle.

    Admissibility follows the owner's design rule: the sagittal curvature
    correction demanded of the WORST heliostat must not exceed the built
    design's, and every beam must land on the flank.  ``correction_ratio`` is
    that comparison, 1.0 = exactly at the cap.
    """
    field = field or load_field_data()
    out = _eval_axicon_raw(tip_mm, angle_deg, field)
    ref = built_axicon(field)
    out["energy_index"] = out["energy_raw"] / ref["energy_raw"]
    out["correction_ratio"] = (out["max_sagittal_correction"]
                               / ref["max_sagittal_correction"])
    out["correction_cap"] = ref["max_sagittal_correction"]
    out["r90_vs_built"] = out["r90_mm"] / ref["r90_mm"] - 1.0
    out["admissible"] = bool(out["feasible"] and out["correction_ratio"] <= 1.0 + 1e-9)
    out["notes"] = []
    if not out["feasible"]:
        out["notes"].append(
            "some beams miss the cone flank: they fall inside the truncated tip "
            f"(r < {TIP_RAD_MM:,.0f} mm) or outside the {RIM_RADIUS_MM:,.0f} mm rim")
    elif out["correction_ratio"] > 1.0 + 1e-9:
        out["notes"].append(
            "the inner heliostats need more sagittal correction than the built "
            "design does -- outside the owner's admissibility rule")
    return out


def eval_cassegrain(rim_z_mm: float, f1_mm: float,
                    field: FieldData | None = None) -> dict:
    """Evaluate a hyperboloid secondary: dish rim height and prime focus F1.

    Feasibility is coverage: F1 has to sit above the rim (or there is nothing
    for the dish to intercept), the conic has to close, and every heliostat's
    beam must land inside the 15 m rim.
    """
    field = field or load_field_data()
    ref = built_axicon(field)
    theta_half = np.deg2rad(0.5 * SUN_FULL_DEG)

    base = {
        "layout": "cassegrain",
        "rim_z_mm": float(rim_z_mm),
        "f1_mm": float(f1_mm),
        "feasible": False,
        "notes": [],
    }

    try:
        des = close_design(float(rim_z_mm), field.R_max, F2_MM,
                           RIM_RADIUS_MM, float(f1_mm))
    except ValueError as exc:
        base["notes"].append(str(exc))
        return base

    o, d = rays_to_f1(field.x, field.y, des)
    t, hit, ok, _ = des.intersect(o, d)
    r_hit = np.hypot(hit[:, 0], hit[:, 1])
    tol = 1e-9
    inside = bool(ok.all() and (r_hit <= des.rim_r * (1.0 + tol)).all())

    base.update({
        "K": float(des.K),
        "R_v_mm": float(abs(des.R_v)),
        "vertex_z_mm": float(des.vertex_z),
        "sag_mm": float(des.sag_mm),
        "a_mm": float(des.a),
        "c_mm": float(des.c),
        "fill_fraction": float(r_hit.max() / des.rim_r) if ok.all() else float("nan"),
        "aperture_fill_f1_mm": float(rim_z_mm * field.R_max
                                     / (field.R_max - RIM_RADIUS_MM)),
    })

    if not inside:
        base["notes"].append(
            f"{int((~ok).sum() + (ok & (r_hit > des.rim_r * (1 + tol))).sum())} of "
            f"{field.n_heliostats} beams clear the {des.rim_r/1000:.0f} m rim -- "
            f"F1 is above the aperture-fill limit "
            f"{base['aperture_fill_f1_mm']:,.0f} mm for this rim height, so the "
            f"outer field spills past the dish")
        return base

    # The relay: the ideal image formed at F1 is re-imaged to F2 with
    # magnification |P->F2| / |P->F1| at the intercept P.
    slant = np.linalg.norm(des.F1[None, :] - o, axis=1)
    u = np.linalg.norm(des.F1[None, :] - hit, axis=1)
    v = np.linalg.norm(des.F2[None, :] - hit, axis=1)
    mag = v / u
    disk_d = np.deg2rad(SUN_FULL_DEG) * slant * mag

    # Obliquity of arrival at the horizontal receiver plane, from the dish.
    stretch = 1.0 / np.cos(np.arctan2(r_hit, hit[:, 2] - des.z2))
    r90 = np.sqrt(0.9) * theta_half * slant * mag * np.sqrt(stretch)

    aim = np.repeat(des.F1[None, :], field.n_heliostats, axis=0)
    cos_annual, w = _cos_annual(field, aim)

    # Every heliostat aims at the same on-axis point, so the aim rays cross the
    # axis exactly at F1 -- no averaging needed, unlike the axicon.
    occ = occlusion_at(des.z1)

    base.update({
        "feasible": True,
        "disk_outer_mm": float(disk_d[np.argmax(field.R)]),
        "disk_max_mm": float(disk_d.max()),
        "magnification_min": float(mag.min()),
        "magnification_max": float(mag.max()),
        "cos_annual": cos_annual,
        "occlusion": occ,
        "mean_aim_height_mm": float(des.z1),
        "energy_raw": cos_annual * occ,
        "energy_index": cos_annual * occ / ref["energy_raw"],
        "r90_mm": float(np.average(r90, weights=w)),
        "r90_inner_mm": float(r90[np.argmin(field.R)]),
        "r90_outer_mm": float(r90[np.argmax(field.R)]),
        "r90_vs_built": float(np.average(r90, weights=w)) / ref["r90_mm"] - 1.0,
        "blocking_note": _blocking_note(des.z1),
    })
    return base


def eval_prime_focus(f1_mm: float, field: FieldData | None = None) -> dict:
    """Evaluate the no-secondary layout: the receiver sits at the prime focus.

    There is no relay and no second reflection, so the energy index carries one
    bounce of mirror reflectivity where the axicon and cassegrain carry two --
    that alone is worth about +11% against the built axicon and is the reason
    this layout wins on energy while losing on where the receiver has to hang.
    """
    field = field or load_field_data()
    ref = built_axicon(field)
    theta_half = np.deg2rad(0.5 * SUN_FULL_DEG)

    aim = np.array([[0.0, 0.0, float(f1_mm)]]).repeat(field.n_heliostats, axis=0)
    slant = np.hypot(field.R, float(f1_mm))

    # Arrival obliquity at a horizontal receiver plane at F1: a beam from radius
    # R arrives atan(R / z1) off vertical, and the image stretches by 1/cos.
    tilt = np.arctan2(field.R, float(f1_mm))
    stretch = 1.0 / np.cos(tilt)
    disk_d = np.deg2rad(SUN_FULL_DEG) * slant
    r90 = np.sqrt(0.9) * theta_half * slant * np.sqrt(stretch)

    cos_annual, w = _cos_annual(field, aim)
    occ = occlusion_at(float(f1_mm))

    # One reflection instead of two.
    bounce_gain = 1.0 / field.reflectivity
    energy_raw = cos_annual * occ * bounce_gain

    return {
        "layout": "prime_focus",
        "f1_mm": float(f1_mm),
        "feasible": True,
        "notes": [],
        "disk_outer_mm": float(disk_d[np.argmax(field.R)]),
        "disk_max_mm": float(disk_d.max()),
        "disk_inner_mm": float(disk_d[np.argmin(field.R)]),
        "stretch_inner": float(stretch[np.argmin(field.R)]),
        "stretch_outer": float(stretch[np.argmax(field.R)]),
        "tilt_inner_deg": float(np.rad2deg(tilt[np.argmin(field.R)])),
        "tilt_outer_deg": float(np.rad2deg(tilt[np.argmax(field.R)])),
        "cos_annual": cos_annual,
        "occlusion": occ,
        "mean_aim_height_mm": float(f1_mm),
        "reflectivity": field.reflectivity,
        "bounce_gain": bounce_gain,
        "energy_raw": energy_raw,
        "energy_index": energy_raw / ref["energy_raw"],
        "r90_mm": float(np.average(r90, weights=w)),
        "r90_inner_mm": float(r90[np.argmin(field.R)]),
        "r90_outer_mm": float(r90[np.argmax(field.R)]),
        "r90_vs_built": float(np.average(r90, weights=w)) / ref["r90_mm"] - 1.0,
        "blocking_note": _blocking_note(float(f1_mm)),
    }


def _blocking_note(f1_mm: float) -> str:
    """Plain-language comparison of this focus height's blocking with the axicon's."""
    delta = f1_mm - BLOCKING_COMPARABLE_F1_MM
    if abs(delta) < 250.0:
        return ("about the same field blocking as the built axicon "
                "(F1 36 m is the measured match)")
    if delta < 0:
        return (f"{abs(delta)/1000:.1f} m below the F1 that matches the axicon's "
                f"blocking -- flatter beams, so more mirror-on-mirror blocking")
    return (f"{delta/1000:.1f} m above the F1 that matches the axicon's blocking "
            f"-- steeper beams, so less mirror-on-mirror blocking")


# ==========================================================================
# One-line honesty footer, so every consumer says the same thing.
# ==========================================================================
HONESTY_FOOTER = (
    "Analytic estimates (validated ±0.6% vs traced runs); spots are "
    "ideal-optics — real traces run larger."
)
