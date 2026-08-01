"""Build the 25-configuration figure model for a given date and time.

What the 25 configurations ARE
------------------------------
Twenty-five HELIOSTATS, all at ONE instant -- not 25 sun positions, not a 5x5
grid of anything.  That is settled by the model file itself rather than by
inference: in ``heliostat_field_model_mcfg.optx`` the sun is carried by

    <single_param id="solaz">   <single_param id="solze">

and a ``single_param`` holds ONE value for the whole model (Quadoa returns NaN
for ``getMulticonfParam(name, i>0)`` on them).  A configuration therefore
*cannot* carry its own sun.  What varies per configuration is exactly the seven
``<param>`` entries, which are exactly one heliostat's placement, pointing and
shape::

    posx  posy       position in the field, mm
    rot_az  rot_el   pointing, deg
    c3  c4  c5       Zernike terms giving the mirror its optical power

Which 25 heliostats is settled too, and not by this script:
``config.toml`` ``[field] n_configs = 25`` ("heliostats in the figure/3D-view
model") plus ``downselect_file = "data/downselected_x,y centers.xlsx"``, which
holds exactly 25 (x, y) pairs.  :func:`beamdown.field.load_or_build_downselect`
matches them back to field rows, and falls back to the farthest-point
downselect if the file is ever missing.  This script calls that function; it
does not have its own opinion.

Why this script exists
----------------------
``beamdown.model_edit.build_figure_model`` was the generator, and it was NEVER
called from anywhere in the repository -- no CLI verb, no GUI action, no test.
It also could not have produced a usable file without a licence seat: it did
the column surgery in Python but then opened a ``QuadoaSession`` to write the
per-config values through ``setMulticonfParam`` and ``saveModelFile``.  The
shipped ``models/figure_model_25cfg.optx`` shows that second half never
happened.  It is byte-for-byte ``heliostat_field_model_mcfg.optx`` plus seven
``val_24`` lines and ``columns="24"`` -> ``"25"``, and nothing else: config 0
holds one stale heliostat at (22300, -60000) and configs 1..24 are all zeros.
A zero heliostat is not "off" -- it is a mirror at the origin pointing at the
horizon, so the file opens fine in the GUI and shows 24 wrong mirrors.  That
superseded function has since been deleted from ``beamdown/model_edit.py``;
only its column-surgery helpers survive there, for the inspect/export path.

So this script does the whole job by text surgery and needs NO licence.  The
values a ``setMulticonfParam`` + ``saveModelFile`` round trip would have left in
the file are simply written into the file directly; they live in plain
``<variable ... value="...">`` attributes and nothing else in the document
refers to them.

Sun position and pointing are NOT reimplemented here: ``beamdown.solar``
provides the former and the configured ``beamdown.secondary`` strategy the
latter, so a figure model and a sweep step at the same instant agree by
construction.

Same guarded-text-surgery contract as ``scripts/build_prime_focus_model.py``:
every edit is anchored, every anchor is asserted before it is used, the base is
never modified, and the output is verified by re-reading it.

Usage::

    python scripts/build_figure_model.py --date 2026-06-21 --hour 12.0
    python scripts/build_figure_model.py --date 2026-06-21 --hour 12.0 --flat
    python scripts/build_figure_model.py --date 2026-06-21 --hour 12.0 --check

``--check`` re-reads the file from disk, recomputes the sun and the pointing
solve from the same ``--date``/``--hour``, and asserts every one of the
25 x 7 + 5 values landed.  It is a pure-Python check and needs no seat; the
licence-gated half (does Quadoa agree there are 25 configurations, do the
writes reach the trace) is ``scripts/verify_figure_model.py``.

The output NEVER overwrites ``models/figure_model_25cfg.optx``.  That file is
the historical artefact this script was written to explain; leave it alone.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

import numpy as np  # noqa: E402

# Reused verbatim rather than reimplemented: Quadoa writes one tag with a
# duplicate attribute, so a strict XML parser rejects the base model as-is.
from build_prime_focus_model import (  # noqa: E402
    StructuralSurprise,
    newline_of,
    xml_check,
)

BASE = REPO / "models" / "heliostat_field_model_mcfg.optx"

# The file the archaeology was done on. Never written to.
PRISTINE = REPO / "models" / "figure_model_25cfg.optx"

# The seven per-configuration parameters, in the model's own order, and the
# five model-wide ones. Asserted exactly -- a base model with a different set
# is a different model and this script must not guess at it.
PARAMS = ["rot_az", "rot_el", "posx", "posy", "c3", "c4", "c5"]
SINGLE_PARAMS = ["solaz", "solze", "sec_height", "rec_offset", "axi_angle"]

# Quadoa's own indentation inside the multiconfig block.
VAR_INDENT = "\t\t\t"
CLOSE_INDENT = "\t\t"

_MULTICONF = re.compile(r'<multiconfig\s+columns="(\d+)"([^>]*)>')
_PARAM_BLOCK = re.compile(
    r'(?P<open><param\s+id="(?P<id>[^"]+)"[^>]*>)(?P<body>.*?)(?P<close></param>)',
    re.DOTALL,
)
_SINGLE_BLOCK = re.compile(
    r'(?P<open><single_param\s+id="(?P<id>[^"]+)"[^>]*>)(?P<body>.*?)'
    r"(?P<close></single_param>)",
    re.DOTALL,
)
_VAR_LINE = re.compile(r'<variable\s+name="val_(\d+)"[^>]*?/>')
_SINGLE_VAR = re.compile(r'<variable\s+name="val"[^>]*?/>')
_VALUE_ATTR = re.compile(r'value="[^"]*"')


def _fmt(value: float) -> str:
    """Serialise like Quadoa does: repr-ish, full precision, no exponent games.

    The base model carries ``value="128.812744044287"`` -- 15 significant
    digits, plain decimal for values of ordinary magnitude and E-notation for
    the Zernikes (``-7.959e-08``).  ``repr`` of a Python float reproduces both
    conventions and always round-trips exactly, which is what ``--check``
    depends on.
    """
    return repr(float(value))


# --------------------------------------------------------------------------
# Reading the base: assert everything before touching anything.
# --------------------------------------------------------------------------
def _multiconfig_span(text: str) -> tuple[int, int]:
    m = _MULTICONF.search(text)
    if m is None:
        raise StructuralSurprise("no <multiconfig columns=...> element")
    close = text.find("</multiconfig>", m.end())
    if close < 0:
        raise StructuralSurprise("<multiconfig> is never closed")
    return m.start(), close + len("</multiconfig>")


def check_base(text: str) -> list[str]:
    """Everything this script assumes about the base, checked out loud."""
    notes: list[str] = []

    ok, fixed, msg = xml_check(text)
    if not ok:
        raise StructuralSurprise(f"base model does not parse: {msg}")
    notes.append(msg)

    start, end = _multiconfig_span(text)
    block = text[start:end]
    m = _MULTICONF.search(text)
    cols = int(m.group(1))
    notes.append(
        f"multiconfig: columns={cols}, "
        f"lines {text[:start].count(chr(10)) + 1}-{text[:end].count(chr(10)) + 1}"
    )

    have_p = re.findall(r'<param id="([^"]*)"', text)
    have_s = re.findall(r'<single_param id="([^"]*)"', text)
    if have_p != PARAMS or have_s != SINGLE_PARAMS:
        raise StructuralSurprise(
            f"multiconfig is params={have_p} single_params={have_s}, expected "
            f"{PARAMS} / {SINGLE_PARAMS} -- this is not the base model this "
            f"script was written against"
        )
    notes.append(f"  {len(have_p)} per-config <param>: {', '.join(have_p)}")
    notes.append(f"  {len(have_s)} model-wide <single_param>: {', '.join(have_s)}")

    # Every <param> must carry a contiguous val_0..val_{cols-1}. A gap would
    # make "column i" mean different things in different parameters, and the
    # rebuild below would silently paper over it.
    for pm in _PARAM_BLOCK.finditer(block):
        idx = [int(i) for i in re.findall(r'<variable\s+name="val_(\d+)"', pm.group("body"))]
        if idx != list(range(cols)):
            raise StructuralSurprise(
                f"<param id={pm.group('id')!r}> has val indices {idx[:5]}... "
                f"({len(idx)} of them), expected a contiguous 0..{cols - 1}"
            )
    notes.append(f"  every <param> carries a contiguous val_0..val_{cols - 1}")

    for sm in _SINGLE_BLOCK.finditer(block):
        n = len(_SINGLE_VAR.findall(sm.group("body")))
        if n != 1:
            raise StructuralSurprise(
                f"<single_param id={sm.group('id')!r}> has {n} "
                f'<variable name="val">, expected exactly 1'
            )
    notes.append('  every <single_param> carries exactly one <variable name="val">')

    # Nothing outside the multiconfig block is touched, so record its size to
    # prove that afterwards.
    notes.append(
        f"  outside the block: {start:,} bytes before, "
        f"{len(text) - end:,} after (must be identical in the output)"
    )
    return notes


# --------------------------------------------------------------------------
# The build.
# --------------------------------------------------------------------------
def _rewrite_param(block_text: str, pid: str, values: list[float]) -> str:
    """Replace one ``<param>``'s variable list with ``len(values)`` entries.

    Existing lines keep every attribute they had -- only ``value`` changes --
    so bounds, ``is_active`` and ``var_id`` survive untouched.  Columns beyond
    what the base had are cloned from the LAST existing line, which is the one
    Quadoa itself most recently wrote in this shape.
    """
    found = {}

    def rewrite(m: re.Match) -> str:
        if m.group("id") != pid:
            return m.group(0)
        body = m.group("body")
        lines = _VAR_LINE.findall(body)
        template_ms = list(_VAR_LINE.finditer(body))
        if not template_ms:
            raise StructuralSurprise(f"<param id={pid!r}> has no variable lines")
        have = len(template_ms)

        out = []
        for i, val in enumerate(values):
            src = template_ms[i] if i < have else template_ms[-1]
            line = src.group(0)
            if i >= have:
                line = re.sub(r'name="val_\d+"', f'name="val_{i}"', line, count=1)
            new, n = _VALUE_ATTR.subn(f'value="{_fmt(val)}"', line, count=1)
            if n != 1:
                raise StructuralSurprise(
                    f"variable line for {pid} val_{i} has no value attribute: {line!r}"
                )
            out.append(VAR_INDENT + new + "\n")

        found[pid] = (have, len(values))
        return m.group("open") + "\n" + "".join(out) + CLOSE_INDENT + m.group("close")

    new_text = _PARAM_BLOCK.sub(rewrite, block_text)
    if pid not in found:
        raise StructuralSurprise(f"<param id={pid!r}> not found in the block")
    return new_text


def _rewrite_single(block_text: str, pid: str, value: float) -> str:
    hit = []

    def rewrite(m: re.Match) -> str:
        if m.group("id") != pid:
            return m.group(0)
        body, n = _SINGLE_VAR.subn(
            lambda v: _VALUE_ATTR.sub(f'value="{_fmt(value)}"', v.group(0), count=1),
            m.group("body"),
            count=1,
        )
        if n != 1:
            raise StructuralSurprise(f"<single_param id={pid!r}>: {n} val variables")
        hit.append(pid)
        return m.group("open") + body + m.group("close")

    new_text = _SINGLE_BLOCK.sub(rewrite, block_text)
    if not hit:
        raise StructuralSurprise(f"<single_param id={pid!r}> not found in the block")
    return new_text


def build(base_text: str, per_config: dict[str, list[float]],
          globals_: dict[str, float]) -> str:
    """Return the base with the multiconfig block rewritten. Nothing else moves."""
    n = len(next(iter(per_config.values())))
    if any(len(v) != n for v in per_config.values()):
        raise ValueError("every per-config parameter needs the same column count")
    if sorted(per_config) != sorted(PARAMS):
        raise ValueError(f"per-config keys {sorted(per_config)} != {sorted(PARAMS)}")
    if sorted(globals_) != sorted(SINGLE_PARAMS):
        raise ValueError(f"global keys {sorted(globals_)} != {sorted(SINGLE_PARAMS)}")

    start, end = _multiconfig_span(base_text)
    block = base_text[start:end]

    for pid in PARAMS:
        block = _rewrite_param(block, pid, per_config[pid])
    for pid in SINGLE_PARAMS:
        block = _rewrite_single(block, pid, globals_[pid])

    # The column count LAST, so a failure above leaves a block that still
    # agrees with its own header rather than a half-grown one.
    block, k = _MULTICONF.subn(
        lambda m: f'<multiconfig columns="{n}"{m.group(2)}>', block, count=1
    )
    if k != 1:
        raise StructuralSurprise("could not rewrite the multiconfig header")

    return base_text[:start] + block + base_text[end:]


# --------------------------------------------------------------------------
# The physics, all of it borrowed.
# --------------------------------------------------------------------------
def solve_instant(cfg, date: _dt.date, hour: float, flat: bool = False):
    """``(TimeStep, indices, provenance, per_config, globals_)`` for one instant."""
    from beamdown import field as F
    from beamdown.secondary import get_strategy
    from beamdown.solar import TimeStep, sun_position

    az, el = sun_position(
        cfg.site.latitude, cfg.site.longitude, cfg.site.timezone,
        date.year, date.month, date.day, float(hour),
    )
    step = TimeStep(date=date, hour=float(hour), solar_az_deg=az, solar_el_deg=el)

    full = F.load_field(cfg)
    idx, provenance = F.load_or_build_downselect(cfg, full)
    sub = full.subset(idx)

    strategy = get_strategy(cfg, flat=flat) if flat else get_strategy(cfg)

    per_config: dict[str, list[float]] = {p: [] for p in PARAMS}
    for i in range(len(sub)):
        x, y = float(sub.x_mm[i]), float(sub.y_mm[i])
        sol = strategy.solve(x, y, az, el, cfg.geometry)
        per_config["posx"].append(x)
        per_config["posy"].append(y)
        per_config["rot_az"].append(float(sol.rot_az_deg))
        per_config["rot_el"].append(float(sol.rot_el_deg))
        per_config["c3"].append(float(sol.c3))
        per_config["c4"].append(float(sol.c4))
        per_config["c5"].append(float(sol.c5))

    # solze is the ZENITH angle -- the model has no elevation parameter. Getting
    # this backwards points every mirror at the complement of the real sun and
    # nothing errors; session.set_sun writes 90 - elevation, and so does this.
    globals_ = {"solaz": float(az), "solze": float(90.0 - el)}
    globals_.update({k: float(v) for k, v in strategy.global_params(cfg.geometry).items()})

    # A non-axicon strategy on the axicon base is a deliberate use, not an
    # accident: the field gets that layout's POINTING while the base keeps its
    # own scenery. Reconcile the single_params explicitly, out loud, and refuse
    # anything that is not one of the two understood mismatches.
    #   - pf_height (prime focus supplies it, this base has no such param):
    #     dropped. Quadoa would silently ignore the write anyway; the aim
    #     height lives in the per-config pointing/Zernikes, nowhere else. There
    #     is no prime-focus detector in this model.
    #   - axi_angle (this base carries it, non-axicon strategies do not):
    #     filled from config. The cone stays where the base put it -- for a
    #     figure/design model it is scenery to replace in the GUI, but a TRACE
    #     of this file still bounces off it. Do not sweep this model.
    for k in sorted(set(globals_) - set(SINGLE_PARAMS)):
        if k != "pf_height":
            raise StructuralSurprise(
                f"strategy {strategy.describe()!r} supplies {k!r}, which this "
                f"base model has no <single_param> for -- Quadoa ignores writes "
                f"to parameters it does not have, silently. Wrong base model?"
            )
        print(f"  note: dropping {k!r} = {globals_.pop(k):g} -- the base has no "
              f"such single_param; the aim height is carried by the pointing")
    for k in sorted(set(SINGLE_PARAMS) - set(globals_)):
        if k != "axi_angle":
            raise StructuralSurprise(
                f"strategy {strategy.describe()!r} does not supply {k!r}, which "
                f"this base model carries as <single_param>. Writing a model "
                f"parameter Quadoa has but Python never sets leaves it at "
                f"whatever the base happened to hold -- refusing to build."
            )
        globals_[k] = float(cfg.geometry.axicon_angle_deg)
        print(f"  note: base keeps its axicon cone (axi_angle {globals_[k]:g} "
              f"from config) -- scenery for GUI design work; do not trace")

    return step, idx, provenance, sub, per_config, globals_, strategy


# --------------------------------------------------------------------------
# Verification: re-read what was written and prove the values landed.
# --------------------------------------------------------------------------
def read_multiconfig(text: str) -> tuple[int, dict[str, list[float]], dict[str, float]]:
    """Parse a built file back into ``(columns, per_config, globals_)``."""
    start, end = _multiconfig_span(text)
    block = text[start:end]
    cols = int(_MULTICONF.search(text).group(1))

    per_config: dict[str, list[float]] = {}
    for m in _PARAM_BLOCK.finditer(block):
        vals = re.findall(
            r'<variable\s+name="val_(\d+)"[^>]*?value="([^"]*)"', m.group("body")
        )
        ordered = sorted((int(i), float(v)) for i, v in vals)
        if [i for i, _ in ordered] != list(range(len(ordered))):
            raise StructuralSurprise(f"{m.group('id')}: non-contiguous val indices")
        per_config[m.group("id")] = [v for _, v in ordered]

    globals_: dict[str, float] = {}
    for m in _SINGLE_BLOCK.finditer(block):
        v = re.search(r'<variable\s+name="val"[^>]*?value="([^"]*)"', m.group("body"))
        if v is None:
            raise StructuralSurprise(f"{m.group('id')}: no val variable")
        globals_[m.group("id")] = float(v.group(1))

    return cols, per_config, globals_


def verify(base_text: str, text: str, per_config: dict[str, list[float]],
           globals_: dict[str, float]) -> list[tuple[bool, str]]:
    checks: list[tuple[bool, str]] = []

    def add(ok, msg):
        checks.append((bool(ok), msg))

    n = len(per_config["posx"])

    ok, fixed, msg = xml_check(text)
    _, base_fixed, _ = xml_check(base_text)
    add(ok, msg)
    if not ok:
        return checks
    add(fixed == base_fixed,
        f"no new malformed tags: {fixed} duplicate-attribute tags, base has "
        f"{base_fixed}")

    cols, got_p, got_g = read_multiconfig(text)
    add(cols == n, f"<multiconfig columns=\"{cols}\"> (base had "
                   f"{int(_MULTICONF.search(base_text).group(1))})")

    add(list(got_p) == PARAMS, f"per-config params unchanged: {list(got_p)}")
    add(list(got_g) == SINGLE_PARAMS, f"single_params unchanged: {list(got_g)}")

    widths = {k: len(v) for k, v in got_p.items()}
    add(set(widths.values()) == {n},
        f"every <param> has exactly {n} columns: {sorted(set(widths.values()))}")

    # The values themselves, exactly -- repr round-trips, so this is ==, not
    # a tolerance. A tolerance here would hide a formatting bug.
    worst = 0.0
    bad = []
    for pid in PARAMS:
        for i, (want, got) in enumerate(zip(per_config[pid], got_p[pid])):
            if got != want:
                bad.append(f"{pid}[{i}] {got!r} != {want!r}")
            worst = max(worst, abs(got - want))
    add(not bad, f"all {n * len(PARAMS)} per-config values round-trip exactly"
                 + (f" -- {len(bad)} bad, first: {bad[0]}" if bad else
                    f" (max |delta| {worst:g})"))

    gbad = [f"{k} {got_g[k]!r} != {v!r}" for k, v in globals_.items() if got_g[k] != v]
    add(not gbad,
        "all 5 single_params round-trip exactly: "
        + ", ".join(f"{k}={globals_[k]:g}" for k in SINGLE_PARAMS)
        + (f" -- BAD: {gbad}" if gbad else ""))

    # The strongest structural check: everything outside the multiconfig block
    # is byte-identical. This script has no business touching geometry,
    # sequences, sources or coatings, so prove it did not.
    b0, b1 = _multiconfig_span(base_text)
    t0, t1 = _multiconfig_span(text)
    add(base_text[:b0] == text[:t0] and base_text[b1:] == text[t1:],
        f"everything outside <multiconfig> is byte-identical to the base "
        f"({t0:,} bytes before, {len(text) - t1:,} after)")

    tc = re.findall(r'tcuid="([^"]*)"', text)
    btc = re.findall(r'tcuid="([^"]*)"', base_text)
    add(tc == btc and len(tc) == len(set(tc)),
        f"tcuids unchanged and unique: {len(tc)} entries")

    # Distinctness: the failure this whole exercise is about is a file whose
    # configurations are all the same (or all zero). Two heliostats in a
    # 25-way max-dispersion downselect are never at the same place.
    xy = list(zip(got_p["posx"], got_p["posy"]))
    add(len(set(xy)) == n,
        f"all {n} configurations hold DISTINCT positions "
        f"({len(set(xy))} unique) -- not the all-zeros file this replaces")
    nonzero = sum(1 for x, y in xy if (x, y) != (0.0, 0.0))
    add(nonzero == n, f"no configuration is left at the origin: {nonzero}/{n} placed")

    return checks


# --------------------------------------------------------------------------
def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--date", required=True,
                   help="local date, YYYY-MM-DD (e.g. 2026-06-21)")
    p.add_argument("--hour", required=True, type=float,
                   help="local clock hour, decimal (e.g. 12.0, 16.5)")
    p.add_argument("--flat", action="store_true",
                   help="flat heliostats: same pointing, Zernikes zeroed")
    p.add_argument("--secondary", default=None,
                   choices=["axicon", "prime_focus", "cassegrain"],
                   help="override config's layout for the POINTING solve only; "
                        "the base model (and its axicon cone) is unchanged "
                        "scenery. prime_focus/cassegrain also need "
                        "--focus-height-mm.")
    p.add_argument("--focus-height-mm", type=float, default=None,
                   help="override [geometry] focus_height_mm: F1 the whole "
                        "field aims and focuses at (e.g. 36000 = 9 m above "
                        "the axicon tip)")
    p.add_argument("--base", type=Path, default=BASE)
    p.add_argument("--out", type=Path, default=None,
                   help="default models/figure_model_<n>cfg_<key>[_flat].optx")
    p.add_argument("--config", type=Path, default=None, help="config.toml to use")
    p.add_argument("--force", action="store_true", help="overwrite an existing output")
    p.add_argument("--check", action="store_true",
                   help="re-read the output from disk and re-assert every value")
    a = p.parse_args(argv)

    from beamdown.config import load_config, validate_layout

    cfg = load_config(a.config)
    if a.secondary is not None:
        # Same override mechanism as verify_prime_focus_model.py: mutate the
        # loaded cfg rather than requiring a config.toml edit (VALUES in that
        # file must never change while a sweep runs). n_mirrors follows the
        # layout; it only matters to validate_layout here -- a figure model
        # carries no throughput.
        object.__setattr__(cfg.optics, "secondary", a.secondary)
        object.__setattr__(cfg.optics, "n_mirrors",
                           1 if a.secondary == "prime_focus" else 2)
        if a.focus_height_mm is not None:
            object.__setattr__(cfg.geometry, "focus_height_mm",
                               float(a.focus_height_mm))
        validate_layout(cfg)
    elif a.focus_height_mm is not None:
        print("--focus-height-mm without --secondary does nothing for the "
              "axicon; refusing so the intent is explicit.")
        return 2
    date = _dt.date.fromisoformat(a.date)

    step, idx, provenance, sub, per_config, globals_, strategy = solve_instant(
        cfg, date, a.hour, flat=a.flat
    )
    n = len(sub)

    layout_tag = ""
    if a.secondary is not None:
        short = {"prime_focus": "pf", "cassegrain": "cass"}.get(a.secondary,
                                                                a.secondary)
        fh = cfg.geometry.focus_height_mm
        layout_tag = f"_{short}" + (f"{fh:.0f}" if fh is not None else "")
    out = a.out or (REPO / "models" /
                    f"figure_model_{n}cfg_{step.key}{layout_tag}"
                    f"{'_flat' if a.flat else ''}.optx")
    out = Path(out)
    if out.resolve() == PRISTINE.resolve():
        print(f"refusing to write {out.name}: that is the original artefact this "
              f"script exists to explain. Choose another --out.")
        return 2

    print(f"instant: {step}")
    print(f"  site lat {cfg.site.latitude}, lon {cfg.site.longitude}, "
          f"tz {cfg.site.timezone:+g}")
    print(f"  solaz {globals_['solaz']:.4f}  solze {globals_['solze']:.4f} "
          f"(= 90 - elevation)")
    print(f"  layout {strategy.describe()!r}, {n} heliostats, downselect "
          f"{provenance}")
    if step.solar_el_deg <= 0.0:
        print(f"  WARNING: the sun is BELOW the horizon (el "
              f"{step.solar_el_deg:.2f} deg). The pointing solve still returns "
              f"numbers; they are not physical.")
    print()

    base_text = a.base.read_text(encoding="utf-8")
    print(f"base: {a.base}  ({len(base_text):,} bytes)")
    for note in check_base(base_text):
        print(f"  {note}")
    print()

    previous = out.read_text(encoding="utf-8") if out.exists() else None
    if previous is not None and not a.force and not a.check:
        print(f"refusing to overwrite {out} -- pass --force")
        return 2

    text = build(base_text, per_config, globals_)
    out.parent.mkdir(parents=True, exist_ok=True)
    eol = newline_of(a.base)
    with open(out, "w", encoding="utf-8", newline=eol) as fh:
        fh.write(text)

    on_disk = out.stat().st_size
    print(f"wrote {out}  ({on_disk:,} bytes on disk, "
          f"{on_disk - a.base.stat().st_size:+,} over the base; line endings "
          f"{'CRLF' if len(eol) == 2 else 'LF'}, as the base)")
    print(f"  {n} configurations x {len(PARAMS)} params + "
          f"{len(SINGLE_PARAMS)} single_params written")
    if previous is not None:
        print(f"  idempotent: {'identical to' if previous == text else 'DIFFERS from'}"
              f" the previous {out.name}")
    print()

    print(f"  {'cfg':>3} {'x (m)':>9} {'y (m)':>9} {'rot_az':>10} {'rot_el':>9} "
          f"{'c3':>11} {'c4':>11} {'c5':>11}")
    for i in range(n):
        print(f"  {i:>3} {per_config['posx'][i] / 1000:>9.2f} "
              f"{per_config['posy'][i] / 1000:>9.2f} "
              f"{per_config['rot_az'][i]:>10.4f} {per_config['rot_el'][i]:>9.4f} "
              f"{per_config['c3'][i]:>11.4g} {per_config['c4'][i]:>11.4g} "
              f"{per_config['c5'][i]:>11.4g}")
    print()

    if a.check:
        # Deliberately re-read from disk and recompute the expectation from the
        # SAME CLI arguments, rather than comparing against the in-memory
        # dictionaries -- otherwise a bug in _fmt or in the writer that also
        # affects the reader would cancel out and pass.
        text = out.read_text(encoding="utf-8")
        _, _, _, _, want_p, want_g, _ = solve_instant(cfg, date, a.hour, flat=a.flat)
        checks = verify(base_text, text, want_p, want_g)
        print("  --check: re-read from disk, expectation recomputed from "
              "--date/--hour")
    else:
        checks = verify(base_text, text, per_config, globals_)

    bad = 0
    for ok, msg in checks:
        print(f"  [{'ok' if ok else 'FAIL'}] {msg}")
        bad += not ok
    print(f"\n  {len(checks) - bad}/{len(checks)} checks passed")
    print("  " + ("PASS -- next: scripts/verify_figure_model.py, once a seat is free"
                  if not bad else "FAIL -- do not open this model"))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
