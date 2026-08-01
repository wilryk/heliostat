"""Build a cassegrain model: swap the axicon cone for the settled hyperboloid.

NO licence. Same guarded-text-surgery contract as build_prime_focus_model.py:
every edit is anchored, every anchor asserted, the base never modified, the
output re-read and checked. What it does to the base model:

  secondary surface
    z: the named ``sec_height`` reference (or whatever literal is there)
       becomes the LITERAL hyperboloid vertex, z = 27,151.783 mm.  The vertex
       is not the axicon tip height, so the parameter wiring cannot stay.
    form: ``<form type="axicon">`` (or a previous asphere) becomes
       ``<form type="asphere">`` with radius = +31,548.867 mm and
       conic K = -6.582109, both frozen (is_active="0").  Positive radius =
       bowl opening up toward F1, vertex below rim -- the convention the
       owner's own GUI edit established and design_cassegrain.py documents.
    aperture: circular, locked, 15,000 mm (the 30 m manufacturability cap).
  receiver surface
    becomes ABSOLUTE z = 7,000.0 mm, frozen.  It was ``rec_offset`` relative
    to the secondary, and the secondary just moved by +151.8 mm.

  CONSEQUENCE, on purpose: the model keeps its ``sec_height``/``rec_offset``
  (and ``axi_angle``) single_params but NO SURFACE references them any more.
  ``QuadoaSession.set_global_geometry`` still writes them every session;
  the writes land in the parameters and move nothing.  The literals above
  are the single source of truth.  Do NOT run an axicon sweep against this
  file -- the cone is gone.

The conic constants are RECOMPUTED here from design_cassegrain.close_design
(rim 15,000 mm at z 30,000, F1 36,000, F2 7,000) rather than pasted, so this
builder can never drift from the certified design.  Design rationale:
scan_prime_focus_height.py (blocking), scan_cassegrain_annual.py (energy).

Usage::

    python scripts/build_cassegrain_model.py            # sweep model
    python scripts/build_cassegrain_model.py --base models/figure_model_25cfg_20260220_0927_pf36000.optx \\
        --out models/figure_model_25cfg_20260220_0927_cass30.optx
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

import numpy as np  # noqa: E402

from build_prime_focus_model import StructuralSurprise, newline_of, xml_check  # noqa: E402
from design_cassegrain import close_design, read_config, read_field_mm  # noqa: E402

# The settled design, and the defaults of --rim-z-mm / --f1-mm. F2 and the rim
# RADIUS are plant, not choice: the receiver is where it is and the 15 m rim is
# the 30 m manufacturability cap, so neither is exposed as a flag.
RIM_Z, F1_Z, F2_Z, RIM_R = 30000.0, 36000.0, 7000.0, 15000.0

_SURF = "<variable name=\"{n}\" value=\"{v}\" is_active=\"0\""


def surf_span(text: str, sid: str) -> tuple[int, int]:
    m = re.search(rf'<surf id="{sid}"[^>]*>', text)
    if m is None:
        raise StructuralSurprise(f'no <surf id="{sid}"> in the base')
    end = text.find("</surf>", m.end())
    if end < 0:
        raise StructuralSurprise(f'<surf id="{sid}"> never closes')
    return m.start(), end + len("</surf>")


def set_var(block: str, name: str, value: str, what: str) -> str:
    """Rewrite one <variable name=...> line: new value, frozen."""
    pat = re.compile(rf'(<variable name="{name}" )value="[^"]*" is_active="[^"]*"')
    block, n = pat.subn(rf'\1value="{value}" is_active="0"', block, count=1)
    if n != 1:
        raise StructuralSurprise(f"{what}: <variable name={name!r}> not found")
    return block


def fmt(v: float) -> str:
    return repr(round(float(v), 9))


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--base", type=Path,
                   default=REPO / "models" / "heliostat_field_model_mcfg.optx")
    p.add_argument("--out", type=Path,
                   default=REPO / "models" / "heliostat_field_cassegrain.optx")
    p.add_argument("--force", action="store_true")
    p.add_argument("--rim-z-mm", type=float, default=RIM_Z,
                   help=f"height of the {RIM_R:,.0f} mm rim circle, i.e. how high "
                        f"the dish's edge hangs (default {RIM_Z:,.0f}, the settled "
                        f"design)")
    p.add_argument("--f1-mm", type=float, default=F1_Z,
                   help=f"prime focus F1: the one point the whole field aims and "
                        f"focuses at, which the dish relays to the receiver "
                        f"(default {F1_Z:,.0f}, the settled design)")
    a = p.parse_args(argv)

    # The design, recomputed from close_design -- never pasted -- so this builder
    # cannot drift from beamdown.design_eval, which is what the GUI's Design tab
    # shows. Non-default rim/F1 are exactly as certified: close_design refuses a
    # geometry that does not close (F1 at or below the rim, degenerate conic).
    cfg = read_config(REPO / "config.toml")
    x_mm, y_mm = read_field_mm(REPO / cfg["positions_file"])
    rim_z, f1_z = float(a.rim_z_mm), float(a.f1_mm)
    try:
        des = close_design(rim_z, float(np.hypot(x_mm, y_mm).max()), F2_Z, RIM_R,
                           f1_z)
    except ValueError as exc:
        print(f"cannot close the design at rim z {rim_z:,.0f}, F1 {f1_z:,.0f}: {exc}")
        return 2
    K = -(des.c / des.a) ** 2
    R_v = des.b2 / des.a
    vertex = des.z_c + des.a
    settled = (rim_z == RIM_Z and f1_z == F1_Z)
    print(f"design: vertex z {vertex:,.3f}  |R_v| {R_v:,.3f}  K {K:.9f}  "
          f"rim {RIM_R:,.0f} at z {rim_z:,.0f}  F1 {f1_z:,.0f}  F2 {F2_Z:,.0f}"
          + ("" if settled else "   [NOT the settled 30,000/36,000 design]"))

    base_text = a.base.read_text(encoding="utf-8")
    ok, _, msg = xml_check(base_text)
    if not ok:
        raise StructuralSurprise(f"base does not parse: {msg}")
    print(f"base: {a.base.name} ({len(base_text):,} bytes; {msg})")

    if a.out.exists() and not a.force:
        print(f"refusing to overwrite {a.out} -- pass --force")
        return 2

    # ---- secondary ------------------------------------------------------
    s0, s1 = surf_span(base_text, "secondary")
    sec = base_text[s0:s1]

    sec = set_var(sec, "z", fmt(vertex), "secondary")

    fm = re.search(r"([ \t]*)<form type=\"(axicon|asphere)\"[^>]*>.*?</form>",
                   sec, re.S)
    if fm is None:
        raise StructuralSurprise("secondary has no axicon/asphere form to replace")
    ind = fm.group(1)
    vind = ind + "\t"
    tail = ('min="-inf" max="inf" bound_type="2" bound_weight="1.0" '
            'bound_tolerance="1e-08" var_id=""')
    form = (
        f'{ind}<form type="asphere" maxorder="0" is_visible="true" '
        f'store_hidden="false" active="true">\n'
        f'{vind}<variable name="radius" value="{fmt(R_v)}" is_active="0" {tail} />\n'
        f'{vind}<variable name="conic" value="{fmt(K)}" is_active="0" {tail} />\n'
        f'{ind}</form>'
    )
    sec = sec[:fm.start()] + form + sec[fm.end():]
    print(f"  secondary: z -> {fmt(vertex)}, form {fm.group(2)!r} -> asphere "
          f"(radius {fmt(R_v)}, conic {fmt(K)}, frozen)")

    # The locked circular aperture IS the 30 m cap; make sure it survives.
    ap = re.search(r'<aperture type="circ"[^>]*>.*?</aperture>', sec, re.S)
    if ap is None:
        ins = sec.rfind("</surf>")
        block = (f'{ind}<aperture type="circ" c_x="0.0" c_y="0.0" '
                 f'is_visible="true" store_hidden="false">\n'
                 f'{vind}<float_ap name="radius" value="{fmt(RIM_R)}" '
                 f'is_locked="true" />\n{ind}</aperture>\n\t\t')
        sec = sec[:ins] + block + sec[ins:]
        print("  secondary: circular aperture was missing -- restored, locked 15,000")
    else:
        apx, n = re.subn(r'(<float_ap name="radius" )value="[^"]*"',
                         rf'\1value="{fmt(RIM_R)}"', ap.group(0), count=1)
        if n != 1:
            raise StructuralSurprise("aperture block has no radius float_ap")
        sec = sec[:ap.start()] + apx + sec[ap.end():]
        print("  secondary: circular aperture kept at 15,000")

    text = base_text[:s0] + sec + base_text[s1:]

    # ---- receiver -------------------------------------------------------
    r0, r1 = surf_span(text, "receiver")
    rec = text[r0:r1]
    rec, n = re.subn(r'(<surf id="receiver" )refe="rel"', r'\1refe="abs"', rec)
    rel_note = "refe rel -> abs, " if n else ""
    rec = set_var(rec, "z", fmt(F2_Z), "receiver")
    text = text[:r0] + rec + text[r1:]
    print(f"  receiver: {rel_note}z -> {fmt(F2_Z)} (absolute, frozen)")

    # ---- write + verify -------------------------------------------------
    eol = newline_of(a.base)
    with open(a.out, "w", encoding="utf-8", newline=eol) as fh:
        fh.write(text)

    back = a.out.read_text(encoding="utf-8")
    ok, _, msg = xml_check(back)
    checks = [
        (ok, msg),
        (f'value="{fmt(vertex)}"' in back, "vertex z landed"),
        (f'value="{fmt(R_v)}"' in back, "radius landed"),
        (f'value="{fmt(K)}"' in back, "conic landed"),
        ('"sec_height"' not in re.search(r'<surf id="secondary".*?</surf>', back, re.S).group(0),
         "secondary no longer references sec_height"),
        ('"rec_offset"' not in re.search(r'<surf id="receiver".*?</surf>', back, re.S).group(0),
         "receiver no longer references rec_offset"),
        ('<form type="axicon"' not in back, "no axicon form remains"),
        (back[:s0] == base_text[:s0], "bytes before the secondary block identical"),
    ]
    bad = 0
    for okc, m in checks:
        print(f"  [{'ok' if okc else 'FAIL'}] {m}")
        bad += not okc
    print(f"\nwrote {a.out} ({a.out.stat().st_size:,} bytes)")
    print("NOTE: sec_height/rec_offset/axi_angle single_params remain in the "
          "file but nothing references them; session writes to them are no-ops."
          "\nNext, when a seat is free: verify with scripts/verify_figure_model.py"
          " (25cfg) or a self_test trace before the first sweep.")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
