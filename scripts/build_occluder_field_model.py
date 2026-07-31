"""Build ``models/heliostat_field_occluders.optx`` from the bulk-sweep model.

The bulk sweep traces one heliostat per configuration and knows nothing about
its neighbours.  This adds a fixed pool of *occluder slots* to that model: 10
on the incoming leg (shading) and 4 on the outgoing leg (blocking), each an
infinite transparent plane carrying a rectangular ``<obscuration>`` the size of
a heliostat facet -- plus one more, ``ax0``, for the secondary's own shadow.

    sun -> ax0 -> so0..so9 -> helio_surf -> bo0..bo3 -> secondary -> receiver

Every slot is placed entirely from multiconfig parameters, so a caller can
position all fifteen at trace time through ``setMulticonfParam`` without ever
touching the file again::

    NAME_x, NAME_y     ground position of the occluding heliostat (mm)
    NAME_az, NAME_el   its azimuth / elevation (deg)  [not on ax0]

They are ``<single_param>`` rather than ``<param>`` entries because the sweep
only ever writes configuration 0; a single_param holds one value for the whole
model in a lone ``<variable name="val">`` child instead of 24 ``val_i`` ones.

Unused slots are *parked*: the defaults put them 200 m away, where the
obscuration -- 2.5 m by 1.5 m -- can intersect nothing.  A parked slot must cost
zero rays, which is the property the verification at the bottom of this file
checks statically and ``traceRays`` checks for real.  See ``PARK_XY`` for why
that distance is 200 m and not the 1000 km you would expect.

Text surgery, not an XML round-trip: Quadoa's own formatting survives, and the
diff against the base model stays readable.

Usage::

    python scripts/build_occluder_field_model.py [--base P] [--out P]
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
BASE = REPO / "models" / "heliostat_field_model_mcfg.optx"
OUT = REPO / "models" / "heliostat_field_occluders.optx"

SHADING = [f"so{i}" for i in range(10)]
BLOCKING = [f"bo{i}" for i in range(4)]
SLOTS = SHADING + BLOCKING

# Every generated uid lives at or above this number.  The base model's own uids
# are single- and double-digit, so nothing can collide even if it grows.
UID0 = 200

# Where a parked slot sits.  200 m: ten times the field radius, so the 2.5 m
# obscuration cannot touch a bundle that never leaves the heliostat's own
# neighbourhood -- and close enough that the tracer stays exact.
#
# Not 1000 km, which is the intuitive choice and is wrong.  A sequential trace
# visits every listed surface, so a slot parked at 1e9 mm sends each ray on a
# 2.6e9 mm round trip before it reaches the mirror, and Quadoa comes back with
# a ray that misses.  Measured on this model, 20k rays, against the unmodified
# base (mean 7440, sd 43):
#
#     park 1e5   +0.06 %      park 5e5  -24.4 %
#     park 2e5   +0.00 %      park 7e5   -4.9 %
#     park 3e5   +0.87 %      park 1e6  -22.5 %
#     park 1e9  -27.4 %
#
# The rays are not stopped by the parked planes -- all 20000 still cross every
# one of them -- they arrive at `helio_surf` displaced enough to miss it. The
# cliff is somewhere between 3e5 and 5e5 and is not monotone past it, so this
# sits a factor 2.5 below the nearest failure rather than just inside it.
PARK_XY = 200000.0

# -- the axicon slot -------------------------------------------------------
# ax0 is the secondary's own shadow on the field, and it is a different animal
# from the fourteen heliostat slots:
#
#   * It is HORIZONTAL.  Every other slot inherits rx=-90/rz=90 on the surface
#     plus the parent's az/el, because it stands upright like the mirror it
#     mimics.  ax0 has no rotation anywhere, so its normal is +z.
#   * It sits at a fixed height, as a literal, not a parameter -- there is only
#     one axicon and it does not move.  `beamdown.occluder_slots.AXICON_PLANE_Z_MM`
#     must agree with this number; that is where the shadow centre is projected.
#   * Its obscuration is a CIRCLE of the axicon rim's radius, not a facet
#     rectangle.  A horizontal circle projected along a fixed direction onto a
#     horizontal plane is a congruent circle, so the shadow is exact.
#
# It also needs no `_az`/`_el`: a circle has no orientation and the plane stays
# horizontal, so two parameters place it completely.
#
# Parking it is free in a way parking the others is not.  Translating an
# infinite *tilted* plane in x/y slides it along its own normal and moves where
# every ray crosses it -- the reason PARK_XY above is 200 m and not 1000 km.
# Translating an infinite *horizontal* plane in x/y leaves the plane itself
# identical; only the obscuration moves.  So ax0 can park as far away as we
# like, and 1 km is chosen simply to be unmistakably off a 90 m field.
AXICON = "ax0"
AX_SLOT = 14            # uid slot index, one past the last heliostat slot
AX_PARK_XY = 1000000.0
AX_PLANE_Z = 13500.0
AX_RADIUS = 15000.0

# A circular obscuration carries its radius in a <float_ap> child, exactly as
# the `secondary` surface's circular *aperture* does -- there is no `radius`
# attribute on the element itself, the way `rect` has `s_x`/`s_y`.
#
# Confirmed two ways against Quadoa 25.09, because a mis-spelled obscuration is
# not an error -- it is a plane that quietly stops nothing:
#
#   * Round-trip.  Loading this file and saving it straight back out re-emits
#     the block below character for character, so the parser kept both the
#     `circ` type and the radius.  (The only thing Quadoa rewrites anywhere in
#     the file is the surfaces' own unset `float_ap name="radius"` clear
#     aperture, 0.0 -> inf, which it does to the base model's surfaces too.)
#   * Functionally.  Sitting on the ray bundle it takes 20,000 rays to 0.
#
# `addObscurationToSurface(row, quadoa.obscurationType.obscurationCirc)` is the
# API route to the same element, if you ever want Quadoa to author one for you.
AX_OBSCURATION = (
    '<obscuration type="circ" c_x="0.0" c_y="0.0" is_visible="true"'
    ' store_hidden="false">\n{i}\t<float_ap name="radius" value="{r}"'
    ' is_locked="true" />\n{i}</obscuration>'
)

# The analysis sequence -- 4th <sequence> element, "sequence 4" in the GUI.
SEQ_INDEX = 3
SEQ_SURFACES = ["sun", "helio_surf", "secondary", "receiver"]

_VAR = ('<variable name="{name}" value="{value}" is_active="0" min="-inf"'
        ' max="inf" bound_type="2" bound_weight="1.0" bound_tolerance="1e-08"'
        ' var_id="" />')


def _var(indent: str, name: str, value) -> str:
    return indent + _VAR.format(name=name, value=value)


def _span_of_pos(text: str, start: int) -> tuple[int, int]:
    """Character span of the ``<pos>`` element beginning at ``start``.

    A depth counter over pos tags, because the blocks nest and a regex for the
    closing tag would stop at the first inner one.
    """
    depth = 0
    for m in re.finditer(r"<pos\b|</pos>", text[start:]):
        depth += 1 if m.group(0) == "<pos" else -1
        if depth == 0:
            return start, start + m.end()
    raise ValueError("unterminated <pos> element")


def occluder_block(name: str, slot: int) -> str:
    """One occluder assembly: outer pos (position), inner pos (angles), surface.

    The split into two nested ``<pos>`` elements is not decoration.  The outer
    one uses ``order="zyx"`` and carries only the translation; the inner one
    uses ``order="xyz"`` and carries only the rotation, so azimuth and elevation
    compose about the occluder's own centre rather than about the field origin.
    That is exactly how ``heli_pos``/``heli_coord_shift`` drive the heliostat
    itself, and it is why ``ry`` takes ``-NAME_el`` rather than ``NAME_el``.
    """
    a_out, a_in = UID0 + 1 + 2 * slot, UID0 + 2 + 2 * slot
    i3, i4, i5, i6 = "\t\t\t", "\t\t\t\t", "\t\t\t\t\t", "\t\t\t\t\t\t"
    L = [
        f'{i3}<pos id="p_{name}" refe="abs" off_axis="true" order="zyx"'
        f' posrot="pos" rot_c_x="0.0" rot_c_y="0.0" rot_c_z="0.0"'
        f' is_visible="true" store_hidden="false" puuid="assembly{a_out}">',
        _var(i4, "z", "0.0"),
        _var(i4, "x", f"{name}_x"),
        _var(i4, "y", f"{name}_y"),
        _var(i4, "rx", "0.0"),
        _var(i4, "ry", "0.0"),
        _var(i4, "rz", "0.0"),
        f'{i4}<pos id="h_{name}" refe="abs" off_axis="true" order="xyz"'
        f' posrot="pos" rot_c_x="0.0" rot_c_y="0.0" rot_c_z="0.0"'
        f' is_visible="true" store_hidden="false" puuid="assembly{a_in}">',
        _var(i5, "z", "0.0"),
        _var(i5, "x", "0.0"),
        _var(i5, "y", "0.0"),
        _var(i5, "rx", "0.0"),
        _var(i5, "ry", f"-{name}_el"),
        _var(i5, "rz", f"{name}_az"),
        f'{i5}<surf id="{name}" refe="abs" off_axis="true" order="xyz"'
        f' posrot="pos" rot_c_x="0.0" rot_c_y="0.0" rot_c_z="0.0"'
        f' is_visible="true" store_hidden="false" puuid="surf{UID0 + slot}"'
        f' suid="surf_uid{UID0 + slot}" ghost="true">',
        _var(i6, "radius", "inf"),
        _var(i6, "z", "0.0"),
        f'{i6}<float_ap name="radius" value="0.0" is_locked="false" />',
        _var(i6, "x", "0.0"),
        _var(i6, "y", "0.0"),
        _var(i6, "rx", "-90"),
        _var(i6, "ry", "0.0"),
        _var(i6, "rz", "90.0"),
        f'{i6}<obscuration type="rect" s_x="2500.0" s_y="1500.0" c_x="0.0"'
        f' c_y="0.0" rot="0.0" is_visible="true" store_hidden="false" />',
        f"{i5}</surf>",
        f"{i4}</pos>",
        f"{i3}</pos>",
    ]
    return "\n".join(L)


def axicon_block() -> str:
    """The ``ax0`` assembly: one flat plane carrying a circular obscuration.

    One ``<pos>``, not the two the heliostat slots use.  The nested pair exists
    only so azimuth and elevation compose about the slot's own centre; with no
    rotation to compose there is nothing for the inner one to do.
    """
    i3, i4, i5, i6 = "\t\t\t", "\t\t\t\t", "\t\t\t\t\t", "\t\t\t\t\t\t"
    a, u = UID0 + 1 + 2 * AX_SLOT, UID0 + AX_SLOT
    L = [
        f'{i3}<pos id="p_{AXICON}" refe="abs" off_axis="true" order="zyx"'
        f' posrot="pos" rot_c_x="0.0" rot_c_y="0.0" rot_c_z="0.0"'
        f' is_visible="true" store_hidden="false" puuid="assembly{a}">',
        # The height is a literal: there is one axicon and it does not move.
        _var(i4, "z", repr(AX_PLANE_Z)),
        _var(i4, "x", f"{AXICON}_x"),
        _var(i4, "y", f"{AXICON}_y"),
        _var(i4, "rx", "0.0"),
        _var(i4, "ry", "0.0"),
        _var(i4, "rz", "0.0"),
        f'{i4}<surf id="{AXICON}" refe="abs" off_axis="true" order="xyz"'
        f' posrot="pos" rot_c_x="0.0" rot_c_y="0.0" rot_c_z="0.0"'
        f' is_visible="true" store_hidden="false" puuid="surf{u}"'
        f' suid="surf_uid{u}" ghost="true">',
        _var(i5, "radius", "inf"),
        _var(i5, "z", "0.0"),
        f'{i5}<float_ap name="radius" value="0.0" is_locked="false" />',
        _var(i5, "x", "0.0"),
        _var(i5, "y", "0.0"),
        # No rotation anywhere -- the plane's normal is +z and stays there.
        _var(i5, "rx", "0.0"),
        _var(i5, "ry", "0.0"),
        _var(i5, "rz", "0.0"),
        i5 + AX_OBSCURATION.format(i=i5, r=repr(AX_RADIUS)),
        f"{i4}</surf>",
        f"{i3}</pos>",
    ]
    return "\n".join(L)


def wrapper_block() -> str:
    """The ``shade_block`` assembly holding all fifteen slots.

    An all-zero container, sibling of ``heli_pos``: it exists so the slots are
    one collapsible group in the GUI tree and so their coordinates stay
    absolute field coordinates rather than heliostat-relative ones.
    """
    i2, i3 = "\t\t", "\t\t\t"
    L = [
        f'{i2}<pos id="shade_block" refe="abs" off_axis="true" order="zyx"'
        f' posrot="pos" rot_c_x="0.0" rot_c_y="0.0" rot_c_z="0.0"'
        f' is_visible="true" store_hidden="false" puuid="assembly{UID0}">',
    ]
    L += [_var(i3, n, "0.0") for n in ("z", "x", "y", "rx", "ry", "rz")]
    L.append(axicon_block())            # first in the sequence, first in the tree
    L += [occluder_block(n, i) for i, n in enumerate(SLOTS)]
    L.append(f"{i2}</pos>")
    return "\n".join(L)


def single_params() -> str:
    """The 58 ``<single_param>`` entries: four per heliostat slot, two for ax0."""
    i2, i3 = "\t\t", "\t\t\t"
    out = []

    def emit(name: str, tc: int, value: float) -> None:
        out.append(f'{i2}<single_param id="{name}" tcuid="tr_{tc}">')
        out.append(_var(i3, "val", repr(value)))
        out.append(f"{i2}</single_param>")

    for slot, name in enumerate(SLOTS):
        for k, (suffix, value) in enumerate(
                (("x", PARK_XY), ("y", PARK_XY), ("az", 0.0), ("el", 0.0))):
            emit(f"{name}_{suffix}", UID0 + 4 * slot + k, value)

    # ax0 takes only x and y -- see AXICON above for why there is no az/el.
    for k, suffix in enumerate(("x", "y")):
        emit(f"{AXICON}_{suffix}", UID0 + 4 * AX_SLOT + k, AX_PARK_XY)
    return "\n".join(out)


def sequence_ref(name: str, slot: int) -> str:
    """One sequence entry for an occluder.

    ``action="transmission"``: the plane must stop the rays its obscuration
    covers without reflecting the rest.
    """
    return (f'\t\t\t\t<surf id="{name}" surf="1" action="transmission"'
            f' diff_order="1" diff_order_2="0" ray_type="ord" refl_mat="env"'
            f' jones="Auto" is_visible="true" store_hidden="false"'
            f' spuid="surfp_uid{UID0 + slot}" hide_to="false"'
            f' hide_from="false" />')


def build(base: Path = BASE, out: Path = OUT) -> str:
    text = base.read_text(encoding="utf-8")

    # -- geometry: shade_block goes in as a sibling of heli_pos --------------
    _, end = _span_of_pos(text, text.index('<pos id="heli_pos"'))
    text = text[:end] + "\n" + wrapper_block() + text[end:]

    # -- multiconfig: the 58 single_params, appended to the block -----------
    close = text.index("\t</multiconfig>")
    text = text[:close] + single_params() + "\n" + text[close:]

    # -- sequence 4 (index 3): wrap the occluders around helio_surf ---------
    seqs = list(re.finditer(r"<sequence\b.*?</sequence>", text, re.S))
    seq = seqs[SEQ_INDEX]
    body = seq.group(0)
    listed = re.findall(r'<surf id="([^"]*)"', body)
    if listed != SEQ_SURFACES:
        raise ValueError(f"sequence {SEQ_INDEX} is {listed}, expected "
                         f"{SEQ_SURFACES} -- refusing to edit the wrong one")

    # ax0 leads: the axicon's shadow is cast before the light reaches anything
    # else, so it belongs between the source and the first neighbour.
    shade = "\n".join([sequence_ref(AXICON, AX_SLOT)]
                      + [sequence_ref(n, i) for i, n in enumerate(SHADING)])
    block = "\n".join(sequence_ref(n, i + len(SHADING))
                      for i, n in enumerate(BLOCKING))
    helio = re.search(r'[ \t]*<surf id="helio_surf"[^>]*/>\n', body)
    body = (body[:helio.start()] + shade + "\n"
            + helio.group(0) + block + "\n" + body[helio.end():])

    # `stopsurface` is an index INTO THE SEQUENCE, so inserting 15 surfaces ahead
    # of the receiver silently re-points it. It named `secondary`; left alone it
    # would name a parked obscuration, and the aperture stop governs how the
    # source samples the pupil -- the trace still runs and quietly loses 20-35%
    # of its rays, on every heliostat, occluded or not.
    ids = re.findall(r'<surf id="([^"]*)"\s+surf=', body)
    head_end = body.index(">") + 1
    stop = ids.index("secondary")
    head, n = re.subn(r'stopsurface="\d+"', f'stopsurface="{stop}"', body[:head_end])
    if n != 1:
        raise ValueError("could not re-point stopsurface on the analysis sequence")
    body = head + body[head_end:]

    text = text[:seq.start()] + body + text[seq.end():]

    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text, encoding="utf-8")
    return text


# --------------------------------------------------------------------------
# Static verification.  Cheap enough to run on every build, and it catches the
# whole class of copy-paste uid collisions that Quadoa reports only as a
# silently mis-drawn model.
# --------------------------------------------------------------------------

def verify(base_text: str, text: str) -> list[tuple[bool, str]]:
    checks: list[tuple[bool, str]] = []

    def n(pat, t):
        return len(re.findall(pat, t))

    # 1. tag balance.  Two <pos> per heliostat slot, one for the wrapper, and
    #    one for ax0 -- which needs no inner rotation pos.
    want_pos = 2 * len(SLOTS) + 1 + 1
    want_surf = len(SLOTS) + 1
    popen, pclose = n(r"<pos\b", text), n(r"</pos>", text)
    dpos = popen - n(r"<pos\b", base_text)
    dsurf = n(r"</surf>", text) - n(r"</surf>", base_text)
    dobsc = n(r"<obscuration\b", text)
    ok = (popen == pclose and dpos == want_pos and dsurf == want_surf
          and dobsc == n(r"</obscuration>", text) + len(SLOTS))
    checks.append((ok, f"tag balance: <pos> {popen} open / {pclose} close "
                       f"(+{dpos} vs base, want {want_pos}); "
                       f"<surf> containers +{dsurf} (want {want_surf}); "
                       f"<obscuration> {dobsc} ({len(SLOTS)} self-closing rect "
                       f"+ {n(r'</obscuration>', text)} circ)"))

    # 2. no duplicate puuid / suid anywhere; no duplicate spuid within a sequence
    dup = []
    for attr in ("puuid", "suid"):
        v = re.findall(attr + r'="([^"]*)"', text)
        dup += [f"{attr}={x}" for x in {y for y in v if v.count(y) > 1}]
    for i, m in enumerate(re.finditer(r"<sequence\b.*?</sequence>", text, re.S)):
        v = re.findall(r'spuid="([^"]*)"', m.group(0))
        dup += [f"seq{i} spuid={x}" for x in {y for y in v if v.count(y) > 1}]
    checks.append((not dup, f"unique ids: {'OK' if not dup else dup}"))

    # 3. every single_param present and referenced exactly once
    want = ([f"{s}_{k}" for s in SLOTS for k in ("x", "y", "az", "el")]
            + [f"{AXICON}_x", f"{AXICON}_y"])
    have = re.findall(r'<single_param id="([^"]*)"', text)
    refs = re.findall(r'<variable name="(?:x|y|ry|rz)" value="(-?[a-z]\w*)"', text)
    missing = [p for p in want if p not in have]
    badrefs = [p for p in want
               if refs.count(p) + refs.count("-" + p) != 1]
    ok = not missing and not badrefs and len(have) == 5 + len(want)
    checks.append((ok, f"single_params: {len(have)} total "
                       f"(5 base + {len(want)} new), missing={missing or 'none'}, "
                       f"mis-referenced={badrefs or 'none'}"))

    # 4. the analysis sequence reads in the right order
    seq = list(re.finditer(r"<sequence\b.*?</sequence>", text, re.S))[SEQ_INDEX]
    listed = re.findall(r'<surf id="([^"]*)"', seq.group(0))
    want_seq = (["sun", AXICON] + SHADING + ["helio_surf"] + BLOCKING
                + ["secondary", "receiver"])
    checks.append((listed == want_seq, f"sequence {SEQ_INDEX}: "
                                       f"{' '.join(listed)}"))

    # 4b. the aperture stop must still be `secondary`, by name and not by luck
    stop = int(re.search(r'stopsurface="(\d+)"', seq.group(0)).group(1))
    at = listed[stop] if stop < len(listed) else "(out of range)"
    checks.append((at == "secondary",
                   f"stopsurface={stop} -> {at}"))

    # the other three sequences must be byte-identical to the base
    a = re.findall(r"<sequence\b.*?</sequence>", base_text, re.S)
    b = re.findall(r"<sequence\b.*?</sequence>", text, re.S)
    same = [i for i in range(len(a)) if i != SEQ_INDEX and a[i] == b[i]]
    checks.append((len(same) == len(a) - 1,
                   f"other sequences untouched: {same}"))
    return checks


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--base", type=Path, default=BASE)
    p.add_argument("--out", type=Path, default=OUT)
    a = p.parse_args(argv)

    base_text = a.base.read_text(encoding="utf-8")
    text = build(a.base, a.out)

    print(f"wrote {a.out}  ({len(text):,} bytes, "
          f"+{len(text) - len(base_text):,} over the base)")
    print(f"  {len(SHADING)} shading slots {SHADING[0]}..{SHADING[-1]}, "
          f"{len(BLOCKING)} blocking slots {BLOCKING[0]}..{BLOCKING[-1]}")
    print(f"  {4 * len(SLOTS)} single_params, all parked at "
          f"x=y={PARK_XY:.0f} mm, az=el=0")
    print(f"  1 axicon slot {AXICON}: horizontal plane at z={AX_PLANE_Z:.0f} mm, "
          f"circular obscuration r={AX_RADIUS:.0f} mm,")
    print(f"    2 single_params {AXICON}_x/{AXICON}_y parked at "
          f"x=y={AX_PARK_XY:.0f} mm")
    print(f"  uids: assembly{UID0}..assembly{UID0 + 1 + 2 * AX_SLOT}, "
          f"surf/surf_uid/surfp_uid {UID0}..{UID0 + AX_SLOT}, "
          f"tcuid tr_{UID0}..tr_{UID0 + 4 * AX_SLOT + 1}")
    print()
    bad = 0
    for ok, msg in verify(base_text, text):
        print(f"  [{'ok' if ok else 'FAIL'}] {msg}")
        bad += not ok
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
