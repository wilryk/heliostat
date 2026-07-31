"""Build ``models/heliostat_field_prime_focus.optx`` from the bulk-sweep model.

The prime-focus layout has no secondary at all: one reflection off ``helio_surf``
onto a horizontal detector facing down, sitting physically at the field's common
aim point ``F1 = (0, 0, focus_height_mm)``.  The shipped model already carries
the surface -- ``prime_focus``, a flat absolute-positioned plane -- and already
carries the path, as its own sequence 0.  What it does not carry is the wiring:
the surface's height is a *literal* ``27000.0`` that no parameter tracks, and the
sequence the analysis stack actually traces (index 3, ``config.toml``
``analysis_seq``) still runs the axicon path.  This script fixes both, on a copy.

Three edits, and nothing else:

1. ``prime_focus``'s ``z`` becomes the new ``<single_param>`` ``pf_height``,
   default 47000.0 -- the same mechanism ``secondary`` uses for ``sec_height``
   and ``receiver`` for ``rec_offset``, so
   :meth:`PrimeFocusStrategy.global_params` can write it once per session and the
   detector cannot drift away from the Python aim point.
2. ``prime_focus``'s ``float_ap radius`` becomes 2500.0, purely so the surface is
   drawn at a sensible size instead of as a zero-radius dot.  It does NOT clip --
   see ``DRAW_RADIUS``.  The detector stays unbounded on purpose: spillage is
   computed in post from the stored rays, against ``[receiver] window_mm``.
3. Sequence 3 becomes ``sun -> helio_surf -> prime_focus``.

Why 47000
---------
Symmetric throw.  The axicon layout's receiver sits 20 m BELOW the axicon vertex
(``sec_height`` 27000 + ``rec_offset`` -20000 = 7000 mm).  The prime-focus
receiver is that same 20 m throw taken the other way, ABOVE the vertex: 27000 +
20000 = 47000 mm.  The two layouts then differ by the secondary and by one
reflection, not by how far the light travels.

Why sequence 3 and not sequence 0
---------------------------------
Sequence 0 already runs this exact path, and it is left untouched.  But
``analysis_seq = 3`` in ``config.toml`` stays valid for every model file this way
-- no new config key, no per-model branch in the trace stack -- and sequence 3's
``setRayDistributionCount1`` semantics (a LITERAL ray count) are the ones this
repository has measured and relies on.  Sequence 0 is documented as literal too,
but it is the never-actually-traced column, so that is an inherited claim rather
than a measured one.  See "Notes on Quadoa" in README.md.

Text surgery, not an XML round-trip, and not the API: Quadoa's Python API cannot
create geometry, and re-emitting the file through ElementTree would reformat
every line and bury a three-line change in a 65,000-byte diff.  Same approach and
same loud-assertions style as ``scripts/build_occluder_field_model.py``.

Usage::

    python scripts/build_prime_focus_model.py [--base P] [--out P] [--force]
"""

from __future__ import annotations

import argparse
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
BASE = REPO / "models" / "heliostat_field_model_mcfg.optx"
OUT = REPO / "models" / "heliostat_field_prime_focus.optx"

DETECTOR = "prime_focus"

# The new single_param.  Named like its siblings (`sec_height`, `rec_offset`)
# and referenced the same way: a geometry <variable>'s `value` attribute holds
# the parameter's *id* as a bare string instead of a number.  That is the whole
# mechanism -- there is no var_id / uid indirection to replicate.  `secondary`
# carries `value="sec_height"` and `receiver` carries `value="rec_offset"` in the
# base model, and the occluder builder created fifty-eight new single_params
# referenced exactly this way, so the pattern is established twice over.
PF_PARAM = "pf_height"
PF_HEIGHT = 47000.0

# tcuid must be unique across the multiconfig block.  The base model uses tr_1..
# tr_12 and nothing else; the occluder builder's generated ones start at tr_200.
# tr_13 can therefore collide with neither, now or after either builder grows.
PF_TCUID = "tr_13"

# Drawing size only.  A surface's `<float_ap name="radius">` is its *clear
# aperture* for display and auto-sizing; what actually stops a ray is a separate
# `<aperture>` child element, which `secondary` has (circ, r=15000) and which
# neither `receiver` nor `prime_focus` has.
#
# Measured, not assumed: `receiver` carries float_ap radius 795.76 mm, and in
# analysis_output/full7, timestep 20260320_0800, 21.5% of the rays that landed on
# it are further out than that -- reaching 2510 mm, i.e. the corners of the
# +/-2000 mm storage window.  A float_ap that let 21.5% of the bundle through
# beyond its own radius is not clipping anything.
#
# 2500 rather than the receiver's own 795.76 because this detector has to look
# like what it is: a window comfortably larger than `[receiver] window_mm` =
# 2000, so nothing about the drawing suggests a boundary the physics does not
# have.  `is_locked="false"` copies `receiver` exactly -- an unlocked float_ap is
# Quadoa's to auto-size, and letting it do so is the point: the number is a hint,
# not a constraint.
DRAW_RADIUS = 2500.0

# What the base model must look like before any of this is safe.
BASE_Z = "27000.0"
BASE_FLOAT_AP = "0.0"
PARAMS = ["rot_az", "rot_el", "posx", "posy", "c3", "c4", "c5"]
SINGLE_PARAMS = ["solaz", "solze", "sec_height", "rec_offset", "axi_angle"]

# The analysis sequence -- 4th <sequence> element, "sequence 4" in the GUI,
# `config.toml` [trace] analysis_seq = 3.
SEQ_INDEX = 3
SEQ_BEFORE = ["sun", "helio_surf", "secondary", "receiver"]
SEQ_AFTER = ["sun", "helio_surf", DETECTOR]

# Sequence 0, Quadoa's own prime-focus sequence.  Not edited -- used as the
# reference the rewritten sequence 3 is checked against, attribute by attribute.
REF_SEQ_INDEX = 0
REF_NAME = "prime_focus"
REF_SEQUID = "seq_uid2"
SEQ_NAME = "secondary_focus"
SEQ_SEQUID = "seq_uid3"

_SINGLE_PARAM = (
    '\t\t<single_param id="{id}" tcuid="{tcuid}">\n'
    '\t\t\t<variable name="val" value="{value}" is_active="0" min="-inf"'
    ' max="inf" bound_type="2" bound_weight="1.0" bound_tolerance="1e-08"'
    ' var_id="" />\n'
    '\t\t</single_param>\n'
)


class StructuralSurprise(RuntimeError):
    """The base model is not what this script was written against."""


# Quadoa writes one tag that strict XML rejects.  `<opti_settings_plot>` carries
# `axis_type="log"` TWICE:
#
#     ... y_scale_g="full" axis_type="log" axis_type="log" is_visible="true" ...
#
# ElementTree refuses the whole document over it ("duplicate attribute", line
# 521 of the base model), and it is in the optimiser's plot settings -- nothing
# to do with geometry, sequences or the multiconfig block.  It is pre-existing,
# this script does not go near it, and `scripts/build_occluder_field_model.py`
# never noticed because its checks are all textual.
#
# Rather than drop the XML check, the duplicate is collapsed before parsing and
# the collapse is COUNTED: the base and the output must need exactly the same
# number of them, so a well-formedness defect introduced by an edit here still
# fails loudly instead of hiding behind Quadoa's own.
_ATTR = re.compile(r'\s+([A-Za-z_][\w.-]*)="[^"]*"')


def _dedup_attrs(text: str) -> tuple[str, int]:
    """Drop repeated attributes within each tag.  Returns (text, tags fixed)."""
    fixed = 0
    out = []
    pos = 0
    for tag in re.finditer(r"<[A-Za-z_][^>]*>", text):
        seen: set[str] = set()
        body = tag.group(0)
        keep, last = [], 0
        for m in _ATTR.finditer(body):
            if m.group(1) in seen:
                keep.append(body[last:m.start()])
                last = m.end()
            seen.add(m.group(1))
        if keep:
            fixed += 1
            out.append(text[pos:tag.start()])
            out.append("".join(keep) + body[last:])
            pos = tag.end()
    out.append(text[pos:])
    return "".join(out), fixed


def xml_check(text: str) -> tuple[bool, int, str]:
    """``(parses?, tags de-duplicated, message)``."""
    clean, fixed = _dedup_attrs(text)
    try:
        root = ET.fromstring(clean)
    except ET.ParseError as exc:
        return False, fixed, f"XML PARSE ERROR: {exc}"
    n = sum(1 for _ in root.iter())
    return True, fixed, (f"XML well-formed: {n} elements "
                         f"({fixed} tag{'s' if fixed != 1 else ''} needed a "
                         f"duplicate attribute collapsed first -- Quadoa's own, "
                         f"see _dedup_attrs)")


def _span_of_surf(text: str, start: int) -> tuple[int, int]:
    """Character span of the ``<surf>`` element beginning at ``start``.

    A depth counter rather than a search for the next ``</surf>``: surfaces are
    not nested in *this* block, but they are elsewhere in the file, and a builder
    that only works because of where it happens to be pointed is a builder that
    breaks silently the first time the model grows.
    """
    depth = 0
    for m in re.finditer(r"<surf\b|</surf>|/>", text[start:]):
        tok = m.group(0)
        if tok == "<surf":
            depth += 1
        elif tok == "</surf>":
            depth -= 1
        elif depth == 1 and text[start + m.start() - 1] != '"':
            # A self-closing <surf ... /> never happens in the geometry block,
            # but if it ever did, the depth counter would run to the end of the
            # file looking for a close tag it will not find.
            continue
        if depth == 0:
            return start, start + m.end()
    raise StructuralSurprise("unterminated <surf> element")


def newline_of(path: Path) -> str:
    r"""The file's own line ending, so writing it back does not change all 532.

    ``Path.write_text`` on Windows translates every ``\n`` to ``\r\n``, and the
    base model is LF throughout -- so the obvious round trip rewrites every line
    in the file and buries a three-line change under a whole-file diff.  (That is
    what happened to ``models/heliostat_field_occluders.optx``, which is CRLF
    while its own base is LF.  Quadoa reads both perfectly well -- that model is
    what the sweep is tracing right now -- so this is about the diff being
    readable, not about the file being valid.)
    """
    return "\r\n" if b"\r\n" in path.read_bytes() else "\n"


def _geometry(text: str) -> tuple[int, int]:
    """Span of the ``<geometry>`` block.

    Needed because ``<surf id="prime_focus"`` matches in two unrelated places:
    the geometry, where it is the surface, and every sequence that lists it,
    where it is a reference.  A search over the whole file finds both and a
    builder that edits "the first one" is one Quadoa re-save away from editing
    the wrong one.
    """
    start = text.index("\t<geometry>")
    return start, text.index("\t</geometry>") + len("\t</geometry>")


def _detector_span(text: str) -> tuple[int, int]:
    g0, g1 = _geometry(text)
    hits = [g0 + m.start()
            for m in re.finditer(r'<surf id="%s" refe=' % DETECTOR, text[g0:g1])]
    if len(hits) != 1:
        raise StructuralSurprise(
            f"expected exactly one <surf id={DETECTOR!r}> in the geometry, "
            f"found {len(hits)} -- refusing to guess which one to move")
    return _span_of_surf(text, hits[0])


def _sequences(text: str) -> list[re.Match]:
    return list(re.finditer(r"<sequence\b.*?</sequence>", text, re.S))


def _normalise_seq(body: str) -> str:
    """Sequence 3's body, rewritten so it should equal sequence 0's exactly.

    Three attributes differ between the two by design and stay different in the
    output: the sequence's name, its uid, and the visibility flag that makes it
    the one the GUI shows.  Everything else -- the surface list, the aperture
    type, the stop index, and the whole ``<source>`` sub-tree with its 3500 mm
    disk, 0.0024 rad spread and 38484.5 W -- must match, and this is how that is
    asserted rather than eyeballed.
    """
    return (body.replace(f'name="{SEQ_NAME}"', f'name="{REF_NAME}"')
                .replace(f'sequid="{SEQ_SEQUID}"', f'sequid="{REF_SEQUID}"')
                .replace('is_visible="true"', 'is_visible="false"'))


# --------------------------------------------------------------------------
# Pre-flight.  Everything this script assumes about the base model, asserted
# before a single byte is changed -- because every one of these assumptions,
# violated, produces a file that loads fine and traces the wrong thing.
# --------------------------------------------------------------------------

def check_base(text: str) -> list[str]:
    notes: list[str] = []

    ok, fixed, msg = xml_check(text)
    if not ok:
        raise StructuralSurprise(f"base model does not parse: {msg}")
    notes.append(msg)

    start, end = _detector_span(text)
    block = text[start:end]
    notes.append(f"{DETECTOR} surface: {end - start} bytes, "
                 f"lines {text[:start].count(chr(10)) + 1}-{text[:end].count(chr(10)) + 1}")

    z = re.search(r'<variable name="z" value="([^"]*)"([^>]*)>', block)
    if z is None or z.group(1) != BASE_Z:
        raise StructuralSurprise(
            f"{DETECTOR}'s z is {z and z.group(1)!r}, expected the literal "
            f"{BASE_Z!r} -- someone has already edited this surface")
    if 'is_active="0"' not in z.group(2) or 'var_id=""' not in z.group(2):
        raise StructuralSurprise(
            f"{DETECTOR}'s z is not an inert literal: {z.group(0)!r}")
    notes.append(f"  z = {BASE_Z} (literal, is_active=0, no var_id)")

    ap = re.search(r'<float_ap name="radius" value="([^"]*)" is_locked="([^"]*)" />',
                   block)
    if ap is None or ap.group(1) != BASE_FLOAT_AP:
        raise StructuralSurprise(
            f"{DETECTOR}'s float_ap is {ap and ap.group(1)!r}, expected "
            f"{BASE_FLOAT_AP!r}")
    notes.append(f"  float_ap radius = {BASE_FLOAT_AP} "
                 f"(is_locked={ap.group(2)!r}) -- unset, drawn as a dot")

    # The whole point of the detector: nothing clips it.  An <aperture> child
    # WOULD, and adding one would silently truncate the spot and make every
    # spillage number in post agree with itself and disagree with reality.
    if "<aperture" in block:
        raise StructuralSurprise(
            f"{DETECTOR} already carries an <aperture> element, which clips "
            f"rays -- an unbounded detector is what this layout needs")
    notes.append("  no <aperture> child: unbounded, nothing is clipped")

    # multiconfig
    mc = re.search(r"<multiconfig\b[^>]*>", text)
    cols = re.search(r'columns="(\d+)"', mc.group(0)).group(1)
    have_p = re.findall(r'<param id="([^"]*)"', text)
    have_s = re.findall(r'<single_param id="([^"]*)"', text)
    if have_p != PARAMS or have_s != SINGLE_PARAMS:
        raise StructuralSurprise(
            f"multiconfig is params={have_p} single_params={have_s}, expected "
            f"{PARAMS} / {SINGLE_PARAMS}")
    if PF_TCUID in re.findall(r'tcuid="([^"]*)"', text):
        raise StructuralSurprise(f"tcuid {PF_TCUID} is already taken")
    notes.append(f"multiconfig: columns={cols}, {len(have_p)} params, "
                 f"{len(have_s)} single_params, {PF_TCUID} free")

    seqs = _sequences(text)
    if len(seqs) != 4:
        raise StructuralSurprise(f"expected 4 sequences, found {len(seqs)}")
    listed = re.findall(r'<surf id="([^"]*)"', seqs[SEQ_INDEX].group(0))
    if listed != SEQ_BEFORE:
        raise StructuralSurprise(
            f"sequence {SEQ_INDEX} is {listed}, expected {SEQ_BEFORE} -- "
            f"refusing to edit the wrong one")
    ref = re.findall(r'<surf id="([^"]*)"', seqs[REF_SEQ_INDEX].group(0))
    if ref != SEQ_AFTER:
        raise StructuralSurprise(
            f"sequence {REF_SEQ_INDEX} is {ref}, expected {SEQ_AFTER} -- it is "
            f"the reference the rewrite is checked against")
    notes.append(f"sequence {SEQ_INDEX}: {' -> '.join(listed)}")
    notes.append(f"sequence {REF_SEQ_INDEX}: {' -> '.join(ref)}  (reference, "
                 f"left untouched)")

    stop = int(re.search(r'stopsurface="(\d+)"', seqs[SEQ_INDEX].group(0)).group(1))
    rstop = int(re.search(r'stopsurface="(\d+)"',
                          seqs[REF_SEQ_INDEX].group(0)).group(1))
    notes.append(f"  stopsurface={stop} -> {listed[stop]!r}   "
                 f"(sequence {REF_SEQ_INDEX}: {rstop} -> {ref[rstop]!r})")
    return notes


# --------------------------------------------------------------------------
# The build.
# --------------------------------------------------------------------------

def build(base_text: str) -> str:
    text = base_text

    # -- 1/2. the detector: height -> parameter, and a visible draw radius ----
    start, end = _detector_span(text)
    block = text[start:end]
    block, n = re.subn(r'(<variable name="z" value=")%s(")' % re.escape(BASE_Z),
                       r"\g<1>%s\g<2>" % PF_PARAM, block)
    if n != 1:
        raise StructuralSurprise(f"replaced {n} z literals, expected 1")
    block, n = re.subn(
        r'(<float_ap name="radius" value=")%s(")' % re.escape(BASE_FLOAT_AP),
        r"\g<1>%s\g<2>" % repr(DRAW_RADIUS), block)
    if n != 1:
        raise StructuralSurprise(f"replaced {n} float_ap radii, expected 1")
    text = text[:start] + block + text[end:]

    # -- 3. the single_param, appended to the multiconfig block ---------------
    # Last, after axi_angle, so the diff is an insertion at one point rather
    # than a reshuffle -- and so the file still reads in the order Quadoa wrote.
    close = text.index("\t</multiconfig>")
    text = (text[:close]
            + _SINGLE_PARAM.format(id=PF_PARAM, tcuid=PF_TCUID,
                                   value=repr(PF_HEIGHT))
            + text[close:])

    # -- 4. sequence 3: sun -> helio_surf -> prime_focus ----------------------
    seq = _sequences(text)[SEQ_INDEX]
    body = seq.group(0)

    # `receiver` goes; `secondary` becomes the detector in place.  Editing the
    # third entry rather than deleting both and writing a fresh one keeps its
    # `spuid` -- which is a slot id within the sequence, not a property of the
    # surface: `surfp_uid4` is the third entry's spuid in all four sequences,
    # whichever surface sits there.
    body, n = re.subn(r'[ \t]*<surf id="receiver"[^>]*/>\n', "", body)
    if n != 1:
        raise StructuralSurprise(f"removed {n} receiver entries, expected 1")
    body, n = re.subn(r'<surf id="secondary"(\s+surf="1")\s+action="reflection"',
                      r'<surf id="%s"\g<1> action="detector"' % DETECTOR, body)
    if n != 1:
        raise StructuralSurprise(f"retargeted {n} secondary entries, expected 1")

    # `stopsurface` is an index INTO THE SEQUENCE.  It was 2 -> `secondary`;
    # dropping a surface after it leaves it at 2 -> `prime_focus`, which is
    # exactly where sequence 0 points its own stop.  Verified by NAME, and
    # rewritten by name if the index ever needs to move, because the aperture
    # stop governs how the source samples the pupil: pointed at the wrong
    # surface the trace still runs and quietly loses part of its bundle.
    ids = re.findall(r'<surf id="([^"]*)"\s+surf=', body)
    head_end = body.index(">") + 1
    want = re.findall(r'<surf id="([^"]*)"',
                      _sequences(text)[REF_SEQ_INDEX].group(0))[
        int(re.search(r'stopsurface="(\d+)"',
                      _sequences(text)[REF_SEQ_INDEX].group(0)).group(1))]
    head, n = re.subn(r'stopsurface="\d+"', f'stopsurface="{ids.index(want)}"',
                      body[:head_end])
    if n != 1:
        raise StructuralSurprise("could not re-point stopsurface")
    body = head + body[head_end:]

    return text[:seq.start()] + body + text[seq.end():]


# --------------------------------------------------------------------------
# Post-flight.  Cheap enough to run on every build.
# --------------------------------------------------------------------------

def verify(base_text: str, text: str) -> list[tuple[bool, str]]:
    checks: list[tuple[bool, str]] = []

    def add(ok, msg):
        checks.append((bool(ok), msg))

    ok, fixed, msg = xml_check(text)
    _, base_fixed, _ = xml_check(base_text)
    add(ok, msg)
    if not ok:
        return checks
    add(fixed == base_fixed,
        f"no new malformed tags: {fixed} duplicate-attribute tags, base has "
        f"{base_fixed}")

    # 1. exactly one detector surface, still
    g0, g1 = _geometry(text)
    n_geo = len(re.findall(r'<surf id="%s" refe=' % DETECTOR, text[g0:g1]))
    n_ref = len(re.findall(r'<surf id="%s" surf=' % DETECTOR, text))
    add(n_geo == 1, f"exactly one {DETECTOR} geometry surface: {n_geo} "
                    f"({n_ref} sequence references: sequences 0 and {SEQ_INDEX})")

    # 2. it is parameterised, drawn, and unbounded
    start, end = _detector_span(text)
    block = text[start:end]
    z = re.search(r'<variable name="z" value="([^"]*)"', block).group(1)
    ap = re.search(r'<float_ap name="radius" value="([^"]*)" is_locked="([^"]*)"',
                   block)
    add(z == PF_PARAM, f"{DETECTOR} z = {z!r} (was {BASE_Z!r})")
    add(ap.group(1) == repr(DRAW_RADIUS),
        f"{DETECTOR} float_ap radius = {ap.group(1)} is_locked={ap.group(2)!r} "
        f"-- drawing only, does not clip")
    add("<aperture" not in block,
        f"{DETECTOR} has no <aperture>: unbounded, spillage stays a post-step")

    # 3. multiconfig grew by exactly the one parameter
    have_p = re.findall(r'<param id="([^"]*)"', text)
    have_s = re.findall(r'<single_param id="([^"]*)"', text)
    cols = re.search(r'<multiconfig[^>]*columns="(\d+)"', text).group(1)
    base_cols = re.search(r'<multiconfig[^>]*columns="(\d+)"', base_text).group(1)
    add(have_p == PARAMS and cols == base_cols,
        f"per-config params unchanged: {have_p}, columns={cols}")
    add(have_s == SINGLE_PARAMS + [PF_PARAM],
        f"single_params: {have_s}")
    val = re.search(r'<single_param id="%s"[^>]*>\s*<variable name="val" '
                    r'value="([^"]*)"' % PF_PARAM, text)
    add(val is not None and float(val.group(1)) == PF_HEIGHT,
        f"{PF_PARAM} = {val and val.group(1)} mm "
        f"(= 27000 vertex + 20000 throw, mirroring rec_offset = -20000)")
    tc = re.findall(r'tcuid="([^"]*)"', text)
    add(len(tc) == len(set(tc)), f"tcuids unique: {len(tc)} entries")

    # 4. the referencing works at all: the name appears as a value exactly once
    refs = re.findall(r'<variable name="[a-z]+" value="%s"' % PF_PARAM, text)
    add(len(refs) == 1, f"{PF_PARAM} referenced by exactly one geometry "
                        f"variable: {len(refs)}")

    # 5. sequence 3
    seqs = _sequences(text)
    body = seqs[SEQ_INDEX].group(0)
    listed = re.findall(r'<surf id="([^"]*)"', body)
    add(listed == SEQ_AFTER, f"sequence {SEQ_INDEX}: {' -> '.join(listed)}")
    stop = int(re.search(r'stopsurface="(\d+)"', body).group(1))
    at = listed[stop] if stop < len(listed) else "(out of range)"
    add(at == DETECTOR, f"sequence {SEQ_INDEX} stopsurface={stop} -> {at!r}")
    actions = dict(re.findall(r'<surf id="([^"]*)"[^>]*action="([^"]*)"', body))
    add(actions == {"sun": "source", "helio_surf": "reflection",
                    DETECTOR: "detector"}, f"actions: {actions}")

    # 6. the strongest one: the rewritten sequence 3 IS sequence 0, modulo the
    #    three attributes that must differ (name, uid, visibility).
    ref = seqs[REF_SEQ_INDEX].group(0)
    add(_normalise_seq(body) == ref,
        "sequence 3 is byte-identical to sequence 0 apart from name/sequid/"
        "is_visible -- same surfaces, same stop, same 3500 mm 38484.5 W source")

    # 7. nothing else moved
    a, b = _sequences(base_text), seqs
    same = [i for i in range(len(a)) if i != SEQ_INDEX and a[i].group(0) == b[i].group(0)]
    add(len(same) == len(a) - 1, f"other sequences untouched: {same}")

    # 8. the two structural checks build_occluder_field_model.py runs, applied
    #    here for parity: this builder adds no geometry, so both must come out
    #    *unchanged* from the base rather than merely self-consistent.
    def n(pat, t):
        return len(re.findall(pat, t))

    balance = [(n(r"<pos\b", text), n(r"</pos>", text)),
               (n(r"<surf\b", text[g0:g1]), n(r"</surf>", text[g0:g1]))]
    same_counts = all(n(p, text) == n(p, base_text)
                      for p in (r"<pos\b", r"</pos>", r"</surf>", r"<obscuration\b"))
    add(balance[0][0] == balance[0][1] and balance[1][0] == balance[1][1]
        and same_counts,
        f"tag balance unchanged: <pos> {balance[0][0]}/{balance[0][1]}, "
        f"geometry <surf> {balance[1][0]}/{balance[1][1]}, all counts equal to "
        f"the base")

    dup = []
    for attr in ("puuid", "suid"):
        v = re.findall(attr + r'="([^"]*)"', text)
        dup += [f"{attr}={x}" for x in {y for y in v if v.count(y) > 1}]
    for i, m in enumerate(seqs):
        v = re.findall(r'spuid="([^"]*)"', m.group(0))
        dup += [f"seq{i} spuid={x}" for x in {y for y in v if v.count(y) > 1}]
    add(not dup, f"unique ids: {'OK' if not dup else dup}")

    # 9. the diff is exactly the intended shape and nothing more
    import difflib
    d = [ln for ln in difflib.unified_diff(base_text.splitlines(),
                                           text.splitlines(), n=0, lineterm="")
         if ln[:1] in "+-" and ln[:3] not in ("+++", "---")]
    added = [ln for ln in d if ln.startswith("+")]
    removed = [ln for ln in d if ln.startswith("-")]
    add(len(added) == 6 and len(removed) == 4,
        f"line diff vs base: {len(added)} added / {len(removed)} removed "
        f"(want 6/4: z, float_ap and the seq-3 surf line each rewritten = 3 "
        f"pairs, +3 single_param lines, -1 receiver entry)")
    touched = sum(w in ln for ln in d
                  for w in ("<variable name=\"z\"", "float_ap", "single_param",
                            "<variable name=\"val\"", "<surf id="))
    add(touched == len(d),
        f"every changed line is a z / float_ap / single_param / sequence-surf "
        f"line: {touched}/{len(d)}")
    return checks


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--base", type=Path, default=BASE)
    p.add_argument("--out", type=Path, default=OUT)
    p.add_argument("--force", action="store_true",
                   help="overwrite an existing output file")
    a = p.parse_args(argv)

    previous = a.out.read_text(encoding="utf-8") if a.out.exists() else None
    if previous is not None and not a.force:
        print(f"refusing to overwrite {a.out} -- pass --force")
        return 2

    base_text = a.base.read_text(encoding="utf-8")
    print(f"base: {a.base}  ({len(base_text):,} bytes)")
    for note in check_base(base_text):
        print(f"  {note}")
    print()

    text = build(base_text)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    eol = newline_of(a.base)
    with open(a.out, "w", encoding="utf-8", newline=eol) as fh:
        fh.write(text)

    on_disk = a.out.stat().st_size
    label = "CRLF" if len(eol) == 2 else "LF"
    print(f"wrote {a.out}  ({on_disk:,} bytes on disk, "
          f"{on_disk - a.base.stat().st_size:+,} over the base; "
          f"line endings {label}, as the base)")
    print(f"  {DETECTOR} z: literal {BASE_Z} -> single_param {PF_PARAM} "
          f"= {PF_HEIGHT:.0f} mm")
    print(f"  {DETECTOR} float_ap radius: {BASE_FLOAT_AP} -> {DRAW_RADIUS:.0f} mm "
          f"(drawing only)")
    print(f"  sequence {SEQ_INDEX}: {' -> '.join(SEQ_BEFORE)}")
    print(f"           -> {' -> '.join(SEQ_AFTER)}")
    if previous is not None:
        print(f"  idempotent: {'identical to' if previous == text else 'DIFFERS from'}"
              f" the previous {a.out.name}")
    print()

    bad = 0
    for ok, msg in verify(base_text, text):
        print(f"  [{'ok' if ok else 'FAIL'}] {msg}")
        bad += not ok
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
