"""Write an .optx with one heliostat plus its shading *and* blocking occluders.

The analytic shading/blocking never entered the ray trace, so it has never been
checked against Quadoa on the same geometry. This builds the model that checks
it: one heliostat, the neighbours that occlude it, and a sequence that meets them
on both legs of the path.

    sun -> [shading occluders] -> helio_surf -> [blocking occluders] -> secondary
        -> receiver

The same physical mirror appears twice, once per leg, as two surfaces sharing a
position and orientation. That is not a duplication error -- a sequential trace
meets each listed surface once, and the beam really does pass that neighbour
twice: on the way in as sunlight, on the way out as the reflected beam. Giving
each leg its own surface keeps the sequence unambiguous rather than relying on
Quadoa's behaviour when an id is referenced twice.

Built by text surgery on a template ``.optx`` rather than through the API,
because the API has no way to create geometry -- its ``add*`` methods attach
features to surfaces that already exist.

Occluders are ``action="transmission"`` surfaces carrying a rectangular
``<obscuration>``: transparent everywhere except the mirror's own footprint, so
they stop rays without reflecting any.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

# Per-configuration multiconfig parameters that describe one heliostat.
HELIOSTAT_PARAMS = ("rot_az", "rot_el", "posx", "posy", "c3", "c4", "c5")


@dataclass
class BuildReport:
    path: Path
    heliostat_id: int
    timestep: str
    shading: list
    blocking: list
    eta_shade: float
    eta_block: float
    eta_secondary: float

    def describe(self) -> str:
        L = [
            f"wrote {self.path}",
            f"  heliostat {self.heliostat_id} at {self.timestep}",
            "",
            "  sequence 4 (index 3):  sun -> "
            + " -> ".join(f"s{o.heliostat_id}" for o in self.shading)
            + " -> helio_surf -> "
            + " -> ".join(f"b{o.heliostat_id}" for o in self.blocking)
            + " -> secondary -> receiver",
            "",
            "  beamdown predicts, for a trace of this model:",
            f"    rays reaching helio_surf, vs the same model without the shading"
            f" obscurations   {self.eta_shade:.4f}",
            f"    rays reaching secondary, vs the same without the blocking"
            f" obscurations     {self.eta_block:.4f}",
        ]
        if self.eta_secondary < 0.999:
            L.append(f"    (the axicon itself also shades this one: "
                     f"eta_secondary {self.eta_secondary:.4f})")
        # Two callers with two shapes of occluder: the bespoke builder carries
        # both effects on one record, the slot planner ranks each leg separately.
        L += ["", "  occluders:",
              f"    {'leg':>6s} {'id':>5s} {'x (mm)':>12s} {'y (mm)':>12s} "
              f"{'rot_az':>9s} {'rot_el':>8s} {'effect':>8s}"]
        for leg, group in (("shade", self.shading), ("block", self.blocking)):
            for o in group:
                x = getattr(o, "posx_mm", None)
                x = o.x_mm if x is None else x
                y = getattr(o, "posy_mm", None)
                y = o.y_mm if y is None else y
                pct = getattr(o, "fraction", None)
                pct = (o.shades_pct if leg == "shade" else o.blocks_pct) \
                    if pct is None else pct * 100.0
                L.append(f"    {leg:>6s} {o.heliostat_id:5d} {x:12.1f} {y:12.1f} "
                         f"{o.rot_az_deg:9.3f} {o.rot_el_deg:8.3f} {pct:7.2f}%")
        return "\n".join(L)


def _span_of_pos(text: str, start: int) -> tuple[int, int]:
    """Character span of the ``<pos>`` element beginning at ``start``.

    A depth counter over pos tags, because the blocks nest two deep and a
    regex for the closing tag would stop at the inner one.
    """
    depth = 0
    for m in re.finditer(r"<pos\b|</pos>", text[start:]):
        depth += 1 if m.group(0) == "<pos" else -1
        if depth == 0:
            return start, start + m.end()
    raise ValueError("unterminated <pos> element")


def _set_value(text: str, tag: str, param_id: str, value: float) -> str:
    """Set ``val_0`` of a multiconfig entry.

    ``<param>`` holds one value per configuration, in ``val_0``, ``val_1``, ...;
    ``<single_param>`` holds one value for the whole model, in a child named just
    ``val``. Same edit either way once the right child name is used -- and the
    difference is the reason a config index means nothing on a single_param.
    """
    m = re.search(r'<%s id="%s"[^>]*>' % (tag, re.escape(param_id)), text)
    if not m:
        raise KeyError(f"no <{tag} id={param_id!r}> in the template")
    child = "val_0" if tag == "param" else "val"
    n = re.compile(r'(<variable name="%s" value=")([^"]*)(")' % child).search(
        text, m.end())
    if not n:
        raise KeyError(f"{param_id} has no {child}")
    return text[:n.start()] + f"{n.group(1)}{value!r}{n.group(3)}" + text[n.end():]


def _set_param(text: str, param_id: str, value: float) -> str:
    return _set_value(text, "param", param_id, value)


def _set_single(text: str, param_id: str, value: float) -> str:
    return _set_value(text, "single_param", param_id, value)


def build_from_slot_model(cfg, summary, heliostat_id: int, timestep: str,
                          out_path: Path | None = None,
                          template: Path | None = None) -> BuildReport:
    """Fill the sweep's occluder-slot model in for one heliostat, to look at.

    Unlike :func:`build`, which clones a bespoke assembly per neighbour, this
    writes the slot values into ``heliostat_field_occluders.optx`` -- the model
    the sweep itself traces. So what you open is exactly the geometry the sweep
    saw for that heliostat, not a reconstruction of it. Unused slots keep their
    parked position a thousand kilometres out and simply do not appear.

    Pure text editing: no Quadoa session, no licence.
    """
    from . import field as field_mod
    from . import occluder_slots as slots_mod
    from . import shading as shading_mod
    from .secondary import get_strategy

    template = Path(template or (cfg.repo_root / "models"
                                 / "heliostat_field_occluders.optx"))
    if not template.exists():
        raise FileNotFoundError(
            f"{template} does not exist. Build it with "
            f"scripts/build_occluder_field_model.py first.")
    text = template.read_text(encoding="utf-8")

    rows = summary[summary.timestep == timestep]
    if not len(rows):
        raise KeyError(f"timestep {timestep!r} is not in this run")
    az = float(rows.solar_az_deg.iloc[0])
    el = float(rows.solar_el_deg.iloc[0])

    fld = field_mod.load_field(cfg)
    strategy = get_strategy(cfg)
    sols = [strategy.solve(float(fld.x_mm[i]), float(fld.y_mm[i]), az, el, cfg.geometry)
            for i in range(len(fld))]
    geoms, aims = shading_mod.build_geometries(fld, sols, cfg)
    neighbours = field_mod.neighbour_pairs(
        fld, shading_mod.search_radius_for(el, cfg.field.mirror_height_mm,
                                           cfg.field.mirror_width_mm))
    plans = slots_mod.plan_field(geoms, aims, fld.ids, neighbours,
                                 shading_mod.sun_vector(az, el))
    i = int(next(k for k, h in enumerate(fld.ids) if int(h) == heliostat_id))
    plan = plans[i]
    sol = sols[i]

    values = {"rot_az": sol.rot_az_deg, "rot_el": sol.rot_el_deg,
              "posx": float(fld.x_mm[i]), "posy": float(fld.y_mm[i]),
              "c3": sol.c3, "c4": sol.c4, "c5": sol.c5}
    for name in HELIOSTAT_PARAMS:
        text = _set_param(text, name, values[name])
    text = _set_single(text, "solaz", az)
    text = _set_single(text, "solze", 90.0 - el)

    for prefix, slots, count in (("so", plan.shading, slots_mod.N_SHADE),
                                 ("bo", plan.blocking, slots_mod.N_BLOCK)):
        for k in range(count):
            name = f"{prefix}{k}"
            if k < len(slots):
                s = slots[k]
                text = _set_single(text, f"{name}_x", s.x_mm)
                text = _set_single(text, f"{name}_y", s.y_mm)
                text = _set_single(text, f"{name}_az", s.rot_az_deg)
                text = _set_single(text, f"{name}_el", s.rot_el_deg)
            else:
                text = _set_single(text, f"{name}_x", slots_mod.PARK_MM)
                text = _set_single(text, f"{name}_y", slots_mod.PARK_MM)

    if out_path is None:
        out_path = cfg.path(f"models/inspect_h{heliostat_id}_{timestep}_slots.optx")
    out_path = Path(out_path)
    out_path.write_text(text, encoding="utf-8")

    step = rows.set_index("heliostat_id")
    r = step.loc[heliostat_id] if heliostat_id in step.index else None
    return BuildReport(
        path=out_path, heliostat_id=heliostat_id, timestep=timestep,
        shading=plan.shading, blocking=plan.blocking,
        eta_shade=float(r.eta_shade) if r is not None else float("nan"),
        eta_block=float(r.eta_block) if r is not None else float("nan"),
        eta_secondary=float(r.get("eta_secondary", 1.0)) if r is not None else 1.0,
    )


def build(cfg, summary, heliostat_id: int, timestep: str, template: Path,
          out_path: Path | None = None) -> BuildReport:
    """Write the model. Pure text editing -- no Quadoa session, no licence."""
    from . import field as field_mod
    from . import occluders as occ_mod
    from .secondary import get_strategy

    template = Path(template)
    text = template.read_text(encoding="utf-8")

    found, totals = occ_mod.occluders_for(cfg, summary, heliostat_id, timestep)
    shading = [o for o in found if o.shades_pct > 0.05]
    blocking = [o for o in found if o.blocks_pct > 0.05]
    if not shading and not blocking:
        raise ValueError(f"heliostat {heliostat_id} at {timestep} is not occluded")

    # -- the heliostat itself ------------------------------------------------
    fld = field_mod.load_field(cfg)
    row = int(next(k for k, i in enumerate(fld.ids) if int(i) == heliostat_id))
    sol = get_strategy(cfg).solve(
        float(fld.x_mm[row]), float(fld.y_mm[row]),
        totals["solar_az_deg"], totals["solar_el_deg"], cfg.geometry)

    values = {"rot_az": sol.rot_az_deg, "rot_el": sol.rot_el_deg,
              "posx": float(fld.x_mm[row]), "posy": float(fld.y_mm[row]),
              "c3": sol.c3, "c4": sol.c4, "c5": sol.c5}
    for name in HELIOSTAT_PARAMS:
        text = _set_param(text, name, values[name])
    # The model takes a zenith angle, not an elevation.
    text = _set_single(text, "solaz", totals["solar_az_deg"])
    text = _set_single(text, "solze", 90.0 - totals["solar_el_deg"])

    # -- clone the occluder assemblies ---------------------------------------
    first = text.index('<pos id="p327"')
    last_start = text.index('<pos id="p399"')
    _, last_end = _span_of_pos(text, last_start)
    tpl_start, tpl_end = _span_of_pos(text, first)
    block_tpl = text[tpl_start:tpl_end]

    def clone(o, prefix: str, n: int) -> str:
        b = block_tpl
        b = b.replace('id="p327"', f'id="p{prefix}{o.heliostat_id}"')
        b = b.replace('id="h327"', f'id="h{prefix}{o.heliostat_id}"')
        b = b.replace('id="s327"', f'id="{prefix}{o.heliostat_id}"')
        b = b.replace('puuid="assembly6"', f'puuid="assembly{100 + 3 * n}"')
        b = b.replace('puuid="assembly7"', f'puuid="assembly{101 + 3 * n}"')
        b = b.replace('puuid="surf6"', f'puuid="surf{100 + n}"')
        b = b.replace('suid="surf_uid6"', f'suid="surf_uid{100 + n}"')
        # Position: the outer <pos> carries x/y, the inner one the rotations.
        head, tail = b.split('<pos id="h', 1)
        head = re.sub(r'(name="x" value=")[^"]*(")', rf'\g<1>{o.posx_mm!r}\g<2>', head)
        head = re.sub(r'(name="y" value=")[^"]*(")', rf'\g<1>{o.posy_mm!r}\g<2>', head)
        tail = re.sub(r'(name="ry" value=")[^"]*(")',
                      rf'\g<1>{-o.rot_el_deg!r}\g<2>', tail, count=1)
        tail = re.sub(r'(name="rz" value=")[^"]*(")',
                      rf'\g<1>{o.rot_az_deg!r}\g<2>', tail, count=1)
        return head + '<pos id="h' + tail

    blocks = [clone(o, "s", n) for n, o in enumerate(shading)]
    blocks += [clone(o, "b", n + len(shading)) for n, o in enumerate(blocking)]
    text = text[:tpl_start] + "\n\t\t\t".join(blocks) + text[last_end:]

    # -- rewire sequence 4 ---------------------------------------------------
    seqs = list(re.finditer(r"<sequence\b.*?</sequence>", text, re.S))
    seq = seqs[3]
    body = seq.group(0)
    ref = re.search(r'<surf id="s327"[^>]*/>', body)
    if ref is None:
        raise ValueError("sequence 4 does not reference the template occluder")
    ref_tpl = ref.group(0)

    def seq_ref(o, prefix, n):
        r = ref_tpl.replace('id="s327"', f'id="{prefix}{o.heliostat_id}"')
        return re.sub(r'spuid="surfp_uid\d+"', f'spuid="surfp_uid{100 + n}"', r)

    shade_refs = "\n\t\t\t\t".join(seq_ref(o, "s", n) for n, o in enumerate(shading))
    block_refs = "\n\t\t\t\t".join(
        seq_ref(o, "b", n + len(shading)) for n, o in enumerate(blocking))

    # Drop the template's three occluder references, then insert both groups
    # around helio_surf.
    new_body = re.sub(r'\s*<surf id="s(327|328|399)"[^>]*/>', "", body)
    helio = re.search(r'\s*<surf id="helio_surf"[^>]*/>', new_body)
    new_body = (new_body[:helio.start()]
                + "\n\t\t\t\t" + shade_refs
                + helio.group(0)
                + "\n\t\t\t\t" + block_refs
                + new_body[helio.end():])
    text = text[:seq.start()] + new_body + text[seq.end():]

    if out_path is None:
        out_path = cfg.path(f"models/inspect_h{heliostat_id}_{timestep}"
                            f"_shade_block.optx")
    out_path = Path(out_path)
    out_path.write_text(text, encoding="utf-8")

    return BuildReport(path=out_path, heliostat_id=heliostat_id, timestep=timestep,
                       shading=shading, blocking=blocking,
                       eta_shade=totals["eta_shade"], eta_block=totals["eta_block"],
                       eta_secondary=totals["eta_secondary"])
