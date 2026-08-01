"""Edit the Quadoa .optx multiconfig block.

The shipped ``heliostat_field_model_mcfg.optx`` has 24 configuration columns;
an export may want more (one per heliostat being inspected).

Why edit the file rather than call the API: checked against the complete
QuadoaCore surface (256 public methods, dumped by ``tools/dump_api.py``). The
only configuration-related methods are

    getNrConfigs  getCurrentConfig  setConfig  nextConfig  previousConfig
    get/setMulticonfParam           get/setCurrentMulticonfParam

-- all of which *select* a configuration or read/write parameters within an
existing one. The six ``add*`` methods all attach features to surfaces
(aperture, coating, form, obscuration, phase, polarization). Nothing creates a
configuration, under that name or under Quadoa's own term "realisation". So the
column count can only be changed in the file.

This is deliberately **text surgery, not an XML round-trip**: re-serialising the
whole document through ElementTree would reorder attributes and rewrite
whitespace across a file Quadoa wrote, risking changes far outside the part we
care about. Here only the multiconfig header and the per-param variable lists are
touched, and the result is verified by loading it back and reading the new column.

The bulk 645-heliostat sweep does not need any of this -- it reuses a single
configuration. The live caller is :func:`beamdown.inspect_model.
export_for_inspection`, which grows a copy of the model so each inspected
heliostat gets its own configuration.

The 25-configuration *figure* model is NOT built here. This module once carried
a ``build_figure_model`` that did the column surgery in Python and then needed a
licence seat to write the per-config values through ``setMulticonfParam``; the
second half never ran anywhere, which is why the shipped
``models/figure_model_25cfg.optx`` has configs 1-24 all zeros. The working,
licence-free generator is ``scripts/build_figure_model.py`` (with
``scripts/verify_figure_model.py`` for the seat-gated half).
"""

from __future__ import annotations

import re
import shutil
from pathlib import Path

_MULTICONF_RE = re.compile(r'(<multiconfig\s+columns=")(\d+)(")')
_PARAM_BLOCK_RE = re.compile(
    r'(<param\s+id="(?P<id>[^"]+)"[^>]*>)(?P<body>.*?)(</param>)', re.DOTALL
)
_VAR_RE = re.compile(r'[ \t]*<variable\s+name="val_(\d+)"[^/]*/>\s*\n')


def current_columns(text: str) -> int:
    m = _MULTICONF_RE.search(text)
    if not m:
        raise ValueError("No <multiconfig columns=...> block found")
    return int(m.group(2))


def expand_multiconfig(src: Path, dst: Path, n_columns: int, backup: bool = True) -> dict:
    """Grow the multiconfig block to ``n_columns``. Returns a change report."""
    src, dst = Path(src), Path(dst)
    text = src.read_text(encoding="utf-8")
    have = current_columns(text)

    if n_columns < have:
        raise ValueError(
            f"Refusing to shrink the multiconfig block ({have} -> {n_columns}); "
            f"that would discard stored per-config values."
        )
    if n_columns == have:
        if src != dst:
            shutil.copyfile(src, dst)
        return {"columns_before": have, "columns_after": have, "params_extended": 0}

    if backup and dst.exists():
        shutil.copyfile(dst, dst.with_suffix(dst.suffix + ".bak"))

    extended = []

    def grow(match: re.Match) -> str:
        head, body, tail = match.group(1), match.group("body"), match.group(4)
        vars_found = _VAR_RE.findall(body)
        if not vars_found:
            return match.group(0)
        highest = max(int(v) for v in vars_found)
        if highest >= n_columns - 1:
            return match.group(0)

        # Clone the last <variable> line so indentation and attributes match.
        last = None
        for m in _VAR_RE.finditer(body):
            last = m
        template = last.group(0)
        additions = []
        for i in range(highest + 1, n_columns):
            line = re.sub(r'name="val_\d+"', f'name="val_{i}"', template)
            line = re.sub(r'value="[^"]*"', 'value="0.0"', line, count=1)
            additions.append(line)

        extended.append(match.group("id"))
        return head + body.rstrip("\n \t") + "\n" + "".join(additions) + tail

    new_text = _PARAM_BLOCK_RE.sub(grow, text)
    new_text = _MULTICONF_RE.sub(lambda m: m.group(1) + str(n_columns) + m.group(3),
                                 new_text, count=1)

    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(new_text, encoding="utf-8")
    return {
        "columns_before": have,
        "columns_after": n_columns,
        "params_extended": len(extended),
        "params": extended,
        "path": str(dst),
    }
