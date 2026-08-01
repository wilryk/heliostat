# Beam-down heliostat field analysis

Traces each heliostat of a beam-down solar concentrator **individually** through a
Quadoa optical model, stores the receiver-plane rays, and answers per-heliostat
questions across time of day and time of year:

- What does one heliostat produce? What do all of them produce together?
- Which heliostat is worst, at a given hour or over a whole year?
- How much does mutual shading and blocking cost?
- How much energy does the field collect in a year?

## Layout

```
config.toml       every tunable constant -- site, geometry, ray budget, dates, storage
beamdown/         the package
data/             heliostat field positions and downselects
models/           Quadoa .optx/.uix models
tests/            verify.py, smoke_sweep.py, test_axicon_parity.py, test_gui.py
tools/            dump_api.py (QuadoaCore API reference)
legacy/           the original scripts, kept for reference and parity testing
analysis_output/  sweep results (generated)
```

## Quick start

```bash
python tests/verify.py --no-quadoa     # 11 stages, no license needed
python tests/test_gui.py               # GUI views, table regression, Trace tab; no license
python -m beamdown info                # time grid, trace count, cost estimate
python -m beamdown sweep --output analysis_output/demo25
python -m beamdown figures --output analysis_output/demo25
python -m beamdown gui --output analysis_output/demo25       # interactive
python -m beamdown compare RUN_A RUN_B --attribute           # is one better?
```

## Comparing two sweeps

`compare` answers "did that change help?" over two runs of the same field. The
headline is **power delivered inside an aperture of radius R, swept over R** --
one radius cannot settle it, because a change that tightens the tail while
shaving the peak wins at small R and loses at large R. If one run dominates at
every radius that is a real result; if the curves cross, the crossing radius is
the number worth quoting.

`--attribute` bands the per-heliostat change by `1/L^2`, the size of the axicon
shape correction. An improvement that really comes from that correction has to
concentrate where `1/L^2` is largest -- a flat trend means something else moved.
Everything is post-processing on stored counts: seconds, no license.

## The GUI

`python -m beamdown gui` opens a desktop explorer over a sweep (Tkinter +
matplotlib, both already present -- nothing to install). One heliostat and one
timestep are selected at a time and every view follows that selection:

| tab | shows |
|---|---|
| Field | the field from directly above, each heliostat its actual rectangle, coloured by any summary column, with the shading and blocking losses washed over it; click to select |
| Spot | the selected heliostat's receiver flux beside the whole-field flux, as an image or as encircled energy |
| Through day | the selected heliostat's metric per hour, over each traced date, above the field total |
| Distribution | that metric across the field, with the selection's rank marked |
| Table | every row for this timestep, sortable; click to select |
| Trace | set up a new sweep: every option, the exact command it will run, and Run/Stop with a live log tail |

Nothing re-traces. DNI, shading/blocking and the aperture are applied to stored
counts, so they respond instantly. **Open in Quadoa** exports an `.optx` with the
selected heliostat's pointing, shape and sun already loaded -- it runs on a
background thread and reports plainly when no license seat is free, which is the
normal case while a sweep is running.

**The aperture is a radius, and it is physical, not a zoom.** Setting it
recomputes power inside it, spillage, and the `power_in_aperture_w` column the
field map can be coloured by, using the same circular mask `beamdown compare`
sweeps. Leave it blank for the whole receiver window.

**The shading/blocking checkbox reaches every view.** The sweep folds shading x
blocking and a 1000 W/m² DNI into `power_w` and `peak_flux_w_m2` before writing
them, so those columns are weighted whatever the checkbox says; the GUI divides
the weights back out when it is off. Without that, a single heliostat's spot
looked unaffected by shading -- a scalar scales every pixel *and* the autoscaled
colour bar, so only the printed numbers can show it -- while the field spot,
being a sum of 645 different weights, visibly changed shape.

**Encircled energy is measured about the receiver axis, not the spot centroid**,
because that is what an aperture intercepts; the summary's `r50_mm`/`r90_mm` are
centroid-referenced and will differ for a spot that has walked off axis.

**The field map is a top-down projection**, so pointing is legible: at low sun the
mirrors stand up and draw as slivers, at noon as near-full rectangles. Two
overlays sit on it, and they are different things:

- **Blue on a mirror** is that mirror's shading loss, and **red** its blocking
  loss -- the numbers the sweep computed, drawn where they belong.
- **The pale wash is the shadow pattern on the ground**, drawn from a 5 m
  `draw_pedestal_height_mm` so shadows fall beside the mirrors rather than
  straddling them. It is thrown clear of the mirror that cast it -- 29 m at 9.7
  deg -- so do not read it as "this mirror is shaded": the mirror actually shaded
  sits at pedestal height, where the offset is zero, so it is the nearest up-sun
  neighbour. Lifting the field for the picture is exact rather than cosmetic,
  because mutual shading is invariant to a height every mirror shares.

The metric ramp is a truncated grey deliberately: a translucent overlay is
invisible over viridis's dark purple and unreadable over its yellow.

**Spot bins** re-histograms the receiver map. Coarsening block-sums the stored
bins (exact, instant, no raw rays needed); refining goes back to the rays, which
costs a whole-file read for the field and is reported in the status bar.
Delivered power is bin-independent and stays fixed to the last digit across every
setting; peak flux is not, and is quoted with its bin size.

A sweep still in progress can be opened read-only; **Reload** picks up timesteps
as they land. `beamdown/explore.py` is the older single-window matplotlib
version, superseded by this and kept only as a fallback.

Run the full field over declination-spaced dates:

```bash
python -m beamdown sweep --all-heliostats --suggest-dates 8 --output analysis_output/full
```

Sweeps are resumable -- rerun the same command and completed timesteps are skipped.

## Setting up a run

Every option is a `beamdown sweep` flag, and none of them requires editing
config.toml. That is deliberate: a sweep already running re-reads config.toml
when it writes its end-of-run report, so an edit made mid-run corrupts that
report, and two people cannot set up two different runs at once.

```bash
python -m beamdown sweep --help
```

| flag | overrides | notes |
|---|---|---|
| `--output DIR` | `storage.root` | where the run lands |
| `--dates D…` / `--suggest-dates N` | `sweep.dates` | `--dates` wins if both are given |
| `--hour-step H` | `sweep.hour_step` | MAXIMUM sample spacing, not a clock grid |
| `--sunrise-margin-min M` | `sweep.sunrise_margin_min` | |
| `--all-heliostats` | -- | all 645 instead of the downselect |
| `--rays N` | `trace.rays_per_heliostat` | total ray budget per heliostat; without `--rays-per-trace` the chunk size is clamped to it |
| `--rays-per-trace N` | `trace.rays_per_trace` | rays per `traceRays` call, so each heliostat costs `ceil(rays / N)` calls |
| `--workers N` | `trace.n_workers` | 1-4; **leave it at 1**, see below |
| `--no-resume` | -- | without it, finished timesteps are skipped |
| `--occluders` | `trace.model_file` | traces occlusion as geometry; picks its own model |
| `--secondary NAME` | `optics.secondary` | `axicon` / `prime_focus` / `cassegrain` |
| `--focus-height-mm Z` | `geometry.focus_height_mm` | required by prime_focus and cassegrain |
| `--rim-height-mm Z` | `geometry.secondary_rim_height_mm` | required by cassegrain |
| `--n-mirrors {1,2}` | `optics.n_mirrors` | 1 for prime_focus, 2 for the others |
| `--model-file PATH` | `trace.model_file` | different layouts need different `.optx` |

So a Cassegrain run needs no config edit at all:

```bash
python -m beamdown sweep --secondary cassegrain --focus-height-mm 24000 \
    --rim-height-mm 20000 --n-mirrors 2 --model-file models/cassegrain.optx \
    --all-heliostats --rays 60000 --workers 1 --output analysis_output/cass1
```

**The overrides are re-applied inside every worker.** A worker process calls
`load_config` itself, from disk, so a value set only on the driver's copy reaches
the report but not the trace. The whole set therefore travels to the workers as
data (`config.apply_overrides`, `sweep._init_worker`) and is recorded in the run's
`manifest.json` under `overrides`, so a stored run can say what it was asked to
do rather than only what config.toml said at the time.

### Ray budget and iteration count

`--rays` sets *how many* rays a heliostat gets; `--rays-per-trace` sets how many
`traceRays` calls it takes to emit them. They are different costs and the second
one used to be reachable only by editing config.toml.

```bash
python -m beamdown sweep --rays 60000 --rays-per-trace 12000 …   # 5 calls each
python -m beamdown sweep --rays 60000 --rays-per-trace 60000 …   # 1 call each
```

The interaction rule:

* `--rays` **alone** keeps the historical clamp — the chunk stays config.toml's
  value capped at the new budget, so a 1,000-ray smoke run does not claim
  60,000-ray chunks.
* **Both** are honoured literally, and `rays_per_trace > rays` is a startup
  error rather than a silent clamp, because the two flags then disagree about
  what was asked for (`config.validate_trace`).
* The chunks always sum **exactly** to `rays_per_heliostat`; a budget that does
  not divide evenly gets a short final chunk
  (`chunk_plan(100000, 30000) == [30000, 30000, 30000, 10000]`). One splitter,
  `config.chunk_plan`, is used by the trace, by `beamdown info` and by the GUI's
  derived call-count label, so none of them can describe a split that will not
  happen.

**How much a ray costs is not known.** Until the override fix below landed,
`--rays` never reached the sweep workers, so *every* stored run in
`analysis_output/` traced config.toml's 120,000 rays whatever it asked for —
full8 included, which asked for 60,000 and whose `manifest.json` still records
60,000 while the trace emitted 120,000. (Watts per ray is
`source_power_w / rays_emitted`, so anything scaling full8's stored counts by its
manifest ray budget reads **2× high** until that field is patched to `120000`;
ray counts and ratios such as `transmission` and `eta_*` are unaffected.)

Consequently there is no measurement anywhere here of how trace time splits
between the fixed per-heliostat cost, the fixed per-`traceRays`-call cost and the
marginal per-ray cost — and therefore no basis for asserting that rays are
expensive *or* that they are cheap. `scripts/probe_ray_cost.py` times all three
separately and prints the ray budget below which halving the rays saves less than
10%. It must not be run while a sweep holds the licence seat, and refuses to
start while any `analysis_output/.*.lock` exists.

The *precision* side is established: the per-heliostat Monte-Carlo noise floor is
0.0067 at 120,000 rays, and field power sums 645 independent traces, so
field-level relative noise is `0.0067/sqrt(645)` = 0.03% at 120,000 rays and
0.04% at 60,000 — both well under the 0.46% energy-integration residual.

**Layout validation re-runs after the overrides land**, so `--secondary
prime_focus` with no focus height fails immediately with the config module's own
message instead of mis-aiming 645 heliostats, and a self-consistent
`--secondary prime_focus --n-mirrors 1` does not warn about whatever config.toml
happens to say.

### The Trace tab

The GUI's **Trace** tab is the same thing with the flags spelled out as widgets:
dates (multi-select, plus free text), time grid, heliostat set, ray budget, rays
per `traceRays` call with a live read-only **→ N traceRays calls per heliostat**
readout beside it, workers, secondary layout with the height entries it needs,
occluders, model file, output directory and resume. It shows the exact
`python -m beamdown sweep …` line it would run, updated on every keystroke, with
a **Copy** button -- so it can be used purely to compose a command for a shell.

**Run** launches that command and refuses in three cases, before taking a licence
seat: any `analysis_output/.*.lock` directory exists (it names the lock and the
pid inside it), the output directory already exists and *resume* is not ticked, or
the form is incomplete. It takes the same `.<name>.lock` the run scripts take,
writes stdout and stderr to `analysis_output/<name>.log` as they do, and the
process is **detached** -- closing the window leaves the sweep running and its
lock in place, which the tab says on its face. **Stop** kills the process tree;
timesteps already written survive, so relaunching with *resume* continues from the
last completed one.

The **workers** spinbox is capped at 1 until explicitly unlocked. Measured:
asking for 2 got `asked for 2, 1 got a license seat`, and the failed second
request then leaked the first, so the next single-worker run found zero seats and
waited ten minutes -- with a modal H0038 dialog to dismiss. Quadoa's tracer
already runs at ~4x parallelism inside one session, so there is little to gain.

## How it works

**Counts are stored, never scaled flux.** The store holds raw receiver x/y
(int16-quantised) plus binned bin-counts. Watts-per-ray, mirror reflectivity,
shading/blocking and DNI are all applied at *read* time, so any of them can be
revised without re-tracing a multi-hour sweep.

**Flux maps add linearly.** The whole-field map is a weighted sum of the
per-heliostat maps (verified to 0.000e+00), so "what do all heliostats produce"
never needs another trace, and per-heliostat efficiencies drop in as weights.

**Optical efficiency depends only on sun direction.** Sun direction at a fixed
site is exactly determined by (declination, hour angle), so two dates with the
same declination give identical optics -- only DNI differs. The ray trace
therefore builds a dimensionless `eta_optical(declination, hour_angle)` surface,
and the annual integral runs over all 8760 hours with a full DNI series. This is
what turns a handful of traced days into a real annual number.

A corollary: **tracing half a year is enough** (December solstice to June
solstice sweeps the full declination range), and the two equinoxes are nearly
redundant -- they differ by 0.14 deg of declination. `energy.suggest_sweep_dates()`
picks dates by declination coverage instead.

**Shading and blocking are computed analytically**, not with Quadoa blocker
geometry -- the same physics for opaque flat rectangles, at a tiny fraction of
the cost, and revisable without re-tracing. At 8.8 deg sun elevation the full
645-heliostat field loses ~44% to mutual shading and blocking, so this is a
first-order effect rather than a refinement.

Two properties of that calculation are worth knowing before arguing with a
number it produces, and `shading.self_check` pins both:

- **Sunlight is collimated; the outgoing beam is not.** Shading uses one sun
  direction for the whole aperture, which is exact. Blocking gives every point on
  the aperture its own direction to the aim point, because across a 5 m mirror
  aiming 60 m away those directions differ by nearly 5 deg. Treating the outgoing
  beam as collimated understated `eta_block` by ~0.4% of field power.
- **Low-sun shading saturates at the nearest up-sun neighbour.** Heliostats at a
  common height with near-equal normals are parallel planes, so the target maps
  rigidly onto its neighbour: row 2 up-sun shadows a strict subset of row 1 and
  adds exactly nothing. Combined with the ray climbing 1.03 m over a 6 m pitch at
  9.7 deg -- 39% of a mirror's vertical extent, which therefore sees straight over
  its neighbour -- this is why sunrise losses are ~38% and not near-total.

**The secondary shades the field too**, and it is not a small effect. The axicon
is a 30 m wide opaque body with its vertex 27 m up and its rim 5.46 m higher,
directly over a field whose innermost heliostats are at 30 m radius. Where its
shadow lands depends sharply on sun elevation:

| sun elevation | shadow band from the axis | heliostats hit | field power |
|---|---|---|---|
| below ~17° | 158–190 m at 9.7° | 0 | thrown clear of the field |
| 25–40° | 27–46 m | up to 30 | **−3.5%** |
| 50–56° | ~20 m | 11–15 | −1.0 to −1.8% |
| above ~65° | inside 12 m | 0 | falls in the central hole |

`shading.SecondaryCone` does an exact ray-cone intersection rather than treating
the axicon as a disc, because vertex and rim differ by 5.46 m — at 10° sun that
displaces the shadow by 31 m, which is most of the way across the inner field.
It enters shading only, never blocking: the secondary is what every heliostat
aims at, so a beam reaching it is the beam arriving.

`eta_shade` is the union of neighbour and secondary shading, not their product —
a patch shaded by both is lost once. `eta_secondary` is reported alongside so the
secondary's share can be read off on its own.

Because none of this entered the ray trace, a correction to the geometry is a
multiply over the summary rather than hours of Quadoa:

```bash
python -m beamdown rescale analysis_output/full5 --apply
```

It re-solves the pointing and aborts if that disagrees with what the sweep
recorded — a mismatch means the model changed, so the traced spots are stale too
and rescaling would paper over it. `raw/` and `flux/` are never touched, and the
original summary is copied to `summary_before_secondary_shading.csv`.

To check any of it against Quadoa, `occluders` lists the neighbours that actually
occlude one heliostat, in the model's own parameterisation:

```bash
python -m beamdown occluders --output analysis_output/full5 --heliostat 326 --timestep 20260320_1800
```

Each is a copy of the heliostat assembly with four numbers changed — `posx`,
`posy`, `rot_az`, `rot_el` — because an occluding neighbour *is* a heliostat.
Drop the Zernike form and make the surface absorbing, or its reflections land on
the receiver. The rectangle axes were verified against the model's transform
chain rather than assumed: composing `Rz(rot_az) Ry(-rot_el)` with `Rz(90)
Rx(-90)` sends local +x to `normalize(z × n)` and local +y to `n × u`, matching
`shading.py` to twelve decimals.

The classical formulation, projecting each rectangle's corners to the ground and
overlapping the polygons (`shading.corner_shadow`), gives bit-identical answers
and is checked against the sampled one. It has a trap: ground shadows cannot tell
which side of the target a neighbour is on, so overlapping every neighbour counts
the down-sun ones too and roughly doubles the apparent loss (0.29 against the
correct 0.62 at 07:00). The ray formulation gets the up-sun test for free.

**The secondary is pluggable** — see the next section.

## Secondary layouts

`beamdown/secondary/` defines a strategy interface returning pointing and shape
coefficients, selected by `[optics] secondary` in `config.toml`. Three layouts
are registered, and the entire sweep, storage, analysis and figure stack is reused
across all three.

| `secondary =` | secondary optic | reflections | aim point | shades the field |
|---|---|---|---|---|
| `axicon` | cone, 20° half-angle | 2 | **per heliostat**, from its radial position | `SecondaryCone`, exact ray-cone test |
| `prime_focus` | none | 1 | `F1 = (0, 0, focus_height_mm)` | nothing |
| `cassegrain` | hyperboloid | 2 | `F1 = (0, 0, focus_height_mm)` | `SecondaryDisc`, circle at `secondary_rim_height_mm` |

### Mirror figure: focused, flat, fixed

Orthogonal to the layout, and composing with all three of them, is what the
heliostat's own surface does. All three write the same three numbers — the
`z3`/`z4`/`z5` coefficients of the one active Zernike form on `helio_surf`,
whose base radius is `inf` — and nothing else. Pointing tracks the sun
identically in all three cases.

| | how the figure is chosen | flag | manifest |
|---|---|---|---|
| focused | re-solved every timestep for that instant's AOI and slant range | (default) | `flat_mirrors: false`, no `fixed_shapes` |
| flat | forced to zero: a plane | `--flat-mirrors` | `flat_mirrors: true` |
| fixed | per-heliostat, from a CSV, frozen for the whole run | `--fixed-shapes CSV` | `fixed_shapes: "<path>"` |

Focused is an idealisation — no ground glass re-figures itself hourly. Flat is
the opposite bound. **Fixed** is the physical case: one figure per mirror,
ground once, pointing still tracking. Build the table with
`scripts/build_fixed_shapes.py`; it writes
`heliostat,x_mm,y_mm,c3,c4,c5` with `#` metadata lines, and heliostats are
matched by position rounded to 1e-3 mm.

```bash
python -m beamdown sweep --secondary prime_focus --focus-height-mm 36000 \
    --n-mirrors 1 --fixed-shapes data/fixed_shapes_pf36000_mean_cos.csv \
    --all-heliostats --workers 1 --output analysis_output/prime_focus_f36_meancos
```

A heliostat missing from the table **stops the run** — it is never a fall-back to
the solved shape, because a run that mixed ground and re-figured mirrors would
report an annual number between the two and look entirely plausible. Flat and
fixed are refused together (argparse for the two flags, `get_strategy` for the
config-says-flat case) rather than one silently winning; use
`--focused-mirrors --fixed-shapes ...` if `config.toml` sets `flat_mirrors = true`.
An **absent** `fixed_shapes` manifest key means the historical behaviour, the
re-figured mirror, so runs written before the option read exactly as they did.

### The shared-focus contract

`prime_focus` and `cassegrain` share **one aim point for the whole field**: every
heliostat aims at and focuses on the single on-axis point `F1 = (0, 0,
focus_height_mm)`. They therefore share one solver
(`secondary/shared_focus.py`) and produce byte-identical `solve()` output for
identical inputs; everything that distinguishes them lives in the *other* seams —
which body sits over the field, how many reflections `n_mirrors` counts, and which
`.optx` is loaded.

That is the substantive contrast with the axicon, which has **no single aim point
at all**: a cone has no focus, so `axicon.py` computes a different aim point for
each heliostat as a function of radial position (`receiver_correction`), pushed out
along that heliostat's own radial direction. One consequence: the axicon rejects a
heliostat at the field origin, which has no defined radial direction, while the
shared-focus layouts handle it fine — it simply looks straight up the axis.

Neither new layout needs a heliostat shape correction. Prime focus has no second
optic to correct for; the Cassegrain's relay is the hyperboloid's own job, since a
hyperboloid is stigmatic between its two foci. So the Python side **never needs the
hyperboloid's conic constants for pointing** — build the surface by hand in Quadoa
and the trace does the rest. Only the axicon contributes an extra sagittal-only
astigmatism, which is why `axicon_shape_correction` exists and the other two pass
zeros to `mirror.to_quadoa_zernike`.

### The aim-extras contract

Wider than the ABC's one abstract method: **every** strategy must set

```python
extras["aim_x_mm"], extras["aim_y_mm"], extras["aim_z_mm"]
```

to the world point this heliostat is aimed at. `shading.build_geometries` reads
them to build the outgoing beam direction that *blocking* is measured along, and
falls back to `(0, 0, secondary_height_mm)` when a key is missing. That fallback
is silent and wrong for every layout whose aim point is not the secondary vertex,
which is all three of them — a strategy that forgets these keys does not fail, it
reports plausible and incorrect blocking.

`SecondaryStrategy.global_params(geometry)` is the other part of the contract: it
names the model-wide `single_param` values `session.set_global_geometry()` writes.
The base set is `sec_height` and `rec_offset`; the axicon adds `axi_angle`, which a
prime-focus or Cassegrain `.optx` has no parameter for. (Writing a parameter a
model does not have is silently ignored by Quadoa, so this is about honesty rather
than about avoiding an error.)

### Config keys per layout

```toml
[geometry]
secondary_height_mm = 27000.0
receiver_offset_mm  = -20000.0
axicon_angle_deg    = 20.0          # axicon only
focus_height_mm     = 24000.0       # REQUIRED for prime_focus and cassegrain
secondary_rim_height_mm = 20000.0   # REQUIRED for cassegrain (shadow-circle height)

[optics]
secondary   = "cassegrain"
n_mirrors   = 2                     # 1 for prime_focus, 2 for axicon/cassegrain
```

`load_config` raises naming the missing key if a shared-focus layout has no
`focus_height_mm` (or a Cassegrain no `secondary_rim_height_mm`). A mismatched
`n_mirrors` is a **warning, not an error, and is never auto-corrected**:
`optics.throughput` is applied when a stored run is *read*, not when it is written,
so silently flipping `n_mirrors` would rescale the reported numbers for every
existing run in `analysis_output/`. Decide which runs that should affect, then edit
it yourself.

`axicon_aperture_radius_mm` (default 15000) is shared: the Cassegrain hyperboloid
is the same 30 m across as the axicon, so it doubles as the disc's radius. The name
keeps its `axicon_` prefix because renaming it would invalidate every stored run's
config copy.

### What each layout's Quadoa model must contain

**`axicon`** — the shipped `models/heliostat_field_model_mcfg.optx`. Sequence
`secondary_focus`, `single_param`s `sec_height` / `rec_offset` / `axi_angle`, and
per-config `posx` / `posy` / `rot_az` / `rot_el` / `c3` / `c4` / `c5`.

**`prime_focus`** — `models/heliostat_field_prime_focus.optx`: one reflection off
`helio_surf` onto a horizontal detector facing down, at `F1`. No secondary surface
in the path, and no `axi_angle` parameter needed. See "The prime-focus model file"
below for how it is built and what to check before trusting it.

**`cassegrain`** — two reflections: `helio_surf` then a hyperboloid secondary,
onto the existing receiver. Build the hyperboloid by hand (vertex radius, conic
constant, position) so that `F1` is its far focus and the receiver its near one;
nothing in Python reads those numbers. `secondary_rim_height_mm` in config must
match the rim height you gave it, since that is what the shading circle is
projected from.

### The prime-focus model file

`models/heliostat_field_prime_focus.optx` is generated, not hand-drawn:

```
python scripts/build_prime_focus_model.py [--force]
```

It copies `heliostat_field_model_mcfg.optx` and changes **three things and
nothing else** — the builder asserts that, by diffing itself against the base and
refusing any change that is not one of these (6 lines added, 4 removed):

| edit | from | to |
|---|---|---|
| `prime_focus` surface `z` | literal `27000.0` | `single_param` **`pf_height`** = `47000.0` |
| `prime_focus` `float_ap radius` | `0.0` (drawn as a dot) | `2500.0`, drawing only |
| sequence 3 | `sun → helio_surf → secondary → receiver` | `sun → helio_surf → prime_focus` |

The `pf_height` parameter is the point of the exercise. It uses the mechanism the
base model already uses twice — a geometry `<variable>` whose `value` attribute
holds a parameter *id* instead of a number, exactly as `secondary` carries
`value="sec_height"` and `receiver` `value="rec_offset"` — so
`PrimeFocusStrategy.global_params` can write `pf_height = focus_height_mm` on
every session and the detector plane and the Python aim point become the same
number rather than two numbers that happen to agree. The failure mode this
prevents is silent: with the height frozen in the file, a mismatched
`focus_height_mm` still traces, still gives a round spot, still puts it on axis,
and only the *size* is wrong.

**Why 47 m.** Symmetric throw. The axicon's receiver sits 20 m **below** the
axicon vertex (`secondary_height_mm` 27000 + `receiver_offset_mm` −20000 = 7000
mm). The prime-focus receiver is the same 20 m taken the other way, **above** the
vertex: 27000 + 20000 = 47000 mm. The two layouts then differ by the secondary and
by one reflection, not by how far the light travels.

**Why sequence 3 rather than switching `analysis_seq` to 0.** Sequence 0 already
runs this path and is left byte-identical; the builder uses it as the *reference*
and asserts the rewritten sequence 3 is identical to it apart from `name`,
`sequid` and `is_visible`. Rewriting 3 keeps `analysis_seq = 3` valid for every
model file — no new config key, no per-layout branch in the trace stack — and
sequence 3's `setRayDistributionCount1` semantics (a literal ray count, see "Notes
on Quadoa") are the ones actually measured here. Sequence 0 is documented as
literal too, but it is the never-traced column, so that is inherited rather than
measured.

**The detector is unbounded, on purpose.** `float_ap radius` is a clear aperture
for *display and auto-sizing*; what stops a ray is a separate `<aperture>` child
element, which `prime_focus` does not have and must not be given. Measured rather
than assumed: `receiver` carries `float_ap radius = 795.76` and in
`analysis_output/full7` 21.5% of the rays that landed on it are further out than
that, reaching 2510 mm. Spillage stays a post-processing step against
`[receiver] window_mm`.

**Before the first prime-focus sweep**, with a licence seat free:

```
python scripts/verify_prime_focus_model.py
```

Six checks that need Quadoa and so could not be run when the file was built:
`pf_height` reads back as `focus_height_mm` (against the stock model the write is
silently ignored — this is the check that catches the wrong `.optx`);
`getSequenceImageSurface(3)` resolves to the `prime_focus` surface; `setRay
DistributionCount1` is still a literal ray count on sequence 3 (request 300, count
300); `self_test` passes; a heliostat aimed at `(0, 0, 47000)` puts its centroid on
the axis, which is exact by construction since the detector plane *is* the aim
plane; and the same heliostat traced flat throws a far larger spot, which is what
proves the `c3`/`c4`/`c5` writes still reach *this file's* `helio_surf`.

## Tracing occlusion as real geometry

The analytic model above is fast and already validated against Quadoa for a
handful of heliostats, but it is still a scalar multiplied onto a clean trace.
`--occluders` swaps in `models/heliostat_field_occluders.optx`, which carries
14 spare heliostat-shaped slots per mirror (10 shading, 4 blocking) plus one
`ax0` slot for the secondary, so occlusion removes rays in the trace itself
instead of being applied afterwards:

```bash
python -m beamdown sweep --occluders --all-heliostats --output analysis_output/full6
python -m beamdown compare analysis_output/full5 analysis_output/full6 --labels scalar traced
```

`beamdown/occluder_slots.py` decides which neighbours occupy those slots. It
never parks a slot out of the way — a tilted occluder parked at 1 km lost up to
34% of rays to a precision cliff between 3×10⁵ and 5×10⁶ mm, non-monotone with
distance, so every slot holds a real neighbour: the strongest occluders first,
then the nearest heliostat still ahead of the mirror if fewer than 14 actually
occlude it, falling back to a synthetic filler only if nothing is ahead at all.
Occlusion is always ranked over the **full 645-heliostat field**, never the
traced subset, because scoring a downselect against its own neighbours
under-reports shading by construction.

**The axicon's shadow is exactly a circle**, which is what makes it cheap to
trace. A horizontal circle projected along a fixed direction onto a horizontal
plane is congruent to itself, just translated — so above about 20° sun
elevation (where the shadow's rim-vs-vertex spread stays inside the rim's own
radius) the traced geometry is one infinite horizontal plane at 13.5 m carrying
a 15 m circular obscuration, slid to wherever the shadow actually falls. That
plane can also park at any distance with zero error, because translating an
infinite horizontal plane changes nothing but the obscuration's position on it
— unlike the tilted heliostat slots, there is no precision cliff to avoid.

`--occluders` alone traces neighbours only and still applies the secondary as a
scalar (this is `full6`); adding the `ax0` slot removes the secondary's shadow
in the trace too, and once every occluder is real geometry **no scalar
efficiency is applied at all** — the summary's power comes purely from rays
that survived the trace. `scripts/verify_occluder_trace.py` checks this against
the analytic prediction heliostat by heliostat before trusting a full sweep:

```bash
python scripts/verify_occluder_trace.py --run analysis_output/full5
```

Both checks came back clean at full-field, full-year scale. Tracing neighbours
alone (`full6` vs the scalar `full5`) moved field power by +0.48% at a 700 mm
aperture, concentrated at low sun (+2.79% below 20° elevation, +0.10% above
60°) — the analytic model's small conservative bias at the angles where
shading is most severe. Then tracing the axicon as geometry too (`full7` vs
`full6`) moved nothing: every aperture radius from 300–1500 mm agreed to within
±0.00–0.18%, no systematic sign, consistent with ray-count noise. That null
result is the point — it means the closed-form circle shadow already used in
the scalar model was exact, not merely close.

### Scalar vs traced occlusion: what was vetted

Tracing occlusion roughly doubles a sweep's cost, so the question is whether the
remaining comparison sweeps can skip it and apply the scalars in post.
`scripts/vet_occlusion_scalars.py` answers that from the three runs already on
disk — no license, nothing re-traced — and writes CSVs, three figures and a
verdict block to `analysis_output/vet_occlusion/`:

```bash
python scripts/vet_occlusion_scalars.py                       # full5 vs full6 vs full7
python scripts/vet_occlusion_scalars.py --runs full6 full7    # one channel at a time
```

The three runs share one 44-step grid, 645 heliostats and 120,000 rays each, and
their stored analytic `eta_shade`/`eta_block`/`eta_secondary` agree to 0.0 — so
`full5→full6` isolates the neighbour channel and `full6→full7` the secondary's.
The script re-derives from each manifest which efficiency that run's `power_w`
already carries and checks it on all 28,380 rows before comparing anything
(`power_w == rays_landed × scale_factor × eff`, agreement 1e-16): `eta_shade ×
eta_block` for `full5`, `eta_secondary` alone for `full6`, exactly 1 for `full7`.
Getting that wrong would move the answer by the very quantity being measured.

**Annualised collected energy inside the 700 mm aperture (monthly DNI) is 0.338%
lower on the scalar path than on the fully traced one** — 9,755.2 against
9,788.2 MWh, with a Monte-Carlo uncertainty of 0.004% on the annual total. The
field total (no aperture) differs by 0.287%. Per traced date the trapezoid
integral differs by 0.22–0.42%. The difference is one-sided and it is not noise:
125σ on the aggregate, against a per-heliostat noise model validated against the
data itself (predicted spread 0.69%, observed 0.71% on heliostats nothing
occludes).

**The cause is the scalar path's product form, not its geometry.** `eta_shade ×
eta_block` charges twice for any patch of mirror that is both shaded and blocked;
`shading.occlusion_efficiency` unions them instead, which is why it exists.
Re-weighting the *same* analytic model in union form cuts the annual gap from
0.338% to 0.114%, and the low-sun bands from +3.28% to +0.63% (el 5–15°) and
+1.25% to +0.20% (el 15–30°). Above 30° the two forms coincide, because nothing
overlaps there. The analytic geometry's own validation is separate and does not
involve a ray trace: `shading.self_check` compares the sampled shading fraction
against a closed-form rectangle overlap (0.3491 vs 0.3502, 0.3540 vs 0.3536,
0.4117 vs 0.4107).

Where the gap sits: +3.28% below 15° elevation, +0.09% above 30°; +7.98% on
heliostats with `eta_shade × eta_block < 0.5` (1.9% of field power), −0.01% on
the 39% of power collected by heliostats nothing occludes. The neighbour channel
is all of it (`full5→full6` +0.484% of aperture power, reproducing the +0.48%
recorded above; the +2.79% "below 20°" recorded above is the unweighted mean of
the per-timestep deltas — power-weighted over the same timesteps it is +2.68%);
the secondary channel is +0.002% ± 0.004%. But *in aggregate*
is the operative phrase for the secondary: on the 1.18% of rows whose mirror the
axicon's shadow rim actually crosses, the one-number-per-heliostat scalar
disagrees with the trace by −30% to +16% (5–95%, sd 17%, against a 1.1% noise
floor). Those mirrors carry 0.75% of the field's power, so it cancels to +0.003%
of the total — and shows up instead as per-timestep scatter, which is what the
±0.18% recorded above actually was.

**Slot overflow is the one regime where the traced path is the approximation.**
The three compared runs never overflow — their lowest sun is 8.78° and their
sweep logs carry no overflow warning. The 12-date grid in `config.toml` does, at
its lowest step at each end of the day: the probe (analytic, no tracing —
`occlusion_efficiency` over all neighbours against the same union over only the
neighbours that fit the 10+4 slots) puts 250 of 645 heliostats over the limit at
el = 1.75° and the traced model 0.80% high there. Nothing overflows by el =
13.93°. Elevations below 5° carry 0.868% of modelled annual energy, so that
regime is worth ~0.007% of the year.

**The distribution is the caveat.** A scalar multiplies the whole spot; traced
occlusion deletes particular rays, so it changes the spot's shape. At the
receiver that change is real but small, and it is measured here against a null
made of heliostats the analytic model says nothing occludes (there the two runs
differ only by their independent rays). Median per occluded heliostat at el =
8.78°: Δr50 −2.49 mm, Δr90 −6.18 mm, Δrms −3.68 mm, aperture fraction +0.257 pp,
against nulls of +0.20, −0.12, +0.28 and −0.078. At el = 78.6° it is Δr90
−0.48 mm and +0.012 pp. Field-summed over 645 heliostats the worst aperture
fraction shift is 0.414 pp, so scalar-vs-traced barely moves spillage.

That last number is a receiver-plane result and nothing more. These runs store
rays **at the receiver only**, so no plane between the secondary and the receiver
was measured, and none can be inferred from them: at the receiver 645
heliostats' occluded edges land in different places and average out, which is
exactly the reason the totals agree. Through focus they need not. **Any
through-focus work must trace occlusion** — this vet does not license scalars
there, and it makes no claim about that regime.

**What the sweep does now (owner decision, 2026-07-31).** A sweep run without
`--occluders` applies the **union** form: `shading.occlusion_efficiency`, stored
per row as `eta_occlusion` and recorded in the manifest as
`"occlusion_form": "union"`. Every run written before that date used the product
form and has no `occlusion_form` key, which readers must take to mean
`"product"` — the two are different numbers from the same geometry, and the
overlap is not recoverable from `eta_shade` and `eta_block` after the fact.
`store.occlusion_weight_columns` is the single place that turns a manifest into
the columns a reader must multiply into the stored counts; the GUI, `compare`,
`figures` and `rescale` all ask it rather than deciding for themselves.
`--occluders` runs are untouched by any of this: their neighbours are in the ray
path and only the secondary stays a scalar (`"occlusion_form": "traced"`).

The GUI's **"Export with shading + blocking geometry"** button
(`build_occluder_model.build_from_slot_model`) writes the selected heliostat
and timestep into a copy of that same traced model, with its occluder slots
filled in and every other slot parked out of sight — a pure text edit, no
license needed — so what opens in Quadoa is the actual geometry the sweep saw
for that heliostat, not a reconstruction of it.

## Notes on Quadoa

- `import quadoa` requires `os.add_dll_directory(r"C:\Program Files\Quadoa")`
  first; `beamdown.session` handles it.
- Licensing is a USB HASP key with a small number of seats. Avoid opening and
  closing sessions in quick succession -- seats take time to be released, and an
  exhausted key raises a modal dialog. `beamdown.session` never retries a
  license failure for that reason, and the worker pool degrades to fewer workers
  rather than failing.
- **Sequence indices are 0-based** (the GUI shows index+1). Measured for
  `heliostat_field_model_mcfg.optx` by `~/.claude/skills/quadoa-python/scripts/probe_model.py`:

  | py index | GUI | image surf | `setRayDistributionCount1` | receiver spot | use |
  |---|---|---|---|---|---|
  | 0 | 1 | 2 | literal ray count | — | — |
  | 1 | 2 | 3 | **grid density** | ±100 mm | 3D views only |
  | 2 | 3 | 2 | **grid density** | ±575 mm (window) | — |
  | 3 | 4 | 3 | literal ray count | ±700 mm | **analysis** |

  Two traps here. `setRayDistributionCount1` is a *grid density* on sequences 1
  and 2 -- `n=200` yields 31,064 rays, `n=1200` over a million -- and a literal
  ray count on 0 and 3. And sequence 1 gives a far smaller receiver spot than
  sequence 3 at the same surface, so it is fine for looking at geometry in the
  GUI but **must not be used for radiometry**.
- `applyChangesAndInitModel` is **not** needed after writing multiconfig
  parameters: `traceRays` updates the model itself (confirmed in the API docs and
  measured against a Monte-Carlo noise floor).
- Quadoa's tracer already runs at ~4x parallelism, so extra worker processes add
  little on an 8-core machine.
- The API can select configurations and read/write their parameters, but cannot
  create one; `beamdown/model_edit.py` adds columns by editing the `.optx`
  (used by the inspect/export path). The 25-config **figure model** is built
  entirely by guarded text surgery and needs **no licence seat**:
  `python scripts/build_figure_model.py --date D --hour H [--flat] --check`;
  `scripts/verify_figure_model.py` is the seat-gated confirmation, and refuses
  to run under any `analysis_output/.*.lock`.

## Configuration

Everything lives in `config.toml`. The values that most affect cost:

| key | effect |
|---|---|
| `trace.rays_per_heliostat` | 120,000 gives ~45,000 rays on the receiver. Its effect on run time is **unmeasured** — see "Ray budget and iteration count"; `scripts/probe_ray_cost.py` is how to find out |
| `trace.rays_per_trace` | `traceRays` calls per heliostat = `ceil(rays_per_heliostat / this)`; 2 at the shipped values |
| `sweep.dates` | more dates = better annual accuracy; see `suggest_sweep_dates` |
| `storage.raw_rays` | `all` / `downselected` / `none` -- raw rays dominate disk use |
| `receiver.window_mm` | must exceed the largest spot; r90 reaches ~1140 mm at low sun |
| `trace.n_workers` | 1-4, license permitting; 1 is near-optimal here |
