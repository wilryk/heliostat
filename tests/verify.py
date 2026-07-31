"""One-shot verification of the beamdown package.

Runs in stages, cheapest and least license-dependent first, so a failure stops
before wasting a Quadoa seat.

    python tests/verify.py            # everything
    python tests/verify.py --no-quadoa   # skip anything needing a license
    python tests/verify.py --api      # also dump the QuadoaCore API docs
"""

from __future__ import annotations

import argparse
import sys
import traceback
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
# The parity stages compare against the original scripts, which now live under
# legacy/ after the repository reorganisation.
sys.path.insert(0, str(REPO / "legacy"))

RESULTS: list[tuple[str, bool, str]] = []


def stage(name):
    def deco(fn):
        def wrapper(*a, **kw):
            print(f"\n{'=' * 70}\n{name}\n{'=' * 70}")
            try:
                ok, note = fn(*a, **kw)
            except Exception as exc:
                traceback.print_exc()
                ok, note = False, f"{type(exc).__name__}: {exc}"
            RESULTS.append((name, ok, note))
            return ok
        return wrapper
    return deco


@stage("1. Config loads and derived quantities are consistent")
def check_config():
    from beamdown.config import load_config

    cfg = load_config()
    # Throughput is asserted as the 0.9^n arithmetic rather than a hardcoded 0.81,
    # because n_mirrors legitimately depends on the layout: prime focus has one
    # reflection, the axicon and Cassegrain have two. Pinning the literal would
    # fail the moment [optics] is switched to prime_focus, for no good reason --
    # what actually needs guarding is that throughput IS reflectivity**n_mirrors
    # and nothing has slipped a fudge factor in.
    expected_throughput = cfg.optics.mirror_reflectivity ** cfg.optics.n_mirrors
    checks = {
        "implied DNI == 1000 W/m2": abs(cfg.source.dni_w_m2 - 1000.0) < 0.01,
        f"throughput == {cfg.optics.mirror_reflectivity}^{cfg.optics.n_mirrors} "
        f"== {expected_throughput:.4g}":
            abs(cfg.optics.throughput - expected_throughput) < 1e-12,
        "receiver height == 7000 mm": abs(cfg.geometry.receiver_height_mm - 7000.0) < 1e-9,
        "source oversizes mirror": cfg.source.aperture_radius_mm > cfg.field.mirror_half_diagonal_mm,
        "chunks sum to ray budget": sum(cfg.trace.chunk_sizes) == cfg.trace.rays_per_heliostat,
        "workers in 1..4": 1 <= cfg.trace.n_workers <= 4,
    }
    for k, v in checks.items():
        print(f"  {'OK ' if v else 'FAIL'}  {k}")
    print(f"  mirror {cfg.field.mirror_width_mm:.0f} x {cfg.field.mirror_height_mm:.0f} mm "
          f"= {cfg.field.mirror_area_m2:.1f} m2, half-diagonal {cfg.field.mirror_half_diagonal_mm:.0f} mm")
    return all(checks.values()), f"{sum(checks.values())}/{len(checks)} checks"


@stage("2. Solar parity with legacy noaa_solar")
def check_solar():
    from beamdown import solar
    from heliostat.noaa_solar import solar_position_calculator as legacy

    lat, lon, tz = -10.0, -52.0, -3
    worst_az = worst_el = 0.0
    n = 0
    for (y, m, d) in [(2026, 3, 20), (2026, 6, 21), (2026, 9, 22), (2026, 12, 21)]:
        for h in np.arange(6.0, 18.01, 0.25):
            oa, oe = legacy(Lat=lat, Lon=lon, TZone=tz, Year=y, Month=m, Day=d, Time=h / 24)
            na, ne = solar.sun_position(lat, lon, tz, y, m, d, h)
            worst_az = max(worst_az, abs(float(oa) - na))
            worst_el = max(worst_el, abs(float(oe) - ne))
            n += 1
    print(f"  {n} comparisons: max |dAz| {worst_az:.2e}, max |dEl| {worst_el:.2e}")
    ok = worst_az < 1e-9 and worst_el < 1e-9
    return ok, f"max dAz {worst_az:.1e}"


@stage("3. Axicon parity with legacy heliostat_shape_solve")
def check_axicon():
    """Legacy parity, with the one deliberate correction held aside.

    The shape now differs from the legacy function by the foreshortening factor
    in ``axicon_shape_correction``. Forcing that factor back to 1 must reproduce
    the legacy values exactly -- which is what makes the departure a single
    physical correction rather than a drift.
    """
    import functools

    from beamdown.config import load_config
    from beamdown import field as F
    from beamdown.secondary import axicon as AX
    from beamdown.secondary import get_strategy
    from heliostat.heliostat_shape_solve import get_heliostat_axicon_shape as legacy

    cfg = load_config()
    fld = F.load_field(cfg)
    # Explicitly "axicon", not cfg.optics.secondary: legacy parity is a property
    # of the axicon strategy and must be checked whatever layout config selects.
    strat = get_strategy("axicon")
    g = cfg.geometry

    def sweep():
        worst, n = np.zeros(5), 0
        for az, el in [(77.9, 50.1), (41.9, 76.7), (294.3, 66.7), (65.6, 3.1), (155.9, 75.2)]:
            for i in range(0, len(fld), 11):
                x, y = float(fld.x_mm[i]), float(fld.y_mm[i])
                old = np.array(legacy(x, y, az, el, g.secondary_height_mm,
                                      g.receiver_height_mm, g.axicon_angle_deg), dtype=float)
                s = strat.solve(x, y, az, el, g)
                new = np.array([s.rot_az_deg, s.rot_el_deg, s.c3, s.c4, s.c5])
                worst = np.maximum(worst, np.abs(old - new))
                n += 1
        return worst, n

    real = AX.axicon_shape_correction
    AX.axicon_shape_correction = functools.partial(real, foreshorten=1.0)
    try:
        flat, n = sweep()
    finally:
        AX.axicon_shape_correction = real
    live, _ = sweep()

    labels = ["rot_az", "rot_el", "c3", "c4", "c5"]
    for lab, f, l in zip(labels, flat, live):
        print(f"  max |delta {lab:6s}| = {f:.3e} with L^2=1,  {l:.3e} as corrected")
    ok = bool((flat < 1e-12).all() and (live[:2] < 1e-12).all() and (live[2:] > 0).all())
    return ok, f"{n} cases, parity {flat.max():.1e}, correction {live[2:].max():.1e}"


# Axicon solve() outputs captured from the code as it stood BEFORE the shared
# mirror-math refactoring (heliostat_orientation / heliostat_shape /
# to_quadoa_zernike moving from axicon.py into mirror.py). Columns:
#
#   x_mm, y_mm, solar_az_deg, solar_el_deg,
#   rot_az_deg, rot_el_deg, c3, c4, c5, aoi_deg, focal_dist_mm
#
# These are asserted bit-for-bit, not to a tolerance, and that is the point: the
# refactoring moved arithmetic between modules and factored the Zernike sign
# conversion into a helper, and the only way to know that was free is to compare
# against values produced before it happened. Do not regenerate this table from
# the current code -- that would make the test tautological. If the axicon physics
# is ever deliberately changed, capture a fresh table from the new code and say so
# in the commit.
_AXICON_GEOMETRY = dict(secondary_height_mm=27000.0, receiver_offset_mm=-20000.0,
                        axicon_angle_deg=20.0)

AXICON_REGRESSION = [
    (30000.0, 0.0, 90.0, 10.0, 0.0, 72.67989092040304, 1.784543137603177e-22, -2.1397021252030787e-06, -2.248270021557599e-06, 62.67989092040304, 60230.167940927604),
    (30000.0, 0.0, 135.0, 35.0, -102.866547578607, 65.03575262178892, 5.232590346648069e-08, -1.5758455644241557e-06, 9.283901153580077e-07, 45.26112320849068, 60230.167940927604),
    (30000.0, 0.0, 180.0, 70.0, -154.32730706858234, 64.32666465385815, -4.3731029461515704e-07, -1.475261900681876e-06, -2.4417775915933826e-07, 24.33947859466284, 60230.167940927604),
    (30000.0, 0.0, 270.0, 45.0, -180.0, 44.820109079596975, -8.915919255884722e-23, -1.454950724954611e-06, -3.63085313293398e-07, 0.17989092040210275, 60230.167940927604),
    (30000.0, 0.0, 0.0, 88.0, 177.1919887420081, 67.28830769706236, 6.013111118814497e-08, -1.48029557689118e-06, -5.263097567326015e-07, 22.69711952136743, 60230.167940927604),
    (0.0, 30000.0, 90.0, 10.0, -35.848470579799006, 35.80138806249052, 1.5631435061635687e-07, -1.4432335513018144e-06, 1.9869502166725274e-07, 41.495814423512094, 60230.167940927604),
    (0.0, 30000.0, 135.0, 35.0, -65.8318791943627, 42.05287691654801, -1.456686865890713e-07, -1.4442073753142393e-06, -2.1981410375353524e-07, 17.69816829175137, 60230.167940927604),
    (0.0, 30000.0, 180.0, 70.0, -90.0, 57.320109079596975, -9.627513166101472e-23, -1.4617325261226177e-06, -4.139947783331158e-07, 12.679890920403025, 60230.167940927604),
    (0.0, 30000.0, 270.0, 45.0, -134.8212351434889, 54.56569804647598, 3.457175094947164e-07, -1.45801311190299e-06, 1.7364214128004508e-07, 30.104076555370618, 60230.167940927604),
    (0.0, 30000.0, 0.0, 88.0, -90.0, 68.32010907959696, -1.1537920567278133e-22, -1.4831903212454369e-06, -5.456874368260213e-07, 23.679890920403032, 60230.167940927604),
    (-45000.0, 0.0, 90.0, 10.0, 0.0, 23.09254352473797, -4.655004695553083e-23, -1.154785798797841e-06, -2.463128433586031e-07, 13.092543524737966, 71682.25509848492),
    (-45000.0, 0.0, 135.0, 35.0, -22.67567487320531, 37.76494472833333, 8.828677931690459e-08, -1.1453725744121647e-06, -9.78347976029614e-08, 18.139751594412868, 71682.25509848492),
    (-45000.0, 0.0, 180.0, 70.0, -22.965138139867133, 60.19152322065841, 3.2247098038913637e-07, -1.1717483083045465e-06, -1.8870664383923365e-07, 28.151860370397454, 71682.25509848492),
    (-45000.0, 0.0, 270.0, 45.0, -4.961204204368009e-14, 85.59254352473793, 1.955369538444707e-21, -1.3222120919155045e-06, -9.434506699389284e-07, 49.40745647526203, 71682.25509848492),
    (-45000.0, 0.0, 0.0, 88.0, 2.475919190313282, 63.06207811527201, -4.2685626244461526e-08, -1.174643737881058e-06, -3.8905435840194544e-07, 26.920220003896855, 71682.25509848492),
    (0.0, -60000.0, 90.0, 10.0, 41.28435746287014, 27.283167457480044, -1.1737371807278219e-07, -9.655786247467331e-07, 2.6300711470297994e-07, 42.49809563525377, 84255.67198602235),
    (0.0, -60000.0, 135.0, 35.0, 26.236225357586168, 59.02709893653628, 2.4436468022101735e-07, -1.0664114130703068e-06, 6.580407873807868e-07, 51.14192879540031, 84255.67198602235),
    (0.0, -60000.0, 180.0, 70.0, 90.0, 70.07585281522704, 1.4606557406907033e-22, -1.0065451747494806e-06, -4.945013130073114e-07, 39.924147184772956, 84255.67198602235),
    (0.0, -60000.0, 270.0, 45.0, 129.27456722202328, 47.274267790547825, -2.3508858017294026e-07, -9.637278154454402e-07, 1.4332178809162725e-07, 34.5979588848767, 84255.67198602235),
    (0.0, -60000.0, 0.0, 88.0, 90.0, 59.07585281522705, 1.6491719865946985e-23, -9.689991724448294e-07, -3.10135395175329e-07, 28.92414718477295, 84255.67198602235),
    (35000.0, 35000.0, 90.0, 10.0, -55.6580407204794, 46.05743593818338, 5.127653446053505e-07, -1.3126362939542097e-06, 9.340503207512243e-07, 59.29907381595013, 75359.02425091648),
    (35000.0, 35000.0, 135.0, 35.0, -90.28754288710606, 44.27378456730055, -1.5911813794050444e-07, -1.0923930260373823e-06, 2.2511482572365105e-07, 35.60456328108742, 75359.02425091648),
    (35000.0, 35000.0, 180.0, 70.0, -122.2553085723951, 53.86226893050895, -1.5102818048804293e-07, -1.0917663501633727e-06, -2.2460684608577807e-07, 21.64768100706252, 75359.02425091648),
    (35000.0, 35000.0, 270.0, 45.0, -155.63997687933752, 41.810017920311, 1.218190679241386e-07, -1.0805266132279704e-06, -9.790739701337715e-08, 17.912329130184, 75359.02425091648),
    (35000.0, 35000.0, 0.0, 88.0, -136.76082541274258, 62.77440163983646, 2.971526417579401e-08, -1.1074825093918709e-06, -3.7616783515449293e-07, 28.62994128902863, 75359.02425091648),
    (-40000.0, 55000.0, 90.0, 10.0, -25.445609169666245, 20.915362716676714, 4.6822829078277124e-08, -8.628965994025585e-07, 9.086124903066353e-09, 26.789294548060447, 91268.28223361552),
    (-40000.0, 55000.0, 135.0, 35.0, -49.662535865262825, 31.3909858128712, -1.3100283079587079e-08, -8.654912451415594e-07, -1.052255951955455e-07, 5.313548085430169, 91268.28223361552),
    (-40000.0, 55000.0, 180.0, 70.0, -63.78942782982399, 49.945203146586245, 9.780203530725274e-08, -8.742673606076854e-07, -1.7948662904312555e-07, 23.5468018723623, 91268.28223361552),
    (-40000.0, 55000.0, 270.0, 45.0, -104.55110014818497, 57.69551314442608, 3.7820632028899874e-07, -9.341900977858411e-07, 3.39817759046437e-07, 46.16345593407874, 91268.28223361552),
    (-40000.0, 55000.0, 0.0, 88.0, -52.60176794921531, 59.61090639363924, -1.9941304953584872e-08, -8.896104082485114e-07, -3.090462000293208e-07, 31.998810052643634, 91268.28223361552),
    (70000.0, -20000.0, 90.0, 10.0, 63.5072660554011, 65.94836360362211, -1.0132272299087805e-06, -1.412020860567293e-06, 1.2721258684423717e-06, 70.26831943437462, 95541.37204459061),
    (70000.0, -20000.0, 135.0, 35.0, -130.34744962670115, 66.7452082836899, -3.786546655409803e-07, -9.966222129589396e-07, 7.074320732232033e-07, 56.412384085894274, 95541.37204459061),
    (70000.0, -20000.0, 180.0, 70.0, -173.6638680365169, 57.89972855833111, -2.766675578989876e-07, -8.494261973119756e-07, -1.54378155839683e-07, 35.30471577694688, 95541.37204459061),
    (70000.0, -20000.0, 270.0, 45.0, 171.07946672860967, 35.9066401638107, -3.317457972075252e-08, -8.228659604292586e-07, -1.0219612595969006e-07, 11.33216660228582, 95541.37204459061),
    (70000.0, -20000.0, 0.0, 88.0, 161.93373844196137, 57.84397484819247, 2.868335840280505e-08, -8.433617921758261e-07, -2.811108556944454e-07, 31.586932298233947, 95541.37204459061),
    (-80000.0, -80000.0, 90.0, 10.0, 22.04743987200361, 15.406042710641094, -1.4025748806377116e-08, -5.731893455358436e-07, 1.4879034495562867e-08, 22.161529553181317, 132910.69428631233),
    (-80000.0, -80000.0, 135.0, 35.0, 4.168752096374046, 35.44985136909887, 1.2949175360853124e-07, -5.905682165147071e-07, 1.5524038806829074e-07, 39.73833426359954, 132910.69428631233),
    (-80000.0, -80000.0, 180.0, 70.0, 26.093424936209956, 59.321406841048, 2.2778686840082194e-07, -6.107872061654443e-07, -1.9380717386208012e-07, 42.994265969182294, 132910.69428631233),
    (-80000.0, -80000.0, 270.0, 45.0, 93.14276930421707, 56.791043156344365, -3.7312138523633044e-07, -6.513462899658459e-07, 2.293697155453335e-07, 52.203891186786365, 132910.69428631233),
    (-80000.0, -80000.0, 0.0, 88.0, 46.453404586405895, 53.56140718853963, -1.5095853252309907e-08, -5.915010665980212e-07, -2.0699080755738978e-07, 35.01234329157888, 132910.69428631233),
    (120000.0, 15000.0, 90.0, 10.0, -71.86451811819356, 75.32669642823457, 7.569639229989504e-07, -1.2144669527147255e-06, 1.3385123318892945e-06, 75.78077688352883, 140323.62463459908),
    (120000.0, 15000.0, 135.0, 35.0, -117.74227236240857, 47.99752773868576, -2.421316752530123e-07, -6.205522659891654e-07, 3.5537609142574873e-07, 53.924381626194936, 140323.62463459908),
    (120000.0, 15000.0, 180.0, 70.0, -154.05636676101943, 49.71582675554827, -1.5780538546087086e-07, -5.574235264409218e-07, -1.0625027458086362e-07, 35.551934645887144, 140323.62463459908),
    (120000.0, 15000.0, 270.0, 45.0, -175.9083527895807, 31.324822963633256, 1.116801089843342e-08, -5.426600498117973e-07, -5.987336407246383e-08, 14.043430024756137, 140323.62463459908),
    (120000.0, 15000.0, 0.0, 88.0, -174.9646027321702, 53.869998557991536, 2.140334159944759e-08, -5.603219791704686e-07, -2.054580797162678e-07, 36.352727111170516, 140323.62463459908),
    (-15000.0, 120000.0, 90.0, 10.0, -40.61888616858797, 18.109319994620602, 2.541004840729463e-08, -5.541137293625995e-07, 1.6789301487998955e-07, 40.13998436866743, 140323.62463459908),
    (-15000.0, 120000.0, 135.0, 35.0, -65.42641464493298, 27.555638078911397, -4.453617605928984e-08, -5.419551122954291e-07, -1.4032155130815764e-08, 18.931572474599154, 140323.62463459908),
    (-15000.0, 120000.0, 180.0, 70.0, -84.75441034206999, 43.8197234767116, 1.3176987842145093e-08, -5.478621337277911e-07, -1.2199133988214364e-07, 26.314173497684425, 140323.62463459908),
    (-15000.0, 120000.0, 270.0, 45.0, -121.8986302320969, 42.15093061526142, 1.9266203129926305e-07, -5.638717075132096e-07, 1.1633413786444586e-07, 41.27566542271004, 140323.62463459908),
    (-15000.0, 120000.0, 0.0, 88.0, -82.60504939403178, 54.76849131863697, -2.8548027490780098e-09, -5.620184813479643e-07, -2.155715234531334e-07, 37.21564814181011, 140323.62463459908),
]


@stage("3b. Axicon solve() unchanged by the shared-mirror-math refactor")
def check_axicon_regression():
    """Bit-level regression against values captured before the refactoring.

    Stage 3 proves the axicon matches the *legacy* solver, which is the stronger
    statement -- but only for ``rot_az``, ``rot_el`` and the L^2=1 shape. This
    stage covers the live, corrected coefficients and the diagnostics too, and it
    is the specific proof that factoring the Zernike sign conversion out of
    ``axicon.py`` into ``mirror.to_quadoa_zernike`` cost nothing.
    """
    from beamdown.config import Geometry
    from beamdown.secondary import get_strategy

    geom = Geometry(**_AXICON_GEOMETRY)
    strat = get_strategy("axicon")

    labels = ["rot_az", "rot_el", "c3", "c4", "c5", "aoi", "focal"]
    worst = [0.0] * len(labels)
    exact = True
    for row in AXICON_REGRESSION:
        x, y, az, el = row[:4]
        sol = strat.solve(x, y, az, el, geom)
        got = (sol.rot_az_deg, sol.rot_el_deg, sol.c3, sol.c4, sol.c5,
               sol.aoi_deg, sol.focal_dist_mm)
        for k, (want, have) in enumerate(zip(row[4:], got)):
            worst[k] = max(worst[k], abs(want - have))
            if want != have:
                exact = False

    for lab, w in zip(labels, worst):
        print(f"  max |delta {lab:6s}| = {w:.3e}")
    print(f"  {len(AXICON_REGRESSION)} cases, bit-identical to the captured table: {exact}")
    ok = exact and all(w <= 1e-12 for w in worst)
    return ok, f"{len(AXICON_REGRESSION)} cases, bit-identical={exact}"


@stage("4. Field loading and downselect")
def check_field():
    from beamdown.config import load_config
    from beamdown import field as F

    cfg = load_config()
    fld = F.load_field(cfg)
    idx, prov = F.load_or_build_downselect(cfg, fld)
    print(f"  field: {len(fld)} heliostats")
    print(f"  downselect: {prov}, n={idx.size}")
    a = F.downselect(fld, 25, "farthest_point")
    b = F.downselect(fld, 25, "farthest_point")
    reproducible = np.array_equal(a, b)
    print(f"  farthest_point reproducible: {reproducible}")
    ok = len(fld) == 645 and idx.size == cfg.field.n_configs and reproducible
    return ok, f"{len(fld)} heliostats, downselect {idx.size}"


@stage("5. Shading geometry self-check")
def check_shading():
    from beamdown import shading

    ok = shading.self_check(verbose=True)
    return ok, "hand-checkable occlusion cases"


@stage("5b. Secondary strategies: prime_focus, cassegrain, and the disc body")
def check_secondary_strategies():
    """The shared-focus layouts, checked without Quadoa.

    Four independent things, in order of how badly a bug in them would hurt:

    1. the aim point really is the single on-axis F1 for every heliostat;
    2. the pointing angles really do bisect sun and aim, reconstructed from
       ``rot_az``/``rot_el`` the way the rest of the code reconstructs them;
    3. the two layouts are byte-for-byte the same solver, and neither inherits the
       axicon's field-origin restriction (while the axicon keeps it);
    4. :class:`~beamdown.shading.SecondaryDisc` occludes what it should, and its
       drawn silhouette agrees with its own occlusion test.
    """
    from beamdown.config import Geometry
    from beamdown.secondary import get_strategy
    from beamdown import shading as S

    FOCUS_Z = 24000.0
    geom = Geometry(secondary_height_mm=27000.0, receiver_offset_mm=-20000.0,
                    axicon_angle_deg=20.0, focus_height_mm=FOCUS_Z,
                    secondary_rim_height_mm=20000.0)
    F1 = np.array([0.0, 0.0, FOCUS_Z])

    pf = get_strategy("prime_focus")
    cg = get_strategy("cassegrain")

    # ~20 field positions spanning the real field's radial range, including the
    # axis itself, and several sun positions from just-above-horizon to overhead.
    positions = [(0.0, 0.0)]
    for r in (30000.0, 45000.0, 70000.0, 100000.0, 140000.0):
        for th in (0.0, 1.1, 2.6, 4.4):
            positions.append((r * np.cos(th), r * np.sin(th)))
    suns = [(88.0, 9.7), (135.0, 35.0), (180.0, 70.0), (270.0, 45.0), (12.0, 84.0)]

    worst_aim = worst_bisect = worst_focal = worst_identical = 0.0
    n = 0
    for (x, y) in positions:
        for (az, el) in suns:
            a = pf.solve(x, y, az, el, geom)
            b = cg.solve(x, y, az, el, geom)

            # (d) identical solvers -- compared exactly, not to a tolerance.
            for ka, kb in ((a.rot_az_deg, b.rot_az_deg), (a.rot_el_deg, b.rot_el_deg),
                           (a.c3, b.c3), (a.c4, b.c4), (a.c5, b.c5),
                           (a.aoi_deg, b.aoi_deg), (a.focal_dist_mm, b.focal_dist_mm)):
                worst_identical = max(worst_identical, abs(ka - kb))
            for key in a.extras:
                worst_identical = max(worst_identical, abs(a.extras[key] - b.extras[key]))

            # (a) the aim extras ARE F1, exactly -- one point for the whole field.
            aim = np.array([a.extras["aim_x_mm"], a.extras["aim_y_mm"],
                            a.extras["aim_z_mm"]])
            worst_aim = max(worst_aim, float(np.abs(aim - F1).max()))

            # (b) the normal implied by rot_az/rot_el bisects sun and mirror->F1.
            # Reconstructed exactly as heliostat_orientation defines it, and as
            # shading.normal_from_angles reads it back.
            mirror = np.array([x, y, 0.0])
            to_sun = S.sun_vector(az, el)
            to_aim = F1 - mirror
            to_aim = to_aim / np.linalg.norm(to_aim)
            want = to_sun + to_aim
            want = want / np.linalg.norm(want)
            got = S.normal_from_angles(a.rot_az_deg, a.rot_el_deg)
            worst_bisect = max(worst_bisect, float(np.abs(want - got).max()))
            # And the defining property itself: equal angles either side.
            worst_bisect = max(worst_bisect,
                               abs(float(got @ to_sun) - float(got @ to_aim)))

            # (c) focal distance is the plain |F1 - mirror|.
            worst_focal = max(worst_focal,
                              abs(a.focal_dist_mm - float(np.linalg.norm(F1 - mirror))))
            n += 1

    print(f"  {n} cases ({len(positions)} positions x {len(suns)} sun positions)")
    print(f"  max |aim - F1|                      = {worst_aim:.3e} mm")
    print(f"  max normal-bisector residual        = {worst_bisect:.3e}")
    print(f"  max |focal_dist - |F1 - mirror||    = {worst_focal:.3e} mm")
    print(f"  max |prime_focus - cassegrain|      = {worst_identical:.3e}")

    checks = {
        "aim point is exactly (0, 0, focus_height)": worst_aim == 0.0,
        "normal bisects sun and mirror->F1": worst_bisect < 1e-9,
        "focal_dist == |F1 - mirror|": worst_focal < 1e-9,
        "prime_focus and cassegrain solve identically": worst_identical == 0.0,
    }

    # (e)/(f) the field origin: well posed for the shared-focus layouts, still
    # rejected by the axicon, which needs a radial direction it does not have.
    try:
        origin = pf.solve(0.0, 0.0, 135.0, 40.0, geom)
        at_origin = abs(origin.focal_dist_mm - FOCUS_Z) < 1e-9
    except Exception as exc:
        print(f"  prime_focus at the origin RAISED {type(exc).__name__}: {exc}")
        at_origin = False
    checks["prime_focus solves at the field origin"] = at_origin

    try:
        get_strategy("axicon").solve(0.0, 0.0, 135.0, 40.0, geom)
        axicon_raises = False
    except ValueError:
        axicon_raises = True
    checks["axicon still raises at the field origin"] = axicon_raises

    # focus_height_mm is not optional for these layouts.
    try:
        pf.solve(30000.0, 0.0, 135.0, 40.0,
                 Geometry(secondary_height_mm=27000.0, receiver_offset_mm=-20000.0,
                          axicon_angle_deg=20.0))
        checks["missing focus_height_mm raises"] = False
    except ValueError as exc:
        checks["missing focus_height_mm raises"] = "focus_height_mm" in str(exc)

    # -- the disc body ------------------------------------------------------
    disc = S.SecondaryDisc(z_mm=20000.0, radius_mm=15000.0)
    checks["disc rim_height_mm == z_mm"] = disc.rim_height_mm == 20000.0

    # Straight up from the axis: hit. Straight down: the disc is behind, not ahead.
    axis = np.array([[0.0, 0.0, 0.0]])
    checks["disc: straight up from axis hits"] = bool(
        disc.occludes(axis, np.array([0.0, 0.0, 1.0]))[0])
    checks["disc: straight down misses"] = not bool(
        disc.occludes(axis, np.array([0.0, 0.0, -1.0]))[0])
    # Horizontal ray never reaches the plane.
    checks["disc: horizontal ray misses"] = not bool(
        disc.occludes(axis, np.array([1.0, 0.0, 0.0]))[0])
    # Just inside / just outside the rim, straight up: the boundary is the radius.
    checks["disc: just inside the rim hits"] = bool(
        disc.occludes(np.array([[14600.0, 0.0, 0.0]]), np.array([0.0, 0.0, 1.0]))[0])
    checks["disc: just outside the rim misses"] = not bool(
        disc.occludes(np.array([[15400.0, 0.0, 0.0]]), np.array([0.0, 0.0, 1.0]))[0])

    # Rays from below at a range of angles, against the closed form: a ray from
    # (px, 0, 0) along to_sun crosses z = 20000 at radius |px - throw|, where
    # throw = z / tan(el) for a due-east sun. Independent of occludes().
    agree = True
    for el in (5.0, 15.0, 30.0, 45.0, 60.0, 80.0):
        to_sun = S.sun_vector(90.0, el)
        throw = disc.z_mm / np.tan(np.deg2rad(el))
        for delta in (-18000.0, -14000.0, -1000.0, 0.0, 9000.0, 16000.0):
            px = -(throw + delta)
            expect = abs(-px - throw) <= disc.radius_mm
            got = bool(disc.occludes(np.array([[px, 0.0, 0.0]]), to_sun)[0])
            if got != expect:
                agree = False
                print(f"    disc mismatch el={el} delta={delta}: {got} != {expect}")
    checks["disc: occlusion matches closed form from below"] = agree

    # The drawn silhouette must be the same set the occlusion test uses. Sample a
    # grid at mirror height, ask occludes(), and ask point-in-polygon of the
    # silhouette; they must label every sample the same way.
    to_sun = S.sun_vector(112.0, 38.0)
    poly = S.disc_shadow(disc, to_sun, ground_z=0.0, n=512)
    g = np.linspace(-90000.0, 90000.0, 61)
    gx, gy = np.meshgrid(g, g, indexing="ij")
    pts = np.column_stack([gx.ravel(), gy.ravel(), np.zeros(gx.size)])
    by_ray = disc.occludes(pts, to_sun)
    centre = poly.mean(axis=0)
    by_poly = (np.hypot(pts[:, 0] - centre[0], pts[:, 1] - centre[1])
               <= disc.radius_mm)
    disagree = int((by_ray != by_poly).sum())
    print(f"  disc silhouette vs occludes(): {disagree} of {len(pts)} samples "
          f"disagree, {int(by_ray.sum())} shaded")
    checks["disc: silhouette agrees with occludes()"] = disagree == 0 and by_ray.any()

    # And that the factory hands out the right body per layout.
    class _Cfg:
        def __init__(self, layout):
            self.geometry = geom

            class _O:
                secondary = layout
            self.optics = _O()

    bodies = {lay: S.secondary_body(_Cfg(lay))
              for lay in ("axicon", "cassegrain", "prime_focus")}
    checks["factory: axicon -> cone"] = isinstance(bodies["axicon"], S.SecondaryCone)
    checks["factory: cassegrain -> disc"] = isinstance(bodies["cassegrain"], S.SecondaryDisc)
    checks["factory: prime_focus -> None"] = bodies["prime_focus"] is None
    checks["prime_focus draws no silhouette"] = len(
        S.secondary_shadow(None, to_sun)) == 0

    for k, v in checks.items():
        print(f"  {'OK ' if v else 'FAIL'}  {k}")
    return all(checks.values()), f"{sum(checks.values())}/{len(checks)} checks, {n} solves"


@stage("5c. Flat heliostats: only the shape changes, and no path can skip it")
def check_flat_mirrors():
    """``[optics] flat_mirrors``, checked on all three layouts at once.

    Four things, in the order a bug in them would hurt:

    1. **Flat really is flat.** ``c3 == c4 == c5 == 0.0`` exactly -- not "small",
       exactly zero, because they are written straight into the model's active
       ``<form type="zernike">`` whose base surface has ``radius = inf``. Zero
       coefficients there leave a plane; anything else leaves a lens.
    2. **Nothing else moved.** ``rot_az``, ``rot_el``, ``aoi``, ``focal_dist``,
       ``cosine_efficiency`` and every ``extras`` key (which includes the
       ``aim_*_mm`` the *blocking* test is measured along) must be BIT-identical
       to the focused solve. The flag is a shape flag; if it also nudged the
       pointing, the flat-vs-focused comparison would be measuring two changes.
    3. **The seam is the only seam.** ``get_strategy`` is the sole constructor of
       a strategy anywhere in the package, and it is given the whole config, so a
       flat run cannot reach a solve() through a code path that forgot. Checked
       statically over the source rather than by trusting the grep of the day.
    4. **The wrapper delegates.** ``global_params`` (which decides what gets
       written into the .optx) must not be swallowed by the wrapper.
    """
    import ast

    from beamdown.config import Config, Geometry, OpticsSpec, load_config
    from beamdown.secondary import FlatHeliostats, get_strategy
    from beamdown.secondary.base import SecondaryStrategy

    geom = Geometry(secondary_height_mm=27000.0, receiver_offset_mm=-20000.0,
                    axicon_angle_deg=20.0, focus_height_mm=24000.0,
                    secondary_rim_height_mm=20000.0)

    # The axicon is undefined at the field origin, so the grid starts off-axis.
    positions = []
    for r in (30000.0, 45000.0, 70000.0, 100000.0, 140000.0):
        for th in (0.0, 1.1, 2.6, 4.4, 5.9):
            positions.append((r * np.cos(th), r * np.sin(th)))
    suns = [(88.0, 9.7), (135.0, 35.0), (180.0, 70.0), (270.0, 45.0),
            (12.0, 84.0), (315.0, 22.5)]

    checks: dict[str, bool] = {}
    n = 0
    worst_shape = 0.0
    unchanged = True
    curved_anywhere = False
    first_bad = ""

    for layout in ("axicon", "prime_focus", "cassegrain"):
        focused = get_strategy(layout)
        flat = get_strategy(layout, flat=True)
        for (x, y) in positions:
            for (az, el) in suns:
                a = focused.solve(x, y, az, el, geom)
                b = flat.solve(x, y, az, el, geom)
                n += 1

                # Exactly zero, not "small": these three go straight into the
                # model's active zernike form, and anything non-zero there is a
                # lens on top of a radius = inf plane.
                worst_shape = max(worst_shape, abs(b.c3), abs(b.c4), abs(b.c5))
                if a.c3 or a.c4 or a.c5:
                    curved_anywhere = True

                # Everything that is NOT the shape, compared with ==, not a
                # tolerance: the flag must be a pure projection of the solution.
                for field_name in ("rot_az_deg", "rot_el_deg", "aoi_deg",
                                   "focal_dist_mm", "cosine_efficiency"):
                    if getattr(a, field_name) != getattr(b, field_name):
                        unchanged = False
                        first_bad = first_bad or f"{layout}.{field_name}"
                if set(a.extras) != set(b.extras):
                    unchanged = False
                    first_bad = first_bad or f"{layout}.extras keys"
                for k in a.extras:
                    if a.extras[k] != b.extras[k]:
                        unchanged = False
                        first_bad = first_bad or f"{layout}.extras[{k}]"

    print(f"  {n} solves ({len(positions)} positions x {len(suns)} suns x 3 layouts)")
    print(f"  max |c3|,|c4|,|c5| with flat on = {worst_shape:.3e} (want exactly 0)")
    print(f"  pointing/diagnostics bit-identical to the focused solve: {unchanged}"
          + (f"  first difference at {first_bad}" if first_bad else ""))
    checks["flat forces c3 = c4 = c5 = 0.0 exactly"] = worst_shape == 0.0
    checks["flat changes nothing but the shape, bit for bit"] = unchanged
    # A control: if the focused solve were also flat the test above would pass
    # vacuously.
    checks["the focused solve is genuinely curved"] = curved_anywhere

    # -- the flag has to arrive from the config, not just from flat=True -------
    def _cfg(layout, flat):
        base = load_config()
        return Config(
            site=base.site, geometry=geom, field=base.field,
            optics=OpticsSpec(secondary=layout, mirror_reflectivity=0.9,
                              n_mirrors=1 if layout == "prime_focus" else 2,
                              flat_mirrors=flat),
            source=base.source, dni=base.dni, trace=base.trace, sweep=base.sweep,
            receiver=base.receiver, storage=base.storage, repo_root=base.repo_root,
        )

    by_config = True
    for layout in ("axicon", "prime_focus", "cassegrain"):
        on = get_strategy(_cfg(layout, True))
        off = get_strategy(_cfg(layout, False))
        by_config &= isinstance(on, FlatHeliostats) and not isinstance(off, FlatHeliostats)
        sol = on.solve(45000.0, 12000.0, 135.0, 35.0, geom)
        by_config &= (sol.c3, sol.c4, sol.c5) == (0.0, 0.0, 0.0)
        # The name stays the layout's, so anything matching on it still matches;
        # describe() is what says "flat".
        by_config &= on.name == layout and "flat" in on.describe()
        # global_params must survive the wrapper: for the axicon that is the
        # cone half-angle, and losing it would leave the .optx unwritten.
        by_config &= on.global_params(geom) == off.global_params(geom)
    checks["a Config with flat_mirrors = true produces flat solutions"] = by_config
    checks["OpticsSpec defaults to focused"] = OpticsSpec(
        secondary="axicon", mirror_reflectivity=0.9, n_mirrors=2).flat_mirrors is False

    # -- (3) no other way to build a strategy ---------------------------------
    #
    # Every solve() in the package goes through an object that get_strategy
    # returned, and get_strategy is given the whole config, so it cannot be
    # unaware of the flag. What could reintroduce the hole is a call site
    # reverting to get_strategy(cfg.optics.secondary) -- passing the NAME drops
    # the flag on the floor and silently keeps the curvature. Parsed, not
    # grepped, so a match inside a string or a comment cannot fool it.
    package = sorted((REPO / "beamdown").rglob("*.py"))
    offenders, n_calls = [], 0
    direct_construction = []
    strategy_classes = {"AxiconStrategy", "PrimeFocusStrategy", "CassegrainStrategy",
                        "SharedFocusStrategy"}
    for path in package:
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call) or not isinstance(node.func, ast.Name):
                continue
            if node.func.id in strategy_classes and path.name not in (
                    "axicon.py", "prime_focus.py", "cassegrain.py"):
                direct_construction.append(f"{path.name}:{node.lineno}")
            if node.func.id != "get_strategy":
                continue
            n_calls += 1
            arg = node.args[0] if node.args else None
            # Accept cfg / self.cfg / cfg.optics. Reject anything ending in
            # .secondary, which is the name-only form the flag cannot ride on.
            if isinstance(arg, ast.Attribute) and arg.attr == "secondary":
                offenders.append(f"{path.relative_to(REPO)}:{node.lineno}")
    print(f"  {n_calls} get_strategy() call(s) across {len(package)} package modules; "
          f"{len(offenders)} pass a bare layout name")
    for bad in offenders:
        print(f"    LEAK  {bad} passes the layout NAME, so flat_mirrors is dropped")
    checks["no package call site passes the layout name alone"] = not offenders
    checks["get_strategy is actually used (>= 8 call sites)"] = n_calls >= 8
    checks["no module instantiates a strategy class directly"] = not direct_construction

    # -- (4) the wrapper is a faithful SecondaryStrategy -----------------------
    wrapped = get_strategy("axicon", flat=True)
    checks["FlatHeliostats is a SecondaryStrategy"] = isinstance(wrapped, SecondaryStrategy)
    checks["flat=False beats a flat config"] = not isinstance(
        get_strategy(_cfg("axicon", True), flat=False), FlatHeliostats)
    checks["flat=True beats a focused config"] = isinstance(
        get_strategy(_cfg("axicon", False), flat=True), FlatHeliostats)
    # A bare name still means "focused", which is what the parity stages rely on.
    checks["a bare layout name is still focused"] = not isinstance(
        get_strategy("axicon"), FlatHeliostats)

    for k, v in checks.items():
        print(f"  {'OK ' if v else 'FAIL'}  {k}")
    return all(checks.values()), f"{sum(checks.values())}/{len(checks)} checks, {n} solves"


@stage("6. Store round-trip (quantisation, index, rebin)")
def check_store():
    import shutil
    import pandas as pd
    from beamdown.config import load_config
    from beamdown.store import RunStore, TimestepResult

    cfg = load_config()
    tmp = REPO / "_verify_store_tmp"
    if tmp.exists():
        shutil.rmtree(tmp)

    rng = np.random.default_rng(0)
    grid = cfg.receiver.grid_size
    n_helio, n_rays = 3, 5000

    store = RunStore(tmp, cfg=cfg, mode="w")
    store.write_manifest({"n_heliostats": n_helio})

    truth, chunks, index = [], [], np.zeros((n_helio, 3), dtype=np.int64)
    counts = np.zeros((n_helio, grid, grid), dtype=np.uint32)
    cursor = 0
    for i in range(n_helio):
        xy = rng.normal(0, 120, size=(n_rays, 2)).astype(np.float32) + i * 50
        truth.append(xy)
        q = RunStore.quantise(xy, cfg.receiver.window_mm)
        chunks.append(q)
        index[i] = (i, cursor, q.shape[0])
        cursor += q.shape[0]
        edges = cfg.receiver.edges
        c, _, _ = np.histogram2d(xy[:, 1], xy[:, 0], bins=[edges, edges])
        counts[i] = c.astype(np.uint32)

    rows = pd.DataFrame({"date": ["2026-03-20"] * n_helio, "hour": [12.0] * n_helio,
                         "heliostat_id": range(n_helio), "power_w": [1.0] * n_helio})
    store.write_timestep(TimestepResult(
        key="20260320_1200", date="2026-03-20", hour=12.0,
        solar_az_deg=41.9, solar_el_deg=76.7,
        heliostat_ids=np.arange(n_helio), rays_emitted=120000,
        counts=counts, rays=np.concatenate(chunks), index=index, rows=rows,
    ))

    reader = RunStore(tmp, cfg=cfg, mode="r")
    back = reader.read_rays("20260320_1200", heliostat_id=1)
    err = np.abs(back - truth[1]).max()
    resolution = cfg.receiver.window_mm / 32767.0
    print(f"  quantisation round-trip max error: {err:.5f} mm (step {resolution:.5f} mm)")
    print(f"  counts sum matches rays: {int(counts[1].sum())} vs {n_rays}")

    rebinned = reader.rebin("20260320_1200", 32, cfg.receiver.window_mm, heliostat_id=1)
    print(f"  rebin to 32x32: shape {rebinned.shape}, total {int(rebinned.sum())}")

    field_flux = reader.field_flux("20260320_1200")
    print(f"  field_flux: shape {field_flux.shape}, peak {field_flux.max():.1f} W/m2")

    ok = err <= resolution and rebinned.shape == (32, 32)
    shutil.rmtree(tmp, ignore_errors=True)
    return ok, f"max quantisation error {err:.5f} mm"


@stage("7. Metrics on a synthetic spot")
def check_metrics():
    from beamdown.config import load_config
    from beamdown import metrics

    cfg = load_config()
    rng = np.random.default_rng(1)
    sigma = 100.0
    xy = rng.normal(0, sigma, size=(200000, 2)).astype(np.float32)

    m = metrics.spot_metrics(xy, 500000, cfg, aperture_radius_mm=300.0)
    # For a 2D Gaussian, r50 = sigma*sqrt(2 ln 2), r90 = sigma*sqrt(2 ln 10)
    exp50 = sigma * np.sqrt(2 * np.log(2))
    exp90 = sigma * np.sqrt(2 * np.log(10))
    print(f"  r50 {m['r50_mm']:.2f} mm (analytic {exp50:.2f})")
    print(f"  r90 {m['r90_mm']:.2f} mm (analytic {exp90:.2f})")
    print(f"  rms {m['rms_radius_mm']:.2f} mm (analytic {sigma*np.sqrt(2):.2f})")
    print(f"  power {m['power_w']:.1f} W, spillage {m['spillage']:.4f}")
    ok = abs(m["r50_mm"] - exp50) < 2.0 and abs(m["r90_mm"] - exp90) < 5.0
    return ok, f"r50 within {abs(m['r50_mm']-exp50):.2f} mm of analytic"


@stage("8. Date selection and annual-energy scaffolding")
def check_energy():
    from beamdown.config import load_config
    from beamdown import energy, solar

    cfg = load_config()
    canon = cfg.sweep.dates
    print(f"  canonical {len(canon)} dates -> "
          f"{energy.distinct_declinations(cfg, canon)} distinct declinations")

    half = energy.suggest_sweep_dates(cfg, 8, branch="ascending", must_include=canon)
    steps = solar.build_time_grid(cfg, half)
    print(f"  suggested  {len(half)} dates -> "
          f"{energy.distinct_declinations(cfg, half)} distinct declinations, "
          f"{len(steps)} timesteps")
    print("   ", ", ".join(str(d) for d in half))

    hy = solar.hours_of_year(cfg, 2026)
    print(f"  hours_of_year: {len(hy)} rows, {(hy.solar_el_deg > 0).sum()} daylight")
    ok = len(hy) == 8760 and energy.distinct_declinations(cfg, half) >= 6
    return ok, f"{len(half)} dates, {energy.distinct_declinations(cfg, half)} declinations"


@stage("8b. Time-grid sampling: uniform daylight window")
def check_time_grid_sampling():
    """build_time_grid now samples uniformly across the daylight window instead
    of snapping to an hour_step clock grid (see beamdown/solar.py). These checks
    guard the arithmetic directly against monkeypatching build_time_grid."""
    import datetime as _dt

    from beamdown.config import load_config
    from beamdown import solar

    checks = {}

    # Worked example from the task: sunrise 07:00, sunset 17:50, 10 min margin,
    # 1h max spacing -> 10.5h span -> ceil(10.5/1.0)=11 intervals -> 12 points.
    rise, set_, margin, step = 7.0, 17.0 + 50.0 / 60.0, 10.0 / 60.0, 1.0
    hours = solar._sample_hours(rise, set_, margin, step)
    diffs = np.diff(hours)
    checks["worked example: 12 points"] = len(hours) == 12
    checks["worked example: first 07:10"] = abs(hours[0] - (7.0 + 10.0 / 60.0)) < 1e-9
    checks["worked example: last 17:40"] = abs(hours[-1] - (17.0 + 40.0 / 60.0)) < 1e-9
    checks["worked example: uniform spacing"] = np.ptp(diffs) < 1e-9
    print(f"  worked example: {len(hours)} points, first {hours[0]:.4f}h "
          f"last {hours[-1]:.4f}h spacing {diffs[0] * 60:.2f} min")

    # Every real configured date: samples stay inside [rise+margin, set-margin],
    # spacing never exceeds hour_step, and hours strictly increase.
    cfg = load_config()
    site = cfg.site
    m = cfg.sweep.sunrise_margin_min / 60.0
    within_window = within_step = increasing = True
    for date in cfg.sweep.dates:
        r, s = solar.sunrise_sunset(site.latitude, site.longitude, site.timezone,
                                     date.year, date.month, date.day)
        h = solar._sample_hours(r, s, m, cfg.sweep.hour_step)
        if len(h) == 0:
            continue
        if h[0] < r + m - 1e-9 or h[-1] > s - m + 1e-9:
            within_window = False
        if len(h) > 1:
            d = np.diff(h)
            increasing = increasing and bool((d > 0).all())
            within_step = within_step and bool((d <= cfg.sweep.hour_step + 1e-9).all())
    checks["configured dates: samples within window"] = within_window
    checks["configured dates: spacing <= hour_step"] = within_step
    checks["configured dates: hours strictly increasing"] = increasing

    ts = solar.TimeStep(date=_dt.date(2026, 3, 20), hour=6.7, solar_az_deg=0.0, solar_el_deg=0.0)
    checks["TimeStep.key minute round-trip"] = ts.key == "20260320_0642"

    for k, v in checks.items():
        print(f"  {'OK ' if v else 'FAIL'}  {k}")
    return all(checks.values()), f"{sum(checks.values())}/{len(checks)} checks"


@stage("9. Quadoa: environment, license seats, and API")
def check_quadoa_env(dump_api=False):
    from beamdown.config import load_config
    from beamdown import session as S

    cfg = load_config()
    procs = S.running_quadoa_processes()
    print(f"  running quadoa/hasp/python processes: {len(procs)}")
    for p in procs[:12]:
        print(f"    {p['name']} (pid {p['pid']})")

    sess = S.QuadoaSession(cfg)
    print("  session opened OK")

    if dump_api:
        text = S.describe_api(cfg.trace.quadoa_folder,
                              ["config", "multiconf", "applyChanges", "init", "trace", "rayPos"])
        (REPO / "quadoa_api_reference.txt").write_text(text, encoding="utf-8")
        print(f"  API docs -> quadoa_api_reference.txt ({len(text.splitlines())} lines)")

    sess.close()
    return True, f"{len(procs)} related processes"


@stage("10. Quadoa: parameters actually reach the trace")
def check_param_writes():
    from beamdown.config import load_config
    from beamdown import session as S
    from beamdown.secondary import get_strategy

    cfg = load_config()
    strat = get_strategy(cfg)
    sess = S.QuadoaSession(cfg)

    print("  -- is applyChangesAndInitModel needed? --")
    needed = sess.check_reinit_needed(strat)

    print("  -- do distinct heliostats give distinct spots? --")
    res = sess.self_test(strat)
    sess.close()
    return res["passed"], f"separation {res['separation_mm']:.1f} mm, reinit_needed={needed}"


@stage("11. Quadoa: concurrent session capacity")
def check_workers(confirm=False):
    import gc

    from beamdown.config import load_config
    from beamdown.sweep import check_worker_capacity

    gc.collect()  # make sure stages 9-10 have released their seats
    cfg = load_config()
    usable = check_worker_capacity(cfg, max_workers=4, confirm=confirm)
    if usable is None:
        return True, "skipped (use --workers to probe)"
    return usable >= 1, f"{usable} concurrent sessions usable"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-quadoa", action="store_true", help="skip license-dependent stages")
    ap.add_argument("--api", action="store_true", help="dump QuadoaCore API docs")
    ap.add_argument("--workers", action="store_true",
                    help="probe concurrent session capacity (pops modal license dialogs)")
    args = ap.parse_args()

    check_config()
    check_solar()
    check_axicon()
    check_axicon_regression()
    check_field()
    check_shading()
    check_secondary_strategies()
    check_flat_mirrors()
    check_store()
    check_metrics()
    check_energy()
    check_time_grid_sampling()

    if not args.no_quadoa:
        if check_quadoa_env(dump_api=args.api):
            check_param_writes()
            check_workers(confirm=args.workers)

    print(f"\n{'=' * 70}\nSUMMARY\n{'=' * 70}")
    for name, ok, note in RESULTS:
        print(f"  {'PASS' if ok else 'FAIL'}  {name}\n        {note}")
    failed = [n for n, ok, _ in RESULTS if not ok]
    print(f"\n{len(RESULTS) - len(failed)}/{len(RESULTS)} stages passed")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
