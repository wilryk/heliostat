"""Beam-down heliostat field analysis.

Traces each heliostat individually through a Quadoa optical model, stores the
receiver-plane rays, and provides per-heliostat analysis across time of day and
time of year.

Nothing in this package imports quadoa at module level -- only
:mod:`beamdown.session` does, and only when a session is actually created. That
keeps config/solar/field/store/metrics/plots usable without a license.
"""

__all__ = ["config", "solar", "field", "store", "shading", "metrics"]
