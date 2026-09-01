"""
Form factor convention: f(q, E) = f0(s) + f1(E) + i*f2(E)
where s = q / (4*pi), f0 is the Thomson term, f1/f2 are the Chantler
anomalous corrections.
"""

import re

import numpy as np
from xraydb import chantler_energies, f0, f1_chantler, f2_chantler

_CHARGE_RE = re.compile(r"[0-9]*[+\-]+$")

def log_q_grid(qMin, qMax, n_points):
    """Log-spaced q grid. qMin must be > 0."""
    if qMin <= 0:
        raise ValueError("log_q_grid: qMin must be > 0")
    if qMax <= qMin:
        raise ValueError("log_q_grid: qMax must be > qMin")
    return np.geomspace(qMin, qMax, n_points, dtype=np.float64)

def compute_form_factors(ions, energy, qvals):
    """Form factors for a batch of ions, one entry per unique ion.

    ions : list[str]
    energy : float, eV
    qvals : ndarray of float64, shape (Q,)

    Returns (ff, log):
        ff:     dict[str, ndarray of complex128], one row per unique ion
                that isn't a dummy site
        log:    list[str], one "TIER ion" line per non-full ion, meant to
                be printed by the caller, not written to a file here
    """
    def __strip_charge(ion):
        """Bare element symbol, e.g. 'Fe3+' -> 'Fe'."""
        return _CHARGE_RE.sub("", ion)

    def __classify(element, energy):
        """Which formula applies to this element at this energy.

        Returns "dummy", "f0_only", or "full".
        """
        try:
            f0(element, np.array([0.0]))
        except Exception:
            return "dummy"

        try:
            valid = chantler_energies(element, emin=0, emax=1e9)
        except Exception:
            return "f0_only"
        if len(valid) == 0:
            return "f0_only"

        emin, emax = min(valid), max(valid)
        if energy < emin or energy > emax:
            return "f0_only"

        return "full"

    def __form_factor(ion, s, energy, tier):
        """Complex form factor for one ion, given its tier from __classify().

        s : ndarray of float64, the Cromer-Mann variable q/(4*pi)
        Returns ndarray of complex128, same shape as s.
        """
        element = __strip_charge(ion)
        try:
            f0_ = np.asarray(f0(ion, s), dtype=np.float64)
        except Exception:
            f0_ = np.asarray(f0(element, s), dtype=np.float64)

        if tier == "f0_only":
            return f0_.astype(np.complex128)

        f1_ = float(np.asarray(f1_chantler(element, energy)).squeeze())
        f2_ = float(np.asarray(f2_chantler(element, energy)).squeeze())
        return (f0_ + f1_ + 1j * f2_).astype(np.complex128)

    s = qvals / (4.0 * np.pi)
    ff = {}
    log = []
    for ion in dict.fromkeys(ions):
        element = __strip_charge(ion)
        tier = __classify(element, energy)
        if tier == "dummy":
            log.append("DUMMY   " + ion)
            continue
        if tier == "f0_only":
            log.append("F0-ONLY " + ion)
        ff[ion] = __form_factor(ion, s, energy, tier)
    return ff, log
