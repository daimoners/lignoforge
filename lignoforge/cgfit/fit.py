"""
Parameter fitting.

* **Bonded terms**: Boltzmann inversion of the mapped reference distributions
  (harmonic wells, see :mod:`lignoforge.cg.boltzmann`).  Specific keys with
  enough samples are fitted on their own; generic keys (bead-type level) are
  always fitted from the pooled samples and serve as fallback.
* **Non-bonded terms**: 12-6 Lennard-Jones per bead-type pair, optimised by
  *relative-entropy minimisation* (Shell, J. Chem. Phys. 2008).  With
  ``U = sum_p [eps_p, sigma_p]``-parametrised pair potentials the gradient of
  the relative entropy is ``beta (<dU/dlambda>_ref - <dU/dlambda>_CG)`` and the
  Hessian is estimated from the frame-to-frame covariance of ``dU/dlambda`` in
  the CG ensemble; a damped, trust-region-limited Newton step is taken in the
  logarithm of the parameters.  Only the pair-distance histograms of the two
  ensembles are needed.
* A **pressure step** rescales all epsilons to bring the CG pressure to the
  atomistic one (the relative-entropy gradient alone does not fix it).
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from lignoforge.cg.boltzmann import KB, _harmonic_from_histogram
from lignoforge.cgfit.model import CGParams
from lignoforge.cgfit.sampling import Stats

K_BOND_MAX, K_ANGLE_MAX = 1.0e5, 1.0e3


# ── bonded ────────────────────────────────────────────────────────────────────

def fit_bonded(st: Stats, temperature: float, min_samples: int = 300,
               ) -> Tuple[Dict[str, dict], Dict[str, dict], dict]:
    kT = KB * temperature
    bonds: Dict[str, dict] = {}
    angles: Dict[str, dict] = {}
    report = {"fitted": {}, "failed": []}
    pooled: Dict[str, List[np.ndarray]] = {}
    arrays = {k: np.concatenate(v) for k, v in st.bonds.items()}
    for k, r in arrays.items():
        pooled.setdefault(k.split(":", 1)[1], []).append(r)
    for key, r in list(arrays.items()) + [(g, np.concatenate(v)) for g, v in pooled.items()]:
        if len(r) < min_samples:
            report["failed"].append(f"bond {key} (n={len(r)})")
            continue
        res = _harmonic_from_histogram(r, lambda c: c ** 2, kT, r.min(), r.max())
        if res is None:
            report["failed"].append(f"bond {key} (no harmonic well)")
            continue
        r0, k = res
        bonds[key] = {"r0": round(float(r0), 4), "k": round(float(min(k, K_BOND_MAX)), 1)}
        report["fitted"][f"bond {key}"] = {**bonds[key], "n": len(r)}
    for key, th in st.angles.items():
        th = np.concatenate(th)
        if len(th) < min_samples:
            report["failed"].append(f"angle {key} (n={len(th)})")
            continue
        res = _harmonic_from_histogram(np.radians(th), np.sin, kT, 0.2, np.pi, bounds=(0.0, np.pi))
        if res is None:
            report["failed"].append(f"angle {key} (no harmonic well)")
            continue
        t0, k = res
        angles[key] = {"theta0": round(float(np.degrees(t0)), 2),
                       "k": round(float(min(k, K_ANGLE_MAX)), 2)}
        report["fitted"][f"angle {key}"] = {**angles[key], "n": len(th)}
    return bonds, angles, report


# ── non-bonded ────────────────────────────────────────────────────────────────

def lj(r: np.ndarray, sigma: float, eps: float, rc: float):
    """Shifted 12-6 LJ: ``u``, ``du/deps``, ``du/dsigma`` on *r* (zero beyond rc)."""
    def raw(rr):
        s6 = (sigma / rr) ** 6
        return 4 * eps * (s6 * s6 - s6), 4 * (s6 * s6 - s6), 4 * eps * (12 * s6 * s6 - 6 * s6) / sigma
    u, de, ds = raw(r)
    u0, de0, ds0 = raw(np.array(rc))
    inside = r < rc
    return (np.where(inside, u - u0, 0.0), np.where(inside, de - de0, 0.0),
            np.where(inside, ds - ds0, 0.0))


def initial_nonbonded(ref_gr: Dict[str, np.ndarray], r: np.ndarray, epsilon: float = 2.0,
                      types: Optional[Sequence[str]] = None) -> Dict[str, dict]:
    """sigma from the first g(r) peak (LJ minimum 2^(1/6) sigma); fixed epsilon."""
    out = {}
    for p, g in ref_gr.items():
        ok = (r > 0.25) & np.isfinite(g)
        peak = r[ok][int(np.argmax(g[ok]))] if ok.any() and np.nanmax(g[ok]) > 0 else 0.6
        out[p] = {"sigma": round(float(peak / 2 ** (1 / 6)), 4), "epsilon": float(epsilon)}
    return out


class Mixing:
    """Map free parameters ``theta`` to per-pair ``(eps, sigma)``.

    ``free``: every pair has its own epsilon and sigma.
    ``lb``  : per-type epsilon/sigma combined with Lorentz-Berthelot rules
              (fewer parameters, and cross terms to other species come for free).
    """

    def __init__(self, pairs: Sequence[str], mode: str = "free"):
        self.pairs, self.mode = list(pairs), mode
        self.types = sorted({t for p in pairs for t in p.split("-")})

    def n(self) -> int:
        return 2 * (len(self.pairs) if self.mode == "free" else len(self.types))

    def theta_from(self, nb: Dict[str, dict]) -> np.ndarray:
        if self.mode == "free":
            return np.array([nb[p][k] for p in self.pairs for k in ("epsilon", "sigma")])
        return np.array([nb[f"{t}-{t}"][k] for t in self.types for k in ("epsilon", "sigma")])

    def pairs_from(self, theta: np.ndarray) -> Dict[str, dict]:
        if self.mode == "free":
            return {p: {"epsilon": float(theta[2 * n]), "sigma": float(theta[2 * n + 1])}
                    for n, p in enumerate(self.pairs)}
        par = {t: (theta[2 * n], theta[2 * n + 1]) for n, t in enumerate(self.types)}
        out = {}
        for p in self.pairs:
            a, b = p.split("-")
            out[p] = {"epsilon": float(np.sqrt(par[a][0] * par[b][0])),
                      "sigma": float(0.5 * (par[a][1] + par[b][1]))}
        return out

    def jacobian(self, theta: np.ndarray) -> np.ndarray:
        """d(pair eps/sigma, flattened per pair) / d theta (numerical)."""
        base = self._flat(theta)
        J = np.zeros((len(base), len(theta)))
        for i in range(len(theta)):
            h = 1e-6 * max(abs(theta[i]), 1e-3)
            t = theta.copy()
            t[i] += h
            J[:, i] = (self._flat(t) - base) / h
        return J

    def _flat(self, theta):
        d = self.pairs_from(theta)
        return np.array([d[p][k] for p in self.pairs for k in ("epsilon", "sigma")])


def _frame_observables(hists: np.ndarray, r: np.ndarray, pairs: Sequence[str],
                       nb: Dict[str, dict], rc: float) -> np.ndarray:
    """(n_frames, 2*n_pairs): per-frame sum over pairs of dU/d(eps, sigma)."""
    X = np.zeros((hists.shape[0], 2 * len(pairs)))
    for n, p in enumerate(pairs):
        _, de, ds = lj(r, nb[p]["sigma"], nb[p]["epsilon"], rc)
        X[:, 2 * n] = hists[:, n, :] @ de
        X[:, 2 * n + 1] = hists[:, n, :] @ ds
    return X


def re_step(ref: Stats, cg: Stats, nb: Dict[str, dict], mixing: Mixing, rc: float,
            temperature: float, trust_eps: float = 0.15, trust_sigma: float = 0.04,
            damping: float = 0.3) -> Tuple[Dict[str, dict], dict]:
    """One damped Newton step of the relative-entropy minimisation."""
    beta = 1.0 / (KB * temperature)
    r = ref.r_centers
    pairs = mixing.pairs
    Xref = _frame_observables(np.array(ref.pair_hist), r, pairs, nb, rc)
    Xcg = _frame_observables(np.array(cg.pair_hist), r, pairs, nb, rc)
    g_pair = beta * (Xref.mean(0) - Xcg.mean(0))
    H_pair = beta ** 2 * np.atleast_2d(np.cov(Xcg, rowvar=False))

    theta = mixing.theta_from(nb)
    J = mixing.jacobian(theta)
    g = J.T @ g_pair
    H = J.T @ H_pair @ J
    # work in ln(theta): d/dln = theta * d/dtheta
    S = np.diag(theta)
    gl, Hl = S @ g, S @ H @ S
    A = Hl + damping * np.diag(np.maximum(np.diag(Hl), 1e-12)) + 1e-12 * np.eye(len(theta))
    try:
        step = -np.linalg.solve(A, gl)
    except np.linalg.LinAlgError:
        step = -gl / np.maximum(np.diag(A), 1e-12)
    limit = np.array([trust_eps, trust_sigma] * (len(theta) // 2))
    step = np.clip(step, -limit, limit)
    new_theta = theta * np.exp(step)
    new = mixing.pairs_from(new_theta)
    for p in new:
        new[p] = {"epsilon": round(new[p]["epsilon"], 5), "sigma": round(new[p]["sigma"], 5)}
    return new, {"grad_norm": float(np.linalg.norm(gl)), "max_step": float(np.abs(step).max())}


def pressure_step(nb: Dict[str, dict], p_cg: float, p_ref: float, gain: float = 0.05,
                  scale: float = 300.0) -> Tuple[Dict[str, dict], float]:
    """Scale every epsilon: pressure too high -> more attraction (and vice versa)."""
    f = float(np.exp(gain * np.tanh((p_cg - p_ref) / scale)))
    return {p: {"epsilon": round(q["epsilon"] * f, 5), "sigma": q["sigma"]}
            for p, q in nb.items()}, f


def rdf_rmse(ref_gr: Dict[str, np.ndarray], cg_gr: Dict[str, np.ndarray], r: np.ndarray,
             weights: Dict[str, float], rmin: float = 0.3) -> Dict[str, float]:
    """RMSE of g(r) per pair above *rmin* and a pair-count-weighted overall value."""
    out, tot, wsum = {}, 0.0, 0.0
    for p in ref_gr:
        m = r > rmin
        e = float(np.sqrt(np.nanmean((ref_gr[p][m] - cg_gr[p][m]) ** 2)))
        out[p] = e
        tot += weights[p] * e
        wsum += weights[p]
    out["overall"] = tot / wsum if wsum else float("nan")
    return out
