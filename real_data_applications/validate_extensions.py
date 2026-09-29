"""Numerical validation tables for the two extensions.

(A) tau>1 white-noise bias of the disequilibrium: exact closed form vs Monte Carlo,
    for the plug-in and distinct-pairs (U) estimators, plus the gap estimator with
    w=(m-1)tau (should be ~0). For continuous white noise the true D is 0, so the
    plug-in/U means ARE the bias.

(B) AR(1) reference disequilibrium for m=3,4,5: the orthant-probability value vs a
    direct large-sample simulation of the AR(1) process (independent check). The
    orthant reduces to the closed-form arcsine expressions for m<=4.

Requires apply_diagnostics.py, fsb_extensions.py.
"""
import numpy as np
from math import factorial
from scipy.signal import lfilter
import apply_diagnostics as A
import fsb_extensions as E


# ---------- A. tau>1 bias: exact vs Monte Carlo ----------
def mc_tau_bias(m, tau, N, B, seed=0):
    M = factorial(m); w = (m - 1) * tau
    rng = np.random.default_rng(seed)
    Dp = np.empty(B); Du = np.empty(B); Dg = np.empty(B)
    for b in range(B):
        s = A.symbols(rng.standard_normal(N), m, tau, policy="index")
        co = A.coordinates(s, M)
        Dp[b] = co["D"]; Du[b] = co["D_U"]; Dg[b] = A.gap_D(s, M, w)
    return Dp, Du, Dg


print("A. tau>1 white-noise bias of D   (N=2000, D_true=0, B=10000 reps)")
print(f"{'m':>2} {'tau':>4} {'n':>5} | {'plug pred':>11} {'plug MC':>11} | "
      f"{'U pred':>11} {'U MC':>11} | {'gap MC':>11}")
N, B = 2000, 10000
for m in (3, 4):
    for tau in (2, 3):
        n = N - (m - 1) * tau
        bp, bu = E.tau_biases(m, tau, n)                      # exact (Fraction)
        Dp, Du, Dg = mc_tau_bias(m, tau, N, B, seed=11)
        se = lambda a: a.std() / np.sqrt(len(a))
        print(f"{m:>2} {tau:>4} {n:>5} | {float(bp):>11.3e} {Dp.mean():>11.3e} | "
              f"{float(bu):>11.3e} {Du.mean():>11.3e} | {Dg.mean():>11.3e}")
        print(f"{'':>2} {'(se)':>4} {'':>5} | {'':>11} {se(Dp):>11.1e} | "
              f"{'':>11} {se(Du):>11.1e} | {se(Dg):>11.1e}")


# ---------- B. AR(1) reference D: orthant vs simulation ----------
def ar1(N, phi, seed):
    rng = np.random.default_rng(seed)
    x = lfilter([1.0], [1.0, -phi], rng.standard_normal(N + 1000))  # burn-in 1000
    return x[1000:]

def sim_D(m, phi, N, seed=3):
    M = factorial(m)
    s = A.symbols(ar1(N, phi, seed), m, 1, policy="index")
    p = np.bincount(s, minlength=M).astype(float); p /= p.sum()
    return ((p - 1 / M) ** 2).sum(), p

print("\nB. AR(1) reference D  (orthant vs simulation, N=5e6, tau=1)")
print(f"{'m':>2} {'phi':>4} | {'D orthant':>12} {'D sim':>12} {'|dD|':>9} {'max|dp|':>9}")
Nsim = 5_000_000
for m in (3, 4, 5):
    for phi in (0.5, 0.9):
        p_o = E.ar_probabilities_orthant(m, phi)
        D_o = E.functionals(p_o)[1]
        D_s, p_s = sim_D(m, phi, Nsim)
        print(f"{m:>2} {phi:>4} | {D_o:>12.6e} {D_s:>12.6e} {abs(D_o-D_s):>9.1e} "
              f"{np.abs(p_o-p_s).max():>9.1e}")
