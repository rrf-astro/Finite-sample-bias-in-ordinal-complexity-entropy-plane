"""Real physiological short series: RR intervals from an ECG record (MIT-BIH
Long-Term ECG Database, record 14134, fs = 128 Hz). RR is measured in integer
samples (the exact physical quantization); in milliseconds it lives on a lattice
of 1000/128 = 7.8125 ms, so ties are frequent. This complements the Forbush case:
here the series carries genuine cardiac structure at unit delay, and the tie
policy, not the overlap correction, is what moves the plane.

Produces figure_rr.(png|pdf) and rr_summary.csv. Reads RealData/RR_14134.csv
(one column RR_samples). Requires apply_diagnostics.py.

Run:  python3 rr_analysis.py
"""
import numpy as np, pandas as pd
from math import factorial
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import apply_diagnostics as A

# Match the manuscript figure style (SuperMongo-inspired, from the main notebook).
WONG = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00"]

def _apply_style():
    plt.rcParams.update({
        "font.family": "serif", "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
        "mathtext.fontset": "stix", "font.size": 9, "axes.labelsize": 10, "axes.titlesize": 9,
        "axes.linewidth": 1.0, "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True, "xtick.major.size": 5, "ytick.major.size": 5,
        "xtick.minor.size": 2.5, "ytick.minor.size": 2.5, "xtick.major.width": 0.9,
        "ytick.major.width": 0.9, "xtick.minor.width": 0.7, "ytick.minor.width": 0.7,
        "legend.frameon": True, "legend.framealpha": 0.90, "legend.facecolor": "white",
        "legend.edgecolor": "0.75", "legend.fontsize": 7.2, "lines.linewidth": 1.35,
        "lines.markersize": 4.2, "pdf.fonttype": 42, "ps.fonttype": 42})

def sm(ax):
    ax.minorticks_on(); ax.tick_params(which="both", direction="in", top=True, right=True)
    for s in ax.spines.values():
        s.set_visible(True); s.set_linewidth(1.0)

HERE = Path(__file__).resolve().parent
FS = 128.0
MS = 1000.0 / FS               # 7.8125 ms per sample
m, M = 4, factorial(4)
REP_START, REP_N = 2000, 2000  # representative window (no gross artifacts)
WIN = 1500                     # window for the across-record tie/policy scan

def clean(x):
    """Winsorize gross non-physiological RR (detector gaps / missed beats) to the
    window median, preserving sequence contiguity. Threshold 40..250 samples
    (about 31 bpm to 192 bpm)."""
    x = x.astype(float).copy(); med = np.median(x)
    x[(x > 250) | (x < 40)] = med
    return x

def load():
    return pd.read_csv(HERE / "RealData" / "RR_14134.csv")["RR_samples"].values.astype(float)

def analyze():
    rr = load()
    rep = clean(rr[REP_START:REP_START + REP_N])
    out = {"rep": rep, "rep_ms": rep * MS}
    # representative window: estimators and policies
    res = {}
    for pol in ("index", "jitter"):
        s = A.symbols(rep, m, 1, policy=pol, rng=np.random.default_rng(0)); co = A.coordinates(s, M)
        res[pol] = dict(H=co["H"], D=co["D"], D_U=co["D_U"], D_gap=A.gap_D(s, M, 3), C_JS=co["C_JS"])
    T, p = A.rank_test(rep, m, 1, policy="index", B_null=9999, seed=3)
    res["tie"] = A.tie_fraction(rep, m, 1); res["rank_p"] = p
    out["rep_res"] = res
    # across-record scan: tie fraction and index-vs-jitter policy gap
    scan = []
    for st in range(0, len(rr) - WIN, WIN):
        x = clean(rr[st:st + WIN]); tf = A.tie_fraction(x, m, 1)
        si = A.coordinates(A.symbols(x, m, 1, policy="index"), M)
        sj = A.coordinates(A.symbols(x, m, 1, policy="jitter", rng=np.random.default_rng(0)), M)
        scan.append((st, tf, si["H"], sj["H"], si["C_JS"], sj["C_JS"]))
    out["scan"] = np.array(scan)
    return out

def figure(out, suffix=""):
    rep_ms = out["rep_ms"]; scan = out["scan"]; res = out["rep_res"]
    tie = 100 * scan[:, 1]; dH = scan[:, 2] - scan[:, 3]          # index - jitter (signed)
    _apply_style()
    fig, ax = plt.subplots(2, 2, figsize=(7.4, 5.6), layout="constrained")
    # (a) RR segment (ms)
    axa = ax[0, 0]; t = np.arange(len(rep_ms))
    axa.plot(t, rep_ms, color=WONG[0], lw=0.7)
    axa.set_xlabel("beat number"); axa.set_ylabel("RR interval (ms)")
    axa.set_title("(a) RR intervals, record 14134 (128 Hz)", loc="left")
    axa.text(0.97, 0.05, f"tie fraction {100*res['tie']:.0f}%  ·  rank $p<10^{{-3}}$",
             transform=axa.transAxes, ha="right", fontsize=7.2, color="0.35")
    # (b) tie fraction across the record
    axb = ax[0, 1]; axb.plot(scan[:, 0] / 1000, tie, color=WONG[0], marker="o", ms=3, lw=0.9)
    axb.axhline(25, color="0.5", lw=0.7, ls=":")
    axb.set_xlabel("record position (thousands of beats)"); axb.set_ylabel("tied ordinal windows (%)")
    axb.set_title("(b) Ties vary with heart-rate stability", loc="left")
    # (c) tie policy sensitivity vs tie fraction
    axc = ax[1, 0]; axc.scatter(tie, dH, s=16, color=WONG[1], edgecolor="k", lw=0.3, zorder=3)
    bslope, a0 = np.polyfit(tie, dH, 1); xx = np.linspace(tie.min(), tie.max(), 50)
    axc.plot(xx, a0 + bslope * xx, "k--", lw=0.9)
    axc.axhline(0, color="0.4", lw=0.7)
    axc.set_xlabel("tied ordinal windows (%)"); axc.set_ylabel(r"$H_{\mathrm{index}} - H_{\mathrm{jitter}}$")
    axc.set_title(r"(c) Tie policy shifts $H$ (corr $=%.2f$)" % np.corrcoef(tie, dH)[0, 1], loc="left")
    # (d) estimator comparison for the representative window
    axd = ax[1, 1]
    vals = [res["index"]["D"], res["index"]["D_U"], res["index"]["D_gap"]]
    axd.bar(range(3), vals, color=[WONG[0], WONG[1], WONG[2]], width=0.6, edgecolor="k", lw=0.6)
    axd.set_xticks(range(3)); axd.set_xticklabels(["plug-in", "distinct-\npairs $U$", "gap"])
    axd.set_ylabel(r"disequilibrium $D$")
    axd.set_title("(d) Overlap correction negligible for strong structure", loc="left")
    axd.set_ylim(0, max(vals) * 1.25)
    for i, v in enumerate(vals):
        axd.text(i, v, f"{v:.4f}", ha="center", va="bottom", fontsize=7.2)
    for a in (axa, axb, axc):
        sm(a)
    axd.tick_params(which="both", direction="in", top=True)
    fig.savefig(HERE / f"figure_rr{suffix}.png", dpi=400, bbox_inches="tight", pad_inches=0.05)
    fig.savefig(HERE / f"figure_rr{suffix}.pdf", bbox_inches="tight", pad_inches=0.05)

def main():
    out = analyze(); res = out["rep_res"]
    print(f"representative window [{REP_START}:{REP_START+REP_N}]  tie={100*res['tie']:.1f}%  rank p={res['rank_p']:.4f}")
    for pol in ("index", "jitter"):
        r = res[pol]; print(f"  {pol:7s} H={r['H']:.4f} D plug={r['D']:.5f} U={r['D_U']:+.5f} gap={r['D_gap']:+.5f} C_JS={r['C_JS']:.5f}")
    scan = out["scan"]; tie = 100 * scan[:, 1]
    print(f"across record: tie {tie.min():.0f}-{tie.max():.0f}%  corr(tie,|dH|)={np.corrcoef(tie, np.abs(scan[:,2]-scan[:,3]))[0,1]:.3f}")
    pd.DataFrame(dict(start=scan[:, 0].astype(int), tie_frac=scan[:, 1],
                      H_index=scan[:, 2], H_jitter=scan[:, 3],
                      CJS_index=scan[:, 4], CJS_jitter=scan[:, 5])).to_csv(HERE / "rr_summary.csv", index=False)
    figure(out); print("saved figure_rr.png/pdf, rr_summary.csv")

if __name__ == "__main__":
    main()
