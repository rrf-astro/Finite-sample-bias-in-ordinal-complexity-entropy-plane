"""Real short-series application: two Forbush decreases, Rome and Oulu neutron
monitors, integer counts recovered from NMDB rates.

Events:
  FD240510  Gannon storm, 2024-05-10/11  (deeper depression)
  FD240324  storm of 2024-03-24/25       (shallower depression)

Produces figure_forbush.(png|pdf) and forbush_summary.csv. Reads the NMDB export
files in RealData/. Requires apply_diagnostics.py in the same folder.

Run:  python3 forbush_analysis.py
"""
import numpy as np, pandas as pd, json
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
DATA = HERE / "RealData"
m, M = 4, factorial(4)
EVENTS = {"FD240510": "Gannon 2024-05-10", "FD240324": "storm 2024-03-24"}
STATIONS = ("ROME", "OULU")
TAUS = [1, 2, 3, 5, 8, 10, 12, 20, 30, 45, 60, 90, 120]
RANK_TAUS = (1, 10, 30, 60)
B_NULL = 9999

def load_nmdb(path):
    """Robust NMDB reader: infers station columns from the STATION header."""
    stations, rows = None, []
    for ln in open(path):
        if ln.startswith("#"):
            if "STATION:" in ln:
                stations = [s.strip() for s in ln.split("STATION:")[1].split(",")]
            continue
        if ln.strip() == "" or ";" not in ln:
            continue
        p = ln.strip().split(";")
        if stations and len(p) == len(stations) + 1:
            rows.append([p[0]] + [float(v) for v in p[1:]])
    df = pd.DataFrame(rows, columns=["dt"] + stations)
    return df.set_index(pd.to_datetime(df.pop("dt"))), stations

def integer_counts(df, st):
    """Rome and Oulu uncorrected rates are exactly (counts per minute)/60."""
    c = np.round(df[st].values * 60)
    assert np.abs(df[st].values * 60 - c).max() < 0.03, f"{st}: no clean integer scale"
    return c.astype(int)

def analyze():
    out = {}
    for ev in EVENTS:
        df, _ = load_nmdb(DATA / f"{ev}_Uncorrected.txt")
        rec = {"index": df.index, "series": {}}
        for st in STATIONS:
            x = integer_counts(df, st)
            curves = {k: [] for k in ("H", "Dplug", "DU", "Dgap", "CJS")}
            for tau in TAUS:
                s = A.symbols(x, m, tau, policy="index"); co = A.coordinates(s, M)
                curves["H"].append(co["H"]); curves["Dplug"].append(co["D"])
                curves["DU"].append(co["D_U"]); curves["Dgap"].append(A.gap_D(s, M, (m-1)*tau))
                curves["CJS"].append(co["C_JS"])
            ranks = {t: A.rank_test(x, m, t, policy="index", B_null=B_NULL, seed=7)[1] for t in RANK_TAUS}
            base = x[:120].mean()
            rec["series"][st] = dict(counts=x, curves=curves, ranks=ranks,
                                     depth=100*(base - x.min())/base)
        out[ev] = rec
    # controlled quantization on Rome of the Gannon event (tie demonstration)
    xg = out["FD240510"]["series"]["ROME"]["counts"]
    H1 = A.coordinates(A.symbols(xg, m, 1, policy="index"), M)["H"]
    D1 = A.coordinates(A.symbols(xg, m, 1, policy="index"), M)["D_U"]
    quant = {"q": [], "tie": [], "dH": [], "dDU": []}
    for q in [1, 2, 5, 10, 20, 40, 80, 160]:
        xq = np.round(xg / q) * q; s = A.symbols(xq, m, 1, policy="index"); co = A.coordinates(s, M)
        quant["q"].append(q); quant["tie"].append(A.tie_fraction(xq, m, 1))
        quant["dH"].append(co["H"] - H1); quant["dDU"].append(co["D_U"] - D1)
    return out, quant

def write_summary(out):
    rows = []
    for ev in EVENTS:
        for st in STATIONS:
            s = out[ev]["series"][st]
            for t in RANK_TAUS:
                i = TAUS.index(t) if t in TAUS else None
                H = s["curves"]["H"][i] if i is not None else np.nan
                Dg = s["curves"]["Dgap"][i] if i is not None else np.nan
                rows.append(dict(event=ev, station=st, depth_pct=round(s["depth"], 1),
                                 tau=t, H=round(H, 4) if i is not None else "",
                                 D_gap=round(Dg, 5) if i is not None else "",
                                 rank_p=round(s["ranks"][t], 4)))
    pd.DataFrame(rows).to_csv(HERE / "forbush_summary.csv", index=False)

def figure(out, quant, suffix=""):
    C = {"ROME": WONG[0], "OULU": WONG[1]}
    LS = {"FD240510": "-", "FD240324": "--"}
    MK = {"ROME": "o", "OULU": "s"}
    _apply_style()
    fig, ax = plt.subplots(2, 2, figsize=(7.4, 5.6), layout="constrained")
    # (a) OULU depression, both events, aligned by elapsed hours
    axa = ax[0, 0]
    for ev in EVENTS:
        idx = out[ev]["index"]; hrs = (idx - idx[0]).total_seconds() / 3600
        v = out[ev]["series"]["OULU"]["counts"].astype(float); pct = 100 * (v / v[:120].mean() - 1)
        axa.plot(hrs, pd.Series(pct).rolling(15, center=True).mean(), color=C["OULU"], ls=LS[ev], lw=1.35,
                 label=f"{EVENTS[ev]} (OULU, {out[ev]['series']['OULU']['depth']:.0f}%)")
    axa.axhline(0, color="0.4", lw=0.7); axa.set_xlabel("hours from event start")
    axa.set_ylabel("count deviation (%)"); axa.set_title("(a) Two Forbush decreases (Oulu)", loc="left")
    axa.legend(loc="lower left")
    # (b) H vs tau, four series, rank significance marked
    axb = ax[0, 1]
    for ev in EVENTS:
        for st in STATIONS:
            s = out[ev]["series"][st]
            axb.plot(TAUS, s["curves"]["H"], color=C[st], ls=LS[ev], lw=1.2, marker=MK[st], ms=3.2,
                     label=f"{st} {'Gan' if ev=='FD240510' else 'Mar'}")
            for t in RANK_TAUS:
                if t in TAUS:
                    p = s["ranks"][t]; Hi = s["curves"]["H"][TAUS.index(t)]
                    axb.plot(t, Hi, marker=MK[st], ms=7, mfc=(C[st] if p < 0.05 else "white"),
                             mec=C[st], mew=1.1, zorder=5)
    axb.set_xscale("log"); axb.set_xlabel(r"embedding delay $\tau$ (min)")
    axb.set_ylabel(r"normalized permutation entropy $H$")
    axb.set_title(r"(b) Emergence with $\tau$ (filled: $p<0.05$)", loc="left")
    axb.legend(loc="lower left", ncol=2)
    # (c) gap-corrected D vs tau, four series
    axc = ax[1, 0]
    for ev in EVENTS:
        for st in STATIONS:
            axc.plot(TAUS, out[ev]["series"][st]["curves"]["Dgap"], color=C[st], ls=LS[ev], lw=1.2,
                     marker=MK[st], ms=3.2, label=f"{st} {'Gan' if ev=='FD240510' else 'Mar'}")
    axc.axhline(0, color="0.4", lw=0.7); axc.set_xscale("log")
    axc.set_xlabel(r"embedding delay $\tau$ (min)"); axc.set_ylabel(r"gap-corrected disequilibrium $D$")
    axc.set_title("(c) Overlap-corrected structure vs event strength", loc="left")
    axc.legend(loc="upper left", ncol=2)
    # (d) controlled quantization (Gannon, Rome)
    axd = ax[1, 1]; tie = 100 * np.array(quant["tie"])
    axd.plot(tie, quant["dH"], color=WONG[0], marker="o", ms=4)
    axd.set_xlabel("tied ordinal windows (%)"); axd.set_ylabel(r"$\Delta H$", color=WONG[0])
    axd.tick_params(axis="y", labelcolor=WONG[0]); axd.axhline(0, color="0.4", lw=0.7)
    axd2 = axd.twinx(); axd2.plot(tie, quant["dDU"], color=WONG[1], ls="--", marker="s", ms=4)
    axd2.set_ylabel(r"$\Delta D_U$", color=WONG[1]); axd2.tick_params(axis="y", labelcolor=WONG[1])
    axd.set_title("(d) Controlled quantization (Rome, Gannon)", loc="left")
    for a in (axa, axb, axc):
        sm(a)
    axd.minorticks_on(); axd.tick_params(which="both", direction="in", top=True)
    fig.savefig(HERE / f"figure_forbush{suffix}.png", dpi=400, bbox_inches="tight", pad_inches=0.05)
    fig.savefig(HERE / f"figure_forbush{suffix}.pdf", bbox_inches="tight", pad_inches=0.05)

def main():
    out, quant = analyze()
    for ev in EVENTS:
        for st in STATIONS:
            s = out[ev]["series"][st]
            print(f"{ev} {st}: depth {s['depth']:.1f}%  rank p", {t: round(s["ranks"][t], 4) for t in RANK_TAUS})
    write_summary(out); figure(out, quant)
    print("saved figure_forbush.png/pdf, forbush_summary.csv")

if __name__ == "__main__":
    main()
