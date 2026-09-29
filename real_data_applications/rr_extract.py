"""Reproducible RR extraction from the WFDB source of record 14134
(MIT-BIH Long-Term ECG Database). Reads RealData/14134.atr (beat annotations)
and writes RealData/RR_14134.csv (one column, RR in integer samples).

Only the header (.hea) and annotation (.atr) files are needed; the signal (.dat)
is not read. Requires the wfdb package. This regenerates the exact series used by
rr_analysis.py.

Run:  python3 rr_extract.py
"""
import numpy as np, pandas as pd, wfdb
from pathlib import Path

REC = Path(__file__).resolve().parent / "RealData" / "14134"
BEAT_SYMBOLS = list("NLReAaJSVEFjnf/Q")   # standard MIT beat labels (excludes rhythm/noise marks)

def main():
    hdr = wfdb.rdheader(str(REC)); fs = hdr.fs
    ann = wfdb.rdann(str(REC), "atr")
    samp = np.asarray(ann.sample); sym = np.asarray(ann.symbol)
    beats = np.isin(sym, BEAT_SYMBOLS)
    rr = np.diff(samp[beats]).astype(int)          # RR in integer samples (exact quantization)
    pd.DataFrame({"RR_samples": rr}).to_csv(REC.parent / "RR_14134.csv", index=False)
    print(f"fs={fs}  beats={int(beats.sum())}  RR intervals={len(rr)}  "
          f"median={int(np.median(rr))} samples ({60*fs/np.median(rr):.1f} bpm)  "
          f"ms lattice step={1000/fs:.4f}")
    print("beat labels present:", sorted(set(sym[beats].tolist())))

if __name__ == "__main__":
    main()
