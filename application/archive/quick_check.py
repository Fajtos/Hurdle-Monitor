"""Feasibility check: fetch candidate openFDA series and report what matters for a hurdle study.
Colab: %env OPENFDA_API_KEY=...  then  !python quick_check.py        (edit PAIRS / TRAIN_END to try others)"""
import subprocess, sys, os, pandas as pd, numpy as np
PAIRS = [("r1", "valsartan", ""), ("r2", "metformin", "LACTIC ACIDOSIS"), ("r3", "lisinopril", "COUGH"), ("r4", "atorvastatin", "MYALGIA"),
         ("r5", "semaglutide", ""), ("r6", "isotretinoin", "DEPRESSION"), ("r7", "ranitidine", "PROSTATE CANCER"), ("r8", "levothyroxine", "ALOPECIA")]
START, END, TRAIN_END = 2012, 2024, "2017-12-31"
os.makedirs("qdata", exist_ok=True)
with open("qpanel.csv", "w") as f:
    f.write("id,drug,reaction,role,note\n"); [f.write(f"{i},{d},{r},target,check\n") for i, d, r in PAIRS]
subprocess.run([sys.executable, "fetch_openfda.py", "qpanel.csv", "qdata", "--start", str(START), "--end", str(END)], check=True)
print("\nRule of thumb: keep a series if the WEEKDAY zero share in training is in [0.15,0.85] and positive days are not all 1 (mean positive count >= 1.3)\n")
for i, d, r in PAIRS:
    p = f"qdata/{i}.csv"
    if not os.path.exists(p): print(i, "no file"); continue
    s = pd.read_csv(p, parse_dates=["date"]); tr = s[s["date"] <= TRAIN_END]; wd = lambda x: x[x["date"].dt.dayofweek < 5]["count"].to_numpy(); we = tr[tr["date"].dt.dayofweek >= 5]["count"].to_numpy()
    w = wd(tr); pos = w[w > 0]; ok = 0.15 <= np.mean(w == 0) <= 0.85 and len(pos) > 30 and pos.mean() >= 1.3
    print(f"{i} {d} {r}: total reports={int(s['count'].sum())}  training weekday zero share={np.mean(w==0):.2f}  weekend mean={we.mean():.2f} vs weekday mean={w.mean():.2f}")
    print(f"   positive weekdays={len(pos)}  mean positive count={pos.mean() if len(pos) else float('nan'):.2f}  share of positives >1={np.mean(pos>1) if len(pos) else float('nan'):.2f}  max={int(s['count'].max())}")
    print("   yearly mean:", s.groupby(s['date'].dt.year)['count'].mean().round(2).to_dict(), "\n   USABLE:", ok)
