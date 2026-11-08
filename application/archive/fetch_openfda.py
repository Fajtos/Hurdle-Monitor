"""Download daily report counts from the openFDA drug adverse-event API (FAERS), one CSV per panel row.
usage: python3 fetch_openfda.py panel.csv outdir [--start 2012] [--end 2024]
panel.csv columns: id,drug,reaction,role,note   (reaction may be empty = all reports for the drug; role = target|control)
Needs internet and a free openFDA API key (https://open.fda.gov/apis/authentication/) in the environment variable OPENFDA_API_KEY (never hard-code or share it).
NOT run in the authoring sandbox (outbound access blocked): test on one row first and inspect the CSV.
Query form (documented openFDA syntax): /drug/event.json?search=patient.drug.openfda.generic_name:"X"+AND+patient.reaction.reactionmeddrapt:"Y"+AND+receivedate:[YYYYMMDD+TO+YYYYMMDD]&count=receivedate&limit=1000"""
import sys, os, csv, json, time, urllib.request, urllib.error, urllib.parse, argparse
import pandas as pd
BASE = "https://api.fda.gov/drug/event.json"
def q(s): return urllib.parse.quote(s, safe="")
def url(drug, reaction, d0, d1, key):
    Q = "%22"                                                      # percent-encoded double quote
    parts = [f'patient.drug.openfda.generic_name:{Q}{q(drug)}{Q}']
    if reaction: parts.append(f'patient.reaction.reactionmeddrapt:{Q}{q(reaction)}{Q}')
    parts.append(f"receivedate:%5B{d0}+TO+{d1}%5D")
    u = f"{BASE}?search={'+AND+'.join(parts)}&count=receivedate&limit=1000"
    return u + (f"&api_key={key}" if key else "")
UAS = ["hurdle-monitor-research", "Mozilla/5.0 (X11; Linux x86_64) research-script", None]   # retried in turn on a 403
def get(u, tries=4):
    ua_i = 0
    for k in range(tries + len(UAS)):
        try:
            h = {"User-Agent": UAS[ua_i]} if UAS[ua_i] else {}
            with urllib.request.urlopen(urllib.request.Request(u, headers=h), timeout=60) as r: return json.load(r)
        except urllib.error.HTTPError as e:
            body = e.read().decode("utf8", "replace")[:300]
            if e.code == 404: return {"results": []}          # openFDA returns 404 when nothing matches
            if e.code == 429: time.sleep(5 * (k + 1)); continue
            if e.code == 403 and ua_i < len(UAS) - 1: ua_i += 1; continue
            print("HTTP", e.code, "for", u, "\nresponse body:", body); raise
        except Exception: time.sleep(2 * (k + 1))
    raise RuntimeError("request failed: " + u)
def main():
    ap = argparse.ArgumentParser(); ap.add_argument("panel"); ap.add_argument("outdir"); ap.add_argument("--start", type=int, default=2012); ap.add_argument("--end", type=int, default=2024); a = ap.parse_args()
    key = os.environ.get("OPENFDA_API_KEY")
    if not key: sys.exit("openFDA now requires an API key. Get a free one at https://open.fda.gov/apis/authentication/ and set it: in Colab run  %env OPENFDA_API_KEY=your_key  ; in a terminal  export OPENFDA_API_KEY=your_key")
    os.makedirs(a.outdir, exist_ok=True)
    for row in csv.DictReader(open(a.panel)):
        recs = []
        for y0 in range(a.start, a.end + 1, 2):                    # <=1000 distinct days per call
            y1 = min(y0 + 1, a.end); js = get(url(row["drug"], row.get("reaction", ""), f"{y0}0101", f"{y1}1231", key)); recs += [(r["time"], r["count"]) for r in js.get("results", [])]; time.sleep(0.3)
        idx = pd.date_range(f"{a.start}-01-01", f"{a.end}-12-31"); s = pd.Series(0, index=idx)
        for t, c in recs: s[pd.Timestamp(t)] = c
        out = os.path.join(a.outdir, f"{row['id']}.csv"); pd.DataFrame({"date": idx, "count": s.values}).to_csv(out, index=False); print(row["id"], "days", len(idx), "zero share", round(float((s.values == 0).mean()), 3), "->", out)
if __name__ == "__main__": main()
