"""Step 6 — FOLH1-locus methylation and the coordinate/annotation audit for Figure 5.

PART A (replication)
  A1  36 locus probes vs FOLH1 expression in the 497 paired TCGA tumours
PART B (independent audit of the sequence-level claims)
  B1  FOLH1 GRCh38 coordinates and strand (Ensembl) vs the manuscript
  B2  are the three scored "CpG-disrupting" positions actually CpG cytosines?
  B3  the nominated gene-body element: position relative to the gene and the TSS
  B4  ChIP-Atlas: is FOLH1 a target of AR / FOXA1 / HOXB13 / EZH2?
"""
import json
import numpy as np
import pandas as pd
import requests
from scipy import stats as st
from statsmodels.stats.multitest import multipletests

from common import DATA, OUT, load_tcga

res = {}
HDR = {"Content-Type": "application/json"}

# ---------------------------------------------------------------- A1
tcga_log, _ = load_tcga()
meth = pd.read_csv(f"{DATA}/folh1_locus_meth.csv", index_col=0)
ann = pd.read_csv(f"{DATA}/folh1_locus_meth_annot.csv", index_col=0)
y = tcga_log.loc["FOLH1"]
y.index = [c[:12] for c in y.index]
y = y[~y.index.duplicated(keep="first")]
mt = [c for c in meth.columns if c.endswith("-01")]
m = meth[mt].copy()
m.columns = [c[:12] for c in mt]
m = m.loc[:, ~m.columns.duplicated()]
common = y.index.intersection(m.columns)
print(f"[A1] paired tumours: {len(common)}   [manuscript 497]")
rows = []
for p in m.index:
    b = m.loc[p, common]
    ok = b.notna()
    if ok.sum() < 50:
        continue
    r, pv = st.spearmanr(b[ok], y[common][ok])
    rows.append(dict(probe=p, pos_hg19=int(ann.loc[p, "chromStart"]), gene=ann.loc[p, "gene"],
                     rho=r, p=pv, n=int(ok.sum())))
mr = pd.DataFrame(rows).sort_values("rho")
mr["padj"] = multipletests(mr["p"], method="fdr_bh")[1]
print(mr.round(4).to_string(index=False))
res["A1"] = mr.round(5).to_dict("records")
mr.to_csv(f"{OUT}/s06_meth_probes.csv", index=False)
for p, exp in [("cg06980460", -0.287), ("cg10006235", -0.238), ("cg25346915", -0.123)]:
    got = float(mr.loc[mr.probe == p, "rho"].iloc[0])
    print(f"[A1] {p}: rho={got:+.3f}  [manuscript {exp:+.3f}]")

# ---------------------------------------------------------------- B1
g = requests.get("https://rest.ensembl.org/lookup/symbol/homo_sapiens/FOLH1",
                 headers=HDR, timeout=90).json()
print(f"\n[B1] Ensembl GRCh38 FOLH1: chr{g['seq_region_name']}:{g['start']:,}-{g['end']:,} "
      f"strand {g['strand']}   [manuscript chr11:49,145,092-49,208,641, minus strand]")
res["B1"] = dict(start=g["start"], end=g["end"], strand=g["strand"],
                 matches_manuscript=(g["start"] == 49145092 and g["end"] == 49208641
                                     and g["strand"] == -1))
TSS = g["end"]

# ---------------------------------------------------------------- B2
print("\n[B2] reference sequence at the three scored CpG-disruption positions")
scored = {"cg06980460": 49208996, "cg25346915": 49208484, "cg10006235": 49204818}
lifted = {"cg06980460": 49208995, "cg25346915": 49208483, "cg10006235": 49204817}
b2 = {}
for probe, pos in scored.items():
    s0, s1 = pos - 6, pos + 6
    seq = requests.get(f"https://rest.ensembl.org/sequence/region/human/11:{s0}-{s1}"
                       f"?coord_system_version=GRCh38",
                       headers={"Content-Type": "text/plain"}, timeout=90).text.strip().upper()
    centre = seq[pos - s0]
    nxt = seq[pos - s0 + 1]
    prv = seq[pos - s0 - 1]
    lift = lifted[probe]
    centre_l = seq[lift - s0]
    nxt_l = seq[lift - s0 + 1]
    is_cpg_scored = (centre == "C" and nxt == "G") or (centre == "G" and prv == "C")
    is_cpg_lifted = (centre_l == "C" and nxt_l == "G")
    print(f"   {probe}: window {s0}-{s1} = {seq}")
    print(f"      scored position {pos}: base '{centre}' (next '{nxt}') -> "
          f"CpG context: {is_cpg_scored};  C>T as modelled: {centre=='C'}")
    print(f"      chain-lifted position {lift}: base '{centre_l}{nxt_l}' -> CpG: {is_cpg_lifted}")
    print(f"      distance from TSS ({TSS:,}) = {TSS - lift:+,} bp "
          f"({'upstream' if lift > TSS else 'inside gene'} on the minus strand)")
    b2[probe] = dict(seq=seq, scored_pos=pos, scored_base=centre, next_base=nxt,
                     cpg_at_scored=bool(is_cpg_scored), lifted_pos=lift,
                     lifted_dinuc=centre_l + nxt_l, cpg_at_lifted=bool(is_cpg_lifted),
                     dist_from_tss=int(TSS - lift))
res["B2"] = b2

# ---------------------------------------------------------------- B3
ENH_S, ENH_E = 49177981, 49178322
print(f"\n[B3] nominated gene-body element chr11:{ENH_S:,}-{ENH_E:,}")
print(f"     inside the gene body: {g['start'] <= ENH_S and ENH_E <= g['end']}")
print(f"     distance from the TSS: {(TSS-ENH_E)/1000:.1f} kb (manuscript figure labels it "
      f"a gene-body DNase element)")
core_s, core_e = 49178125, 49178193
print(f"     'sensitive core' chr11:{core_s:,}-{core_e:,} = {core_e-core_s} bp, "
      f"inside the element: {ENH_S <= core_s and core_e <= ENH_E}")
res["B3"] = dict(in_gene_body=bool(g["start"] <= ENH_S and ENH_E <= g["end"]),
                 kb_from_tss=float((TSS - ENH_E) / 1000))

# ---------------------------------------------------------------- B4
print("\n[B4] ChIP-Atlas precomputed target genes (hg38, +/-10 kb TSS)")
b4 = {}
for tf in ["AR", "FOXA1", "HOXB13", "EZH2"]:
    try:
        r = requests.get(f"https://chip-atlas.dbcls.jp/data/hg38/target/{tf}.10.tsv", timeout=120)
        if r.status_code == 200:
            hits = [l for l in r.text.split("\n") if l.startswith("FOLH1\t") or "\tFOLH1\t" in l
                    or l.split("\t")[0] == "FOLH1"]
            n = r.text.count("\n")
            print(f"   {tf}: {n} target genes; FOLH1 rows {len(hits)}")
            if hits:
                print(f"      {hits[0][:160]}")
            b4[tf] = dict(n_targets=n, folh1_rows=len(hits), first=hits[0][:200] if hits else None)
        else:
            print(f"   {tf}: HTTP {r.status_code}")
            b4[tf] = dict(error=f"HTTP {r.status_code}")
    except Exception as e:
        print(f"   {tf}: FAILED {type(e).__name__}")
        b4[tf] = dict(error=type(e).__name__)
res["B4"] = b4

json.dump(res, open(f"{OUT}/s06_results.json", "w"), indent=1)
print("\nsaved s06 outputs")
