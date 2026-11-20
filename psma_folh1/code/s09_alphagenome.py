"""Step 9 — AlphaGenome replication (manuscript Results section 6 / Figure 5).

Run with a stage argument:  landscape | variants | ism | ism_control | allele | window | contact | loci

PART A (replication, same recipe as the authors' notebook, AlphaGenome v1 model)
  A1  prostate DNase landscape over the FOLH1 locus and the nominated regulatory peaks
  A2  variant effects: peak SNVs + the three CpG-disrupting variants
  A3  4-bp in-silico mutagenesis scan across the nominated gene-body element
  A4  contact maps (H1-hESC) and the four nominated target loci
PART B (independent validation)
  B1  ALLELE SPECIFICITY: is the +0.29 at cg25346915 really CpG disruption, or a position effect?
  B2  BACKGROUND ISM: the same 4-bp scan over a matched control region with no DNase peak
  B3  WINDOW SENSITIVITY: does the effect size depend on the prediction window?
  B4  DETERMINISM: re-scoring the same variant
"""
import json
import os
import sys
import time

import numpy as np
import pandas as pd
import requests

KEY_PATH = os.environ.get("ALPHAGENOME_KEY_FILE",
                          "/tmp/claude-0/-home-user-cloude/a7bfcead-6f81-5243-ba6f-8839a5b9ecee"
                          "/scratchpad/.alphagenome_key")
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "out")
CACHE = os.path.join(OUT, "ag_cache")
os.makedirs(CACHE, exist_ok=True)

from alphagenome.models.v1 import dna_client
from alphagenome.data import genome

MODEL = dna_client.create(api_key=open(KEY_PATH).read().strip())

# manuscript / Ensembl GRCh38 coordinates
GENE_S, GENE_E = 49_145_092, 49_208_641          # FOLH1, minus strand; TSS at GENE_E
LOCUS = genome.Interval(chromosome="chr11", start=49_107_464, end=49_238_536)     # 128 kb
COMMON_WIN = genome.Interval(chromosome="chr11", start=48_950_000, end=49_474_288)  # 512 kb
ENH_S, ENH_E = 49_177_981, 49_178_322            # nominated gene-body element
PEAKS = {"promoter_peak": (49_207_936, 49_208_734),
         "genebody_peak1": (ENH_S, ENH_E),
         "upstream_20.7kb": (49_229_272, 49_229_430),
         "intron1_peak": (49_206_722, 49_206_871)}
SUMMITS = {"promoter_peak": 49_208_588, "genebody_peak1": 49_178_154,
           "upstream_20.7kb": 49_229_326, "intron1_peak": 49_206_797}
CPG = {"cg06980460": (49_208_996, "C", "T"),     # promoter, 354 bp upstream of the TSS
       "cg25346915": (49_208_484, "C", "T"),     # 0.2 kb from the TSS
       "cg10006235": (49_204_818, "C", "T")}     # 3.8 kb into intron 1


def score_variant(pos, ref, alt, win=COMMON_WIN, label=""):
    """The authors' scoring function: FOLH1 = minus-strand RNA-seq summed over the gene body."""
    var = genome.Variant(chromosome="chr11", position=int(pos),
                         reference_bases=ref, alternate_bases=alt)
    o = MODEL.predict_variant(
        interval=win, variant=var,
        requested_outputs=[dna_client.OutputType.RNA_SEQ, dna_client.OutputType.DNASE],
        ontology_terms=["UBERON:0002367"])
    i0, i1 = GENE_S - win.start, GENE_E - win.start
    ref_expr = float(o.reference.rna_seq.values[i0:i1, 1].sum())
    alt_expr = float(o.alternate.rna_seq.values[i0:i1, 1].sum())
    vi = int(pos) - win.start
    ref_dn = float(o.reference.dnase.values[vi - 250:vi + 250, 0].mean())
    alt_dn = float(o.alternate.dnase.values[vi - 250:vi + 250, 0].mean())
    return dict(label=label, pos=int(pos), ref=ref, alt=alt,
                window=f"{win.start}-{win.end}", ref_expr=ref_expr, alt_expr=alt_expr,
                expr_log2fc=float(np.log2((alt_expr + 1e-3) / (ref_expr + 1e-3))),
                dnase_delta=alt_dn - ref_dn)


def run_many(jobs, path, sleep=0.3):
    """jobs = list of dicts(pos, ref, alt, label[, win]); cached to `path`."""
    done = {}
    if os.path.exists(path):
        done = {r["label"]: r for r in pd.read_csv(path).to_dict("records")}
    rows = list(done.values())
    for j in jobs:
        if j["label"] in done:
            continue
        for attempt in range(4):
            try:
                r = score_variant(j["pos"], j["ref"], j["alt"],
                                  j.get("win", COMMON_WIN), j["label"])
                rows.append(r)
                print(f"  {j['label']:44s} log2FC={r['expr_log2fc']:+.4f}  "
                      f"dDNase={r['dnase_delta']:+.4f}", flush=True)
                break
            except Exception as e:
                print(f"  {j['label']}: attempt {attempt} {type(e).__name__} {str(e)[:90]}",
                      flush=True)
                time.sleep(4)
        pd.DataFrame(rows).to_csv(path, index=False)
        time.sleep(sleep)
    return pd.DataFrame(rows)


def ensembl_snvs(region, summit, n=3):
    s, e = region
    for attempt in range(6):
        try:
            r = requests.get(
                f"https://rest.ensembl.org/overlap/region/human/11:{s}-{e}?feature=variation",
                headers={"Content-Type": "application/json"}, timeout=60)
            if r.status_code == 200 and r.text.strip().startswith("["):
                vs = [v for v in r.json() if len("/".join(v.get("alleles", []))) == 3]
                vs = sorted(vs, key=lambda v: abs(v["start"] - summit))[:n]
                return [(v["id"], int(v["start"]), *"/".join(v["alleles"]).split("/")) for v in vs]
        except Exception:
            pass
        time.sleep(6)
    return []


# ======================================================================== A1
def stage_landscape():
    out = MODEL.predict_interval(
        interval=LOCUS,
        requested_outputs=[dna_client.OutputType.DNASE, dna_client.OutputType.RNA_SEQ,
                           dna_client.OutputType.CAGE],
        ontology_terms=["UBERON:0002367"])
    dnase = out.dnase.values[:, 0]
    pos = np.arange(LOCUS.start, LOCUS.end)
    sm = np.convolve(dnase, np.ones(129) / 129, mode="same")
    thr = np.percentile(sm, 95)
    segs, in_seg, s0 = [], False, 0
    for i, a in enumerate(sm > thr):
        if a and not in_seg:
            s0, in_seg = i, True
        elif not a and in_seg:
            segs.append((pos[s0], pos[i - 1], sm[s0:i].max(), pos[s0:i][np.argmax(sm[s0:i])]))
            in_seg = False
    segs = pd.DataFrame(segs, columns=["start", "end", "peak_height", "peak_pos"])
    segs["width"] = segs["end"] - segs["start"]

    def annot(r):
        if r["start"] <= GENE_E <= r["end"]:
            return "promoter/TSS"
        if r["start"] > GENE_E:
            return f"upstream(TSS)+{(r['start']-GENE_E)/1000:.1f}kb"
        if r["end"] < GENE_S:
            return f"downstream-{(GENE_S-r['end'])/1000:.1f}kb"
        return f"gene body ({(GENE_E-r['end'])/1000:.1f} kb from TSS)"

    segs["location"] = segs.apply(annot, axis=1)
    segs = segs.sort_values("peak_height", ascending=False)
    print(f"[A1] predicted prostate DNase peaks in the locus: {len(segs)}")
    print(segs.head(12).to_string(index=False))
    segs.to_csv(f"{OUT}/s09_dnase_peaks.csv", index=False)
    np.save(f"{CACHE}/dnase_smoothed.npy", sm)
    np.save(f"{CACHE}/rna_minus.npy", out.rna_seq.values[:, 1])
    print("\n[A1] RNA_SEQ tracks:", list(out.rna_seq.metadata["name"]))
    print("[A1] CAGE tracks:", list(out.cage.metadata["name"]))


# ======================================================================== A2
def stage_variants():
    jobs = []
    for name, region in PEAKS.items():
        for rsid, p, ref, alt in ensembl_snvs(region, SUMMITS[name]):
            jobs.append(dict(label=f"{name}:{rsid}", pos=p, ref=ref, alt=alt))
    for probe, (p, ref, alt) in CPG.items():
        jobs.append(dict(label=f"{probe}_CpG>{ref}>{alt}", pos=p, ref=ref, alt=alt))
    print(f"[A2] scoring {len(jobs)} variants in the common 512 kb window")
    df = run_many(jobs, f"{OUT}/s09_variant_effects.csv")
    print(df[["label", "pos", "ref", "alt", "expr_log2fc", "dnase_delta"]]
          .sort_values("expr_log2fc").round(4).to_string(index=False))


# ======================================================================== A3
def stage_ism(region=(ENH_S, ENH_E), step=4, tag="enhancer"):
    s, e = region
    seq = None
    for attempt in range(8):
        try:
            r = requests.get(f"https://rest.ensembl.org/sequence/region/human/11:{s}-{e-1}"
                             f"?coord_system_version=GRCh38",
                             headers={"Content-Type": "text/plain"}, timeout=60)
            if r.status_code == 200 and len(r.text.strip()) == e - s:
                seq = r.text.strip().upper()
                break
        except Exception:
            pass
        time.sleep(8)
    assert seq, "reference sequence fetch failed"
    print(f"[A3/{tag}] reference sequence {len(seq)} nt")
    comp = {"A": "C", "C": "G", "G": "T", "T": "A"}
    jobs = []
    for p in range(s, e, step):
        ref = seq[p - s]
        jobs.append(dict(label=f"ism_{tag}_{p}", pos=p, ref=ref, alt=comp[ref]))
    print(f"[A3/{tag}] {len(jobs)} positions")
    df = run_many(jobs, f"{OUT}/s09_ism_{tag}.csv")
    print(df["expr_log2fc"].describe().round(4).to_string())
    top = df.reindex(df["expr_log2fc"].abs().sort_values(ascending=False).index).head(10)
    print(top[["pos", "expr_log2fc", "dnase_delta"]].round(4).to_string(index=False))


# ======================================================================== B1
def stage_allele():
    """Is the cg25346915 effect CpG-disruption-specific, or a position effect?"""
    jobs = []
    for probe, (p, ref, _) in CPG.items():
        for alt in "ACGT":
            if alt == ref:
                continue
            jobs.append(dict(label=f"{probe}_{ref}>{alt}", pos=p, ref=ref, alt=alt))
    # also the two flanking positions of cg25346915, all alternate alleles
    p0 = CPG["cg25346915"][0]
    seq = requests.get(f"https://rest.ensembl.org/sequence/region/human/11:{p0-3}-{p0+3}"
                       f"?coord_system_version=GRCh38",
                       headers={"Content-Type": "text/plain"}, timeout=60).text.strip().upper()
    print(f"[B1] context around cg25346915 ({p0}): {seq}")
    for off in (-2, -1, 1, 2):
        ref = seq[3 + off]
        for alt in "ACGT":
            if alt == ref:
                continue
            jobs.append(dict(label=f"cg25346915_offset{off:+d}_{ref}>{alt}",
                             pos=p0 + off, ref=ref, alt=alt))
    print(f"[B1] scoring {len(jobs)} allele/position variants")
    df = run_many(jobs, f"{OUT}/s09_allele_specificity.csv")
    print(df[["label", "pos", "ref", "alt", "expr_log2fc"]].sort_values("pos")
          .round(4).to_string(index=False))


# ======================================================================== B3
def stage_window():
    """Does the headline effect depend on the prediction window?"""
    wins = {
        "128kb_locus": LOCUS,
        "512kb_common": COMMON_WIN,
        "512kb_shifted_+50kb": genome.Interval("chr11", 49_000_000, 49_524_288),
        "1Mb": genome.Interval("chr11", 48_700_000, 49_748_576),
    }
    jobs = []
    for wname, w in wins.items():
        for probe in ("cg25346915", "cg06980460"):
            p, ref, alt = CPG[probe]
            jobs.append(dict(label=f"{probe}@{wname}", pos=p, ref=ref, alt=alt, win=w))
        jobs.append(dict(label=f"genebody_summit@{wname}", pos=SUMMITS["genebody_peak1"],
                         ref=None, alt=None, win=w))
    jobs = [j for j in jobs if j["ref"] is not None]
    print(f"[B3] scoring {len(jobs)} window variants")
    df = run_many(jobs, f"{OUT}/s09_window_sensitivity.csv")
    print(df[["label", "ref_expr", "alt_expr", "expr_log2fc"]].round(4).to_string(index=False))


# ======================================================================== B4
def stage_determinism():
    p, ref, alt = CPG["cg25346915"]
    vals = []
    for i in range(3):
        r = score_variant(p, ref, alt, COMMON_WIN, f"repeat{i}")
        vals.append(r["expr_log2fc"])
        print(f"  repeat {i}: ref_expr={r['ref_expr']:.3f} alt_expr={r['alt_expr']:.3f} "
              f"log2FC={r['expr_log2fc']:+.6f}")
    print(f"[B4] spread across repeats = {max(vals)-min(vals):.2e}")
    json.dump({"repeats": vals}, open(f"{OUT}/s09_determinism.json", "w"), indent=1)


# ======================================================================== A4
def stage_contact():
    cm_win = genome.Interval("chr11", 48_877_856, 49_402_144)
    out = MODEL.predict_interval(interval=cm_win,
                                 requested_outputs=[dna_client.OutputType.CONTACT_MAPS],
                                 ontology_terms=["EFO:0003042"])
    cm_all = out.contact_maps.values
    nb = cm_all.shape[0]
    res = (cm_win.end - cm_win.start) / nb
    b = lambda p: int((p - cm_win.start) / res)
    enh_b, tss_b, g3_b = b(49_178_150), b(GENE_E), b(GENE_S)
    rows = []
    for tr in range(cm_all.shape[2]):
        cm = cm_all[:, :, tr]
        d = abs(tss_b - enh_b)
        obs = cm[enh_b, tss_b]
        bg = np.array([cm[i, i + d] for i in range(max(0, enh_b - 40), min(nb - d, enh_b + 40))])
        name = out.contact_maps.metadata.iloc[tr].get("Assay title", f"track{tr}")
        rows.append(dict(track=str(name), enh_promoter=float(obs), bg_mean=float(bg.mean()),
                         bg_sd=float(bg.std()), enrichment=float(obs / bg.mean()),
                         z=float((obs - bg.mean()) / (bg.std() + 1e-9))))
        print(f"  {name}: enh-promoter={obs:.3f} bg={bg.mean():.3f}+/-{bg.std():.3f} "
              f"enrichment={obs/bg.mean():.2f}x z={(obs-bg.mean())/(bg.std()+1e-9):.2f}")
    pd.DataFrame(rows).to_csv(f"{OUT}/s09_contact_maps.csv", index=False)


def stage_loci():
    rows = []
    for sym in ["TXNRD1", "SLC16A3", "SLC7A5", "PLAUR"]:
        g = None
        for attempt in range(5):
            try:
                r = requests.get(f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{sym}",
                                 headers={"Content-Type": "application/json"}, timeout=45)
                if r.status_code == 200:
                    g = r.json()
                    break
            except Exception:
                pass
            time.sleep(4)
        if not g:
            print(f"  {sym}: lookup failed")
            continue
        centre = (g["start"] + g["end"]) // 2
        win = genome.Interval(f"chr{g['seq_region_name']}", centre - 65536, centre + 65536)
        o = MODEL.predict_interval(
            interval=win,
            requested_outputs=[dna_client.OutputType.DNASE, dna_client.OutputType.RNA_SEQ],
            ontology_terms=["UBERON:0002367"])
        i0, i1 = g["start"] - win.start, g["end"] - win.start
        si = 0 if g["strand"] == 1 else 1
        expr = float(o.rna_seq.values[i0:i1, si].sum())
        tss = g["start"] if g["strand"] == 1 else g["end"]
        ti = tss - win.start
        prom = float(o.dnase.values[max(0, ti - 5000):ti + 5000, 0].max())
        rows.append(dict(gene=sym, chrom=g["seq_region_name"], start=g["start"], end=g["end"],
                         strand=g["strand"], expr_genebody=expr, promoter_dnase_max=prom))
        print(f"  {sym}: gene-body RNA={expr:.1f}, promoter DNase max={prom:.3f}")
    pd.DataFrame(rows).to_csv(f"{OUT}/s09_target_loci.csv", index=False)


STAGES = {"landscape": stage_landscape, "variants": stage_variants,
          "ism": lambda: stage_ism(),
          "ism_control": lambda: stage_ism((49_190_000, 49_190_341), 4, "control"),
          "allele": stage_allele, "window": stage_window,
          "determinism": stage_determinism, "contact": stage_contact, "loci": stage_loci}

if __name__ == "__main__":
    for name in sys.argv[1:]:
        print(f"\n{'='*22} {name} {'='*22}", flush=True)
        STAGES[name]()
