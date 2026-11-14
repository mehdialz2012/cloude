"""Phase I - variant identification, coordinate verification and annotation.

Independently verifies, from primary public sources and without relying on the
source manuscript:

  1. GRCh38 coordinates of the three dbSNP variants (NCBI dbSNP + Ensembl).
  2. Reference and alternate alleles, with explicit strand orientation.
  3. That each stated REF base matches the GRCh38 reference sequence.
  4. The genomic relationship of each variant to GHSR (TSS distance,
     transcripts, VEP consequence).
  5. Overlap with Ensembl Regulatory Build features.
  6. Protein-coding and non-coding genes in the 1 Mb scoring window.
  7. Population allele frequencies (gnomAD v4 / 1000 Genomes via Ensembl).
  8. The genomic location of the manuscript's PCR amplicon, by mapping the
     published primer pair onto GRCh38.

All raw API responses are preserved under data/raw/phase1/.

Run:  python -m ghsr_ag2.phase1_annotation
"""

from __future__ import annotations

import json
import pathlib
import time
from typing import Any

import pandas as pd
import requests

from . import config, provenance

ENSEMBL = "https://rest.ensembl.org"
DBSNP = "https://api.ncbi.nlm.nih.gov/variation/v0/refsnp"
RAW = config.DATA_RAW / "phase1"
RAW.mkdir(parents=True, exist_ok=True)

# A generic client name; no user or project identifying information is sent.
HEADERS = {
    "Content-Type": "application/json",
    "Accept": "application/json",
    "User-Agent": "alphagenome-ghsr-analysis/1.0",
}

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def _get(url: str, cache_name: str, params: dict | None = None) -> Any:
    """GET with on-disk caching of the raw response."""
    path = RAW / f"{cache_name}.json"
    if path.exists():
        return json.loads(path.read_text())
    for attempt in range(4):
        try:
            r = requests.get(url, headers=HEADERS, params=params, timeout=120)
            r.raise_for_status()
            payload = r.json()
            path.write_text(json.dumps(payload, indent=1))
            return payload
        except Exception as exc:  # noqa: BLE001 - retry transient failures
            if attempt == 3:
                raise RuntimeError(f"failed to GET {url}: {exc}") from exc
            time.sleep(2 * (attempt + 1))
    raise AssertionError("unreachable")


# --------------------------------------------------------------------------
# Individual retrieval steps
# --------------------------------------------------------------------------
def fetch_dbsnp(rsid: str) -> dict[str, Any]:
    """Return the primary top-level placement (PTLP) alleles from dbSNP."""
    payload = _get(f"{DBSNP}/{rsid.removeprefix('rs')}", f"dbsnp_{rsid}")
    snap = payload.get("primary_snapshot_data", {})
    out: dict[str, Any] = {
        "rsid": f"rs{payload.get('refsnp_id')}",
        "dbsnp_last_update": payload.get("last_update_date"),
        "dbsnp_build": payload.get("last_update_build_id"),
        "variant_type": snap.get("variant_type"),
        "placements": [],
    }
    for pl in snap.get("placements_with_allele", []):
        if not pl.get("is_ptlp"):
            continue
        assemblies = [
            a.get("assembly_name")
            for a in pl.get("placement_annot", {}).get(
                "seq_id_traits_by_assembly", []
            )
        ]
        for al in pl.get("alleles", []):
            spdi = al.get("allele", {}).get("spdi", {})
            deleted = spdi.get("deleted_sequence")
            inserted = spdi.get("inserted_sequence")
            if deleted == inserted:
                continue  # the reference "allele" record
            out["placements"].append(
                {
                    "seq_id": spdi.get("seq_id"),
                    "assemblies": assemblies,
                    # SPDI positions are 0-based; convert to 1-based.
                    "position_1based": int(spdi.get("position")) + 1,
                    "ref": deleted,
                    "alt": inserted,
                    "hgvs": al.get("hgvs"),
                }
            )
    return out


def fetch_ensembl_variation(rsid: str) -> dict[str, Any]:
    payload = _get(
        f"{ENSEMBL}/variation/human/{rsid}",
        f"ensembl_variation_{rsid}",
        {"content-type": "application/json", "pops": "1"},
    )
    mappings = [
        m
        for m in payload.get("mappings", [])
        if m.get("assembly_name") == config.ENSEMBL_ASSEMBLY
    ]
    return {
        "rsid": payload.get("name"),
        "synonyms": payload.get("synonyms", []),
        "minor_allele": payload.get("minor_allele"),
        "MAF": payload.get("MAF"),
        "most_severe_consequence_variation_api": payload.get(
            "most_severe_consequence"
        ),
        "mappings": [
            {
                "chrom": m.get("seq_region_name"),
                "start": m.get("start"),
                "end": m.get("end"),
                "allele_string": m.get("allele_string"),
                "strand": m.get("strand"),
            }
            for m in mappings
        ],
        "populations": payload.get("populations", []),
    }


def fetch_vep(rsid: str) -> dict[str, Any]:
    payload = _get(f"{ENSEMBL}/vep/human/id/{rsid}", f"vep_{rsid}")
    rec = payload[0]
    tx = []
    seen = set()
    for t in rec.get("transcript_consequences", []):
        k = (t.get("transcript_id"), tuple(t.get("consequence_terms", [])))
        if k in seen:
            continue
        seen.add(k)
        tx.append(
            {
                "gene_symbol": t.get("gene_symbol"),
                "gene_id": t.get("gene_id"),
                "transcript_id": t.get("transcript_id"),
                "biotype": t.get("biotype"),
                "consequence_terms": t.get("consequence_terms"),
                "distance": t.get("distance"),
                "strand": t.get("strand"),
            }
        )
    return {
        "rsid": rsid,
        "allele_string": rec.get("allele_string"),
        "most_severe_consequence": rec.get("most_severe_consequence"),
        "transcript_consequences": tx,
        "regulatory_feature_consequences": rec.get(
            "regulatory_feature_consequences", []
        ),
        "motif_feature_consequences": rec.get("motif_feature_consequences", []),
    }


def fetch_ghsr_gene() -> dict[str, Any]:
    payload = _get(
        f"{ENSEMBL}/lookup/symbol/homo_sapiens/{config.GHSR_SYMBOL}",
        "ensembl_ghsr_gene",
        {"expand": "1", "content-type": "application/json"},
    )
    transcripts = []
    for t in payload.get("Transcript", []):
        tss = t["start"] if t["strand"] == 1 else t["end"]
        transcripts.append(
            {
                "transcript_id": t["id"],
                "name": t.get("display_name"),
                "biotype": t["biotype"],
                "start": t["start"],
                "end": t["end"],
                "strand": t["strand"],
                "tss": tss,
                "is_canonical": bool(t.get("is_canonical")),
            }
        )
    return {
        "gene_id": payload["id"],
        "symbol": payload.get("display_name"),
        "description": payload.get("description"),
        "chrom": payload["seq_region_name"],
        "start": payload["start"],
        "end": payload["end"],
        "strand": payload["strand"],
        "biotype": payload["biotype"],
        "assembly": payload.get("assembly_name"),
        "transcripts": transcripts,
    }


def fetch_sequence(start: int, end: int, tag: str) -> str:
    payload = _get(
        f"{ENSEMBL}/sequence/region/human/"
        f"{config.ENSEMBL_CHROM}:{start}..{end}:1",
        f"ensembl_sequence_{tag}",
        {
            "content-type": "application/json",
            "coord_system_version": config.ENSEMBL_ASSEMBLY,
        },
    )
    return payload["seq"].upper()


def fetch_regulatory_features(start: int, end: int) -> list[dict[str, Any]]:
    payload = _get(
        f"{ENSEMBL}/overlap/region/human/"
        f"{config.ENSEMBL_CHROM}:{start}-{end}",
        f"ensembl_regulatory_{start}_{end}",
        {"feature": "regulatory", "content-type": "application/json"},
    )
    return sorted(payload, key=lambda f: f["start"])


def fetch_genes_in_window(start: int, end: int) -> list[dict[str, Any]]:
    payload = _get(
        f"{ENSEMBL}/overlap/region/human/"
        f"{config.ENSEMBL_CHROM}:{start}-{end}",
        f"ensembl_genes_{start}_{end}",
        {"feature": "gene", "content-type": "application/json"},
    )
    return payload


# --------------------------------------------------------------------------
# Derived tables
# --------------------------------------------------------------------------
def _pop_freq(populations: list[dict], allele: str, pop: str) -> float | None:
    for p in populations:
        if p.get("population") == pop and p.get("allele") == allele:
            return float(p["frequency"])
    return None


def build_variant_annotation() -> tuple[pd.DataFrame, dict[str, Any]]:
    gene = fetch_ghsr_gene()
    canonical = next(
        t for t in gene["transcripts"] if t["transcript_id"] == config.GHSR_CANONICAL_TX
    )
    ghsr_tss = canonical["tss"]

    # Reference sequence spanning the amplicon plus flanks, for REF checks and
    # primer mapping.
    seq_start, seq_end = 172_449_100, 172_449_700
    seq = fetch_sequence(seq_start, seq_end, f"{seq_start}_{seq_end}")

    rows = []
    per_rsid: dict[str, Any] = {}
    for spec in config.VARIANTS:
        if spec.rsid not in per_rsid:
            per_rsid[spec.rsid] = {
                "dbsnp": fetch_dbsnp(spec.rsid),
                "ensembl": fetch_ensembl_variation(spec.rsid),
                "vep": fetch_vep(spec.rsid),
            }
        dbsnp = per_rsid[spec.rsid]["dbsnp"]
        ens = per_rsid[spec.rsid]["ensembl"]
        vep = per_rsid[spec.rsid]["vep"]

        # --- coordinate agreement between the two primary sources ----------
        dbsnp_positions = {p["position_1based"] for p in dbsnp["placements"]}
        ens_positions = {m["start"] for m in ens["mappings"]}
        coord_agree = (
            dbsnp_positions == ens_positions == {spec.position}
        )

        # --- REF base check against the GRCh38 reference sequence ----------
        observed_ref = seq[spec.position - seq_start]
        ref_matches = observed_ref == spec.ref

        # --- allele presence in dbSNP --------------------------------------
        dbsnp_alts = {
            p["alt"] for p in dbsnp["placements"] if p["ref"] == spec.ref
        }
        alt_in_dbsnp = spec.alt in dbsnp_alts

        # --- Ensembl allele string / strand --------------------------------
        mapping = ens["mappings"][0] if ens["mappings"] else {}
        ens_alleles = mapping.get("allele_string", "")
        ens_strand = mapping.get("strand")

        # --- strand reconciliation with the manuscript ---------------------
        m_ref, m_alt = (
            spec.manuscript_notation.split(">")
            if ">" in spec.manuscript_notation
            else ("", "")
        )
        plus_matches_manuscript = (m_ref, m_alt) == (spec.ref, spec.alt)
        minus_matches_manuscript = (m_ref, m_alt) == (
            revcomp(spec.ref),
            revcomp(spec.alt),
        )
        if plus_matches_manuscript:
            orientation = "plus strand (reference-genome orientation)"
        elif minus_matches_manuscript:
            orientation = "minus strand (GHSR gene orientation)"
        else:
            orientation = "neither strand as written (REF/ALT inverted)"

        # --- distance to GHSR TSS (gene is on the minus strand) ------------
        # Positive distance = upstream of the TSS in the direction of
        # transcription.
        tss_distance_upstream = spec.position - ghsr_tss

        # --- VEP consequence for GHSR --------------------------------------
        ghsr_tc = [
            t
            for t in vep["transcript_consequences"]
            if t.get("gene_symbol") == config.GHSR_SYMBOL
            and t.get("transcript_id") == config.GHSR_CANONICAL_TX
        ]
        consequence = (
            ",".join(sorted(set(ghsr_tc[0]["consequence_terms"])))
            if ghsr_tc
            else "none reported"
        )
        vep_distance = ghsr_tc[0]["distance"] if ghsr_tc else None

        # --- regulatory feature overlap ------------------------------------
        regs = fetch_regulatory_features(172_439_000, 172_459_000)
        overlapping = [
            f"{f['id']}({f.get('description')})"
            for f in regs
            if f["start"] <= spec.position <= f["end"]
        ]
        nearest = min(
            regs,
            key=lambda f: 0
            if f["start"] <= spec.position <= f["end"]
            else min(
                abs(f["start"] - spec.position), abs(f["end"] - spec.position)
            ),
        )
        nearest_dist = (
            0
            if nearest["start"] <= spec.position <= nearest["end"]
            else min(
                abs(nearest["start"] - spec.position),
                abs(nearest["end"] - spec.position),
            )
        )

        pops = ens["populations"]
        rows.append(
            {
                "variant_key": spec.key,
                "rsid": spec.rsid,
                "manuscript_label": spec.manuscript_label,
                "chrom": config.CHROM,
                "position_grch38": spec.position,
                "ref_plus_strand": spec.ref,
                "alt_plus_strand": spec.alt,
                "is_primary_configuration": spec.is_primary,
                "hgvs_grch38": f"NC_000003.12:g.{spec.position}{spec.ref}>{spec.alt}",
                "dbsnp_build": dbsnp["dbsnp_build"],
                "dbsnp_last_update": dbsnp["dbsnp_last_update"],
                "dbsnp_variant_type": dbsnp["variant_type"],
                "dbsnp_all_alts": ",".join(sorted(dbsnp_alts)),
                "alt_present_in_dbsnp": alt_in_dbsnp,
                "ensembl_allele_string": ens_alleles,
                "ensembl_strand": ens_strand,
                "coordinates_agree_dbsnp_ensembl": coord_agree,
                "grch38_ref_base_observed": observed_ref,
                "ref_matches_reference_genome": ref_matches,
                "manuscript_notation": spec.manuscript_notation,
                "manuscript_orientation_resolved": orientation,
                "ensembl_minor_allele": ens["minor_allele"],
                "ensembl_MAF": ens["MAF"],
                "af_gnomadg_all": _pop_freq(pops, spec.alt, "gnomADg:ALL"),
                "af_gnomadg_mid": _pop_freq(pops, spec.alt, "gnomADg:mid"),
                "af_gnomadg_nfe": _pop_freq(pops, spec.alt, "gnomADg:nfe"),
                "af_gnomadg_sas": _pop_freq(pops, spec.alt, "gnomADg:sas"),
                "af_gnomadg_eas": _pop_freq(pops, spec.alt, "gnomADg:eas"),
                "af_gnomadg_afr": _pop_freq(pops, spec.alt, "gnomADg:afr"),
                "af_1000g_all": _pop_freq(
                    pops, spec.alt, "1000GENOMES:phase_3:ALL"
                ),
                "vep_most_severe_consequence": vep["most_severe_consequence"],
                "vep_consequence_ghsr_canonical": consequence,
                "vep_distance_to_ghsr_bp": vep_distance,
                "ghsr_gene_id": gene["gene_id"],
                "ghsr_strand": gene["strand"],
                "ghsr_canonical_tss": ghsr_tss,
                "bp_upstream_of_ghsr_tss": tss_distance_upstream,
                "inside_ghsr_gene_body": gene["start"]
                <= spec.position
                <= gene["end"],
                "overlapping_regulatory_features": ";".join(overlapping)
                or "none",
                "nearest_regulatory_feature": f"{nearest['id']}"
                f"({nearest.get('description')})",
                "nearest_regulatory_feature_distance_bp": nearest_dist,
                "inside_manuscript_amplicon": config.AMPLICON_START
                <= spec.position
                <= config.AMPLICON_END,
                "provenance_note": spec.note,
            }
        )

    extra = {
        "ghsr_gene": gene,
        "sequence_window": {"start": seq_start, "end": seq_end},
        "regulatory_features": fetch_regulatory_features(
            172_439_000, 172_459_000
        ),
    }
    return pd.DataFrame(rows), extra


def build_amplicon_check(seq_start: int = 172_449_100) -> dict[str, Any]:
    """Map the published primer pair onto GRCh38 and verify product size."""
    seq = fetch_sequence(seq_start, 172_449_700, f"{seq_start}_172449700")
    fwd, rev = config.PRIMER_FORWARD, config.PRIMER_REVERSE

    def locate(primer: str) -> dict[str, Any]:
        i = seq.find(primer)
        if i >= 0:
            return {
                "match": "plus_strand",
                "start": seq_start + i,
                "end": seq_start + i + len(primer) - 1,
            }
        j = seq.find(revcomp(primer))
        if j >= 0:
            return {
                "match": "minus_strand",
                "start": seq_start + j,
                "end": seq_start + j + len(primer) - 1,
            }
        return {"match": "not_found", "start": None, "end": None}

    f_loc, r_loc = locate(fwd), locate(rev)
    coords = [
        c
        for loc in (f_loc, r_loc)
        for c in (loc["start"], loc["end"])
        if c is not None
    ]
    amp_start, amp_end = (min(coords), max(coords)) if coords else (None, None)
    size = (amp_end - amp_start + 1) if coords else None
    return {
        "forward_primer": fwd,
        "forward_primer_location": f_loc,
        "reverse_primer": rev,
        "reverse_primer_location": r_loc,
        "mapped_amplicon_start": amp_start,
        "mapped_amplicon_end": amp_end,
        "mapped_amplicon_size_bp": size,
        "reported_amplicon_size_bp": config.REPORTED_AMPLICON_BP,
        "size_matches_manuscript": size == config.REPORTED_AMPLICON_BP,
        "variants_inside_amplicon": [
            v.rsid
            for v in config.PRIMARY_VARIANTS
            if amp_start is not None and amp_start <= v.position <= amp_end
        ],
    }


def build_gene_window_table(
    center: int, half_width: int = 524_288
) -> pd.DataFrame:
    genes = fetch_genes_in_window(center - half_width, center + half_width)
    rows = []
    for g in genes:
        tss = g["start"] if g["strand"] == 1 else g["end"]
        inside = g["start"] <= center <= g["end"]
        rows.append(
            {
                "gene_id": g["id"],
                "gene_symbol": g.get("external_name") or g["id"],
                "biotype": g.get("biotype"),
                "start": g["start"],
                "end": g["end"],
                "strand": g["strand"],
                "tss": tss,
                "signed_distance_variant_to_tss_bp": tss - center,
                "abs_distance_to_tss_bp": abs(tss - center),
                "variant_inside_gene_body": inside,
            }
        )
    df = pd.DataFrame(rows).sort_values("abs_distance_to_tss_bp")
    return df.reset_index(drop=True)


# --------------------------------------------------------------------------
def main() -> None:
    print("Phase I: variant annotation and coordinate verification")
    annot, extra = build_variant_annotation()
    amplicon = build_amplicon_check()
    genes = build_gene_window_table(config.VARIANTS[0].position)

    annot_path = config.TABLES / "table_S1_variant_annotation.csv"
    annot.to_csv(annot_path, index=False)
    genes_path = config.TABLES / "table_S2_genes_in_1mb_window.csv"
    genes.to_csv(genes_path, index=False)

    (config.DATA_PROCESSED / "phase1_amplicon_mapping.json").write_text(
        json.dumps(amplicon, indent=2)
    )
    (config.DATA_PROCESSED / "phase1_ghsr_gene.json").write_text(
        json.dumps(extra["ghsr_gene"], indent=2)
    )
    (config.DATA_PROCESSED / "phase1_regulatory_features.json").write_text(
        json.dumps(extra["regulatory_features"], indent=2)
    )

    # ---- verification summary --------------------------------------------
    checks = {
        "all_ref_alleles_match_grch38": bool(
            annot["ref_matches_reference_genome"].all()
        ),
        "all_coordinates_agree_dbsnp_vs_ensembl": bool(
            annot["coordinates_agree_dbsnp_ensembl"].all()
        ),
        "all_alts_present_in_dbsnp": bool(annot["alt_present_in_dbsnp"].all()),
        "all_variants_inside_manuscript_amplicon": bool(
            annot["inside_manuscript_amplicon"].all()
        ),
        "amplicon_size_matches_manuscript": amplicon["size_matches_manuscript"],
        "no_variant_inside_ghsr_gene_body": bool(
            (~annot["inside_ghsr_gene_body"]).all()
        ),
        "no_variant_overlaps_ensembl_regulatory_feature": bool(
            (annot["overlapping_regulatory_features"] == "none").all()
        ),
    }
    (config.DATA_PROCESSED / "phase1_checks.json").write_text(
        json.dumps(checks, indent=2)
    )

    print(annot[
        [
            "variant_key",
            "position_grch38",
            "ref_plus_strand",
            "alt_plus_strand",
            "ref_matches_reference_genome",
            "bp_upstream_of_ghsr_tss",
            "vep_consequence_ghsr_canonical",
            "manuscript_orientation_resolved",
        ]
    ].to_string(index=False))
    print("\nAmplicon mapping:")
    print(json.dumps(amplicon, indent=2))
    print("\nVerification checks:")
    print(json.dumps(checks, indent=2))
    print(f"\nWrote {annot_path}")
    print(f"Wrote {genes_path}")

    provenance.write(
        "phase1",
        {
            "databases": {
                "NCBI dbSNP": {
                    "endpoint": DBSNP,
                    "build": str(annot["dbsnp_build"].iloc[0]),
                },
                "Ensembl REST": {
                    "endpoint": ENSEMBL,
                    "assembly": config.ENSEMBL_ASSEMBLY,
                },
            },
            "checks": checks,
            "raw_response_dir": str(RAW),
        },
    )


if __name__ == "__main__":
    main()
