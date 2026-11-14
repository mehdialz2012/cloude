"""Phase III-b - REF vs ALT base-resolution track predictions.

`score_variant` returns aggregated scores; this module calls `predict_variant`
to retrieve the underlying base-resolution REF and ALT prediction tracks around
the GHSR locus. These are needed for:

  * Figure 2 (REF vs ALT profiles at the variant and at the GHSR promoter);
  * an independent check that the aggregated scores in Phase III are consistent
    with the raw tracks they are computed from (Phase IV).

A 100 kb context is used so that the GHSR gene body, its TSS and the three
variants all fall inside the predicted window at a resolution that can be
plotted, while keeping the stored arrays small.

Run:  python -m ghsr_ag2.phase3b_tracks
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
from alphagenome.data import genome
from alphagenome.models.v2 import dna_client, dna_output

from . import config, provenance

RAW = config.DATA_RAW / "phase3b"
RAW.mkdir(parents=True, exist_ok=True)

# Output types retrieved as base-resolution tracks. Restricted to the
# mechanistically plausible primary modalities to keep array sizes manageable.
TRACK_OUTPUTS = (
    dna_output.OutputType.RNA_SEQ,
    dna_output.OutputType.CAGE,
    dna_output.OutputType.DNASE,
    dna_output.OutputType.ATAC,
    dna_output.OutputType.CHIP_TF,
    dna_output.OutputType.CHIP_HISTONE,
)

SEQUENCE_LENGTH_KEY = "SEQUENCE_LENGTH_100KB"

# Conditions plotted in Figure 2: the best-supported representative of each
# biological panel, chosen before inspecting the effect sizes.
FIGURE_CONDITIONS = (
    "UBERON:0000007",             # pituitary gland - highest predicted GHSR
    "UBERON:0001898",             # hypothalamus
    "UBERON:0000955",             # brain (fully observed)
    "UBERON:0002481_CL:0000062",  # osteoblast
    "CL:0000137",                 # osteocyte
    "UBERON:0000945",             # stomach (fully observed reference)
)


def extract_profiles(
    client: dna_client.DnaClient,
    spec: config.VariantSpec,
    condition_curie: str,
) -> dict[str, object]:
    variant = genome.Variant(
        chromosome=config.CHROM,
        position=spec.position,
        reference_bases=spec.ref,
        alternate_bases=spec.alt,
        name=f"{spec.rsid}_{spec.ref}>{spec.alt}",
    )
    interval = variant.reference_interval.resize(
        dna_client.SUPPORTED_SEQUENCE_LENGTHS[SEQUENCE_LENGTH_KEY]
    )
    out = client.predict_variant(
        interval=interval,
        variant=variant,
        requested_outputs=TRACK_OUTPUTS,
        condition=condition_curie,
    )
    record: dict[str, object] = {
        "variant_key": spec.key,
        "rsid": spec.rsid,
        "condition_curie": condition_curie,
        "condition_label": config.CONDITION_BY_CURIE[condition_curie].label,
        "interval": str(interval),
        "interval_start": interval.start,
        "interval_end": interval.end,
        "sequence_length_key": SEQUENCE_LENGTH_KEY,
        "outputs": {},
    }
    for ot in TRACK_OUTPUTS:
        ref_td = out.reference.get(ot)
        alt_td = out.alternate.get(ot)
        if ref_td is None or alt_td is None:
            record["outputs"][ot.name] = {"available": False}
            continue
        ref = np.asarray(ref_td.values, dtype=np.float32)
        alt = np.asarray(alt_td.values, dtype=np.float32)
        meta = ref_td.metadata.reset_index(drop=True)
        np.savez_compressed(
            RAW / f"tracks_{spec.key}_{condition_curie.replace(':', '-')}_{ot.name}.npz",
            ref=ref,
            alt=alt,
            interval_start=np.int64(interval.start),
            interval_end=np.int64(interval.end),
            resolution=np.int64(ref_td.resolution),
        )
        meta.to_csv(
            RAW
            / f"trackmeta_{spec.key}_{condition_curie.replace(':', '-')}_{ot.name}.csv",
            index=False,
        )
        record["outputs"][ot.name] = {
            "available": True,
            "n_positions": int(ref.shape[0]),
            "n_tracks": int(ref.shape[1]),
            "resolution_bp": int(ref_td.resolution),
            "max_abs_alt_minus_ref": float(np.nanmax(np.abs(alt - ref))),
            "n_observed_tracks": int(
                (meta["imputation_mode"] == "observed").sum()
            )
            if "imputation_mode" in meta.columns
            else None,
        }
    return record


def main() -> None:
    print("Phase III-b: REF vs ALT base-resolution track predictions")
    client = dna_client.create(config.api_key(), timeout=120)
    records = []
    for spec in config.VARIANTS:
        for curie in FIGURE_CONDITIONS:
            rec = extract_profiles(client, spec, curie)
            records.append(rec)
            avail = {
                k: v.get("max_abs_alt_minus_ref")
                for k, v in rec["outputs"].items()
                if v.get("available")
            }
            print(
                f"[phase3b] {spec.key:<18} {rec['condition_label']:<24} "
                + " ".join(
                    f"{k}:{(v if v is not None else float('nan')):.4g}"
                    for k, v in avail.items()
                )
            )
    (RAW / "phase3b_index.json").write_text(json.dumps(records, indent=2))
    pd.json_normalize(records).to_csv(
        config.DATA_PROCESSED / "phase3b_track_index.csv", index=False
    )
    provenance.write(
        "phase3b",
        {
            "sequence_length_key": SEQUENCE_LENGTH_KEY,
            "requested_outputs": [o.name for o in TRACK_OUTPUTS],
            "figure_conditions": list(FIGURE_CONDITIONS),
            "n_records": len(records),
            "raw_output_dir": str(RAW),
        },
    )
    print(f"[phase3b] wrote {len(records)} track records -> {RAW}")


if __name__ == "__main__":
    main()
