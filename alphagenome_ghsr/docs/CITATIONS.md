# Software, model and database citations

Execution date of the reported analysis: **2026-11-14** (UTC).
Genome reference assembly: **GRCh38/hg38** throughout.

---

## 1. Model

**AlphaGenome 2** (Google DeepMind), accessed through the official Python SDK
`alphagenome` version **2.0.0** (released 2026-11-11), module
`alphagenome.models.v2`, served by the public API endpoint
`gdmscience.googleapis.com:443`.

* Model version requested: default checkpoint (`ModelVersion.ALL_FOLDS`;
  the client sends no explicit `model_version`, so the server's current
  default model is used).
* Organism: `HOMO_SAPIENS`.
* SDK repository: <https://github.com/google-deepmind/alphagenome>
* SDK documentation: <https://www.alphagenomedocs.com/>
* SDK author / copyright: Google LLC. SDK licence: Apache-2.0.
* SDK summary string (from package metadata): "A Python SDK for interacting
  and visualizing genomic models."

> **Action required before submission.** The canonical literature citation for
> the AlphaGenome 2 model should be taken from the official documentation
> ("How to cite") at <https://www.alphagenomedocs.com/>. That host was **not
> reachable from the analysis environment**, so the exact title, author list,
> venue and DOI of the AlphaGenome 2 reference could not be verified here and
> **has deliberately not been reproduced from memory**. The primary methods
> reference for the AlphaGenome model family is the AlphaGenome preprint
> (Avsec et al., Google DeepMind, bioRxiv, 2025); its precise bibliographic
> details and DOI, and whether a separate AlphaGenome 2 reference now exists,
> must be confirmed against the official documentation. Do not cite this
> analysis as establishing either reference.

Because the API serves a hosted model, the exact served checkpoint cannot be
pinned from the client. The SDK version, the API host, the request parameters
and the execution date are recorded in
`data/raw/provenance_phase3_main.json` and are the reproducibility anchor for
these results.

### Scorers used

All 22 officially recommended human variant scorers were requested for every
variant/condition pair, via
`alphagenome.models.v2.variant_scorers.get_recommended_scorers(
Organism.HOMO_SAPIENS, cell_conditioning=False)`. The exact scorer list is
recorded verbatim in `data/raw/provenance_phase3_main.json`.

---

## 2. Databases

| Resource | Version / build | Endpoint | Accessed |
|---|---|---|---|
| NCBI dbSNP | build **157** (records last updated 2024-10-31) | `https://api.ncbi.nlm.nih.gov/variation/v0/refsnp` | 2026-11-14 |
| Ensembl REST API | GRCh38 | `https://rest.ensembl.org` | 2026-11-14 |
| Ensembl VEP (via REST) | GRCh38 | `https://rest.ensembl.org/vep/human/id` | 2026-11-14 |
| Ensembl Regulatory Build (via REST `overlap/region`) | GRCh38 | `https://rest.ensembl.org/overlap/region/human` | 2026-11-14 |
| gnomAD genomes allele frequencies (served through Ensembl `pops=1`) | gnomADg, as served by Ensembl on the access date | `https://rest.ensembl.org/variation/human` | 2026-11-14 |
| 1000 Genomes phase 3 allele frequencies (served through Ensembl) | phase 3 | `https://rest.ensembl.org/variation/human` | 2026-11-14 |
| GENCODE / Ensembl gene and transcript annotation (served through Ensembl) | GRCh38 | `https://rest.ensembl.org/lookup`, `overlap/region` | 2026-11-14 |

Every raw response is preserved verbatim under `data/raw/phase1/`, so the
annotation layer can be audited exactly as retrieved. The Ensembl release
number served on the access date is not returned by the endpoints used; the
preserved responses are the authoritative record.

Suggested citations for these resources (standard references, to be formatted
to the target journal's style):

* **dbSNP** - Sherry ST, Ward MH, Kholodov M, Baker J, Phan L, Smigielski EM,
  Sirotkin K. dbSNP: the NCBI database of genetic variation. *Nucleic Acids
  Research* 2001;29(1):308-311.
* **Ensembl** - the current Ensembl database paper for the release served on
  the access date (Ensembl publishes an annual database update in *Nucleic
  Acids Research*).
* **Ensembl VEP** - McLaren W, Gil L, Hunt SE, Riat HS, Ritchie GRS,
  Thormann A, Flicek P, Cunningham F. The Ensembl Variant Effect Predictor.
  *Genome Biology* 2016;17:122.
* **Ensembl Regulatory Build** - Zerbino DR, Wilder SP, Johnson N, Juettemann T,
  Flicek PR. The Ensembl Regulatory Build. *Genome Biology* 2015;16:56.
* **gnomAD** - the current gnomAD flagship paper for the version served
  (Chen S, Francioli LC, Goodrich JK, et al., *Nature* 2024;625:92-100 for
  gnomAD v4).
* **1000 Genomes Project** - The 1000 Genomes Project Consortium. A global
  reference for human genetic variation. *Nature* 2015;526:68-74.

---

## 3. Software

| Package | Version |
|---|---|
| Python | 3.13.18 |
| alphagenome | 2.0.0 |
| anndata | 0.13.4 |
| grpcio | 1.85.1 |
| matplotlib | 3.11.2 |
| numpy | 2.5.3 |
| pandas | 3.0.5 |
| pyarrow | 26.0.1 |
| requests | 2.34.2 |
| scipy | 1.18.1 |
| seaborn | 0.13.2 |
| nbformat | 5.11.1 |
| nbconvert | 7.17.2 |
| ipykernel | 7.4.0 |
| jupyter-core | 5.9.1 |
| python-docx | 1.2.0 |

Platform: Linux x86_64 (kernel 6.18.44). Full records, regenerated on every
run, are in `data/raw/provenance_*.json`.

Suggested citations:

* **NumPy** - Harris CR, Millman KJ, van der Walt SJ, et al. Array programming
  with NumPy. *Nature* 2020;585:357-362.
* **SciPy** - Virtanen P, Gommers R, Oliphant TE, et al. SciPy 1.0:
  fundamental algorithms for scientific computing in Python. *Nature Methods*
  2020;17:261-272.
* **pandas** - McKinney W. Data structures for statistical computing in
  Python. *Proceedings of the 9th Python in Science Conference* 2010;51-56.
* **Matplotlib** - Hunter JD. Matplotlib: a 2D graphics environment.
  *Computing in Science & Engineering* 2007;9(3):90-95.
* **AnnData** - Virshup I, Rybakov S, Theis FJ, Angerer P, Wolf FA. anndata:
  annotated data. *bioRxiv* 2021.12.16.473007.

---

## 4. Source study

Saeidi K, Khamisipour G, Motamed N, Farhadi A, Alizadeh M, Akbarzadeh S.
Ghrelin Receptor Gene Polymorphisms, Cognitive Screening Status, and
Osteoporosis-Related Parameters: An Exploratory Case-Control Study.
Manuscript, file `ARTICLE__revised-v2.docx`.

The cohort from which that study drew its participants:

* Shafiee G, Ostovar A, Heshmat R, et al. Bushehr Elderly Health (BEH)
  programme: study protocol and design of musculoskeletal system and cognitive
  function (stage II). *BMJ Open* 2017;7(8):e013606.

Level 2 and Level 3 evidence in `results/tables/table_3_integrated_evidence.csv`
is transcribed from that manuscript (not recomputed); the transcription,
with the manuscript table or section each value comes from, is in
`data/processed/manuscript_findings_transcription.json`.

---

## 5. Data and code availability statement (suggested wording)

> The computational analysis pipeline, all raw AlphaGenome 2 model outputs at
> per-track resolution, all derived tables and all figure source data are
> available at <repository URL>. Variant annotation used public records only
> (NCBI dbSNP build 157; Ensembl GRCh38). No participant-level genetic,
> clinical, biochemical or densitometric data were transmitted to any external
> service: the AlphaGenome API received only public GRCh38 variant coordinates
> and public ontology identifiers. Access to the AlphaGenome 2 API requires a
> key obtained from Google DeepMind; the key is read from an environment
> variable and is not contained in the repository.
