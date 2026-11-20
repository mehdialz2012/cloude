"""Extract the FOLH1-locus HM450 probes (hg19 chr11:48.9-49.25 Mb) from the TCGA-PRAD
Xena methylation matrix, mirroring the authors' probe window."""
import gzip
import numpy as np
import pandas as pd
from common import DATA

pm = pd.read_csv(f"{DATA}/illumina450_probemap.tsv", sep="\t")
pm.columns = [c.strip("#") for c in pm.columns]
print("probemap cols:", pm.columns.tolist())
idc, chrc, startc = pm.columns[0], "chrom", "chromStart"
sel = pm[(pm[chrc].astype(str).isin(["chr11", "11"])) &
         (pm[startc] >= 48_900_000) & (pm[startc] <= 49_250_000)]
print(f"probes in chr11:48.90-49.25Mb: {len(sel)}")
want = set(sel[idc])

rows, header = [], None
with gzip.open(f"{DATA}/tcga_prad_meth450.gz", "rt") as fh:
    header = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        pid = line.split("\t", 1)[0]
        if pid in want:
            rows.append(line.rstrip("\n").split("\t"))
print(f"rows pulled: {len(rows)}")
m = pd.DataFrame(rows).set_index(0)
m.columns = header[1:]
m = m.apply(pd.to_numeric, errors="coerce")
m.to_csv(f"{DATA}/folh1_locus_meth.csv")
ann = sel.set_index(idc)[[chrc, startc, "gene"]]
ann.to_csv(f"{DATA}/folh1_locus_meth_annot.csv")
print(m.shape, "saved")
print(m.iloc[:3, :4])
