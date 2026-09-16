"""
Reconstructs tully_merge.R's original input (tully_data/genus_counts.csv,
which no longer exists on this VM) from the new combined genus table.

That old file was a Tully-only extract of the same genus_wide table
postprocessing/01_build_microbiome_matrices.ipynb produces, with ANI_ID/
SeqID columns split out of SampleID (e.g. "Tully__10732_S4" -> ANI_ID
"10732", SeqID "S4"). Run after 01_build_microbiome_matrices.ipynb;
writes output/tully_genus_counts.csv (gitignored).
"""

import pandas as pd

df = pd.read_csv("../output/genus_counts_no_controls_515F_806R.csv")
tully = df[df["SampleID"].str.startswith("Tully__")].copy()

extracted = tully["SampleID"].str.extract(r"^Tully__(?P<ANI_ID>[^_]+)_(?P<SeqID>S\d+)$")
tully["ANI_ID"] = extracted["ANI_ID"]
tully["SeqID"] = extracted["SeqID"]

assert tully["ANI_ID"].notna().all() and tully["SeqID"].notna().all(), \
    "some Tully SampleIDs didn't match the expected Tully__<ANI_ID>_<SeqID> pattern"

genus_cols = [c for c in tully.columns if c not in ("SampleID", "ANI_ID", "SeqID")]
tully = tully[["SampleID", "ANI_ID", "SeqID"] + genus_cols]

tully.to_csv("output/tully_genus_counts.csv", index=False)
print(tully.shape, "written to output/tully_genus_counts.csv")
