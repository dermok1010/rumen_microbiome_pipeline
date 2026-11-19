"""Write the GG2 2024.09 backbone taxonomy truncated at genus (species rank dropped) as a plain TSV
for `qiime tools import`. Needed because training on all 28,054 species-level taxonomy strings was
killed for memory (exit 137) on the 16 GB VM; at genus level there are about 9,000 classes.
Only the genus call is used downstream, and classify-sklearn's confidence at genus is the summed
posterior of the species under that genus, so the genus calls should be close to a species-level
model's (priors per genus are the same; the split of k-mer profiles across species is lost).
Reads the .qza as a zip (data only). Run with python -I."""
import sys, zipfile
from pathlib import Path

ref, out = Path(sys.argv[1]), Path(sys.argv[2])
with zipfile.ZipFile(ref / "2024.09.backbone.tax.qza") as z:
    name = [n for n in z.namelist() if n.endswith("data/taxonomy.tsv")][0]
    lines = z.read(name).decode().splitlines()
n = 0
with open(out, "w") as f:
    f.write("Feature ID\tTaxon\n")
    for ln in lines[1:]:
        fid, tax = ln.split("\t")[:2]
        ranks = [r.strip() for r in tax.split(";")]
        f.write(fid + "\t" + "; ".join(r for r in ranks if not r.startswith("s__")) + "\n"); n += 1
print(n, "records ->", out)
