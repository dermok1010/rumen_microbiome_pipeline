"""The kernel prediction test (09) at ASV level, on the paper's shared 250 bp V4 core.

This is the comparison the 16S vs RE-RRS paper relies on across primer sets (loco_16s_v4core):
features are V4-core sequences, so a V4-trained model reads V3-V4 test samples sequence for
sequence, with no classifier in between. Same training samples, test cows, covariates and kernel
code as 09; only the features change (cores in >= 5% of training samples). Test sources are each
cow's autumn V4 sample (T2) and its V3-V4 library of the same rumen sample (the one with more reads).

Run from the repo root with ~/envs/micro_prediction/bin/python, after 10_asv_core_tables.py.
Outputs: outputs/region_comparison_v3v4_v4/prediction_kernel_asv/ (kernel_folds/, results/, summary.csv).
"""
import sys
from importlib import import_module
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
k09 = import_module("09_dairy_kernel_prediction")
pred07 = k09.pred07
CORE = pred07.OUT / "asv_core"


def read(name):
    return pd.read_csv(CORE / f"{name}.tsv", sep="\t", index_col=0, dtype={"sample": str})


def main():
    pheno = pred07.phenotypes()
    _, ptr = pred07.training_set(pheno)
    _, pte = pred07.test_sources(pheno)
    train = read("train_v4").loc[ptr["sample"]]
    sheet = pd.read_csv(pred07.OUT / "samples.csv", dtype={"cow": str, "sample": str})
    v4 = read("crt23_v4")
    t2 = {s: c for s, c in zip(sheet["sample"], sheet["cow"]) if s.startswith("CRT23__T2.")}
    v4_t2 = v4.loc[list(t2)].set_axis([t2[s] for s in t2])
    v34 = read("bsi_v3v4")
    main34 = sheet[sheet.region == "V3V4"].sort_values("reads", ascending=False).drop_duplicates("cow")
    v34 = v34.loc[main34["sample"]].set_axis(main34["cow"].tolist())
    src = {"V4_T2": v4_t2.loc[pte.cow], "V3V4": v34.loc[pte.cow]}
    k09.run_kernel(train, ptr, src, pte, pred07.OUT / "prediction_kernel_asv", "V4-core ASVs")


if __name__ == "__main__":
    main()
