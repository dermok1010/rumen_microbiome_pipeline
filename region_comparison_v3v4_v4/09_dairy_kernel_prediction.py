"""The 07 prediction test refitted with the 16S vs RE-RRS paper's own kernel model.

07 used ridge and RF; the paper's main microbiome model is the REML kernel model ("mrm"), so this
reruns the same question with that model's code, unchanged, from micro_prediction_NZ:
  ~/microbiome_prediction/03_export_kernels.py  kernel_blocks(), full_rank_design()
  ~/microbiome_prediction/common.py             ASVTransform, CovariateTransform
  ~/microbiome_prediction/04_mixed_model.R      the REML fit and scoring
Training set, test cows, test profiles and covariates are 07's (its loaders are reused). Features
follow the paper's genus configs instead of 07's rule: CLR (pseudocount 0.5) of genera present in
>= 5% of training samples. No per-dataset centring: the kernel is centred on the training mean, and
shifting a test source by a constant only shifts its predictions, so it changes the bias, not r.

Each test source is written as one "split" (V4_T2, V4_T1, V3V4, V3V4_core) with the same training
set, then 04_mixed_model.R fits microbiome, covariates and combined for each.

Run from the repo root with ~/envs/micro_prediction/bin/python, after 07's inputs exist.
Outputs: outputs/region_comparison_v3v4_v4/prediction_kernel/ (kernel_folds/, results/, summary.csv).
"""
import subprocess
import sys
from importlib import import_module
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
MP = Path.home() / "microbiome_prediction"
sys.path[:0] = [str(HERE), str(MP)]
pred07 = import_module("07_dairy_prediction")
kern = import_module("03_export_kernels")
from common import ASVTransform, CovariateTransform  # noqa: E402

OUT = pred07.OUT / "prediction_kernel"
PREVALENCE, PSEUDO = 0.05, 0.5  # the paper's genus configs (re_rrs_comparison/configs/*_16s_genus.json)
CFG = {"numeric_covariates": pred07.COVS, "categorical_covariates": []}


def main():
    pheno = pred07.phenotypes()
    gtr, ptr = pred07.training_set(pheno)
    src, pte = pred07.test_sources(pheno)
    run_kernel(gtr.drop(columns="Unassigned_genus", errors="ignore"), ptr, src, pte, OUT, "genera")


def run_kernel(train_counts, ptr, src, pte, out_dir, unit):
    """Fit the paper's kernel model once per test source.

    train_counts: training samples x features, rows in ptr order; src: test source -> cows x features
    counts, rows in pte order; features missing from a test source count as zero."""
    if out_dir.exists():
        raise FileExistsError(f"{out_dir} already exists; move it aside to rerun")
    folds, results = out_dir / "kernel_folds", out_dir / "results"
    folds.mkdir(parents=True)
    genera = train_counts.columns
    rtr = train_counts.to_numpy(float)
    tx = ASVTransform(PREVALENCE, PSEUDO).fit(rtr)
    cx = CovariateTransform(CFG).fit(ptr)
    print(f"train: {len(ptr)} samples, {int(tx.keep.sum())} {unit} at >= {PREVALENCE:.0%} prevalence; "
          f"test: {len(pte)} CRT23 cows")
    kept = genera[tx.keep]
    for name, counts in src.items():
        c = counts.reindex(columns=genera, fill_value=0)
        print(f"  {name}: {100 * c[kept].to_numpy().sum() / counts.to_numpy().sum():.1f}% of test reads "
              f"in the retained {unit}")
    train = pd.DataFrame({"SampleID": ptr["sample"], "EID": "train_" + ptr["sample"].astype(str),
                          "study": ptr["dataset"], "y": ptr["y"]})
    for name, counts in src.items():
        out = folds / name
        out.mkdir()
        rte = counts.reindex(columns=genera, fill_value=0).to_numpy(float)
        ktt, kvt = kern.kernel_blocks(tx.transform(rtr), tx.transform(rte))
        a = np.column_stack([np.ones(len(ptr)), cx.transform(ptr)])
        b = np.column_stack([np.ones(len(pte)), cx.transform(pte)])
        a, b, _ = kern.full_rank_design(a, b)
        for fname, value in [("K_train", ktt), ("K_test_train", kvt), ("X_train", a), ("X_test", b)]:
            np.savetxt(out / f"{fname}.csv.gz", value, delimiter=",", fmt="%.17g")
        train.to_csv(out / "train.csv", index=False)
        pd.DataFrame({"SampleID": [f"{name}.{c}" for c in pte.cow], "EID": "CRT23_" + pte.cow,
                      "study": "CRT23", "y_obs": pte.y.to_numpy()}).to_csv(out / "test.csv", index=False)
        print(name, "kernel exported:", ktt.shape, kvt.shape, flush=True)
    subprocess.run(["Rscript", str(MP / "04_mixed_model.R"), str(folds), str(results)], check=True,
                   stdout=subprocess.DEVNULL)

    rows, ref = [], {}
    for name in src:
        for pset in ["microbiome", "combined", "covariates"]:
            d = results / "mrm" / pset / name
            m = pd.read_csv(d / "metrics.csv").iloc[0]
            vc = pd.read_csv(d / "variance_components.csv").iloc[0]
            p = pd.read_csv(d / "predictions.csv").y_pred.to_numpy()
            if name == "V4_T2":
                ref[pset] = p
            rows.append({"source": name, "predictor_set": pset, "n": int(m.n), "r_obs": round(m.COR, 3),
                         "bias": round(m.bias_pred_minus_obs, 1), "rmse": round(m.RMSE, 1),
                         "m2": round(vc.m2, 3),
                         "r_vs_V4_T2": None if name == "V4_T2" else round(np.corrcoef(p, ref[pset])[0, 1], 3),
                         "shift_vs_V4_T2": None if name == "V4_T2" else round((p - ref[pset]).mean(), 1)})
    summary = pd.DataFrame(rows)
    summary.to_csv(out_dir / "summary.csv", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
