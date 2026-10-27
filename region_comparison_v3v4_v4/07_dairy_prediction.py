"""Does a methane model trained on V4 data still work when fed V3-V4 data from the same cows?

Train on the other dairy trials' V4 genus profiles (HM24, clover22: GreenFeed CH4, g/day), then
predict the CRT23 cows' Autumn methane from four profiles of each cow:
  V4_T2       V4, Autumn sample: the reference
  V4_T1       V4, Summer sample of the same cow: how much predictions move with biology alone
  V3V4        BSI's V3-V4 sample, its own primer-matched classifier
  V3V4_core   the same reads trimmed to the V4 core and classified with the V4 classifier
BSI's aliquot is closest to T2 in all 39 cows (timepoint_matching.csv), so V3V4 vs V4_T2 is
mostly the region/lab difference.

Features: CLR of genera with mean relative abundance >= 0.01% in training (pseudocount 0.5).
Transforms: "raw"; "centred", where each dataset's mean CLR profile is subtracted (each training
trial separately, and each test source across its own cows: unsupervised, no methane used);
"ranked" (RF only), per-feature percentile ranks within each dataset.
Reported per model x transform x source: r and bias against observed methane, and agreement
(r, mean difference) with the V4_T2 predictions.

Run from the repo root with ~/envs/micro_prediction/bin/python, after 06_dairy_classify.sh.
"""
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import RidgeCV
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "outputs" / "region_comparison_v3v4_v4"
RES = OUT / "prediction"
META = Path.home() / "hpc_incoming" / "dairy_meta"
V4_TABLE = REPO / "results" / "rumen_combined_feature-table.tsv"
MIN_REL, PSEUDO, SEED = 1e-4, 0.5, 42
COVS = ["bwt", "yield", "parity"]


def genus_of(taxon):
    for rank in str(taxon).split(";"):
        rank = rank.strip()
        if rank.startswith("g__") and len(rank) > 3:
            return rank[3:]
    return "Unassigned_genus"


def training_genus():
    cols = pd.read_csv(V4_TABLE, sep="\t", skiprows=1, nrows=0).columns
    sel = [c for c in cols if c.startswith(("EN00011687__Dairy_Trt", "clover22__"))]
    t = pd.read_csv(V4_TABLE, sep="\t", skiprows=1, usecols=[cols[0]] + sel, index_col=0)
    tax = pd.concat([pd.read_csv(OUT / f, sep="\t", index_col=0)["Taxon"] for f in
                     ["taxonomy_V4_gtdb_r214.tsv", "taxonomy_V4_dairy_train_gtdb_r214.tsv"]])
    tax = tax[~tax.index.duplicated()]
    t = t[(t > 0).any(axis=1)]
    missing = (~t.index.isin(tax.index)).sum()
    assert missing == 0, f"{missing} training ASVs without taxonomy"
    return t.groupby(t.index.map(lambda a: genus_of(tax[a]))).sum().T


def phenotypes():
    c = pd.read_csv(META / "CRT23_cow_data.csv")
    c = pd.DataFrame({"key": "CRT23|" + c.farm_number.astype(str) + "|" + c.run, "y": c.ch4adj,
                      "bwt": c.BWT, "yield": c.daily_yield, "parity": c.parity})
    h = pd.read_csv(META / "HM24_cow_data.csv")
    h = pd.DataFrame({"key": "HM24|" + h.farm_name.astype(str) + "|" + h.run.astype(str), "y": h.ch4,
                      "bwt": h.bwt, "yield": h.daily_yield, "parity": h.parity})
    v = pd.read_csv(META / "clover22_cow_data.csv")
    v = pd.DataFrame({"key": "clover22|" + v.Farm_Number.astype(str) + "|" + v.season, "y": v.CH4adj,
                      "bwt": v.BWT, "yield": v.DAILY_YIELD, "parity": v.PARITY})
    return pd.concat([c, h, v]).set_index("key")


def training_set(pheno):
    g = training_genus()
    bridge = pd.read_csv(META / "HM24 Microbe sample ID.csv")
    hm = dict(zip(bridge.sample_id, "HM24|" + bridge.farm_name.astype(str) + "|" + bridge.run.astype(str)))
    season = {"R1": "Spring", "R2": "Summer", "R3": "Autumn"}

    def key(s):
        if s in hm:
            return hm[s]
        m = re.match(r"clover22__(R\d)\.(\d+)", s)
        return f"clover22|{m.group(2)}|{season[m.group(1)]}"

    keys = g.index.map(key)
    ok = keys.isin(pheno.index)
    g, keys = g[ok], keys[ok]
    p = pheno.loc[keys].reset_index(drop=True)
    p["dataset"] = [k.split("|")[0] for k in keys]
    p["sample"] = g.index
    p = p.dropna(subset=["y"] + COVS)
    return g.loc[p["sample"]].reset_index(drop=True), p.reset_index(drop=True)


def test_sources(pheno):
    """Genus counts per CRT23 cow for each source, and the cows' Autumn phenotype."""
    v4 = pd.read_csv(OUT / "genus_V4_gtdb_r214.tsv", sep="\t", index_col=0)
    cows = pd.read_csv(OUT / "samples.csv", dtype={"cow": str})
    cow = dict(zip(cows["sample"], cows["cow"]))

    def by_cow(t):  # V3-V4 replicate library preps are summed per cow
        return t.groupby(t.index.astype(str).map(lambda s: cow[s])).sum()

    v3v4 = by_cow(pd.read_csv(OUT / "genus_V3V4_gtdb_r214.tsv", sep="\t", index_col=0))
    core = by_cow(pd.read_csv(OUT / "genus_V3V4_core_gtdb_r214.tsv", sep="\t", index_col=0))
    src = {"V4_T2": by_cow(v4[v4.index.str.startswith("CRT23__T2.")]),
           "V4_T1": by_cow(v4[v4.index.str.startswith("CRT23__T1.")]),
           "V3V4": v3v4, "V3V4_core": core}
    keys = pd.Index([f"CRT23|{c}|Autumn" for c in v3v4.index])
    have = keys.isin(pheno.index)
    cows_ok = [c for c, h in zip(v3v4.index, have) if h and all(c in s.index for s in src.values())]
    p = pheno.loc[[f"CRT23|{c}|Autumn" for c in cows_ok]].reset_index(drop=True)
    p["cow"] = cows_ok
    return {k: s.loc[cows_ok] for k, s in src.items()}, p


def clr(counts, genera):
    x = counts.reindex(columns=genera, fill_value=0).to_numpy(float) + PSEUDO
    lx = np.log(x)
    return lx - lx.mean(axis=1, keepdims=True)


def transform(x, groups, how):
    x = x.copy()
    for g in np.unique(groups):
        m = groups == g
        if how == "centred":
            x[m] -= x[m].mean(axis=0)
        elif how == "ranked":
            x[m] = stats.rankdata(x[m], axis=0) / m.sum()
    return x


def models():
    return {"ridge": lambda: make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(-2, 5, 30))),
            "rf": lambda: RandomForestRegressor(n_estimators=1000, min_samples_leaf=3, max_features=0.33,
                                                n_jobs=4, random_state=SEED)}


def main():
    RES.mkdir(exist_ok=True)
    pheno = phenotypes()
    gtr, ptr = training_set(pheno)
    rel = gtr.div(gtr.sum(axis=1), axis=0)
    genera = [g for g in rel.columns[rel.mean() >= MIN_REL] if g != "Unassigned_genus"]
    src, pte = test_sources(pheno)
    print(f"train: {len(ptr)} samples ({ptr.dataset.value_counts().to_dict()}), {len(genera)} genera; "
          f"test: {len(pte)} CRT23 cows, Autumn CH4 mean {pte.y.mean():.0f} g/d (SD {pte.y.std():.0f})")
    print("test genera present in training set (% of reads): " + ", ".join(
        f"{k} {100 * s.reindex(columns=genera, fill_value=0).to_numpy().sum() / s.to_numpy().sum():.1f}"
        for k, s in src.items()))

    xtr_clr = clr(gtr, genera)
    ctr = ptr[COVS].to_numpy(float)
    rows, preds = [], []
    for how in ["raw", "centred", "ranked"]:
        xtr = transform(xtr_clr, ptr.dataset.to_numpy(), how)
        for name, make in models().items():
            if how == "ranked" and name == "ridge":
                continue
            for pset in ["microbiome", "combined"]:
                a = xtr if pset == "microbiome" else np.column_stack([xtr, ctr])
                m = make().fit(a, ptr.y)
                for s, counts in src.items():
                    xte = transform(clr(counts, genera), np.zeros(len(counts)), how)
                    b = xte if pset == "microbiome" else np.column_stack([xte, pte[COVS].to_numpy(float)])
                    preds.append(pd.DataFrame({"model": name, "transform": how, "predictor_set": pset,
                                               "source": s, "cow": pte.cow, "y": pte.y, "y_pred": m.predict(b)}))
    # covariates alone: the baseline any microbiome model has to beat; the same for every source
    for name, make in models().items():
        m = make().fit(ctr, ptr.y)
        preds.append(pd.DataFrame({"model": name, "transform": "-", "predictor_set": "covariates",
                                   "source": "-", "cow": pte.cow, "y": pte.y,
                                   "y_pred": m.predict(pte[COVS].to_numpy(float))}))
    preds = pd.concat(preds, ignore_index=True)
    preds.to_csv(RES / "predictions.csv", index=False)

    ref = preds[preds.source == "V4_T2"].set_index(["model", "transform", "predictor_set", "cow"]).y_pred
    for (mo, tr, ps, s), d in preds.groupby(["model", "transform", "predictor_set", "source"], sort=False):
        r, p = stats.pearsonr(d.y_pred, d.y)
        row = {"model": mo, "transform": tr, "predictor_set": ps, "source": s, "n": len(d),
               "r_obs": round(r, 3), "p_obs": round(p, 3), "bias": round((d.y_pred - d.y).mean(), 1),
               "rmse": round(np.sqrt(((d.y_pred - d.y) ** 2).mean()), 1)}
        if s not in ("-", "V4_T2"):
            rp = ref.loc[(mo, tr, ps)].reindex(d.cow).to_numpy()
            row["r_vs_V4_T2"] = round(stats.pearsonr(d.y_pred, rp)[0], 3)
            row["shift_vs_V4_T2"] = round((d.y_pred.to_numpy() - rp).mean(), 1)
        rows.append(row)
    summary = pd.DataFrame(rows)
    summary.to_csv(RES / "summary.csv", index=False)
    with pd.option_context("display.width", 200, "display.max_rows", 200):
        print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
