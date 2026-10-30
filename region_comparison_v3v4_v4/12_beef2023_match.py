"""Are rumen-core's unlabelled Beef_2023 V4 samples the same rumen samples as BSI's Grange V3-V4 beef?

Nothing links Beef_2023's sample numbers to animal tags, so this matches profiles instead. Every
sample is cut to the shared V4 core (10_asv_core_tables.py's rule), and each BSI V3-V4 library is
matched to its nearest V4 sample (Bray-Curtis on relative abundance, and Aitchison on CLR of cores
in >= 10% of samples).

Calibration: BSI's Moorepark dairy libraries against all V4 dairy samples (CRT23 T1/T2 plus the
HM24/clover22 herds as decoys). Their true partner is known (same cow), so the share of libraries
whose nearest V4 sample is their own cow shows how well matching works across primer sets.
Test: BSI's Grange beef libraries against all V4 beef (Beef_2023, plus Beef_115/117 and Tully as
decoys). Same-sample pairs should show as Grange libraries whose nearest V4 sample is a Beef_2023
sample, one-to-one (mutual nearest neighbours), at distances like the dairy true pairs. Herd-level
closeness alone would also fit Grange and Beef_2023 merely sharing a farm or diet.

Writes outputs/region_comparison_v3v4_v4/beef2023_match/{nearest.csv,summary.csv}.
Run from the repo root with ~/envs/rumen_microbiome_pipeline/bin/python (needs h5py).
"""
import io
import re
import sys
import zipfile
from importlib import import_module
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(Path.home() / "microbiome_prediction" / "re_rrs_comparison")]
gt = import_module("02_genus_tables")
nzi = import_module("build_nz_inputs")
act = import_module("10_asv_core_tables")
OUT = gt.OUT / "beef2023_match"
V4_GROUPS = {"Beef_2023": "Beef_2023", "Beef_115": "Beef_115", "Beef_117": "Beef_117", "Tully__": "Tully",
             "CRT23__": "CRT23", "clover22__": "dairy_other", "EN00011687__Dairy_Trt": "dairy_other"}


def v4_group(sample):
    return next((g for k, g in V4_GROUPS.items() if k in sample), None)


def bsi_cattle():
    """BSI V3-V4 cattle libraries (ASVs x libraries) and library -> (herd, animal tag)."""
    info = pd.read_csv(gt.CMP / "sampleInfo.csv", dtype=str)
    cattle = info[info.Species == "CATTLE"]
    herd = {lab: ("Grange" if c.startswith("Grange") else "Moorepark", tag)
            for lab, c, tag in zip(cattle.LabSampleID, cattle.COMMENT.fillna(""), cattle.UIDTag)}
    z = zipfile.ZipFile(gt.V3V4_TABLE)
    name = next(n for n in z.namelist() if n.endswith("feature-table.biom"))
    with h5py.File(io.BytesIO(z.read(name)), "r") as f:
        obs = [x.decode() for x in f["observation/ids"][:]]
        smp = [x.decode() for x in f["sample/ids"][:]]
        data, ind, ptr = f["sample/matrix/data"][:], f["sample/matrix/indices"][:], f["sample/matrix/indptr"][:]
    rows, meta = {}, {}
    for i, s in enumerate(smp):  # replicate library preps carry a trailing letter; keep them
        base = re.sub(r"[a-z]+$", "", s)
        if base in herd:
            rows[s] = pd.Series(data[ptr[i]:ptr[i + 1]], index=[obs[j] for j in ind[ptr[i]:ptr[i + 1]]])
            meta[s] = herd[base]
    return pd.DataFrame(rows).fillna(0), meta


def main():
    OUT.mkdir(exist_ok=True)
    cols = pd.read_csv(gt.V4_TABLE, sep="\t", skiprows=1, nrows=0).columns
    sel = [c for c in cols if v4_group(c)]
    v4 = pd.read_csv(gt.V4_TABLE, sep="\t", skiprows=1, usecols=[cols[0]] + sel, index_col=0)
    v4 = v4[(v4 > 0).any(axis=1)]
    v34, meta = bsi_cattle()
    v34 = v34[(v34 > 0).any(axis=1)]
    # The VM's rep-seqs are from the 2026-09-15 run; the table is the 2026-10-01 run, which redid the
    # three Tully flowcells, so some (mostly Tully) ASVs have no sequence here and are dropped.
    seqs = act.fasta(nzi.IE_SEQS)
    has = v4.index.isin(list(seqs))
    kept = (v4[has].sum() / v4.sum()).groupby(lambda s: v4_group(s))
    print("V4 reads with a known sequence, median by group:", kept.median().round(3).to_dict())
    v4 = v4[has]
    a, _, _ = act.core_table(v4, seqs, act.v4_core)
    b, _, _ = act.core_table(v34, act.fasta(nzi.NZ_SEQS), nzi.nz_core)
    # libraries and samples below 5,000 core reads are too thin to match on
    a, b = a[a.sum(axis=1) >= 5000], b[b.sum(axis=1) >= 5000]
    cores = a.columns.union(b.columns)
    a, b = a.reindex(columns=cores, fill_value=0), b.reindex(columns=cores, fill_value=0)
    ra, rb = a.div(a.sum(axis=1), axis=0), b.div(b.sum(axis=1), axis=0)
    common = cores[(pd.concat([a, b]) > 0).mean() >= 0.10]

    def clr(m):
        x = np.log(m[common].to_numpy(float) + 0.5)
        return x - x.mean(axis=1, keepdims=True)

    # rows of both sum to 1, so Bray-Curtis is half the L1 distance; numpy only (no scipy in this env)
    xa, xb, ca, cb = ra.to_numpy(), rb.to_numpy(), clr(a), clr(b)
    dist = {"bray_curtis": np.array([np.abs(xa - r).sum(axis=1) / 2 for r in xb]),
            "aitchison": np.array([np.sqrt(((ca - r) ** 2).sum(axis=1)) for r in cb])}
    groups = np.array([v4_group(s) for s in a.index])
    rows = []
    for metric, D in dist.items():
        for herd, pool in [("Moorepark", ["CRT23", "dairy_other"]),
                           ("Grange", ["Beef_2023", "Beef_115", "Beef_117", "Tully"])]:
            cand = np.isin(groups, pool)
            libs = [i for i, s in enumerate(b.index) if meta[s][0] == herd]
            sub = D[np.ix_(libs, np.where(cand)[0])]
            cand_ids = a.index[cand]
            nearest_lib = {}  # V4 sample -> its nearest BSI library, for mutual-nearest checks
            for j in range(sub.shape[1]):
                nearest_lib[cand_ids[j]] = b.index[libs[int(np.argmin(sub[:, j]))]]
            for k, i in enumerate(libs):
                order = np.argsort(sub[k])
                nn, d1, d2 = cand_ids[order[0]], sub[k, order[0]], sub[k, order[1]]
                tag = meta[b.index[i]][1]
                own = [j for j, s in enumerate(cand_ids) if v4_group(s) == "CRT23" and re.search(rf"\.{tag}$", s)]
                rows.append({"metric": metric, "herd": herd, "library": b.index[i], "tag": tag,
                             "nearest": nn, "nearest_group": v4_group(nn), "d_nearest": d1,
                             "d_second": d2, "mutual": nearest_lib[nn] == b.index[i],
                             "own_rank": (int(min(np.where(np.isin(order, own))[0])) + 1) if own else None,
                             "d_own": float(sub[k, own].min()) if own else None})
    near = pd.DataFrame(rows)
    near.to_csv(OUT / "nearest.csv", index=False)
    summ = []
    for (metric, herd), g in near.groupby(["metric", "herd"]):
        r = {"metric": metric, "herd": herd, "libraries": len(g),
             "median_d_nearest": round(g.d_nearest.median(), 3), "mutual_nn": int(g.mutual.sum())}
        r.update({f"nn_{k}": int(v) for k, v in g.nearest_group.value_counts().items()})
        if herd == "Moorepark":
            r.update({"own_cow_is_nearest": int((g.own_rank == 1).sum()),
                      "own_cow_in_top5": int((g.own_rank <= 5).sum()),
                      "median_d_own": round(g.d_own.median(), 3)})
        else:
            m = g[(g.nearest_group == "Beef_2023") & g.mutual]
            r.update({"mutual_nn_beef2023": len(m), "beef2023_samples_claimed": m.nearest.nunique()})
        summ.append(r)
    summ = pd.DataFrame(summ)
    summ.to_csv(OUT / "summary.csv", index=False)
    print(f"V4 samples: {pd.Series(groups).value_counts().to_dict()}; BSI libraries: "
          f"{pd.Series([meta[s][0] for s in b.index]).value_counts().to_dict()}; "
          f"cores {len(cores)}, CLR on {len(common)}")
    print(summ.T.to_string())


if __name__ == "__main__":
    main()
