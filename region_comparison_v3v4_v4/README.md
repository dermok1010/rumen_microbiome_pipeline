# V3-V4 vs V4 16S on the same dairy cows

**Question.** Are V3-V4 (341F/805R) 16S results directly comparable with this pipeline's V4
(515F/806R) results? Asked by Dermot, 2026-10-27, as a check before treating the two amplicons as
interchangeable.

**Animals.** 39 Moorepark dairy cows from the CRT23 trial that appear in two 16S datasets:
- **V3-V4**: BSI's `comparison_16s_vs_rerrs_2024` cattle (43 rows, 39 cows), sequenced by BSI and
  processed by Dermot on the Teagasc HPC on 30 Aug 2026 with this pipeline's Docker image
  (QIIME2 2024.2; DADA2 trim_left 17/21, trunc_len 280/260).
- **V4**: `CRT23__T1.<cow>` and `CRT23__T2.<cow>` in this pipeline's combined run (Summer and
  Autumn 2023; 82 samples, 41 cows).

Whether BSI's aliquot is the same rumen sample as T1 or T2 is not recorded (BSI's cattle
SAMPLEDATE is blank); the analysis checks which timepoint each cow's V3-V4 profile is closest to.

**Taxonomy.** Both regions classified on GTDB r214 with naive-Bayes classifiers trained on the
matching primer extract (`classifier/build/gtdb_r214/02_extract_train.sh`, `REGION=V4` / `V3V4`),
so taxonomy differences reflect the region, not the database.

## Steps

    bash region_comparison_v3v4_v4/01_classify.sh     # GTDB r214 taxonomy for both regions' ASVs
    python region_comparison_v3v4_v4/02_genus_tables.py

Outputs: `outputs/region_comparison_v3v4_v4/` (gitignored).

## Prediction test (2026-10-27)

    bash region_comparison_v3v4_v4/06_dairy_classify.sh   # GTDB r214 for the training herds' V4 ASVs
    python region_comparison_v3v4_v4/07_dairy_prediction.py

Ridge and RF on genus CLR trained on HM24 + clover22 V4 (153 samples, GreenFeed CH4), predicting
Autumn CH4 of 35 CRT23 cows from their V4 T2, V4 T1, V3-V4 and V3-V4-trimmed-to-V4-core profiles.
Outputs in `outputs/region_comparison_v3v4_v4/prediction/`. Result: V3-V4 inputs gave the same or
slightly higher r than V4 (microbiome only 0.34-0.46 vs 0.12-0.34; differences within noise at
n = 35), and V3-V4 predictions agreed with the V4 T2 ones more closely than the same cow's V4 T1
predictions did. The region's main effect was a level shift of about -20 to -26 g/day in raw-CLR
ridge, which per-dataset centring removes. Microbiome-only skill across herds is weak even on V4,
so the test can only rule out a large accuracy loss, not a small one.

## Prediction test with the paper's kernel model (2026-10-29)

    ~/envs/micro_prediction/bin/python region_comparison_v3v4_v4/09_dairy_kernel_prediction.py

The 07 question refitted with the 16S vs RE-RRS paper's REML kernel model, using its code unchanged
from `~/microbiome_prediction` (`03_export_kernels.py` kernel and design, `common.py` CLR and
covariate transforms, `04_mixed_model.R` fit). Same 153 training samples, 35 test cows and four test
profiles as 07; features per the paper's genus configs (genera in >= 5% of training samples, 385;
CLR, pseudocount 0.5); no per-dataset centring. Pearson r, microbiome only: V4 T2 0.30, V4 T1 0.19,
V3-V4 0.24, V3-V4 core 0.31; combined 0.61 / 0.59 / 0.59 / 0.60; covariates alone 0.39. V3-V4 input
shifts predictions by -11.5 g/day on average against V4 T2, about the same as the cows' own summer
V4 sample (-9.7). This replaces ridge/RF as the result reported to the team. Outputs in
`outputs/region_comparison_v3v4_v4/prediction_kernel/`.

## Prediction test at ASV level on the shared V4 core (2026-10-29)

    ~/envs/rumen_microbiome_pipeline/bin/python region_comparison_v3v4_v4/10_asv_core_tables.py
    ~/envs/micro_prediction/bin/python region_comparison_v3v4_v4/11_dairy_kernel_asv.py

The test the paper relies on across primer sets: features are the first 250 bp after 515F
(`loco_16s_v4core`'s rule, `build_nz_inputs.py`), so the V4-trained kernel reads V3-V4 test samples
sequence for sequence. V4 ASVs from some dairy flow cells (CRT23, part of HM24/clover22) keep 8-12
bases of the 515F 3' end; 10 removes them before cutting, which the paper's `ie_core` does not need
for the Irish sheep. 2,746 cores in >= 5% of training samples, holding 93.8% of the test cows' V4
reads and 89.3% of their V3-V4 reads. Pearson r, same 35 cows, V4 vs V3-V4 test set: microbiome
0.42 vs 0.34, combined 0.62 vs 0.62; V3-V4 shifts microbiome predictions by -2.3 g/day. This is the
result reported to the team; the genus-level version (09) is kept for reference. Outputs in
`outputs/region_comparison_v3v4_v4/{asv_core,prediction_kernel_asv}/`.

## Depth check (2026-10-29)

    micromamba run -p ~/envs/micro_prediction Rscript region_comparison_v3v4_v4/08_depth_check.R

BSI's V3-V4 libraries are about 2.6x deeper than the CRT23 V4 samples (median 111k vs 44k reads on
the shared genera; zero cells 3.3% vs 16.7%). Rebuilding the Aitchison PCoA with every sample
rarefied to 18,237 reads, or on proportions, still splits the primer sets completely on PC1
(27.8-28.2% and 32.1%, against 26.5% on raw counts), so the separation is not a depth artefact.
Output: `outputs/region_comparison_v3v4_v4/depth_check.csv`.

## Beef_2023 = Grange Beef_Diet_PRG: profile matching (2026-10-30)

    ~/envs/rumen_microbiome_pipeline/bin/python region_comparison_v3v4_v4/12_beef2023_match.py

rumen-core's 70 `Beef_2023` V4 samples have no animal-ID key. Every BSI V3-V4 cattle library and
every V4 cattle sample was cut to the shared V4 core (10's rule) and each library matched to its
nearest V4 sample. Calibration on the 43 Moorepark dairy libraries against 240 V4 dairy samples:
own cow nearest for 39 (Aitchison) and 30 (Bray-Curtis). Result:

- All 75 Grange `Beef_Diet_PRG` libraries have a `Beef_2023` sample as nearest; all 70 `Beef_2023`
  samples are the mutual nearest neighbour of a different Diet PRG library (Bray-Curtis), and
  Aitchison picks the same partner for all 70.
- Pair distances match known same-cow pairs: Bray-Curtis median 0.35 (dairy own cow 0.37), against
  0.66 for the 46 Grange `Beef_Treatment` GC/MSS libraries, which have no V4 partner and match
  mostly Tully (32) or loosely Beef_2023 (14).

So `Beef_2023` (`EN00011686` 27-30 + `EN00011687` 35-40 samples) is very likely the same rumen
sampling of the Grange Diet PRG suckler cows (75 cows, 3-digit numbers, born 2015-2021), 70 of 75
sequenced by Macrogen. Not settled: same extract or same sampling day vs a nearby one (the dairy
true pairs also span a few weeks). The one-to-one pairs in `nearest.csv` (metric bray_curtis,
herd Grange, mutual, nearest_group Beef_2023) give a candidate sample-to-cow key; confirm with
the Macrogen sheet (asked of Wasim, 30 Oct) before joining phenotypes on it. The VM's rep-seqs
are from the 2026-09-15 run, so a few Tully ASVs from the 2026-10-01 rerun are dropped (Tully
keeps 99.8% of reads; it is only a decoy here). This would give beef a same-animal V4 vs V3-V4
pairing like the dairy cows, once Grange phenotypes exist. Outputs:
`outputs/region_comparison_v3v4_v4/beef2023_match/{nearest.csv,summary.csv}`.
