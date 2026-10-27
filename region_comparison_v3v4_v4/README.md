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
