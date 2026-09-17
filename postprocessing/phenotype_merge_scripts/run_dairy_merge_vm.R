####################################################################
# run_dairy_merge_vm.R
#
# VM-adapted successor to charlie_1.R (never run -- see postprocessing/
# README.md's "phenotype_merge_scripts/" section: all 4 of its inputs
# were missing from this VM until the user uploaded them 2026-09-17 to
# ~/hpc_incoming/dairy_meta/).
#
# Merges the dairy portion of the 17-flowcell 16S genus table
# (postprocessing/output/genus_counts_no_controls_515F_806R.csv) against
# animal phenotype/breeding-value data. Four separate dairy
# sub-populations, four separate join schemes, four output files:
#
#   1. EN00011687 batch, cohorts Dairy_Trt1/2/3 (92 16S samples, "HM24"
#      methane trial) -- joined directly on SampleID against the bridge
#      file HM24 Microbe sample ID.csv (sample_id -> farm_name/run),
#      itself then joined to HM24_cow_data.csv (methane/milk phenotypes)
#      on (farm_name, run). Checked 2026-09-17: bridge sample_ids match
#      the genus table's EN00011687 Dairy_Trt* SampleIDs exactly (92/92).
#      87/92 also match phenotype data on the full (farm_name, run) key --
#      the other 5 (farm 9289 run 1, farm 308 run 3: animal has no rows
#      in HM24_cow_data at all; farm 1329 run 3 x2, farm 75 run 2: that
#      run number doesn't exist for that animal, likely a run-labelling
#      slip in the bridge file) are left unmatched (NA phenotype) rather
#      than guessed at.
#   2. EN00010710 batch ("Heifers", 72 16S samples) -- SampleID suffix
#      (e.g. EN00010710__3064, stripped of a trailing _1/_2 replicate-run
#      suffix where present) matched directly against HF_24_EBI_Profile.xlsx's
#      "FB" column. Checked 2026-09-17: all 69 unique base tags matched
#      (100%); 3 of them (3097, 3108, 3164) have two 16S replicate runs
#      each. Note this file carries EBI/genomic breeding values (Milk,
#      Fert, Carbon, Calv, Beef, Health, ...), not methane phenotypes --
#      it's a different kind of data from the HM24/Tully/sheep methane
#      merges, not a gap in this script.
#
# CRT23_cow_data.csv and clover22_cow_data.csv (also in that upload) ARE
# both part of the 17-flowcell run too -- an earlier pass wrongly reported
# zero overlap because it only checked the EN0001xxxx-numbered dairy
# batches and missed that CRT23/clover22 are their own top-level SampleID
# prefixes in the genus table (CRT23__T{1,2}.<farm_number>,
# clover22__R{1,2,3}.<Farm_Number> -- 2 and 3 within-trial runs per farm,
# matching CRT23_cow_data.csv's 2-level `run` column (Summer/Autumn) and
# clover22_cow_data.csv's 3-level `season` column (Spring/Summer/Autumn)
# exactly). Farm-level match: CRT23 37/41 farms, clover22 22/22 farms.
#
#   3. CRT23 (41 farms x 2 runs = 82 16S samples) -- T1/T2 mapped to
#      Summer/Autumn **chronologically, per the user's confirmation
#      2026-09-17** (no on-VM manifest has real collection dates to
#      verify this independently -- flagged and confirmed, not assumed
#      silently). Joined on (farm_number, run).
#   4. clover22 (22 farms x 3 runs = 66 16S samples) -- R1/R2/R3 mapped to
#      Spring/Summer/Autumn the same way. Joined on (Farm_Number, season).
#
# CRT23 also matches strongly (39/43 by POTTLE_ID/UIDTag) to the RE-RRS
# "Moorepark Dairy" cohort in
# ~/microbiome_prediction/data/re_rrs/comparison_16s_vs_rerrs_2024/ (see
# that repo's docs/re_rrs_data.md) -- i.e. most of these cows now have
# both 16S (here) and RE-RRS data.
####################################################################

library(dplyr)
library(readr)
library(readxl)

genus <- read_csv("output/genus_counts_no_controls_515F_806R.csv",
                   show_col_types = FALSE)

# ------------------------------------------------------------
# 1. HM24 (EN00011687 Dairy_Trt1/2/3) x HM24_cow_data.csv
# ------------------------------------------------------------

hm24_dairy_genus <- genus %>%
  filter(grepl("^EN00011687__Dairy_Trt", SampleID))
cat("HM24 dairy 16S samples:", nrow(hm24_dairy_genus), "\n")

bridge <- read_csv("/home/dermodkkelly/hpc_incoming/dairy_meta/HM24 Microbe sample ID.csv",
                    show_col_types = FALSE) %>%
  mutate(farm_name = as.character(farm_name), run = as.character(run))

stopifnot(setequal(bridge$sample_id, hm24_dairy_genus$SampleID))

hm24_pheno <- read_csv("/home/dermodkkelly/hpc_incoming/dairy_meta/HM24_cow_data.csv",
                        show_col_types = FALSE) %>%
  mutate(farm_name = as.character(farm_name), run = as.character(run))

cat("bridge rows:", nrow(bridge), "\n")
cat("bridge (farm_name, run) keys matched in HM24_cow_data:",
    sum(paste(bridge$farm_name, bridge$run) %in%
          paste(hm24_pheno$farm_name, hm24_pheno$run)),
    "of", n_distinct(paste(bridge$farm_name, bridge$run)), "\n")

hm24_full <- hm24_dairy_genus %>%
  left_join(bridge %>% select(sample_id, batch, cohort, farm_name, run),
            by = c("SampleID" = "sample_id")) %>%
  left_join(hm24_pheno, by = c("farm_name", "run"))

cat("HM24 merged rows:", nrow(hm24_full), "\n")
cat("HM24 rows with a phenotype match (non-NA ch4):",
    sum(!is.na(hm24_full$ch4)), "\n")

write_csv(hm24_full, "output/dairy_hm24_full_clean.csv")

# ------------------------------------------------------------
# 2. Heifers (EN00010710) x HF_24_EBI_Profile.xlsx
# ------------------------------------------------------------

heifer_genus <- genus %>%
  filter(grepl("^EN00010710__", SampleID)) %>%
  mutate(FB = sub("^EN00010710__", "", SampleID),
         FB = sub("_[12]$", "", FB))  # strip replicate-run suffix
cat("Heifer (EN00010710) 16S samples:", nrow(heifer_genus), "\n")
cat("Unique heifer FB tags:", n_distinct(heifer_genus$FB), "\n")

ebi <- read_excel("/home/dermodkkelly/hpc_incoming/dairy_meta/HF_24_EBI_Profile.xlsx") %>%
  mutate(FB = as.character(FB))

cat("EBI profile animals:", n_distinct(ebi$FB), "\n")
cat("Heifer FB tags matched in EBI profile:",
    sum(unique(heifer_genus$FB) %in% ebi$FB), "of",
    n_distinct(heifer_genus$FB), "\n")

heifer_full <- heifer_genus %>%
  left_join(ebi, by = "FB")

cat("Heifer merged rows:", nrow(heifer_full), "\n")
cat("Heifer rows with an EBI match:", sum(!is.na(heifer_full$EBI)), "\n")

write_csv(heifer_full, "output/dairy_heifers_ebi_clean.csv")

# ------------------------------------------------------------
# 3. CRT23 (T1/T2 -> Summer/Autumn) x CRT23_cow_data.csv
# ------------------------------------------------------------

run_map_crt23 <- c(T1 = "Summer", T2 = "Autumn")

crt23_genus <- genus %>%
  filter(grepl("^CRT23__", SampleID)) %>%
  mutate(
    run_tag     = sub("^CRT23__(T[12])\\..*$", "\\1", SampleID),
    farm_number = sub("^CRT23__T[12]\\.", "", SampleID),
    run         = unname(run_map_crt23[run_tag])
  )
cat("CRT23 16S samples:", nrow(crt23_genus), "\n")

crt23_pheno <- read_csv("/home/dermodkkelly/hpc_incoming/dairy_meta/CRT23_cow_data.csv",
                         show_col_types = FALSE) %>%
  mutate(farm_number = as.character(farm_number))

cat("CRT23 (farm_number, run) matched:",
    sum(paste(crt23_genus$farm_number, crt23_genus$run) %in%
          paste(crt23_pheno$farm_number, crt23_pheno$run)),
    "of", nrow(crt23_genus), "\n")

crt23_full <- crt23_genus %>%
  left_join(crt23_pheno, by = c("farm_number", "run"))

cat("CRT23 merged rows:", nrow(crt23_full), "\n")
cat("CRT23 rows with a phenotype match (non-NA ch4mean):",
    sum(!is.na(crt23_full$ch4mean)), "\n")

write_csv(crt23_full, "output/dairy_crt23_full_clean.csv")

# ------------------------------------------------------------
# 4. clover22 (R1/R2/R3 -> Spring/Summer/Autumn) x clover22_cow_data.csv
# ------------------------------------------------------------

run_map_clover22 <- c(R1 = "Spring", R2 = "Summer", R3 = "Autumn")

clover22_genus <- genus %>%
  filter(grepl("^clover22__", SampleID)) %>%
  mutate(
    run_tag     = sub("^clover22__(R[123])\\..*$", "\\1", SampleID),
    Farm_Number = sub("^clover22__R[123]\\.", "", SampleID),
    season      = unname(run_map_clover22[run_tag])
  )
cat("clover22 16S samples:", nrow(clover22_genus), "\n")

clover22_pheno <- read_csv("/home/dermodkkelly/hpc_incoming/dairy_meta/clover22_cow_data.csv",
                            show_col_types = FALSE) %>%
  mutate(Farm_Number = as.character(Farm_Number))

cat("clover22 (Farm_Number, season) matched:",
    sum(paste(clover22_genus$Farm_Number, clover22_genus$season) %in%
          paste(clover22_pheno$Farm_Number, clover22_pheno$season)),
    "of", nrow(clover22_genus), "\n")

clover22_full <- clover22_genus %>%
  left_join(clover22_pheno, by = c("Farm_Number", "season"))

cat("clover22 merged rows:", nrow(clover22_full), "\n")
cat("clover22 rows with a phenotype match (non-NA CH4mean):",
    sum(!is.na(clover22_full$CH4mean)), "\n")

write_csv(clover22_full, "output/dairy_clover22_full_clean.csv")
