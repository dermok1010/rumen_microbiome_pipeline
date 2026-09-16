####################################################################
# 02_merge_sheep_methane.R
# sheep_all (EID + samp_date)  ->  nearest-date methane measurement
####################################################################

library(dplyr)
library(readr)

setwd("/home/dermot.kelly/Dermot_analysis/Phd/Paper_2/merging/")

clean_id <- function(x) gsub("[^A-Za-z0-9]", "", as.character(x))

# --- inputs ---
sheep_all <- read_csv("sheep_data/sheep_sample_to_EID.csv",
                      col_types = cols(EID = col_character(),
                                       .default = col_guess()))
methane <- read.csv(
  "/home/dermot.kelly/Dermot_analysis/Phd/PAC_data_pipeline/data/PACfile_ani_id.csv"
) %>%
  mutate(EID = clean_id(EID),
         meth_date = as.Date(date, format = "%d/%m/%Y"))

cat("sheep samples:", nrow(sheep_all),
    "| resolved:", sum(!is.na(sheep_all$EID)), "\n")
cat("methane records:", nrow(methane),
    "| unique animals:", n_distinct(methane$EID), "\n\n")

# --- how many of our animals appear in the methane file at all? ---
sheep_eids <- unique(sheep_all$EID[!is.na(sheep_all$EID)])
cat("sheep animals with any methane record:",
    sum(sheep_eids %in% methane$EID), "of", length(sheep_eids), "\n\n")

# --- nearest-date match, one methane record per microbiome sample ---
sheep_methane <- sheep_all %>%
  filter(!is.na(EID), !is.na(samp_date)) %>%
  inner_join(methane, by = "EID", suffix = c("_micro", "_meth")) %>%
  mutate(date_diff = abs(as.numeric(samp_date - meth_date))) %>%
  group_by(sample_id) %>%
  slice_min(date_diff, n = 1, with_ties = FALSE) %>%
  ungroup()

cat("samples with a methane match:", nrow(sheep_methane),
    "| unique sample_ids:", n_distinct(sheep_methane$sample_id), "\n")
print(summary(sheep_methane$date_diff))


# --- how the match quality distributes, by cohort ---
sheep_methane %>%
  group_by(cohort) %>%
  summarise(n = n(),
            same_day  = sum(date_diff == 0),
            within_7  = sum(date_diff <= 7),
            within_30 = sum(date_diff <= 30),
            median_diff = median(date_diff),
            max_diff    = max(date_diff),
            .groups = "drop") %>%
  print()



sheep_final <- sheep_methane %>% filter(date_diff <= 3)

cat("within 3 days:", nrow(sheep_final),
    "| excluded:", nrow(sheep_methane) - nrow(sheep_final), "\n")

sheep_final %>% count(cohort)

# what got dropped
sheep_methane %>%
  filter(date_diff > 3) %>%
  select(sample_id, cohort, EID, samp_date, meth_date, date_diff) %>%
  arrange(desc(date_diff)) %>%
  print()



clr <- read.csv(
  "/home/dermot.kelly/Dermot_analysis/Phd/Paper_2/rumen_microbiome_pipeline/full_combined_results/CLR_genus_only.csv",
  check.names = FALSE
)

genus_counts <- read.csv(
  "/home/dermot.kelly/Dermot_analysis/Phd/Paper_2/rumen_microbiome_pipeline/full_combined_results/genus_counts_no_controls.csv",
  check.names = FALSE
)

sheep_clr <- sheep_final %>%
  rename(SampleID = sample_id) %>%
  inner_join(clr, by = "SampleID")

sheep_counts <- sheep_final %>%
  rename(SampleID = sample_id) %>%
  inner_join(genus_counts, by = "SampleID")

cat("CLR set:   ", nrow(sheep_clr),  "rows,", ncol(sheep_clr),  "cols\n")
cat("Counts set:", nrow(sheep_counts), "rows,", ncol(sheep_counts), "cols\n")

write_csv(sheep_clr,    "../network_analysis/data/sheep_clr_methane.csv")
write_csv(sheep_counts, "../network_analysis/data/sheep_counts_methane.csv")



library(readr); library(tibble)

# --- ASV table: ASVs = rows, samples = columns ---
asv <- read_tsv(
  "/home/dermot.kelly/Dermot_analysis/Phd/Paper_2/rumen_microbiome_pipeline/full_combined_results/rumen_combined_feature-table.tsv",
  skip = 1, show_col_types = FALSE)          # skip '# Constructed from biom' line
colnames(asv)[1] <- "FeatureID"

# sheep sample IDs in the final (within-3-day) set
sheep_ids <- sheep_final$sample_id
present   <- intersect(sheep_ids, colnames(asv))
cat("sheep samples found in ASV table:", length(present), "of", length(sheep_ids), "\n")

# subset ASV columns to sheep, drop ASVs absent from all sheep samples
sheep_asv_wide <- asv %>% select(FeatureID, all_of(present))
sheep_asv_wide <- sheep_asv_wide[rowSums(sheep_asv_wide[,-1]) > 0, ]
cat("ASVs present in sheep:", nrow(sheep_asv_wide), "\n")

# --- transpose to samples-as-rows, then merge with methane/metadata ---
sheep_asv_long <- sheep_asv_wide %>%
  column_to_rownames("FeatureID") %>%
  t() %>% as.data.frame() %>%
  rownames_to_column("SampleID")

sheep_asv_methane <- sheep_final %>%
  rename(SampleID = sample_id) %>%
  inner_join(sheep_asv_long, by = "SampleID")

cat("ASV+methane set:", nrow(sheep_asv_methane), "rows,",
    ncol(sheep_asv_methane), "cols\n")

# --- write outputs ---
write_tsv(sheep_asv_wide,      "../network_analysis/data/sheep_asv_table.tsv")           # ASVs x samples (raw)
write_csv(sheep_asv_methane,   "../network_analysis/data/sheep_asv_methane.csv")         # samples x (meta+methane+ASVs)





####################################################################
# CLR-transform sheep ASV count data
####################################################################

# sheep_asv_long currently contains:
# SampleID | ASV_1 | ASV_2 | ... | ASV_n

asv_count_matrix <- sheep_asv_long %>%
  column_to_rownames("SampleID") %>%
  as.matrix()

storage.mode(asv_count_matrix) <- "numeric"

cat(
  "\nRaw sheep ASV matrix:",
  nrow(asv_count_matrix), "samples x",
  ncol(asv_count_matrix), "ASVs\n"
)

# ------------------------------------------------------------
# Filter rare ASVs
# Retain ASVs present in at least 5% of sheep samples
# ------------------------------------------------------------

min_prevalence <- 0.05

asv_prevalence <- colMeans(asv_count_matrix > 0)

keep_asv <- asv_prevalence >= min_prevalence

asv_count_filtered <- asv_count_matrix[
  ,
  keep_asv,
  drop = FALSE
]

cat(
  "ASVs retained at",
  min_prevalence * 100,
  "% prevalence:",
  ncol(asv_count_filtered),
  "of",
  ncol(asv_count_matrix),
  "\n"
)

# ------------------------------------------------------------
# CLR transformation
# ------------------------------------------------------------

pseudocount <- 0.5

asv_log <- log(asv_count_filtered + pseudocount)

asv_clr_matrix <- sweep(
  asv_log,
  MARGIN = 1,
  STATS = rowMeans(asv_log),
  FUN = "-"
)

# CLR values within every sample should have mean approximately zero
clr_row_means <- rowMeans(asv_clr_matrix)

cat(
  "Maximum absolute CLR row mean:",
  max(abs(clr_row_means)),
  "\n"
)

# ------------------------------------------------------------
# Convert back to a data frame and attach metadata/methane
# ------------------------------------------------------------

sheep_asv_clr <- as.data.frame(
  asv_clr_matrix,
  check.names = FALSE
) %>%
  rownames_to_column("SampleID")

sheep_asv_clr_methane <- sheep_final %>%
  rename(SampleID = sample_id) %>%
  inner_join(sheep_asv_clr, by = "SampleID")

cat(
  "CLR ASV+methane set:",
  nrow(sheep_asv_clr_methane), "rows,",
  ncol(sheep_asv_clr_methane), "columns\n"
)

# Check for missing or non-finite transformed values
clr_feature_names <- colnames(asv_clr_matrix)

cat(
  "Missing CLR values:",
  sum(is.na(sheep_asv_clr_methane[, clr_feature_names])),
  "\n"
)

cat(
  "Non-finite CLR values:",
  sum(!is.finite(as.matrix(
    sheep_asv_clr_methane[, clr_feature_names]
  ))),
  "\n"
)

# ------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------

write_csv(
  sheep_asv_clr,
  "../network_analysis/data/sheep_asv_clr.csv"
)

write_csv(
  sheep_asv_clr_methane,
  "../network_analysis/data/sheep_asv_clr_methane.csv"
)

cat("\nSaved:\n")
cat("  sheep_asv_clr.csv\n")
cat("  sheep_asv_clr_methane.csv\n")