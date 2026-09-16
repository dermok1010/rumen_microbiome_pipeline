

library(dplyr)
library(readr)

setwd("/home/dermot.kelly/Dermot_analysis/Phd/Paper_2/merging/dairy_metadata/")  # or wherever you've put them

# --- load ---
bridge  <- read_csv("HM24 Microbe sample ID.csv", show_col_types = FALSE)
hm24    <- read_csv("HM24_cow_data.csv", show_col_types = FALSE)
crt23   <- read_csv("CRT23_cow_data.csv", show_col_types = FALSE)
clover22<- read_csv("clover22_cow_data.csv", show_col_types = FALSE)

# --- the farm IDs the bridge needs to resolve ---
bridge_farms <- unique(as.character(bridge$farm_name))
cat("Bridge references", length(bridge_farms), "unique farm_names\n\n")

# --- farm ID columns in each cow-data sheet (note different column names) ---
hm24_farms    <- unique(as.character(hm24$farm_name))
crt23_farms   <- unique(as.character(crt23$farm_number))
clover22_farms<- unique(as.character(clover22$Farm_Number))

# --- for each bridge farm, which sheets contain it? ---
match_check <- data.frame(farm = bridge_farms) %>%
  mutate(
    in_HM24     = farm %in% hm24_farms,
    in_CRT23    = farm %in% crt23_farms,
    in_clover22 = farm %in% clover22_farms,
    n_sheets    = in_HM24 + in_CRT23 + in_clover22
  )

# summary
cat("How many sheets does each bridge farm match?\n")
print(table(match_check$n_sheets))

cat("\nFarms matching MULTIPLE sheets (ambiguous):\n")
match_check %>% filter(n_sheets > 1) %>% print()

cat("\nFarms matching NO sheet (unresolved):\n")
match_check %>% filter(n_sheets == 0) %>% print()



# what run values are in each?
bridge %>% count(run)
hm24 %>% count(run)

# for one cow, do the runs correspond?
bridge %>% filter(farm_name == 606) %>% select(sample_id, cohort, farm_name, run) %>% arrange(run)
hm24   %>% filter(farm_name == 606) %>% select(farm_name, run, run_date, ch4) %>% arrange(run)

dairy_hm24 <- bridge %>%
  left_join(hm24, by = c("farm_name", "run"))

# verify it stayed 1:1 (no explosion, no warning)
nrow(dairy_hm24)              # expect 92
sum(is.na(dairy_hm24$ch4))    # how many samples got no methane match


dairy_hm24 %>%
  filter(is.na(ch4)) %>%
  select(sample_id, cohort, farm_name, run)


clr <- read.csv("/home/dermot.kelly/Dermot_analysis/Phd/Paper_2/rumen_microbiome_pipeline/full_combined_results/CLR_genus_only.csv", check.names = FALSE)

dairy_hm24_full <- dairy_hm24 %>%
  rename(SampleID = sample_id) %>%
  inner_join(clr, by = "SampleID")

nrow(dairy_hm24_full)    # expect 92 (all HM24 samples are in the CLR matrix)

print(head(dairy_hm24_full))


write_csv(dairy_hm24_full, "/home/dermot.kelly/Dermot_analysis/Phd/Paper_2/network_analysis/lactating_cows.csv")





# ============================================================
# Prepare raw genus counts for dairy SPIEC-EASI
# ============================================================
# ============================================================
# Extract dairy genus counts for SPIEC-EASI
# ============================================================

genus_counts_path <- paste0(
  "/home/dermot.kelly/Dermot_analysis/Phd/Paper_2/",
  "rumen_microbiome_pipeline/full_combined_results/",
  "genus_counts_no_controls.csv"
)

genus_counts <- read.csv(
  genus_counts_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("Genus count table dimensions:\n")
print(dim(genus_counts))

cat("\nFirst columns:\n")
print(names(genus_counts)[1:10])

print(genus_counts[1:5, 1:6])


# ============================================================
# Extract dairy genus counts
# ============================================================

dairy_sample_ids <- as.character(dairy_hm24_full$SampleID)

dairy_genus_counts_df <- genus_counts %>%
  filter(SampleID %in% dairy_sample_ids) %>%
  arrange(match(SampleID, dairy_sample_ids))

cat("Dairy samples recovered:\n")
print(nrow(dairy_genus_counts_df))

missing_dairy_ids <- setdiff(
  dairy_sample_ids,
  dairy_genus_counts_df$SampleID
)

cat("Missing dairy samples:\n")
print(missing_dairy_ids)

stopifnot(
  nrow(dairy_genus_counts_df) == 92,
  identical(
    dairy_genus_counts_df$SampleID,
    dairy_sample_ids
  )
)

dairy_genus_count_mat <- dairy_genus_counts_df %>%
  select(-SampleID) %>%
  mutate(
    across(
      everything(),
      ~ suppressWarnings(as.numeric(.x))
    )
  ) %>%
  as.matrix()

rownames(dairy_genus_count_mat) <- dairy_genus_counts_df$SampleID
storage.mode(dairy_genus_count_mat) <- "numeric"

cat("\nInitial dairy genus count matrix:\n")
print(dim(dairy_genus_count_mat))

# ============================================================
# Save dairy raw genus counts
# ============================================================

dairy_genus_counts_path <- paste0(
  "/home/dermot.kelly/Dermot_analysis/Phd/Paper_2/",
  "network_analysis/data/dairy_genus_counts.csv"
)

write_csv(
  dairy_genus_counts_df,
  dairy_genus_counts_path
)

cat(
  "\nSaved dairy genus counts to:\n",
  dairy_genus_counts_path,
  "\n"
)









# ============================================================
# Dairy: phenotypes + FULL unfiltered raw genus counts
# ============================================================

# 1. full raw counts for all biological samples (519 genera, 1459 samples)
genus_counts <- read.csv(
  "/home/dermot.kelly/Dermot_analysis/Phd/Paper_2/rumen_microbiome_pipeline/full_combined_results/genus_counts_no_controls.csv",
  check.names = FALSE, stringsAsFactors = FALSE
)
cat("Full counts:", dim(genus_counts), "\n")

# 2. dairy phenotypes only (no genus columns) — rebuild from the bridge join
dairy_pheno <- bridge %>%
  left_join(hm24, by = c("farm_name", "run")) %>%
  rename(SampleID = sample_id)

cat("Dairy phenotype rows:", nrow(dairy_pheno), "| cols:", ncol(dairy_pheno), "\n")

# 3. join: phenotypes first, then counts
dairy_counts_pheno <- dairy_pheno %>%
  inner_join(genus_counts, by = "SampleID")

cat("Result:", dim(dairy_counts_pheno), "\n")
cat("  phenotype cols:", ncol(dairy_pheno), "| genus cols:", ncol(genus_counts)-1, "\n")

# sanity checks
stopifnot(nrow(dairy_counts_pheno) == 92)
stopifnot(!any(grepl("\\.x$|\\.y$", names(dairy_counts_pheno))))  # no name collisions

write_csv(
  dairy_counts_pheno,
  "/home/dermot.kelly/Dermot_analysis/Phd/Paper_2/network_analysis/data/dairy_counts_with_phenotypes.csv"
)