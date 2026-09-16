####################################################################
# run_tully_merge.R
#
# VM-adapted rerun of tully_merge.R (from
# postprocessing/phenotype_merge_scripts/) against the new 17-flowcell
# combined genus table.
#
# tully_merge.R's original input, tully_data/genus_counts.csv, no longer
# exists on this VM. It was a Tully-only extract of the same genus_wide
# table 01_build_microbiome_matrices.ipynb produces, with ANI_ID/SeqID
# columns split out of SampleID (e.g. "Tully__10732_S4" -> ANI_ID
# "10732", SeqID "S4"). Reconstructed here as
# postprocessing/phenotype_merge_scripts/output/tully_genus_counts.csv
# by filtering the new genus_counts_no_controls_515F_806R.csv to
# SampleID starting with "Tully__" and parsing that same pattern (all
# 126 Tully rows matched it cleanly -- see git history for the
# derivation).
#
# Otherwise unchanged from tully_merge.R except: library(dplyr) moved
# to the top (the original called n_distinct()/count() at line 22
# before loading dplyr at line 67 -- only worked interactively because
# dplyr was already attached from a prior script in the same session),
# and the three I/O paths below.
####################################################################

library(dplyr)

data <- read.csv(
  "output/tully_genus_counts.csv",
  check.names = FALSE)

# what ID columns are there?
names(data)[1:6]
head(data[, 1:4])



beef_pheno <- read.csv(
  "/home/dermodkkelly/hpc_incoming/beef_meta/final with VFA pH.csv",
  check.names = TRUE)   # let R repair the empty/duplicate names

beef_pheno$short_tag <- as.character(beef_pheno$short_tag)
# standardise join keys
data$ANI_ID       <- as.character(data$ANI_ID)
beef_pheno$short_tag <- as.character(beef_pheno$short_tag)

cat("tully microbiome animals:", n_distinct(data$ANI_ID), "\n")
cat("phenotype animals:", n_distinct(beef_pheno$short_tag), "\n")
cat("overlap:", sum(unique(data$ANI_ID) %in% beef_pheno$short_tag), "\n")

# which tully animals have no phenotype match?
setdiff(unique(data$ANI_ID), beef_pheno$short_tag)


# duplicates in the phenotype file?
beef_pheno %>% count(short_tag) %>% filter(n > 1) %>% arrange(desc(n))


# ------------------------------------------------------------
# 1. Check repeated animal IDs in the Tully microbiome data
# ------------------------------------------------------------

data %>%
  count(ANI_ID) %>%
  filter(n > 1) %>%
  arrange(desc(n))

# inspect the actual duplicated samples
data %>%
  filter(ANI_ID %in% (
    data %>%
      count(ANI_ID) %>%
      filter(n > 1) %>%
      pull(ANI_ID)
  )) %>%
  select(SampleID, ANI_ID, SeqID)

# animals duplicated in microbiome
dup_micro <- data %>%
  count(ANI_ID) %>%
  filter(n > 1) %>%
  pull(ANI_ID)

# animals duplicated in phenotype data
dup_pheno <- beef_pheno %>%
  filter(!is.na(short_tag)) %>%
  count(short_tag) %>%
  filter(n > 1) %>%
  pull(short_tag)

# drop anything ambiguous in either dataset
drop_ids <- union(dup_micro, dup_pheno)

drop_ids
# keep only unambiguous microbiome and phenotype records
data_clean <- data %>%
  filter(!ANI_ID %in% drop_ids)

pheno_clean <- beef_pheno %>%
  filter(
    !is.na(short_tag),
    !short_tag %in% drop_ids
  )

# merge
beef_full <- data_clean %>%
  left_join(
    pheno_clean,
    by = c("ANI_ID" = "short_tag")
  )

# checks
cat("microbiome rows:", nrow(data_clean), "\n")
cat("unique microbiome animals:", n_distinct(data_clean$ANI_ID), "\n")
cat("beef_full rows:", nrow(beef_full), "\n")
cat("unique animals:", n_distinct(beef_full$ANI_ID), "\n")
cat("with CH4:", sum(!is.na(beef_full$CH4__g_d_)), "\n")

# should be empty
beef_full %>%
  count(ANI_ID) %>%
  filter(n > 1)

# check unmatched phenotypes
beef_full %>%
  filter(is.na(CH4__g_d_)) %>%
  select(SampleID, ANI_ID, SeqID)



write.csv(
  beef_full,
  "output/beef_full_clean.csv",
  row.names = FALSE
)
