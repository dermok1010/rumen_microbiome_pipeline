

data <- read.csv(
  "/home/dermot.kelly/Dermot_analysis/Phd/Paper_2/microbiome_ml/tully_data/genus_counts.csv",
  check.names = FALSE)

# what ID columns are there?
names(data)[1:6]
head(data[, 1:4])



beef_pheno <- read.csv(
  "/home/dermot.kelly/Dermot_analysis/Phd/Paper_2/merging/beef_meta/final with VFA pH.csv",
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

intersect(c("20906","50105","60106"), unique(data$ANI_ID))
beef_pheno %>%
  filter(short_tag %in% c("20906","50105","60106")) %>%
  select(short_tag, Year, Intake_Id, CH4__g_d_, Total_DMI) %>%
  arrange(short_tag, Year)


# for the non-duplicated tully animals, which year is their phenotype?
tully_anis <- unique(data$ANI_ID)

beef_pheno %>%
  filter(short_tag %in% tully_anis) %>%
  count(short_tag) %>%
  filter(n == 1) %>%                       # the unambiguous ones
  left_join(beef_pheno %>% select(short_tag, Year), by = "short_tag") %>%
  count(Year)


# --- merge tully beef microbiome with phenotypes, dropping year-ambiguous animals ---
ambiguous <- c("50105", "60106")

beef_full <- data %>%
  filter(!ANI_ID %in% ambiguous) %>%
  left_join(beef_pheno, by = c("ANI_ID" = "short_tag"))

# checks
cat("beef microbiome rows (after drop):", nrow(data) - sum(unique(data$ANI_ID) %in% ambiguous), "\n")
cat("beef_full rows:", nrow(beef_full), "\n")   # should be 122, no inflation
cat("with CH4_g_d:", sum(!is.na(beef_full$CH4__g_d_)), "\n")

# confirm no duplication crept in
beef_full %>% count(ANI_ID) %>% filter(n > 1)   # should be empty


library(dplyr)

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
  "/home/dermot.kelly/Dermot_analysis/Phd/Paper_2/microbiome_ml/tully_data/beef_full_clean.csv",
  row.names = FALSE
)