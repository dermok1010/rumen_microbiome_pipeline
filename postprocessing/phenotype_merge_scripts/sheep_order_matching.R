####################################################################
# 01_map_sheep_samples_to_EID.R
#
# Maps sequenced sheep rumen microbiome samples to animal EIDs.
#
#   FASTQ sample name  --(Macrogen order sheet)-->  tube ID
#   tube ID            --(study sample sheet)   -->  EID + sampling date
#
# Covers four sources:
#   Sheep_CT24    -> "CT 2024" sheet          (key: Sample ID)
#   Sheep_NZAC    -> "INZAC 2023" sheet       (key: Sample ID, + VID where ambiguous)
#   Sheep_Haggot  -> "GF Hoggetts 24" sheet   (key: Rumen Sample No)
#   CTmicro       -> Macrogen Sample List Feb 2025 (key: VID)
#
# Sheep_CT25 is excluded: its study of origin is not yet identified.
#
# Output: sheep_all (one row per sequenced sample, with EID + samp_date)
#         unresolved_sheep_samples.csv (documented exclusions)
####################################################################

library(dplyr)
library(readr)
library(readxl)
library(stringr)
library(purrr)

setwd("/home/dermot.kelly/Dermot_analysis/Phd/Paper_2/merging/")


# ------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------
ORDER_SHEET   <- "sheep_data/sheep_order_sheet_2026.csv"
STUDY_SHEETS  <- "Methane Predict - Microbiome Samples for Grange.xlsx"
CTMICRO_SHEET <- "/data/Genetics/primary/R1240_microbiome/CT_microbiome_data/Macrogen Sample List Feb 2025.csv"
INVENTORY     <- paste0("/home/dermot.kelly/Dermot_analysis/Phd/Paper_2/",
                        "rumen_microbiome_pipeline/manifest_files/combined_run/",
                        "full_sample_inventory.tsv")

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------
clean_id <- function(x) gsub("[^A-Za-z0-9]", "", as.character(x))

parse_date_flex <- function(x) {
  if (inherits(x, "Date"))    return(x)
  if (inherits(x, "POSIXct")) return(as.Date(x))
  as.Date(as.character(x), format = "%d/%m/%Y")
}

# ==================================================================
# 1. Load inputs
# ==================================================================
sample_matching <- read_tsv(INVENTORY, show_col_types = FALSE)

order_sheet <- read_csv(ORDER_SHEET, show_col_types = FALSE)

sheets <- list(
  CT2024    = read_excel(STUDY_SHEETS, sheet = "CT 2024"),
  INZAC2023 = read_excel(STUDY_SHEETS, sheet = "INZAC 2023"),
  GF_Hogg   = read_excel(STUDY_SHEETS, sheet = "GF Hoggetts 24"),
  CT25      = read_excel("sheep_data/INZAC_April2025_RumenSampling.xlsx", sheet = "Sheet1")
)

ctmicro_sheet <- read.csv(CTMICRO_SHEET) %>%
  select(-any_of(c("X", "X.1", "X.2")))



# ==================================================================
# 2. Bridge: Macrogen sample name -> tube ID
# ==================================================================
# Two NZAC rows carry "sampleNo-VID" tube IDs (e.g. "35-28193") because the
# study sheet has two animals under that sample number; the VID disambiguates.
bridge <- order_sheet %>%
  select(macrogen_sample_id = 1, tube_id = 2) %>%
  mutate(
    macrogen_sample_id = str_trim(macrogen_sample_id),
    tube_id            = str_trim(tube_id),
    is_control = str_detect(macrogen_sample_id, "^\\(\\+\\)|^\\(-\\)"),
    cohort     = if_else(is_control, "CONTROL",
                         str_remove(macrogen_sample_id, "_\\d+$")),
    seq_number = suppressWarnings(as.integer(str_extract(macrogen_sample_id, "\\d+$"))),
    tube_num   = suppressWarnings(as.integer(tube_id)),
    tube_vid   = suppressWarnings(as.integer(str_extract(tube_id, "(?<=-)\\d+$"))),
    tube_num   = if_else(is.na(tube_num),
                         suppressWarnings(as.integer(str_extract(tube_id, "^\\d+(?=-)"))),
                         tube_num)
  )

# ==================================================================
# 3. Link sequencing samples to tube IDs
#    (Macrogen truncated "Sheep_NZAC2003" to "Sheep_NZAC" in the FASTQ names)
# ==================================================================
seq_side <- sample_matching %>%
  filter(str_detect(cohort, "^Sheep")) %>%
  mutate(number = as.integer(number))

bridge_j <- bridge %>%
  filter(!is_control) %>%
  mutate(cohort_h = str_replace(cohort, "^Sheep_NZAC2003$", "Sheep_NZAC")) %>%
  select(cohort_h, seq_number, tube_id, tube_num, tube_vid)

linked <- seq_side %>%
  left_join(bridge_j, by = c("cohort" = "cohort_h", "number" = "seq_number"))

stopifnot(sum(is.na(linked$tube_id)) == 0)   # every sequenced sample must be in the order sheet

# ==================================================================
# 4. Resolve tube ID -> EID + sampling date, per cohort
# ==================================================================

## 4a. Cohorts with unambiguous sample numbering ---------------------
build_cohort <- function(co, sheet_name, key_col) {
  df      <- sheets[[sheet_name]]
  datecol <- intersect(c("Date Sampled", "Date", "Sample Date"), names(df))[1]
  
  lookup <- df %>%
    transmute(tube_num  = suppressWarnings(as.integer(.data[[key_col]])),
              EID,
              samp_date = .data[[datecol]]) %>%
    filter(!is.na(tube_num))
  
  stopifnot(!any(duplicated(lookup$tube_num)))   # guard against many-to-many
  
  linked %>%
    filter(cohort == co) %>%
    left_join(lookup, by = "tube_num") %>%
    mutate(source_sheet = sheet_name,
           resolution   = if_else(is.na(EID), "UNRESOLVED_no_sheet_row", "by_tube_number"))
}

ct24_resolved   <- build_cohort("Sheep_CT24",   "CT2024",  "Sample ID")
haggot_resolved <- build_cohort("Sheep_Haggot", "GF_Hogg", "Rumen Sample No")
ct25_resolved   <- build_cohort("Sheep_CT25",   "CT25",    "Sample No.")

## 4b. NZAC: sample numbers 35 and 36 each map to two animals --------
inzac_lookup <- sheets$INZAC2023 %>%
  transmute(tube_num  = suppressWarnings(as.integer(`Sample ID`)),
            vid_num   = suppressWarnings(as.integer(VID)),
            EID,
            samp_date = `Date Sampled`) %>%
  filter(!is.na(tube_num))

ambiguous_nums <- inzac_lookup %>% count(tube_num) %>% filter(n > 1) %>% pull(tube_num)

nzac_resolved <- linked %>%
  filter(cohort == "Sheep_NZAC") %>%
  # VID takes precedence where the order sheet supplied one
  left_join(inzac_lookup %>% select(vid_num, EID, samp_date),
            by = c("tube_vid" = "vid_num")) %>%
  rename(EID_vid = EID, date_vid = samp_date) %>%
  left_join(inzac_lookup %>% filter(!tube_num %in% ambiguous_nums) %>%
              select(tube_num, EID, samp_date),
            by = "tube_num") %>%
  mutate(
    EID          = coalesce(EID_vid, EID),
    samp_date    = coalesce(date_vid, samp_date),
    source_sheet = "INZAC2023",
    resolution   = case_when(
      !is.na(EID_vid)              ~ "by_VID",
      !is.na(EID)                  ~ "by_tube_number",
      tube_num %in% ambiguous_nums ~ "UNRESOLVED_two_animals_share_sample_number",
      TRUE                         ~ "UNRESOLVED_no_sheet_row"
    )
  ) %>%
  select(-EID_vid, -date_vid)

## 4c. CTmicro: FASTQ stem is a zero-padded VID ----------------------
ctmicro_resolved <- sample_matching %>%
  filter(batch == "CTmicro") %>%
  mutate(number   = as.character(number),
         cohort   = "CTmicro",
         tube_num = suppressWarnings(as.integer(str_remove(sample_id, "^CTmicro__")))) %>%
  left_join(
    ctmicro_sheet %>%
      transmute(tube_num = suppressWarnings(as.integer(VID)),
                EID, samp_date = Date, Group),
    by = "tube_num"
  ) %>%
  mutate(tube_id      = as.character(tube_num),
         source_sheet = "MacrogenSampleListFeb2025",
         resolution   = if_else(is.na(EID), "UNRESOLVED_no_sheet_row", "by_VID"))

# ==================================================================
# 5. Combine and standardise
# ==================================================================
sheep_all <- bind_rows(
  ct24_resolved    %>% mutate(number = as.character(number), samp_date = parse_date_flex(samp_date)),
  nzac_resolved    %>% mutate(number = as.character(number), samp_date = parse_date_flex(samp_date)),
  haggot_resolved  %>% mutate(number = as.character(number), samp_date = parse_date_flex(samp_date)),
  ct25_resolved    %>% mutate(number = as.character(number), samp_date = parse_date_flex(samp_date)),
  ctmicro_resolved %>% mutate(samp_date = parse_date_flex(samp_date))
) %>% 
  mutate(EID = na_if(clean_id(EID), "NA"))

# ==================================================================
# 6. Summary and documented exclusions
# ==================================================================
sheep_all %>%
  group_by(cohort) %>%
  summarise(n = n(),
            resolved  = sum(!is.na(EID)),
            with_date = sum(!is.na(samp_date)),
            .groups = "drop") %>%
  print()

cat("\nTOTAL:", nrow(sheep_all),
    "| resolved:", sum(!is.na(sheep_all$EID)),
    "| unresolved:", sum(is.na(sheep_all$EID)), "\n\n")

print(count(sheep_all, resolution))


#write_csv(unresolved, "sheep_data/unresolved_sheep_samples.csv")
write_csv(sheep_all,  "sheep_data/sheep_sample_to_EID.csv")
