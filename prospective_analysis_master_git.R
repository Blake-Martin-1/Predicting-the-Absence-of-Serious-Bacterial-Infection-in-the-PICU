### Script to evaluate the prospectively collected inputs and predictions of the RF and PLR models. Data is from 9/1/2023 to 4/06/24

# Call setup script to import needed filepaths
source("/home/martinbl/sbi_blake/setup_10yr.R")

# Load additional scripts
library(lubridate)
library(data.table)
library(pedbp)
library(DescTools)
library(stringi)
library(stringr)
library(readxl)
library(devtools)
library(rriskDistributions)
library(mltools)
library(tidyr)
library(data.table)
library(rsample)
library(parsnip)
library(ranger)
library(pROC)
library(caret)
library(iml)

# Helpful function to view vector as a df
# This function turns a vector into a vertical one column dataframe for ease of visualization
get_df <- function(the_vec){
  View(data.frame("var" = the_vec))
}

# Function to return a dataframe in which you only see the columns with names that have a given string within them
get_cols <- function(df, temp_string){
  temp_loc <- str_detect(string = colnames(df), pattern = temp_string)
  return(df[, temp_loc])
}

# Function to return a dataframe with the names of all of the columns
see_col_names <- function(df){
  temp_names <- colnames(df)
  new_df <- data.frame("col_names" = temp_names)
  return(new_df)
}

# Function to return a data frame that has all of the classes of the columns within a dataframe
show_classes <- function(df){
  class_vec <- sapply(X = df, FUN = class)
  temp_df <- data.frame("col_name" = class_vec)
  return(temp_df)
}

# Function to change a given POSIXct object to the America/Denver time zone, e.g. UTC value would lose 6 hours
to_mdt <- function(df, my_col){
  df[[paste(my_col)]] <- with_tz(df[[paste(my_col)]], tzone = "America/Denver")
  return(df)
}

# change to POSIXct format
to_posix <- function(df, my_col){
  df[[paste(my_col)]] <- ymd_hms(as.character(df[[paste(my_col)]]), tz = "America/Denver")
  return(df)
}



# Load the prospective data:
epic_7_23  <- read_csv(
  file = "/phi/sbi/prospective_data/Prospective/model_output/model_output_Sep2023.csv",
  col_types = cols(
    ICU_START_INSTANT = col_character(),
    VALID_START_INSTANT = col_character(),
    VALID_END_INSTANT = col_character(),
    INSTANT_UTC_DTTM = col_character(),
    ORIGIN = col_character(),
    IS_INTERVENTIONAL_RADIOLOGY = col_character(),
    DBP_TS_IS_COMPLETE = col_skip(),
    DBP_TS_LEN = col_skip(),
    FIO2_TS_IS_COMPLETE = col_skip(),
    FIO2_TS_LEN = col_skip(),
    HR_TS_IS_COMPLETE = col_skip(),
    HR_TS_LEN = col_skip(),
    O2SAT_TS_IS_COMPLETE = col_skip(),
    O2SAT_TS_LEN = col_skip(),
    RR_TS_IS_COMPLETE = col_skip(),
    RR_TS_LEN = col_skip(),
    SBP_TS_IS_COMPLETE = col_skip(),
    SBP_TS_LEN = col_skip()
  )
)

epic_9_23 <- read_csv(
  file = "/phi/sbi/prospective_data/Prospective/model_output/model_output_Oct2023.csv",
  col_types = cols(
    ICU_START_INSTANT = col_character(),
    VALID_START_INSTANT = col_character(),
    VALID_END_INSTANT = col_character(),
    INSTANT_UTC_DTTM = col_character(),
    ORIGIN = col_character(),
    IS_INTERVENTIONAL_RADIOLOGY = col_character(),
    DBP_TS_IS_COMPLETE = col_skip(),
    DBP_TS_LEN = col_skip(),
    FIO2_TS_IS_COMPLETE = col_skip(),
    FIO2_TS_LEN = col_skip(),
    HR_TS_IS_COMPLETE = col_skip(),
    HR_TS_LEN = col_skip(),
    O2SAT_TS_IS_COMPLETE = col_skip(),
    O2SAT_TS_LEN = col_skip(),
    RR_TS_IS_COMPLETE = col_skip(),
    RR_TS_LEN = col_skip(),
    SBP_TS_IS_COMPLETE = col_skip(),
    SBP_TS_LEN = col_skip()
  )
)

epic_10_23 <- read_csv(
  file = "/phi/sbi/prospective_data/Prospective/model_output/model_output_Nov2023.csv",
  col_types = cols(
    ICU_START_INSTANT = col_character(),
    VALID_START_INSTANT = col_character(),
    VALID_END_INSTANT = col_character(),
    INSTANT_UTC_DTTM = col_character(),
    ORIGIN = col_character(),
    IS_INTERVENTIONAL_RADIOLOGY = col_character(),
    DBP_TS_IS_COMPLETE = col_skip(),
    DBP_TS_LEN = col_skip(),
    FIO2_TS_IS_COMPLETE = col_skip(),
    FIO2_TS_LEN = col_skip(),
    HR_TS_IS_COMPLETE = col_skip(),
    HR_TS_LEN = col_skip(),
    O2SAT_TS_IS_COMPLETE = col_skip(),
    O2SAT_TS_LEN = col_skip(),
    RR_TS_IS_COMPLETE = col_skip(),
    RR_TS_LEN = col_skip(),
    SBP_TS_IS_COMPLETE = col_skip(),
    SBP_TS_LEN = col_skip()
  )
)

epic_11_23 <- read_csv(
  file = "/phi/sbi/prospective_data/Prospective/model_output/model_output_Dec2023.csv",
  col_types = cols(
    ICU_START_INSTANT = col_character(),
    VALID_START_INSTANT = col_character(),
    VALID_END_INSTANT = col_character(),
    INSTANT_UTC_DTTM = col_character(),
    ORIGIN = col_character(),
    IS_INTERVENTIONAL_RADIOLOGY = col_character(),
    DBP_TS_IS_COMPLETE = col_skip(),
    DBP_TS_LEN = col_skip(),
    FIO2_TS_IS_COMPLETE = col_skip(),
    FIO2_TS_LEN = col_skip(),
    HR_TS_IS_COMPLETE = col_skip(),
    HR_TS_LEN = col_skip(),
    O2SAT_TS_IS_COMPLETE = col_skip(),
    O2SAT_TS_LEN = col_skip(),
    RR_TS_IS_COMPLETE = col_skip(),
    RR_TS_LEN = col_skip(),
    SBP_TS_IS_COMPLETE = col_skip(),
    SBP_TS_LEN = col_skip()
  )
)

epic_12_23  <- read_csv(
  file = "/phi/sbi/prospective_data/Prospective/model_output/model_output_Jan2024.csv",
  col_types = cols(
    ICU_START_INSTANT = col_character(),
    VALID_START_INSTANT = col_character(),
    VALID_END_INSTANT = col_character(),
    INSTANT_UTC_DTTM = col_character(),
    ORIGIN = col_character(),
    IS_INTERVENTIONAL_RADIOLOGY = col_character(),
    DBP_TS_IS_COMPLETE = col_skip(),
    DBP_TS_LEN = col_skip(),
    FIO2_TS_IS_COMPLETE = col_skip(),
    FIO2_TS_LEN = col_skip(),
    HR_TS_IS_COMPLETE = col_skip(),
    HR_TS_LEN = col_skip(),
    O2SAT_TS_IS_COMPLETE = col_skip(),
    O2SAT_TS_LEN = col_skip(),
    RR_TS_IS_COMPLETE = col_skip(),
    RR_TS_LEN = col_skip(),
    SBP_TS_IS_COMPLETE = col_skip(),
    SBP_TS_LEN = col_skip()
  )
)

epic_1_24  <- read_csv(
  file = "/phi/sbi/prospective_data/Prospective/model_output/model_output_Feb2024.csv",
  col_types = cols(
    ICU_START_INSTANT = col_character(),
    VALID_START_INSTANT = col_character(),
    VALID_END_INSTANT = col_character(),
    INSTANT_UTC_DTTM = col_character(),
    ORIGIN = col_character(),
    IS_INTERVENTIONAL_RADIOLOGY = col_character(),
    DBP_TS_IS_COMPLETE = col_skip(),
    DBP_TS_LEN = col_skip(),
    FIO2_TS_IS_COMPLETE = col_skip(),
    FIO2_TS_LEN = col_skip(),
    HR_TS_IS_COMPLETE = col_skip(),
    HR_TS_LEN = col_skip(),
    O2SAT_TS_IS_COMPLETE = col_skip(),
    O2SAT_TS_LEN = col_skip(),
    RR_TS_IS_COMPLETE = col_skip(),
    RR_TS_LEN = col_skip(),
    SBP_TS_IS_COMPLETE = col_skip(),
    SBP_TS_LEN = col_skip()
  )
)

epic_2_24  <- read_csv(
  file = "/phi/sbi/prospective_data/Prospective/model_output/model_output_actual_feb2024.csv",
  col_types = cols(
    ICU_START_INSTANT = col_character(),
    VALID_START_INSTANT = col_character(),
    VALID_END_INSTANT = col_character(),
    INSTANT_UTC_DTTM = col_character(),
    ORIGIN = col_character(),
    IS_INTERVENTIONAL_RADIOLOGY = col_character(),
    DBP_TS_IS_COMPLETE = col_skip(),
    DBP_TS_LEN = col_skip(),
    FIO2_TS_IS_COMPLETE = col_skip(),
    FIO2_TS_LEN = col_skip(),
    HR_TS_IS_COMPLETE = col_skip(),
    HR_TS_LEN = col_skip(),
    O2SAT_TS_IS_COMPLETE = col_skip(),
    O2SAT_TS_LEN = col_skip(),
    RR_TS_IS_COMPLETE = col_skip(),
    RR_TS_LEN = col_skip(),
    SBP_TS_IS_COMPLETE = col_skip(),
    SBP_TS_LEN = col_skip()
  )
)

epic_3_24  <- read_csv(
  file = "/phi/sbi/prospective_data/Prospective/model_output/model_output_March2024.csv",
  col_types = cols(
    ICU_START_INSTANT = col_character(),
    VALID_START_INSTANT = col_character(),
    VALID_END_INSTANT = col_character(),
    INSTANT_UTC_DTTM = col_character(),
    ORIGIN = col_character(),
    IS_INTERVENTIONAL_RADIOLOGY = col_character(),
    DBP_TS_IS_COMPLETE = col_skip(),
    DBP_TS_LEN = col_skip(),
    FIO2_TS_IS_COMPLETE = col_skip(),
    FIO2_TS_LEN = col_skip(),
    HR_TS_IS_COMPLETE = col_skip(),
    HR_TS_LEN = col_skip(),
    O2SAT_TS_IS_COMPLETE = col_skip(),
    O2SAT_TS_LEN = col_skip(),
    RR_TS_IS_COMPLETE = col_skip(),
    RR_TS_LEN = col_skip(),
    SBP_TS_IS_COMPLETE = col_skip(),
    SBP_TS_LEN = col_skip()
  )
)

epic_4_24  <- read_csv(
  file = "/phi/sbi/prospective_data/Prospective/model_output/model_output_April2024.csv",
  col_types = cols(
    ICU_START_INSTANT = col_character(),
    VALID_START_INSTANT = col_character(),
    VALID_END_INSTANT = col_character(),
    INSTANT_UTC_DTTM = col_character(),
    ORIGIN = col_character(),
    IS_INTERVENTIONAL_RADIOLOGY = col_character(),
    DBP_TS_IS_COMPLETE = col_skip(),
    DBP_TS_LEN = col_skip(),
    FIO2_TS_IS_COMPLETE = col_skip(),
    FIO2_TS_LEN = col_skip(),
    HR_TS_IS_COMPLETE = col_skip(),
    HR_TS_LEN = col_skip(),
    O2SAT_TS_IS_COMPLETE = col_skip(),
    O2SAT_TS_LEN = col_skip(),
    RR_TS_IS_COMPLETE = col_skip(),
    RR_TS_LEN = col_skip(),
    SBP_TS_IS_COMPLETE = col_skip(),
    SBP_TS_LEN = col_skip()
  )
)

epic_5_24  <- read_csv(
  file = "/phi/sbi/prospective_data/Prospective/model_output/model_output_May2024.csv",
  col_types = cols(
    ICU_START_INSTANT = col_character(),
    VALID_START_INSTANT = col_character(),
    VALID_END_INSTANT = col_character(),
    INSTANT_UTC_DTTM = col_character(),
    ORIGIN = col_character(),
    IS_INTERVENTIONAL_RADIOLOGY = col_character(),
    DBP_TS_IS_COMPLETE = col_skip(),
    DBP_TS_LEN = col_skip(),
    FIO2_TS_IS_COMPLETE = col_skip(),
    FIO2_TS_LEN = col_skip(),
    HR_TS_IS_COMPLETE = col_skip(),
    HR_TS_LEN = col_skip(),
    O2SAT_TS_IS_COMPLETE = col_skip(),
    O2SAT_TS_LEN = col_skip(),
    RR_TS_IS_COMPLETE = col_skip(),
    RR_TS_LEN = col_skip(),
    SBP_TS_IS_COMPLETE = col_skip(),
    SBP_TS_LEN = col_skip()
  )
)

epic_6_24  <- read_csv(
  file = "/phi/sbi/prospective_data/Prospective/model_output/model_output_June2024.csv",
  col_types = cols(
    ICU_START_INSTANT = col_character(),
    VALID_START_INSTANT = col_character(),
    VALID_END_INSTANT = col_character(),
    INSTANT_UTC_DTTM = col_character(),
    ORIGIN = col_character(),
    IS_INTERVENTIONAL_RADIOLOGY = col_character(),
    DBP_TS_IS_COMPLETE = col_skip(),
    DBP_TS_LEN = col_skip(),
    FIO2_TS_IS_COMPLETE = col_skip(),
    FIO2_TS_LEN = col_skip(),
    HR_TS_IS_COMPLETE = col_skip(),
    HR_TS_LEN = col_skip(),
    O2SAT_TS_IS_COMPLETE = col_skip(),
    O2SAT_TS_LEN = col_skip(),
    RR_TS_IS_COMPLETE = col_skip(),
    RR_TS_LEN = col_skip(),
    SBP_TS_IS_COMPLETE = col_skip(),
    SBP_TS_LEN = col_skip()
  )
)

epic_7_24  <- read_csv(
  file = "/phi/sbi/prospective_data/Prospective/model_output/model_output_July2024.csv",
  col_types = cols(
    ICU_START_INSTANT = col_character(),
    VALID_START_INSTANT = col_character(),
    VALID_END_INSTANT = col_character(),
    INSTANT_UTC_DTTM = col_character(),
    ORIGIN = col_character(),
    IS_INTERVENTIONAL_RADIOLOGY = col_character(),
    DBP_TS_IS_COMPLETE = col_skip(),
    DBP_TS_LEN = col_skip(),
    FIO2_TS_IS_COMPLETE = col_skip(),
    FIO2_TS_LEN = col_skip(),
    HR_TS_IS_COMPLETE = col_skip(),
    HR_TS_LEN = col_skip(),
    O2SAT_TS_IS_COMPLETE = col_skip(),
    O2SAT_TS_LEN = col_skip(),
    RR_TS_IS_COMPLETE = col_skip(),
    RR_TS_LEN = col_skip(),
    SBP_TS_IS_COMPLETE = col_skip(),
    SBP_TS_LEN = col_skip()
  )
)

epic_8_24  <- read_csv(
  file = "/phi/sbi/prospective_data/Prospective/model_output/model_output_August2024.csv",
  col_types = cols(
    ICU_START_INSTANT = col_character(),
    VALID_START_INSTANT = col_character(),
    VALID_END_INSTANT = col_character(),
    INSTANT_UTC_DTTM = col_character(),
    ORIGIN = col_character(),
    IS_INTERVENTIONAL_RADIOLOGY = col_character(),
    DBP_TS_IS_COMPLETE = col_skip(),
    DBP_TS_LEN = col_skip(),
    FIO2_TS_IS_COMPLETE = col_skip(),
    FIO2_TS_LEN = col_skip(),
    HR_TS_IS_COMPLETE = col_skip(),
    HR_TS_LEN = col_skip(),
    O2SAT_TS_IS_COMPLETE = col_skip(),
    O2SAT_TS_LEN = col_skip(),
    RR_TS_IS_COMPLETE = col_skip(),
    RR_TS_LEN = col_skip(),
    SBP_TS_IS_COMPLETE = col_skip(),
    SBP_TS_LEN = col_skip()
  )
)

# Fix a few issues with datasets to ensure can bind the rows (i.e. data class)

epic_pros_raw <- bind_rows(epic_9_23, epic_10_23, epic_11_23, epic_12_23, epic_1_24, epic_2_24, epic_3_24, epic_4_24, epic_5_24, epic_6_24,
                           epic_7_24, epic_8_24)
epic_pros_data <- epic_pros_raw[, !grepl("_TS$", names(epic_pros_raw))] # Remove columns with long, incomprehensible strings

# # Call separate script to extract all of the individual lab and vital sign values, I think optional?
# source(file = "/phi/sbi/sbi_blake/extract_values.R")


# Rename certain columns
pros_df <- epic_pros_data %>% rename(score_time = INSTANT_UTC_DTTM, model_score = TOTAL_SCORE, model_type = ACUITY_SYSTEM_ID)

pros_df$model_type[pros_df$model_type == "100133"] <- "LR_no_abx"
pros_df$model_type[pros_df$model_type == "100134"] <- "LR_yes_abx"
pros_df$model_type[pros_df$model_type == "100148"] <- "RF_no_abx"
pros_df$model_type[pros_df$model_type == "100149"] <- "RF_yes_abx"

# Lower case of column names
colnames(pros_df) <- str_to_lower(colnames(pros_df))

# Fix relevant date/times
pros_df <- pros_df %>% mutate(score_time = as.POSIXct(score_time, tz = "UTC", format = "%Y-%m-%d %H:%M:%S"))
pros_df$score_time <- with_tz(time = pros_df$score_time, tzone = "America/Denver")

# Rename ICU admission time and fix the timezone
pros_df <- pros_df %>% mutate(picu_adm_date_time = as.POSIXct(icu_start_instant, format = "%Y-%m-%d %H:%M:%S"))
pros_df <- pros_df %>% mutate(picu_adm_date_time = force_tz(picu_adm_date_time, tzone = "America/Denver"))


# Filter for correct study start time
study_start_time <- as.POSIXct("2023-09-01 00:00:01", tz = "America/Denver")
pros_df <- pros_df %>% filter(score_time >= study_start_time)

# Load in the PICU admit/transfer dates and times to identify the correct identifier for each icu admission for each child
adt_raw <- read_csv(file = "/phi/sbi/prospective_data/Prospective/adt_export_pros_100724.csv")
adt_df <- adt_raw
adt_df$intime <- as.POSIXct(adt_df$intime, format = "%m/%d/%y %H:%M")
adt_df$outtime <- as.POSIXct(adt_df$outtime, format = "%m/%d/%y %H:%M")

adt_df <- adt_df %>% mutate(study_id = paste0(pat_enc_csn_id, "_", cicu_event_count))

# Convert to data.tables
pros_dt <- as.data.table(pros_df)
adt_dt <- as.data.table(adt_df)

# Create interval columns for score_time (required by foverlaps)
pros_dt[, score_start := score_time]
pros_dt[, score_end := score_time]

# Ensure all datetime columns have the same timezone
tz_target <- "America/Denver"

pros_dt[, score_time := as.POSIXct(score_time, tz = tz_target)]
pros_dt[, score_start := as.POSIXct(score_start, tz = tz_target)]
pros_dt[, score_end := as.POSIXct(score_end, tz = tz_target)]

adt_dt[, intime := as.POSIXct(intime, tz = tz_target)]
adt_dt[, outtime := as.POSIXct(outtime, tz = tz_target)]

# Set keys for adt_dt
setkey(adt_dt, pat_enc_csn_id, intime, outtime)

# Perform non-equi join
pros_dt_with_study_id <- foverlaps(
  x = pros_dt,  # now has interval columns
  y = adt_dt[, .(pat_enc_csn_id, intime, outtime, study_id)],
  by.x = c("pat_enc_csn_id", "score_start", "score_end"),
  by.y = c("pat_enc_csn_id", "intime", "outtime"),
  type = "within",
  nomatch = NA
)

# Optional: drop the interval columns if you don't want to keep them
pros_dt_with_study_id[, c("score_start", "score_end") := NULL]

pros_dt_with_study_id <- pros_dt_with_study_id %>%
  mutate(hours_since_picu_adm = as.numeric(difftime(score_time, picu_adm_date_time, units = "hours")))

# Filter out predictions for patients without in/out times (wrong age, admitted before study start date, etc)
pros_dt_with_study_id <- pros_dt_with_study_id %>% filter(!is.na(intime))
pros_df <- pros_dt_with_study_id


# Reorder
pros_df <- pros_df %>% relocate(study_id, pat_enc_csn_id, picu_adm_date_time, score_time, hours_since_picu_adm) %>% arrange(study_id, picu_adm_date_time, score_time)


########### Determine Presence of SBI ###########
# Determination of culture negative sepsis with lactate > 2 and SBP min < 5th%ile
pros_cx_neg <- pros_df %>% mutate(cx_neg_sepsis = ifelse(sbp_min < 5 & lactate_max > 2, yes = 1, no = 0)) # note this gives anyone with the lactate sbp criteria, but still need BCx sent


# Determination of pneumonia presence
# Read in the csv with csn's of patients with VPS dx of pna
vps_pna_df_sep <- read_xlsx(path = "/phi/sbi/prospective_data/Prospective/model_output/9_1_23_to_12_31_23_PNA_rerun.xlsx")
vps_pna_jan_mar <- read_xlsx(path = "/phi/sbi/prospective_data/Prospective/model_output/jan_2024_to_mar_2024_pna.xlsx")
vps_pna_apr_jun <- read_xlsx(path = "/phi/sbi/prospective_data/Prospective/model_output/apr_2024_to_june_2024_pna.xlsx")
vps_pna_jul_sep <- read_xlsx(path = "/phi/sbi/prospective_data/Prospective/model_output/7_1_24_to_9_30_24_PNA.xlsx")

# Fix names and classes of columns, and then bind all vps df's together
sept <- vps_pna_df_sep %>% dplyr::select(MRN, ICUAdmDateTime, AccountNum, `Present On Admission`)
colnames(sept) <- c("mrn", "picu_adm_date_time", "hsp_account_id", "pna_on_admit")

jan <- vps_pna_jan_mar %>% dplyr::select(MRN, ICUAdmDateTime, AccountNum, `Present On Admission`)
colnames(jan) <- c("mrn", "picu_adm_date_time", "hsp_account_id", "pna_on_admit")

apr <- vps_pna_apr_jun %>% dplyr::select(MRN, ICUAdmDateTime, AccountNum, `Present On Admission`)
colnames(apr) <- c("mrn", "picu_adm_date_time", "hsp_account_id", "pna_on_admit")

july <- vps_pna_jul_sep %>% dplyr::select(MRN, ICUAdmDateTime, AccountNum, `Present On Admission`)
colnames(july) <- c("mrn", "picu_adm_date_time", "hsp_account_id", "pna_on_admit")

pna_all <- bind_rows(sept, jan, apr, july)

# Fix all of the pna timezone info
pna_all <- pna_all %>% mutate(picu_adm_date_time = force_tz(picu_adm_date_time, tzone = "America/Denver"))


# Load in demographic info to help match csn and har
demo_ref <- read_csv(file = "/phi/sbi/prospective_data/Prospective/model_output/demog_export_pros_092024.csv")
id_df <- demo_ref %>% dplyr::select(hsp_account_id, pat_enc_csn_id, pat_mrn_id) %>% distinct()

# Make har and mrn numeric to match pros dfs
pna_all$hsp_account_id <- as.numeric(pna_all$hsp_account_id)
pna_all$mrn <- as.numeric(pna_all$mrn)

# Filter for pneumonias present on admission
pna_admit <- pna_all %>% filter(pna_on_admit == "Yes") %>% distinct() # from 1,004 pna rows to 625 pneumonias

# Ensure demo_ref matches the name demo_df used in below code
demo_df <- demo_ref



# ###############################################################################
# Want to identify children with pneumonia present on admission for that PICU stay, so bind by csn and picu_admit date to assign pneumonias
csn_har <- id_df %>% dplyr::select(pat_enc_csn_id, hsp_account_id) %>% distinct()

# Add csn to the pna df
pna_adm_csn <- pna_admit %>% left_join(csn_har, by = "hsp_account_id")

# Those with na for the csn will be too old or too young so filter out
pna_adm_csn <- pna_adm_csn %>% dplyr::filter(!is.na(pat_enc_csn_id)) # down to 565 pneumonias

# Now need to make sure that the picu admit dates match up so that I assign the pneumonias to the right rows of pros
## ── 1.  Make interval columns for pros_cx_neg  ──────────────────────────
setDT(pros_cx_neg)
pros_cx_neg[ , `:=`(
  window_start = picu_adm_date_time - lubridate::dhours(1),
  window_end   = picu_adm_date_time + lubridate::dhours(1)
)]

## ── 2.  Prepare pna_adm_csn as zero-length intervals  ───────────────────
setDT(pna_adm_csn)
pna_adm_csn[ , `:=`(pna_start = picu_adm_date_time,
                    pna_end   = picu_adm_date_time)]

## ── 3.  Key the tables for an interval overlap join  ────────────────────
setkey(pros_cx_neg, pat_enc_csn_id, window_start, window_end)
setkey(pna_adm_csn, pat_enc_csn_id, pna_start,  pna_end)

## ── 4.  Initialise flag to 0  ───────────────────────────────────────────
pros_cx_neg[ , pna_1_0 := 0L]

## ── 5.  Mark rows that have *any* overlap (±1 h) with a pneumonia row  ──
pros_cx_neg[
  pna_adm_csn,                                   # i-table
  on      = .(pat_enc_csn_id,
              window_start <= pna_start,         # interval overlap logic
              window_end   >= pna_end),
  pna_1_0 := 1L,
  nomatch = 0L                                   # only rows that overlap
]

## ── 6.  Clean up helper columns if you like  ────────────────────────────
pros_cx_neg[ , c("window_start", "window_end") := NULL]





# Assign the pneumonia label to these patients
pros_pna <- pros_cx_neg #62212 total predictions with pna outcome (2.6% of predictions in patients with pna), 448 study IDs with pna

# Fix cx neg sepsis in setting of pna
pros_pna$cx_neg_sepsis[pros_pna$pna_1_0 == 1] <- 0

# Run code to identify bacterial infections from microbiologic specimens, note that micro_df is the ultimate df produced by below script call
source(file = "/home/martinbl/sbi_blake/prospective_validation_files/micro_tidying_ongoing_pros_validation.R")

###### Note that the micro file (and prob other files) don't have outcome data for kids admitted to the regular hospital prior to 9/1/2023.
# micro_df is all positive encounters by study id with 598 positive specimens, 157 study id values with positive cultures out of 2625 picu encounters = 6%


## Identify study id's with a micro sbi
pros_all <- pros_pna %>% mutate(micro_sbi_1_0 = ifelse(study_id %in% micro_df$study_id, yes = 1, no = 0))

# Clarify that cnss has to occur in setting of drawing a blood culture, drops cnss rate from 7% to 3.7%
pros_all <- pros_all %>% mutate(cnss_true = ifelse(cx_neg_sepsis == 1 & study_id %in% bcx_mrn_hosp_adm$study_id, yes = 1, no = 0))

# Now identify sbi overall with combo of pneumonia, cnss, and micro sbi's
pros_all <- pros_all %>% mutate(sbi_present = ifelse((pna_1_0 + cnss_true + micro_sbi_1_0) > 0, yes = 1, no = 0)) # 26.9% of all predictions have an SBI present in 24hr period around PICU admit

# Identify study_id level rate of SBI
# BY prediction
sbi_rate_by_prediction <- round(sum(pros_all$sbi_present) / nrow(pros_all) * 100, digits = 2)

# By study_id
n_with_sbi <- pros_all %>%
  group_by(study_id) %>%                 # one group per PICU stay
  summarise(has_sbi = any(sbi_present == 1), .groups = "drop") %>%
  filter(has_sbi) %>%                    # keep stays with ≥1 “1”
  nrow()

sbi_rate_by_study_id <- round(n_with_sbi / n_distinct(pros_all$study_id) * 100, digits = 2) # 22.6% of all PICU admissions



# Now filter to appropriate time period based on available data
study_end <- as.POSIXct("2024-08-31 23:59:59", tz = "America/Denver")
study_start_time <- as.POSIXct("2023-09-01 00:00:01", tz = "America/Denver")
pros_full_data <- pros_all # store full dataset in case needed
pros_all <- pros_all %>% filter(score_time <= study_end) %>% filter (score_time >= study_start_time)# filter for end of study period

# Identify CDiff patients
cdiff_test_info <- micro %>% filter(organism_name == "c_diff")

# Get list of distinct patient ids, csn's, and PICU admit date/times so that we can filter for CDiff tests after PICU admission
id_csn_adm <- pros_all %>% dplyr::select(study_id, pat_enc_csn_id, picu_adm_date_time) %>% distinct()


###### Now call the script to load in the antimicrobials and determine antibiotic administration prior to each prediction
source(file = "/phi/sbi/sbi_blake/abx_pros_validation.R")


# Filter out predictions that occur >24 hours after picu admission
pros_all<- pros_all %>% mutate(t_diff = difftime(score_time, picu_adm_date_time, units = "hours"))
pros_all <- pros_all %>% filter(t_diff <= 24) # reduces predictions to 249586



# Recalculate number of encounters
n_enc <- n_distinct(pros_all$study_id) #2794 study encounters
n_pts <- nrow(pros_all %>% dplyr::select(pat_enc_csn_id) %>% distinct()) #2,628 unique hospitalizations
n_pred_pre <- nrow(pros_all) #249,586 predictions


# Quality control to see the number of unique encounters and those with and without SBI
n_enc <- n_distinct(pros_all$study_id) #
n_enc # show number of unique study encounters

# Fix it so that sbi_present is 1 if cx_neg_sepsis is 1 at any point during the PICU encounter
# Step 1: Identify study_ids where cx_neg_sepsis is 1 in any row
study_ids_with_cx_neg_sepsis <- pros_all %>%
  group_by(study_id) %>%
  summarise(max_cx_neg_sepsis = max(cnss_true)) %>%
  filter(max_cx_neg_sepsis == 1) %>%
  pull(study_id)

# Create another column to show if cx_neg_sepsis was every present (instead of just becoming 1 once)
pros_all <- pros_all %>%
  mutate("ever_cx_neg_sepsis" = if_else(study_id %in% study_ids_with_cx_neg_sepsis, 1, 0))

pros_all <- pros_all %>% mutate(sbi_present = ifelse(sbi_present == 1 | ever_cx_neg_sepsis == 1, yes = 1, no = sbi_present))

# Add in race, ethnicity, language, insurance data from demo_slim after adding hsp_account_id via id_df,
pros_all <- pros_all %>% left_join(id_df, by = "pat_enc_csn_id")

# Remove those without hsp_account ids because they represent patients admitted before study period and others who didn't move through PICU normally (e.g. to CICU, wrong age)
pros_all_rows_no_har <- pros_all
pros_all <- pros_all %>% filter(!is.na(hsp_account_id)) #249576 rows

pros_all <- pros_all %>% rename(old_race = race) %>% left_join(demo_slim %>% dplyr::select(hsp_account_id, race, ethnicity, language, insurance_type) %>% distinct(), by = "hsp_account_id")

# Fix ethnicity and other SDOH categories
pros_all$ethnicity[pros_all$ethnicity == "Decline to Answer"] <- "Other_Or_Unknown"
pros_all$ethnicity[pros_all$ethnicity == "Not Reported"] <- "Other_Or_Unknown"
pros_all$ethnicity[pros_all$ethnicity == "Unknown"] <- "Other_Or_Unknown"
pros_all$ethnicity[is.na(pros_all$ethnicity)] <- "Other_Or_Unknown"

pros_all$race[is.na(pros_all$race)] <- "Other_Or_Unknown"
pros_all$race[pros_all$race == "Not Reported"] <- "Other_Or_Unknown"
pros_all$race[pros_all$race == "Decline to Answer"] <- "Other_Or_Unknown"
pros_all$race[pros_all$race == "Other"] <- "Other_Or_Unknown"
pros_all$race[pros_all$race == "Unknown"] <- "Other_Or_Unknown"

pros_all$language[!(pros_all$language %in% c("English", "Spanish"))] <- "Other_Or_Unknown"
pros_all$insurance_type[pros_all$insurance_type == "Other" | is.na(pros_all$insurance_type)] <- "Other_Or_Unknown"

# Establish slim dataframe with unique study id and other identifiers
adm_and_csn <- pros_all %>% dplyr::select(study_id, pat_enc_csn_id, hsp_account_id, pat_mrn_id, picu_adm_date_time) %>% distinct()

# # Call the AKI script to identify patients with delayed AKI
# source(file = "/home/martinbl/sbi_blake/prospective_validation_files/aki_pros.R")

# # Call the outcomes script to evaluate AKI and C.Diff in patients with and without antibiotic exposure
# source(file = "/home/martinbl/sbi_blake/prospective_validation_files/aae_pros.R")



# ## Calculate duration of antibiotics after sbi negative prediction
# mean((pros_all %>% filter(sbi_present == 0))$abx_duration_after_score) #mean duration of abx for all SBI negative children is 2.22 days.
# mean((pros_all %>% filter(sbi_present == 0) %>% filter(abx_duration_after_score > 0))$abx_duration_after_score) # when they get antibiotics the mean duration is 7.19


# Now calculate number of encounters with SBI and without SBI
n_sbi_with <- pros_all %>%
  group_by(study_id) %>%
  summarise(has_sbi = any(sbi_present == 1), .groups = "drop") %>%
  filter(has_sbi) %>%
  nrow()

n_sbi_without <- pros_all %>%
  group_by(study_id) %>%
  summarise(has_sbi = any(sbi_present == 1), .groups = "drop") %>%
  filter(!has_sbi) %>%
  nrow()


# QC: Output the results
n_sbi_with #702
n_sbi_without #2403
n_sbi_with + n_sbi_without  # This should now equal 3105 which is the number of study_ids that are unique
n_enc <- nrow(pros_all %>% dplyr::select(study_id) %>% distinct()); n_enc # overall 22.6% of picu encounters have an SBI, n_enc = 3104

# Ensure that if a patient has a pneumonia, micro, or cnss SBI that that study_id has all rows with pna present
pros_all <- pros_all %>%
  group_by(study_id) %>%
  mutate(pna_1_0 = if_else(any(pna_1_0 == 1), 1, 0)) %>%
  ungroup()


pros_all <- pros_all %>%
  group_by(study_id) %>%
  mutate(micro_sbi_1_0 = if_else(any(micro_sbi_1_0 == 1), 1, 0)) %>%
  ungroup()

pros_all <- pros_all %>%
  group_by(study_id) %>%
  mutate(ever_cx_neg_sepsis = if_else(any(ever_cx_neg_sepsis == 1), 1, 0)) %>%
  ungroup()

pros_all <- pros_all %>%
  group_by(study_id) %>%
  mutate(sbi_present = if_else(any(sbi_present == 1), 1, 0)) %>%
  ungroup()


# Store / recall pros_all as of 4.29.25
# write.csv(x = pros_all, file = "/phi/sbi/prospective_data/Prospective/pros_all_4_29_25.csv")
# write.csv(x = pros_all, file = "/phi/sbi/prospective_data/Prospective/pros_all_5_13_25.csv") # has antibiotic duration data
# write.csv(x = pros_all, file = "/phi/sbi/prospective_data/Prospective/pros_all_6_20_25.csv") # has antibiotic duration data, data thtur 8/31/25
# write.csv(x = pros_all, file = "/phi/sbi/prospective_data/Prospective/pros_all_6_21_25.csv") # has antibiotic duration data, data through 8/31/25, includes March 2024 data

# pros_all <- read.csv(file = "/phi/sbi/prospective_data/Prospective/pros_all_6_21_25.csv")

# Eval the 2x2 matrix of antibiotic use and SBI presence through 3/31/2024

j_sbi_abx <- pros_all %>% filter(picu_adm_date_time < as.Date("2024-04-01")) %>% dplyr::select(study_id, sbi_present, abx_exp)


j_sbi_abx_summary <- j_sbi_abx %>%
  group_by(study_id) %>%
  summarize(
    sbi_present = max(sbi_present),  # or use `unique()` if you'd prefer to check
    abx_exp = max(abx_exp),
    .groups = "drop"
  )

confusionMatrix(data = as.factor(j_sbi_abx_summary$abx_exp), reference = as.factor(j_sbi_abx_summary$sbi_present))

# Below 2x2 for SBI presence and abx exposure (need to redo with only PICU antibiotics)
# Confusion Matrix and Statistics
#
#               Reference
# Prediction    0     1
#       0       1072  48
#       1       513   419
# We are giving antibiotics to 513 patients with no SBI meaning that 513 / (1072 + 513) * 100 = 32.4% of patients without SBI are getting antibiotics
# Conversely, we're missing 10.3% of SBIs in first 24 hours

# Analyze patients with SBI that didn't get treated
sbi_not_treated <- pros_all %>% filter(sbi_present == 1, abx_exp == 0) %>% filter(picu_adm_date_time < as.Date("2024-04-01"))
sbi_not_treated_slim <- sbi_not_treated %>% dplyr::select(study_id, pat_enc_csn_id, abx_exp, sbi_present, ever_cx_neg_sepsis, micro_sbi_1_0, pna_1_0, picu_adm_date_time,
                                                          age, is_female, malignancy_pccc, pccc, ethnicity) %>% distinct()

# Add row number for aid with including the predictions / results later
pros_old <- pros_all # ensure scores and model types are stored in case needed later

# Select only one models rows to recreate what data a model would be seeing
pros_all <- pros_all %>% filter(model_type == "RF_yes_abx")

pros_all <- pros_all %>% dplyr::select(-model_type, -model_score, -error_id, -pat_id, -line, -time_elapsed, -abx_log_last_rdi_utc, -abx_rf_last_rdi_utc,
                                       -noabx_log_last_rdi_utc, -noabx_rf_last_rdi_utc)
if("X" %in% colnames(pros_all)){pros_all <- pros_all %>% dplyr::select(-X)}
pros_all <- pros_all %>% dplyr::select(-score_time) %>% mutate(hours_since_picu_adm = round(hours_since_picu_adm, digits = 0))
pros_all <- pros_all %>% dplyr::select(-start, -end, -rowid, -t_diff)


nrow(pros_all) #63,173
pros_all <- pros_all %>% distinct()
nrow(pros_all) #62,814


pros_all <- pros_all %>% mutate(rowid = row_number())

# Determine proportion of unnecessary antibiotics and median days of therapy received
# 1. Filter to patients without SBI
no_sbi <- pros_all %>%
  filter(sbi_present == 0)

# 2. Identify those who received antibiotics in the PICU
abx_in_picu <- no_sbi %>%
  group_by(study_id) %>%
  arrange(hours_since_picu_adm) %>%
  mutate(
    abx_started_in_picu = any(abx_exp == 0 & lead(abx_exp, default = last(abx_exp)) == 1),
    continued_abx = any(abx_exp == 1 & abx_duration_after_score > 0)
  ) %>%
  summarise(
    abx_in_picu = any(abx_started_in_picu | continued_abx)
  ) %>%
  filter(abx_in_picu)

# 3. Calculate proportion of all SBI-negative patients who received antibiotics in the PICU
num_no_sbi <- n_distinct(no_sbi$study_id)
num_abx_in_picu <- nrow(abx_in_picu)
prop_abx_in_picu <- num_abx_in_picu / num_no_sbi # 19.8% of all SBI-negative patients got abx in first 24 hrs

# 4. Median days of antibiotics administered (using first value per patient)
abx_days <- no_sbi %>%
  filter(study_id %in% abx_in_picu$study_id, abx_duration_after_score > 0) %>%
  group_by(study_id) %>%
  arrange(hours_since_picu_adm) %>%
  slice(1) %>%
  summarise(first_abx_days = abx_duration_after_score)

median_abx_days <- median(abx_days$first_abx_days, na.rm = TRUE)

# Output
list(
  proportion_of_no_sbi_patients_with_abx = prop_abx_in_picu,
  median_abx_days_among_no_sbi_patients_who_received_abx = median_abx_days
) # median 2.5 days antibiotics

# Plot sample run-chart for proportion of admissions who got unnecessary antibiotics
# First need to ensure we have the patients with SBI abx data
abx_in_picu_all <- pros_all %>%
  group_by(study_id) %>%
  arrange(hours_since_picu_adm) %>%
  mutate(
    abx_started_in_picu = any(abx_exp == 0 & lead(abx_exp, default = last(abx_exp)) == 1),
    continued_abx = any(abx_exp == 1 & abx_duration_after_score > 0)
  ) %>%
  summarise(
    abx_in_picu = any(abx_started_in_picu | continued_abx)
  ) %>%
  filter(abx_in_picu)

ids_sbi_date <- pros_all %>% dplyr::select(study_id, picu_adm_date_time, sbi_present) %>% distinct()
ids_sbi_date <- left_join(ids_sbi_date, abx_in_picu_all, by = "study_id")
ids_sbi_date$abx_in_picu[is.na(ids_sbi_date$abx_in_picu)] <- FALSE
no_sbi_yes_abx <- ids_sbi_date %>% filter(sbi_present == 0) %>% distinct()

# Step 1: Fix time/date and prepare monthly counts
no_sbi_yes_abx$picu_adm_date_time <- ymd_hms(no_sbi_yes_abx$picu_adm_date_time, tz = "America/Denver")
monthly_counts <- no_sbi_yes_abx %>%
  filter(sbi_present == 0, abx_in_picu == TRUE) %>%
  mutate(month = as.Date(floor_date(picu_adm_date_time, unit = "month"))) %>%
  group_by(month) %>%
  summarise(num_patients = n()) %>%
  ungroup()

# Step 2: Plot run chart
ggplot(monthly_counts, aes(x = month, y = num_patients)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "steelblue", size = 2) +
  labs(
    title = "Monthly Count of Unnecessary Antibiotic Use in PICU",
    x = "PICU Admission Month",
    y = "Number of Patients (No SBI + Antibiotics)"
  ) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
  theme_minimal(base_size = 14) +
  theme(
    text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold")
  )

### Now look at proportion of SBI negative patients getting antibiotics instead of raw count
# Step 1: Add month column (Date class)
no_sbi_yes_abx <- no_sbi_yes_abx %>%
  mutate(month = as.Date(floor_date(picu_adm_date_time, unit = "month")))

# Step 2: Calculate monthly numerator (unnecessary antibiotics)
numerator_df <- no_sbi_yes_abx %>%
  filter(sbi_present == 0, abx_in_picu == TRUE) %>%
  group_by(month) %>%
  summarise(num_unnecessary_abx = n())

# Step 3: Calculate monthly denominator (all PICU patients that month)
denominator_df <- no_sbi_yes_abx %>%
  group_by(month) %>%
  summarise(total_patients = n())

# Step 4: Join numerator and denominator, calculate proportion
prop_df <- left_join(numerator_df, denominator_df, by = "month") %>%
  mutate(prop_unnecessary_abx = num_unnecessary_abx / total_patients)

# Step 5: Plot
ggplot(prop_df, aes(x = month, y = prop_unnecessary_abx)) +
  geom_line(color = "firebrick", linewidth = 1) +
  geom_point(color = "firebrick", size = 2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
  labs(
    title = "Monthly Proportion of PICU Patients Without SBI Receiving Antibiotics",
    x = "PICU Admission Month",
    y = "Proportion of Patients (No SBI + Antibiotics)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold")
  )


# Make QI Chart instead
library(qicharts2)

# Ensure prop_df has:
# - x: time (month)
# - n: denominator (total_patients)
# - y: numerator (num_unnecessary_abx)

# Reconstruct if needed:
p_chart_data <- left_join(numerator_df, denominator_df, by = "month") %>%
  mutate(
    month = as.Date(month),
    prop_unnecessary_abx = num_unnecessary_abx / total_patients
  )

# Step 3: Create the p-chart
qic(
  x = month,
  y = num_unnecessary_abx,
  n = total_patients,
  data = p_chart_data,
  chart = "p",
  y.percent = TRUE,
  x.format = "%b %Y",
  title = "p-Chart: Monthly Proportion of SBI-Negative PICU Patients Receiving Antibiotics",
  xlab = "PICU Admission Month",
  ylab = "Proportion of Patients (No SBI + Antibiotics)"
)



#### Now will look at proportion of children with SBI who got abx in first 24 hours

# Step 1: Use the same cleaned dataset, ensure month column exists
sbi_yes_df <- ids_sbi_date %>%
  filter(sbi_present == 1) %>%
  mutate(
    picu_adm_date_time = ymd_hms(picu_adm_date_time, tz = "America/Denver")) %>%
    mutate(month = as.Date(floor_date(picu_adm_date_time, unit = "month")))


# Step 2: Numerator = SBI-positive patients who received antibiotics in PICU
numerator_df <- sbi_yes_df %>%
  filter(abx_in_picu == TRUE) %>%
  group_by(month) %>%
  summarise(num_with_abx = n())

# Step 3: Denominator = All SBI-positive patients
denominator_df <- sbi_yes_df %>%
  group_by(month) %>%
  summarise(total_sbi_patients = n())

# Step 4: Join and calculate proportion
sbi_prop_df <- left_join(numerator_df, denominator_df, by = "month") %>%
  mutate(prop_with_abx = num_with_abx / total_sbi_patients)

# Step 5: Plot
ggplot(sbi_prop_df, aes(x = month, y = prop_with_abx)) +
  geom_line(color = "darkgreen", linewidth = 1) +
  geom_point(color = "darkgreen", size = 2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
  labs(
    title = "Monthly Proportion of SBI-Positive PICU Patients Receiving Antibiotics",
    x = "PICU Admission Month",
    y = "Proportion of SBI Patients Treated with Antibiotics"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold")
  )

# Step 3: Create the p-chart
p_chart_data_sbi <- left_join(numerator_df, denominator_df, by = "month") %>%
  mutate(
    month = as.Date(month),
    prop_abx = num_with_abx / total_sbi_patients
  )

qic(
  x = month,
  y = num_with_abx,
  n = total_sbi_patients,
  data = p_chart_data_sbi,
  chart = "p",
  y.percent = TRUE,
  x.format = "%b %Y",
  title = "p-Chart: Monthly Proportion of SBI+ PICU Patients Receiving Antibiotics",
  xlab = "PICU Admission Month",
  ylab = "Proportion of Patients (SBI Present + Antibiotics)"
)


#### Now will determine what % of patients with SBI had abx started and how quickly ####
# 1. Filter to patients with SBI
sbi_patients <- pros_all %>%
  filter(sbi_present == 1)

# 2. For each study_id, detect:
#    a. If abx_exp goes from 0 to 1 (abx started in PICU)
#    b. If any abx_duration_after_score > 0 while abx_exp == 1 (continued in PICU)
sbi_summary <- sbi_patients %>%
  group_by(study_id) %>%
  arrange(hours_since_picu_adm) %>%
  mutate(
    abx_started_in_picu = any(abx_exp == 0 & dplyr::lead(abx_exp, default = last(abx_exp)) == 1),
    continued_abx = any(abx_exp == 1 & abx_duration_after_score > 0)
  ) %>%
  summarise(
    abx_in_picu = any(abx_started_in_picu | continued_abx)
  )

# 3. Proportion of SBI patients who received antibiotics in the PICU
total_sbi_patients <- n_distinct(sbi_patients$study_id)
sbi_with_abx_in_picu <- sbi_summary %>% filter(abx_in_picu)
num_sbi_with_abx <- nrow(sbi_with_abx_in_picu)
prop_sbi_with_abx <- num_sbi_with_abx / total_sbi_patients

# 4. Median number of hours from PICU admission to abx start (among those who STARTED in PICU)
abx_start_times <- sbi_patients %>%
  filter(study_id %in% sbi_with_abx_in_picu$study_id) %>%
  group_by(study_id) %>%
  arrange(hours_since_picu_adm) %>%
  mutate(abx_shift = abx_exp - dplyr::lag(abx_exp, default = 0)) %>%
  filter(abx_shift == 1) %>%  # captures the first time abx changes from 0 to 1
  slice(1) %>%                # only keep first change per patient
  ungroup()

median_hours_to_abx_start <- median(abx_start_times$hours_since_picu_adm, na.rm = TRUE)
quantile(abx_start_times$hours_since_picu_adm, na.rm = TRUE) # median 1 (0-4)

# Output
list(
  proportion_of_sbi_patients_with_abx_in_picu = prop_sbi_with_abx, # 60.6% of SBI patients got abx in the PICU. in 1st 24 hours
  median_hours_until_abx_started_in_picu = median_hours_to_abx_start # started median 1 hour after PICU admission
)

quantile(abx_start_times$hours_since_picu_adm, na.rm = TRUE)

# Only 60.6% of patients with SBI received antibiotics during that first 24 hours in the PICU. What? Will investigate
# Filter sbi_patients to only those who did NOT receive antibiotics in the PICU
sbi_no_abx_in_picu <- sbi_patients %>% filter(picu_adm_date_time < as.Date("2024-04-01")) %>%
  filter(!(study_id %in% sbi_with_abx_in_picu$study_id)) %>% dplyr::select(study_id, sbi_present, micro_sbi_1_0, pna_1_0, ever_cx_neg_sepsis) %>% distinct()

n_sbi_no_tx <- sbi_no_abx_in_picu %>% dplyr::select(study_id) %>% n_distinct() # 188 patients not treated

sbi_no_tx_summary <- data.frame("Total" = n_sbi_no_tx, "PNA" = sum(sbi_no_abx_in_picu$pna_1_0), "Micro" = sum(sbi_no_abx_in_picu$micro_sbi_1_0),
                                "CNSS" = sum(sbi_no_abx_in_picu$ever_cx_neg_sepsis))
View(sbi_no_tx_summary)
# 175 / 246 = 71.1% are pneumonias






### Train new models using prospective data ###
# slim down to just the study_id, predictor info, and sbi outcome for model training

model_slim <- pros_all %>% dplyr::select(-pat_enc_csn_id, -picu_adm_date_time, -hours_since_picu_adm, -intime, -outtime,
                                                 -contact_creation_to_admission_delta, -height, -height_date, -origin, -valid_start_instant, -valid_end_instant, -cx_neg_sepsis,
                                         -pna_1_0, -bcx_sent, -micro_sbi_1_0, -cnss_true, -ever_cx_neg_sepsis, -hsp_account_id, -pat_mrn_id, -race, -ethnicity,
                                         -language, -insurance_type, -icu_start_instant, -old_race, -prior_unit, -weight, -weight_date, -is_interventional_radiology,
                                         -abx_duration_after_score
                                         )

raw_names <- colnames(model_slim)[str_ends(colnames(model_slim), "_raw")]

model_slim <- model_slim %>% dplyr::select(-all_of(raw_names))

# Re-order
model_slim <- model_slim %>% relocate(study_id, sbi_present, everything())

# Remove duplicate rows
model_slim <- model_slim %>% distinct()


# Now some QC to see IDs where SBI-present changes:
# IDs with mixed SBI labels
mixed_ids <- model_slim %>%
  group_by(study_id) %>%
  filter(n_distinct(sbi_present) > 1) %>%      # TRUE if both 0 and 1 seen
  ungroup() # currently no rows

## Remove bad columns
predictors <- setdiff(names(model_slim), c("study_id", "sbi_present", "SBI", "abx_exp", "rowid"))

nzv <- nearZeroVar(model_slim[ , predictors])
if(length(nzv) > 0) {
  predictors <- predictors[-nzv]
} # gets rid of 28 predictors

bad_cols <- names(model_slim)[
  sapply(model_slim[ , predictors], function(x) any(is.na(x) | is.infinite(x)))
] # no bad columns


set.seed(2025)

## Collapse to one row per admission **first** if you only need admission-level labels
admission_df <- model_slim %>%                 # getting one study id per encounter
  group_by(study_id) %>% dplyr::slice(1) %>% ungroup()

## Create a vector of unique IDs
all_ids <- unique(admission_df$study_id)

## Shuffle the IDs
shuffled_ids <- sample(all_ids)

## Determine sizes
n_total <- length(shuffled_ids)
n_train <- floor(0.60 * n_total)
n_valid <- floor(0.20 * n_total)

## Split into train / validation / test
train_ids <- shuffled_ids[1:n_train]
valid_ids <- shuffled_ids[(n_train + 1):(n_train + n_valid)]
test_ids  <- shuffled_ids[(n_train + n_valid + 1):n_total]

## Build the dataframes
train_df <- model_slim %>% filter(study_id %in% train_ids)
valid_df <- model_slim %>% filter(study_id %in% valid_ids)
test_df  <- model_slim %>% filter(study_id %in% test_ids)


set.seed(2025)

train_df$SBI <- factor(ifelse(train_df$sbi_present == 1, "pos", "neg"), levels = c("pos","neg"))
valid_df$SBI <- factor(ifelse(valid_df$sbi_present == 1, "pos", "neg"), levels = c("pos","neg"))
test_df$SBI <- factor(ifelse(test_df$sbi_present == 1, "pos", "neg"), levels = c("pos","neg"))

# Make sure to create factors
train_df  <- train_df  %>% mutate(across(where(is.character), as.factor))
valid_df <- valid_df%>% mutate(across(where(is.character), as.factor))
test_df<- test_df%>% mutate(across(where(is.character), as.factor))



# Train and tune
set.seed(2025)

library(doParallel)
registerDoParallel(cores = parallel::detectCores() - 1)  # leave 1 core free

# Train initial model
## 1. Build patient-level folds ----------------------------
# one vector of the grouping variable (length == nrow(train_df))
group_var  <- train_df$study_id

# groupKFold() returns a *list*; each element is the rows that go into the
# *training* portion of that fold (the held-out rows are taken as the complement)
folds <- groupKFold(group_var, k = 5)           # same k as your CV

## 2. Tell caret to use them -------------------------------
ctrl <- trainControl(method        = "cv",
                     number        = 5,         # ignored when you give `index`
                     summaryFunction = twoClassSummary,
                     classProbs    = TRUE,
                     verboseIter   = TRUE,
                     index         = folds)     # <-- key line

## 3. Fit the model ----------------------------------------
rf_grid <- expand.grid(
  mtry            = c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 69),
  splitrule       = c("gini", "extratrees"),
  min.node.size   = c(1:2)
)


rf_init <- train(SBI ~ .,
                 tuneGrid = rf_grid,
                 data       = train_df[, c(predictors, "SBI")],
                 method     = "ranger",
                 num.trees  = 1000,
                 importance = "impurity",
                 respect.unordered.factors = TRUE,
                 na.action  = na.omit,
                 trControl  = ctrl) # the ntree = 100 results are: 15, extratrees, 1, AUC 0.80



# rf_tune <- train(SBI ~ .,
#                  tuneGrid = rf_grid,
#                  data       = train_df[, c(predictors, "SBI")],
#                  method     = "ranger",
#                  num.trees  = 100,
#                  importance = "none",
#                  respect.unordered.factors = TRUE,
#                  na.action  = na.omit,
#                  trControl  = ctrl)  # yields mtry 34, extratrees and mns of 2

rf_tune <- rf_init

# # Only using the identified optimal values
# rf_grid <- expand.grid(
#   mtry            = c(23),
#   splitrule       = c("extratrees"),
#   min.node.size   = c(3)
# )
# rf_tune <- train(SBI ~ ., tuneGrid = rf_grid, data = (train_df %>% dplyr::select(-study_id, -sbi_present)) , method = "ranger", na.action = na.omit, importance = 'none',
#                  respect.unordered.factors = TRUE, num.trees = 100,
#                  trControl = trainControl(method = "cv", number = 5, summaryFunction = twoClassSummary, classProbs = TRUE,
#                                           verboseIter = TRUE))


### Optional code to plot top 40 predictors with scaled variable importance
# 1. pull raw importance from the caret model
imp_tbl <- varImp(rf_init, scale = FALSE)$importance %>%
  tibble::rownames_to_column("predictor") %>%   # keep predictor names
  rename(raw_imp = Overall)

# 2. linearly rescale to 0–100 (100 = most important)
imp_tbl <- imp_tbl %>%
  mutate(scaled_imp = 100 * raw_imp / max(raw_imp, na.rm = TRUE))

# 3. select the top 40 predictors
top40 <- imp_tbl %>%
  arrange(desc(scaled_imp)) %>%
  slice_head(n = 40)

# 4. horizontal bar plot
ggplot(top40,
       aes(x = reorder(predictor, scaled_imp), y = scaled_imp)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Random-forest variable importance (top 40)",
       x = NULL,
       y = "Scaled importance (0–100)") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank())



# Apply to training and test sets
rf_pred_prob_train <- predict(rf_tune, train_df %>% dplyr::select(-study_id, -sbi_present), type = "prob")[, "pos"]
rf_pred_prob_val <- predict(rf_tune, valid_df %>% dplyr::select(-study_id, -sbi_present), type = "prob")[, "pos"]
rf_pred_prob_test <- predict(rf_tune, test_df %>% dplyr::select(-study_id, -sbi_present), type = "prob")[, "pos"]

colAUC(rf_pred_prob_val, valid_df$SBI, plotROC = TRUE) # determine AUROC in test set = 0.82
colAUC(rf_pred_prob_test, test_df$SBI, plotROC = TRUE) # determine AUROC in test set = 0.85

# Identify best threshold for NPV in validation set
# Loop over thresholds
thresholds <- seq(0, 1, by = 0.01)
npv_values <- sapply(thresholds, function(thresh) {
  preds_val <- ifelse(rf_pred_prob_val > thresh, "pos", "neg")
  preds_val <- factor(preds_val, levels = c("neg", "pos"))

  # Ensure actuals are character-levels matching the confusion matrix access
  actuals <- valid_df$SBI


  cm <- table(preds_val, actuals)

  tn <- cm["neg", "neg"]
  fn <- cm["neg", "pos"]

  npv <- tn / (tn + fn)
  return(npv)
})

# Likely doesn't make sense to look at NPV thresholds in training set
# npv_values_train <- sapply(thresholds, function(thresh) {
#   preds <- ifelse(rf_pred_prob_train > thresh, "pos", "neg")
#   cm <- table(preds, train_df$SBI)
#   tn <- cm["neg", "neg"]
#   fn <- cm["neg", "pos"]
#   npv <- tn / (tn + fn)
#   return(npv)
# })

View(data.frame("threshold" = thresholds, "npv" = npv_values)) #yields 0.12 as highest threshold that still yields 96% NPV in training set


plot(x = thresholds, y = npv_values)

my_thresh <- 0.12 #identified as highest threshold with 96% NPV in validation set, now test in test set
preds_test <- ifelse(rf_pred_prob_test > my_thresh, "pos", "neg")
cm <- table(preds_test, test_df$SBI) # identifies 2810 negative rows = 38.6% of negative rows
tn <- cm["neg", "neg"]
fn <- cm["neg", "pos"]
npv <- tn / (tn + fn)
npv #0.95 in test set

############## Repeat this process to identify PPV > 0.8 in the validation set and then test in test set
# Step 1: Loop through thresholds to calculate PPV on validation set
thresholds <- seq(0, 1, by = 0.01)
ppv_values <- sapply(thresholds, function(thresh) {
  preds_val <- ifelse(rf_pred_prob_val > thresh, "pos", "neg")
  preds_val <- factor(preds_val, levels = c("neg", "pos"))

  actuals <- valid_df$SBI

  cm <- table(preds_val, actuals)

  tp <- cm["pos", "pos"]
  fp <- cm["pos", "neg"]

  if ((tp + fp) == 0) return(NA)  # Avoid division by zero

  ppv <- tp / (tp + fp)
  return(ppv)
})

# View threshold vs PPV
ppv_df <- data.frame(threshold = thresholds, ppv = ppv_values)
View(ppv_df)

# Optional: plot PPV curve
plot(x = thresholds, y = ppv_values, type = "l", main = "PPV vs. Threshold", xlab = "Threshold", ylab = "PPV")

# Step 2: Choose highest threshold with PPV > 0.8
valid_thresh_ppv <- min(ppv_df$threshold[ppv_df$ppv > 0.8], na.rm = TRUE) # yields a threshold of 0.66 to achieve PPV > 0.8 in validation set

# Step 3: Evaluate that threshold in the test set
preds_test_ppv <- ifelse(rf_pred_prob_test > valid_thresh_ppv, "pos", "neg")
preds_test_ppv <- factor(preds_test_ppv, levels = c("neg", "pos"))
actuals_test_ppv <- test_df$SBI

cm_ppv <- table(preds_test_ppv, actuals_test_ppv) # identifies 491 of 2741 = 18% of the true positives

tp_ppv <- cm_ppv["pos", "pos"]
fp_ppv <- cm_ppv["pos", "neg"]

ppv_test <- tp_ppv / (tp_ppv + fp_ppv)
ppv_test # PPV of 0.90 in test set with threshold of 0.66

levels(test_df$SBI)

#### Create calibration plot ####
test_df <- test_df %>%
  mutate(pred_prob = rf_pred_prob_test)

# Create bins (e.g., deciles of predicted probability)
test_df <- test_df %>%
  mutate(prob_bin = cut(pred_prob,
                        breaks = seq(0, 1, by = 0.1),
                        include.lowest = TRUE,
                        right = FALSE,
                        labels = paste0(seq(0, 90, by = 10), "-", seq(10, 100, by = 10), "%")))

# Summarize observed SBI rate in each bin
bin_summary <- test_df %>%
  group_by(prob_bin) %>%
  summarize(mean_pred_prob = mean(pred_prob),
            observed_sbi_rate = mean(sbi_present),
            count = n(),
            .groups = "drop")

# Plot the calibration curve
ggplot(bin_summary, aes(x = mean_pred_prob, y = observed_sbi_rate)) +
  geom_point(size = 3) +
  geom_line() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  labs(title = "Calibration Plot: Predicted vs. Observed SBI Rate",
       x = "Mean Predicted SBI Probability (per bin)",
       y = "Observed Proportion with SBI") +
  theme_minimal()




just_outcomes <- pros_all %>% dplyr::select(rowid, study_id, hours_since_picu_adm, sbi_present, pna_1_0, micro_sbi_1_0, ever_cx_neg_sepsis, abx_exp, abx_duration_after_score) %>% distinct()
test_w_preds <- test_df %>% left_join(just_outcomes, by = "rowid")
test_w_preds <- test_w_preds %>% bind_cols(data.frame("prediction" = preds_test))
test_w_preds <- test_w_preds %>% mutate("pred_prob" = rf_pred_prob_test)
test_w_preds <- test_w_preds %>% relocate(rowid, study_id.x, study_id.y, hours_since_picu_adm, sbi_present.x, sbi_present.y, prediction, pred_prob, abx_duration_after_score, pna_1_0, micro_sbi_1_0, ever_cx_neg_sepsis)

# remove extra columns and rename
test_w_preds <- test_w_preds %>% dplyr::select(-study_id.y, -sbi_present.x) %>% rename(study_id = study_id.x, sbi_present = sbi_present.y)

# Plot some trajectories
# 1. Filter to 0–24 hours (if needed) and sample 10 random patients
plot_seed <- 2024
set.seed(plot_seed)  # for reproducibility, 37 = best, other good options 2029, 16
sample_ids <- test_w_preds %>%
  filter(hours_since_picu_adm >= 0 & hours_since_picu_adm <= 24) %>% dplyr::select(study_id) %>%
  distinct(study_id) %>%
  slice_sample(n = 8) %>%
  pull(study_id)

# 2. Subset data to just these patients
plot_data <- test_w_preds %>%
  filter(study_id %in% sample_ids, hours_since_picu_adm <= 24)

# 3. Plot
ggplot(plot_data, aes(x = hours_since_picu_adm, y = pred_prob, group = study_id, color = as.factor(sbi_present))) +
  geom_line(linewidth = 1) +

  # Add horizontal lines
  geom_hline(yintercept = 0.66, color = "black", linetype = "solid") +
  geom_hline(yintercept = 0.12, color = "black", linetype = "solid") +

  # Add labels for the threshold lines
  annotate("text", x = 24, y = 0.66 + 0.005, label = "Prob = 0.66", hjust = 1, vjust = 0, size = 5, fontface = "bold") +
  annotate("text", x = 24, y = 0.12 + 0.005, label = "Prob = 0.12", hjust = 1, vjust = 0, size = 5, fontface = "bold") +

  # Color by SBI status
  scale_color_manual(
    values = c("0" = "blue", "1" = "red"),
    labels = c("No SBI", "SBI"),
    name = "SBI Status"
  ) +

  labs(
    title = "SBI Prediction Trajectories Over First 24 Hours in PICU",
    x = "Hours Since PICU Admission",
    y = "Predicted Probability of SBI"
  ) +

  # Bold and enlarge all text
  theme_minimal(base_size = 16) +  # base font size
  theme(
    legend.position = "bottom",
    text = element_text(face = "bold"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 16, face = "bold")
  )


# Plot of SHAP values for a random SBI patient predicted to be SBI negative at 6 hour mark
# Step 1: Filter for eligible patients
eligible_patients <- test_w_preds %>%
  filter(sbi_present == 0,
         hours_since_picu_adm == 6,
         pred_prob < 0.12, abx_exp.x == 1)

# Step 2: Randomly select one study_id
set.seed(2028)
selected_study_id <- sample(unique(eligible_patients$study_id), 1)

# Step 3: Extract full feature vector for that patient at hour 6
patient_row <- test_w_preds %>%
  filter(study_id == selected_study_id, hours_since_picu_adm == 6)

# Step 4: Prepare data for SHAP
# Assuming your model was trained using a set of predictors named `predictors`
# and your rf_tune model was trained with formula interface like: SBI ~ .
# If you have `predictors` stored as a vector, use that.
patient_features <- patient_row %>% dplyr::select(all_of(predictors))

# Step 5: Create iml Predictor object
predictor_iml <- Predictor$new(
  model = rf_tune,
  data = train_df %>% dplyr::select(dplyr::all_of(predictors)),  # use your training features here
  y = train_df$SBI,                                # binary outcome
  type = "prob",                                   # model outputs probabilities
  class = "pos"                                        # class index for positive class (SBI = 1)
)

# Step 6: Compute SHAP values
shap <- Shapley$new(predictor_iml, x.interest = patient_features)

# Step 7: Get and plot top 10 SHAP values
shap_values <- shap$results %>%
  arrange(desc(abs(phi))) %>%
  slice(1:10)

# Plot
ggplot(shap_values, aes(x = reorder(feature.value, phi), y = phi, fill = phi > 0)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  labs(
    title = "Top 10 SHAP Values for a Low-Risk, SBI-Negative Patient",
    x = "Feature and Value",
    y = "SHAP Value"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14, face = "bold"),
    text = element_text(face = "bold")
  )




# Plot Distribution of SBI scores by SBI status
library(ggplot2)

ggplot(test_w_preds, aes(x = pred_prob, fill = as.factor(sbi_present))) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(
    values = c("0" = "blue", "1" = "red"),
    labels = c("No SBI", "SBI"),
    name = "SBI Status"
  ) +
  labs(
    title = "Distribution of Predicted SBI Probabilities",
    x = "Predicted Probability of SBI",
    y = "Density"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    text = element_text(face = "bold"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14, face = "bold"),
    legend.position = "top"
  )





true_negs <- test_w_preds %>% left_join( (pros_all %>% dplyr::select(study_id, picu_adm_date_time) %>% distinct()), by = "study_id") %>%
  filter(sbi_present == 0) %>% filter(prediction == "neg") %>% filter(picu_adm_date_time < as.Date("2024-04-01"))
true_negs$abx_duration_after_score[true_negs$abx_duration_after_score > 7] <- 7

true_negs$abx_duration_after_score[true_negs$abx_duration_after_score == 0] <- NA
median(true_negs$abx_duration_after_score, na.rm = T) # median 3.0 days antibiotics
mean(true_negs$abx_duration_after_score, na.rm = T) # mean 3.1 antibiotic days for those who get abx
stats::sd(true_negs$abx_duration_after_score, na.rm = T) # sd 2.0
quantile(true_negs$abx_duration_after_score, na.rm = T)
# median 3.0 and IQR 1.4-4.5 days


# Determine clinical relevance of missed sbis
n_test_predictions <- nrow(test_w_preds) # 12,697
n_test_patients <- test_w_preds %>% dplyr::select(study_id) %>% n_distinct() # 622 patients
n_neg_predictions <- sum(test_w_preds$prediction == "neg") # 4881
n_wrong_neg_predictions <- nrow(test_w_preds %>% filter(prediction == "neg", sbi_present == 1)) #249

fn_pred_rows <- test_w_preds %>% filter(prediction == "neg", sbi_present == 1)
n_wrong_patients <- fn_pred_rows %>% dplyr::select(study_id) %>% n_distinct() #29 encounters of 622

n_predictions_wrong_pna <- nrow(fn_pred_rows %>% filter(pna_1_0 == 1)) #148 wrong predictions in pna patients
n_patients_wrong_pna <- nrow(fn_pred_rows %>% filter(pna_1_0 == 1) %>% dplyr::select(study_id) %>% distinct()) #18 patients with pneumonia & FN pred at some point

n_predictions_wrong_cnss <- nrow(fn_pred_rows %>% filter(ever_cx_neg_sepsis == 1)) #7 wrong predictions in culture negative sepsis patients
n_patients_wrong_cnss <- nrow(fn_pred_rows %>% filter(ever_cx_neg_sepsis == 1) %>% dplyr::select(study_id) %>% distinct()) # 5 patients with cnss & FN pred at some point

n_predictions_wrong_micro <- nrow(fn_pred_rows %>% filter(micro_sbi_1_0 == 1)) # 99 wrong predictions in culture negative sepsis patients
n_patients_wrong_micro <- nrow(fn_pred_rows %>% filter(micro_sbi_1_0 == 1) %>% dplyr::select(study_id) %>% distinct()) # 9 patients with micro confirmed infxn & FN pred at some point


################## Now redo most of the above but with only the top 40 variables from original abx_exp model
predictors_ae <- c("age", "hematocrit_mean", "lactate_max", "dbp_last", "fio2_last", "hr_last", "sbp_last", "temp_last", "los_before_icu_days", "dbp_max",
                   "hr_max", "sbp_max", "temp_max", "dbp_mean", "fio2_mean", "hr_mean", "rr_mean", "o2sat_mean", "sbp_mean", "temp_mean",
                   "dbp_median", "fio2_median", "hr_median", "rr_median", "sbp_median", "temp_median", "dbp_min", "rr_min", "hr_min", "o2sat_min",
                   "sbp_min", "o2sat_count", "preicu_6thfloor", "preicu_8thfloor", "preicu_9thfloor", "preicu_er", "preicu_hem_onc_bmt", "preicu_nocer",
                   "preicu_nocinpatient", "preicu_nonchco", "preicu_or", "preicu_othericu", "preicu_outpatient", "preicu_procedure_center", "dbp_slope", "hr_slope",
                   "rr_slope", "o2sat_slope", "sbp_slope", "temp_slope")



# Train and tune
set.seed(2025)

library(doParallel)
registerDoParallel(cores = parallel::detectCores() - 1)  # leave 1 core free

# Train initial model
## 1. Build patient-level folds ----------------------------
# one vector of the grouping variable (length == nrow(train_df))
group_var  <- train_df$study_id

# groupKFold() returns a *list*; each element is the rows that go into the
# *training* portion of that fold (the held-out rows are taken as the complement)
folds <- groupKFold(group_var, k = 5)           # same k as your CV

## 2. Tell caret to use them -------------------------------
ctrl <- trainControl(method        = "cv",
                     number        = 5,         # ignored when you give `index`
                     summaryFunction = twoClassSummary,
                     classProbs    = TRUE,
                     verboseIter   = TRUE,
                     index         = folds)     # <-- key line

## 3. Fit the model ----------------------------------------
rf_grid_ae <- expand.grid(
  mtry            = c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50),
  splitrule       = c("gini", "extratrees"),
  min.node.size   = c(1:2)
)

rf_init_ae_40 <- train(SBI ~ .,
                 data       = train_df[, c(predictors_ae, "SBI")],
                 method     = "ranger",
                 num.trees  = 1000,
                 importance = "impurity",
                 respect.unordered.factors = TRUE,
                 na.action  = na.omit,
                 trControl  = ctrl) #2 gini and 1 are best, AUC 0.775

rf_grid_ae <- expand.grid(
  mtry            = c(1:3, 25:27, 51:69),
  splitrule       = c("gini", "extratrees"),
  min.node.size   = c(1:4)
)

# If already did grid search
# rf_tune_ae_40 <- rf_init_ae_40

rf_tune_ae_40 <- train(SBI ~ .,
                 tuneGrid = rf_grid_ae,
                 data       = train_df[, c(predictors, "SBI")],
                 method     = "ranger",
                 num.trees  = 1000,
                 importance = "none",
                 respect.unordered.factors = TRUE,
                 na.action  = na.omit,
                 trControl  = ctrl)  # yields mtry 34, extratrees and mns of 2

# # Only using the identified optimal values
# rf_grid <- expand.grid(
#   mtry            = c(23),
#   splitrule       = c("extratrees"),
#   min.node.size   = c(3)
# )
# rf_tune <- train(SBI ~ ., tuneGrid = rf_grid, data = (train_df %>% dplyr::select(-study_id, -sbi_present)) , method = "ranger", na.action = na.omit, importance = 'none',
#                  respect.unordered.factors = TRUE, num.trees = 100,
#                  trControl = trainControl(method = "cv", number = 5, summaryFunction = twoClassSummary, classProbs = TRUE,
#                                           verboseIter = TRUE))

# Apply to training and test sets
rf_pred_prob_train_ae_40 <- predict(rf_tune_ae_40, train_df %>% dplyr::select(-study_id, -sbi_present), type = "prob")[, "pos"]
rf_pred_prob_val_ae_40 <- predict(rf_tune_ae_40, valid_df %>% dplyr::select(-study_id, -sbi_present), type = "prob")[, "pos"]
rf_pred_prob_test_ae_40 <- predict(rf_tune_ae_40, test_df %>% dplyr::select(-study_id, -sbi_present), type = "prob")[, "pos"]

colAUC(rf_pred_prob_val_ae_40, valid_df$SBI, plotROC = TRUE) # determine AUROC in test set = 0.76
colAUC(rf_pred_prob_test_ae_40, test_df$SBI, plotROC = TRUE) # determine AUROC in test set = 0.781

# Identify best threshold for NPV in validation set
# Loop over thresholds
thresholds_ae_40 <- seq(0, 1, by = 0.01)
npv_values_ae_40 <- sapply(thresholds, function(thresh) {
  preds_val_ae_40 <- ifelse(rf_pred_prob_val_ae_40 > thresh, "pos", "neg")
  preds_val_ae_40 <- factor(preds_val_ae_40, levels = c("neg", "pos"))

  # Ensure actuals are character-levels matching the confusion matrix access
  actuals <- valid_df$SBI


  cm_ae_40 <- table(preds_val_ae_40, actuals)

  tn <- cm_ae_40["neg", "neg"]
  fn <- cm_ae_40["neg", "pos"]

  npv <- tn / (tn + fn)
  return(npv)
})

# Likely doesn't make sense to look at NPV thresholds in training set
# npv_values_train <- sapply(thresholds, function(thresh) {
#   preds <- ifelse(rf_pred_prob_train > thresh, "pos", "neg")
#   cm <- table(preds, train_df$SBI)
#   tn <- cm["neg", "neg"]
#   fn <- cm["neg", "pos"]
#   npv <- tn / (tn + fn)
#   return(npv)
# })

View(data.frame("threshold" = thresholds, "npv" = npv_values_ae_40)) #yields 0.10 as highest threshold that still yields 95% NPV in training set

plot(x = thresholds_ae_40, y = npv_values_ae_40)

my_thresh_ae_40 <- 0.1 #identified as highest threshold with 95% NPV in validation set, now test in test set
preds_test_ae_40 <- ifelse(rf_pred_prob_test_ae_40 > my_thresh_ae_40, "pos", "neg")
cm_ae_40 <- table(preds_test_ae_40, test_df$SBI) # identifies 1718 / 7282 = 24% of total TNs
tn <- cm_ae_40["neg", "neg"]
fn <- cm_ae_40["neg", "pos"]
npv_ae_40 <- tn / (tn + fn)
npv_ae_40 #0.93 in test set



######################## Now repeat process with a logistic regression model ############################
## 1. Fit the model (no tuning parameters for plain GLM)
set.seed(2025)
logit_model <- train(
  SBI ~ .,
  data      = train_df[, c(predictors, "SBI")],
  method    = "glm",
  family    = binomial,            # logistic link
  na.action = na.omit,
  trControl = ctrl                 # AUC 0.78 in training set
)

## 2. Predicted probabilities
logit_pred_prob_train <- predict(
  logit_model,
  newdata = train_df %>% dplyr::select(-study_id, -sbi_present),
  type    = "prob")[ , "pos"]

logit_pred_prob_val <- predict(
  logit_model,
  newdata = valid_df %>% dplyr::select(-study_id, -sbi_present),
  type    = "prob")[ , "pos"]

logit_pred_prob_test <- predict(
  logit_model,
  newdata = test_df  %>% dplyr::select(-study_id, -sbi_present),
  type    = "prob")[ , "pos"]

## 3. AUROC in validation & test sets
colAUC(logit_pred_prob_val,  valid_df$SBI,  plotROC = TRUE)   # AUROC 0.77
colAUC(logit_pred_prob_test, test_df$SBI,  plotROC = TRUE)    # AUROC 0.80

## 4. Pick the highest threshold that still gives ≥95 % NPV in the *validation* set
thresholds <- seq(0, 1, by = 0.01)

# make “actuals” a factor whose levels match the confusion-matrix access later
actuals_val_lr <- valid_df$SBI

npv_values_lr <- sapply(thresholds, function(thresh) {
  preds_val_lr <- factor(ifelse(logit_pred_prob_val > thresh, "pos", "neg"),
                      levels = c("neg", "pos"))
  cm_lr <- table(preds_val_lr, actuals_val_lr)

  tn <- cm_lr["neg", "neg"]
  fn <- cm_lr["neg", "pos"]
  tn / (tn + fn)            # NPV
})

# inspect or plot if you like
View(data.frame(threshold = thresholds, npv = npv_values_lr))

best_thresh_lr <- max(thresholds[npv_values_lr >= 0.95])  # highest threshold with NPV ≥ 0.95
best_thresh_lr <- 0.04

## 5. Apply that threshold to the *test* set
actuals_test <- test_df$SBI
preds_test   <- factor(ifelse(logit_pred_prob_test > best_thresh_lr, "pos", "neg"),
                       levels = c("neg", "pos"))
cm_test <- table(preds_test, actuals_test)

tn <- cm_test["neg", "neg"]
fn <- cm_test["neg", "pos"]
npv_test <- tn / (tn + fn) # NPY 0.98 and identifies 16 percent of SBI-negative children

cat(sprintf("Best threshold from validation = %.2f\n", best_thresh))
cat(sprintf("NPV in test set at that threshold = %.2f\n", npv_test)) # NPV = 0.98
























##### Now create separate dataframes for the 4 models using the right abx exposure cohort
no_abx_rf <- pros_all %>% filter(abx_exp == 0) %>% filter(model_type == "RF_no_abx")
yes_abx_rf <- pros_all %>% filter(abx_exp == 1) %>% filter(model_type == "RF_yes_abx")
no_abx_lr <- pros_all %>% filter(abx_exp == 0) %>% filter(model_type == "LR_no_abx")
yes_abx_lr <- pros_all %>% filter(abx_exp == 1) %>% filter(model_type == "LR_yes_abx")

# Establish cut-point
cut_point_yes_abx <- 10
cut_point_no_abx <- 5

# Compute confusion matrices using appropriate cut-points
yes_abx_lr <- yes_abx_lr %>% mutate(pred_1_0 = ifelse(model_score <= cut_point_yes_abx, yes = 0, no = 1))
yes_abx_rf <- yes_abx_rf %>% mutate(pred_1_0 = ifelse(model_score <= cut_point_yes_abx, yes = 0, no = 1))

no_abx_lr <- no_abx_lr %>% mutate(pred_1_0 = ifelse(model_score <= cut_point_no_abx, yes = 0, no = 1))
no_abx_rf <- no_abx_rf %>% mutate(pred_1_0 = ifelse(model_score <= cut_point_no_abx, yes = 0, no = 1))


# Create the confusion matrices using the confusionMatrix function
cm_yes_abx_lr <- confusionMatrix(as.factor(yes_abx_lr$pred_1_0), as.factor(yes_abx_lr$sbi_present), positive = '1')
cm_yes_abx_rf <- confusionMatrix(as.factor(yes_abx_rf$pred_1_0), as.factor(yes_abx_rf$sbi_present), positive = '1')
cm_no_abx_lr <- confusionMatrix(as.factor(no_abx_lr$pred_1_0), as.factor(no_abx_lr$sbi_present), positive = '1')
cm_no_abx_rf <- confusionMatrix(as.factor(no_abx_rf$pred_1_0), as.factor(no_abx_rf$sbi_present), positive = '1')

# Print out the new CM's
cm_yes_abx_lr
cm_yes_abx_rf
cm_no_abx_lr
cm_no_abx_rf

# Call the create_models_val script to create models from the original training set (RF)
source(file = "/home/martinbl/sbi_blake/create_models_for_val.R")

# Use rf-model from retro data (testing set) to get the sbi probabilities
sbi_prob_retrospective <- exp_test_df_au$prob_of_sbi_by_model * 100 #1731 rows in the abx unexposed test set

# 2. Take the prospective probabilities from your dataset
sbi_prob_prospective <- no_abx_rf$model_score

# 3. Combine into a single data frame for plotting
plot_df <- data.frame(
  prob  = c(sbi_prob_retrospective, sbi_prob_prospective),
  group = c(
    rep("Retrospective", length(sbi_prob_retrospective)),
    rep("Prospective",   length(sbi_prob_prospective))
  )
)

# 4. Plot the distributions side by side (density plot)
ggplot(plot_df, aes(x = prob, color = group, fill = group)) +
  geom_density(alpha = 0.3) +
  scale_color_manual(
    name = "Cohort",  # Legend title
    values = c("Retrospective" = "red", "Prospective" = "blue")
  ) +
  scale_fill_manual(
    name = "Cohort",  # Legend title
    values = c("Retrospective" = "red", "Prospective" = "blue")
  ) +
  labs(
    title = "Distribution of Random Forest Predicted Probabilities: Abx Unexposed",
    x     = "Predicted Probability of SBI",
    y     = "Density"
  ) +
  geom_vline(xintercept = 5, color = "black", linetype = "dashed") +
  scale_x_continuous(limits = c(0, 100)) +
  theme_minimal() +
  theme(
    # Make all text larger and bold
    text              = element_text(size = 14, face = "bold"),
    axis.title        = element_text(size = 14, face = "bold"),
    axis.text         = element_text(size = 12, face = "bold"),
    legend.title      = element_text(size = 14, face = "bold"),
    legend.text       = element_text(size = 12, face = "bold"),
    plot.title        = element_text(size = 16, face = "bold")
  )


###### Now repeat using the abx exposed cohort
# Use rf-model from retro data to get the sbi probabilities
sbi_prob_retrospective_abx <- exp_test_df_ae$prob_of_sbi_by_model * 100

# 2. Take the prospective probabilities from your dataset
sbi_prob_prospective_abx <- yes_abx_rf$model_score

# 3. Combine into a single data frame for plotting
plot_df_abx <- data.frame(
  prob  = c(sbi_prob_retrospective_abx, sbi_prob_prospective_abx),
  group = c(
    rep("Retrospective", length(sbi_prob_retrospective_abx)),
    rep("Prospective",   length(sbi_prob_prospective_abx))
  )
)

# 4. Plot the distributions side by side (density plot)
ggplot(plot_df_abx, aes(x = prob, color = group, fill = group)) +
  geom_density(alpha = 0.3) +
  scale_color_manual(
    name = "Cohort",  # Legend title
    values = c("Retrospective" = "red", "Prospective" = "blue")
  ) +
  scale_fill_manual(
    name = "Cohort",  # Legend title
    values = c("Retrospective" = "red", "Prospective" = "blue")
  ) +
  labs(
    title = "Distribution of Random Forest Predicted Probabilities: Abx Unexposed",
    x     = "Predicted Probability of SBI",
    y     = "Density"
  ) +
  geom_vline(xintercept = 10, color = "black", linetype = "dashed") +
  theme_minimal() +
  theme(
    # Make all text larger and bold
    text              = element_text(size = 14, face = "bold"),
    axis.title        = element_text(size = 14, face = "bold"),
    axis.text         = element_text(size = 12, face = "bold"),
    legend.title      = element_text(size = 14, face = "bold"),
    legend.text       = element_text(size = 12, face = "bold"),
    plot.title        = element_text(size = 16, face = "bold")
  )








# Plot density distributions

ggplot() +
  geom_density(data = no_abx_lr, aes(x = model_score, fill = "Logistic Regression"), alpha = 0.5) +
  geom_density(data = no_abx_rf, aes(x = model_score, fill = "Random Forest"), alpha = 0.5) +
  scale_fill_manual(values = c("Logistic Regression" = "blue", "Random Forest" = "green"), name = "Model Type") +
  geom_vline(aes(xintercept = 5), color = "red", linetype = "dashed", size = 1) +
  labs(title = "Distribution of Model Scores: Antibiotic Unexposed",
       x = "Model Score",
       y = "Density") +
  annotate("text", x = 0.055, y = 0.09, label = "Training Data Cut-Point", color = "red", hjust = -0.2, fontface = "bold", size = 5) +  # Adjusted label positioning
  theme_minimal() +
  theme(text = element_text(face = "bold", size = 16),  # Apply bold and larger text to all labels and titles
        plot.title = element_text(size = 20, face = "bold"),  # Make title larger and bold
        legend.title = element_text(face = "bold", size = 16),  # Bold and larger legend titles
        legend.text = element_text(face = "bold", size = 14)) +  # Bold and larger legend text
  ylim(0, 0.1)  # Set y-axis limits


ggplot() +
  geom_density(data = yes_abx_lr, aes(x = model_score, fill = "Logistic Regression"), alpha = 0.5) +
  geom_density(data = yes_abx_rf, aes(x = model_score, fill = "Random Forest"), alpha = 0.5) +
  scale_fill_manual(values = c("Logistic Regression" = "blue", "Random Forest" = "green"), name = "Model Type") +
  geom_vline(aes(xintercept = 10), color = "red", linetype = "dashed", size = 1) +
  labs(title = "Distribution of Model Scores: Antibiotic Exposed",
       x = "Model Score",
       y = "Density") +
  annotate("text", x = 0.055, y = 0.09, label = "Training Data Cut-Point", color = "red", hjust = -0.35, fontface = "bold", size = 5) +  # Adjusted label positioning
  theme_minimal() +
  theme(text = element_text(face = "bold", size = 16),  # Apply bold and larger text to all labels and titles
        plot.title = element_text(size = 20, face = "bold"),  # Make title larger and bold
        legend.title = element_text(face = "bold", size = 16),  # Bold and larger legend titles
        legend.text = element_text(face = "bold", size = 14)) +  # Bold and larger legend text
  ylim(0, 0.1)  # Set y-axis limits


########### Now start computing other metrics e.g. AUROC, AUPRC, etc ##########################
compute_metrics <- function(df, threshold = 0.5) {
  # Ensure required packages are installed
  if (!requireNamespace("pROC", quietly = TRUE)) {
    install.packages("pROC")
  }
  if (!requireNamespace("PRROC", quietly = TRUE)) {
    install.packages("PRROC")
  }

  # Load packages
  library(pROC)
  library(PRROC)

  # Extract true labels and model scores
  y_true <- df$sbi_present
  y_score <- df$model_score

  # Ensure y_true is binary numeric (0/1)
  y_true <- as.numeric(as.character(y_true))

  # Check if y_true has at least two levels
  if (length(unique(y_true)) < 2) {
    # If there is only one level, return NA for AUROC and AUPRC
    auroc <- NA
    auprc <- NA
  } else {
    # Compute AUROC
    roc_obj <- roc(y_true, y_score)
    auroc <- auc(roc_obj)

    # Compute AUPRC
    pr <- pr.curve(scores.class0 = y_score[y_true == 1],
                   scores.class1 = y_score[y_true == 0],
                   curve = FALSE)
    auprc <- pr$auc.integral
  }

  # Generate binary predictions based on the threshold
  y_pred <- ifelse(y_score >= threshold, 1, 0)

  # Calculate confusion matrix components
  TP <- sum(y_pred == 1 & y_true == 1)
  TN <- sum(y_pred == 0 & y_true == 0)
  FP <- sum(y_pred == 1 & y_true == 0)
  FN <- sum(y_pred == 0 & y_true == 1)

  # Compute performance metrics with checks to avoid division by zero
  sensitivity <- ifelse((TP + FN) > 0, TP / (TP + FN), NA)
  specificity <- ifelse((TN + FP) > 0, TN / (TN + FP), NA)
  ppv <- ifelse((TP + FP) > 0, TP / (TP + FP), NA)
  npv <- ifelse((TN + FN) > 0, TN / (TN + FN), NA)
  nnh <- ifelse(FN > 0, TN / FN, Inf)  # Avoid division by zero

  # Compile metrics into a list
  metrics <- list(
    AUROC = auroc,
    AUPRC = auprc,
    Sensitivity = sensitivity,
    Specificity = specificity,
    PPV = ppv,
    NPV = npv,
    Number_Needed_to_Harm = nnh
  )

  return(metrics)
}


# Calculate the metrics for the models
compute_metrics(df = no_abx_lr, cut_point_no_abx)
compute_metrics(df = no_abx_rf, cut_point_no_abx)

compute_metrics(df = yes_abx_lr, cut_point_yes_abx)
compute_metrics(df = yes_abx_rf, cut_point_yes_abx)

#Call below for prospective result subgroup analysis
source(file = "/home/martinbl/sbi_blake/prospective_validation_files/pros_subgroup_analysis.R")

# Call below to create new random forrest and logistic regression model
source(file = "/home/martinbl/sbi_blake/prospective_validation_files/pros_new_model.R")



######### Quality Control code on pros_all #############
# Plot HR values for different subgroups.

# Load necessary libraries
library(ggplot2)

# Create the plot
ggplot(pros_all, aes(x = hr_median, fill = factor(preicu_or))) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Distribution of hr_median Stratified by preicu_or",
    x = "Median Heart Rate",
    y = "Density",
    fill = "Pre ICU OR"
  ) +
  scale_fill_manual(values = c("blue", "red"), labels = c("Pre ICU OR = 0", "Pre ICU OR = 1")) +
  theme_minimal()
 # above shows a chart with the non-OR patients having a higher proportion of large median heart rate values

# Now will evaluate distribution of HR median values <= 2 hours after PICU admission and those after.
pros_qc <- pros_all
pros_qc <- pros_qc %>% mutate(t_after_picu_adm = as.numeric(difftime(time1 = score_time, time2 = picu_adm_date_time, units = "hours")))
pros_qc <- pros_qc %>% mutate(above_2 = ifelse(t_after_picu_adm > 2, yes = 1, no = 0))

# Create the plot
ggplot(pros_qc, aes(x = hr_median, fill = factor(above_2))) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Distribution of hr_median Stratified by above_2",
    x = "Median Heart Rate",
    y = "Density",
    fill = "Above 2 Hours"
  ) +
  scale_fill_manual(values = c("blue", "green"), labels = c("<= 2 hours", "> 2 hours")) +
  theme_minimal() # distribution is basically the same though a slightly lower high HR median peak for the <=2 hour subgroup

# Compare distribution of hr_median values between retro and pros datasets when only one value per patient (around 2 hours after
# PICU admit) is retained
pros_one_row_per_pt <- pros_qc %>% filter(t_after_picu_adm > 2) %>% filter(t_after_picu_adm < 3.5)
pros_one_row_per_pt <- pros_one_row_per_pt %>%group_by(study_id) %>% filter (t_after_picu_adm == min(t_after_picu_adm)) %>% ungroup()
pros_one_row_per_pt <- pros_one_row_per_pt %>% dplyr::select(study_id, pat_enc_csn_id, hsp_account_id, age, hr_median, score_time) %>% distinct()

# First load retrospective dataset
retro_df <- read.fst(path = "~/sbi_blake/jan_25_23_model_data_df.fst")
retro_df <- retro_df %>% rename(age = age_years) # fix age for later comparison
retro_df <- retro_df %>% rename(hr_median = median_hr)



# Create the plot to compare hr median values between the two
# Add a column to differentiate the datasets
retro_df <- retro_df %>% mutate(dataset = "Retrospective")
pros_one_row_per_pt <- pros_one_row_per_pt %>% mutate(dataset = "Prospective")

# Combine the two datasets
combined_df <- bind_rows(retro_df, pros_one_row_per_pt)

# Create the plot
ggplot(combined_df, aes(x = hr_median, fill = dataset)) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Distribution of hr_median: Retrospective vs Prospective",
    x = "Median Heart Rate",
    y = "Density",
    fill = "Dataset"
  ) +
  scale_fill_manual(values = c("blue", "green")) +
  theme_minimal()

# Create the plot of HR median distributions for the whole prospective datadset and retrospective dataset
pros_hr_all <- pros_all %>% mutate(dataset = "Prospective")
pros_hr_all <- pros_hr_all %>% dplyr::select(hr_median, dataset)
retro_df <- retro_df %>% mutate(dataset = "Retrospective")
combined_hr <- bind_rows(retro_df, pros_hr_all)

ggplot(combined_hr, aes(x = hr_median, fill = dataset)) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Distribution of hr_median: Retrospective vs Prospective",
    x = "Median Heart Rate",
    y = "Density",
    fill = "Dataset"
  ) +
  scale_fill_manual(values = c("blue", "green")) +
  theme_minimal() +
  theme(
    text              = element_text(size = 14, face = "bold"),  # Base text
    axis.title        = element_text(size = 14, face = "bold"),  # Axis titles
    axis.text         = element_text(size = 12, face = "bold"),  # Axis tick labels
    legend.title      = element_text(size = 14, face = "bold"),  # Legend title
    legend.text       = element_text(size = 12, face = "bold"),  # Legend labels
    plot.title        = element_text(size = 16, face = "bold")   # Plot title
  )
# Ok so apparently that HR issue was not really a huge issue. There's still a larger spike in the retrospective HR values at higher
# percentiles, but the general shape is the same. Will need to look into the other variables to see if distributions similar...




# Now look at the age distribution between the retro and prospective groups
# Add a column to differentiate the datasets
pros_qc <- pros_qc %>% mutate(dataset = "Prospective") %>% dplyr::select(age, dataset)
retro_df <- retro_df %>% mutate(dataset = "Retrospective") %>% dplyr::select(age, dataset)

# Combine the two datasets
combined_df <- bind_rows(pros_qc, retro_df)

# Create the plot
ggplot(combined_df, aes(x = age, fill = dataset)) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Overlay of Age Distributions: Prospective vs Retrospective",
    x = "Age (years)",
    y = "Density",
    fill = "Dataset"
  ) +
  scale_fill_manual(values = c("blue", "green")) +
  theme_minimal() +
  theme(
    text              = element_text(size = 14, face = "bold"),  # Base text
    axis.title        = element_text(size = 14, face = "bold"),  # Axis titles
    axis.text         = element_text(size = 12, face = "bold"),  # Axis tick labels
    legend.title      = element_text(size = 14, face = "bold"),  # Legend title
    legend.text       = element_text(size = 12, face = "bold"),  # Legend labels
    plot.title        = element_text(size = 16, face = "bold")   # Plot title
  )
# There are generally similar shapes with a peak around 1yo, however this peak is certainly larger in the prospective dataset





#### 10-17-24, hr_median is fine, need to check other variables... ####
# Write code to compare distributions of retrospective data to two different sets of prospective data:
# 1. All predictions, 2. Only those predictions around the 2-hour post PICU admission mark

retro_df <- read.fst(path = "~/sbi_blake/jan_25_23_model_data_df.fst")
retro_df <- retro_df %>% rename(age = age_years) # fix age for later comparison
retro_df <- retro_df %>% rename(hr_median = median_hr)

# Start by transforming column names from retro to match prospective data
# List of vital sign suffixes to search for
suffixes <- c("_dbp", "_fio2", "_rr", "_sat", "_sbp", "_hr", "_temp")

# Function to rename columns based on specified suffixes
rename_vital_sign_columns <- function(col_name) {
  # Loop through each suffix and check if the column name ends with the suffix
  for (suffix in suffixes) {
    if (str_ends(col_name, suffix)) {
      # Extract the vital sign abbreviation (e.g., hr, dbp)
      vital_sign <- str_extract(col_name, suffix)
      # Remove the suffix from the original column name
      remaining_name <- str_remove(col_name, suffix)
      # Create the new column name by moving the vital sign abbreviation to the front
      return(paste0(str_remove(vital_sign, "_"), "_", remaining_name))
    }
  }
  # If no suffix match, return the original column name
  return(col_name)
}

# Apply the renaming function to all column names in the dataframe
colnames(retro_df) <- sapply(colnames(retro_df), rename_vital_sign_columns)

## Now fix the sat_ columns to be o2sat_
# Function to rename columns starting with "sat_" to "o2sat_"
rename_sat_columns <- function(col_name) {
  if (str_starts(col_name, "sat_")) {
    # Replace "sat_" with "o2sat_"
    return(str_replace(col_name, "^sat_", "o2sat_"))
  }
  # Return the original name if it doesn't start with "sat_"
  return(col_name)
}

# Apply the renaming function to all column names in the dataframe
colnames(retro_df) <- sapply(colnames(retro_df), rename_sat_columns)

## Replace _number with _count in column names
# Function to replace "_number" with "_count" in column names
replace_number_with_count <- function(col_name) {
  # Replace occurrences of "_number" with "_count"
  return(str_replace(col_name, "_number", "_count"))
}

# Apply the function to all column names in the dataframe
colnames(retro_df) <- sapply(colnames(retro_df), replace_number_with_count)


# Fix other one-off columns
retro_df <- retro_df %>% rename(is_female = sex,
                                albumin_mean = albumin_blood,
                                ast_present = ast_blood_pres,
                                bacteria_urine_mean = bacteria_urine,
                                bands_perc_mean = bands_perc_blood,
                                bands_perc_present = bands_perc_blood_pres,
                                base_excess_present = base_excess_blood_pres,
                                bun_mean = bun_blood,
                                chloride_mean = chloride_blood,
                                crp_mean = crp_blood,
                                fibrinogen_mean = fibrinogen,
                                gcs_total_min = min_tot_gcs,
                                gcs_verbal_min = min_verbal,
                                hematocrit_mean = hematocrit_blood,
                                hemoglobin_present = hemoglobin_blood_pres,
                                lactate_max = lactate_blood,
                                lactate_present = lactate_blood_pres,
                                ldh_total_mean = ldh_total_serum,
                                leukocytes_urine_mean = leukocytes_urine,
                                leukocytes_urine_present = leukocytes_urine_pres,
                                lipase_mean = lipase_blood,
                                malignancy_pccc = cancer,
                                mcv_mean = mcv_blood,
                                monocytes_perc_mean = monocytes_perc_blood,
                                nitrate_urine_mean = nitrite_urine,
                                origin = pre_icu,
                                pccc = ccc,
                                pco2_mean = pco2_gas_blood,
                                ph_gas_blood_mean = ph_gas_blood,
                                ph_urine_mean = ph_urine,
                                platelet_count_mean = platelet_count_blood,
                                po2_arterial_present = pao2_gas_blood_pres,
                                po2_venous_mean = po2_venous,
                                po2_venous_present = po2_venous_pres,
                                respiratory_support_any_positive = invasive_o2,
                                salicylates_present = salicylates_serum_pres,
                                scheduled_admit = sched_adm,
                                sodium_present = sodium_blood_pres,
                                uric_acid_present = uric_acid_blood_pres,
                                wbc_csf_present = wbc_csf_pres,
                                wbc_urine_mean = wbcs_urine,
                                weight_present = wt_pres
                                )

# Find the column names that are common between retro_df and pros_all
common_columns <- intersect(colnames(retro_df), colnames(pros_all))

# Create a dataframe comparing the classes of the common columns
class_comparison_df <- data.frame(
  col_name = common_columns,                     # Common column names
  retro_class = sapply(retro_df[common_columns], class),  # Class in retro_df
  pros_class = sapply(pros_all[common_columns], class)    # Class in pros_all
)

# Below line can be used to select which lines / predictors to plot comparisons of
cols_to_comp <- class_comparison_df$col_name[1:nrow(class_comparison_df)]



# Loop through each column and generate the plots
for (col in cols_to_comp) {
  # Create two separate data frames for retro_df and pros_all for the current column
  retro_data <- data.frame(value = retro_df[[col]], dataset = "retro_df")
  pros_data  <- data.frame(value = pros_all[[col]], dataset = "pros_all")

  # Combine both data frames
  plot_data <- rbind(retro_data, pros_data)

  # If there's nothing in plot_data, just break
  if (is.null(plot_data)) {
    break
  }

  # Check if the current column is numeric or categorical
  if (is.numeric(plot_data$value)) {
    # --- Numeric Data: Use density plot ---
    p <- ggplot(plot_data, aes(x = value, fill = dataset)) +
      geom_density(alpha = 0.5) +
      labs(
        title = paste("Distribution of", col),
        x     = col,
        y     = "Density"
      ) +
      theme_minimal() +
      theme(
        text              = element_text(size = 14, face = "bold"),  # Base text
        axis.title        = element_text(size = 14, face = "bold"),  # Axis titles
        axis.text         = element_text(size = 12, face = "bold"),  # Axis tick labels
        legend.title      = element_text(size = 14, face = "bold"),  # Legend title
        legend.text       = element_text(size = 12, face = "bold"),  # Legend labels
        plot.title        = element_text(size = 16, face = "bold")   # Plot title
      )
  } else {
    # --- Categorical/Binary Data: Use grouped bar chart of % ---
    # 1. Summarize data: count how many rows fall into each category per dataset
    plot_data_summarized <- plot_data %>%
      group_by(dataset, value) %>%
      summarise(n = n(), .groups = "drop") %>%
      # 2. Calculate percent within each dataset
      group_by(dataset) %>%
      mutate(percent = (n / sum(n)) * 100) %>%
      ungroup()

    # 3. Plot as a bar chart with percent on the y-axis
    p <- ggplot(plot_data_summarized, aes(x = value, y = percent, fill = dataset)) +
      geom_col(position = "dodge") +
      labs(
        title = paste("Distribution of", col),
        x     = col,
        y     = "Percent"
      ) +
      theme_minimal() +
      theme(
        text              = element_text(size = 14, face = "bold"),  # Base text
        axis.title        = element_text(size = 14, face = "bold"),  # Axis titles
        axis.text         = element_text(size = 12, face = "bold"),  # Axis tick labels
        legend.title      = element_text(size = 14, face = "bold"),  # Legend title
        legend.text       = element_text(size = 12, face = "bold"),  # Legend labels
        plot.title        = element_text(size = 16, face = "bold")   # Plot title
      )
  }

  # Print the plot
  print(p)
}


# Now produce charts on the SBi outcomes and percentages:
#add column to retro_df to show if there was a micro diagnosed sbi
retro_orig <- retro_df
retro_df <- retro_df %>% mutate(micro_sbi_1_0 = as.numeric(retro_df$sbi) - 1)

# Rename the retro_df sbi_outcomes
retro_df <- retro_df %>% rename(pna_1_0 = sbi_pneumonia, cx_neg_sepsis = sbi_cx_neg_sepsis)


# 1) Combine the relevant columns from retro_df and pros_all into one long data frame
df_retro <- data.frame(
  dataset        = "Retrospective",
  pna_1_0        = retro_df$pna_1_0,
  cx_neg_sepsis  = retro_df$cx_neg_sepsis,
  micro_sbi_1_0  = retro_df$micro_sbi_1_0
)

df_pros <- data.frame(
  dataset        = "Prospective",
  pna_1_0        = pros_all$pna_1_0,
  cx_neg_sepsis  = pros_all$cx_neg_sepsis,
  micro_sbi_1_0  = pros_all$micro_sbi_1_0
)

# Combine them row-wise
df_combined <- rbind(df_retro, df_pros)

# Convert from "wide" to "long" format, creating a single 'outcome' column
df_long <- df_combined %>%
  pivot_longer(
    cols      = c("pna_1_0", "cx_neg_sepsis", "micro_sbi_1_0"),
    names_to  = "outcome",
    values_to = "value"
  )

# 2) Calculate the percent of rows with value == 1 for each (dataset, outcome) pair
df_summary <- df_long %>%
  group_by(dataset, outcome) %>%
  summarize(
    count_ones = sum(value == 1, na.rm = TRUE),
    n_total    = n(),
    percent_1  = (count_ones / n_total) * 100,
    .groups    = "drop"
  )

# 3) Plot a grouped bar chart of these percentages
ggplot(df_summary, aes(x = outcome, y = percent_1, fill = dataset)) +
  geom_col(position = position_dodge()) +
  labs(
    title = "Percent of Patients WIth SBI Subtypes",
    x     = "SBI Type",
    y     = "Percent of Patients"
  ) +
  theme_minimal() +
  theme(
    text              = element_text(size = 14, face = "bold"),  # Base text
    axis.title        = element_text(size = 14, face = "bold"),  # Axis titles
    axis.text         = element_text(size = 12, face = "bold"),  # Axis tick labels
    legend.title      = element_text(size = 14, face = "bold"),  # Legend title
    legend.text       = element_text(size = 12, face = "bold"),  # Legend labels
    plot.title        = element_text(size = 16, face = "bold")   # Plot title
  )



# Now need to look at calibration across different SBI probability deciles
# First change probability to decimal form
pros_data <- pros_all %>% mutate(model_score = pros_all$model_score / 100)

# Create dfs for each model and abx exposure group
pros_data_lr_no_abx <- pros_data %>% filter(model_type == "LR_no_abx")
pros_data_lr_yes_abx <- pros_data %>% filter(model_type == "LR_yes_abx")
pros_data_rf_no_abx <- pros_data %>% filter(model_type == "RF_no_abx")
pros_data_rf_yes_abx <- pros_data %>% filter(model_type == "RF_yes_abx")


# 1) Bin the predicted probabilities into deciles (ntile = 10)
cal_data_lr_no <- pros_data_lr_no_abx %>%
  mutate(prob_bin = ntile(model_score, 10)) %>%
  group_by(prob_bin) %>%
  summarize(
    mean_pred = mean(model_score, na.rm = TRUE),    # Average predicted probability in this bin
    obs_rate  = mean(sbi_present, na.rm = TRUE),    # Observed proportion of SBI=1 in this bin
    .groups   = "drop"
  )

# 2) Plot the calibration curve
ggplot(cal_data_lr_no, aes(x = mean_pred, y = obs_rate)) +
  geom_point(size = 3) +         # Plot points for each bin
  geom_line(size = 1) +          # Connect the points
  geom_abline(
    slope = 1, intercept = 0,
    linetype = "dashed", color = "red"
  ) +
  labs(
    title = "Calibration Plot for LR Model (Abx Unexposed)",
    x = "Mean Predicted Probability (0–1)",
    y = "Observed SBI Rate (0–1)"
  ) +
  theme_minimal() +
  theme(
    text              = element_text(size = 14, face = "bold"),  # Base text
    axis.title        = element_text(size = 14, face = "bold"),  # Axis titles
    axis.text         = element_text(size = 12, face = "bold"),  # Axis tick labels
    legend.title      = element_text(size = 14, face = "bold"),  # Legend title
    legend.text       = element_text(size = 12, face = "bold"),  # Legend labels
    plot.title        = element_text(size = 16, face = "bold")   # Plot title
  )



# 1) Bin the predicted probabilities into deciles (ntile = 10)
cal_data_lr_yes <- pros_data_lr_yes_abx %>%
  mutate(prob_bin = ntile(model_score, 10)) %>%
  group_by(prob_bin) %>%
  summarize(
    mean_pred = mean(model_score, na.rm = TRUE),    # Average predicted probability in this bin
    obs_rate  = mean(sbi_present, na.rm = TRUE),    # Observed proportion of SBI=1 in this bin
    .groups   = "drop"
  )

# 2) Plot the calibration curve
ggplot(cal_data_lr_yes, aes(x = mean_pred, y = obs_rate)) +
  geom_point(size = 3) +         # Plot points for each bin
  geom_line(size = 1) +          # Connect the points
  geom_abline(
    slope = 1, intercept = 0,
    linetype = "dashed", color = "red"
  ) +                            # Ideal calibration line
  labs(
    title = "Calibration Plot for LR Model (Abx Exposed)",
    x = "Mean Predicted Probability (0–100)",
    y = "Observed SBI Rate (0–1)"
  ) +
  theme_minimal() +
  theme(
    text              = element_text(size = 14, face = "bold"),  # Base text
    axis.title        = element_text(size = 14, face = "bold"),  # Axis titles
    axis.text         = element_text(size = 12, face = "bold"),  # Axis tick labels
    legend.title      = element_text(size = 14, face = "bold"),  # Legend title
    legend.text       = element_text(size = 12, face = "bold"),  # Legend labels
    plot.title        = element_text(size = 16, face = "bold")   # Plot title
  )







# 1) Bin the predicted probabilities into deciles (ntile = 10)
cal_data_rf_no <- pros_data_rf_no_abx %>%
  mutate(prob_bin = ntile(model_score, 10)) %>%
  group_by(prob_bin) %>%
  summarize(
    mean_pred = mean(model_score, na.rm = TRUE),    # Average predicted probability in this bin
    obs_rate  = mean(sbi_present, na.rm = TRUE),    # Observed proportion of SBI=1 in this bin
    .groups   = "drop"
  )

# 2) Plot the calibration curve
ggplot(cal_data_rf_no, aes(x = mean_pred, y = obs_rate)) +
  geom_point(size = 3) +         # Plot points for each bin
  geom_line(size = 1) +          # Connect the points
  geom_abline(
    slope = 1, intercept = 0,
    linetype = "dashed", color = "red"
  ) +                            # Ideal calibration line
  labs(
    title = "Calibration Plot for RF Model (Abx Unexposed)",
    x = "Mean Predicted Probability (0–100)",
    y = "Observed SBI Rate (0–1)"
  ) +
  theme_minimal() +
  theme(
    text              = element_text(size = 14, face = "bold"),  # Base text
    axis.title        = element_text(size = 14, face = "bold"),  # Axis titles
    axis.text         = element_text(size = 12, face = "bold"),  # Axis tick labels
    legend.title      = element_text(size = 14, face = "bold"),  # Legend title
    legend.text       = element_text(size = 12, face = "bold"),  # Legend labels
    plot.title        = element_text(size = 16, face = "bold")   # Plot title
  )



# 1) Bin the predicted probabilities into deciles (ntile = 10)
cal_data_rf_yes <- pros_data_rf_yes_abx %>%
  mutate(prob_bin = ntile(model_score, 10)) %>%
  group_by(prob_bin) %>%
  summarize(
    mean_pred = mean(model_score, na.rm = TRUE),    # Average predicted probability in this bin
    obs_rate  = mean(sbi_present, na.rm = TRUE),    # Observed proportion of SBI=1 in this bin
    .groups   = "drop"
  )

# 2) Plot the calibration curve
ggplot(cal_data_rf_yes, aes(x = mean_pred, y = obs_rate)) +
  geom_point(size = 3) +         # Plot points for each bin
  geom_line(size = 1) +          # Connect the points
  geom_abline(
    slope = 1, intercept = 0,
    linetype = "dashed", color = "red"
  ) +
  scale_x_continuous(limits = c(0, 1)) +   # Force x-axis from 0–1
  scale_y_continuous(limits = c(0, 1)) +   # Force y-axis from 0–1
  labs(
    title = "Calibration Plot for RF Model (Abx Exposed)",
    x = "Mean Predicted Probability (0–100)",
    y = "Observed SBI Rate (0–1)"
  ) +
  theme_minimal() +
  theme(
    text              = element_text(size = 14, face = "bold"),  # Base text
    axis.title        = element_text(size = 14, face = "bold"),  # Axis titles
    axis.text         = element_text(size = 12, face = "bold"),  # Axis tick labels
    legend.title      = element_text(size = 14, face = "bold"),  # Legend title
    legend.text       = element_text(size = 12, face = "bold"),  # Legend labels
    plot.title        = element_text(size = 16, face = "bold")   # Plot title
  )






cal_data <- pros_data_lr_yes_abx %>%
  # 1) Use cut() to bin scores by 0.1 increments between 0 and 1
  mutate(
    prob_bin = cut(
      model_score,
      breaks = seq(0, 1, by = 0.1),
      include.lowest = TRUE,   # Ensures 0 is included in first bin
      right = FALSE            # If right=FALSE, intervals are [0,0.1), [0.1,0.2), etc.
    )
  ) %>%
  group_by(prob_bin) %>%
  summarize(
    mean_pred = mean(model_score, na.rm = TRUE),    # Average predicted probability in this bin
    obs_rate  = mean(sbi_present, na.rm = TRUE),    # Observed proportion of SBI=1 in this bin
    .groups   = "drop"
  )

# 2) Plot the calibration curve
ggplot(cal_data, aes(x = mean_pred, y = obs_rate)) +
  geom_point(size = 3) +         # Plot points for each bin
  geom_line(size = 1) +          # Connect the points
  geom_abline(
    slope = 1, intercept = 0,
    linetype = "dashed", color = "red"
  ) +
  labs(
    title = "LR Model Yes Abx Exposure",
    x = "Mean Predicted Probability (0–1)",
    y = "Observed SBI Rate (0–1)"
  ) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() +
  theme(
    text              = element_text(size = 14, face = "bold"),
    axis.title        = element_text(size = 14, face = "bold"),
    axis.text         = element_text(size = 12, face = "bold"),
    legend.title      = element_text(size = 14, face = "bold"),
    legend.text       = element_text(size = 12, face = "bold"),
    plot.title        = element_text(size = 16, face = "bold")
  )


