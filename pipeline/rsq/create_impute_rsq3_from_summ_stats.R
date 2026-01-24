#!/usr/bin/env Rscript
## Chirag Patel
## create a large dataset to estimate R2, impute, and estimate the R2
## combination of create_mv_data and impute_and_r2.R
## updated  4/13/24
##
## usage:
##   Rscript create_impute_rsq.R -p PHENOTYPE_NAME -s /path/to/adjusted_meta_2.csv
##
## requires:
##   - a meta CSV (formerly adjusted_meta_2.Rdata) with required columns (see below)
##   - a NHANES SQLite database at path_to_nhanes

suppressPackageStartupMessages({
  library(getopt)
  library(tidyverse)
  library(DBI)
  library(devtools)
  library(mice)
  library(finalfit)
  library(splines)
})

load_all("..")

MISSING_PCT_THRESHOLD <- 40 # default (not currently used downstream)
ADDTC <- FALSE

spec <- matrix(c(
  "pvarname",     "p", 1, "character", # CFDDS
  "missing_pct",  "m", 2, "double", # 40
  "meta_rds",     "s", 1, "character",
  "nhanes_sqlite", "d", 1, "character"
), byrow = TRUE, ncol = 4)

opt <- getopt(spec)

pvarname_to_query <- opt$pvarname
meta_rds_path <- opt$meta_rds

# Allow override of DB path; otherwise default to prior hard-coded location
path_to_nhanes <- if (!is.null(opt$nhanes_sqlite) && opt$nhanes_sqlite != "") {
  opt$nhanes_sqlite
} else {
  "../../db/nhanes_031725.sqlite"
}

if (is.null(pvarname_to_query) || pvarname_to_query == "") {
  stop("You must provide -p / --pvarname (phenotype variable name), e.g., -p BMXHT")
}
if (is.null(meta_rds_path) || meta_rds_path == "") {
  stop("You must provide -s / --meta_rds pointing to the adjusted_meta_2 rds file.")
}
if (!file.exists(meta_rds_path)) {
  stop("meta_csv file not found: ", meta_rds_path)
}
if (!file.exists(path_to_nhanes)) {
  stop("NHANES sqlite file not found: ", path_to_nhanes)
}

# Read meta/summary stats from CSV (replaces load(.Rdata))
adjusted_meta_2 <- readr::read_rds(meta_rds_path)

# Validate required columns exist
required_cols <- c(
  "pvarname", "evarname", "p.value", "nobs",
  "rsq_adjusted_base_diff", "model_number", "term_name"
)
missing_cols <- setdiff(required_cols, colnames(adjusted_meta_2))
if (length(missing_cols) > 0) {
  stop("meta_rds is missing required columns: ", paste(missing_cols, collapse = ", "))
}

summ_stats <- adjusted_meta_2 |>
  filter(
    pvarname == pvarname_to_query,
    p.value < .05,
    nobs >= 2,
    model_number == 2
  ) |>
  filter(term_name == "expo" | term_name == "expo1") |>
  select(pvarname, evarname, p.value, nobs, rsq_adjusted_base_diff) |>
  rename(rsq_adj = rsq_adjusted_base_diff)

cat("assembling the data\n")
nhanes_con <- DBI::dbConnect(RSQLite::SQLite(), dbname = path_to_nhanes)

demo_table_for_begin_year <- function(nhanes_con, yr) {
  demo_name <- "DEMO"
  if (yr == 1999) {
    demo_name <- "DEMO"
  } else if (yr == 2001) {
    demo_name <- "DEMO_B"
  } else if (yr == 2003) {
    demo_name <- "DEMO_C"
  } else if (yr == 2005) {
    demo_name <- "DEMO_D"
  } else if (yr == 2007) {
    demo_name <- "DEMO_E"
  } else if (yr == 2009) {
    demo_name <- "DEMO_F"
  } else if (yr == 2011) {
    demo_name <- "DEMO_G"
  } else if (yr == 2013) {
    demo_name <- "DEMO_H"
  } else if (yr == 2015) {
    demo_name <- "DEMO_I"
  } else if (yr == 2017) {
    demo_name <- "DEMO_J"
  } else {
    stop("Error: The year provided is not valid. Use one of 1999, 2001, 2003, 2005, 2007, 2009, 2011, 2013, 2015, or 2017.")
  }
  dplyr::tbl(nhanes_con, demo_name)
}

get_phenotype_across_surveys <- function(nhanes_con, pvarname) {
  variable_names_epcf <- tbl(nhanes_con, "variable_names_epcf")
  dfiles <- variable_names_epcf |>
    filter(Variable.Name == pvarname) |>
    select(Variable.Name, Data.File.Name, Begin.Year) |>
    collect()
  
  if (nrow(dfiles) == 0) {
    stop("Phenotype variable not found in variable_names_epcf: ", pvarname)
  }
  
  pheno_data <- vector("list", length = nrow(dfiles))
  for (i in seq_len(nrow(dfiles))) {
    tble_name <- dfiles |> slice(i) |> pull(Data.File.Name)
    yr <- dfiles |> slice(i) |> pull(Begin.Year)
    vname <- dfiles |> slice(i) |> pull(Variable.Name)
    
    demo <- demo_table_for_begin_year(nhanes_con, yr)
    
    pheno_data[[i]] <- demo |>
      left_join(
        tbl(nhanes_con, tble_name) |>
          select(SEQN, all_of(vname)),
        by = "SEQN"
      ) |>
      collect()
  }
  
  bind_rows(pheno_data)
}

get_exposure_across_surveys <- function(nhanes_con, evarname) {
  etables <- get_table_names_for_varname(nhanes_con, varname = evarname)
  
  if (nrow(etables) == 0) {
    warning("No tables found for exposure variable: ", evarname)
    return(tibble(SEQN = integer(), !!evarname := numeric()))
  }
  
  etable_list <- vector("list", length = nrow(etables))
  for (ii in seq_len(nrow(etables))) {
    dname <- etables |> slice(ii) |> pull(Data.File.Name)
    etable_list[[ii]] <- tbl(nhanes_con, dname) |>
      select(SEQN, all_of(evarname)) |>
      collect()
  }
  bind_rows(etable_list)
}

get_exposures_across_surveys <- function(con, evar_table, pheno_data) {
  for (ii in seq_len(nrow(evar_table))) {
    evar <- evar_table |> slice(ii) |> pull(evarname)
    e_to_merge <- get_exposure_across_surveys(con, evar)
    pheno_data <- pheno_data |> left_join(e_to_merge, by = "SEQN")
  }
  pheno_data
}

ptable <- get_phenotype_across_surveys(nhanes_con, pvarname_to_query)

# Exposures to pull (top determined later; here we pull all candidates in summ_stats)
summ_stats_variables <- summ_stats |> select(evarname)
if (ADDTC) {
  tc_to_add <- tibble(evarname = "LBXTC")
  summ_stats_variables <- summ_stats_variables |> bind_rows(tc_to_add)
}

p_mv <- get_exposures_across_surveys(nhanes_con, summ_stats_variables, pheno_data = ptable)

DBI::dbDisconnect(nhanes_con)

cols_to_keep <- c(
  "SEQN", "SDMVPSU", "WTMEC2YR", "SDMVSTRA", "INDFMPIR", "SDDSRVYR", "RIAGENDR", "RIDAGEYR",
  "ETHNICITY_NONHISPANICWHITE", "ETHNICITY_MEXICAN", "ETHNICITY_OTHERHISPANIC", "ETHNICITY_OTHER",
  "ETHNICITY_NONHISPANICBLACK", "AGE_SQUARED",
  "EDUCATION_LESS9", "EDUCATION_9_11", "EDUCATION_HSGRAD", "EDUCATION_AA", "EDUCATION_COLLEGEGRAD",
  "BORN_INUSA", "URXUCR_adj", "DRXTKCAL_adj", "DR1TKCAL_adj", "DR2TKCAL_adj",
  "BMXBMI_adj", "BMXHT_adj"
)

cols_to_keep <- intersect(colnames(p_mv), cols_to_keep)
cols_to_keep <- c(cols_to_keep, summ_stats_variables |> pull(evarname), pvarname_to_query)

p_mv <- p_mv |> select(all_of(cols_to_keep))

big_data <- list(big_data = p_mv, summ_stats = summ_stats, phenotype = pvarname_to_query)

################# impute and R2
cat("now imputing and estimating R2\n")
big_missing_gl <- missing_glimpse(big_data$big_data)

bd <- big_data$big_data
summ_stat <- big_data$summ_stats
phenotype <- big_data$phenotype

## select n=20 highest by R2
M_TOP <- 20
exposures_selected <- summ_stat |> arrange(desc(rsq_adj)) |> slice_head(n = M_TOP)

if (nrow(exposures_selected) == 0) {
  message("no exposures above the threshold of missingness to impute, quitting")
  quit("no", 0)
}

cols <- c(
  "SDMVPSU", "WTMEC2YR", "INDFMPIR", "SDDSRVYR", "SDMVSTRA",
  "ETHNICITY_NONHISPANICWHITE", "ETHNICITY_MEXICAN", "ETHNICITY_OTHERHISPANIC", "ETHNICITY_OTHER",
  "ETHNICITY_NONHISPANICBLACK", "AGE_SQUARED", "RIAGENDR", "RIDAGEYR",
  "EDUCATION_LESS9", "EDUCATION_9_11", "EDUCATION_HSGRAD", "EDUCATION_AA", "EDUCATION_COLLEGEGRAD",
  "BORN_INUSA"
)

if (ADDTC) {
  cols <- c(cols, "LBXTC")
}

bd_exposure_tibble <- bd |>
  select(all_of(c("SEQN", phenotype, cols, exposures_selected$evarname))) |>
  rename(pheno = all_of(phenotype)) |>
  filter(!is.na(pheno)) |>
  filter(!is.na(INDFMPIR))

bd_m <- mice(bd_exposure_tibble, m = 10, print = FALSE)

base_formula <- "I(scale(pheno)) ~ RIAGENDR + RIDAGEYR + AGE_SQUARED + ETHNICITY_NONHISPANICBLACK + ETHNICITY_OTHER + ETHNICITY_MEXICAN + ETHNICITY_OTHERHISPANIC + EDUCATION_LESS9 + EDUCATION_9_11 + EDUCATION_HSGRAD"
if (ADDTC) {
  base_formula <- sprintf("%s + %s", base_formula, "ns(log(LBXTC), df = 3)")
}

extended_formula <- sprintf("%s + %s", base_formula, paste(exposures_selected$evarname, collapse = "+"))

fit1 <- with(bd_m, lm(formula(base_formula)))
fit2 <- with(bd_m, lm(formula(extended_formula)))

r2_1 <- pool.r.squared(fit1)
r2_2 <- pool.r.squared(fit2)

pooled_results <- pool(fit2)

# Multivariate R2
mv_diff <- r2_2[1] - r2_1[1]
print(mv_diff)

uni_r2 <- exposures_selected |> summarize(s = sum(rsq_adj, na.rm = TRUE)) |> pull(s)

to_save <- list(
  rsq_summary = tibble(
    phenotype = phenotype,
    rsq_exposures = mv_diff,
    rsq_base = r2_1[1],
    rsq_uni = uni_r2,
    number_exposures_in_model = nrow(exposures_selected)
  ),
  exposures_selected = exposures_selected,
  fit1 = r2_1,
  fit2 = r2_2,
  pooled_fit = pooled_results
)

dir.create("./out", showWarnings = FALSE, recursive = TRUE)
file_out <- sprintf("./out/%s_imp_R2_cirsq3_from_summ_stats.rds", phenotype)
write_rds(to_save, file = file_out)

cat("saved: ", file_out, "\n")