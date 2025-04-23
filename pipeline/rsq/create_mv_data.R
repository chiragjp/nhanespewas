## Chirag Patel
## create a large dataset to estimate R2
## 05/15/23
## updated  4/13/24

## nutrition and biomarkers of exposure that are measured in ~10k
## usage: Rscript create_mv_data.R -p PHENOTYPE_NAME
## requires a summary_stats and a nhanes individual level database

library(getopt)
library(tidyverse)
library(DBI)
library(devtools)
load_all('..')

spec <- matrix(c(
  'pvarname', 'p', 1, "character"
), byrow=TRUE, ncol=4);
opt <- getopt(spec)

pvarname_to_query <- opt$pvarname

path_to_nhanes <- '../../db/nhanes_031725.sqlite'
path_to_summary_stats <- '../../db/pe_summary_stats_01_2025.sqlite'


nhanes_con <- DBI::dbConnect(RSQLite::SQLite(), dbname=path_to_nhanes)
summary_stats_con <- DBI::dbConnect(RSQLite::SQLite(), dbname=path_to_summary_stats)

summ_stats <- tbl(summary_stats_con, "adjusted_meta") |> filter(pvarname == pvarname_to_query,
                                                                p.value < 5e-4, nobs > 2, # FDR significant
                                                                model_number == 2) |> filter(term_name == 'expo' | term_name == 'expo1') |> select(pvarname, evarname, p.value, nobs)

rsq_1 <- tbl(summary_stats_con, "rsq_overall") |> filter(model_number == 2, aggregate_base_model ==0, phenotype==pvarname_to_query) |> select(exposure, rsq)
rsq_2 <- tbl(summary_stats_con, "rsq_overall") |> filter(model_number == 2, aggregate_base_model ==1, phenotype==pvarname_to_query) |> select(exposure, rsq) |> rename(rsq_base=rsq)

rsq <- rsq_1 |>left_join(rsq_2, by="exposure") |> mutate(rsq_adj = rsq-rsq_base)
summ_stats <- summ_stats |> left_join(rsq, by=c("evarname"="exposure")) |> collect()

DBI::dbDisconnect(summary_stats_con); remove(summary_stats_con)

demo_table_for_begin_year <- function(nhanes_con, yr) {
  demo_name <- "DEMO"
  # Map the beginning year to the corresponding DEMO table name
  if(yr == 1999) {
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
    stop("Error: The year provided is not valid. Use one of 2001, 2003, 2005, 2007, 2009, 2011, 2013, 2015, or 2017.")
  }

  # Return the table from the connection corresponding to the determined demo_name
  return(dplyr::tbl(nhanes_con, demo_name))
}



get_phenotype_across_surveys <- function(nhanes_con, pvarname) {
  pvarname_to_query <- pvarname
  variable_names_epcf <- tbl(nhanes_con, "variable_names_epcf")
  dfiles <- variable_names_epcf |> filter(Variable.Name == pvarname_to_query) |> select(Variable.Name, Data.File.Name, Begin.Year) |> collect()
  n_files <- dfiles |> count() |> pull()
  pheno_data <- vector("list", length = nrow(dfiles))
  for(i in 1:n_files) {
    tble_name <- dfiles |> slice(i) |> pull(Data.File.Name)
    yr <- dfiles |> slice(i) |> pull(Begin.Year)
    vname <- dfiles |> slice(i) |> pull(Variable.Name)
    demo <- demo_table_for_begin_year(nhanes_con, yr)
    pheno_data[[i]] <- demo |> left_join(tbl(nhanes_con, tble_name) |>
                                           select(SEQN,all_of(vname)), by="SEQN") |> collect()
  }

  pheno_data <- pheno_data |> bind_rows()

}


get_exposure_across_surveys <- function(nhanes_con, evarname) {
  etables <- get_table_names_for_varname(nhanes_con, varname = evarname)
  etable_list <- vector("list", length=nrow(etables))
  for(ii in 1:nrow(etables)) {
    dname <- etables |> slice(ii) |> pull(Data.File.Name)
    etable_list[[ii]] <- tbl(nhanes_con, dname) |> select("SEQN", all_of(evarname)) |> collect()
  }

  etable_list <- bind_rows(etable_list)

}

get_exposures_across_surveys <- function(con, evar_table, pheno_data) {
  # evar_table is just a table that i a list of evarname
  for(ii in 1:nrow(summ_stats)) {
    evar <- summ_stats |> slice(ii) |> pull(evarname)
    e_to_merge <- get_exposure_across_surveys(con, evar)
    pheno_data <- pheno_data |> left_join(e_to_merge, by="SEQN")
  }
  pheno_data
}


ptable <- get_phenotype_across_surveys(nhanes_con, pvarname_to_query)
p_mv <- get_exposures_across_surveys(nhanes_con, summ_stats |> select(evarname), pheno_data=ptable)

DBI::dbDisconnect(nhanes_con)

cols_to_keep <- c("SEQN", "SDMVPSU", "WTMEC2YR", "INDFMPIR", "SDDSRVYR","RIAGENDR","RIDAGEYR",
                  "ETHNICITY_NONHISPANICWHITE", "ETHNICITY_MEXICAN", "ETHNICITY_OTHERHISPANIC", "ETHNICITY_OTHER", "ETHNICITY_NONHISPANICBLACK", "AGE_SQUARED",
                  "EDUCATION_LESS9",  "EDUCATION_9_11", "EDUCATION_HSGRAD", "EDUCATION_AA", "EDUCATION_COLLEGEGRAD", "BORN_INUSA", "URXUCR_adj",
                  "DR1TKCAL_adj", "DR2TKCAL_adj", "BMXBMI_adj", "BMXHT_adj")

cols_to_keep <- c(cols_to_keep, summ_stats |> pull(evarname), pvarname_to_query)
p_mv <- p_mv |> select(all_of(cols_to_keep))
file_out <- sprintf('%s_mv.rds', pvarname_to_query)

to_out <- list(big_data=p_mv, summ_stats=summ_stats, phenotype=pvarname_to_query)
write_rds(to_out, file_out)

