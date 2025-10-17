## Chirag Patel
## create a large dataset to estimate R2, impute, and estimate the R2
## combination of create_mv_data and impute_and_r2.R
## updated  4/13/24


## usage: Rscript create_impute_rsq.R -p PHENOTYPE_NAME
## requires a summary_stats and a nhanes individual level database

library(getopt)
library(tidyverse)
library(DBI)
library(devtools)
library(mice)
library(finalfit)
library(splines)

load_all('..')

MISSING_PCT_THRESHOLD <- 40 # default
ADDTC <- F
spec <- matrix(c(
  'pvarname', 'p', 1, "character",
  'missing_pct', 'm', 2, "double"
), byrow=TRUE, ncol=4);
opt <- getopt(spec)

pvarname_to_query <- opt$pvarname
#pvarname_to_query <- 'BMXHT'

### first create the data
path_to_nhanes <- '../../db/nhanes_031725.sqlite'

load("../../rmd/adjusted_meta_2.Rdata")

#pvarname_to_query <- "GrimAgeMort"

summ_stats <- adjusted_meta_2 |> filter(pvarname == pvarname_to_query,
       #p.value < 5e-4, nobs >= 2, # FDR significant
       p.value < 5e-2, nobs >= 2, 
       model_number == 2) |> filter(term_name == 'expo' | term_name == 'expo1') |> select(pvarname, evarname, p.value, nobs, rsq_adjusted_base_diff) |> rename(rsq_adj=rsq_adjusted_base_diff)


cat("assembling the data\n")
### first create the data
path_to_nhanes <- '../../db/nhanes_031725.sqlite'
nhanes_con <- DBI::dbConnect(RSQLite::SQLite(), dbname=path_to_nhanes)


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
  for(ii in 1:nrow(evar_table)) {
    evar <- evar_table |> slice(ii) |> pull(evarname)
    e_to_merge <- get_exposure_across_surveys(con, evar)
    pheno_data <- pheno_data |> left_join(e_to_merge, by="SEQN")
  }
  pheno_data
}


ptable <- get_phenotype_across_surveys(nhanes_con, pvarname_to_query)
# add in LBXTC if flag is set
summ_stats_variables <-summ_stats |> select(evarname)
if(ADDTC) {
  tc_to_add <- tibble(evarname="LBXTC")
  summ_stats_variables <- summ_stats |> select(evarname) |> rbind(tc_to_add)
}


p_mv <- get_exposures_across_surveys(nhanes_con,summ_stats_variables , pheno_data=ptable)


DBI::dbDisconnect(nhanes_con)

cols_to_keep <- c("SEQN", "SDMVPSU", "WTMEC2YR", "SDMVSTRA", "INDFMPIR", "SDDSRVYR","RIAGENDR","RIDAGEYR",
                  "ETHNICITY_NONHISPANICWHITE", "ETHNICITY_MEXICAN", "ETHNICITY_OTHERHISPANIC", "ETHNICITY_OTHER", "ETHNICITY_NONHISPANICBLACK", "AGE_SQUARED",
                  "EDUCATION_LESS9",  "EDUCATION_9_11", "EDUCATION_HSGRAD", "EDUCATION_AA", "EDUCATION_COLLEGEGRAD", "BORN_INUSA", "URXUCR_adj", "DRXTKCAL_adj", "DR1TKCAL_adj", "DR2TKCAL_adj",
                   "BMXBMI_adj", "BMXHT_adj")

cols_to_keep <- intersect(colnames(p_mv), cols_to_keep)

cols_to_keep <- c(cols_to_keep, summ_stats_variables |> pull(evarname), pvarname_to_query)

p_mv <- p_mv |> select(all_of(cols_to_keep))
file_out <- sprintf('%s_mv.rds', pvarname_to_query)

big_data <- list(big_data=p_mv, summ_stats=summ_stats, phenotype=pvarname_to_query)




################# impute and R2
cat("now imputing and estimating R2\n")
big_missing_gl <- missing_glimpse(big_data$big_data)
bd <- big_data$big_data
summ_stat <- big_data$summ_stats
phenotype <- big_data$phenotype

## select n=10 highest by R2
M_TOP <- 20
exposures_selected <- summ_stat |> arrange(-rsq_adj) |> slice_head(n=M_TOP)


if(nrow(exposures_selected) == 0) {
  message("no exposures above the threshold of missingness to impute, qutting")
  quit("no", 0)
}

cols <- c("SDMVPSU", "WTMEC2YR", "INDFMPIR", "SDDSRVYR", "SDMVSTRA",
          "ETHNICITY_NONHISPANICWHITE", "ETHNICITY_MEXICAN", "ETHNICITY_OTHERHISPANIC",  "ETHNICITY_OTHER", "ETHNICITY_NONHISPANICBLACK", "AGE_SQUARED", "RIAGENDR", "RIDAGEYR",
          "EDUCATION_LESS9",  "EDUCATION_9_11", "EDUCATION_HSGRAD", "EDUCATION_AA", "EDUCATION_COLLEGEGRAD", "BORN_INUSA")

if(ADDTC) {
  cols <- c(cols, "LBXTC")
}

bd_exposure_tibble <- bd |> select(all_of(c("SEQN", phenotype, cols, exposures_selected$evarname)))
bd_exposure_tibble <- bd_exposure_tibble |> rename(pheno=all_of(phenotype))
bd_exposure_tibble <- bd_exposure_tibble |> filter(!is.na(pheno)) |> filter(!is.na(INDFMPIR))

bd_m <- mice(bd_exposure_tibble, m = 10, print=F)

base_formula <- "I(scale(pheno)) ~ RIAGENDR + RIDAGEYR + AGE_SQUARED + ETHNICITY_NONHISPANICBLACK + ETHNICITY_OTHER + ETHNICITY_MEXICAN + ETHNICITY_OTHERHISPANIC + EDUCATION_LESS9 + EDUCATION_9_11 + EDUCATION_HSGRAD"
if(ADDTC) {
  base_formula <- sprintf("%s + %s", base_formula, "ns(log(LBXTC), df = 3)")
}


extended_formula <- sprintf('%s + %s', base_formula, paste(exposures_selected$evarname, collapse='+'))

fit1 <- with(bd_m, lm(formula(base_formula)))
fit2 <- with(bd_m, lm(formula(extended_formula)))

r2_1 <- pool.r.squared(fit1)
r2_2 <- pool.r.squared(fit2)

pooled_results <- pool(fit2)

# Multivariate R2
mv_diff <- r2_2[1] - r2_1[1]
print(mv_diff)


uni_r2 <- exposures_selected |> summarize(s=sum(rsq_adj)) |> pull(s)


to_save <- list(rsq_summary=tibble(phenotype=phenotype, rsq_exposures=mv_diff, rsq_base=r2_1[1], rsq_uni=uni_r2, number_exposures_in_model=nrow(exposures_selected)),
                exposures_selected=exposures_selected, fit1=r2_1, fit2=r2_2, pooled_fit=pooled_results)

file_out <- sprintf("./out/%s_imp_R2.rds", phenotype)
write_rds(to_save, file=file_out)