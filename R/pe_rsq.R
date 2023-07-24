# Chirag Patel
# estimate R2 from ExWAS significant variables
# 5/22/23



# Source File and Variable Definitions
#
# The adjustmentVariables are the covariates used for adjustment in the statistical models.
#
# The surveyVariables are additional variables that might be used in the survey design, such as sampling weights and stratification variables.
#
# search_colnames_for_weights function
# This function takes in a table (tab) and returns the name of the column that is assumed to represent the sampling weight. The function first looks for any column names that start with 'WT'. If multiple 'WT' column names are found, it applies a series of rules to select the most appropriate weight column.
#
# get_mv_expo_pheno_tables_for_big_table function
# This function retrieves and merges multiple tables of exposure and phenotype data from a database, based on the series that the phenotype belongs to.
#
# get_individual_level_table function
# This function retrieves the individual-level data for a given phenotype variable name from a summary statistics database and an NHANES database. It then adjusts the data based on a number of criteria and returns a list with the final set of exposure variables, the selected weight to use, and the selected data.
#
# unweighted_r2 function
# This function calculates the R-squared value for two linear models. The first model includes only the adjustment variables as predictors, and the second model includes both the adjustment variables and the exposure variables. This function does not use survey weights.
#
# svy_weighted_r2 function
# This function calculates the R-squared value for two linear models as in the unweighted_r2 function, but it takes into account the survey design by using the survey package in R to calculate a weighted R-squared.

# source('quantpe.R')
#adjustmentVariables <- c("RIDAGEYR", "AGE_SQUARED",
#                         "RIAGENDR",
#                         "INDFMPIR",
#                         "EDUCATION_LESS9","EDUCATION_9_11","EDUCATION_AA","EDUCATION_COLLEGEGRAD",
                         #"BORN_INUSA",
#                         "ETHNICITY_MEXICAN", "ETHNICITY_OTHERHISPANIC","ETHNICITY_OTHER", "ETHNICITY_NONHISPANICBLACK")

#surveyVariables <- c('WTMEC2YR', 'WTMEC4YR', 'SDMVPSU', 'SDMVSTRA', 'SDDSRVYR', 'SDDSRVYR')
#surveyVariables <- c('WTMEC2YR', 'SDMVPSU', 'SDMVSTRA', 'SDDSRVYR', 'SDDSRVYR')



#' Retrieve multivariate exposure and phenotype tables for a large table
#'
#' This function retrieves exposure and phenotype tables based on the provided table names,
#' performs several operations including joining and filtering, and returns a list containing
#' the final processed tables and relevant information.
#'
#' @param con A DBI connection object to the database.
#' @param pheno_table_name The name of the phenotype table.
#' @param expo_table_names A vector containing the names of the exposure tables.
#'
#' @return A list containing the final merged table (`merged_tab`), the exposure table (`e_table`),
#' the phenotype table (`p_table`), the series name (`series`), and a tibble containing information
#' about the table weights (`table_weights`).
#' @export
#' @examples
#' \dontrun{
#' result <- get_mv_expo_pheno_tables_for_big_table(con, "pheno_table", c("expo_table1", "expo_table2"))
#' }
#' @importFrom DBI dbGetQuery
#' @importFrom dplyr tbl filter pull collect select full_join inner_join left_join mutate_if
#' @importFrom tidyr tibble
get_mv_expo_pheno_tables_for_big_table <- function(con, pheno_table_name, expo_table_names) {
  table_names <- tbl(con, "table_names_epcf")
  seriesName <- table_names |> filter(Data.File.Name == pheno_table_name) |> pull(series)
  log_info("Series of phenotype is { seriesName } ")
  demo_table_name <- table_names |> filter(component == 'DEMO', series==seriesName) |> collect() |> pull(Data.File.Name)
  log_info("Demographics table is { demo_table_name } ")
  demo <- tbl(con, demo_table_name) |> collect()

  e_table <- NULL
  table_weight_names <- tibble()
  for(ii in 1:length(expo_table_names)) {
    e_table_lcl <- tbl(con, expo_table_names[ii]) |> collect()
    log_info("Exposure table {expo_table_names[ii]} has {e_table_lcl |> count() |> pull(n) } rows")
    if(ii == 1) {
      e_table <- e_table_lcl
      nme <- search_colnames_for_weights(e_table_lcl)
      n <- length(unique(e_table_lcl$SEQN))
      table_weight_names <- rbind(table_weight_names, tibble(table_name=expo_table_names[ii], weight_name=nme, sample_size=n))
      next;
    }
    prev_colnames <- colnames(e_table)
    nme <- search_colnames_for_weights(e_table_lcl)
    n <- length(unique(e_table_lcl$SEQN))
    table_weight_names <- rbind(table_weight_names, tibble(table_name=expo_table_names[ii], weight_name=nme, sample_size=n))

    cols_to_keep <- c(setdiff(colnames(e_table_lcl), prev_colnames), "SEQN")
    e_table <- e_table |> full_join(e_table_lcl |> select(all_of(cols_to_keep)), by="SEQN")
  }

  #log_info("New exposure table has {e_table |> collect() |> count() |> pull(n) } rows")

  p_table <- tbl(con, pheno_table_name) |> collect()
  nme <-  search_colnames_for_weights(p_table)
  n <- length(unique(p_table$SEQN))
  table_weight_names <- rbind(table_weight_names, tibble(table_name=expo_table_names[ii], weight_name=nme, sample_size=n))
  log_info("Pheno table {pheno_table_name} has {p_table |> count() |> pull(n) } rows")
  small_tab <- demo |> inner_join(p_table, by="SEQN") |> left_join(e_table, by="SEQN")
  log_info("Merged table has { small_tab |> count() |> pull(n) } rows")

  small_tab <- small_tab |> mutate_if(is.not.numeric, as.numeric)
  return(list(merged_tab=small_tab, e_table=e_table, p_table=p_table, series=seriesName, table_weights=table_weight_names))
}

#' Retrieve individual level data table based on a set of a variables
#'
#' Assumes that findings have already been associated with a phenotype: |> filter(sig_levels == 'Bonf.<0.05')
#'
#' This function queries individual level data from a provided database connection
#' based on specific parameters and performs several operations including joining,
#' filtering, summarizing, and selecting of data, before returning the final processed data.
#'
#' @param summary_stats_con A DBI connection object to the summary statistics database.
#' @param nhanes_con A DBI connection object to the NHANES database.
#' @param pvarname_to_query The variable name to be queried from the database.
#'
#' @return A list containing the final selected variables (`evars`), the selected weight to use (`weight_to_use`),
#' and the selected data (`selected_data`).
#' @export
#' @examples
#' \dontrun{
#' result <- get_individual_level_table(summary_stats_con, nhanes_con, "varname")
#' }
#' @importFrom DBI dbGetQuery
#' @importFrom dplyr tbl filter select inner_join group_by summarize collect rowwise mutate ungroup pull bind_rows arrange
#' @importFrom tidyr unnest_wider gather
get_individual_level_table <- function(summary_stats_con, nhanes_con, pvarname_to_query) {
  adjusted_meta_2 <- tbl(summary_stats_con, "adjusted_meta_2")
  expos_wide <- tbl(summary_stats_con, "expos_wide")
  evars_for_p <- adjusted_meta_2 |> filter(pvarname == pvarname_to_query) |> filter(sig_levels == 'Bonf.<0.05')  # evariables
  table_list <- expos_wide |> filter(pvarname == pvarname_to_query) # |>  filter(term == 'expo')
  table_list <- table_list |> inner_join(evars_for_p |> select(evarname)) |> select(evarname, pvarname, exposure_table_name, ecategory,pcategory, phenotype_table_name, nobs_adjusted)


  ep_tables <- table_list |> collect() |> select(exposure_table_name, phenotype_table_name) |> unique() |> group_by(phenotype_table_name) |> summarize(exposure_table_names=list(exposure_table_name))
  ep_tables <- ep_tables|> rowwise() |>
    mutate(m_struct=list(get_mv_expo_pheno_tables_for_big_table(nhanes_con, phenotype_table_name, exposure_table_names)))

  ep_tables <- ep_tables |> ungroup() |> unnest_wider(m_struct)

  weights <- ep_tables |> pull(table_weights) |> bind_rows()
  uniq_weights <- weights |> group_by(weight_name) |> summarize(combined_sample_size=sum(sample_size))

  ss_per_variable <- table_list |> group_by(evarname) |> summarize(total_n=sum(nobs_adjusted))
  the_big_one <- ep_tables |> pull(merged_tab) |> bind_rows()

  big_samples <- ss_per_variable |> filter(total_n >= 10000)

  big_cc <- the_big_one |> select(
    SEQN, all_of(pvarname_to_query), all_of(surveyVariables), all_of(uniq_weights$weight_name),
    all_of(adjustmentVariables), all_of( big_samples |> pull(evarname))) |> filter(!is.na(get(pvarname_to_query)))


  missingness <- big_cc |>
    summarise_all(~(sum(is.na(.))/length(.)*100)) |> gather(variable, missing_perc) |> arrange(missing_perc)
  selected_vars <- missingness |> filter(missing_perc < 30) |> pull(variable)
  selected_data <- big_cc |> select(all_of(selected_vars)) |> stats::na.omit()

  selected_weights <- intersect(selected_vars, uniq_weights$weight_name)
  weight_to_use <- uniq_weights |> filter(weight_name %in% selected_weights) |> slice_min(combined_sample_size)

  log_info("MV table has {selected_data |> count() |> pull(n) } rows")
  evars_final <- setdiff(selected_vars, c(pvarname_to_query, "SEQN", selected_weights, surveyVariables, adjustmentVariables))
  log_info("MV table has {evars_final |> length() } number of E")
  log_info("MV table SDDSRVYR: {unique(selected_data$SDDSRVYR) }")
  list(evars=evars_final, weight_to_use=weight_to_use, selected_data=selected_data)
}

search_colnames_for_weights <- function(tab) {
  weight_names <- colnames(tab)[grep("^WT", colnames(tab))]
  weight_name <- 'WTMEC2YR' ## default demographic
  if(length(weight_names) > 1) {
    if('WTDRD1' %in% weight_names) { ## dietary variable
      weight_name <- 'WTDRD1'
    } else if(any(grepl('2YR$', weight_names))) {
      weight_name <- weight_names[grep('2YR$', weight_names)]
    } else {
      weight_name <- weight_names[1]
    }
  }
  if(length(weight_names) ==1) {
    weight_name <- weight_names[1]
  }
  return(weight_name)
}

is.not.numeric <- function(x, ...) {
  !is.numeric(x, ...)
}


unweighted_r2 <- function(pvarname, evarnames, adjustment_variables, data) {
  baseformula <- as.formula(sprintf("%s ~ %s", pvarname, paste(adjustment_variables, collapse="+")))
  baseevars <- addToBase(baseformula, adjustingVariables = evarnames)
  gl1 <- glance(lm(baseformula , selected_data))
  gl2 <- glance(lm(baseevars, selected_data))
  return(list(base=gl1, mve=gl2))
}

svy_weighted_r2 <- function(pvarname, evarnames, adjustment_variables, dat, weight_name="WTMEC2YR") {
  dat <- dat |> mutate(wt = !!sym(weight_name)) |> mutate(wt = wt / length(dat$SDDSRVYR)) # assume 1999 and 2001 are not being analyzed

  dsn <- svydesign(ids=~SDMVPSU, strata=~SDMVSTRA, weights=~wt, nest=T, data = dat |> filter(wt > 0))
  baseformula <- as.formula(sprintf("%s ~ %s", pvarname, paste(adjustment_variables, collapse="+")))
  baseevars <- addToBase(baseformula, adjustingVariables = evarnames)
  modbase <- svyglm(baseformula, dsn)
  modevars <- svyglm(baseevars, dsn)
  r2base <- svyrsquared(modbase)
  r2evars <- svyrsquared(modevars)
  return(list(base=r2base, mve=r2evars))
}
