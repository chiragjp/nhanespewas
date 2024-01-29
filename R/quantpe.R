## Chirag J Patel
## 07/14/23
## quantpe.R
## machinery to associate a real-valued outcome with a categorical or real-valued exposure


surveyVariables <- c('WTMEC2YR', 'WTMEC4YR', 'SDMVPSU', 'SDMVSTRA', 'SDDSRVYR')

seriesBeginYearMap <- function(seriesName) {
    seriesBeginYear <- dplyr::case_when(
    seriesName == 'A' ~ 1999,
    seriesName == 'B' ~ 2001,
    seriesName == 'C' ~ 2003,
    seriesName == 'D' ~ 2005,
    seriesName == 'E' ~ 2007,
    seriesName == 'F' ~ 2009,
    seriesName == 'G' ~ 2011,
    seriesName == 'H' ~ 2013,
    seriesName == 'I' ~ 2015,
    seriesName == 'J' ~ 2017,
    seriesName == 'K' ~ 2019,
    seriesName == 'L' ~ 2021,
    TRUE ~ NA_real_
  )
  seriesBeginYear
}


adjustment_variables <- function() {
  varnames <- c("RIDAGEYR", "RIAGENDR", "INDFMPIR", "ETHNICITY_MEXICAN",
                "ETHNICITY_OTHERHISPANIC", "ETHNICITY_OTHER",
                "ETHNICITY_NONHISPANICBLACK", "ETHNICITY_NONHISPANICWHITE",
                "EDUCATION_LESS9", "EDUCATION_9_11", "EDUCATION_HSGRAD",
                "EDUCATION_AA", "EDUCATION_COLLEGEGRAD")
}

get_path_to_extdata_database<- function() {
  path_to_db <- system.file("extdata", "nhanes_pewas_a-d.sqlite", package = "nhanespewas")
}


#' Connect to the PEWAS NHANES database
#'
#' This function establishes a connection to the PEWAS NHANES database .
#'
#' @param path_to_data The file path to the participant-level NHANES data. Each table in the database corresponds to a table in the NHANES.
#'
#' @return A DBI connection object to the NHANES database.
#' @export
#' @examples
#' \dontrun{
#' con <- connect_pewas_data("/path/to/database.sqlite")
#' }
#' @importFrom DBI dbConnect
#' @importFrom RSQLite SQLite
connect_pewas_data <- function(path_to_data = NULL) {

  ## check if the database exists
  if (is.null(path_to_data)) {
    path_to_data <- get_path_to_extdata_database()
  } else if(!file.exists(path_to_data)) {
    stop("nhanes pewas database does not exist")
  }
  con <- DBI::dbConnect(RSQLite::SQLite(), dbname=path_to_data)
  # check version date and type of database (e.g., summary stats or NHANES raw)
  if(DBI::dbExistsTable(con, "description")) {
    desc <- dplyr::tbl(con, "description")
    return(con)
  } else{
    stop("nhanes pewas database is not in correct format")
  }

}

#' Disconnect from a NHANES database
#'
#' This function disconnects from a database given a DBI connection object.
#'
#' @param dbConn A NHANES DBI connection object.
#'
#' @return Invisible NULL. The connection will be completely severed.
#' @export
#' @examples
#' \dontrun{
#' disconnect_pewas_data(dbConn)
#' }
#' @importFrom DBI dbDisconnect
disconnect_pewas_data <- function(dbConn) {
  DBI::dbDisconnect(dbConn)
}

#' Retrieve and merge tables from a database based on the series name, exposure variable name, and phenotype variable name
#'
#' @param pvarname Character. The phenotype variable name.
#' @param evarname Character. The exposure variable name.
#' @param seriesName Character. The series name.
#' @param con The database connection.
#' @param pheno_table_name Character. The optional phenotype table name. If NULL (default), the function will attempt to find the table based on the `pvarname`.
#' @param expo_table_name Character. The optional exposure table name. If NULL (default), the function will attempt to find the table based on the `evarname`.
#'
#' @return A list containing the merged table, the exposure table, the phenotype table, and the series name.
#' @examples
#' \dontrun{
#' conn <- connect_pewas_data(...)
#' tables <- get_tables("LBXGLU", "LBXGTC", "C", conn, "L10AM_C", "L45VIT_C")
#' }
#' @export
get_tables <- function(pvarname, evarname, seriesName, con, pheno_table_name=NULL, expo_table_name=NULL) {
  variables <- dplyr::tbl(con, "variable_names_epcf")
  table_names <- dplyr::tbl(con, "table_names_epcf")
  series_begin_year <- seriesBeginYearMap(seriesName)
  demo_table_name <- table_names |> dplyr::filter(series == seriesName & component == 'DEMO') |> dplyr::collect()
  demo <- dplyr::tbl(con, demo_table_name$Data.File.Name)
  ## get table for the phenotype
  logger::log_info("Getting tables for {pvarname} ~ {evarname}")

  p_table <- tibble::tibble()
  if(!is.null(pheno_table_name)) {
    #p_table_name <- p_table_name |> filter(Data.File.Name == pheno_table_name)
    p_table <- tbl(con, pheno_table_name)
  } else{
    p_table_name <- variables |> dplyr::filter(Variable.Name == pvarname & Begin.Year == series_begin_year) |> dplyr::collect()
    pheno_table_name <- p_table_name$Data.File.Name
    p_table <- dplyr::tbl(con, pheno_table_name)
  }

  logger::log_info("Pheno table {pheno_table_name} has {p_table |> count() |> pull(n) } rows")
  ## get table for exposure
  e_table <- tibble::tibble()
  if(!is.null(expo_table_name)) {
    e_table <- dplyr::tbl(con, expo_table_name)
  } else {
    e_table_name <- variables |> dplyr::filter(Variable.Name == evarname & Begin.Year == series_begin_year) |> dplyr::collect()
    expo_table_name <- e_table_name$Data.File.Name
    e_table <- dplyr::tbl(con, expo_table_name)
  }

  if(expo_table_name == pheno_table_name) {
    p_table <- p_table |> dplyr::select(SEQN, pvarname)
  }

  etable_nr <- e_table |> dplyr::count() |> dplyr::pull(n)
  logger::log_info("Exposure table {expo_table_name} has {etable_nr} rows")
  small_tab <- demo |> dplyr::inner_join(p_table, by="SEQN") |> dplyr::inner_join(e_table, by="SEQN") |> dplyr::collect()
  logger::log_info("Merged table has {small_tab |> count() |> pull(n) } rows")
  return(list(merged_tab=small_tab, e_table=e_table, p_table=p_table, series=seriesName))
}


#' Retrieve and merge demographic, exposure, and phenotype tables from a database
#'
#' @param con The database connection.
#' @param pheno_table_name Character. The phenotype table name.
#' @param expo_table_name Character. The exposure table name.
#'
#' @return A list containing the merged table, the exposure table, the phenotype table, and the series name.
#' @examples
#' \dontrun{
#' conn <- connect_pewas_data()
#' results <- get_expo_pheno_tables(conn, pheno_table_name="L10AM_C", expo_table_name="L45VIT_C")
#' }
#' @export
get_expo_pheno_tables <- function(con, pheno_table_name, expo_table_name) {
  table_names <- dplyr::tbl(con, "table_names_epcf")
  seriesName <- table_names |> dplyr::filter(Data.File.Name == pheno_table_name) |> dplyr::pull(series)
  logger::log_info("Series of phenotype is { seriesName } ")
  demo_table_name <- table_names |> dplyr::filter(component == 'DEMO', series==seriesName) |> dplyr::collect() |> dplyr::pull(Data.File.Name)
  logger::log_info("Demographics table is { demo_table_name } ")
  demo <- dplyr::tbl(con, demo_table_name)
  e_table <- dplyr::tbl(con, expo_table_name)
  e_table_nrow <- e_table |> collect() |> nrow()
  logger::log_info("Exposure table {expo_table_name} has { e_table_nrow } rows")
  p_table <- dplyr::tbl(con, pheno_table_name)
  if(pheno_table_name == expo_table_name) { # hack to preserve data structure
    p_table <- p_table |> dplyr::select(SEQN)
  }
  p_table_nrow <- p_table |> collect() |> nrow()
  logger::log_info("Pheno table {pheno_table_name} has {p_table_nrow } rows")
  small_tab <- demo |> dplyr::inner_join(p_table, by="SEQN") |> dplyr::inner_join(e_table, by="SEQN") |> dplyr::collect()
  small_tab_nrow <- small_tab |> collect() |> nrow()
  logger::log_info("Merged table has { small_tab_nrow } rows")
  return(list(merged_tab=small_tab, e_table=e_table, p_table=p_table, series=seriesName))
}


#' Retrieve and merge demographic, multi-variable exposure, and phenotype tables from a database
#'
#' This function retrieves and merges tables from a database based on provided table names. Specifically, it merges demographic, multiple exposure variables and phenotype tables. The function logs the series of the phenotype, the name of the demographics table, the number of rows in each exposure and phenotype table, and finally the number of rows in the merged table.
#'
#' @param con The database connection.
#' @param pheno_table_name Character. The name of the phenotype table.
#' @param expo_table_names Character vector. The names of the exposure tables.
#'
#' @return A list containing the merged table, the final exposure table after joining all exposure tables, the phenotype table, and the series name.
#' @examples
#' \dontrun{
#' conn <- connect_pewas_data()
#' results <- get_mv_expo_pheno_tables(conn, "L10AM_C", c("L45VIT_C", "L06COT_C"))
#' }
get_mv_expo_pheno_tables <- function(con, pheno_table_name, expo_table_names) {
  table_names <- dplyr::tbl(con, "table_names_epcf")
  seriesName <- table_names |> dplyr::filter(Data.File.Name == pheno_table_name) |> dplyr::pull(series)
  logger::log_info("Series of phenotype is { seriesName } ")
  demo_table_name <- table_names |> dplyr::filter(component == 'DEMO', series==seriesName) |> dplyr::collect() |> dplyr::pull(Data.File.Name)
  logger::log_info("Demographics table is { demo_table_name } ")
  demo <- dplyr::tbl(con, demo_table_name) |> dplyr::collect()

  e_table <- NULL
  for(ii in 1:length(expo_table_names)) {
    e_table_lcl <- dplyr::tbl(con, expo_table_names[ii]) |> dplyr::collect()
    logger::log_info("Exposure table {expo_table_names[ii]} has {e_table_lcl |> count() |> pull(n) } rows")
    if(ii == 1) {
      e_table <- e_table_lcl
      next;
    }
    prev_colnames <- colnames(e_table)
    cols_to_keep <- c(setdiff(colnames(e_table_lcl), prev_colnames), "SEQN")
    e_table <- e_table |> dplyr::inner_join(e_table_lcl |> dplyr::select(tidyselect::all_of(cols_to_keep)), by="SEQN")
  }

  #log_info("New exposure table has {e_table |> collect() |> count() |> pull(n) } rows")

  p_table <- dplyr::tbl(con, pheno_table_name) |> dplyr::collect()

  logger::log_info("Pheno table {pheno_table_name} has {p_table |> count() |> pull(n) } rows")
  small_tab <- demo |> dplyr::inner_join(p_table, by="SEQN") |> dplyr::inner_join(e_table, by="SEQN")
  logger::log_info("Merged table has { small_tab |> count() |> pull(n) } rows")

  return(list(merged_tab=small_tab, e_table=e_table, p_table=p_table, series=seriesName))
}



#' Identify the appropriate weight column from the exposure and phenotype tables
#'
#' This function identifies the appropriate weight column from the exposure (e_table) and phenotype (p_table) tables, contained within the list object produced by the `get_tables()` function. The weight column is determined based on the following criteria:
#' 1. If a column named "WTDRD1" exists, it is selected.
#' 2. If not, and a column ending in "2YR" exists, that is selected.
#' 3. If neither of the above conditions are met, the first weight column is selected.
#' If there is no weight column in either table, a default weight column ("WTMEC2YR") is used.
#' Once the weight column is identified, a new column "wt" is added to the merged table (merged_tab), which is equivalent to the weight column.
#'
#' @param get_tables_obj A list object produced by the `get_tables()` function. Contains the exposure table (e_table), phenotype table (p_table), merged table (merged_tab), and series name.
#' @return The input list object with an additional column "wt" in the merged table (merged_tab) corresponding to the identified weight column.
#' @examples
#' \dontrun{
#' conn <- connect_pewas_data()
#' get_tables_results <- get_tables("LBXGLU", "LBXGTC", "C", conn, "L10AM_C", "L45VIT_C")
#' weighted_tables <- figure_out_weight(get_tables_results)
#' }
#' @export
figure_out_weight <- function(get_tables_obj) {
  ##  need to select weights that are not the bootstrap weights
  weight_name_demo <- 'WTMEC2YR'
  e_table <- get_tables_obj$e_table
  p_table <- get_tables_obj$p_table

  weight_name_e <- colnames(e_table)[grep("^WT", colnames(e_table))]
  weight_name_p <- colnames(p_table)[grep("^WT", colnames(p_table))]

  ## filter here for a single weight
  ## if e_weight is empty, move on
  ## if e_weight  has dietary weight (wtdrd1), select that (vs. the 2 day weight)
  ## else if e_weight has 2yr weight select that
  ## else select the first weight

  if(length(weight_name_e) > 1) {
    if('WTDRD1' %in% weight_name_e) { ## dietary variable
      weight_name_e <- 'WTDRD1'
    } else if(any(grepl('2YR$', weight_name_e))) {
      weight_name_e <- weight_name_e[grep('2YR$', weight_name_e)]
    } else {
      weight_name_e <- weight_name_e[1]
    }
  }

  if(length(weight_name_p) > 1) {
    if('WTDRD1' %in% weight_name_p) {
      weight_name_p <- 'WTDRD1'
    } else if(any(grepl('2YR$', weight_name_p))) {
      weight_name_p <- weight_name_p[grep('2YR$', weight_name_p)]
    } else {
      weight_name_p <- weight_name_p[1]
    }
  }

  logger::log_info("p weight name: { weight_name_p }")
  logger::log_info("e weight name: { weight_name_e }")
  etable_nrows <- e_table |> dplyr::count() |> pull(n)
  ptable_nrows <- p_table |> dplyr::count() |> pull(n)

  if(rlang::is_empty(weight_name_e) & rlang::is_empty(weight_name_p)) {
    logger::log_info("no weights in e or p table")
    weight_name <- weight_name_demo
  } else if(!rlang::is_empty(weight_name_e) & !rlang::is_empty(weight_name_p)) {
    logger::log_info("weights in both e and p table")

    if(etable_nrows < ptable_nrows) {
      weight_name <-  weight_name_e
    } else {
      weight_name <-  weight_name_p
    }

    if(weight_name_e == weight_name_p) {
      weight_name <- sprintf('%s.y', weight_name)
    }

  } else if(!rlang::is_empty(weight_name_p)) {
    logger::log_info("weight in p table")
    weight_name <- weight_name_p
  } else if (!rlang::is_empty(weight_name_e)) {
    logger::log_info("weight in e table")
    weight_name <- weight_name_e
  } else {
    weight_name <- weight_name_demo
  }
  logger::log_info("final weight name: { weight_name }")
  get_tables_obj$merged_tab <- get_tables_obj$merged_tab |> dplyr::mutate(wt = !!as.name(weight_name))
  get_tables_obj
}



log10_xform_variable <- function(x) {
  # are there zeros?
  if(any(x==0, na.rm = T)) {
    x[x==0 & !is.na(x)] <- sqrt(min(x[x>0 & !is.na(x)]))  # replace with the sqrt of the minimum
  }
  log10(x)
}


comment_from_exponame <- function(exposure_name) {
  corevar <- substr(exposure_name, 4, 100)
  firsttwo <- substr(exposure_name, 1, 2)
  lcName <- sprintf("%sD%sLC", firsttwo, corevar)
  lcName
}

create_svydesign <- function(data) {
  options(survey.lonely.psu="adjust")
  survey::svydesign(ids=~SDMVPSU, strata=~SDMVSTRA, weights=~wt, nest=T, data=data)
}

addToBase <- function(base_formula, adjustingVariables) {
  form <- base_formula
  if(length(adjustingVariables)) {
    addStr <- stats::as.formula(sprintf('~ . + %s', paste(adjustingVariables, collapse='+')))
    form <- stats::update.formula(base_formula, addStr)
  }
  return(form)
}



#' Calculate R-squared and adjusted R-squared for a survey-weighted linear model
#'
#' This function calculates the R-squared and adjusted R-squared values for a
#' linear model based on an analysis object which contains a survey design,
#' residuals, and a response variable. The calculation is conducted using
#' survey-weighted means, totals, and residuals.
#'
#' @param analysisObj A list containing the components necessary for the calculation.
#' It should have the following components:
#'   * survey.design: A survey design object.
#'   * residuals: A numeric vector of residuals from the model.
#'   * y: A numeric vector of the response variable.
#'
#' @return A list with two components:
#'   * rsq: R-squared of the linear model.
#'   * adj.r2: Adjusted R-squared of the linear model.
#'
#' @examples
#' \dontrun{
#' analysisObj <- list(survey.design = dsn, residuals = resid, y = y)
#' results <- svyrsquared(analysisObj)
#' }
svyrsquared <- function(analysisObj) {
  ## returns rsquared and adjusted rsquared for a linear model
  dsn <- analysisObj$survey.design
  resid <- analysisObj$residuals
  yy <- analysisObj$y
  dsn <- stats::update(dsn, y=yy)
  ymean <- survey::svymean(~y, dsn)[1]
  dsn <- stats::update(dsn, ydiff=(yy-ymean)^2)
  dsn <- stats::update(dsn, residual=resid^2)
  SSE <- survey::svytotal(~residual, dsn)
  SST <- survey::svytotal(~ydiff, dsn)
  rsquared <- 1 - (SSE[1] /SST[1])
  ## compute adjusted R squared
  N <- length(yy)
  P <- length(analysisObj$coefficients)
  adjrsquared <- 1-( (N-1)*(1-rsquared)/(N-P) )
  return(list(rsq=rsquared, adj.r2=adjrsquared))
}


#' Run a model with options for scaling and quantile-based categorization of exposure and phenotype variables
#'
#' This function runs a survey-weighted generalized linear model based on the given formula and design,
#' with options for scaling the exposure and phenotype variables and for categorizing the exposure variable
#' based on quantiles or given levels.
#'
#' @param formu A formula specifying the model to be run.
#' @param dsn A survey design object.
#' @param scale_expo Logical, if TRUE, the exposure variable is scaled by its standard deviation.
#' @param scale_pheno Logical, if TRUE, the phenotype variable is scaled by its standard deviation.
#' @param quantile_expo Numeric vector, quantiles to use for categorizing the exposure variable.
#' If not NULL, this overwrites the scale_expo argument.
#' @param expo_levels A vector of levels to categorize the exposure variable.
#' @param save_dsn Logical, if TRUE, the desgin object for each model will be saved
#' If not NULL, this overwrites the quantile_expo and scale_expo arguments.
#'
#' @return A list containing:
#'   * model: The model object (if the parameter is set to T)
#'   * glanced: A one-row data frame of model-level statistics from glance().
#'   * tidied: A data frame of term-level statistics from tidy().
#'   * r2: R-squared and adjusted R-squared from svyrsquared().
#'   * scale_pheno: Whether the phenotype variable was scaled.
#'   * scale_expo: Whether the exposure variable was scaled.
#'   * q_cut_points: Cut points used for quantile-based categorization of the exposure variable.
#'   * quantiles: The input quantiles for categorizing the exposure variable.
#'   * expo_levels: The input levels for categorizing the exposure variable.
#'
#' @examples
#' \dontrun{
#' result <- run_model(formu, dsn, scale_expo=T, scale_pheno=F, quantile_expo=c(0.25, 0.5, 0.75))
#' }
#' @export
run_model <- function(formu, dsn, scale_expo=T, scale_pheno=F, quantile_expo=NULL, expo_levels=NULL, save_svymodel=F) {
  cut_points <- NA
  if(scale_expo) {
    logger::log_info("Exposure as scaled")
    s <- sqrt(mean(survey::svyvar(~expo, dsn, na.rm = T)))
    mn <- mean(survey::svymean(~expo, dsn, na.rm=T))
    dsn <- stats::update(dsn, expo=(expo-mn)/s)
  } else if(!is.null(quantile_expo)) {
    logger::log_info("Exposure as a cut variable by quantiles")
    q <- survey::svyquantile(~expo, quantiles=quantile_expo, design = dsn)
    cut_points <- coef(q)
    labs <- paste("q", 1:(length(cut_points)-1), sep = "")
    dsn <- stats::update(dsn, expo=cut(expo, breaks=cut_points, labels = labs))
  } else if(!is.null(expo_levels)) {
    logger::log_info("Exposure as a categorical variable with levels { paste(expo_levels, collapse=';') }")
    dsn <- stats::update(dsn, expo=factor(expo, levels=expo_levels))
  }

  if(scale_pheno) {
    s <- sqrt(mean(survey::svyvar(~pheno, dsn, na.rm = T)))
    mn <- mean(survey::svymean(~pheno, dsn, na.rm=T))
    dsn <- stats::update(dsn, pheno=(pheno-mn)/s)
  }

  mod <- survey::svyglm(formu,dsn)
  ti <- broom::tidy(mod)
  gl <- broom::glance(mod)
  r2 <- svyrsquared(mod)
  list(glanced=gl, tidied=ti, r2=r2,
       scale_pheno=scale_pheno, scale_expo=scale_expo,
       q_cut_points=cut_points, quantiles=quantile_expo, expo_levels=expo_levels, model=ifelse(save_svymodel, mod, NA)) # not saving dsn
}

run_mv_model <- function(formu, dsn, scale_pheno=F) {

  if(scale_pheno) {
    s <- sqrt(mean(survey::svyvar(~pheno, dsn, na.rm = T)))
    mn <- mean(survey::svymean(~pheno, dsn, na.rm=T))
    dsn <- stats::update(dsn, pheno=(pheno-mn)/s)
  }

  mod <- survey::svyglm(formu,dsn)
  ti <- broom::tidy(mod)
  gl <- broom::glance(mod)
  r2 <- svyrsquared(mod)
  list(model=mod, r2=r2, # tidy and glanced do not work for large models (degree of freedom will be negative - see ?svyglm) - use summary(model, df.resid=Inf)
       dsn=dsn)
}


name_and_xform_pheno_expo <- function(pheno, exposure, table_object, logxform_p=T, logxform_e=T) {

  if(logxform_p) {
    table_object$merged_tab <- table_object$merged_tab |> dplyr::mutate(pheno=log10_xform_variable(!!as.name(pheno)))
  } else {
    table_object$merged_tab <- table_object$merged_tab |> dplyr::mutate(pheno=!!as.name(pheno))
  }
  if(logxform_e) {
    ##
    table_object$merged_tab <- table_object$merged_tab |> dplyr::mutate(expo=(log10_xform_variable(!!as.name(exposure))))
  } else {
    table_object$merged_tab <- table_object$merged_tab |> dplyr::mutate(expo=(!!as.name(exposure)))
  }

  table_object
}

name_and_xform_pheno <- function(pheno, table_object, logxform_p=T) {

  if(logxform_p) {
    table_object$merged_tab <- table_object$merged_tab |> dplyr::mutate(pheno=log10_xform_variable(!!as.name(pheno)))
  } else {
    table_object$merged_tab <- table_object$merged_tab |> dplyr::mutate(pheno=!!as.name(pheno))
  }
  table_object
}



#' Execute an association between a P and E variable
#' This function performs a series of operations including fetching tables,
#' applying weight adjustments, performing log transformations, creating survey designs,
#' running models, and computing demographic breakdowns.
#'
#' @param pheno A variable related to phenotype
#' @param exposure A variable related to exposure
#' @param adjustment_variables A vector of variables to adjust for in the models. They have to be in the demo table.
#' @param series A variable related to series
#' @param con Connection object for the NHANES database
#' @param logxform_p Logical, if TRUE, a log transformation is applied to pheno variable. Default is TRUE.
#' @param logxform_e Logical, if TRUE, a log transformation is applied to exposure variable. Default is TRUE.
#' @param scale_e Logical, if TRUE, exposure variable is scaled. Default is TRUE.
#' @param scale_p Logical, if TRUE, pheno variable is scaled. Default is FALSE.
#' @param pheno_table_name Optional, name of the phenotype table
#' @param expo_table_name Optional, name of the exposure table
#' @param quantile_expo Optional, if specified, exposure variable will be cut by these quantiles
#' @param exposure_levels Optional, if specified, exposure variable will be treated as a categorical variable with these levels
#'
#' @return A list containing: the svydesign object, unweighted number of observations, phenotype, exposure, series, unadjusted model, adjusted model, base model, and demographic breakdown.
#'
#' @examples
#' \dontrun{
#' connection <- connect_pewas_data()
#' adjustments <- adjustment_variables()
#' pe_result <- pe(pheno="BMXBMI", exposure="LBXGTC", adjustment_variables=adjustments, series="C",
#' con=connection, logxform_p=T, logxform_e=T, scale_e=T, scale_p=F)
#' pe_result$adjusted_model$tidied
#' }
#'
#' @export

pe <- function(pheno, exposure, adjustment_variables, series, con,
               logxform_p=T, logxform_e=T, scale_e=T, scale_p=F,
               pheno_table_name=NULL, expo_table_name=NULL,
               quantile_expo=NULL, exposure_levels=NULL) {
  ## adjustment_variables has to be in the demo table

  ## get tables
  tab_obj <- get_tables(pheno, exposure, series, con, pheno_table_name, expo_table_name)
  ## weight
  tab_obj <- figure_out_weight(tab_obj)
  ## do logxforms and rename variables "expo" and "pheno"
  tab_obj <- name_and_xform_pheno_expo(pheno, exposure, tab_obj, logxform_p, logxform_e)
  ## create svydesign
  dat <- tab_obj$merged_tab |> dplyr::filter(!is.na(wt), wt > 0, !is.na(expo), !is.na(pheno), dplyr::if_all(tidyselect::all_of(adjustment_variables), ~!is.na(.)))
  dsn <- create_svydesign(dat)
  ## run models
  baseformula <- stats::as.formula("pheno ~ expo")
  baseadjusted <- addToBase(baseformula, adjustingVariables = adjustment_variables)
  basebase <- stats::as.formula(sprintf("pheno ~ %s", paste(adjustment_variables, collapse="+")))
  ##
  unadjusted_mod <- run_model(baseformula, dsn, scale_expo = scale_e, scale_pheno = scale_p,  quantile_expo=quantile_expo, expo_levels =  exposure_levels)
  adjusted_mod <- run_model(baseadjusted, dsn, scale_expo = scale_e, scale_pheno = scale_p, quantile_expo=quantile_expo, expo_levels = exposure_levels)
  base_mod <- run_model(basebase, dsn, scale_expo = scale_e, scale_pheno = scale_p, quantile_expo=quantile_expo, expo_levels = exposure_levels)
  n <- dsn |> nrow()
  ## demographic breakdown
  demo_break_tbl <- demographic_breakdown(dsn)
  ## return mods
  list(dsn=dsn, unweighted_n=n, phenotype=pheno, exposure=exposure, series=series, unadjusted_model=unadjusted_mod, adjusted_model=adjusted_mod, base_mod=base_mod, demo_breakdown=demo_break_tbl)
}


##

#' PE by Table
#'
#' This function operates similarly to the 'pe' function, but takes a pre-processed table object as input instead of separate phenotype and exposure tables. It is designed to help analysts scale up the pe associations by sepearting the "get_tables" procedure from runnin an association.
#'
#' @param tab_obj A pre-processed table object including both phenotype and exposure information
#' @param pvar A variable related to phenotype
#' @param evar A variable related to exposure
#' @param adjustment_variables A vector of variables to adjust for in the models
#' @param logxform_p Logical, if TRUE, a log transformation is applied to phenotype variable. Default is TRUE.
#' @param logxform_e Logical, if TRUE, a log transformation is applied to exposure variable. Default is TRUE.
#' @param scale_e Logical, if TRUE, exposure variable is scaled. Default is TRUE.
#' @param scale_p Logical, if TRUE, phenotype variable is scaled. Default is FALSE.
#' @param quantile_expo Optional, if specified, exposure variable will be cut by these quantiles
#' @param exposure_levels Optional, if specified, exposure variable will be treated as a categorical variable with these levels
#'
#' @return A list containing: the svydesign object, log transformation status for phenotype and exposure, scaling status for phenotype and exposure, unweighted number of observations, phenotype, series, exposure, unadjusted model, adjusted model, base model, and demographic breakdown.
#'
#' @examples
#' \dontrun{
#' con <- connect_pewas_data()
#' table_object <- get_tables("LBXGLU", "LBXGTC", "C", con, "L10AM_C", "L45VIT_C")
#' pe_result <- pe_by_table(tab_obj=table_object, pvar="LBXGLU", evar="LBXGTC",
#' adjustment_variables=c("RIDAGEYR", "RIAGENDR"),
#' logxform_p=T, logxform_e=T, scale_e=T, scale_p=F,
#' quantile_expo=NULL, exposure_levels=NULL)
#' }
#'
#' @export

pe_by_table <- function(tab_obj, pvar, evar,
                  adjustment_variables,
                  logxform_p=T, logxform_e=T, scale_e=T, scale_p=F,
                  quantile_expo=NULL, exposure_levels=NULL) {
  pheno <- pvar
  exposure <- evar
  tab_obj <- name_and_xform_pheno_expo(pheno, exposure, tab_obj, logxform_p, logxform_e)
  ## create svydesign
  dat <- tab_obj$merged_tab |> dplyr::filter(!is.na(wt), wt > 0, !is.na(expo), !is.na(pheno), dplyr::if_all(tidyselect::all_of(adjustment_variables), ~!is.na(.)))
  dsn <- create_svydesign(dat)
  ## run models
  baseformula <- stats::as.formula("pheno ~ expo")
  baseadjusted <- addToBase(baseformula, adjustingVariables = adjustment_variables)
  basebase <- stats::as.formula(sprintf("pheno ~ %s", paste(adjustment_variables, collapse="+")))
  demo_break_tbl <- demographic_breakdown(dsn)
  ##
  unadjusted_mod <- run_model(baseformula, dsn, scale_expo = scale_e, scale_pheno = scale_p,  quantile_expo=quantile_expo, expo_levels =  exposure_levels)
  adjusted_mod <- run_model(baseadjusted, dsn, scale_expo = scale_e, scale_pheno = scale_p, quantile_expo=quantile_expo, expo_levels = exposure_levels)
  base_mod <- run_model(basebase, dsn, scale_expo = scale_e, scale_pheno = scale_p, quantile_expo=quantile_expo, expo_levels = exposure_levels)
  n <- dsn |> nrow()
  ## return mods
  list(dsn=dsn, log_p = logxform_p, log_e = logxform_e, scaled_p = scale_p, scaled_e=scale_e, unweighted_n=n, phenotype=pheno, series=tab_obj$series, exposure=exposure, unadjusted_model=unadjusted_mod, adjusted_model=adjusted_mod, base_mod=base_mod, demographic_breakdown=demo_break_tbl)
}


#' PE by Table Flexible Adjustments
#'
#' This function operates similarly to the 'pe' function, but takes a pre-processed table object as input instead of separate phenotype and exposure tables. It is designed to help analysts scale up the pe associations by sepearting the "get_tables" procedure from runnin an association.
#' It also takes a tibble for the adjustment variables to flexibly parameterize adjustments
#' @param tab_obj A pre-processed table object including both phenotype and exposure information
#' @param pvar A variable related to phenotype
#' @param evar A variable related to exposure
#' @param adjustment_variables A tibble of variables to adjust for in the models (to consider multiple scenarios)
#' @param logxform_p Logical, if TRUE, a log transformation is applied to phenotype variable. Default is TRUE.
#' @param logxform_e Logical, if TRUE, a log transformation is applied to exposure variable. Default is TRUE.
#' @param scale_e Logical, if TRUE, exposure variable is scaled. Default is TRUE.
#' @param scale_p Logical, if TRUE, phenotype variable is scaled. Default is FALSE.
#' @param quantile_expo Optional, if specified, exposure variable will be cut by these quantiles
#' @param exposure_levels Optional, if specified, exposure variable will be treated as a categorical variable with these levels
#'
#' @return A list containing: the svydesign object, log transformation status for phenotype and exposure, scaling status for phenotype and exposure, unweighted number of observations, phenotype, series, exposure, models (per adjustment) and demographic breakdown.
#'
#' @examples
#' \dontrun{
#' con <- connect_pewas_data()
#' table_object <- get_tables("LBXGLU", "LBXGTC", "C", con, "L10AM_C", "L45VIT_C")
#' pe_result <- pe_by_table_flex_adjust(tab_obj=table_object, pvar="LBXGLU", evar="LBXGTC",
#' adjustment_variables=c("RIDAGEYR", "RIAGENDR"), # tibble
#' logxform_p=T, logxform_e=T, scale_e=T, scale_p=F,
#' quantile_expo=NULL, exposure_levels=NULL)
#' }
#'
#' @export
pe_by_table_flex_adjust <- function(tab_obj, pvar, evar,
                        adjustment_variables, ## tibble, indexed by scenario and list of adjustment variables
                        logxform_p=T, logxform_e=T, scale_e=T, scale_p=F,
                        quantile_expo=NULL, exposure_levels=NULL) {

  pheno <- pvar
  exposure <- evar
  tab_obj <- name_and_xform_pheno_expo(pheno, exposure, tab_obj, logxform_p, logxform_e)
  ## create svydesign

  potential_adjusters <- setdiff(unique(adjustment_variables$variables), NA)

  dat <- tab_obj$merged_tab |> dplyr::filter(!is.na(wt), wt > 0, !is.na(expo), !is.na(pheno), dplyr::if_all(tidyselect::all_of(potential_adjusters), ~!is.na(.)))
  dsn <- create_svydesign(dat)
  demo_break_tbl <- demographic_breakdown(dsn)

  ## run models
  baseformula <- stats::as.formula("pheno ~ expo")
  uniq_model_scenarios <- unique(adjustment_variables$scenario)
  models <- vector(mode="list", length=length(uniq_model_scenarios))
  logger::log_info("total models { length(uniq_model_scenarios)} " )
  for(mod_num in 1:length(uniq_model_scenarios)) {
    scene <- uniq_model_scenarios[mod_num]
    logger::log_info("running model {scene}" )
    adjust_variables_for_scene <- adjustment_variables |> filter(scenario==scene) |> dplyr::pull(variables)
    baseadjusted <- NA
    if(length(adjust_variables_for_scene) == 1 & is.na(adjust_variables_for_scene[1])) {
      baseadjusted <-baseformula
    } else {
      baseadjusted <- addToBase(baseformula, adjust_variables_for_scene)
    }
    models[[mod_num]] <- run_model(baseadjusted,dsn, scale_expo = scale_e, scale_pheno = scale_p,  quantile_expo=quantile_expo, expo_levels =  exposure_levels)
  }
  n <- dsn |> nrow()
  ## return mods
  list(log_p = logxform_p, log_e = logxform_e, scaled_p = scale_p, scaled_e=scale_e, unweighted_n=n, phenotype=pheno, series=tab_obj$series, exposure=exposure, models=models, adjustment_variables=adjustment_variables, demographic_breakdown=demo_break_tbl)
}


#' Provide a demographic breakdown for a survey design object
#'
#' This function calculates the mean for several hard-coded demographic variables in the survey design object.
#' The demographic variables are: RIDAGEYR, RIAGENDR, INDFMPIR, ETHNICITY_MEXICAN, ETHNICITY_OTHERHISPANIC,
#' ETHNICITY_OTHER, ETHNICITY_NONHISPANICBLACK, ETHNICITY_NONHISPANICWHITE, EDUCATION_LESS9, EDUCATION_9_11,
#' EDUCATION_HSGRAD, EDUCATION_AA, and EDUCATION_COLLEGEGRAD.
#'
#' @param svy_dsn A survey design object dervied from the nhanes_pewas::get_*_tables() functions
#' @return A tibble containing the mean for each demographic variable in the survey design object.
#' @export
demographic_breakdown <- function(svy_dsn) {
  ## for a merged table, provide the summary of the demovars3

  varnames <- adjustment_variables()
  # Check for presence of variables
  dataset_vars <- names(svy_dsn$variables)
  missing_vars <- setdiff(varnames, dataset_vars)

  if(length(missing_vars) > 0) {
    stop(paste("The following variables are missing: ", paste(missing_vars, collapse = ", ")))
  }

  mn_obj <- survey::svymean(~RIDAGEYR+I(RIAGENDR-1)+INDFMPIR+ETHNICITY_MEXICAN+ETHNICITY_OTHERHISPANIC+ETHNICITY_OTHER+ETHNICITY_NONHISPANICBLACK+ETHNICITY_NONHISPANICWHITE+EDUCATION_LESS9+EDUCATION_9_11+EDUCATION_HSGRAD+EDUCATION_AA+EDUCATION_COLLEGEGRAD, design=svy_dsn)
  varnames <- mn_obj |> names()
  varnames[grepl("RIAGENDR", varnames)] <- "RIAGENDR"
  ret <- mn_obj |> tibble::as_tibble() |> dplyr::mutate(varname=varnames)
  ret
}





