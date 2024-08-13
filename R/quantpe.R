## Chirag J Patel
## 07/14/23
## pe.R --> consider renaming to pe.R
## machinery to associate phenotypes with a categorical or real-valued exposure

check_e_data_type <- function(varname, con) {
  ret <- list(vartype="continuous", varlevels=NULL)

  if(grepl('CNT$', varname)) {
    return(list(vartype="continuous-rank", varlevels=NULL))
  }

  if(grepl("^PAQ", varname)) {
    return(list(vartype="continuous", varlevels=NULL))
  }

  elvl <- tbl(con, 'e_variable_levels') |>  filter(Variable.Name == varname, !is.na(values)) |> pull(values) |> unique()

  if(length(elvl) == 1) {
    return(list(vartype="continuous", varlevels=NULL))
  } else if(any(elvl < 1 & elvl > 0) | any(round(elvl) != elvl)) {
    return(list(vartype="continuous-rank", varlevels=sort(elvl)))
  } else if(all(round(elvl) == elvl)) {
    return(list(vartype="categorical", varlevels=sort(elvl)))
  }
  return(ret)
}
xysvydesign <- function(get_tables_obj) {
  create_svydesign(get_tables_obj$merged_tab)
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
    logger::log_info("Phenotype {pheno} being logged")
    table_object$merged_tab <- table_object$merged_tab |> dplyr::mutate(pheno=log10_xform_variable(!!as.name(pheno)))
  } else {
    logger::log_info("Phenotype {pheno} being kept on the natural scale")
    table_object$merged_tab <- table_object$merged_tab |> dplyr::mutate(pheno=!!as.name(pheno))
  }

  if(logxform_e) {
    ##
    logger::log_info("Exposure {exposure} being logged")
    table_object$merged_tab <- table_object$merged_tab |> dplyr::mutate(expo=(log10_xform_variable(!!as.name(exposure))))
  } else {
    logger::log_info("Exposure {exposure} being kept on the natural scale")
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


#' P-E association by a single survey
#'
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
pe_by_survey_series <- function(pheno, exposure, adjustment_variables, series, con,
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


#' Phenotype and Exposure Analysis
#'
#' This function performs analysis of phenotype and exposure with adjustment variables. It handles the retrieval of necessary data, weighting, and running of models, and provides a demographic breakdown.
#'
#' @param pheno A character string specifying the phenotype variable.
#' @param exposure A character string specifying the exposure variable.
#' @param adjustment_variables A character vector of adjustment variables that must be present in the demographic table.
#' @param con A database connection object.
#' @param series A character vector specifying the series to consider (default is NULL).
#' @param logxform_p Logical, whether to log-transform the phenotype variable (default is TRUE).
#' @param logxform_e Logical, whether to log-transform the exposure variable (default is TRUE).
#' @param scale_e Logical, whether to scale the exposure variable (default is TRUE).
#' @param scale_p Logical, whether to scale the phenotype variable (default is FALSE).
#' @param pheno_table_name Optional character string specifying the phenotype table name (default is NULL).
#' @param expo_table_name Optional character string specifying the exposure table name (default is NULL).
#' @param quantile_expo Optional vector specifying quantiles for the exposure variable (default is NULL).
#' @param exposure_levels Optional vector specifying levels for the exposure variable (default is NULL).
#' @return A list containing the following elements:
#' \itemize{
#'   \item `dsn`: The survey design object.
#'   \item `unweighted_n`: The number of unweighted observations.
#'   \item `phenotype`: The phenotype variable.
#'   \item `exposure`: The exposure variable.
#'   \item `series`: The series information.
#'   \item `unadjusted_model`: The unadjusted model.
#'   \item `adjusted_model`: The adjusted model.
#'   \item `base_mod`: The base model adjusted for the adjustment variables.
#'   \item `demo_breakdown`: A table with demographic breakdown of the survey design.
#' }
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item If a single series is provided, calls `pe_by_survey_series` for analysis.
#'   \item Retrieves table names for the phenotype and exposure variables.
#'   \item Ensures the phenotype and exposure variables are collected in the same survey.
#'   \item Gets table names for each series and binds the tables.
#'   \item Handles weights by calling `figure_out_multiyear_weight`.
#'   \item Transforms and names the phenotype and exposure variables if specified.
#'   \item Filters data and creates a survey design object.
#'   \item Performs demographic breakdown of the survey design.
#'   \item Runs unadjusted, adjusted, and base models for the analysis.
#' }
#' @examples
#' \dontrun{
#' con <- DBI::dbConnect(RSQLite::SQLite(), dbname = "nhanes.db")
#' adjustment_vars <- c("AGE", "SEX")
#' result <- pe(pheno = "BMI", exposure = "SMOKING", adjustment_variables = adjustment_vars, con = con)
#' }
#' @export
pe <- function(pheno, exposure, adjustment_variables,con, series=NULL,
               logxform_p=T, logxform_e=T, scale_e=T, scale_p=F,
               pheno_table_name=NULL, expo_table_name=NULL,
               quantile_expo=NULL, exposure_levels=NULL) {
  ## adjustment_variables has to be in the demo table
  if(is.null(series) & length(series) == 1) {
    return(pe_by_survey_series(pheno,exposure, adjustment_variables, con, series, logxform_p, logxform_e, scale_e, scale_p, pheno_table_name, expo_table_name, quantile_expo, exposure_levels))
  }
  ## get tables
  ptables <- get_table_names_for_varname(con, varname = pheno, series) |> rename(p_name = Data.File.Name)
  etables <- get_table_names_for_varname(con, varname = exposure, series) |> rename(e_name = Data.File.Name)
  table_set <- ptables |> inner_join(etables, by = "Begin.Year")
  if(nrow(table_set) == 0) {
    stop("Y and X variables not collected in the same survey")
  }

  ## get table names for each series
  tab_obj <- get_x_y_tables_as_list(con,table_set$p_name,table_set$e_name)

  ## weight
  tab_obj <- figure_out_multiyear_weight(tab_obj)
  ## then merge tables

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



#' Phenotype and Exposure Analysis under different modeling scenarios
#'
#' This function performs flexible adjustment for phenotype and exposure analysis using specified adjustment variables. It retrieves necessary data, handles weighting, and runs models based on the specified scenarios.
#'
#' @param pheno A character string specifying the phenotype variable.
#' @param exposure A character string specifying the exposure variable.
#' @param adjustment_variables A data frame containing adjustment variables and scenarios. It should have columns `variables` and `scenario`.
#' @param con A database connection object.
#' @param series A character vector specifying the series to consider (default is NULL).
#' @param logxform_p Logical, whether to log-transform the phenotype variable (default is TRUE).
#' @param logxform_e Logical, whether to log-transform the exposure variable (default is TRUE).
#' @param scale_e Logical, whether to scale the exposure variable (default is TRUE).
#' @param scale_p Logical, whether to scale the phenotype variable (default is FALSE).
#' @param pheno_table_name Optional character string specifying the phenotype table name (default is NULL).
#' @param expo_table_name Optional character string specifying the exposure table name (default is NULL).
#' @param quantile_expo Optional vector specifying quantiles for the exposure variable (default is NULL).
#' @param exposure_levels Optional vector specifying levels for the exposure variable (default is NULL).
#' @return A list containing the following elements:
#' \itemize{
#'   \item `log_p`: Logical, whether the phenotype was log-transformed.
#'   \item `log_e`: Logical, whether the exposure was log-transformed.
#'   \item `scaled_p`: Logical, whether the phenotype was scaled.
#'   \item `scaled_e`: Logical, whether the exposure was scaled.
#'   \item `unweighted_n`: The number of unweighted observations.
#'   \item `phenotype`: The phenotype variable.
#'   \item `series`: The series information from the data.
#'   \item `exposure`: The exposure variable.
#'   \item `models`: A list of models run based on the scenarios.
#'   \item `base_models`: A list of base models run for comparison.
#'   \item `adjustment_variables`: The adjustment variables used in the models, a tibble: see adjustment_models
#'   \item `demographic_breakdown`: A table with demographic breakdown of the survey design.
#' }
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Retrieves table names for the phenotype and exposure variables.
#'   \item Ensures the phenotype and exposure variables are collected in the same survey.
#'   \item Gets table names for each series and binds the tables.
#'   \item Handles weights by calling `figure_out_multiyear_weight`.
#'   \item Transforms and names the phenotype and exposure variables if specified.
#'   \item Filters data and creates a survey design object.
#'   \item Performs demographic breakdown of the survey design.
#'   \item Runs models based on the specified scenarios and adjustment variables.
#' }
#' @examples
#' \dontrun{
#' con <- DBI::dbConnect(RSQLite::SQLite(), dbname = "nhanes.db")
#' adjustment_vars <- data.frame(variables = c("AGE", "SEX"), scenario = c("A", "A"))
#' result <- pe_flex_adjust(pheno = "BMXBMI", exposure = "LBXCOT", adjustment_variables = adjustment_vars, con = con)
#' }
#' @export
pe_flex_adjust <- function(pheno, exposure, adjustment_variables,con, series=NULL,
                           logxform_p=T, logxform_e=T, scale_e=T, scale_p=F,
                           pheno_table_name=NULL, expo_table_name=NULL,
                           quantile_expo=NULL, exposure_levels=NULL) {

  ptables <- get_table_names_for_varname(con, varname = pheno, series) |> rename(p_name = Data.File.Name)
  etables <- get_table_names_for_varname(con, varname = exposure, series) |> rename(e_name = Data.File.Name)
  table_set <- ptables |> inner_join(etables, by = "Begin.Year")
  if(nrow(table_set) == 0) {
    stop("Y and X variables not collected in the same survey")
  }


  ## get table names for each series
  tab_obj <- get_x_y_tables_as_list(con,table_set$p_name,table_set$e_name)

  ## weight
  tab_obj <- figure_out_multiyear_weight(tab_obj)

  tab_obj <- name_and_xform_pheno_expo(pheno, exposure, tab_obj, logxform_p, logxform_e)

  potential_adjusters <- setdiff(unique(adjustment_variables$variables), NA)
  dat <- tab_obj$merged_tab |> dplyr::filter(!is.na(wt), wt > 0, !is.na(expo), !is.na(pheno), dplyr::if_all(tidyselect::all_of(potential_adjusters), ~!is.na(.)))
  dsn <- create_svydesign(dat)
  demo_break_tbl <- demographic_breakdown(dsn)

  ## run models
  baseformula <- stats::as.formula("pheno ~ expo")
  basebase <- stats::as.formula("pheno ~ 1")
  uniq_model_scenarios <- unique(adjustment_variables$scenario)
  models <- vector(mode="list", length=length(uniq_model_scenarios))
  base_models <- vector(mode="list", length=length(uniq_model_scenarios))
  logger::log_info("total models { length(uniq_model_scenarios)} " )
  for(mod_num in 1:length(uniq_model_scenarios)) {
    scene <- uniq_model_scenarios[mod_num]
    logger::log_info("running model {scene}" )
    adjust_variables_for_scene <- adjustment_variables |> filter(scenario==scene) |> dplyr::pull(variables)
    baseadjusted <- NA
    if(length(adjust_variables_for_scene) == 1 & is.na(adjust_variables_for_scene[1])) {
      baseadjusted <- baseformula
    } else {
      baseadjusted <- addToBase(baseformula, adjust_variables_for_scene)
      basebaseadjusted <- addToBase(basebase, adjust_variables_for_scene)
      base_models[[mod_num]] <- run_model(basebaseadjusted,dsn, scale_expo = scale_e, scale_pheno = scale_p,  quantile_expo=quantile_expo, expo_levels =  exposure_levels)
    }
    models[[mod_num]] <- run_model(baseadjusted,dsn, scale_expo = scale_e, scale_pheno = scale_p,  quantile_expo=quantile_expo, expo_levels =  exposure_levels)
  }
  n <- dsn |> nrow()
  list(log_p = logxform_p, log_e = logxform_e, scaled_p = scale_p, scaled_e=scale_e, unweighted_n=n, phenotype=pheno, series=tab_obj$series, exposure=exposure, models=models, base_models=base_models, adjustment_variables=adjustment_variables, demographic_breakdown=demo_break_tbl)
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
                                    logxform_p=F, logxform_e=T, scale_e=T, scale_p=T,
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
  basebase <- stats::as.formula("pheno ~ 1")
  uniq_model_scenarios <- unique(adjustment_variables$scenario)
  models <- vector(mode="list", length=length(uniq_model_scenarios))
  base_models <- vector(mode="list", length=length(uniq_model_scenarios))
  logger::log_info("total models { length(uniq_model_scenarios)} " )
  for(mod_num in 1:length(uniq_model_scenarios)) {
    scene <- uniq_model_scenarios[mod_num]
    logger::log_info("running model {scene}" )
    adjust_variables_for_scene <- adjustment_variables |> filter(scenario==scene) |> dplyr::pull(variables)
    baseadjusted <- NA
    if(length(adjust_variables_for_scene) == 1 & is.na(adjust_variables_for_scene[1])) {
      baseadjusted <- baseformula
    } else {
      baseadjusted <- addToBase(baseformula, adjust_variables_for_scene)
      basebaseadjusted <- addToBase(basebase, adjust_variables_for_scene)
      base_models[[mod_num]] <- run_model(basebaseadjusted,dsn, scale_expo = scale_e, scale_pheno = scale_p,  quantile_expo=quantile_expo, expo_levels =  exposure_levels)
    }
    models[[mod_num]] <- run_model(baseadjusted,dsn, scale_expo = scale_e, scale_pheno = scale_p,  quantile_expo=quantile_expo, expo_levels =  exposure_levels)
  }
  n <- dsn |> nrow()
  list(log_p = logxform_p, log_e = logxform_e, scaled_p = scale_p, scaled_e=scale_e, unweighted_n=n, phenotype=pheno, series=tab_obj$series, exposure=exposure, models=models, base_models=base_models, adjustment_variables=adjustment_variables, demographic_breakdown=demo_break_tbl)
}


#' Provide a demographic breakdown for a survey design object
#'
#' This function calculates the mean for several hard-coded demographic variables in the survey design object.
#' The demographic variables are: RIDAGEYR, RIAGENDR, INDFMPIR, ETHNICITY_MEXICAN, ETHNICITY_OTHERHISPANIC,
#' ETHNICITY_OTHER, ETHNICITY_NONHISPANICBLACK, ETHNICITY_NONHISPANICWHITE, EDUCATION_LESS9, EDUCATION_9_11,
#' EDUCATION_HSGRAD, EDUCATION_AA, and EDUCATION_COLLEGEGRAD.
#'
#' @param svy_dsn A survey design object derived from the nhanespewas::get_*_tables() functions
#' @return A tibble containing the mean for each demographic variable in the survey design object.
#' @export
demographic_breakdown <- function(svy_dsn) {
  ## for a merged table, provide the summary of the demovars
  demographic_variables <- function() {
    varnames <- c("RIDAGEYR", "RIAGENDR", "INDFMPIR", "ETHNICITY_MEXICAN",
                  "ETHNICITY_OTHERHISPANIC", "ETHNICITY_OTHER",
                  "ETHNICITY_NONHISPANICBLACK", "ETHNICITY_NONHISPANICWHITE",
                  "EDUCATION_LESS9", "EDUCATION_9_11", "EDUCATION_HSGRAD",
                  "EDUCATION_AA", "EDUCATION_COLLEGEGRAD")
  }

  varnames <- demographic_variables()
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




#' XY by Table Flexible Adjustments
#'
#' This function operates similarly to the 'pe' function, but takes a pre-processed table object as input instead of separate phenotype and exposure tables. It is designed to help analysts scale up the pe associations by sepearting the "get_tables" procedure from runnin an association.
#' It also takes a tibble for the adjustment variables to flexibly parameterize adjustments
#' @param tab_obj A pre-processed table object including both phenotype and exposure information
#' @param yvar A variable related to phenotype
#' @param xvar A variable related to exposure
#' @param adjustment_variables A tibble of variables to adjust for in the models (to consider multiple scenarios)
#' @param logxform_y Logical, if TRUE, a log transformation is applied to phenotype variable. Default is TRUE.
#' @param logxform_x Logical, if TRUE, a log transformation is applied to exposure variable. Default is TRUE.
#' @param scale_y Logical, if TRUE, exposure variable is scaled. Default is TRUE.
#' @param scale_y Logical, if TRUE, phenotype variable is scaled. Default is FALSE.

#' @return A list containing: the svydesign object, log transformation status for phenotype and exposure, scaling status for phenotype and exposure, unweighted number of observations, phenotype, series, exposure, models (per adjustment) and demographic breakdown.
#'
#' @examples
#' \dontrun{
#' con <- connect_pewas_data()
#' table_object <- get_tables("LBXGLU", "LBXGTC", "C", con, "L10AM_C", "L45VIT_C")
#' pe_result <- xy_by_table_flex_adjust(tab_obj=table_object, yvar="LBXGLU", xvar="LBXGTC",
#' adjustment_variables=c("RIDAGEYR", "RIAGENDR"), # tibble
#' }
#'
#' @export

xy_by_table_flex_adjust <- function(tab_obj, yvar, xvar,
                                    adjustment_variables, ## tibble, indexed by scenario and list of adjustment variables
                                    logxform_x=T, logxform_y=T, scale_x=T, scale_y=T,
                                    quantile_expo=NULL, exposure_levels=NULL) {

  pheno <- yvar
  exposure <- xvar
  tab_obj <- name_and_xform_pheno_expo(pheno, exposure, tab_obj, logxform_y, logxform_x)
  ## create svydesign
  potential_adjusters <- setdiff(unique(adjustment_variables$variables), NA)

  dat <- tab_obj$merged_tab |> dplyr::filter(!is.na(wt), wt > 0, !is.na(expo), !is.na(pheno), dplyr::if_all(tidyselect::all_of(potential_adjusters), ~!is.na(.)))
  dsn <- create_svydesign(dat)
  demo_break_tbl <- demographic_breakdown(dsn)

  ## run models
  baseformula <- stats::as.formula("pheno ~ expo")
  basebase <- stats::as.formula("pheno ~ 1")
  uniq_model_scenarios <- unique(adjustment_variables$scenario)
  models <- vector(mode="list", length=length(uniq_model_scenarios))
  base_models <- vector(mode="list", length=length(uniq_model_scenarios))
  logger::log_info("total models { length(uniq_model_scenarios)} " )
  for(mod_num in 1:length(uniq_model_scenarios)) {
    scene <- uniq_model_scenarios[mod_num]
    logger::log_info("running model {scene}" )
    adjust_variables_for_scene <- adjustment_variables |> filter(scenario==scene) |> dplyr::pull(variables)
    baseadjusted <- NA
    if(length(adjust_variables_for_scene) == 1 & is.na(adjust_variables_for_scene[1])) {
      baseadjusted <- baseformula
    } else {
      baseadjusted <- addToBase(baseformula, adjust_variables_for_scene)
      basebaseadjusted <- addToBase(basebase, adjust_variables_for_scene)
      base_models[[mod_num]] <- run_model(basebaseadjusted,dsn, scale_expo = scale_x, scale_pheno = scale_y,quantile_expo=quantile_expo, expo_levels =  exposure_levels)
    }
    models[[mod_num]] <- run_model(baseadjusted,dsn, scale_expo = scale_x, scale_pheno = scale_y, quantile_expo=quantile_expo, expo_levels =  exposure_levels)
  }
  n <- dsn |> nrow()
  list(log_y = logxform_y, log_x = logxform_x, scaled_y = scale_y, scaled_x=scale_x, unweighted_n=n, yvar=yvar, series=tab_obj$series, xvar=xvar, models=models, base_models=base_models, adjustment_variables=adjustment_variables, demographic_breakdown=demo_break_tbl)
}



check_e_data_type <- function(varname, con) {
  ret <- list(vartype="continuous", varlevels=NULL)

  if(grepl('CNT$', varname)) {
    return(list(vartype="continuous-rank", varlevels=NULL))
  }

  if(grepl("^PAQ", varname)) {
    return(list(vartype="continuous", varlevels=NULL))
  }

  elvl <- tbl(con, 'e_variable_levels') |>  filter(Variable.Name == varname, !is.na(values)) |> pull(values) |> unique()

  if(length(elvl) == 1) {
    return(list(vartype="continuous", varlevels=NULL))
  } else if(any(elvl < 1 & elvl > 0) | any(round(elvl) != elvl)) {
    return(list(vartype="continuous-rank", varlevels=sort(elvl)))
  } else if(all(round(elvl) == elvl)) {
    return(list(vartype="categorical", varlevels=sort(elvl)))
  }
  return(ret)
}
