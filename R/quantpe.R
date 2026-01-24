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

# considers interaction
aceOfBase <- function(base_formula, adjustingVariables, interact_with = NULL, interactor_variable="expo") {
  form <- addToBase(base_formula, adjustingVariables)
  if (!is.null(interact_with) && length(interact_with) > 0) {
    missing_vars <- setdiff(interact_with, adjustingVariables)
    if (length(missing_vars) > 0) {
      stop(sprintf("All interaction variables must be in adjustingVariables. Missing: %s",
                   paste(missing_vars, collapse = ", ")))
    }
    interaction_terms <- paste(sprintf("expo*%s", interact_with), collapse = " + ")
    form <- stats::update.formula(form, paste("~ . +", interaction_terms))
  }
  return(form)
}



#' Calculate R-squared and adjusted R-squared for a survey-weighted linear model
#'
#' This function calculates the R-squared and adjusted R-squared values for a
#' linear model based on an analysis object which containssave   a survey design,
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


svyregtermtest_for_formula <- function(mod,regTermTestForm) {
  # does the default for svyregression, a tidy wrapper for survey::regTermTest
  regTest <-survey::regTermTest(mod, regTermTestForm)
  tibble(
    p.value=as.numeric(regTest$p),
    f.test = as.numeric(regTest$Ftest),
    df=as.numeric(regTest$df),
    ddf = as.numeric(regTest$ddf),
    test.terms = paste( regTest$test.terms, collapse="+")
  )
}

#' Fit a survey‐weighted GLM with optional scaling, categorization, and joint term testing
#'
#' This function fits a survey‐weighted generalized linear model (`svyglm`) using the
#' provided formula and design object, with flexible options to:
#' - scale or quantile‐categorize the exposure variable,
#' - scale the outcome (phenotype) variable via mean/SD, CLR, or inverse‐variance transform,
#' - run a joint term test (Wald test) on a specified set of terms,
#' - optionally return the full model object.
#'
#' @param formu
#'   A two‐sided [`formula`] describing the regression model, e.g. `pheno ~ expo + covariate1`.
#' @param dsn
#'   A [`survey.design`] object created by `survey::svydesign()` or similar.
#' @param scale_expo
#'   Logical (default `TRUE`): if `TRUE`, center and scale `expo` to mean 0, SD 1.
#'   Ignored if `quantile_expo` or `expo_levels` is non‐`NULL`.
#' @param scale_pheno
#'   Logical (default `FALSE`): if `TRUE`, scale the outcome variable (`pheno`).
#'   The method is chosen by `scale_type`.
#' @param quantile_expo
#'   Optional numeric vector of quantiles (e.g. `c(0.25, 0.5, 0.75)`).
#'   If not `NULL`, `expo` is converted to a factor with breaks at those survey‐weighted quantiles.
#' @param expo_levels
#'   Optional character or factor vector of levels for `expo`.
#'   If not `NULL`, `expo` is converted to a factor with those levels.
#'   Overrides `quantile_expo` if both are provided.
#' @param scale_type
#'   Integer (1, 2, or 3; default `1`) selecting the outcome scaling method when `scale_pheno = TRUE`:
#'   \itemize{
#'     \item `1` = standard (center & scale by mean/SD),
#'     \item `2` = centered log‐ratio (CLR),
#'     \item `3` = inverse‐variance (rank‐normal) transform.
#'   }
#' @param regTermTestForm
#'   Optional one‐sided [`formula`] (e.g. `~ expo + expo:age`) specifying terms to jointly test
#'   via `survey::regTermTest()`.  If provided, a tidy summary of that test is returned.
#' @param save_svymodel
#'   Logical (default `FALSE`): if `TRUE`, the full `svyglm` object is returned in the output list
#'   under `$model`; otherwise `$model` is set to `NA`.
#'
#' @return
#' A named list with components:
#' \describe{
#'   \item{glanced}{A one‐row tibble from `broom::glance()` of model‐level fit statistics.}
#'   \item{tidied}{A tibble from `broom::tidy()` of term‐level coefficients and tests.}
#'   \item{r2}{A tibble with R² measures from `svyrsquared()`.}
#'   \item{scale_pheno}{Logical, echoing the `scale_pheno` input.}
#'   \item{scale_expo}{Logical, echoing the `scale_expo` input.}
#'   \item{q_cut_points}{Numeric vector of cut‐points if `quantile_expo` was used; otherwise `NA`.}
#'   \item{quantiles}{Echo of `quantile_expo`.}
#'   \item{expo_levels}{Echo of `expo_levels`.}
#'   \item{reg_term_tidied}{If `regTermTestForm` provided, a one‐row tibble with joint‐term test
#'     results (`p.value`, `f.test`, `df`, `ddf`, `test.terms`); otherwise `NULL`.}
#'   \item{model}{The full `svyglm` object if `save_svymodel = TRUE`, else `NA`.}
#' }
#'
#' @examples
#' \dontrun{
#' library(survey)
#'
#' # Prepare survey design `dsn`...
#' form <- pheno ~ expo + age + sex
#' result <- run_model(
#'   formu            = form,
#'   dsn              = my_design,
#'   scale_expo       = TRUE,
#'   scale_pheno      = TRUE,
#'   quantile_expo    = c(0.25, 0.5, 0.75),
#'   expo_levels      = NULL,
#'   scale_type       = 1,
#'   regTermTestForm  = ~ expo + expo:age,
#'   save_svymodel    = FALSE
#' )
#'
#' # Access tidied coefficients:
#' result$tidied
#' # Joint Wald test of expo terms:
#' result$reg_term_tidied
#' }
#'
#' @export
run_model <- function(formu, dsn, scale_expo=T, scale_pheno=F, quantile_expo=NULL, expo_levels=NULL,  scale_type=1, regTermTestForm=NULL, save_svymodel=F) { # scale type 1 is mean and SD; scale type 2 == CLR scale type 3 is IVT

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

  if(scale_pheno & scale_type == 1) {
    logger::log_info("Phenotype scaled by mean and SD")
    s <- sqrt(mean(survey::svyvar(~pheno, dsn, na.rm = T)))
    mn <- mean(survey::svymean(~pheno, dsn, na.rm=T))
    dsn <- stats::update(dsn, pheno=(pheno-mn)/s)
  } else if(scale_pheno & scale_type == 2) {
    # centered log ratio transformation
    logger::log_info("Phenotype scaled using CLR")
    dsn <- stats::update(dsn, logpheno=log(pheno+.25))
    mn <- mean(survey::svymean(~logpheno, dsn, na.rm=T))
    dsn <- stats::update(dsn, pheno=logpheno-mn)
  } else if(scale_pheno & scale_type == 3) {
    logger::log_info("Phenotype scaled using RankNorm")
    ## inverse variance transform
    #dsn <- RNOmni::rankNormal
    pheno_raw <- dsn$variables$pheno
    # Create an empty vector (with NA values) of the same length
    pheno_trans <- rep(NA, length(pheno_raw))
    # Apply the Rank Normal transformation to non-missing values only
    non_missing <- !is.na(pheno_raw)
    pheno_trans[non_missing] <- RNOmni::RankNorm(pheno_raw[non_missing])
    # Update the survey design object with the transformed phenotype
    dsn <- stats::update(dsn, pheno = pheno_trans)
  }


  mod <- survey::svyglm(formu,dsn)
  ti <- broom::tidy(mod)
  gl <- broom::glance(mod)
  r2 <- svyrsquared(mod)
  regTermTest <- NULL
  if(!is.null(regTermTestForm)) {
    regTermTest <- svyregtermtest_for_formula(mod, regTermTestForm)
  }

  list(glanced=gl, tidied=ti, r2=r2,
       scale_pheno=scale_pheno, scale_expo=scale_expo,
       q_cut_points=cut_points, quantiles=quantile_expo, expo_levels=expo_levels,
       reg_term_tidied=regTermTest, model= if (save_svymodel) mod else NA) # not saving dsn
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
    logger::log_info("Phenotype {pheno} being kept on the natural scale (not being logged)")
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


#' Phenotype and Exposure Association Across Surveys
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
#' adjustment_vars <- c("RIDAGEYR", "RIAGENDR")
#' result <- pe(pheno = "BMXBMI", exposure = "LBXCOT", adjustment_variables = adjustment_vars, con = con)
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



#' Phenotype and Exposure Flexible Adjustment and Interaction Modeling
#'
#' This function performs phenotype–exposure analysis under multiple modeling scenarios,
#' including optional log-transformation, scaling, and interaction terms. It retrieves
#' data from the database, merges phenotype and exposure tables, handles survey weights,
#' applies transformations, filters missingness, constructs a survey design, computes
#' demographic summaries, and fits models for each scenario defined in
#' `adjustment_variables`.
#'
#' @param pheno
#'   A string naming the phenotype (outcome) variable.
#' @param exposure
#'   A string naming the exposure (predictor) variable.
#' @param adjustment_variables
#'   A data frame with columns `variables` (character) and `scenario` (character),
#'   defining which covariates to adjust for in each scenario.
#' @param con
#'   A DBI database connection object used to look up and fetch tables.
#' @param series
#'   Optional character vector of survey series to restrict table lookup;
#'   `NULL` (default) means “use all available series.”
#' @param logxform_p
#'   Logical; if `TRUE` (default), the phenotype is log-transformed before modeling.
#' @param logxform_e
#'   Logical; if `TRUE` (default), the exposure is log-transformed before modeling.
#' @param scale_e
#'   Logical; if `TRUE` (default), the exposure is scaled to mean 0 and SD 1.
#' @param scale_p
#'   Logical; if `TRUE` (default `FALSE`), the phenotype is scaled to mean 0 and SD 1.
#' @param pheno_table_name
#'   Optional string to override the default phenotype table lookup (default `NULL`).
#' @param expo_table_name
#'   Optional string to override the default exposure table lookup (default `NULL`).
#' @param quantile_expo
#'   Optional numeric vector of quantiles (e.g. `c(0.1,0.5,0.9)`) at which to evaluate
#'   exposure effects in `run_model()` (default `NULL`).
#' @param exposure_levels
#'   Optional vector of factor levels for the exposure variable in categorical analyses
#'   (default `NULL`).
#' @param scale_type
#'   Integer code for scaling method (default `1`):
#'   \itemize{
#'     \item `1` = standard (mean/SD),
#'     \item `2` = centered log‐ratio (CLR),
#'     \item `3` = inverse‐variance transformation (IVT).
#'   }
#' @param interact_with
#'   Optional character vector of adjustment variables to include interaction terms
#'   with `expo` (must be a subset of `adjustment_variables$variables`; default `NULL`).
#' @param ...
#'   Additional arguments passed on to [`run_model()`], e.g. model‐specific options.
#'
#' @return
#' A named list with elements:
#' \describe{
#'   \item{log_p}{Logical indicating whether the phenotype was log-transformed.}
#'   \item{log_e}{Logical indicating whether the exposure was log-transformed.}
#'   \item{scaled_p}{Logical indicating whether the phenotype was scaled.}
#'   \item{scaled_e}{Logical indicating whether the exposure was scaled.}
#'   \item{unweighted_n}{Integer: count of unweighted observations in the survey design.}
#'   \item{phenotype}{Character: name of the phenotype variable.}
#'   \item{series}{Data frame of series metadata used for table lookup.}
#'   \item{exposure}{Character: name of the exposure variable.}
#'   \item{models}{List of fitted survey‐weighted models for each scenario (with optional interactions).}
#'   \item{base_models}{List of baseline (intercept‐only or base‐adjusted) models for comparison.}
#'   \item{adjustment_variables}{The input data frame of adjustment variables and scenarios.}
#'   \item{demographic_breakdown}{Data frame of demographic summaries from the survey design.}
#' }
#'
#' @details
#' For each unique `scenario` in `adjustment_variables`, the function:
#' \enumerate{
#'   \item Merges phenotype (`pheno`) and exposure (`exposure`) tables by year.
#'   \item Determines survey weights via `figure_out_multiyear_weight()`.
#'   \item Applies log‐transform and/or scaling to both `pheno` and `exposure` if requested.
#'   \item Filters out rows with missing weights or missing values in the outcome, exposure,
#'         or any adjustment covariates.
#'   \item Constructs a survey design object and computes a demographic breakdown.
#'   \item Uses `aceOfBase()` to build the formula—adding covariates and optional
#'         `expo:covariate` interactions—and then fits models via `run_model()`.
#' }
#'
#' @examples
#' \dontrun{
#' con <- DBI::dbConnect(RSQLite::SQLite(), dbname = "nhanes.db")
#'
#' adj_vars <- data.frame(
#'   variables = c("RIDAGEYR", "RIAGENDR", "INDFMPIR"),
#'   scenario  = c("S1",      "S1",      "S1")
#' )
#'
#' result <- pe_flex_adjust(
#'   pheno                = "LBXCRP",
#'   exposure             = "LBXCOT",
#'   adjustment_variables = adj_vars,
#'   con                  = con,
#'   series               = c("2013-2014", "2015-2016"),
#'   logxform_p           = TRUE,
#'   logxform_e           = TRUE,
#'   scale_e              = TRUE,
#'   scale_p              = FALSE,
#'   quantile_expo        = c(0.1, 0.5, 0.9),
#'   exposure_levels      = c("low", "med", "high"),
#'   scale_type           = 2,
#'   interact_with        = "RIDAGEYR"
#' )
#' }
#'
#' @export
pe_flex_adjust <- function(pheno, exposure, adjustment_variables,con, series=NULL,
                           logxform_p=T, logxform_e=T, scale_e=T, scale_p=F,
                           pheno_table_name=NULL, expo_table_name=NULL,
                           quantile_expo=NULL, exposure_levels=NULL, scale_type=1, interact_with=NULL, ...) {

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

  potential_adjusters <- all.vars(as.formula(sprintf("~%s", paste(potential_adjusters, collapse="+"))))

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
    interact_with_local_model <- interact_with
    regtermformula <- NULL
    if(length(adjust_variables_for_scene) == 1 & is.na(adjust_variables_for_scene[1])) {
      baseadjusted <- baseformula
    } else {
      baseadjusted <- addToBase(baseformula, adjust_variables_for_scene)
      basebaseadjusted <- addToBase(basebase, adjust_variables_for_scene)
      if(!is.null(interact_with)) {
        interact_with_local_model <- intersect(adjust_variables_for_scene, interact_with)
        if(length(interact_with_local_model)>0) {
          term_labels      <- c("expo", paste0("expo:", interact_with_local_model))
          regtermformula <- stats::as.formula(paste("~", paste(term_labels, collapse = " + ")))
          baseadjusted <- aceOfBase(baseformula, adjust_variables_for_scene, interact_with=interact_with_local_model)
          basebaseadjusted <- aceOfBase(basebase, adjust_variables_for_scene)
        }
      }

      base_models[[mod_num]] <- run_model(basebaseadjusted,dsn, scale_expo = scale_e, scale_pheno = scale_p,  quantile_expo=quantile_expo, expo_levels =  exposure_levels, scale_type=scale_type, ...)
    }
    models[[mod_num]] <- run_model(baseadjusted,dsn, scale_expo = scale_e, scale_pheno = scale_p,  quantile_expo=quantile_expo, expo_levels =  exposure_levels,scale_type=scale_type, regTermTestForm=regtermformula, ... )
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

  run_model_safely <- safely(run_model)
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
      base_models[[mod_num]] <- run_model_safely(basebaseadjusted,dsn, scale_expo = scale_e, scale_pheno = scale_p,  quantile_expo=quantile_expo, expo_levels =  exposure_levels)
    }

    models[[mod_num]] <- run_model_safely(baseadjusted,dsn, scale_expo = scale_e, scale_pheno = scale_p,  quantile_expo=quantile_expo, expo_levels =  exposure_levels)
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




#' Estimate exposure ~ demographics (survey-weighted) for a single exposure
#'
#' This is a lightweight helper that reuses your existing “guts”:
#'   get_tables() / figure_out_*weight() / name_and_xform_* / create_svydesign() / run_model()
#'
#' It fits: expo ~ <demographic covariates> (optionally including survey year)
#' and returns a tidy coefficient table plus model fit stats.
#'
#' Notes:
#' - The “outcome” in run_model() is `pheno`, so we map `pheno := expo` and set `expo := 1`
#'   (intercept-only predictor) to make the model estimate E[exposure | demographics].
#' - We force `scale_expo = FALSE` because the predictor is constant (1). Scaling the predictor
#'   would break (svyvar of a constant is 0). The outcome scaling can be toggled via `scale_exposure`.
#'
#' @param exposure Character. NHANES exposure variable name (e.g., "LBXGTC").
#' @param con DBI connection to NHANES sqlite.
#' @param series Optional character vector or NULL. If provided and length==1, uses pe_by_survey_series path.
#' @param adjustment_variables Character vector of demographic covariates.
#' @param logxform_exposure Logical. If TRUE, log10-transform exposure using your log10_xform_variable().
#' @param scale_exposure Logical. If TRUE, scale the exposure outcome (pheno) in run_model().
#' @param include_survey_year Logical. If TRUE, add as.factor(SDDSRVYR) to adjustment_variables (when present).
#' @param expo_table_name Optional override for exposure table name.
#' @param save_svymodel Logical. If TRUE, return the svyglm object from run_model().
#'
#' @return A list with: dsn, unweighted_n, exposure, formula, glanced, tidied, r2, model (optional)
#' @export
expo_on_demographics <- function(exposure,
                                 con,
                                 series = NULL,
                                 adjustment_variables = c(
                                   "RIDAGEYR", "AGE_SQUARED", "RIAGENDR", "INDFMPIR",
                                   "EDUCATION_LESS9", "EDUCATION_9_11", "EDUCATION_AA", "EDUCATION_COLLEGEGRAD",
                                   "ETHNICITY_MEXICAN", "ETHNICITY_OTHERHISPANIC", "ETHNICITY_OTHER", "ETHNICITY_NONHISPANICBLACK"
                                 ),
                                 logxform_exposure = TRUE,
                                 scale_exposure = FALSE,
                                 include_survey_year = FALSE,
                                 expo_table_name = NULL,
                                 save_svymodel = FALSE) {

  # Optionally add survey year as a factor term (only if column exists later)
  if (isTRUE(include_survey_year)) {
    adjustment_variables <- c(adjustment_variables, "as.factor(SDDSRVYR)")
  }

  # ---------- Pull and weight data ----------
  if (!is.null(series) && length(series) == 1) {
    # single-series path using get_tables()
    tab_obj <- get_tables(pheno = exposure, exposure = exposure, series = series, con = con,
                          pheno_table_name = expo_table_name, expo_table_name = expo_table_name)
    tab_obj <- figure_out_weight(tab_obj)
  } else {
    # multi-series path using get_x_y_tables_as_list()
    etables <- get_table_names_for_varname(con, varname = exposure, series) |> rename(e_name = Data.File.Name)
    if (nrow(etables) == 0) stop("Exposure variable not found in get_table_names_for_varname(): ", exposure)

    tab_obj <- get_x_y_tables_as_list(con, etables$e_name, etables$e_name)
    tab_obj <- figure_out_multiyear_weight(tab_obj)
  }

  # ---------- Create `expo` column (the exposure), then map to pheno ----------
  # Use your log10 transform helper (handles zeros)
  if (isTRUE(logxform_exposure)) {
    tab_obj$merged_tab <- tab_obj$merged_tab |>
      dplyr::mutate(expo = log10_xform_variable(!!as.name(exposure)))
  } else {
    tab_obj$merged_tab <- tab_obj$merged_tab |>
      dplyr::mutate(expo = !!as.name(exposure))
  }

  # Invert roles so that run_model() fits: pheno ~ demographics (+ intercept predictor)
  tab_obj$merged_tab <- tab_obj$merged_tab |>
    dplyr::mutate(pheno = expo, expo = 1)

  # Ensure adjustment variables exist (strip function wrappers for checking)
  # e.g. "as.factor(SDDSRVYR)" -> "SDDSRVYR"
  adjust_raw_vars <- all.vars(stats::as.formula(sprintf("~ %s", paste(adjustment_variables, collapse = " + "))))

  dat <- tab_obj$merged_tab |>
    dplyr::filter(
      !is.na(wt), wt > 0,
      !is.na(pheno),
      dplyr::if_all(tidyselect::all_of(adjust_raw_vars), ~ !is.na(.))
    )

  if (nrow(dat) == 0) stop("No complete-case rows after filtering (wt/pheno/adjusters).")

  dsn <- create_svydesign(dat)

  # ---------- Fit model: pheno ~ demographics ----------
  rhs <- paste(adjustment_variables, collapse = " + ")
  form <- stats::as.formula(sprintf("pheno ~ %s", rhs))

  # IMPORTANT: scale_expo must be FALSE because predictor expo==1 is constant
  fit <- run_model(
    formu = form,
    dsn = dsn,
    scale_expo = FALSE,
    scale_pheno = isTRUE(scale_exposure),
    quantile_expo = NULL,
    expo_levels = NULL,
    scale_type = 1,
    regTermTestForm = NULL,
    save_svymodel = save_svymodel
  )

  list(
    dsn = dsn,
    unweighted_n = nrow(dsn),
    exposure = exposure,
    logxform_exposure = logxform_exposure,
    scale_exposure = scale_exposure,
    series = tab_obj$series,
    formula = form,
    glanced = fit$glanced,
    tidied = fit$tidied,
    r2 = fit$r2,
    model = fit$model
  )
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
