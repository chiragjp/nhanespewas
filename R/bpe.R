## Chirag J Patel
## 12/24/24
## bpe.R
## machinery to binary phenotypes with a categorical or real-valued exposure


#' Prepare and Transform Phenotype and Exposure for Logistic Modeling
#'
#' This function assigns the phenotype and exposure variables in a table object and optionally applies a log10 transformation to the exposure variable.
#'
#' @param pheno Character; the name of the phenotype variable.
#' @param exposure Character; the name of the exposure variable.
#' @param table_object A list containing a data frame (`merged_tab`) where the phenotype and exposure variables are present.
#' @param logxform_e Logical; if `TRUE`, the exposure variable is log10-transformed. Default is `TRUE`.
#'
#' @return A modified version of the input `table_object` with the following updates:
#' \item{merged_tab}{The data frame updated with `pheno` and `expo` columns. The `expo` column will be log10-transformed if `logxform_e = TRUE`.}
#'
#' @details
#' The function modifies the input `table_object` by creating a `pheno` column corresponding to the phenotype variable and an `expo` column for the exposure variable. If `logxform_e` is set to `TRUE`, the exposure variable is log10-transformed; otherwise, it is kept on its natural scale. Logging statements provide information on the transformation applied.
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#'
#' # Example table object
#' table_object <- list(
#'   merged_tab = data.frame(
#'     diabetes_status = c(0, 1, 1, 0),
#'     pollution_exposure = c(10, 20, 15, 30)
#'   )
#' )
#'
#' # Assign phenotype and log-transform exposure
#' result <- name_and_xform_logistic_pheno_expo(
#'   pheno = "diabetes_status",
#'   exposure = "pollution_exposure",
#'   table_object = table_object,
#'   logxform_e = TRUE
#' )
#'
#' # View transformed data
#' result$merged_tab
#' }
#'
#' @importFrom dplyr mutate
#' @importFrom logger log_info
#' @export
name_and_xform_logistic_pheno_expo <- function(pheno, exposure, table_object, logxform_e=T) {

  table_object$merged_tab <- table_object$merged_tab |> dplyr::mutate(pheno=!!as.name(pheno))

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


#' Flexible Logistic Regression for Exposure and Adjustment Scenarios
#'
#' This function fits logistic regression models for a phenotype and exposure, allowing flexible adjustment for various scenarios, transformations, and scaling. It integrates survey design handling and demographic breakdowns.
#'
#' @param pheno Character; the name of the phenotype variable.
#' @param exposure Character; the name of the exposure variable.
#' @param adjustment_variables Data frame with two columns: `scenario` (indicating adjustment scenarios) and `variables` (adjustment variables for each scenario).
#' @param con A database connection object to query tables for variables.
#' @param series Character vector (optional); the series of datasets to consider (e.g., survey years).
#' @param logxform_e Logical; if `TRUE`, the exposure variable is log-transformed. Default is `TRUE`.
#' @param scale_e Logical; if `TRUE`, the exposure variable is scaled (z-score transformation). Default is `TRUE`.
#' @param pheno_table_name Character (optional); explicit table name for the phenotype variable. Default is `NULL`.
#' @param expo_table_name Character (optional); explicit table name for the exposure variable. Default is `NULL`.
#' @param quantile_expo Numeric vector of quantiles for discretizing the exposure variable. Default is `NULL`.
#' @param exposure_levels Character vector specifying levels for the exposure variable as a categorical variable. Default is `NULL`.
#'
#' @return A list containing:
#' \item{log_e}{Logical; whether the exposure variable was log-transformed.}
#' \item{scaled_e}{Logical; whether the exposure variable was scaled.}
#' \item{unweighted_n}{Integer; the number of unweighted observations used in the analysis.}
#' \item{phenotype}{The phenotype variable name.}
#' \item{series}{The series of datasets considered.}
#' \item{exposure}{The exposure variable name.}
#' \item{models}{A list of logistic regression models for each adjustment scenario.}
#' \item{base_models}{A list of base models (without exposure) for each adjustment scenario.}
#' \item{adjustment_variables}{The adjustment variables used in the models.}
#' \item{demographic_breakdown}{A demographic breakdown table of the survey design.}
#'
#' @details
#' This function facilitates logistic regression modeling across multiple scenarios for exposure and adjustment variables. It integrates:
#' - Survey design creation and multi-year weights.
#' - Flexible handling of phenotype and exposure transformations (log-transformation, scaling, quantile discretization, or categorization).
#' - Adjustment variable scenarios defined in a structured data frame.
#'
#' The function uses the **survey** package for survey design and regression modeling, and it depends on utility functions like `get_table_names_for_varname`, `get_x_y_tables_as_list`, and `create_svydesign` for data preparation.
#'
#' @examples
#' \dontrun{
#' # Example: Database connection (mock example)
#' con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
#'
#' # Example adjustment variables
#' adjustment_vars <- data.frame(
#'   scenario = c("scenario1"),
#'   variables = c("RIAGENDR")
#' )
#'
#' # Run flexible logistic regression
#' result <- logistic_e_flex_adjust(
#'   pheno = "diabetes",
#'   exposure = "air_pollution",
#'   adjustment_variables = adjustment_vars,
#'   con = con,
#'   series = c("E"),
#'   logxform_e = TRUE,
#'   scale_e = TRUE
#' )
#'
#' # Access demographic breakdown
#' result$demographic_breakdown
#'
#' # Access models
#' result$models
#' }
#'
#' @importFrom dplyr filter if_all rename setdiff inner_join pull
#' @importFrom stats as.formula
#' @importFrom logger log_info
#' @importFrom survey svyquantile svyglm
#' @export
logistic_e_flex_adjust <- function(pheno, exposure, adjustment_variables,con, series=NULL,
                                   logxform_e=T, scale_e=T,
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

  tab_obj <- name_and_xform_logistic_pheno_expo(pheno, exposure, tab_obj, logxform_e)

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
      base_models[[mod_num]] <- run_logistic_model(basebaseadjusted,dsn, scale_expo = scale_e,  quantile_expo=quantile_expo, expo_levels =  exposure_levels)
    }
    models[[mod_num]] <- run_logistic_model(baseadjusted,dsn, scale_expo = scale_e, quantile_expo=quantile_expo, expo_levels =  exposure_levels)
  }
  n <- dsn |> nrow()
  list(log_e = logxform_e, scaled_e=scale_e, unweighted_n=n, phenotype=pheno, series=tab_obj$series, exposure=exposure, models=models, base_models=base_models, adjustment_variables=adjustment_variables, demographic_breakdown=demo_break_tbl)
}


#' Run Logistic Model on Survey Data
#'
#' This function fits a logistic regression model using survey data, allowing for the scaling, quantile transformation, or categorical recoding of an exposure variable.
#'
#' @param formu A formula object specifying the logistic regression model to fit.
#' @param dsn A survey design object created using functions like `svydesign` from the **survey** package.
#' @param scale_expo Logical; if `TRUE`, the exposure variable will be scaled (z-score transformation). Default is `TRUE`.
#' @param quantile_expo Numeric vector of quantiles (e.g., `c(0, 0.25, 0.5, 0.75, 1)`) to discretize the exposure variable into quantile bins. Default is `NULL`.
#' @param expo_levels Character vector of specific levels to treat the exposure as a categorical variable. Default is `NULL`.
#' @param save_svymodel Logical; if `TRUE`, the full survey model object is included in the output. Default is `FALSE`.
#'
#' @return A list with the following elements:
#' \item{glanced}{A tibble with model summary statistics produced by `broom::glance`.}
#' \item{tidied}{A tibble with coefficient estimates and related statistics produced by `broom::tidy`.}
#' \item{scale_expo}{Logical, indicating whether the exposure variable was scaled.}
#' \item{q_cut_points}{Numeric vector of cut points used if the exposure was discretized into quantiles.}
#' \item{quantiles}{The input vector of quantiles used for discretization (if applicable).}
#' \item{expo_levels}{The input levels of the categorical exposure variable (if applicable).}
#' \item{model}{The full survey model object, if `save_svymodel = TRUE`; otherwise, `NA`.}
#'
#' @details
#' This function is designed to run logistic regression models on survey data, providing flexibility for handling exposure variables. It supports:
#' - Scaling the exposure variable (z-score transformation).
#' - Discretizing the exposure into quantile bins.
#' - Treating the exposure as a categorical variable with specified levels.
#'
#' The function utilizes the **survey** package for modeling, with `quasibinomial` as the default family for logistic regression.
#'
#' @examples
#' \dontrun{
#' library(survey)
#' library(broom)
#'
#' # Example survey data
#' data(api)
#' dstrat <- svydesign(id = ~1, strata = ~stype, weights = ~pw, data = apistrat, fpc = ~fpc)
#'
#' # Run logistic model with scaling
#' result <- run_logistic_model(api00 ~ ell + meals, dstrat, scale_expo = TRUE)
#'
#' # Access results
#' result$glanced
#' result$tidied
#' }
#'
#' @importFrom survey svyglm svymean svyvar svyquantile
#' @importFrom stats update
#' @importFrom broom glance tidy
#' @export
run_logistic_model <- function(formu, dsn, scale_expo=T, quantile_expo=NULL, expo_levels=NULL, save_svymodel=F) {
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


  mod <- survey::svyglm(formu,dsn, family=quasibinomial())
  ti <- broom::tidy(mod)
  gl <- broom::glance(mod)
  list(glanced=gl, tidied=ti,
       scale_expo=scale_expo,
       q_cut_points=cut_points, quantiles=quantile_expo, expo_levels=expo_levels, model=ifelse(save_svymodel, mod, NA)) # not saving dsn
}
