## create a data structure to accommodate multiple scenarios and domains of the exposome
adjustmentVariables <- c("RIDAGEYR", "AGE_SQUARED",
                         "RIAGENDR",
                         "INDFMPIR",
                         "EDUCATION_LESS9","EDUCATION_9_11","EDUCATION_AA","EDUCATION_COLLEGEGRAD",
                         "ETHNICITY_MEXICAN", "ETHNICITY_OTHERHISPANIC","ETHNICITY_OTHER", "ETHNICITY_NONHISPANICBLACK")

SURVEY_YEAR_VARIABLE <- "SDDSRVYR"

adjustment_models <- rbind(
  tibble(scenario="base", variables=NA),
  tibble(scenario="age_sex_ethnicity_income_education", variables=adjustmentVariables), #5
  tibble(scenario="age_sex", variables=c("RIDAGEYR", "AGE_SQUARED", "RIAGENDR")), #3
  tibble(scenario="age", variables=c("RIDAGEYR", "AGE_SQUARED")), #2
  tibble(scenario="sex", variables=c("RIAGENDR")), #1
  tibble(scenario="age_sex_ethnicity", variables=c("RIDAGEYR", "AGE_SQUARED", "RIAGENDR", "ETHNICITY_MEXICAN", "ETHNICITY_OTHERHISPANIC","ETHNICITY_OTHER", "ETHNICITY_NONHISPANICBLACK")),
  tibble(scenario="age_sex_income_education", variables=c("RIDAGEYR", "AGE_SQUARED","RIAGENDR", "INDFMPIR", "EDUCATION_LESS9","EDUCATION_9_11","EDUCATION_AA","EDUCATION_COLLEGEGRAD")), #4
  tibble(scenario="income_education", variables=c("INDFMPIR", "EDUCATION_LESS9","EDUCATION_9_11","EDUCATION_AA","EDUCATION_COLLEGEGRAD")),
  tibble(scenario="ethnicity", variables=c("ETHNICITY_MEXICAN", "ETHNICITY_OTHERHISPANIC","ETHNICITY_OTHER", "ETHNICITY_NONHISPANICBLACK"))
) |> mutate(domain = 'default')

adjustment_models_diet_x <- rbind(
  tibble(scenario="base", variables=c("DRXTKCAL_adj")),
  tibble(scenario="age_sex_ethnicity_income_education", variables=c(adjustmentVariables, "DRXTKCAL_adj")),
  tibble(scenario="age_sex", variables=c("RIDAGEYR", "AGE_SQUARED", "RIAGENDR", "DRXTKCAL_adj")),
  tibble(scenario="age", variables=c("RIDAGEYR", "AGE_SQUARED", "DRXTKCAL_adj")),
  tibble(scenario="sex", variables=c("RIAGENDR", "DRXTKCAL_adj")),
  tibble(scenario="age_sex_ethnicity", variables=c("RIDAGEYR", "AGE_SQUARED", "RIAGENDR", "ETHNICITY_MEXICAN", "ETHNICITY_OTHERHISPANIC","ETHNICITY_OTHER", "ETHNICITY_NONHISPANICBLACK", "DRXTKCAL_adj")),
  tibble(scenario="age_sex_income_education", variables=c("RIDAGEYR", "AGE_SQUARED","RIAGENDR", "INDFMPIR", "EDUCATION_LESS9","EDUCATION_9_11","EDUCATION_AA","EDUCATION_COLLEGEGRAD", "DRXTKCAL_adj")),
  tibble(scenario="income_education", variables=c("INDFMPIR", "EDUCATION_LESS9","EDUCATION_9_11","EDUCATION_AA","EDUCATION_COLLEGEGRAD", "DRXTKCAL_adj")),
  tibble(scenario="ethnicity", variables=c("ETHNICITY_MEXICAN", "ETHNICITY_OTHERHISPANIC","ETHNICITY_OTHER", "ETHNICITY_NONHISPANICBLACK", "DRXTKCAL_adj"))
) |> mutate(domain = 'DRX') # dietary recall

adjustment_models_diet_1 <- rbind(
  tibble(scenario="base", variables=c("DR1TKCAL_adj")),
  tibble(scenario="age_sex_ethnicity_income_education", variables=c(adjustmentVariables, "DR1TKCAL_adj")),
  tibble(scenario="age_sex", variables=c("RIDAGEYR", "AGE_SQUARED", "RIAGENDR", "DR1TKCAL_adj")),
  tibble(scenario="age", variables=c("RIDAGEYR", "AGE_SQUARED", "DR1TKCAL_adj")),
  tibble(scenario="sex", variables=c("RIAGENDR", "DR1TKCAL_adj")),
  tibble(scenario="age_sex_ethnicity", variables=c("RIDAGEYR", "AGE_SQUARED", "RIAGENDR", "ETHNICITY_MEXICAN", "ETHNICITY_OTHERHISPANIC","ETHNICITY_OTHER", "ETHNICITY_NONHISPANICBLACK", "DR1TKCAL_adj")),
  tibble(scenario="age_sex_income_education", variables=c("RIDAGEYR", "AGE_SQUARED","RIAGENDR", "INDFMPIR", "EDUCATION_LESS9","EDUCATION_9_11","EDUCATION_AA","EDUCATION_COLLEGEGRAD", "DR1TKCAL_adj")),
  tibble(scenario="income_education", variables=c("INDFMPIR", "EDUCATION_LESS9","EDUCATION_9_11","EDUCATION_AA","EDUCATION_COLLEGEGRAD", "DR1TKCAL_adj")),
  tibble(scenario="ethnicity", variables=c("ETHNICITY_MEXICAN", "ETHNICITY_OTHERHISPANIC","ETHNICITY_OTHER", "ETHNICITY_NONHISPANICBLACK", "DR1TKCAL_adj"))
) |> mutate(domain = 'DR1') # dietary recall


adjustment_models_diet_2 <- rbind(
  tibble(scenario="base", variables=c("DR2TKCAL_adj")),
  tibble(scenario="age_sex_ethnicity_income_education", variables=c(adjustmentVariables, "DR2TKCAL_adj")),
  tibble(scenario="age_sex", variables=c("RIDAGEYR", "AGE_SQUARED", "RIAGENDR", "DR2TKCAL_adj")),
  tibble(scenario="age", variables=c("RIDAGEYR", "AGE_SQUARED", "DR2TKCAL_adj")),
  tibble(scenario="sex", variables=c("RIAGENDR", "DR2TKCAL_adj")),
  tibble(scenario="age_sex_ethnicity", variables=c("RIDAGEYR", "AGE_SQUARED", "RIAGENDR", "ETHNICITY_MEXICAN", "ETHNICITY_OTHERHISPANIC","ETHNICITY_OTHER", "ETHNICITY_NONHISPANICBLACK", "DR2TKCAL_adj")),
  tibble(scenario="age_sex_income_education", variables=c("RIDAGEYR", "AGE_SQUARED","RIAGENDR", "INDFMPIR", "EDUCATION_LESS9","EDUCATION_9_11","EDUCATION_AA","EDUCATION_COLLEGEGRAD", "DR2TKCAL_adj")),
  tibble(scenario="income_education", variables=c("INDFMPIR", "EDUCATION_LESS9","EDUCATION_9_11","EDUCATION_AA","EDUCATION_COLLEGEGRAD", "DR2TKCAL_adj")),
  tibble(scenario="ethnicity", variables=c("ETHNICITY_MEXICAN", "ETHNICITY_OTHERHISPANIC","ETHNICITY_OTHER", "ETHNICITY_NONHISPANICBLACK", "DR2TKCAL_adj"))
) |> mutate(domain = 'DR2') # dietary recall

adjustment_models_ucr <- rbind(
  tibble(scenario="base", variables=c("URXUCR_adj")),
  tibble(scenario="age_sex_ethnicity_income_education", variables=c(adjustmentVariables, "URXUCR_adj")),
  tibble(scenario="age_sex", variables=c("RIDAGEYR", "AGE_SQUARED", "RIAGENDR", "URXUCR_adj")),
  tibble(scenario="age", variables=c("RIDAGEYR", "AGE_SQUARED", "URXUCR_adj")),
  tibble(scenario="sex", variables=c("RIAGENDR", "URXUCR_adj")),
  tibble(scenario="age_sex_ethnicity", variables=c("RIDAGEYR", "AGE_SQUARED", "RIAGENDR", "ETHNICITY_MEXICAN", "ETHNICITY_OTHERHISPANIC","ETHNICITY_OTHER", "ETHNICITY_NONHISPANICBLACK", "URXUCR_adj")),
  tibble(scenario="age_sex_income_education", variables=c("RIDAGEYR", "AGE_SQUARED","RIAGENDR", "INDFMPIR", "EDUCATION_LESS9","EDUCATION_9_11","EDUCATION_AA","EDUCATION_COLLEGEGRAD", "URXUCR_adj")),
  tibble(scenario="income_education", variables=c("INDFMPIR", "EDUCATION_LESS9","EDUCATION_9_11","EDUCATION_AA","EDUCATION_COLLEGEGRAD", "URXUCR_adj")),
  tibble(scenario="ethnicity", variables=c("ETHNICITY_MEXICAN", "ETHNICITY_OTHERHISPANIC","ETHNICITY_OTHER", "ETHNICITY_NONHISPANICBLACK", "URXUCR_adj"))
) |> mutate(domain = 'URX') # urine exposure



#' Determine Adjustment Scenario for a Given Exposure and Phenotype
#'
#' This function determines the appropriate adjustment scenario based on the names of the exposure and phenotype variables. It returns a predefined set of adjustment models specific to the type of variable provided.
#'
#' @param evarname A character string specifying the name of the exposure variable.
#' @param pvarname A character string specifying the name of the phenotype variable.
#' @return The function returns one of the predefined adjustment model sets based on the variable name patterns:
#' \itemize{
#'   \item `adjustment_models_diet_x`: If the exposure variable starts with `"DRX"`, indicating diet-related data adjusted by total calories.
#'   \item `adjustment_models_diet_1`: If the exposure variable starts with `"DR1"`, indicating first-day diet recall data adjusted by total calories.
#'   \item `adjustment_models_diet_2`: If the exposure variable starts with `"DR2"`, indicating second-day diet recall data adjusted by total calories.
#'   \item `adjustment_models_ucr`: If the exposure variable starts with `"UR"` and the phenotype variable does not start with `"UR"`, indicating urine data adjusted by creatinine.
#'   \item `adjustment_models`: A default set of adjustment models if none of the above conditions are met.
#' }
#' @details
#' The function uses the first few characters of the exposure and phenotype variable names to determine the appropriate adjustment scenario. It logs the chosen adjustment type using `log_info`.
#' @examples
#' \dontrun{
#' evarname <- "DR1TKCAL"
#' pvarname <- "BMXBMI"
#' adjustment_model <- adjustment_scenario_for_variable(evarname, pvarname)
#' print(adjustment_model)
#' }
#' @export

adjustment_scenario_for_variable <- function(evarname, pvarname) {
  first_three <- substr(evarname, 1, 3)
  first_two <- substr(evarname, 1, 2)
  first_two_p <- substr(pvarname, 1, 2)
  if(first_three == 'DRX') {
    log_info("Adjusting by total calories")
    return(adjustment_models_diet_x)
  } else if(first_three == 'DR1') {
    log_info("Adjusting by total calories")
    return(adjustment_models_diet_1)
  } else if(first_three == 'DR2') {
    log_info("Adjusting by total calories")
    return(adjustment_models_diet_2)
  } else if(first_two == 'UR' & first_two_p != 'UR') {
    log_info("Adjusting by creatinine")
    return(adjustment_models_ucr)
  }else {
    log_info("Default adjustments")
    return(adjustment_models)
  }
}


#' Add a Covariate to an Adjustment Scenario
#'
#' This function adds a specified covariate to a given adjustment scenario within an adjustment model.
#'
#' @param adjustment_model A data frame or tibble representing the adjustment model,
#' containing at least the columns `domain`, `scenario`, and `variables`.
#' @param scenario_name A character string specifying the name of the adjustment scenario
#' to which the covariate should be added.
#' @param variable_name A character string specifying the name of the covariate (variable) to add.
#'
#' @return A modified version of the `adjustment_model` with the new covariate added to the specified scenario.
#'
#' @details
#' The function assumes that the input `adjustment_model` contains a `domain` column, and the
#' domain value for the new entry is extracted from the first row of the input model. The new
#' scenario and variable are appended as a new row to the `adjustment_model`.
#'
#' @examples
#' # Example adjustment model
#' adjustment_model <- tibble::tibble(
#'   domain = "example_domain",
#'   scenario = c("age_sex", "age_sex"),
#'   variables = c("age", "sex")
#' )
#'
#' # Add a new covariate
#' updated_model <- add_covariate_to_scenario(
#'   adjustment_model,
#'   scenario_name = "age",
#'   variable_name = "as.factor(SDDSRVYR)"
#' )
#'
#' print(updated_model)
#'
#' @export
add_covariate_to_scenario <- function(adjustment_model, scenario_name, variable_name) {
  domain_name <- adjustment_model |> first() |> pull("domain")
  adjustment_model <- adjustment_model |> rbind(tibble(domain=domain_name, scenario=scenario_name, variables=variable_name))
  adjustment_model
}

#' Add Survey Year to Adjustment Scenarios
#'
#' This function updates an adjustment model for exposures by adding the survey year
#' (`SDDSRVYR`) as a covariate across multiple adjustment scenarios. It applies the
#' `nhanespewas::add_covariate_to_scenario()` function to incorporate the survey year
#' in various combinations of demographic and socioeconomic variables.
#'
#' @param adjustment_model_for_e A list or data structure representing the adjustment
#' model for exposures. It is expected to be compatible with the
#' `nhanespewas::add_covariate_to_scenario()` function.
#'
#' @return The updated adjustment model for exposures, with the survey year added as
#' a covariate across various adjustment scenarios.
#'
#' @details
#' The following adjustment scenarios are updated by adding the survey year as
#' a covariate:
#' - `age_sex_ethnicity_income_education`
#' - `age_sex`
#' - `age`
#' - `sex`
#' - `age_sex_ethnicity`
#' - `age_sex_income_education`
#' - `income_education`
#' - `ethnicity`
#'
#' The survey year is added as a factor using `as.factor(SDDSRVYR)`.
#'
#' @examples
#' # Assuming `adjustment_model_for_e` is predefined and compatible:
#' updated_model <- add_survey_year_to_adjustment_scenarios(adjustment_model_for_e)
#'
#' @export
add_survey_year_to_adjustment_scenarios <- function(adjustment_model_for_e) {
  adjustment_model_for_e <- nhanespewas::add_covariate_to_scenario(adjustment_model_for_e,
                                                                   "age_sex_ethnicity_income_education",
                                                                   variable_name = "as.factor(SDDSRVYR)")

  adjustment_model_for_e <- nhanespewas::add_covariate_to_scenario(adjustment_model_for_e,
                                                                   "age_sex",
                                                                   variable_name = "as.factor(SDDSRVYR)")
  adjustment_model_for_e <- nhanespewas::add_covariate_to_scenario(adjustment_model_for_e,
                                                                   "age",
                                                                   variable_name = "as.factor(SDDSRVYR)")

  adjustment_model_for_e <- nhanespewas::add_covariate_to_scenario(adjustment_model_for_e,
                                                                   "sex",
                                                                   variable_name = "as.factor(SDDSRVYR)")

  adjustment_model_for_e <- nhanespewas::add_covariate_to_scenario(adjustment_model_for_e,
                                                                   "age_sex_ethnicity",
                                                                   variable_name = "as.factor(SDDSRVYR)")
  adjustment_model_for_e <- nhanespewas::add_covariate_to_scenario(adjustment_model_for_e,
                                                                   "age_sex_income_education",
                                                                   variable_name = "as.factor(SDDSRVYR)")
  adjustment_model_for_e <- nhanespewas::add_covariate_to_scenario(adjustment_model_for_e,
                                                                   "income_education",
                                                                   variable_name = "as.factor(SDDSRVYR)")

  adjustment_model_for_e <- nhanespewas::add_covariate_to_scenario(adjustment_model_for_e,
                                                                   "ethnicity",
                                                                   variable_name = "as.factor(SDDSRVYR)")
  adjustment_model_for_e

}





