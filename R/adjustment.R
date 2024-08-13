## create a data structure to accommodate multiple scenarios and domains of the exposome
adjustmentVariables <- c("RIDAGEYR", "AGE_SQUARED",
                         "RIAGENDR",
                         "INDFMPIR",
                         "EDUCATION_LESS9","EDUCATION_9_11","EDUCATION_AA","EDUCATION_COLLEGEGRAD",
                         "ETHNICITY_MEXICAN", "ETHNICITY_OTHERHISPANIC","ETHNICITY_OTHER", "ETHNICITY_NONHISPANICBLACK")


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
