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

