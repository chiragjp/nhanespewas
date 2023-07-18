## insert covariates in demo tables
## 


library(tidyverse)
library(DBI)
con <- DBI::dbConnect(RSQLite::SQLite(), dbname='nhanes_122322.sqlite')


ethnicity <- function(demo_table) {
  demo_table |> mutate(ETHNICITY_NONHISPANICWHITE=ifelse(RIDRETH1 == 3, 1, 0), 
                       ETHNICITY_MEXICAN=ifelse(RIDRETH1 == 1, 1, 0),
                       ETHNICITY_OTHERHISPANIC=ifelse(RIDRETH1 == 2, 1, 0),
                       ETHNICITY_OTHER=ifelse(RIDRETH1 == 5, 1, 0),
                       ETHNICITY_NONHISPANICBLACK=ifelse(RIDRETH1 == 4, 1, 0)
                      )
}



education <- function(demo_table) {
  #demo_table |> mutate(AGE_SQUARED=RIDAGEYR*RIDRAGEYR)
  #DMDEDUC2 for adults over 20
  #DMDEDUC3 for kids under 20
  # demo_table |> mutate(
  #   EDUCATION_LESS9 = ifelse(DMDEDUC2==1 | (!is.na(DMDEDUC3) & (DMDEDUC3 < 9 | DMDEDUC3 == 55 | DMDEDUC3 == 66)), 1, 0),
  #   EDUCATION_9_11 = ifelse(DMDEDUC2==2 | (!is.na(DMDEDUC3) & DMDEDUC3 >= 9 & DMDEDUC3 <= 11), 1, 0),
  #   EDUCATION_HSGRAD = ifelse(DMDEDUC2==3 | (!is.na(DMDEDUC3) & DMDEDUC3 > 12 & DMDEDUC3 <= 15), 1, 0),
  #   EDUCATION_AA = ifelse(DMDEDUC2==4, 1, 0),
  #   EDUCATION_COLLEGEGRAD = ifelse(DMDEDUC2==5, 1, 0)
  # )
  
  demo_table <- demo_table |>  mutate(EDUCATION_LESS9 = case_when(
    DMDEDUC2 == 1 ~ 1,
    (DMDEDUC3 < 9 | DMDEDUC3 == 55 | DMDEDUC3 == 66) ~ 1, 
    !is.na(DMDEDUC2) & DMDEDUC2 < 7 ~ 0,
    (DMDEDUC3 >= 9 | DMDEDUC3 != 55 | DMDEDUC3 != 66) ~ 0,
    TRUE ~ NA_real_
  ))
  
  demo_table <- demo_table |> mutate(EDUCATION_9_11 = case_when(
    DMDEDUC2 == 2 ~ 1,
    (DMDEDUC3 >= 9 & DMDEDUC3 <= 12)  ~ 1, 
    !is.na(DMDEDUC2) & DMDEDUC2 < 7 ~ 0,
    (!(DMDEDUC3 >= 9 & DMDEDUC3 <= 12) & DMDEDUC3 <= 66)  ~ 0,
    TRUE ~ NA_real_
  ))
  
  demo_table <- demo_table |> mutate(EDUCATION_HSGRAD = case_when(
    DMDEDUC2 == 3 ~ 1,
    (DMDEDUC3 == 13 | DMDEDUC3 == 14) ~ 1, 
    !is.na(DMDEDUC2) & DMDEDUC2 < 7 ~ 0,
    DMDEDUC3 != 13 & DMDEDUC3 != 14 & DMDEDUC3 < 77 ~ 0,
    TRUE ~ NA_real_
  ))
  
  demo_table <- demo_table |> mutate(EDUCATION_AA = case_when(
    DMDEDUC2 == 4 ~ 1,
    !is.na(DMDEDUC2) & DMDEDUC2 < 7 ~ 0,
    !is.na(DMDEDUC3) & DMDEDUC3 < 77 ~ 0,
    TRUE ~ NA_real_
  ))
  
  demo_table <- demo_table |> mutate(EDUCATION_COLLEGEGRAD = case_when(
    DMDEDUC2 == 5 ~ 1,
    !is.na(DMDEDUC2) & DMDEDUC2 < 7 ~ 0,
    !is.na(DMDEDUC3) & DMDEDUC3 < 77 ~ 0,
    TRUE ~ NA_real_
  ))
  
}


age_squared <- function(demo_table) {
  # age^2
  demo_table |> mutate(AGE_SQUARED=RIDAGEYR*RIDAGEYR)
}

born_in_us_or_elsewhere <- function(demo_table) {
  
  # DMDBORN
  #1	Born in 50 US States or Washington, DC	9408	9408	
  #2	Born in Mexico	993	10401	
  #3	Born Elsewhere	629	11030	
  #7	Refused	4	11034	
  #9	Don't know	0	11034
  # DMDBORN2
  #1	Born in 50 US States or Washington, DC	8399	8399	
  #2	Born in Mexico	730	9129	
  #4	Born in Other Spanish Speaking Country	583	9712	
  #5	Born in Other Non-Spanish Speaking Country	431	10143	
  #7	Refused	5	10148	
  #9	Don't Know	1	10149
  
  #DMDBORN4
  #1	Born in 50 US states or Washington, DC	7668	7668	
  #2	Others	2083	9751	
  #77	Refused	2	9753	
  #99	Don't Know	3	9756	
  #.	Missing	0	9756	
  
  if("DMDBORN" %in% colnames(demo_table)) {
    return(demo_table |> mutate(BORN_INUSA = case_when(DMDBORN == 1 ~ 1, 
                                                       DMDBORN > 1 & DMDBORN < 7 ~ 0,
                                                       TRUE ~ NA_real_
                                                       ) ) )
  }
  
  if("DMDBORN4" %in% colnames(demo_table)) {
    return(demo_table |> mutate(BORN_INUSA = case_when(DMDBORN4 == 1 ~ 1, 
                                                       DMDBORN4 > 1 & DMDBORN4 < 7 ~ 0,
                                                       TRUE ~ NA_real_
    ) ) )
  } 
  
  if("DMDBORN2" %in% colnames(demo_table)) {
    return(demo_table |> mutate(BORN_INUSA = case_when(DMDBORN2 == 1 ~ 1, 
                                                       DMDBORN2 > 1 & DMDBORN2 < 7 ~ 0,
                                                       TRUE ~ NA_real_
    ) ) )
  } 
  
  demo_table
  
}




table_descriptions <- tbl(con, "table_names_epcf") |> filter(epcf == 'c') |> collect()
tables <- vector(mode="list", length=nrow(table_descriptions))
for(ii in 1:nrow(table_descriptions)) {
  cat(sprintf("%s\n", table_descriptions$Data.File.Name[ii]))
  demo_tab <- tbl(con, table_descriptions$Data.File.Name[ii]) |> collect()
  tables[[ii]] <- demo_tab |> ethnicity() |> age_squared() |> education() |> born_in_us_or_elsewhere()
  DBI::dbWriteTable(con, table_descriptions$Data.File.Name[ii], tables[[ii]], overwrite=T, append=F)
}

adjustment_variables <- c("RIDAGEYR", "AGE_SQUARED", 
                          "RIAGENDR", 
                          "INDFMPIR",
                          "EDUCATION_LESS9","EDUCATION_9_11","EDUCATION_AA","EDUCATION_COLLEGEGRAD",
                          "BORN_INUSA",
                          "ETHNICITY_MEXICAN", "ETHNICITY_OTHERHISPANIC","ETHNICITY_OTHER", "ETHNICITY_NONHISPANICBLACK")