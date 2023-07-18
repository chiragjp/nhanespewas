# Chirag J Patel
# 12/23/2022

## select E and P variable tables from the sqlite table

library(tidyverse)
library(DBI)

con <- dbConnect(RSQLite::SQLite(), "nhanes_122322.sqlite")

CLEAN_SLATE  <- F
series_from_table_name <- function(table_name_string_vec) {
  map_chr(
    str_split(table_name_string_vec, "\\_"),
    function(split_list) {
      if(length(split_list)==1) {
        return("A")
      } else{
        return(split_list[[length(split_list)]])
      }
    }
  )
}

if(CLEAN_SLATE) {
  tables_list <- tbl(con, "table_names") %>% as_tibble()
  tables_list <- tables_list %>% mutate(series=series_from_table_name(Data.File.Name)) %>% arrange(Data.File.Name) %>% write_csv(file="select_vars.csv")
  # now categorize and move to next line, setting this flag to F  
}
# categorize these as epcf
tables_list <- read_csv(file = "select_vars.csv")

## get data dictionaries for each
tables_list <- tables_list %>% filter(!is.na(epcf))
for(i in 1:nrow(tables_list)) {
  cat(sprintf("%s\n", tables_list$Data.File.Name[i]))
  tabl_var <- nhanesTableVars(tables_list$component[i], tables_list$Data.File.Name[i], nchar = 75) 
  #dbWriteTable(con, "variable_names", tabl_var, append=T)
}



# Phenotypes:
# Exam variables - adiposity, height, weight,blood pressure, pulse rate
# Lab variables - biochemistry

# Exposome
# metals, pesticides, allergy, 
# DIET, total with adjustment for total calorie

# Covariates - baseline
# age, sex, poverty, education
# DEMO table

# Inclusion criteria
# filter: disease status (MCQ)

dbDisconnect(con)