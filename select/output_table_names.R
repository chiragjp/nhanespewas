# Chirag J Patel
# 12/23/2022

## outputs the table names in order to annotate variables needed for E-P

library(tidyverse)
library(DBI)

con <- dbConnect(RSQLite::SQLite(), "nhanes_122322.sqlite")

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


tables_list <- tbl(con, "table_names") %>% as_tibble()
tables_list <- tables_list %>% mutate(series=series_from_table_name(Data.File.Name)) %>% arrange(Data.File.Name) %>% write_csv(file="select_tables.csv")
  # now categorize and move to next line, setting this flag to F  

# categorize these as epcf
tables_list <- read_csv(file = "select_tables.csv")
dbDisconnect(con)


# Phenotypes - 'p'
# Exam variables - adiposity, height, weight,blood pressure, pulse rate
# Lab variables - biochemistry

# Exposome - 'e'
# metals, pesticides, allergy, 
# DIET, total with adjustment for total calorie

# Covariates - baseline - 'c'
# age, sex, poverty, education
# DEMO table

# Inclusion criteria - 'f'
# filter: disease status (MCQ)

