

library(tidyverse)
library(DBI)

# 
# categorical levels for E variables
dbname <- "nhanes_122322.sqlite"
con  <- dbConnect(RSQLite::SQLite(), dbname)

# go through each E variable and collect the levels

variable_information <- tbl(con, "variable_names_epcf")
table_information <- tbl(con, "table_names_epcf") # c - control; e - exposure, f-filter; p - phenotype
variable_information_joined <- variable_information |> left_join(table_information |> select(Data.File.Name , epcf), by="Data.File.Name")
selected_variables <- read_csv("./select_variables_1.csv") |> filter(select == 1)
variable_information_selected <- variable_information_joined |> collect() |> filter(Variable.Name %in% selected_variables$Variable.Name)
e_tables <- variable_information_selected |> filter(epcf == 'e') |> pull(Data.File.Name) |> unique()


THRESHOLD <- 20
categorical_variables_in_table <- function(con, table_name) {
  tab <- tbl(con, table_name)
  stackjack <- tab |> collect() |>  select_if(is.numeric) |> map(function(x) {
    levs <- unique(x)
    if(length(levs) <= THRESHOLD) {
      return(levs)
    }
    return(0)
  }) |> stack()
  
}

categorical_variables_in_tables <- function(con, tables) {
  stacks <- vector(mode="list", length(tables))
  for(tbl in 1:length(tables))  {
    stacks[[tbl]] <- categorical_variables_in_table(con, tables[tbl]) |> rename(Variable.Name = ind) |> mutate(Data.File.Name = tables[tbl])
  }
  stacks |> bind_rows()
}

## e_tables
in_db <- dbListTables(con)
e_tables <- intersect(in_db, e_tables)
e_categorical_levels <- categorical_variables_in_tables(con, e_tables)

#the_questionable_categoricals  <- e_categorical_levels |> filter(values > 0, values != round(values))


## convert to rank and do not analyze as categorical 
## if number of discrete categories is only 1, ?

dbWriteTable(con, "e_variable_levels", e_categorical_levels, overwrite=T)
demo_tables <- table_information |> filter(component == 'DEMO') |> pull(Data.File.Name)
demo_levels <- categorical_variables_in_tables(con, demo_tables)
dbWriteTable(con, "demo_variable_levels", demo_levels, overwrite=T)

