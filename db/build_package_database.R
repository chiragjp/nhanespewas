## build a .sqlite NHANES database for packaging
## should be < 2GB to also commit to github
## put in extdata

library(DBI)
library(tidyverse)
series_to_save <- c("C")
to_con <- DBI::dbConnect(RSQLite::SQLite(), dbname='./nhanes_pewas_test.sqlite')
from_con <- DBI::dbConnect(RSQLite::SQLite(), dbname='./nhanes_122322.sqlite')

tables_in_db <- dbListTables(from_con)
table_names_epcf <- tbl(from_con, "table_names_epcf")

housekeeping <- c('demo_variable_levels', 'e_variable_levels', 'table_names_epcf', 'variable_names_epcf', 'table_names')

## select those in series A, B, C, D
#series_to_save <- c("A", "B", "C", "D")

table_names_epcf |> filter(series %in% series_to_save) |> collect() |> nrow()
table_names_epcf |> collect() |> nrow()


read_write_table <- function(from_c, to_c, table_name) {
  tab <- tbl(from_c, table_name) |> collect()
  dbWriteTable(to_c, table_name, tab)
}

table_names_to_copy <- table_names_epcf |> filter(series %in% series_to_save) |> collect()
#for(ii in 1:nrow(table_names_to_copy)) {
#  print(table_names_to_copy$Data.File.Name[ii])
#  read_write_table(from_con, to_con, table_names_to_copy$Data.File.Name[ii])
#}

## also write in the housekeeping tables
for(ii in 1:length(housekeeping)) {
  print(housekeeping[ii])
  read_write_table(from_con, to_con, housekeeping[ii])
}


dbDisconnect(from_con)
dbDisconnect(to_con)


