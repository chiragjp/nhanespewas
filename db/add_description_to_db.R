# script to add adds a "description" table to the database

#!/usr/local/bin/Rscript

library(getopt)
library(DBI)
library(tidyverse)

spec = matrix(c(
  'path', 'p', 1, "character",
  'notes', 'n', 1, "character",
  'is_summary', 's', 0, "logical"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

DATA_TYPE <- "individual-level"
if(!is.null(opt$is_summary)) {
  DATA_TYPE <- "summary stats"
}

PATH_TO_DB <- opt$path #'../inst/extdata/nhanes_pewas_a-d.sqlite'
NOTES <- opt$notes #"nhanes A, B, C, D"

print(PATH_TO_DB)
print(NOTES)
print(DATA_TYPE)

con <- DBI::dbConnect(RSQLite::SQLite(), dbname=PATH_TO_DB)
table_name <- 'description'
description_table <- tibble(release_date=as.character(Sys.Date()), release_package="nhanespewas", data_type=DATA_TYPE, notes=NOTES)
dbWriteTable(con, table_name, description_table, overwrite=TRUE)
dbDisconnect(con)
