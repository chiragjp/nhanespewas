# Chirag J Patel
# 12/23/2022
## load all of the tables dictionaries from a .rds into a sqlite database
## see download_dictionary_for_select_variables.R
## db_load_dictionary.R

library(tidyverse)
library(DBI)
library(getopt)

spec <- matrix(c(
  'dictionary_rds', 'd', 1, "character",
  'sqlite_name', 'o', 1, "character"
), byrow=TRUE, ncol=4);
opt <- getopt(spec)

dbname <- "nhanes_122322.sqlite"
dictionary_rds <- "./data/dictionary_DEMO.rds"
if(!is.null(opt$dictionary_rds)) {
  dictionary_rds <- opt$dictionary_rds
  dbname <- opt$sqlite_name
}

cat(sprintf("reading %s\n", dictionary_rds))

dict <- read_rds(dictionary_rds)
con <- dbConnect(RSQLite::SQLite(), dbname)

dbWriteTable(con, "variable_names_epcf", dict$dictionary, append=T)
dbWriteTable(con, "table_names_epcf", dict$tables %>% select(-c(Data.File.Description)), append=T)

dbDisconnect(con)