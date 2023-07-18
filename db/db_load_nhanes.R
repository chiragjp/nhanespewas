# Chirag J Patel
# 12/23/2022
## load all of the tables from a .rds into a sqlite database
## see download_nhanesa.R

library(tidyverse)
library(DBI)
library(nhanesA)
library(getopt)

spec <- matrix(c(
  'component', 'c', 1, "character",
  'table_filename_rds', 't', 1, "character",
  'sqlite_name', 'o', 1, "character"
), byrow=TRUE, ncol=4);
opt <- getopt(spec)

component <- "DEMO"
tabs <- read_rds('./nhanes_DEMO_all_1222222.rds')
dbname <- "nhanes.sqlite"
if(!is.null(opt$component)) {
  component <- opt$component
  tabs <- read_rds(opt$table_filename_rds)
  dbname <- opt$sqlite_name
}

con <- dbConnect(RSQLite::SQLite(), dbname)

frames_to_db <- function(con, frames_for_years, overwrite=T, append=F) {
  for(y in 1:length(frames_for_years$years)) {
    cat(sprintf('%s\n', frames_for_years$years[y]))
    frms <- frames_for_years$tables[[y]]
    table_names <- names(frms)
    for(f in 1:length(frms)) {
      cat(sprintf('%s\n', table_names[f]))
      if(is.null(frms[[f]])) {
        next;
      }
      dbWriteTable(con, table_names[f], frms[[f]], overwrite=overwrite, append=append)
    }
  }
} 

frames_to_db(con, tabs)
table_names <- map(tabs$years, ~nhanesTables(component, .)) %>% bind_rows() %>% mutate(component = component)

dbWriteTable(con, "table_names", table_names, append=T)

dbDisconnect(con)


