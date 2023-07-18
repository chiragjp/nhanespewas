## pull dictionaries for component

library(tidyverse)
library(nhanesA)
library(getopt)

spec <- matrix(c(
  'component_module', 'c', 1, "character"
), byrow=TRUE, ncol=4);
opt <- getopt(spec)


component_module <- opt$component_module 

## get data dictionaries for each table in select_tables.csv
tables_list <- read_csv(file = "select_tables.csv")
tables_list <- tables_list %>% filter(!is.na(epcf), component == component_module) 
dictionary <- vector("list", nrow(tables_list))
cat(sprintf("Dictionaries for %i number of tables\n", nrow(tables_list)))
for(i in 1:nrow(tables_list)) {
  cat(sprintf("%i,%s\n", i, tables_list$Data.File.Name[i]))
  dictionary[[i]] <- nhanesTableVars(tables_list$component[i], tables_list$Data.File.Name[i], nchar = 75, details = T) 
}

dictionary <- dictionary %>% bind_rows()

fileout <- sprintf("dictionary_%s.rds", component_module)
dictionary_struct <- list(
  dictionary = dictionary,
  component = component_module,
  tables=tables_list
)
saveRDS(dictionary_struct, file=fileout)