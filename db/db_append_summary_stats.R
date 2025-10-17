# load summary stats to sqlite for the multiple adjustment scenarios
# appending new data

library(tidyverse)
library(DBI)

con <- DBI::dbConnect(RSQLite::SQLite(), dbname='./nhanes_012324.sqlite')
varnames <- tbl(con, "variable_names_epcf") |> collect()
expo_levels <- tbl(con, "e_variable_levels") |> collect()
tablenames <- tbl(con, "table_names_epcf") |> collect()
dbDisconnect(con)

path_to_summary <- '../pe_summary_0824/'
load(file.path(path_to_summary, 'mage', 'mage_gathered.Rdata'))
#load('../pe_summary_0824/mage/mage_gathered.Rdata')

domains <- read_csv('../select/variable_domains_ep_2.csv')

#con <- DBI::dbConnect(RSQLite::SQLite(), dbname='./pe_summary_stats_08_2024.sqlite')
con <- DBI::dbConnect(RSQLite::SQLite(), dbname='./pe_summary_stats_01_2025.sqlite')
#dbWriteTable(con, "variable_names_epcf", varnames, overwrite=T, append=F)
#dbWriteTable(con, "table_names_epcf", tablenames, overwrite=T, append=F)
dbWriteTable(con, 'pe', pe, overwrite=F, append=T)
dbWriteTable(con, 'glanced', glanced, overwrite=F, append=T)
dbWriteTable(con, 'rsq', rsq, overwrite=F, append=T)
dbWriteTable(con, 'variable_domain', domains, overwrite=T, append=F)


remove(list=ls())
