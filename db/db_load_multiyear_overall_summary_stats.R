
# load in the new overall summary statistics from a multi-year analysis

library(DBI)
library(tidyverse)

#path_to_summary <- '../pe_summary_020424/'
#load(file.path(path_to_summary, 'pe_summary_022524.Rdata'))
#load('../pe_summary_0824/gathered3.Rdata')
load('../pe_summary_0824/gathered_overall_0824.Rdata')

con <- DBI::dbConnect(RSQLite::SQLite(), dbname='./pe_summary_stats_08_2024.sqlite')

dbWriteTable(con, 'pe_overall', pe, overwrite=T, append=F)
dbWriteTable(con, 'glanced_overall', glanced, overwrite=T, append=F)
dbWriteTable(con, 'rsq_overall', rsq, overwrite=T, append=F)

dbDisconnect(con)
