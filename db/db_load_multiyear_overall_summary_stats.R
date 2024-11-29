
# load in the new overall summary statistics from a multi-year analysis

library(DBI)
library(tidyverse)

#load('../pe_summary_0824/gathered_overall_0824.Rdata')
#load('../pe_summary_0824/gathered_quantile_102124.Rdata')
#load('../pe_summary_0824/gathered_quantile_noscale_p_102124.Rdata')
#load('../pe_summary_0824/gathered_quantile_102224.Rdata')
load('../pe_summary_0824/gathered_ptile_no_scale_p_102324.Rdata')


#con <- DBI::dbConnect(RSQLite::SQLite(), dbname='./pe_summary_stats_08_2024.sqlite')
#con <- DBI::dbConnect(RSQLite::SQLite(), dbname='./pe_summary_stats_quantile_10_2024.sqlite')
con <- DBI::dbConnect(RSQLite::SQLite(), dbname='./pe_summary_stats_percentile_noscale_p_10_2024.sqlite')

dbWriteTable(con, 'pe_overall', pe, overwrite=T, append=F)
dbWriteTable(con, 'glanced_overall', glanced, overwrite=T, append=F)
dbWriteTable(con, 'rsq_overall', rsq, overwrite=T, append=F)

dbDisconnect(con)

remove(list=ls())
