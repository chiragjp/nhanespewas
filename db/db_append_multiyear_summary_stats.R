
# load in the new overall summary statistics from a multi-year analysis
library(DBI)
library(tidyverse)

#load('../pe_summary_0824/gathered_overall_0824.Rdata')
#load('../pe_summary_0824/gathered_quantile_102124.Rdata')
#load('../pe_summary_0824/gathered_quantile_noscale_p_102124.Rdata')
#load('../pe_summary_0824/gathered_quantile_102224.Rdata')
#load('../pe_summary_0824/gathered_ptile_no_scale_p_102324.Rdata')
#load("../pe_summary_0125/gathered_overall_linear_0125.Rdata")
#load("../pe_summary_0125/gathered_quantile_scale_p_0125_v2.Rdata")
#load("../pe_summary_0125/gathered_quantile_scale_p_0125.Rdata")
load("../pe_summary_0824/mage/mage_gathered.Rdata")

full_data_con <- DBI::dbConnect(RSQLite::SQLite(), dbname='./nhanes_031725.sqlite')
varnames <- tbl(full_data_con, "variable_names_epcf") |> collect()


#con <- DBI::dbConnect(RSQLite::SQLite(), dbname='./pe_summary_stats_08_2024.sqlite')
#con <- DBI::dbConnect(RSQLite::SQLite(), dbname='./pe_summary_stats_quantile_10_2024.sqlite')
#con <- DBI::dbConnect(RSQLite::SQLite(), dbname='./pe_summary_stats_percentile_noscale_p_10_2024.sqlite')
con <- DBI::dbConnect(RSQLite::SQLite(), dbname='./pe_summary_stats_01_2025.sqlite')
#con <- DBI::dbConnect(RSQLite::SQLite(), dbname='./pe_summary_stats_percentile_noscale_p_01_2025_v2.sqlite')
#con <- DBI::dbConnect(RSQLite::SQLite(), dbname='./pe_summary_stats_percentile_scale_p_01_2025_v2.sqlite')
dbWriteTable(con, 'pe_overall', pe, overwrite=F, append=T)
dbWriteTable(con, 'glanced_overall', glanced, overwrite=F, append=T)
dbWriteTable(con, 'rsq_overall', rsq, overwrite=F, append=T)
dbWriteTable(con, "variable_names_epcf", varnames, overwrite = T, append=F)
dbDisconnect(con)

remove(list=ls())
