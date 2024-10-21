library(tidyverse)
library(DBI)
## this loads the demographic "breakdown" to the summary stats table

#path_to_summary_stats_db <- '../pe_summary_stats.sqlite'
#path_to_summary_stats_db <- '../pe_summary_stats.sqlite'
summary_stats_con <- DBI::dbConnect(RSQLite::SQLite(), dbname='pe_summary_stats_08_2024.sqlite')
#summary_stats_con <- DBI::dbConnect(RSQLite::SQLite(), dbname=path_to_summary_stats_db)

path_to_summary <- '../pe_summary_060623/'
load(file.path(path_to_summary, 'demographic_breakdown.Rdata'))
dbWriteTable(summary_stats_con, "demographic_breakdown", demo_breakdowns, overwrite=T, append=F)
dbDisconnect(summary_stats_con)
