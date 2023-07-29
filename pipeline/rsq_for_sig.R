## Chirag Patel
## multivariate additive association between lifestyle exposures and a phenotype for a survey
## 05/15/23

library(furrr)
#library(this.path)
#setwd(this.dir())
#source('pe_rsq.R')
library(devtools)
load_all("..")
source('db_paths.R')
#path_to_nhanes <- '../nhanes_122322.sqlite'
#path_to_summary_stats <- '../pe_summary_stats.sqlite'

nhanes_con <- DBI::dbConnect(RSQLite::SQLite(), dbname=path_to_nhanes)
summary_stats_con <- DBI::dbConnect(RSQLite::SQLite(), dbname=path_to_summary_stats)
adjusted_meta_2 <- tbl(summary_stats_con, "adjusted_meta_2")
evars_for_ps <- adjusted_meta_2 |> filter(sig_levels == 'Bonf.<0.05')  |> group_by(pvarname) |> count() |> collect()
#evars_for_ps <- adjusted_meta_2 |> filter(sig_levels == 'BY<0.05')  |> group_by(pvarname) |> count() |> collect()
evars_for_ps <- evars_for_ps |> filter(!(pvarname %in% c('LBXRBFSI','LBXFOLSI','LBDIRNSI','LBXINSI', 'LBXMMASI')))

plan(multisession)
r2s <- future_map_dfr(evars_for_ps$pvarname, function(pvarname) {
  nhanes_con <- DBI::dbConnect(RSQLite::SQLite(), dbname=path_to_nhanes)
  summary_stats_con <- DBI::dbConnect(RSQLite::SQLite(), dbname=path_to_summary_stats)
  combined_dat <- get_individual_level_table(summary_stats_con, nhanes_con, pvarname)
  dbDisconnect(nhanes_con)
  dbDisconnect(summary_stats_con)
  svy_r2 <- svy_weighted_r2(pvarname,combined_dat$evars, adjustmentVariables, dat=combined_dat$selected_data, weight_name = combined_dat$weight_to_use$weight_name)
  tibble(pvarname=pvarname,
         n_evars=length(combined_dat$evars),
         n=combined_dat$selected_data |> nrow(),
         base_rsq=svy_r2$base$rsq,base_adj_rsq=svy_r2$base$adj.r2,mve_rsq = svy_r2$mve$rsq,mve_adj_rsq=svy_r2$mve$adj.r2)
}, .options = furrr_options(seed = TRUE))


write_rds(r2s, file='../pe_summary_060623/mvr2.rds')

dbWriteTable(summary_stats_con, 'mvr2', r2s, overwrite=T)
dbDisconnect(nhanes_con)
dbDisconnect(summary_stats_con)

