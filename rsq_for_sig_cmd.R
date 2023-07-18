## Chirag Patel
## multivariate additive association between lifestyle exposures and a phenotype for a survey
## 05/15/23

## nutrition and biomarkers of exposure that are measured in ~10k
## usage: Rscript pe_rsq.R -p PHENOTYPE_NAME
## requires a summary_stats and a nhanes individual level database
## see: rsq.R

library(getopt)
library(this.path)
setwd(this.dir())
source('pe_rsq.R')
path_to_nhanes <- './nhanes_122322.sqlite'
path_to_summary_stats <- './pe_summary_stats.sqlite'

spec <- matrix(c(
  'pvarname', 'p', 1, "character"
), byrow=TRUE, ncol=4);
opt <- getopt(spec)

pvarname_to_query <- opt$pvarname

nhanes_con <- DBI::dbConnect(RSQLite::SQLite(), dbname=path_to_nhanes)
summary_stats_con <- DBI::dbConnect(RSQLite::SQLite(), dbname=path_to_summary_stats)



combined_dat <- get_individual_level_table(summary_stats_con, nhanes_con, pvarname_to_query)
log_info("E variable: for {pvarname_to_query}: {combined_dat$evars}")
log_info("Weight: for {pvarname_to_query}: {combined_dat$weight_to_use$weight_name}")
log_info("Sample size: for {pvarname_to_query}: {combined_dat$selected_data |> nrow() }")
svy_r2 <- svy_weighted_r2(pvarname_to_query,combined_dat$evars, adjustmentVariables, dat=combined_dat$selected_data, weight_name = combined_dat$weight_to_use$weight_name)


log_info("R2 base model for {pvarname_to_query}: {round(svy_r2$base$rsq, 3)}")
log_info("R2 exposure model for {pvarname_to_query}: {round(svy_r2$mve$rsq, 3)}")
log_info("R2 exposure model for {pvarname_to_query} attributable to E: {round(svy_r2$mve$rsq-svy_r2$base$rsq, 5)}")

dbDisconnect(nhanes_con)
dbDisconnect(summary_stats_con)

fileout <- sprintf('%s_r2.rds', pvarname_to_query)
x <- list(evarnames=combined_dat$evars, pvarname=pvarname_to_query, r2=svy_r2)
write_rds(x, file = fileout)
log_info("DONE")