## insert diet specific covariates in demo tables (total caloric intake)
##

library(tidyverse)
library(DBI)
con <- DBI::dbConnect(RSQLite::SQLite(), dbname='nhanes_012324.sqlite')


total_calorie <- function(demo_table, series) {
  ## if from series a or b, use DRXTKCAL, else DR1TKCAL
  ## search columns for TKCAL
  d_covariates <- NULL
  if(series == 'A') {
    d_covariates <- tbl(con, 'DRXTOT') |> select(c("SEQN", "DRXTKCAL")) |> collect() |> rename(DRXTKCAL_adj=DRXTKCAL)
  } else if(series == 'B') {
    d_covariates <- tbl(con, 'DRXTOT_B') |> select(c("SEQN", "DRXTKCAL")) |> collect() |> rename(DRXTKCAL_adj = DRXTKCAL)
  } else {
    d_covariates_day1 <- tbl(con, sprintf("DR1TOT_%s", series)) |> select("SEQN", "DR1TKCAL") |> collect() |> rename(DR1TKCAL_adj = DR1TKCAL)
    d_covariates_day2 <- tbl(con, sprintf("DR2TOT_%s", series)) |> select("SEQN", "DR2TKCAL") |> collect() |> rename(DR2TKCAL_adj = DR2TKCAL)
    d_covariates <- d_covariates_day1 |> full_join(d_covariates_day2, by="SEQN")
  }
  demo_table |> left_join(d_covariates, by="SEQN")
}

table_descriptions <- tbl(con, "table_names_epcf") |> filter(epcf == 'c') |> collect()
tables <- vector(mode="list", length=nrow(table_descriptions))
for(ii in 1:nrow(table_descriptions)) {
  cat(sprintf("%s\n", table_descriptions$Data.File.Name[ii]))
  demo_tab <- tbl(con, table_descriptions$Data.File.Name[ii]) |> collect()
  ## now get diet table
  series <- table_descriptions$series[ii]
  tables[[ii]] <- demo_tab |> total_calorie(series)
  DBI::dbWriteTable(con, table_descriptions$Data.File.Name[ii], tables[[ii]], overwrite=T, append=F)
}
