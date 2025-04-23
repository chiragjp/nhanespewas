## load height

library(tidyverse)
library(DBI)
con <- DBI::dbConnect(RSQLite::SQLite(), dbname='nhanes_031725.sqlite')

ht_join <- function(demo_table, series) {
  d_covariates <- NULL
  if(series == 'A') {
    d_covariates <- tbl(con, 'BMX') |> select(c("SEQN", "BMXHT")) |> collect() |> rename(BMXHT_adj=BMXHT)
  } else {
    d_covariates <- tbl(con, sprintf("BMX_%s", series)) |> select("SEQN", "BMXHT") |> collect() |> rename(BMXHT_adj = BMXHT)
  }
  demo_table |> left_join(d_covariates, by="SEQN")
}

table_descriptions <- tbl(con, "table_names_epcf") |> filter(epcf == 'c') |> collect() ## demo tables
tables <- vector(mode="list", length=nrow(table_descriptions))
for(ii in 1:nrow(table_descriptions)) {
  cat(sprintf("%s\n", table_descriptions$Data.File.Name[ii]))
  demo_tab <- tbl(con, table_descriptions$Data.File.Name[ii]) |> collect()
  ## now get diet table
  series <- table_descriptions$series[ii]
  tables[[ii]] <- demo_tab |> ht_join(series)
  DBI::dbWriteTable(con, table_descriptions$Data.File.Name[ii], tables[[ii]], overwrite=T, append=F)
}


