## insert BMI as BMIADJ
##

library(tidyverse)
library(DBI)
con <- DBI::dbConnect(RSQLite::SQLite(), dbname='nhanes_031725.sqlite')


bmi_join <- function(demo_table, series) {
  ## if from series a or b, use DRXTKCAL, else DR1TKCAL
  ## search columns for TKCAL
  d_covariates <- NULL
  if(series == 'A') {
    d_covariates <- tbl(con, 'BMX') |> select(c("SEQN", "BMXBMI")) |> collect() |> rename(BMXBMI_adj=BMXBMI)
  } else {
    d_covariates <- tbl(con, sprintf("BMX_%s", series)) |> select("SEQN", "BMXBMI") |> collect() |> rename(BMXBMI_adj = BMXBMI)
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
  tables[[ii]] <- demo_tab |> bmi_join(series)
  DBI::dbWriteTable(con, table_descriptions$Data.File.Name[ii], tables[[ii]], overwrite=T, append=F)
}


