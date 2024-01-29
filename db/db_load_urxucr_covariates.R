## insert urine creatinine covariates in demo tables (URXUCR)
##

library(tidyverse)
library(DBI)
con <- DBI::dbConnect(RSQLite::SQLite(), dbname='nhanes_012324.sqlite')


urinary_creatinine <- function(demo_table, series) {
  d_covariates <- NULL
  tab_name <- NULL
  if(series == 'A') {
    tab_name <- "LAB16"
  } else if(series == 'B') {
    tab_name <-  'L16_B'
  } else if(series == 'C') {
    tab_name <- 'L16_C'
  } else {
    tab_name <-sprintf("ALB_CR_%s", series)
  }
  d_covariates <- tbl(con, tab_name) |> select(c("SEQN", "URXUCR")) |> collect() |> rename(URXUCR_adj=URXUCR)
  demo_table |> left_join(d_covariates, by="SEQN")
}

table_descriptions <- tbl(con, "table_names_epcf") |> filter(epcf == 'c') |> collect()
tables <- vector(mode="list", length=nrow(table_descriptions))
for(ii in 1:nrow(table_descriptions)) {
  cat(sprintf("%s\n", table_descriptions$Data.File.Name[ii]))
  demo_tab <- tbl(con, table_descriptions$Data.File.Name[ii]) |> collect()
  series <- table_descriptions$series[ii]
  tables[[ii]] <- demo_tab |> urinary_creatinine(series)
  DBI::dbWriteTable(con, table_descriptions$Data.File.Name[ii], tables[[ii]], overwrite=T, append=F)
}
