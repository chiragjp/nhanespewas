
library(tidyverse)
library(DBI)

con <- DBI::dbConnect(RSQLite::SQLite(), dbname='./nhanes_012324.sqlite')
varnames <- tbl(con, "variable_names_epcf") |> collect()
tablenames <- tbl(con, "table_names_epcf") |> collect()
dbDisconnect(con)

path_to_summary <- '../pp_summary_032124'
load(file.path(path_to_summary, 'pp_summary_032024.Rdata'))

domains <- read_csv('../select/variable_domains_ep_2.csv')

con <- DBI::dbConnect(RSQLite::SQLite(), dbname='./pp_summary_stats_03_2024.sqlite')
dbWriteTable(con, "variable_names_epcf", varnames, overwrite=T, append=F)
dbWriteTable(con, "table_names_epcf", tablenames, overwrite=T, append=F)
dbWriteTable(con, 'pe', pe, overwrite=T, append=F)
dbWriteTable(con, 'glanced', glanced, overwrite=T, append=F)
dbWriteTable(con, 'rsq', rsq, overwrite=T, append=F)
dbWriteTable(con, 'variable_domain', domains, overwrite=T, append=F)

adjusted_meta <- rbind(
  read_rds(file.path(path_to_summary, 'pp_meta_model_1.rds')) |> select(-meta_model) |> ungroup() |> rename(term_name=term) |> unnest(tidied) |> unnest(glanced),
  read_rds(file.path(path_to_summary, 'pp_meta_model_2.rds')) |> select(-meta_model) |> ungroup() |> rename(term_name=term) |> unnest(tidied) |> unnest(glanced)
)



to_remove <- c(
  "LBXFOLSI",
  "LBXINSI",
  "LBXMCHSI",
  "LBXMMASI",
  "LBXRBFSI"
)

adjusted_meta <- adjusted_meta |> filter(!(yvar %in% to_remove))
adjusted_meta <- adjusted_meta |> filter(!(xvar %in% to_remove))
adjusted_meta <- adjusted_meta |> rename(xvarname = xvar, yvarname = yvar)

dbWriteTable(con, 'adjusted_meta', adjusted_meta, overwrite=T, append=F)
dbDisconnect(con)

remove(list=ls())
