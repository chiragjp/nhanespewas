# load summary stats to sqlite
library(tidyverse)
library(DBI)

con <- DBI::dbConnect(RSQLite::SQLite(), dbname='./nhanes_122322.sqlite')
varnames <- tbl(con, "variable_names_epcf") |> collect()
expo_levels <- tbl(con, "e_variable_levels") |> collect()
tablenames <- tbl(con, "table_names_epcf") |> collect()
dbDisconnect(con)

#path_to_summary <- './pe_summary_051223/'

path_to_summary <- './pe_summary_060623/'

load(file.path(path_to_summary, 'pe_out_pscale.Rdata'))

domains <- read_csv('./select/variable_domains_ep_2.csv')

con <- DBI::dbConnect(RSQLite::SQLite(), dbname='./pe_summary_stats.sqlite')
dbWriteTable(con, "variable_names_epcf", varnames, overwrite=T, append=F)
dbWriteTable(con, "table_names_epcf", tablenames, overwrite=T, append=F)
dbWriteTable(con, 'pe', pe, overwrite=T, append=F)
dbWriteTable(con, 'glanced', glanced, overwrite=T, append=F)
dbWriteTable(con, 'variable_domain', domains, overwrite=T, append=F)

adjusted_meta <- rbind(
  read_rds(file.path(path_to_summary, 'pe_meta_adjusted_continuous.rds')) |> mutate(vartype="continuous"),
  read_rds(file.path(path_to_summary, 'pe_meta_adjusted_continuous-rank.rds')) |> mutate(vartype="continuous-rank"),
  read_rds(file.path(path_to_summary, 'pe_meta_adjusted_categorical.rds')) |> mutate(vartype="categorical")
)

adjusted_meta <- adjusted_meta |> select(-error) |> rename(expo_name=term) |> mutate(expo_name=ifelse(is.na(expo_name), "expo", expo_name))

unadjusted_meta <- rbind(
  read_rds(file.path(path_to_summary, 'pe_meta_unadjusted_continuous.rds')) |> mutate(vartype="continuous"),
  read_rds(file.path(path_to_summary, 'pe_meta_unadjusted_continuous-rank.rds')) |> mutate(vartype="continuous-rank"),
  read_rds(file.path(path_to_summary, 'pe_meta_unadjusted_categorical.rds')) |> mutate(vartype="categorical"),
)
  
unadjusted_meta <- unadjusted_meta |> select(-error) |> rename(expo_name=term) |> mutate(expo_name=ifelse(is.na(expo_name), "expo", expo_name))


to_remove <- c(
  "LBXFOLSI",
  "LBXINSI",
  "LBXMCHSI",
  "LBXMMASI",
  "LBXRBFSI"
)

adjusted_meta <- adjusted_meta |> filter(!(pvarname %in% to_remove))
unadjusted_meta <- unadjusted_meta |> filter(!(pvarname %in% to_remove))

dbWriteTable(con, 'adjusted_meta', adjusted_meta |> select(-meta_model) |> unnest(c(tidied, glanced)) |> ungroup(), overwrite=T, append=F)
dbWriteTable(con, 'unadjusted_meta', unadjusted_meta |> select(-meta_model) |> unnest(c(tidied, glanced)) |> ungroup(), overwrite=T, append=F)

dbDisconnect(con)
remove(list=ls())
