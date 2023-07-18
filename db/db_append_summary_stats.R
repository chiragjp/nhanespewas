## append summary stats to sqlite
library(tidyverse)
library(DBI)

con <- DBI::dbConnect(RSQLite::SQLite(), dbname='./pe_summary_stats.sqlite')
#vn <- tbl(con, "variable_names_epcf")
#tn <- tbl(con, "table_names_epcf")
#domains <- read_csv('./select/variable_domains_ep_2.csv')


path_to_summary <- './pe_summary_telo'
load(file.path(path_to_summary, 'pe_out_pscale.Rdata'))
dbWriteTable(con, 'pe', pe, overwrite=FALSE, append=TRUE)
dbWriteTable(con, 'glanced', glanced, overwrite=FALSE, append=TRUE)


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

dbWriteTable(con, 'adjusted_meta', adjusted_meta |> select(-meta_model) |> unnest(c(tidied, glanced)) |> ungroup(), overwrite=FALSE, append=TRUE)
dbWriteTable(con, 'unadjusted_meta', unadjusted_meta |> select(-meta_model) |> unnest(c(tidied, glanced)) |> ungroup(), overwrite=FALSE, append=TRUE)

dbDisconnect(con)
remove(list=ls())




