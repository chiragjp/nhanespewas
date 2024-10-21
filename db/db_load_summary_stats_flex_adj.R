# load summary stats to sqlite for the multiple adjustment scenarios
library(tidyverse)
library(DBI)

con <- DBI::dbConnect(RSQLite::SQLite(), dbname='./nhanes_012324.sqlite')
varnames <- tbl(con, "variable_names_epcf") |> collect()
expo_levels <- tbl(con, "e_variable_levels") |> collect()
tablenames <- tbl(con, "table_names_epcf") |> collect()
dbDisconnect(con)

path_to_summary <- '../pe_summary_0824/'
load(file.path(path_to_summary, 'gathered_by_series_0824.Rdata'))

domains <- read_csv('../select/variable_domains_ep_2.csv')

con <- DBI::dbConnect(RSQLite::SQLite(), dbname='./pe_summary_stats_08_2024.sqlite')
dbWriteTable(con, "variable_names_epcf", varnames, overwrite=T, append=F)
dbWriteTable(con, "table_names_epcf", tablenames, overwrite=T, append=F)
dbWriteTable(con, 'pe', pe, overwrite=T, append=F)
dbWriteTable(con, 'glanced', glanced, overwrite=T, append=F)
dbWriteTable(con, 'rsq', rsq, overwrite=T, append=F)
dbWriteTable(con, 'variable_domain', domains, overwrite=T, append=F)

adjusted_meta <- rbind(
  read_rds('../pe_summary_020424/pe_meta_model_1.rds') |> select(-meta_model) |> ungroup() |> rename(term_name=term) |> unnest(tidied) |> unnest(glanced),
  read_rds('../pe_summary_020424/pe_meta_model_2.rds')|> select(-meta_model) |> ungroup() |> rename(term_name=term) |> unnest(tidied) |> unnest(glanced),
  read_rds('../pe_summary_020424/pe_meta_model_3.rds') |> select(-meta_model) |> ungroup() |> rename(term_name=term) |> unnest(tidied) |> unnest(glanced),
  read_rds('../pe_summary_020424/pe_meta_model_4.rds')|> select(-meta_model) |> ungroup() |> rename(term_name=term) |> unnest(tidied) |> unnest(glanced),
  read_rds('../pe_summary_020424/pe_meta_model_5.rds')|> select(-meta_model) |> ungroup() |> rename(term_name=term) |> unnest(tidied) |> unnest(glanced),
  read_rds('../pe_summary_020424/pe_meta_model_6.rds')|> select(-meta_model) |> ungroup() |> rename(term_name=term) |> unnest(tidied) |> unnest(glanced),
  read_rds('../pe_summary_020424/pe_meta_model_7.rds')|> select(-meta_model) |> ungroup() |> rename(term_name=term) |> unnest(tidied) |> unnest(glanced),
  read_rds('../pe_summary_020424/pe_meta_model_8.rds')|> select(-meta_model) |> ungroup() |> rename(term_name=term) |> unnest(tidied) |> unnest(glanced),
  read_rds('../pe_summary_020424/pe_meta_model_9.rds')|> select(-meta_model) |> ungroup() |> rename(term_name=term) |> unnest(tidied) |> unnest(glanced)
)

adjusted_meta <- adjusted_meta |> rename(expo_name=term) |> mutate(expo_name=ifelse(is.na(expo_name), "expo", expo_name))


to_remove <- c(
  "LBXFOLSI",
  "LBXINSI",
  "LBXMCHSI",
  "LBXMMASI",
  "LBXRBFSI"
)

adjusted_meta <- adjusted_meta |> filter(!(phenotype %in% to_remove))
adjusted_meta <- adjusted_meta |> rename(pvarname = phenotype, evarname = exposure)

dbWriteTable(con, 'adjusted_meta', adjusted_meta, overwrite=T, append=F)
dbDisconnect(con)

remove(list=ls())
