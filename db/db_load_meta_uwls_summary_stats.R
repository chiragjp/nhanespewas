library(tidyverse)
library(DBI)


#con <- DBI::dbConnect(RSQLite::SQLite(), dbname='./pe_summary_stats_02_2024-v2.sqlite')
con <- DBI::dbConnect(RSQLite::SQLite(), dbname='./pe_summary_stats_08_2024.sqlite')

adjusted_meta <- rbind(
  read_rds('../pe_summary_0824/pe_meta_model_uwls_1.rds')  |> ungroup() |> rename(term_name=term),
  read_rds('../pe_summary_0824/pe_meta_model_uwls_2.rds') |> ungroup() |> rename(term_name=term),
  read_rds('../pe_summary_0824/pe_meta_model_uwls_3.rds')  |> ungroup() |> rename(term_name=term),
  read_rds('../pe_summary_0824/pe_meta_model_uwls_4.rds') |> ungroup() |> rename(term_name=term),
  read_rds('../pe_summary_0824/pe_meta_model_uwls_5.rds') |> ungroup() |> rename(term_name=term),
  read_rds('../pe_summary_0824/pe_meta_model_uwls_6.rds') |> ungroup() |> rename(term_name=term),
  read_rds('../pe_summary_0824/pe_meta_model_uwls_7.rds') |> ungroup() |> rename(term_name=term),
  read_rds('../pe_summary_0824/pe_meta_model_uwls_8.rds')  |> ungroup() |> rename(term_name=term),
  read_rds('../pe_summary_0824/pe_meta_model_uwls_9.rds')  |> ungroup() |> rename(term_name=term)
)

adjusted_meta <- adjusted_meta |> select(-error) |> rename(expo_name=term_name) |> mutate(expo_name=ifelse(is.na(expo_name), "expo", expo_name))

to_remove <- c(
  "LBXFOLSI",
  "LBXINSI",
  "LBXMCHSI",
  "LBXMMASI",
  "LBXRBFSI"
)

adjusted_meta <- adjusted_meta |> filter(!(phenotype %in% to_remove))
adjusted_meta <- adjusted_meta |> rename(pvarname = phenotype, evarname = exposure)

dbWriteTable(con, 'adjusted_meta_uwls', adjusted_meta, overwrite=T, append=F)
dbDisconnect(con)

remove(list=ls())
