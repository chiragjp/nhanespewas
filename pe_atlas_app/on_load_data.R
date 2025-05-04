## Chirag
## create on load tables to save time

library(dbplyr)
library(tidyverse)
library(DBI)
library(pool)
library(fst)
pool <- dbPool(drv = RSQLite::SQLite(), dbname='./pe_shiny_2.sqlite')
##

p_variables <- tbl(pool, 'adjusted_meta_2') |> group_by(pvarname) |> count() |> left_join(tbl(pool, "p_variable_domain"), by=c("pvarname"="Variable.Name")) |> collect()
p_variables <- p_variables |> mutate(cat_subcat=ifelse(!is.na(psubcategory), paste(pcategory, psubcategory, sep="-"), pcategory ))
p_variables <- p_variables |> mutate(pvardesc_selector = sprintf("%s-(%s)", pvardesc, pvarname))

e_variables <- tbl(pool, 'adjusted_meta_2') |> group_by(evarname) |> count() |> left_join(tbl(pool, "e_variable_domain"), by=c("evarname"="Variable.Name")) |> collect()
e_variables <- e_variables |> mutate(cat_subcat=ifelse(!is.na(esubcategory), paste(ecategory, esubcategory, sep="-"), ecategory ))
e_variables <- e_variables |> mutate(evardesc_selector = sprintf("%s-(%s)", evardesc, evarname))
e_category_strs <- e_variables |> ungroup() |> select(ecategory, esubcategory) |> group_by(ecategory, esubcategory) |> count()
e_category_strs <- e_category_strs |> ungroup() |> mutate(cat_subcat=ifelse(!is.na(esubcategory), paste(ecategory, esubcategory, sep="-"), ecategory ))


#purrr::walk(DBI::dbListTables(pool), \(tbl) {
#  df <- DBI::dbReadTable(pool)
#  write_fst(df, file.path(paste0(tbl, ".fst")), compress = 50)
#})

tables <- dbListTables(pool)
# 4. Loop over each table: read it and write to .fst
for(tbl in tables) {
  message("Dumping table: ", tbl)
  df <- dbReadTable(pool, tbl)
  write_fst(
    df,
    path    = file.path(paste0(tbl, ".fst")),
    compress = 50        # adjust 0–100: higher = smaller file, slightly slower read
  )
}

poolClose(pool)

write_fst(p_variables, "p_variables.fst")
write_fst(e_variables, "e_variables.fst")
write_fst(e_category_strs, "e_category_strs.fst")


