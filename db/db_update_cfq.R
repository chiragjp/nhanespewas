## Chirag
## 7/27/25
## preapre the CFQ data

library(tidyverse)
library(DBI)

## copy to db
path_to_db <- './nhanes_031725.sqlite'
con <- DBI::dbConnect(RSQLite::SQLite(), dbname=path_to_db)
table_description <- tbl(con, "table_names_epcf")
variable_description <- tbl(con, "variable_names_epcf")

## now write to the table_description and variable_description
table_descs <- tibble(
  Data.File.Name = c("CFQ", "CFQ_B", "CFQ_G", "CFQ_H"),
  component = "QUESTIONNAIRE",
  series = c("A", "B", "G", "H"),
  epcf = 'p'
)


cfq <- tbl(con, "CFQ") |> collect() |> mutate(CFDDS=CFDRIGHT)
cfq_b <- tbl(con, "CFQ_B" ) |> collect() |> mutate(CFDDS=CFDRIGHT)


new_table_description <- table_description |> collect() |> rbind(table_descs)

variable_descs_ab <- tibble(
  Variable.Name = c("CFDRIGHT", "CFDRIGHT",  "CFDDS", "CFDDS"),
  Variable.Description= "Score: number correct",
  Data.File.Name = c("CFQ", "CFQ_B", "CFQ", "CFQ_B"),
  Data.File.Description = "Cognitive Functioning",
  Begin.Year = c(1999, 2001, 1999, 2001),
  EndYear = c(2000, 2002, 2000, 2002),
  Component = "Questionnaire",
  Use.Constraints = "None"
)

variable_descs_gh <- tibble(
  Variable.Name = "CFDDS",
  Variable.Description= "Digit Symbol: Score",
  Data.File.Name = c("CFQ_G", "CFQ_H"),
  Data.File.Description = "Cognitive Functioning",
  Begin.Year = c(2011, 2013),
  EndYear = c(2012, 2014),
  Component = "Questionnaire",
  Use.Constraints = "None"
)

new_variable_description <-variable_description |> collect() |>rbind(variable_descs_ab, variable_descs_gh)

dbWriteTable(con, "CFQ", cfq, append=F, overwrite=T)
dbWriteTable(con, "CFQ_B", cfq_b, append=F, overwrite=T)
dbWriteTable(con, "table_names_epcf", new_table_description, append=F, overwrite=T)
dbWriteTable(con, "variable_names_epcf", new_variable_description, append=F, overwrite=T)

dbDisconnect(con)


rbind(variable_descs_ab, variable_descs_gh) |> filter(Variable.Name == 'CFDDS') |> write_csv("../select/select_cfq_pheno_variables.csv")
