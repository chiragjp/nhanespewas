## Chirag
## 3/31/23
## load in the alpha diversity to the database

library(tidyverse)
library(DBI)

alpha_rb <- read_table('./data/oral_microbiome/dada2rb-alpha.txt', na='.')
alpha_rsv <- read_table('./data/oral_microbiome/dada2rsv-alpha.txt', na='.')
alphas <- alpha_rb |> left_join(alpha_rsv)

## create mean and SD alpha_rb and mean and SD alpha_rsv

# RB_ObservedOTUs
# RB_FaPhyloDiv
# RB_ShanWienDiv
# RB_InverseSimpson

# RSV_ObservedOTUs
# RSV_FaPhyloDiv
# RSV_ShanWienDiv
# RSV_InverseSimpson

start_strings <- c(
  "RB_ObservedOTUs",
  "RB_FaPhyloDiv",
  "RB_ShanWienDiv",
  "RB_InverseSimpson",
  "RSV_ObservedOTUs",
  "RSV_FaPhyloDiv",
  "RSV_ShanWienDiv",
  "RSV_InverseSimpson"
)


alphas <- alphas |> rowwise()
for(strn in start_strings) {
  print(strn)
  new_name_mean <- str_glue("{strn}Mean")
  new_name_sd <- str_glue("{strn}SD")
  alphas <- alphas |>  mutate("{new_name_mean}" := mean(c_across(starts_with(strn)), na.rm=T),
                                "{new_name_sd}" := sd(c_across(starts_with(strn)), na.rm=T))  
}

alphas <- alphas |> ungroup()

## copy to db
path_to_db <- './nhanes_122322.sqlite'
con <- DBI::dbConnect(RSQLite::SQLite(), dbname=path_to_db)
table_description <- tbl(con, "table_names_epcf")
variable_description <- tbl(con, "variable_names_epcf")
## get demo 2009 and 2011
# F and G
demo_f <- tbl(con, "DEMO_F") |> select(SEQN) |> collect()
demo_g <- tbl(con, "DEMO_G") |> select(SEQN) |> collect()

# Table Names: ALPHADIVERSITY_G, ALPHADIVERSITY_F
alpha_f <- demo_f |> inner_join(alphas, by="SEQN")
alpha_g <- demo_g |> inner_join(alphas, by="SEQN")

## insert into the con table
dbWriteTable(con, "ALPHADIVERSITY_F" , alpha_f, overwrite=T, append=F)
dbWriteTable(con, "ALPHADIVERSITY_G" , alpha_g, overwrite=T, append=F)

## now write to the table_description and variable_description
table_descs <- tibble(
  Data.File.Name = c("ALPHADIVERSITY_F", "ALPHADIVERSITY_G"),
  component = "LABORATORY",
  series = c("F", "G"),
  epcf = 'p'
)

new_table_description <- table_description |> collect() |> rbind(table_descs)
variable_names <- colnames(alphas)

variable_descs_f <- tibble(
  Variable.Name = variable_names,
  Variable.Description= variable_names,
  Data.File.Name = 'ALPHADIVERSITY_F',
  Data.File.Description = "Oral Microbiome Alpha Diversity",
  Begin.Year = 2009,
  EndYear = 2010,
  Component = "Laboratory",
  Use.Constraints = "None"
)

variable_descs_g <- tibble(
  Variable.Name = variable_names,
  Variable.Description= variable_names,
  Data.File.Name = 'ALPHADIVERSITY_G',
  Data.File.Description = "Oral Microbiome Alpha Diversity",
  Begin.Year = 2011,
  EndYear = 2012,
  Component = "Laboratory", 
  Use.Constraints = "None"
)

new_variable_description <- variable_description |> collect() |> rbind(variable_descs_f, variable_descs_g)

dbWriteTable(con, "table_names_epcf", new_table_description, append=F, overwrite=T)
dbWriteTable(con, "variable_names_epcf", new_variable_description, append=F, overwrite=T)

dbDisconnect(con)

variable_descs <- rbind(variable_descs_f, variable_descs_g)

variable_descs |> filter(grepl("Mean", Variable.Name)) |> rbind(
  variable_descs |> filter(grepl("SD", Variable.Name))
) |> write_csv(file="./select/select_ubiome_pheno_variables.csv")

