## microbiome data
## this was released on 11/27/24
## chirag 
## 11/28/24

library(tidyverse)
library(DBI)


directory_to_ubiome <- '~/Dropbox (HMS)/projects/nhanespewas/download/ubiome/'
count_files <- dir(directory_to_ubiome, pattern = "*count")
relative_files <- dir(directory_to_ubiome, pattern = "*relative")
annotate_files <- dir(directory_to_ubiome, pattern = "*annotate")
## rb and rsv files

files <- rbind(
  tibble(filename=count_files, type="count"),
  tibble(filename=relative_files, type="relative")
)


annotation_rb <- read_tsv(file.path(directory_to_ubiome, "dada2rb-taxonomy-annotate.txt"))
annotation_rsv <- read_tsv(file.path(directory_to_ubiome, "dada2rsv-taxonomy-annotate.txt"))

annotation_rb <- annotation_rb |> rename(taxonomy_silva_123 = `Taxonomy in SILVA v123`)
annotation_rsv <- annotation_rsv |> rename(taxonomy_silva_123 = `Taxonomy in SILVA v123`)


## new path here
path_to_db <- './nhanes_112824.sqlite'
con <- DBI::dbConnect(RSQLite::SQLite(), dbname=path_to_db)
table_description <- tbl(con, "table_names_epcf")
variable_description <- tbl(con, "variable_names_epcf")

## get demo 2009 and 2011
# F and G
demo_f <- tbl(con, "DEMO_F") |> select(SEQN) |> collect()
demo_g <- tbl(con, "DEMO_G") |> select(SEQN) |> collect()

# merge and load
#table_list <- vector("list", length=length(files))
table_info <- vector("list", nrow(files))
for(ii in 1:nrow(files)) {
  #print(files[[ii]])
  f <- files |> slice(ii) 
  fname <- f |> pull(filename)
  vartype <- f |> pull(type)
  tab <- read_tsv(file.path(directory_to_ubiome, fname))
  tab <- tab |> rename_with(~paste0(.x, "_", vartype), !starts_with("SEQN"))
  tab_f <- tab |> inner_join(demo_f)
  tab_g <- tab |> inner_join(demo_g)
  table_name_f <- sprintf("%s_F", toupper(sub(".txt","" , gsub("-", "_", fname))))
  table_name_g <- sprintf("%s_G", toupper(sub(".txt","" , gsub("-", "_", fname))))
  table_info[[ii]] <- rbind(
    tibble(table_name=table_name_g, col=colnames(tab_g), series = "G", ubiome_type=vartype),
    tibble(table_name=table_name_f, col=colnames(tab_f), series= "F", ubiome_type=vartype))
  print(table_name_f)
  dbWriteTable(con, table_name_f, tab_f, overwrite=T,append=F)
  print(table_name_g)
  dbWriteTable(con, table_name_g, tab_g, overwrite=T,append=F)
}

table_info <- bind_rows(table_info) |> filter(col != "SEQN")

table_descs <- rbind(
  tibble(
    Data.File.Name = table_info |> filter(series == 'G') |> pull(table_name) |> unique(),
    component = "LABORATORY",
    series = c("G"),
    epcf = 'p'
  ), 
  tibble(
    Data.File.Name = table_info |> filter(series == 'F') |> pull(table_name) |> unique(),
    component = "LABORATORY",
    series = c("F"),
    epcf = 'p'
  )
)

new_table_description <- table_description |> select(-5) |> collect() |> rbind(table_descs)

variable_descs <- table_info |> 
  rename(Data.File.Name = table_name, Variable.Name = col) |> 
  mutate(Data.File.Description = Data.File.Name, Component="Laboratory", Use.Constraints = "None") |>
  mutate(Begin.Year = case_when(
    series == "F" ~ 2009,
    series == "G" ~ 2011,
    TRUE ~ NA_integer_
  )) |> 
  mutate(EndYear = Begin.Year + 1)


variable_descs <- variable_descs |> mutate(k1=str_split_i(Variable.Name, "_", 1), 
                                           k2=str_split_i(Variable.Name, "_", 2),
                                           Orig.Name=paste0(k1, "_", k2)) |> select(-c(k1, k2))
variable_descs <- variable_descs |> left_join(annotation_rb |> rbind(annotation_rsv), by=c("Orig.Name"="Name"))
variable_descs <- variable_descs |> rename(Variable.Description = taxonomy_silva_123)

new_variable_description <- variable_description |> collect()  |> rbind(variable_descs |> select (-series, -Orig.Name, -ubiome_type))

dbWriteTable(con, "table_names_epcf", new_table_description, append=F, overwrite=T)
dbWriteTable(con, "variable_names_epcf", new_variable_description, append=F, overwrite=T)

dbDisconnect(con)


write_csv(variable_descs |> select (-series, -Orig.Name, -ubiome_type), file='../select/select_ubiome_count_pheno_variables.csv')