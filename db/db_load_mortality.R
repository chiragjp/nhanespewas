# Chirag J Patel
# 2/4/2024
# load mortality data

library(tidyverse)
library(DBI)

dbname <- "nhanes_012324.sqlite"

mortality_data <- read_rds("../download/mortality_data/mortality_data_1999-2018_2019.rds")
con <- dbConnect(RSQLite::SQLite(), dbname)

dbWriteTable(con, "mortality", mortality_data, append=F, overwrite=T)

dbDisconnect(con)
