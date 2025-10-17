## DNAM Age process data file
# August 30 2025

library(tidyverse)
library(sas7bdat)
sasfile <- sas7bdat::read.sas7bdat('./DNAm Age/dnmepi.sas7bdat')

demo <- read_rds('./data/nhanes_DEMO_all_1222222.rds')
demo_a <- demo$tables[[1]]$DEMO |> select(SEQN, SDDSRVYR)
demo_b <- demo$tables[[2]]$DEMO_B |> select(SEQN, SDDSRVYR)

dnam_a <- sasfile |> inner_join(demo_a) |> write_rds("./DNAm Age/DNAM_A.rds")
dnam_b <- sasfile |> inner_join(demo_b) |> write_rds("./DNAm Age/DNAM_B.rds")

