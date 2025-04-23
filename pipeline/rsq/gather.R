
library(tidyverse)

path_to_rds <- "./out/m_50/"
rds_files <- dir(path_to_rds, pattern = "rds")

read_rds_out <- function(rds_file) {
  jack <- read_rds(rds_file)
  jack$rsq_summary
}

read_rds_summary <- function(rds_file) {
  jack <- read_rds(rds_file)
  jack$exposures_selected
}


rsq_summary <- map_df(file.path(path_to_rds,rds_files), read_rds_out)
rsq_summary <- rsq_summary |> unnest(cols=c(rsq_base))
exposure_selected <- map_df(file.path(path_to_rds,rds_files), read_rds_summary)

save(rsq_summary, exposure_selected, file='mvrsq_50.Rdata')
