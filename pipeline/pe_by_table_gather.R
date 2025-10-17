# Gather all .rds files in a folder ('out')

library(tidyverse)
library(getopt)
spec <- matrix(c(
  'path_to_gathered', 'o', 1, "character",
  'input_directory', 'i', 1, "character"
), byrow=TRUE, ncol=4);
opt <- getopt(spec)


path_to_out <- opt$input_directory
rds_files <- dir(path_to_out, pattern = '*.rds')
pe <- rds_files |> map_df(function(x) {
  dat <- read_rds(file.path(path_to_out, x))
  dat$pe_tidied
})

glanced <- rds_files |> map_df(function(x) {
  dat <- read_rds(file.path(path_to_out, x))
  dat$pe_glanced
})

rsq <- rds_files |> map_df(function(x) {
  dat <- read_rds(file.path(path_to_out, x))
  dat$rsq
})

reg_term_test_tidied <- rds_files |> map_df(function(x) {
  dat <- read_rds(file.path(path_to_out, x))
  dat$reg_term_test_tidied
})

file_out <- opt$path_to_gathered
save(pe, glanced, rsq, reg_term_test_tidied, file = file_out)

