# Chirag
# script to create script for slurm (HMS O2) for arbitrary X on Y analysis
# 3/1/24


library(tidyverse)
source('db_paths.R')
#path_to_nhanes <- './nhanes_122322.sqlite'

n_groups <- 100
ss_file <- '../select/sample_size_pp_category_032024.csv'
use_sbatch <- FALSE
path_out <- '../out/'

cmd_path <- '../pipeline/pp_by_tables_adj.R'
#cmd_path <- './demographic_makeup_by_tables.R'
path_scripts_out <- '../scripts'
sample_size_per_pair <- read_csv(ss_file)
sample_size_threshold <- 500
pairs <- sample_size_per_pair |> filter(n>=sample_size_threshold) |> select(x_table_name, y_table_name) |> unique()


## cut up pairs and put it into a .sh file in preparation for SLURM
## split by 200 each?
index <- rep(1:n_groups, floor(nrow(pairs)/n_groups))
index <- c(index, 1:( nrow(pairs)-length(index) ))
pairs <- pairs |> mutate(index=index)


cat_pairs <- function(table_pairs, sh_file_out) {
  cat("",file=sh_file_out)
  table_pairs |> transpose() |> walk(function(x) {
    cmd <- sprintf('Rscript %s -y %s -x %s -l %s -i %s -o %s \n', file.path(cmd_path), x$y_table_name, x$x_table_name, file.path('.', ss_file), file.path('.', path_to_nhanes), file.path('.', path_out))
    cat(cmd, file = sh_file_out, append = T)
  })
  
}

for(ii in 1:n_groups) {
  script_file_name <- file.path(path_scripts_out, sprintf('%i.sh', ii))
  cat_pairs(pairs |> filter(index == ii), script_file_name)
}


## create master bash file to submit to the sbatch queue
## exclude compute-f-17-*
## compute-f-17-[09-16]
main_out <- file.path(path_scripts_out, "srun_me.sh")
cat("",file=main_out)
for(ii in 1:n_groups) {
  if(use_sbatch) {
    cat(sprintf("sbatch --exclude=compute-f-17-[09-16] -o %i.out -p short -t 02:59:00 -c 1 --mem=2000 --wrap=\"./%i.sh\"\n", ii, ii), file=main_out, append = T)
  } else {
    cat(sprintf("./%i.sh\n", ii), file=main_out, append = T)
  }
}

