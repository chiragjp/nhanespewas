# Chirag
# script to create script for slurm (HMS O2)
# 4/1/23


library(tidyverse)
source('db_paths.R')
#path_to_nhanes <- './nhanes_122322.sqlite'

n_groups <- 100
#ss_file <- './select/sample_size_pe_040523.csv'  #opt$sample_size_pairs_list_file
#ss_file <- './select/sample_size_pe_category_041823.csv'  #opt$sample_size_pairs_list_file
#ss_file <- './select/sample_size_pe_category_051323.csv'  #opt$sample_size_pairs_list_file
#ss_file <- './select/sample_size_pe_category_telo_052823.csv'
ss_file <- '../select/sample_size_pe_category_060623.csv'
#ss_file <- './select/sample_size_alphape.csv'
use_sbatch <- FALSE
path_out <- '../out/'
#cmd_path <- './pe_by_tables.R'
cmd_path <- './demographic_makeup_by_tables.R'
path_scripts_out <- '../scripts'
sample_size_per_pair <- read_csv(ss_file)
sample_size_threshold <- 500
pairs <- sample_size_per_pair |> filter(n>=sample_size_threshold) |> select(e_table_name, p_table_name) |> unique()



## cut up pairs and put it into a .sh file in preparation for SLURM
## split by 200 each?
index <- rep(1:n_groups, floor(nrow(pairs)/n_groups))
index <- c(index, 1:( nrow(pairs)-length(index) ))
pairs <- pairs |> mutate(index=index)


cat_pairs <- function(table_pairs, sh_file_out) {
  cat("",file=sh_file_out)
  table_pairs |> transpose() |> walk(function(x) {
    cmd <- sprintf('Rscript %s -p %s -e %s -l %s -i %s -o %s \n', file.path('..', cmd_path), x$p_table_name, x$e_table_name, file.path('.', ss_file), file.path('.', path_to_nhanes), file.path('.', path_out))
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

