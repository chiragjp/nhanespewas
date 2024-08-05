
#exwas_script_cmd.R

path_out <- '.'
library(getopt)
library(tidyverse)
library(logger)
source('db_paths.R')

path_out <- '../out'
to_do <- read_csv(ss_file) |> group_by(pvarname) |> count()
path_scripts_out <- '../scripts'
cmd_path <- '../pipeline/exwas.R'
use_sbatch <- TRUE

main_out <- file.path(path_scripts_out, "srun_me.sh")
cat("",file=main_out)
for(ii in 1:nrow(to_do)) {
  pvar <- to_do |> slice(ii) |> pull(pvarname)
  rcmd <- sprintf("Rscript %s -p %s -l %s -i %s -o %s", file.path(cmd_path) , pvar, ss_file, path_to_nhanes, path_out)
  
  if(use_sbatch) {
    cat(sprintf("sbatch --exclude=compute-f-17-[09-16] -o %s%s.out -p short -t 10:59:00 -c 1 --mem=2000 --wrap=\"%s\"\n", file.path(path_out), pvar, rcmd),
        file=main_out, append = T)
  } else {
    cat(sprintf("%s\n", rcmd), file=main_out, append = T)
  }
}

