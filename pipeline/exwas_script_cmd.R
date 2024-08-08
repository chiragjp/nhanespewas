
#script to create "srun_me.sh" to execute the batch jobs
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
use_input_exposure_files <- TRUE # this is to filter the exposure inputs to parallelize the job: see 

to_do <- read_csv(ss_file) |> filter(pvarname == phenotype) |> group_by(evarname) |> summarize(total_n=sum(n), num_surveys = n()) |> mutate(pvarname = phenotype)
to_do <- to_do |> filter(num_surveys >=2, total_n >= sample_size_threshold)


main_out <- file.path(path_scripts_out, "srun_me.sh") 
cat("",file=main_out)
for(ii in 1:nrow(to_do)) {
  pvar <- to_do |> slice(ii) |> pull(pvarname)
  rcmd <- sprintf("Rscript %s -p %s -l %s -i %s -o %s", file.path(cmd_path) , pvar, ss_file, path_to_nhanes, path_out)
  
  if(use_sbatch) {
    cat(sprintf("sbatch --exclude=compute-f-17-[09-16] -o %s/%s.out -p short -t 11:59:00 -c 1 --mem=10G --wrap=\"%s\"\n", file.path(path_out), pvar, rcmd),
        file=main_out, append = T)
  } else {
    cat(sprintf("%s\n", rcmd), file=main_out, append = T)
  }
}

