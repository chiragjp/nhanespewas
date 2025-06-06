
#script to create "srun_me.sh" to execute the batch jobs for ubiome associations with phenotypes, such as BMI
#see exwas_script_cmd.R


library(getopt)
library(tidyverse)
library(logger)
source('db_paths.R')

path_out <- '../out'
ss_file <- '../select/sample_size_ubiome_y_112424.csv'
#path_to_db <-   '../db/nhanes_012324.sqlite' # '../nhanes_122322.sqlite'
path_to_nhanes <-'../db/nhanes_112824.sqlite'


to_do <- read_csv(ss_file) |> filter(grepl("relative$", evarname)) |> group_by(pvarname) |> count() |> ungroup()
to_do <- to_do |> filter(!grepl("^DX", pvarname)) # remove dexa


path_scripts_out <- '../scripts'
cmd_path <- '../pipeline/ubiome_phewas.R'
use_sbatch <- TRUE
use_input_exposure_files <- TRUE # this is to filter the exposure inputs to parallelize the job: see exwas_exposure_aux_files.R
num_input_exposure_files <- 5
use_quantile <- FALSE


main_out <- file.path(path_scripts_out, "srun_me.sh")
#print(nrow(to_do))

cat("",file=main_out)
for(ii in 1:nrow(to_do)) {
  pvar <- to_do |> slice(ii) |> pull(pvarname)
  rcmds <- vector("list", 1)
  if(use_input_exposure_files) {
    rcmds <- vector("list", num_input_exposure_files)
    for(i in 1:num_input_exposure_files) {
      e_file <- sprintf("%i.csv", i)
      rcmds[[i]] <- sprintf("Rscript %s -p %s -l %s -i %s -o %s -e %s -q %i", file.path(cmd_path) , pvar, ss_file, path_to_nhanes, path_out, file.path(path_scripts_out, e_file), use_quantile)
    }
  } else {
    rcmds[[1]] <- sprintf("Rscript %s -p %s -l %s -i %s -o %s -q %i", file.path(cmd_path) , pvar, ss_file, path_to_nhanes, path_out, use_quantile)
  }

  if(use_sbatch) {
    for(i in 1:length(rcmds)) {
      cat(sprintf("sbatch --exclude=compute-f-17-[09-16] -o %s/%s_%i.out -p short -t 11:59:00 -c 1 --mem=10G --wrap=\"%s\"\n", file.path(path_out), pvar, i, rcmds[i]),
          file=main_out, append = T)
    }
  } else {
    for(i in 1:length(rcmds)) {
      cat(sprintf("%s\n", rcmds[[i]]), file=main_out, append = T)
    }
  }
}

