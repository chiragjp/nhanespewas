
#script to create "srun_me.sh" to execute the batch jobs
#exwas_script_cmd.R

path_out <- '.'
library(getopt)
library(tidyverse)
library(logger)
source('db_paths.R')

path_out <- '../out'
#ss_file <- '../select/sample_size_pe_category_060623.csv'
ss_file <- '../select/sample_size_pe_category_0824.csv'
to_do <- read_csv(ss_file) |> group_by(pvarname) |> count() |> ungroup()
path_scripts_out <- '../scripts'
cmd_path <- '../pipeline/exwas.R'
use_sbatch <- T
use_input_exposure_files <- T # this is to filter the exposure inputs to parallelize the job: see exwas_exposure_aux_files.R
num_input_exposure_files <- 10
use_quantile <- 0
interact_with <- "RIDAGEYR"
adjustment_scenario <- "age_sex_ethnicity_income_education"

main_out <- file.path(path_scripts_out, "srun_me.sh")
#print(nrow(to_do))

sprintf_interact_param <- function(cmd, interact_with) {
  if(is.null(interact_with)) {
    return(cmd)
  }
  cmd <- sprintf(" %s -w %s", cmd, interact_with)
}

sprintf_adjustment_param <- function(cmd, scenario) {
  if(is.null(scenario)) {
    return(cmd)
  }
  cmd <- sprintf("%s -s %s", cmd, scenario)
}


cat("",file=main_out)
for(ii in 1:nrow(to_do)) {
  pvar <- to_do |> slice(ii) |> pull(pvarname)
  rcmds <- vector("list", 1)
  if(use_input_exposure_files) {
    rcmds <- vector("list", num_input_exposure_files)
    for(i in 1:num_input_exposure_files) {
      e_file <- sprintf("%i.csv", i)
      rcmds[[i]] <- sprintf("Rscript %s -p %s -l %s -i %s -o %s -e %s -q %i", file.path(cmd_path) , pvar, ss_file, path_to_nhanes, path_out, file.path(path_scripts_out, e_file), use_quantile) |>
        sprintf_interact_param(interact_with=interact_with) |> sprintf_adjustment_param(adjustment_scenario)
    }
  } else {
    rcmds[[1]] <- sprintf("Rscript %s -p %s -l %s -i %s -o %s -q %i", file.path(cmd_path) , pvar, ss_file, path_to_nhanes, path_out, use_quantile) |>
      sprintf_interact_param(interact_with=interact_with) |> sprintf_adjustment_param(adjustment_scenario)
  }

  if(use_sbatch) {
    for(i in 1:length(rcmds)) {
      cat(sprintf("sbatch -o %s/%s_%i.out -p short -t 11:59:00 -c 1 --mem=10G --wrap=\"%s\"\n", file.path(path_out), pvar, i, rcmds[[i]]),
          file=main_out, append = T)
      #cat(sprintf("sbatch --exclude=compute-f-17-[09-16] -o %s/%s_%i.out -p short -t 11:59:00 -c 1 --mem=10G --wrap=\"%s\"\n", file.path(path_out), pvar, i, rcmds[i]),
      #    file=main_out, append = T)
    }
  } else {
    for(i in 1:length(rcmds)) {
      cat(sprintf("%s\n", rcmds[[i]]), file=main_out, append = T)
    }
  }
}

