# Chirag
# December 2024
# ubiome_outcome.R: perform ubiome-wide study associating clinical outcome with the microbiome
# modified for ubiome relative abundance (12/24)

library(getopt)
library(tidyverse)
library(logger)
library(tools)
devtools::load_all("..")

TEST <- FALSE
spec <- matrix(c(
  'phenotype', 'p', 1, "character",
  'sample_size_pairs_list_file', 'l', 1, "character", # file that lists the sample sizes for each pair
  'path_to_db', 'i', 1, "character",
  'path_out', 'o', 1, "character",
  'exposures', 'e', 2, "character",
  'use_quantile', 'q', 2, "numeric"
), byrow=TRUE, ncol=4);
opt <- getopt(spec)

sample_size_threshold <- 500



########### debug stuff
exposure <- "RSV_genus29_relative"
phenotype <- "TOOTH_DECAY_OHAROCDT"
ss_file <- '../select/sample_size_ubiome_y_112424.csv'
path_to_db <-'../db/nhanes_112824.sqlite'
path_out <- '.'
use_quantile <- FALSE
if(!TEST) {
  phenotype <- opt$phenotype
  ss_file <- opt$sample_size_pairs_list_file
  path_to_db <- opt$path_to_db
  path_out <- opt$path_out
  use_quantile <- ifelse(is.null(opt$use_quantile), 0, opt$use_quantile)
}


############### end DEBUG

con <- DBI::dbConnect(RSQLite::SQLite(), dbname=path_to_db)
#to_do <- read_csv(ss_file) |> filter(pvarname == phenotype) |> group_by(evarname) |> summarize(total_n=sum(n), num_surveys = n()) |> mutate(pvarname = phenotype)
#to_do <- to_do |> filter(num_surveys >=2, total_n >= sample_size_threshold) # no need to threshold
to_do <- read_csv(ss_file) |>
  filter(pvarname == 'BMXBMI', grepl("genus", evarname), grepl("relative$", evarname), grepl("^RB", evarname)) |> group_by(evarname) |>
  summarize(total_n=sum(n), num_surveys = n()) |> mutate(pvarname = phenotype)

exposure_name_sample <- (to_do |> slice_head() |> pull(evarname))[1]

ptables <- get_table_names_for_varname(con, varname = phenotype) |> rename(p_name = Data.File.Name)
etables <- get_table_names_for_varname(con, varname = exposure_name_sample) |> rename(e_name = Data.File.Name)

table_set <- ptables |> inner_join(etables, by = "Begin.Year")
tab_obj <- get_x_y_tables_as_list(con,table_set$p_name,table_set$e_name)
## weight
tab_obj <- figure_out_multiyear_weight(tab_obj)

if(TEST) {
  to_do <- to_do |> filter(evarname == exposure)
}

if(!is.null(opt$exposures)) {
  exposures <- read_csv(opt$exposures, col_names=F)
  to_do <- to_do |> inner_join(exposures |> rename(evarname=X1), by="evarname")
}

pe_safely <- safely(nhanespewas::logistic_e_flex_adjust_by_table)
#pe_safely <- (nhanespewas::logistic_e_flex_adjust_by_table)
tidied <- vector("list", length = nrow(to_do))
glanced <- vector("list", length = nrow(to_do))

log_info("Process ID: {Sys.getpid()}")
log_info("Number of Exposures: {nrow(to_do)}")

if(nrow(to_do) == 0) {
  done <- dbDisconnect(con)
  log_info("0 pairs, quitting")
  log_info("Done for the phenotype { phenotype } ")
  stop("No pairs to execute")
}

adjustment_model_for_e <- adjustment_models |> filter(scenario == "base" | scenario == "age_sex_ethnicity_income_education")# limit the adjustment scenarios

N <- nrow(to_do)
models <- vector("list", length=N)

for(ii in 1:N) {
  rw <- to_do |> slice(ii)
  log_info("{ii} out of {nrow(to_do)}; expo: {rw$evarname}; pheno: {rw$pvarname} ")
  mod <- NULL
  if(use_quantile==1) {
    mod <- pe_safely(tab_obj,rw$pvarname, rw$evarname, adjustment_model_for_e, logxform_e=F, scale_e=F,
                     quantile_expo=c(0, .25, .5, .75, 1), exposure_levels=NULL)
  } else {
    mod <- pe_safely(tab_obj,rw$pvarname, rw$evarname, adjustment_model_for_e, logxform_e=F, scale_e=T,
                     quantile_expo=NULL, exposure_levels=NULL)
  }
  models[[ii]] <- mod ## each iteration contains a pair, and a data struct of all the models
}



tidied_result <- function(models, aggregate_base_model = F) {
  final_result <- map_dfr(models, ~ {
    # Check for the presence of an error
    if (is.null(.x$error)) {
      # If no error, proceed to extract and bind tibbles
      lcl <- .x$result
      m_index <- 1:length(.x$result$models)
      to_do <- .x$result$models
      if(aggregate_base_model) {
        to_do <- .x$result$base_models
      }

      map2_dfr(to_do, m_index, ~ { ## loops over all adjustment scenarios
        if(is.null(.x)) {
          return(tibble())
        }
        tidied_with_key <- mutate(.x$tidied,
                                  model_number = .y,
                                  series = paste(lcl$series, collapse=";"),
                                  exposure = lcl$exposure,
                                  phenotype= lcl$phenotype,
                                  #log_p= lcl$log_p,
                                  log_e = lcl$log_e,
                                  #scaled_p = lcl$scaled_p,
                                  scaled_e = lcl$scaled_e
        )
        return(tidied_with_key)
      })
    } else {
      # Return an empty tibble if there's an error (or handle it differently as needed)
      tibble()
    }
  })
  final_result
}

glanced_result <-function(models, aggregate_base_model=F) {
  final_result <- map_dfr(models, ~ {
    # Check for the presence of an error
    if (is.null(.x$error)) {
      # If no error, proceed to extract and bind tibbles
      lcl <- .x$result
      m_index <- 1:length(.x$result$models)

      to_do <- .x$result$models
      if(aggregate_base_model) {
        to_do <- .x$result$base_models
      }

      map2_dfr(to_do, m_index, ~ { ## loops over all adjustment scenarios
        if(is.null(.x)) {
          return(tibble())
        }
        tidied_with_key <- mutate(.x$glanced,
                                  model_number = .y,
                                  series = paste(lcl$series, collapse=";"),
                                  exposure = lcl$exposure,
                                  phenotype= lcl$phenotype,
                                  #log_p= lcl$log_p,
                                  log_e = lcl$log_e,
                                  #scaled_p = lcl$scaled_p,
                                  scaled_e = lcl$scaled_e
        )
        return(tidied_with_key)
      })
    } else {
      # Return an empty tibble if there's an error (or handle it differently as needed)
      tibble()
    }
  })
  final_result
}

r2_result <-function(models, aggregate_base_model=F) {
  final_result <- map_dfr(models, ~ {
    # Check for the presence of an error
    if (is.null(.x$error)) {
      # If no error, proceed to extract and bind tibbles
      lcl <- .x$result
      m_index <- 1:length(.x$result$models)

      to_do <- .x$result$models
      if(aggregate_base_model) {
        to_do <- .x$result$base_models
      }

      map2_dfr(to_do, m_index, ~ { ## loops over all adjustment scenarios
        if(is.null(.x)) {
          return(tibble())
        }
        tidied_with_key <- mutate(.x$r2 |> as_tibble(),
                                  model_number = .y,
                                  series = paste(lcl$series, collapse=";"),
                                  exposure = lcl$exposure,
                                  phenotype= lcl$phenotype,
                                  #log_p= lcl$log_p,
                                  log_e = lcl$log_e,
                                  #scaled_p = lcl$scaled_p,
                                  scaled_e = lcl$scaled_e
        )
        return(tidied_with_key)
      })
    } else {
      # Return an empty tibble if there's an error (or handle it differently as needed)
      tibble()
    }
  })
  final_result
}


tidied <- rbind(
  tidied_result(models, F) |> mutate(aggregate_base_model = F),
  tidied_result(models, T) |> mutate(aggregate_base_model = T)
) |> mutate(phenotype = phenotype)
glanced <- rbind(
  glanced_result(models, F) |> mutate(aggregate_base_model = F),
  glanced_result(models, T) |> mutate(aggregate_base_model = T)
) |> mutate(phenotype = phenotype )

#rsq <- rbind(
#  r2_result(models, F) |> mutate(aggregate_base_model = F),
#  r2_result(models, T) |> mutate(aggregate_base_model = T)
#) |> mutate(phenotype=phenotype )

outstruct <- NULL
if(TEST) {
  outstruct <- list(pe_tidied=tidied,pe_glanced=glanced, rsq=NULL,
                    modls=models)
} else {
  outstruct <- list(pe_tidied=tidied,pe_glanced=glanced, rsq=NULL, modls=models)
}

done <- dbDisconnect(con)

log_info("Done with ubiome-wide association for { phenotype }")

outfile_name <- NULL
if(is.null(opt$exposures)) {
  outfile_name <- sprintf('%s.rds', phenotype)
} else {
  exp_name <- file_path_sans_ext(basename(opt$exposures))
  outfile_name <- sprintf('%s_%s.rds', phenotype, exp_name)
}


saveRDS(outstruct, file=file.path(path_out, outfile_name))





