# Chirag
# August 2024
# exwas.R: perform ExWAS for a phenotype


#Description:
  #   This script performs a Phenotype-Exposure-Wide Association Study (ExWAS) for a
  #   specified phenotype using NHANES data stored in a SQLite database. It reads a list
  #   of (phenotype, exposure) pairs that meet minimum sample size criteria, runs survey-
  #   weighted regression models across different adjustment scenarios, and outputs
  #   tidy, glance, and R² results for each exposure–phenotype combination. Results
  #   are saved as an RDS file for downstream meta‐analysis or visualization.
  #
  # Usage:
  #   Rscript exwas.R
  #     --phenotype           <string>   Name of the phenotype variable to analyze
  #     --sample_size_pairs_list_file <string>  CSV file listing sample sizes for each (phenotype, exposure) pair
  #     --path_to_db          <string>   File path to the NHANES SQLite database
  #     --path_out            <string>   Directory where output RDS will be saved
  #     [--exposures          <string>   Optional CSV file (one exposure per line) to restrict the exposures processed]
  #     [--use_quantile       <0 or 1>    Optional flag: if 1, use quantile-based exposure modeling; default is 0]
  #
  # Arguments:
  #   --phenotype
  #       The NHANES variable name of the phenotype to analyze (e.g., "LBXRDW").
  #
  #   --sample_size_pairs_list_file
  #       Path to a CSV that contains columns: pvarname, evarname, n, etc. This file is
  #       used to filter all (phenotype,exposure) pairs by sample size (minimum 500)
  #       and number of surveys (≥2). Each row corresponds to one cohort’s n for a pair. Optional.
  #
  #   --path_to_db
  #       Path to an existing SQLite database containing cleaned NHANES tables. The script
  #       uses DBI and RSQLite to connect and pull data for each exposure and phenotype.
  #
  #   --path_out
  #       Directory where the script will save the output RDS. The filename will be either
  #       "<phenotype>.rds" or "<phenotype>_<exposures_basename>.rds” if --exposures is given.
  #
  #   --exposures (optional)
  #       Path to a one‐column CSV file (no header) listing specific exposure variable names
  #       to include. If provided, the script will only run ExWAS on those exposures that
  #       appear both in the sample‐size file and in this list.
  #
  #   --use_quantile (optional)
  #       Numeric flag (0 or 1). If 1, continuous exposures are modeled using quantile bins
  #       (quartiles and deciles defined by QUANTILES = c(0, .25, .5, .75, .9, 1)). Default is 0.
  #
  # Details:
  #   1. The script begins by parsing command‐line options via getopt().
  #   2. It sets a minimum sample size threshold (default 500) and filters out any (phenotype,
  #      exposure) pairs with fewer than 500 total participants or fewer than 2 survey waves. If this file is not supplied, the script will use the package provided file
  #   3. If --exposures is supplied, only those exposures will be retained for analysis.
  #   4. For each remaining (phenotype, exposure) pair, the script:
  #       a. Determines exposure data type via nhanespewas::check_e_data_type():
  #          - "continuous": log‐transform, scale, or quantile‐bin according to --use_quantile
  #          - "categorical": treat as factor with levels returned by check_e_data_type()
  #          - "continuous-rank": scale as numeric (no log transform)
  #       b. Builds an adjustment scenario (covariates to include) using
  #          nhanespewas::adjustment_scenario_for_variable() + survey year.
  #       c. Calls nhanespewas::pe_flex_adjust() (wrapped in purrr::safely()) with:
  #          - phenotype name
  #          - exposure name
  #          - adjustment formula list
  #          - database connection
  #          - options: log_transform phenotype or exposure, scale flags, quantile bins, scale_type
  #       d. Stores the returned model object (or error) in a list “models”.
  #   5. After looping over all exposures, the script assembles three tables:
  #       • pe_tidied: one row per model (tidy regression coefficients and SEs)
  #       • pe_glanced: one row per model (glance summary: p‐value, df, AIC, etc.)
  #       • rsq: one row per model (variance explained by each exposure, R²)
  #      Each record includes metadata columns: model_number, series (survey waves),
  #      exposure, phenotype, log_p, log_e, scaled_p, scaled_e, and a flag “aggregate_base_model”
  #      indicating whether the base (minimally adjusted) model was used.
  #   6. Finally, the script disconnects from the SQLite DB and saves a single RDS object
  #      containing: list(pe_tidied, pe_glanced, rsq, modls=models). This RDS file is named
  #      "<phenotype>.rds" or "<phenotype>_<exposures_basename>.rds" under --path_out.
  #
  # Dependencies:
  #   - R (>=3.6)
  #   - Packages: getopt, tidyverse, logger, tools, nhanespewas, DBI, RSQLite, purrr
  #
  # Examples:
  #   # Run ExWAS for phenotype LBXRDW using sample-size file and full NHANES DB:
  #   Rscript exwas.R \
  #     --phenotype LBXRDW \
  #     --sample_size_pairs_list_file ../select/sample_size_pe_category_0824.csv \
  #     --path_to_db ../db/nhanes_031725.sqlite \
  #     --path_out ./results
  #
  #   # Restrict to a subset of exposures listed in exposures_to_run.csv, apply quantile bins:
  #   Rscript exwas.R \
  #     --phenotype LBXRDW \
  #     --sample_size_pairs_list_file ../select/sample_size_pe_category_0824.csv \
  #     --path_to_db ../db/nhanes_031725.sqlite \
  #     --path_out ./results \
  #     --exposures exposures_to_run.csv \
  #     --use_quantile 1
  #
  # Notes:
  #   • The file “sample_size_pe_category_0824.csv” must include columns: pvarname, evarname, n, etc.
  #   • Users should ensure the SQLite database includes the necessary NHANES tables and indexes used by nhanespewas.
  #   • If no eligible (phenotype, exposure) pairs remain (after filtering), the script logs “0 pairs, quitting” and exits.
  #   • To debug or test on a single exposure, set TEST <- TRUE near the top; this overrides getopt().
  #
  # License:
  #   MIT License (see LICENSE file in the nhanespewas repository)
  #
  # --------------------------------------------------------------



library(getopt)
library(tidyverse)
library(logger)
library(nhanespewas)
library(tools)
#devtools::load_all(".")

TEST <- F
spec <- matrix(c(
  'phenotype', 'p', 1, "character",
  'sample_size_pairs_list_file', 'l', 2, "character", # file that lists the sample sizes for each pair
  'path_to_db', 'i', 1, "character",
  'path_out', 'o', 1, "character",
  'exposures', 'e', 2, "character",
  'adjustment_scenario', 's', 2, "character",
  'interact_with','w', 2, "character",
  'use_quantile', 'q', 2, "numeric"
), byrow=TRUE, ncol=4);
opt <- getopt(spec)

sample_size_threshold <- 500

########### debug stuff
exposure <- 'LBXBCD'
phenotype <- 'LBXRDW'
#ss_file <- '../select/sample_size_pe_category_0824.csv'
ss_file <- system.file("extdata", "sample_size_pe_category_0824.csv", package = "nhanespewas")
path_to_db <-'../db/nhanes_031725.sqlite'
path_out <- '.'
use_quantile <- 0
#adjustment_scenario <- "age_sex_ethnicity_income_education"
adjustment_scenario <- NULL
interact_with <- NULL
#interact_with <- "RIDAGEYR"

if(!TEST) {
  phenotype <- opt$phenotype
  if(is.null(opt$sample_size_pairs_list_file)) {
    ss_file <- system.file("extdata", "sample_size_pe_category_0824.csv", package = "nhanespewas")
  } else {
    ss_file <- opt$sample_size_pairs_list_file
  }
  path_to_db <- opt$path_to_db
  path_out <- opt$path_out
  use_quantile <- ifelse(is.null(opt$use_quantile), 0, opt$use_quantile)
  adjustment_scenario <- opt$adjustment_scenario
  interact_with <-opt$interact_with
}

QUANTILES <- c(0, .25, .5, .75, .9, 1)
SCALE_P <- T
SCALE_TYPE <- 1 # Scaled by SD
############### end DEBUG

con <- DBI::dbConnect(RSQLite::SQLite(), dbname=path_to_db)
to_do <- read_csv(ss_file) |> filter(pvarname == phenotype) |> group_by(evarname) |> summarize(total_n=sum(n), num_surveys = n()) |> mutate(pvarname = phenotype)
to_do <- to_do |> filter(num_surveys >=2, total_n >= sample_size_threshold)

if(TEST) {
  to_do <- to_do |> filter(evarname == exposure)
}

if(!is.null(opt$exposures)) {
  exposures <- read_csv(opt$exposures, col_names=F)
  to_do <- to_do |> inner_join(exposures |> rename(evarname=X1), by="evarname")
}

pe_safely <- safely(nhanespewas::pe_flex_adjust)
tidied <- vector("list", length = nrow(to_do))
glanced <- vector("list", length = nrow(to_do))

log_info("Process ID: {Sys.getpid()}")
log_info("Number of Exposures: {nrow(to_do)}")

if(nrow(to_do) == 0) {
  done <- DBI::dbDisconnect(con)
  log_info("0 pairs, quitting")
  log_info("Done for the phenotype { phenotype } ")
  stop("No pairs to execute")
}


N <- nrow(to_do)
models <- vector("list", length=N)

for(ii in 1:N) {
  rw <- to_do |> slice(ii)
  log_info("{ii} out of {nrow(to_do)}; expo: {rw$evarname}; pheno: {rw$pvarname} ")

  ## how to transform the phenotype and the exposure - log the phenotype

  e_levels <- nhanespewas::check_e_data_type(rw$evarname, con)
  adjustment_model_for_e <- nhanespewas::adjustment_scenario_for_variable(rw$evarname, rw$pvarname)

  ## add in the SDDRVYR in the non-base models
  adjustment_model_for_e <- nhanespewas:::add_survey_year_to_adjustment_scenarios(adjustment_model_for_e)

  if(!is.null(adjustment_scenario)) {
    adjustment_model_for_e <- adjustment_model_for_e |> filter(scenario == adjustment_scenario)
  }

  log_info("{ii} e_levels { e_levels$vartype } ")
  mod <- NULL
  if(e_levels$vartype == 'continuous') {
    log_info("{ii} continuous { rw$evarname } ")
    # pe
    if(use_quantile==1) {
      mod <- pe_safely(rw$pvarname, rw$evarname, adjustment_model_for_e, con,
                       logxform_p=F, logxform_e=F, scale_e=F, scale_p=SCALE_P,
                       quantile_expo=QUANTILES, exposure_levels=NULL,scale_type=SCALE_TYPE, interact_with=interact_with)
    } else {
      mod <- pe_safely(rw$pvarname, rw$evarname, adjustment_model_for_e, con,
                     logxform_p=F, logxform_e=T, scale_e=T, scale_p=SCALE_P,
                     quantile_expo=NULL, exposure_levels=NULL, scale_type=SCALE_TYPE, interact_with=interact_with)
    }

  } else if(e_levels$vartype == 'categorical') {
    log_info("{ii} categorizing { rw$evarname } ")
    mod <- pe_safely(rw$pvarname, rw$evarname, adjustment_model_for_e, con,
                     logxform_p=F, logxform_e=F, scale_e=F, scale_p=SCALE_P,
                     quantile_expo=NULL, exposure_levels=e_levels$varlevels, scale_type=SCALE_TYPE, interact_with=interact_with)

  } else if(e_levels$vartype == 'continuous-rank') {
    log_info("{ii} as is { rw$evarname } ")
    mod <- pe_safely(rw$pvarname, rw$evarname, adjustment_model_for_e, con,
                     logxform_p=F, logxform_e=F, scale_e=T, scale_p=SCALE_P,
                     quantile_expo=NULL, exposure_levels=NULL, scale_type=SCALE_TYPE, interact_with=interact_with)

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
                                  log_p= lcl$log_p,
                                  log_e = lcl$log_e,
                                  scaled_p = lcl$scaled_p,
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
                                  log_p= lcl$log_p,
                                  log_e = lcl$log_e,
                                  scaled_p = lcl$scaled_p,
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
                                  log_p= lcl$log_p,
                                  log_e = lcl$log_e,
                                  scaled_p = lcl$scaled_p,
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

regterm_result <- function(models) {
  map_dfr(models, ~ {
    # skip errored runs
    if (!is.null(.x$error)) return(tibble())

    lcl <- .x$result
    # pick the right set: main‐models vs. base‐models
    to_do <- lcl$models

    map2_dfr(to_do, seq_along(to_do), ~ {
      mod_item <- .x
      idx      <- .y

      if (is.null(mod_item) || is.null(mod_item$reg_term_tidied)) {
        return(tibble())
      }

      mod_item$reg_term_tidied %>%
        mutate(
          model_number = idx,
          series       = paste(lcl$series, collapse = ";"),
          exposure     = lcl$exposure,
          phenotype    = lcl$phenotype,
          log_p        = lcl$log_p,
          log_e        = lcl$log_e,
          scaled_p     = lcl$scaled_p,
          scaled_e     = lcl$scaled_e
        )
    })
  })
}

tidied <- rbind(
  tidied_result(models, F) |> mutate(aggregate_base_model = F),
  tidied_result(models, T) |> mutate(aggregate_base_model = T)
) |> mutate(phenotype = phenotype)
glanced <- rbind(
  glanced_result(models, F) |> mutate(aggregate_base_model = F),
  glanced_result(models, T) |> mutate(aggregate_base_model = T)
) |> mutate(phenotype = phenotype )

rsq <- rbind(
  r2_result(models, F) |> mutate(aggregate_base_model = F),
  r2_result(models, T) |> mutate(aggregate_base_model = T)
) |> mutate(phenotype=phenotype )

reg_term_test_tidied <- regterm_result(models) |> mutate(aggregate_base_model=F, phenotype=phenotype)

outstruct <- NULL
if(TEST) {
  outstruct <- list(pe_tidied=tidied,pe_glanced=glanced, rsq=rsq, reg_term_test_tidied=reg_term_test_tidied,modls=models)
} else {
  outstruct <- list(pe_tidied=tidied,pe_glanced=glanced, rsq=rsq,reg_term_test_tidied=reg_term_test_tidied,modls=models)
}

done <- DBI::dbDisconnect(con)

log_info("Done with ExWAS for { phenotype }")

outfile_name <- NULL
if(is.null(opt$exposures)) {
  outfile_name <- sprintf('%s.rds', phenotype)
} else {
  exp_name <- file_path_sans_ext(basename(opt$exposures))
  outfile_name <- sprintf('%s_%s.rds', phenotype, exp_name)
}


saveRDS(outstruct, file=file.path(path_out, outfile_name))





