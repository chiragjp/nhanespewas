# Chirag
# August 2024
# exwas.R: perform ExWAS for a phenotype

library(getopt)
library(tidyverse)
library(logger)
library(tools)
devtools::load_all("..")
#Error in check_e_data_type(rw$evarname) : object 'elvl' not found
TEST <- F
spec <- matrix(c(
  'phenotype', 'p', 1, "character",
  'sample_size_pairs_list_file', 'l', 1, "character", # file that lists the sample sizes for each pair
  'path_to_db', 'i', 1, "character",
  'path_out', 'o', 1, "character",
  'exposures', 'e', 2, "character"
), byrow=TRUE, ncol=4);
opt <- getopt(spec)

sample_size_threshold <- 500


########### debug stuff
phenotype <- "BMXBMI"
#ss_file <- './select/sample_size_pe.csv'  #opt$sample_size_pairs_list_file
ss_file <- '../select/sample_size_pe_category_060623.csv'
path_to_db <-   '../db/nhanes_012324.sqlite' # '../nhanes_122322.sqlite'
path_out <- '.'

if(!TEST) {
  phenotype <- opt$phenotype
  ss_file <- opt$sample_size_pairs_list_file
  path_to_db <- opt$path_to_db
  path_out <- opt$path_out
}

############### end DEBUG



con <- DBI::dbConnect(RSQLite::SQLite(), dbname=path_to_db)

to_do <- read_csv(ss_file) |> filter(pvarname == phenotype) |> group_by(evarname) |> summarize(total_n=sum(n), num_surveys = n()) |> mutate(pvarname = phenotype)
to_do <- to_do |> filter(num_surveys >=2, total_n >= sample_size_threshold)

if(!is.null(opt$exposures)) {
  exposures <- read_csv(opt$exposures, col_names=F)
  to_do <- to_do |> inner_join(exposures |> rename(evarname=X1), by="evarname")
}

pe_safely <- safely(nhanespewas::pe_flex_adjust)
tidied <- vector("list", length = nrow(to_do))
glanced <- vector("list", length = nrow(to_do))


check_e_data_type <- function(varname, con) {
  ret <- list(vartype="continuous", varlevels=NULL)

  if(grepl('CNT$', varname)) {
    return(list(vartype="continuous-rank", varlevels=NULL))
  }

  if(grepl("^PAQ", varname)) {
    return(list(vartype="continuous", varlevels=NULL))
  }

  elvl <- tbl(con, 'e_variable_levels') |>  filter(Variable.Name == varname, !is.na(values)) |> pull(values) |> unique()

  if(length(elvl) == 1) {
    return(list(vartype="continuous", varlevels=NULL))
  } else if(any(elvl < 1 & elvl > 0) | any(round(elvl) != elvl)) {
    return(list(vartype="continuous-rank", varlevels=sort(elvl)))
  } else if(all(round(elvl) == elvl)) {
    return(list(vartype="categorical", varlevels=sort(elvl)))
  }
  return(ret)
}

adjustment_scenario_for_variable <- function(evarname, pvarname) {
  first_three <- substr(evarname, 1, 3)
  first_two <- substr(evarname, 1, 2)
  first_two_p <- substr(pvarname, 1, 2)
  if(first_three == 'DRX') {
    log_info("Adjusting by total calories")
    return(adjustment_models_diet_x)
  } else if(first_three == 'DR1') {
    log_info("Adjusting by total calories")
    return(adjustment_models_diet_1)
  } else if(first_three == 'DR2') {
    log_info("Adjusting by total calories")
    return(adjustment_models_diet_2)
  } else if(first_two == 'UR' & first_two_p != 'UR') {
    log_info("Adjusting by creatinine")
    return(adjustment_models_ucr)
  }else {
    log_info("Default adjustments")
    return(adjustment_models)
  }

}


log_info("Process ID: {Sys.getpid()}")
log_info("Number of Exposures: {nrow(to_do)}")

if(nrow(to_do) == 0) {
  done <- dbDisconnect(con)
  log_info("0 pairs, quitting")
  log_info("Done for the phenotype { phenotype } ")
  stop("No pairs to execute")
}


N <- nrow(to_do)
models <- vector("list", length=N)
if(TEST) { N = 3; }
for(ii in 1:N) {
  rw <- to_do |> slice(ii)
  log_info("{ii} out of {nrow(to_do)}; expo: {rw$evarname}; pheno: {rw$pvarname} ")

  ## how to transform the phenotype and the exposure - log the phenotype

  e_levels <- check_e_data_type(rw$evarname, con)
  adjustment_model_for_e <- adjustment_scenario_for_variable(rw$evarname, rw$pvarname)
  #print(adjustment_model_for_e)
  log_info("{ii} e_levels { e_levels$vartype } ")
  mod <- NULL
  if(e_levels$vartype == 'continuous') {
    log_info("{ii} continuous { rw$evarname } ")
    # pe
    #pe <- function(pheno, exposure, adjustment_variables,con, series=NULL,
    #               logxform_p=T, logxform_e=T, scale_e=T, scale_p=F,
    #               pheno_table_name=NULL, expo_table_name=NULL,
    #               quantile_expo=NULL, exposure_levels=NULL) {
    mod <- pe_safely(rw$pvarname, rw$evarname, adjustment_model_for_e, con,
                     logxform_p=T, logxform_e=T, scale_e=T, scale_p=T,
                     quantile_expo=NULL, exposure_levels=NULL)

  } else if(e_levels$vartype == 'categorical') {
    log_info("{ii} categorizing { rw$evarname } ")
    mod <- pe_safely(rw$pvarname, rw$evarname, adjustment_model_for_e, con,
                     logxform_p=F, logxform_e=F, scale_e=F, scale_p=T,
                     quantile_expo=NULL, exposure_levels=e_levels$varlevels)

  } else if(e_levels$vartype == 'continuous-rank') {
    log_info("{ii} as is { rw$evarname } ")
    mod <- pe_safely(rw$pvarname, rw$evarname, adjustment_model_for_e, con,
                     logxform_p=F, logxform_e=F, scale_e=T, scale_p=T,
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

outstruct <- NULL
if(TEST) {
  outstruct <- list(pe_tidied=tidied,pe_glanced=glanced, rsq=rsq,
                    modls=models)
} else {
  outstruct <- list(pe_tidied=tidied,pe_glanced=glanced, rsq=rsq, modls=models)
}

done <- dbDisconnect(con)

log_info("Done with ExWAS for { phenotype }")

outfile_name <- NULL
if(is.null(opt$exposures)) {
  outfile_name <- sprintf('%s.rds', phenotype)
} else {
  exp_name <- file_path_sans_ext(basename(opt$exposures))
  outfile_name <- sprintf('%s_%s.rds', phenotype, exp_name)
}


saveRDS(outstruct, file=file.path(path_out, outfile_name))





