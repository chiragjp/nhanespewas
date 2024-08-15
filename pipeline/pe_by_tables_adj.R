
# inputs a table of NHANES 'phenotype' table (e.g., BMX) and 'exposure' (e.g. DRXTOT) and associates each in one table with the other
# requires a file that pairs viable phenotype-exposures to associate
library(getopt)
library(tidyverse)
library(logger)
#source('db_paths.R')
devtools::load_all("..")

TEST <- T
spec <- matrix(c(
  'phenotype_table', 'p', 1, "character",
  'exposure_table', 'e', 1, "character",
  'sample_size_pairs_list_file', 'l', 1, "character", # file that lists the sample sizes for each pair
  'path_to_db', 'i', 1, "character",
  'path_out', 'o', 1, "character"
), byrow=TRUE, ncol=4);
opt <- getopt(spec)

sample_size_threshold <- 500


########### debug stuff
#phenotype_table <- 'BPX'
#phenotype_table <- 'BMX_E'
#phenotype_table <- 'BIOPRO_F'
#exposure_table <- 'DRXTOT'
#exposure_table <- 'L02HBS'
#exposure_table <- 'PHPYPA'
#exposure_table <- 'LAB06'
#exposure_table <- 'VOCWB_F'
#exposure_table <- 'PAQ_E'
#exposure_table <- 'DR1TOT_E'
#exposure_table <- 'PAQ'
#phenotype_table <- 'BMX_E'
#exposure_table <- 'DS2TOT_E'
#phenotype_table <- 'BMX_D'
#exposure_table <- 'VITAEC_D'

phenotype_table <- "ALB_CR_F" # "ALB_CR_E"
exposure_table <- "HPVSER_F" #HPVSER_E" #


#ss_file <- './select/sample_size_pe.csv'  #opt$sample_size_pairs_list_file
ss_file <- '../select/sample_size_pe_category_060623.csv'
path_to_db <-   '../db/nhanes_012324.sqlite' # '../nhanes_122322.sqlite'
path_out <- '.'

if(!TEST) {
  phenotype_table <- opt$phenotype_table
  exposure_table <- opt$exposure_table
  ss_file <- opt$sample_size_pairs_list_file
  path_to_db <- opt$path_to_db
  path_out <- opt$path_out
}

############### end DEBUG



con <- DBI::dbConnect(RSQLite::SQLite(), dbname=path_to_db)
to_do <- read_csv(ss_file) |> filter(e_table_name == exposure_table, p_table_name == phenotype_table, n >= sample_size_threshold)
m_table <- get_expo_pheno_tables(con, phenotype_table, exposure_table) |> figure_out_weight()
expo_levels <- tbl(con, 'e_variable_levels') |> filter(Data.File.Name == exposure_table) |> collect()
pe_safely <- safely(nhanespewas::pe_by_table_flex_adjust)
tidied <- vector("list", length = nrow(to_do))
glanced <- vector("list", length = nrow(to_do))

check_if_e_categorical <- function(varname) {
  ## old function
  elvl <- expo_levels |> filter(Variable.Name == varname)
  if(nrow(elvl) > 1) {
    return(elvl |> pull(values) |> na.omit() |> sort())
  }
  return(NULL)
}

check_e_data_type <- function(varname) {
  ret <- list(vartype="continuous", varlevels=NULL)
  elvl <- expo_levels |> filter(Variable.Name == varname, !is.na(values)) |> pull(values)

  if(grepl('CNT$', varname)) {
    return(list(vartype="continuous-rank", varlevels=elvl))
  }

  if(grepl("^PAQ", varname)) {
    return(list(vartype="continuous", varlevels=NULL))
  }

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
log_info("num pairs: {nrow(to_do)}")

if(nrow(to_do) == 0) {
  dbDisconnect(con)
  log_info("0 pairs, quitting")
  log_info("Done with PxE: { phenotype_table } x { exposure_table }")
  stop("No pairs to execute")
}


N <- nrow(to_do)
models <- vector("list", length=N)
for(ii in 1:N) {
  rw <- to_do |> slice(ii)
  log_info("{ii} out of {nrow(to_do)}; expo: {rw$evarname}; pheno: {rw$pvarname} ")

  ## how to transform the phenotype and the exposure - log the phenotype

  e_levels <- check_e_data_type(rw$evarname)
  adjustment_model_for_e <- adjustment_scenario_for_variable(rw$evarname, rw$pvarname)
  #print(adjustment_model_for_e)
  log_info("{ii} e_levels { e_levels$vartype } ")
  mod <- NULL
  if(e_levels$vartype == 'continuous') {
    log_info("{ii} continuous { rw$evarname } ")
    mod <- pe_safely(m_table, rw$pvarname, rw$evarname, adjustment_model_for_e,
                     logxform_p=F, logxform_e=T, scale_e=T, scale_p=T,
                     quantile_expo=NULL, exposure_levels=NULL)

  } else if(e_levels$vartype == 'categorical') {
    log_info("{ii} categorizing { rw$evarname } ")
    mod <- pe_safely(m_table, rw$pvarname, rw$evarname, adjustment_model_for_e,
                     logxform_p=F, logxform_e=F, scale_e=F, scale_p=T,
                     quantile_expo=NULL, exposure_levels=e_levels$varlevels)

  } else if(e_levels$vartype == 'continuous-rank') {
    log_info("{ii} as is { rw$evarname } ")
    mod <- pe_safely(m_table, rw$pvarname, rw$evarname, adjustment_model_for_e,
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
                                  series = lcl$series,
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
                                  series = lcl$series,
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
                                  series = lcl$series,
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
) |> mutate(exposure_table_name = exposure_table, phenotype_table_name = phenotype_table )
glanced <- rbind(
  glanced_result(models, F) |> mutate(aggregate_base_model = F),
  glanced_result(models, T) |> mutate(aggregate_base_model = T)
) |> mutate(exposure_table_name = exposure_table, phenotype_table_name = phenotype_table )

rsq <- rbind(
  r2_result(models, F) |> mutate(aggregate_base_model = F),
  r2_result(models, T) |> mutate(aggregate_base_model = T)
) |> mutate(exposure_table_name = exposure_table, phenotype_table_name = phenotype_table )

outstruct <- NULL
if(TEST) {
  outstruct <- list(pe_tidied=tidied,pe_glanced=glanced, rsq=rsq,
                    modls=models)
} else {
  outstruct <- list(pe_tidied=tidied,pe_glanced=glanced, rsq=rsq, modls=models)
}


log_info("Done with PxE: { phenotype_table } x { exposure_table }")
outfile_name <- sprintf('%s_%s.rds', phenotype_table, exposure_table)

saveRDS(outstruct, file=file.path(path_out, outfile_name))





