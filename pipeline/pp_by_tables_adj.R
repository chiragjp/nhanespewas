
# inputs a table of NHANES 'pheno' table (e.g., DRXTOT) and 'pheno' (e.g. DRXTOT) and associates each in one table with the other
# requires a file that pairs viable phenotype-exposures to associate
library(getopt)
library(tidyverse)
library(logger)
devtools::load_all("..")

TEST <- F
spec <- matrix(c(
  'y_table', 'y', 1, "character",
  'x_table', 'x', 1, "character",
  'sample_size_pairs_list_file', 'l', 1, "character", # file that lists the sample sizes for each pair
  'path_to_db', 'i', 1, "character",
  'path_out', 'o', 1, "character"
), byrow=TRUE, ncol=4);
opt <- getopt(spec)

sample_size_threshold <- 500


########### debug stuff
#x_table <- 'BMX_D'
#y_table <- 'BMX_D'

ss_file <- '../select/sample_size_pp_category_032024.csv'
path_to_db <-   '../db/nhanes_012324.sqlite' # '../nhanes_122322.sqlite'
path_out <- '.'

if(!TEST) {
  y_table <- opt$y_table
  x_table <- opt$x_table
  ss_file <- opt$sample_size_pairs_list_file
  path_to_db <- opt$path_to_db
  path_out <- opt$path_out
}

############### end DEBUG

con <- DBI::dbConnect(RSQLite::SQLite(), dbname=path_to_db)
to_do <- read_csv(ss_file) |> filter(x_table_name == x_table, y_table_name == y_table, n >= sample_size_threshold)
m_table <- get_x_y_tables(con, x_table, y_table) 
m_table$p_table <- m_table$table1
m_table$e_table <- m_table$table2 
m_table <- m_table |> figure_out_weight()
ppcor_safely <- safely(xy_by_table_flex_adjust)
tidied <- vector("list", length = nrow(to_do))
glanced <- vector("list", length = nrow(to_do))


pp_adjustment_models <- adjustment_models |> filter(
  scenario == 'base' |
  scenario == 'age_sex_ethnicity_income_education' 
)

log_info("Process ID: {Sys.getpid()}")
log_info("num pairs: {nrow(to_do)}")

if(nrow(to_do) == 0) {
  log_info("0 pairs, quitting")
  log_info("Done with PxP: { y_table } x { x_table }")
  stop("No pairs to execute")
}


N <- nrow(to_do)
models <- vector("list", length=N)
for(ii in 1:N) {
  rw <- to_do |> slice(ii)
  log_info("{ii} out of {nrow(to_do)}; y: {rw$yvarname}; x: {rw$xvarname} ")
  mod <- NULL
  mod <- ppcor_safely(m_table, rw$yvarname, rw$xvarname, pp_adjustment_models, logxform_x=T, logxform_y=T, scale_x=T, scale_y=T)
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
                                  xvar = lcl$xvar,
                                  yvar= lcl$yvar,
                                  log_y= lcl$log_y,
                                  log_x = lcl$log_x,
                                  scaled_y = lcl$scaled_y,
                                  scaled_x = lcl$scaled_x
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
                                  xvar = lcl$xvar,
                                  yvar= lcl$yvar,
                                  log_y= lcl$log_y,
                                  log_x = lcl$log_x,
                                  scaled_y = lcl$scaled_y,
                                  scaled_x = lcl$scaled_x
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
                                  xvar = lcl$xvar,
                                  yvar= lcl$yvar,
                                  log_y= lcl$log_y,
                                  log_x = lcl$log_x,
                                  scaled_y = lcl$scaled_y,
                                  scaled_x = lcl$scaled_x
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
) |> mutate(y_table_name = y_table, x_table_name = x_table )
glanced <- rbind(
  glanced_result(models, F) |> mutate(aggregate_base_model = F),
  glanced_result(models, T) |> mutate(aggregate_base_model = T)
) |> mutate(y_table_name = y_table, x_table_name = x_table )

rsq <- rbind(
  r2_result(models, F) |> mutate(aggregate_base_model = F),
  r2_result(models, T) |> mutate(aggregate_base_model = T)
) |> mutate(y_table_name = y_table, x_table_name = x_table )

outstruct <- NULL
if(TEST) {
  outstruct <- list(pe_tidied=tidied,pe_glanced=glanced, rsq=rsq,
                    modls=models)
} else {
  outstruct <- list(pe_tidied=tidied,pe_glanced=glanced, rsq=rsq, modls=models)
}


log_info("Done with PxP: { y_table } x { x_table }")
outfile_name <- sprintf('%s_%s.rds', y_table, x_table)

saveRDS(outstruct, file=file.path(path_out, outfile_name))





