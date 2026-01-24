#!/usr/bin/env Rscript
# Chirag
# Exposure ~ demographics “ExWAS” (E-DWAS)
#
# Runs survey-weighted regressions with the exposure as the outcome:
#   exposure ~ demographics
#
# - Continuous exposures: log10-transform (with zero-handling) + scale (SD units) by default
# - Continuous-rank exposures: scale only
# - Categorical exposures: runs one-vs-reference binaries for each non-reference level (quasibinomial)
#
# Output mirrors your exwas.R style: tidied, glanced, rsq, plus raw models list.

suppressPackageStartupMessages({
  library(getopt)
  library(tidyverse)
  library(logger)
  library(DBI)
  library(RSQLite)
  library(purrr)
  library(survey)
  library(broom)
  library(nhanespewas)
})

devtools::load_all(".")
TEST <- FALSE

spec <- matrix(c(
  "path_to_db",        "i", 1, "character",
  "path_out",          "o", 1, "character",
  "exposures",         "e", 2, "character",  # optional CSV (1 column, no header) of exposures to run
  "include_survey_year","y", 2, "numeric"   # 0/1; default 1
), byrow = TRUE, ncol = 4)

opt <- getopt(spec)

# -------------------------
# Defaults / debug
# -------------------------
path_to_db <- "../db/nhanes_031725.sqlite"
path_out <- "."
include_survey_year <- 1
to_do <- tibble(evarname = c("LBXGTC", "LBXHPE", "LBXBEC"))

if (!TEST) {
  path_to_db <- opt$path_to_db
  path_out <- opt$path_out
  include_survey_year <- ifelse(is.null(opt$include_survey_year), 1, opt$include_survey_year)
  to_do <- read_csv(opt$exposures)
}


# -------------------------
# Helpers
# -------------------------

demo_table_name_for_begin_year <- function(yr) {
  if (yr == 1999) return("DEMO")
  if (yr == 2001) return("DEMO_B")
  if (yr == 2003) return("DEMO_C")
  if (yr == 2005) return("DEMO_D")
  if (yr == 2007) return("DEMO_E")
  if (yr == 2009) return("DEMO_F")
  if (yr == 2011) return("DEMO_G")
  if (yr == 2013) return("DEMO_H")
  if (yr == 2015) return("DEMO_I")
  if (yr == 2017) return("DEMO_J")
  stop("Unsupported Begin.Year: ", yr)
}

d_safely <- safely(expo_on_demographics)
con <- DBI::dbConnect(RSQLite::SQLite(), dbname=path_to_db)

N <- nrow(to_do)
models <- vector("list", length=N)

for(ii in 1:N) {
  rw <- to_do |> slice(ii)
  evar <- rw$evarname
  log_info("{ii} out of {nrow(to_do)}; expo: {evar} ")
  e_levels <- nhanespewas::check_e_data_type(evar, con)
  mod <- NULL

  if(e_levels$vartype == 'continuous') {
    mod <- d_safely(evar, con, logxform_exposure = T, scale_exposure = T, include_survey_year=T)
  }
  else if(e_levels$vartype == 'categorical') {
    mod <- d_safely(evar, con, logxform_exposure = F, scale_exposure = T, include_survey_year = T)
  }
  else if(e_levels$vartype == 'continuous-rank') {
    mod <- d_safely(evar, con, logxform_exposure = F, scale_exposure = T, include_survey_year = T)
  }
  models[[ii]] <- mod
}


# helper: accept either plain objects or purrr::safely outputs
.unwrap_safely <- function(x) {
  if (is.list(x) && ("error" %in% names(x)) && ("result" %in% names(x))) {
    if (!is.null(x$error) || is.null(x$result)) return(NULL)
    return(x$result)
  }
  x
}

tidied_result <- function(models) {
  map_dfr(models, function(x) {
    lcl <- .unwrap_safely(x)
    if (is.null(lcl) || is.null(lcl$tidied)) return(tibble())

    out <- lcl$tidied |>
      mutate(
        series   = if (!is.null(lcl$series)) paste(lcl$series, collapse = ";") else NA_character_,
        exposure = lcl$exposure %||% NA_character_,
        log_e    = lcl$logxform_exposure    %||% NA,
        scaled_e = lcl$scale_exposure %||% NA
      )

    out
  })
}

glanced_result <- function(models) {
  map_dfr(models, function(x) {
    lcl <- .unwrap_safely(x)
    if (is.null(lcl) || is.null(lcl$glanced)) return(tibble())

    lcl$glanced |>
      mutate(
        series   = if (!is.null(lcl$series)) paste(lcl$series, collapse = ";") else NA_character_,
        exposure = lcl$exposure %||% NA_character_,
        log_e    = lcl$logxform_exposure    %||% NA,
        scaled_e = lcl$scale_exposure %||% NA
      )
  })
}

r2_result <- function(models) {
  map_dfr(models, function(x) {
    lcl <- .unwrap_safely(x)
    if (is.null(lcl) || is.null(lcl$r2)) return(tibble())

    as_tibble(lcl$r2) |>
      mutate(
        series   = if (!is.null(lcl$series)) paste(lcl$series, collapse = ";") else NA_character_,
        exposure = lcl$exposure %||% NA_character_,
        log_e    = lcl$logxform_exposure    %||% NA,
        scaled_e = lcl$scale_exposure %||% NA
      )
  })
}

# small infix helper (optional)
`%||%` <- function(a, b) if (!is.null(a)) a else b


tidied  <- tidied_result(models)
glanced <- glanced_result(models)
rsq     <- r2_result(models)

outstruct <- NULL
if(TEST) {
  outstruct <- list(pe_tidied=tidied,pe_glanced=glanced, rsq=rsq,modls=models)
} else {
  outstruct <- list(pe_tidied=tidied,pe_glanced=glanced, rsq=rsq,modls=NULL)
}


done <- DBI::dbDisconnect(con)
log_info("Done with demographic associations")
base <- basename(opt$exposures)
name <- tools::file_path_sans_ext(base)
outfile_name <- sprintf("dwas_%s.rds", name)
saveRDS(outstruct, file=file.path(path_out, outfile_name))
