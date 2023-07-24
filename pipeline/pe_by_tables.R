
# inputs a table of NHANES 'phenotype' table (e.g., BMX) and 'exposure' (e.g. DRXTOT) and associates each in one table with the other
# requires a file that pairs viable phenotype-exposures to associate
library(getopt)
#source('quantpe.R')
library(devtools)
source('db_paths.R')
load_all("..")

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

phenotype_table <- 'BPX'
#phenotype_table <- 'BMX_E'
#phenotype_table <- 'BIOPRO_F'
#exposure_table <- 'L02HBS'
#exposure_table <- 'PHPYPA'
#exposure_table <- 'LAB06'
#exposure_table <- 'VOCWB_F'
#exposure_table <- 'PAQ_E'
#exposure_table <- 'DS1TOT_E'
exposure_table <- 'PAQ'
#phenotype_table <- 'BMX_E'
#exposure_table <- 'DS2TOT_E'


#ss_file <- './select/sample_size_pe.csv'  #opt$sample_size_pairs_list_file
ss_file <- '../select/sample_size_pe_category_041823.csv'
path_to_db <- path_to_nhanes # '../nhanes_122322.sqlite'
path_out <- '../out'

if(!TEST) {
  phenotype_table <- opt$phenotype_table
  exposure_table <- opt$exposure_table
  ss_file <- opt$sample_size_pairs_list_file
  path_to_db <- opt$path_to_db
  path_out <- opt$path_out
}

#adjustmentVariables <- c("RIDAGEYR", "RIAGENDR", "INDFMPIR")
adjustmentVariables <- c("RIDAGEYR", "AGE_SQUARED",
                          "RIAGENDR",
                          "INDFMPIR",
                          "EDUCATION_LESS9","EDUCATION_9_11","EDUCATION_AA","EDUCATION_COLLEGEGRAD",
                          #"BORN_INUSA",
                          "ETHNICITY_MEXICAN", "ETHNICITY_OTHERHISPANIC","ETHNICITY_OTHER", "ETHNICITY_NONHISPANICBLACK")
con <- DBI::dbConnect(RSQLite::SQLite(), dbname=path_to_db)


to_do <- read_csv(ss_file) |> filter(e_table_name == exposure_table, p_table_name == phenotype_table, n >= sample_size_threshold)

m_table <- get_expo_pheno_tables(con, phenotype_table, exposure_table) |> figure_out_weight()

expo_levels <- tbl(con, 'e_variable_levels') |> filter(Data.File.Name == exposure_table) |> collect()

pe_safely <- safely(pe_by_table)
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

check_e_data_type <- function(varname)

  {
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


log_info("Process ID: {Sys.getpid()}")
log_info("num pairs: {nrow(to_do)}")

if(nrow(to_do) == 0) {
  log_info("0 pairs, quitting")
  log_info("Done with PxE: { phenotype_table } x { exposure_table }")
  stop("No pairs to execute")
}

for(ii in 1:nrow(to_do)) {
  rw <- to_do |> slice(ii)
  log_info("{ii} out of {nrow(to_do)}; expo: {rw$evarname}; pheno: {rw$pvarname} ")

  ## how to transform the phenotype and the exposure - log the phenotype

  #e_levels <- check_if_e_categorical(rw$evarname) ### categorical, rank, or none
  e_levels <- check_e_data_type(rw$evarname)

  log_info("{ii} e_levels { e_levels$vartype } ")
  mod <- NULL
  if(e_levels$vartype == 'continuous') {
    log_info("{ii} doing { rw$evarname } ")
    mod <- pe_safely(m_table, rw$pvarname, rw$evarname, adjustmentVariables,
                       logxform_p=F, logxform_e=T, scale_e=T, scale_p=T,
                       quantile_expo=NULL, exposure_levels=NULL)

  } else if(e_levels$vartype == 'categorical') {
    log_info("{ii} categorizing { rw$evarname } ")
    mod <- pe_safely(m_table, rw$pvarname, rw$evarname, adjustmentVariables,
                       logxform_p=F, logxform_e=F, scale_e=F, scale_p=T,
                       quantile_expo=NULL, exposure_levels=e_levels$varlevels)

  } else if(e_levels$vartype == 'continuous-rank') {
    log_info("{ii} as is { rw$evarname } ")
    mod <- pe_safely(m_table, rw$pvarname, rw$evarname, adjustmentVariables,
                     logxform_p=F, logxform_e=F, scale_e=T, scale_p=T,
                     quantile_expo=NULL, exposure_levels=NULL)

  }

  if(is.null(mod$result)) {
    tidied[[ii]] <- NULL
    glanced[[ii]] <- NULL
    log_error('{ii}: { mod$error }')
    next;
  }

  mod <- mod$result

  tidy_save <-  rbind(
    mod$adjusted$tidied |> mutate(evarname = rw$evarname, pvarname = rw$pvarname, model_type='adjusted', vartype=e_levels$vartype),
    mod$unadjusted$tidied |> mutate(evarname = rw$evarname, pvarname = rw$pvarname, model_type='unadjusted',vartype=e_levels$vartype),
    mod$base$tidied |> mutate(evarname = rw$evarname, pvarname = rw$pvarname, model_type='base',vartype=e_levels$vartype)
  )

  glance_save <-  rbind(
    mod$adjusted$glanced |> mutate(evarname = rw$evarname, pvarname = rw$pvarname, model_type='adjusted') |> cbind(mod$adjusted_model$r2 |> as_tibble()),
    mod$unadjusted$glanced |> mutate(evarname = rw$evarname, pvarname = rw$pvarname, model_type='unadjusted') |> cbind(mod$unadjusted_model$r2 |> as_tibble()),
    mod$base$glanced |> mutate(evarname = rw$evarname, pvarname = rw$pvarname, model_type='base') |> cbind(mod$base$r2 |> as_tibble())
  )
  tidied[[ii]] <- tidy_save
  glanced[[ii]] <- glance_save
}

outstruct <- list(pe_tidied=tidied |> bind_rows() |> mutate(exposure_table_name = exposure_table, phenotype_table_name=phenotype_table),
                  pe_glanced=glanced |> bind_rows() |> mutate(exposure_table_name = exposure_table, phenotype_table_name=phenotype_table))


log_info("Done with PxE: { phenotype_table } x { exposure_table }")
outfile_name <- sprintf('%s_%s.rds', phenotype_table, exposure_table)

saveRDS(outstruct, file=file.path(path_out, outfile_name))





