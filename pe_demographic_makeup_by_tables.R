# inputs a table of NHANES 'phenotype' table (e.g., BMX) and 'exposure' (e.g. DRXTOT) and estimates the demographic makeup
# the "demographic" makeup is the percent female/race/education, average poverty, average age,


library(getopt)
library(this.path)
setwd(this.dir())
source('quantpe.R')

TEST <- F
spec <- matrix(c(
  'phenotype_table', 'p', 1, "character",
  'exposure_table', 'e', 1, "character",
  'sample_size_pairs_list_file', 'l', 1, "character", # file that lists the sample sizes for each pair
  'path_to_db', 'i', 1, "character",
  'path_out', 'o', 1, "character"
), byrow=TRUE, ncol=4);
opt <- getopt(spec)

sample_size_threshold <- 500
adjustmentVariables <- c("RIDAGEYR", "AGE_SQUARED", 
                         "RIAGENDR", 
                         "INDFMPIR",
                         "EDUCATION_LESS9","EDUCATION_9_11","EDUCATION_AA","EDUCATION_COLLEGEGRAD",
                         "ETHNICITY_MEXICAN", "ETHNICITY_OTHERHISPANIC","ETHNICITY_OTHER", "ETHNICITY_NONHISPANICBLACK")
#phenotype_table <- 'BPX'
#phenotype_table <- 'BMX_E'
#phenotype_table <- 'BIOPRO_F'
#exposure_table <- 'L02HBS'
#exposure_table <- 'PHPYPA'
#exposure_table <- 'LAB06'
#exposure_table <- 'VOCWB_F'
#exposure_table <- 'PAQ_E'
#exposure_table <- 'DS1TOT_E'
#exposure_table <- 'PAQ'
#phenotype_table <- 'BMX_E'
#exposure_table <- 'DS2TOT_E'


#ss_file <- './select/sample_size_pe.csv'  #opt$sample_size_pairs_list_file
ss_file <- './select/sample_size_pe_category_060623.csv'
path_to_db <- './nhanes_122322.sqlite'
path_out <- './out'

if(!TEST) {
  phenotype_table <- opt$phenotype_table
  exposure_table <- opt$exposure_table
  ss_file <- opt$sample_size_pairs_list_file
  path_to_db <- opt$path_to_db
  path_out <- opt$path_out
}

con <- DBI::dbConnect(RSQLite::SQLite(), dbname=path_to_db)
to_do <- read_csv(ss_file) |> filter(e_table_name == exposure_table, p_table_name == phenotype_table, n >= sample_size_threshold)

m_table <- get_expo_pheno_tables(con, phenotype_table, exposure_table) |> figure_out_weight()

log_info("num pairs: {nrow(to_do)}")

if(nrow(to_do) == 0) {
  log_info("0 pairs, quitting")
  stop("No pairs to execute")
}


demo_breakdowns <- vector("list", length = nrow(to_do))
for(ii in 1:nrow(to_do)) {
  rw <- to_do |> slice(ii)
  log_info("{ii} out of {nrow(to_do)}; expo: {rw$evarname}; pheno: {rw$pvarname} ")
  pheno <-  rw$pvarname
  expo <- rw$evarname
  
  dat <- m_table$merged_tab |> filter(!is.na(wt), wt > 0, !is.na(pheno), !is.na(expo), if_all(all_of(adjustmentVariables), ~!is.na(.)))
  dsn <- create_svydesign(dat)  
  demo_breakdowns[[ii]] <- demographic_breakdown(dsn) |> mutate(varname_expo=expo, varname_pheno=pheno)
}
demo_breakdowns <- demo_breakdowns |> bind_rows() |> mutate(phenotable=phenotype_table, expotable=exposure_table)
outfile_name <- sprintf('%s_%s_demographic.rds', phenotype_table, exposure_table)
saveRDS(demo_breakdowns, file=file.path(path_out, outfile_name))

