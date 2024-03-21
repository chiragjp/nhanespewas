## Get the sample sizes per correlation and survey year

library(devtools)
load_all("..")

library(corrr)
library(tidyverse)
library(future)
library(logger)
source('db_paths.R')

con <- DBI::dbConnect(RSQLite::SQLite(), dbname=path_to_nhanes)
x_variables <- read_csv('../select/select_expo_variables_3.csv')
y_variables <- read_csv('../select/select_expo_variables_3.csv')

variable_information_selected <- rbind(x_variables, y_variables)

merged_table <- function(con, tab1, tab2) {
  if(tab1 == tab2) {
    return(tbl(con, tab1) |> collect())
  }
  tbl(con, tab1) |> inner_join(tbl(con,tab2), by = "SEQN") |> collect()
}

sample_size_for_table <- function(m_table, xvars, yvars) {
  sample_size_pair <- m_table |> pair_n() |> as_cordf()
  vars <- intersect(c(xvars, yvars), colnames(m_table))
  sample_size_pair <- sample_size_pair |> dice(all_of(vars)) |> stretch(na.rm=F)
  # Check if 'r' column exists, if not, create 'n' with NA values
  if (!"r" %in% colnames(sample_size_pair)) {
    #sample_size_pair$n <- NA
    return(NULL)
  }
  sample_size_pair <- sample_size_pair |> rename(n = r)
  sample_size_pair <- sample_size_pair |> filter(x %in% xvars, y %in% yvars) |> rename(xvarname=x, yvarname=y)
  sample_size_pair
}

sample_size_for_xy_tables <- function(con, xTable, yTable, variable_information, epcf_x = 'e', epcf_y='p') {
  log_info("Working on { xTable } x {yTable}")
  m_table <- merged_table(con, yTable, xTable)  |> collect()
  xvars <- variable_information |> filter(Data.File.Name == xTable, epcf == epcf_x) |> select(Variable.Name) |> pull()
  yvars <- variable_information |> filter(Data.File.Name == yTable, epcf == epcf_y) |> select(Variable.Name) |> pull()
  tab <-sample_size_for_table(m_table, xvars,yvars)
  if(is.null(tab)) {
    return(tab)
  }
  return(tab |> mutate(x_table_name = xTable, y_table_name = yTable))
}

xy_sample_sizes <- function(con, x_vars, y_vars, x_epcf='p', y_epcf='p') {
  ## per survey year
  x_tables <- x_vars |> select(Data.File.Name) |> pull() |> unique()
  y_tables <- y_vars |> select(Data.File.Name) |> pull() |> unique()
  var_information <- rbind(x_vars, y_vars)
  log_info("Num e_tables { length(x_tables) }")
  log_info("Num p_tables { length(y_tables) }")
  for_x <- vector("list", length(x_tables))
  for(x in 1:length(x_tables)) {
    x_tab_name <- x_tables[x]
    if((x %% 10) == 0) {
      log_info("Num x done { x }")
    }
    for_x[[x]] <- map(y_tables, ~sample_size_for_xy_tables(con, x_tab_name, .x, var_information, epcf_x = x_epcf, epcf_y = y_epcf))
  }
  for_x |> bind_rows()
}




log_info("Sample size for survey:")
log_info("A")
samp_size_a <-  xy_sample_sizes(con, x_variables |> filter(Begin.Year == 1999), y_variables |> filter(Begin.Year == 1999), x_epcf='e', y_epcf='e')
log_info("B")
samp_size_b <- xy_sample_sizes(con, x_variables |> filter(Begin.Year == 2001), y_variables |> filter(Begin.Year == 2001), x_epcf='e', y_epcf='e')
log_info("C")
samp_size_c <- xy_sample_sizes(con, x_variables |> filter(Begin.Year == 2003), y_variables |> filter(Begin.Year == 2003), x_epcf='e', y_epcf='e')
log_info("D")
samp_size_d <- xy_sample_sizes(con, x_variables |> filter(Begin.Year == 2005), y_variables |> filter(Begin.Year == 2005), x_epcf='e', y_epcf='e')
log_info("E")
samp_size_e <- xy_sample_sizes(con, x_variables |> filter(Begin.Year == 2007), y_variables |> filter(Begin.Year == 2007), x_epcf='e', y_epcf='e')
log_info("F")
samp_size_f <- xy_sample_sizes(con, x_variables |> filter(Begin.Year == 2009), y_variables |> filter(Begin.Year == 2009), x_epcf='e', y_epcf='e')
log_info("G")
samp_size_g <- xy_sample_sizes(con, x_variables |> filter(Begin.Year == 2011), y_variables |> filter(Begin.Year == 2011), x_epcf='e', y_epcf='e')
log_info("H")
samp_size_h <- xy_sample_sizes(con, x_variables |> filter(Begin.Year == 2013), y_variables |> filter(Begin.Year == 2013), x_epcf='e', y_epcf='e')
log_info("I")
samp_size_i <- xy_sample_sizes(con, x_variables |> filter(Begin.Year == 2015), y_variables |> filter(Begin.Year == 2015), x_epcf='e', y_epcf='e')
log_info("J")
samp_size_j <- xy_sample_sizes(con, x_variables |> filter(Begin.Year == 2017), y_variables |> filter(Begin.Year == 2017), x_epcf='e', y_epcf='e')

samp_sizes <- rbind(
  samp_size_a, samp_size_b, samp_size_c, samp_size_d, samp_size_e, samp_size_f, samp_size_g, samp_size_h, samp_size_i, samp_size_j
)

#samp_sizes |> write_csv('./select/sample_size_pe_category_041823.csv')
#samp_sizes |> write_csv('./select/sample_size_pe_category_051323.csv')
samp_sizes |> write_csv('../select/sample_size_ee_category_032024.csv')
#rbind(samp_size_a, samp_size_b) |> filter(pvarname == 'TeloMean' | pvarname == 'TeloStd') |> write_csv('./select/sample_size_pe_category_telo_052823.csv')

## combine preexisting sample sizes
#read_csv("./select/sample_size_pe_category_051323.csv") |> rbind(read_csv('./select/sample_size_pe_category_telo_052823.csv')) |> write_csv('./select/sample_size_pe_category_060623.csv')
#write_csv('./select/sample_size_pe_category_060623.csv')




