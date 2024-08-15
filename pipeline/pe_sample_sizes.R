## Get the sample sizes per correlation and survey year

library(getopt)
library(corrr)
source('db_paths.R')
library(DBI)
library(tidyverse)
library(logger)
spec <- matrix(c(
  'series', 's', 1, "character"
), byrow=TRUE, ncol=4);
opt <- getopt(spec)

series <- opt$series

con <- DBI::dbConnect(RSQLite::SQLite(), dbname=path_to_nhanes)
e_variables <- read_csv('../select/select_expo_variables_3.csv')
p_variables <- read_csv('../select/select_pheno_variables_3.csv')

#p_variables <- read_csv("./select/select_ubiome_pheno_variables.csv") |> mutate(epcf='p')

variable_information_selected <- rbind(e_variables, p_variables)

merged_table <- function(con, tab1, tab2) {
  if(tab1 == tab2) {
    return(tbl(con, tab1) |> collect())
  }
  tbl(con, tab1) |> inner_join(tbl(con,tab2), by = "SEQN") |> collect()
}

merged_table_dplyr <- function(con, tab1, tab2) {
  if(tab1 == tab2) {
    return(tbl(con, tab1))
  }
  tbl(con, tab1) |> inner_join(tbl(con,tab2), by = "SEQN")
}

sample_size_for_table <- function(m_table, evars, pvars) {
  sample_size_pair <- m_table |> pair_n() |> as_cordf()
  vars <- intersect(c(evars, pvars), colnames(m_table))
  sample_size_pair <- sample_size_pair |> dice(all_of(vars)) |> stretch(na.rm=F) |> rename(n=r)
  sample_size_pair <- sample_size_pair |> filter(x %in% evars, y %in% pvars) |> rename(evarname=x, pvarname=y)
  sample_size_pair
}


sample_size_for_table_dplyr <- function(m_table, evars, pvars) {
  ## to complete
  nn <- vector("list",length = length(evars)*length(pvars))
  i <- 1
  for(evari in 1:length(evars)) {
    for(pvari in 1:length(pvars)) {
      #m_table |> select(pvars[pvari], evars[evari])
      col1 <- pvars[pvari]
      col2 <- evars[evari]
      n <- m_table |> filter(!is.na(!!sym(col1)), !is.na(!!sym(col2))) |>
        summarize(nn = n()) |> pull(nn)
      nn[[i]] <- tibble(pvarname=col1, evarname=col2, n=n)
      i <- i + 1
    }
  }
  # return a tibble
  nn |> bind_rows()
}

sample_size_for_ep_tables <- function(con, eTable, pTable, variable_information) {
log_info("Working on { eTable } x {pTable}")
 m_table <- merged_table(con, pTable, eTable)  |> collect()
 evars <- variable_information |> filter(Data.File.Name == eTable, epcf == 'e') |> select(Variable.Name) |> pull()
 pvars <- variable_information |> filter(Data.File.Name == pTable, epcf == 'p') |> select(Variable.Name) |> pull()
 sample_size_for_table(m_table, evars,pvars) |> mutate(e_table_name = eTable, p_table_name = pTable)
}


sample_size_for_ep_tables_dplyr <- function(con, eTable, pTable, variable_information) {
  log_info("Working on { eTable } x {pTable}")
  m_table <- merged_table_dplyr(con, pTable, eTable)
  evars <- variable_information |> filter(Data.File.Name == eTable, epcf == 'e') |> select(Variable.Name) |> pull()
  pvars <- variable_information |> filter(Data.File.Name == pTable, epcf == 'p') |> select(Variable.Name) |> pull()
  sample_size_for_table_dplyr(m_table, evars,pvars) |> mutate(e_table_name = eTable, p_table_name = pTable)
}

ep_sample_sizes <- function(con, e_vars, p_vars) {
  ## per survey year
  e_tables <- e_vars |> select(Data.File.Name) |> pull() |> unique()
  p_tables <- p_vars |> select(Data.File.Name) |> pull() |> unique()
  var_information <- rbind(e_vars, p_vars)
  log_info("Num e_tables { length(e_tables) }")
  log_info("Num p_tables { length(p_tables) }")
  for_e <- vector("list", length(e_tables))
  for(e in 1:length(e_tables)) {
    e_tab_name <- e_tables[e]
    if((e %% 10) == 0) {
      log_info("Num e done { e }")
    }
    for_e[[e]] <- map(p_tables, ~sample_size_for_ep_tables_dplyr(con, e_tab_name, .x, var_information))
  }
  for_e |> bind_rows()
}


log_info("Sample size for survey:")
if(series == "A") {
  log_info("A")
  samp_size <- ep_sample_sizes(con, e_variables |> filter(Begin.Year == 1999), p_variables |> filter(Begin.Year == 1999))
} else if(series == "B") {
  log_info("B")
  samp_size <- ep_sample_sizes(con, e_variables |> filter(Begin.Year == 2001), p_variables |> filter(Begin.Year == 2001))
} else if(series == "C") {
  log_info("C")
  samp_size <- ep_sample_sizes(con, e_variables |> filter(Begin.Year == 2003), p_variables |> filter(Begin.Year == 2003))
} else if(series == "D") {
  log_info("D")
  samp_size <- ep_sample_sizes(con, e_variables |> filter(Begin.Year == 2005), p_variables |> filter(Begin.Year == 2005))
} else if(series == "E") {
  log_info("E")
  samp_size <- ep_sample_sizes(con, e_variables |> filter(Begin.Year == 2007), p_variables |> filter(Begin.Year == 2007))
} else if(series == "F") {
  log_info("F")
  samp_size <- ep_sample_sizes(con, e_variables |> filter(Begin.Year == 2009), p_variables |> filter(Begin.Year == 2009))
} else if(series == "G") {
  log_info("G")
  samp_size <- ep_sample_sizes(con, e_variables |> filter(Begin.Year == 2011), p_variables |> filter(Begin.Year == 2011))
} else if(series == "H") {
  log_info("H")
  samp_size <- ep_sample_sizes(con, e_variables |> filter(Begin.Year == 2013), p_variables |> filter(Begin.Year == 2013))
} else if(series == "I") {
  log_info("I")
  samp_size <- ep_sample_sizes(con, e_variables |> filter(Begin.Year == 2015), p_variables |> filter(Begin.Year == 2015))
} else if(series == "J") {
  log_info("J")
  samp_size <- ep_sample_sizes(con, e_variables |> filter(Begin.Year == 2017), p_variables |> filter(Begin.Year == 2017))
} else {
  exit("specify a valid series!")
}

samp_size |> write_csv(sprintf('sample_size_pe_category_%s_0824.csv', series))


#log_info("A")
#samp_size_a <- ep_sample_sizes(con, e_variables |> filter(Begin.Year == 1999), p_variables |> filter(Begin.Year == 1999))
#log_info("B")
#samp_size_b <- ep_sample_sizes(con, e_variables |> filter(Begin.Year == 2001), p_variables |> filter(Begin.Year == 2001))
#log_info("C")
#samp_size_c <- ep_sample_sizes(con, e_variables |> filter(Begin.Year == 2003), p_variables |> filter(Begin.Year == 2003))
#log_info("D")
#samp_size_d <- ep_sample_sizes(con, e_variables |> filter(Begin.Year == 2005), p_variables |> filter(Begin.Year == 2005))
#log_info("E")
#samp_size_e <- ep_sample_sizes(con, e_variables |> filter(Begin.Year == 2007), p_variables |> filter(Begin.Year == 2007))
#log_info("F")
#samp_size_f <- ep_sample_sizes(con, e_variables |> filter(Begin.Year == 2009), p_variables |> filter(Begin.Year == 2009))
#log_info("G")
#samp_size_g <- ep_sample_sizes(con, e_variables |> filter(Begin.Year == 2011), p_variables |> filter(Begin.Year == 2011))
#log_info("H")
#samp_size_h <- ep_sample_sizes(con, e_variables |> filter(Begin.Year == 2013), p_variables |> filter(Begin.Year == 2013))
#log_info("I")
#amp_size_i <- ep_sample_sizes(con, e_variables |> filter(Begin.Year == 2015), p_variables |> filter(Begin.Year == 2015))
#log_info("J")
#samp_size_j <- ep_sample_sizes(con, e_variables |> filter(Begin.Year == 2017), p_variables |> filter(Begin.Year == 2017))

#samp_sizes <- rbind(
#  samp_size_a, samp_size_b, samp_size_c, samp_size_d, samp_size_e, samp_size_f, samp_size_g, samp_size_h, samp_size_i, samp_size_j
#)

#samp_sizes |> write_csv('./select/sample_size_pe_category_0824.csv')



#samp_sizes |> write_csv('./select/sample_size_pe_category_041823.csv')
#samp_sizes |> write_csv('./select/sample_size_pe_category_051323.csv')
#rbind(samp_size_a, samp_size_b) |> filter(pvarname == 'TeloMean' | pvarname == 'TeloStd') |> write_csv('./select/sample_size_pe_category_telo_052823.csv')

## combine preexisting sample sizes
#read_csv("./select/sample_size_pe_category_051323.csv") |> rbind(read_csv('./select/sample_size_pe_category_telo_052823.csv')) |> write_csv('./select/sample_size_pe_category_060623.csv')

#write_csv("./sle")



dbDisconnect(con)
