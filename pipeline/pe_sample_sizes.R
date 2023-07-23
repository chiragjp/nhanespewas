## Get the sample sizes per correlation and survey year

source('quantpe.R')
library(corrr)
library(future)

con <- DBI::dbConnect(RSQLite::SQLite(), dbname='nhanes_122322.sqlite')
e_variables <- read_csv('./select/select_expo_variables_3.csv')
p_variables <- read_csv('./select/select_pheno_variables_3.csv')

#p_variables <- read_csv("./select/select_ubiome_pheno_variables.csv") |> mutate(epcf='p')

variable_information_selected <- rbind(e_variables, p_variables)

merged_table <- function(con, tab1, tab2) {
  if(tab1 == tab2) {
    return(tbl(con, tab1) |> collect())
  }
  tbl(con, tab1) |> inner_join(tbl(con,tab2), by = "SEQN") |> collect()
}

sample_size_for_table <- function(m_table, evars, pvars) {
  sample_size_pair <- m_table |> pair_n() |> as_cordf()
  vars <- intersect(c(evars, pvars), colnames(m_table))
  sample_size_pair <- sample_size_pair |> dice(all_of(vars)) |> stretch(na.rm=F) |> rename(n=r)
  sample_size_pair <- sample_size_pair |> filter(x %in% evars, y %in% pvars) |> rename(evarname=x, pvarname=y)
  sample_size_pair
} 

sample_size_for_ep_tables <- function(con, eTable, pTable, variable_information) {
  log_info("Working on { eTable } x {pTable}")
 m_table <- merged_table(con, pTable, eTable)  |> collect()
 evars <- variable_information |> filter(Data.File.Name == eTable, epcf == 'e') |> select(Variable.Name) |> pull()
 pvars <- variable_information |> filter(Data.File.Name == pTable, epcf == 'p') |> select(Variable.Name) |> pull()
 sample_size_for_table(m_table, evars,pvars) |> mutate(e_table_name = eTable, p_table_name = pTable)
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
    for_e[[e]] <- map(p_tables, ~sample_size_for_ep_tables(con, e_tab_name, .x, var_information))
  }
  for_e |> bind_rows()
}




log_info("Sample size for survey:")
log_info("A")
samp_size_a %<-% ep_sample_sizes(con, e_variables |> filter(Begin.Year == 1999), p_variables |> filter(Begin.Year == 1999))
log_info("B")
samp_size_b %<-% ep_sample_sizes(con, e_variables |> filter(Begin.Year == 2001), p_variables |> filter(Begin.Year == 2001))
log_info("C")
samp_size_c <- ep_sample_sizes(con, e_variables |> filter(Begin.Year == 2003), p_variables |> filter(Begin.Year == 2003))
log_info("D")
samp_size_d <- ep_sample_sizes(con, e_variables |> filter(Begin.Year == 2005), p_variables |> filter(Begin.Year == 2005))
log_info("E")
samp_size_e <- ep_sample_sizes(con, e_variables |> filter(Begin.Year == 2007), p_variables |> filter(Begin.Year == 2007))
log_info("F")
samp_size_f <- ep_sample_sizes(con, e_variables |> filter(Begin.Year == 2009), p_variables |> filter(Begin.Year == 2009))
log_info("G")
samp_size_g <- ep_sample_sizes(con, e_variables |> filter(Begin.Year == 2011), p_variables |> filter(Begin.Year == 2011))
log_info("H")
samp_size_h <- ep_sample_sizes(con, e_variables |> filter(Begin.Year == 2013), p_variables |> filter(Begin.Year == 2013))
log_info("I")
samp_size_i <- ep_sample_sizes(con, e_variables |> filter(Begin.Year == 2015), p_variables |> filter(Begin.Year == 2015))
log_info("J")
samp_size_j <- ep_sample_sizes(con, e_variables |> filter(Begin.Year == 2017), p_variables |> filter(Begin.Year == 2017))

samp_sizes <- rbind(
  samp_size_a, samp_size_b, samp_size_c, samp_size_d, samp_size_e, samp_size_f, samp_size_g, samp_size_h, samp_size_i, samp_size_j
)

#samp_sizes |> write_csv('./select/sample_size_pe_category_041823.csv')
#samp_sizes |> write_csv('./select/sample_size_pe_category_051323.csv')
#rbind(samp_size_a, samp_size_b) |> filter(pvarname == 'TeloMean' | pvarname == 'TeloStd') |> write_csv('./select/sample_size_pe_category_telo_052823.csv')

## combine preexisting sample sizes
read_csv("./select/sample_size_pe_category_051323.csv") |> rbind(read_csv('./select/sample_size_pe_category_telo_052823.csv')) |> write_csv('./select/sample_size_pe_category_060623.csv')




