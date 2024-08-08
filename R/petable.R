## Chirag J Patel
## 07/14/23
## petable.R
## machinery to associate a real-valued outcome with a categorical or real-valued exposure


surveyVariables <- c('WTMEC2YR', 'WTMEC4YR', 'SDMVPSU', 'SDMVSTRA', 'SDDSRVYR')

seriesBeginYearMap <- function(seriesName) {
    seriesBeginYear <- dplyr::case_when(
    seriesName == 'A' ~ 1999,
    seriesName == 'B' ~ 2001,
    seriesName == 'C' ~ 2003,
    seriesName == 'D' ~ 2005,
    seriesName == 'E' ~ 2007,
    seriesName == 'F' ~ 2009,
    seriesName == 'G' ~ 2011,
    seriesName == 'H' ~ 2013,
    seriesName == 'I' ~ 2015,
    seriesName == 'J' ~ 2017,
    seriesName == 'K' ~ 2019,
    seriesName == 'L' ~ 2021,
    TRUE ~ NA_real_
  )
  seriesBeginYear
}



get_path_to_extdata_database<- function() {
  path_to_db <- system.file("extdata", "nhanes_pewas_a-d.sqlite", package = "nhanespewas")
}


#' Connect to the PEWAS NHANES database
#'
#' This function establishes a connection to the PEWAS NHANES database .
#'
#' @param path_to_data The file path to the participant-level NHANES data. Each table in the database corresponds to a table in the NHANES.
#'
#' @return A DBI connection object to the NHANES database.
#' @export
#' @examples
#' \dontrun{
#' con <- connect_pewas_data("/path/to/database.sqlite")
#' }
#' @importFrom DBI dbConnect
#' @importFrom RSQLite SQLite
connect_pewas_data <- function(path_to_data = NULL) {

  ## check if the database exists
  if (is.null(path_to_data)) {
    path_to_data <- get_path_to_extdata_database()
  } else if(!file.exists(path_to_data)) {
    stop("nhanes pewas database does not exist")
  }
  con <- DBI::dbConnect(RSQLite::SQLite(), dbname=path_to_data)
  # check version date and type of database (e.g., summary stats or NHANES raw)
  if(DBI::dbExistsTable(con, "description")) {
    desc <- dplyr::tbl(con, "description")
    return(con)
  } else{
    stop("nhanes pewas database is not in correct format")
  }

}

#' Disconnect from a NHANES database
#'
#' This function disconnects from a database given a DBI connection object.
#'
#' @param dbConn A NHANES DBI connection object.
#'
#' @return Invisible NULL. The connection will be completely severed.
#' @export
#' @examples
#' \dontrun{
#' disconnect_pewas_data(dbConn)
#' }
#' @importFrom DBI dbDisconnect
disconnect_pewas_data <- function(dbConn) {
  DBI::dbDisconnect(dbConn)
}

#' Retrieve and merge tables from a database based on the series name, exposure variable name, and phenotype variable name
#'
#' @param pvarname Character. The phenotype variable name.
#' @param evarname Character. The exposure variable name.
#' @param seriesName Character. The series name.
#' @param con The database connection.
#' @param pheno_table_name Character. The optional phenotype table name. If NULL (default), the function will attempt to find the table based on the `pvarname`.
#' @param expo_table_name Character. The optional exposure table name. If NULL (default), the function will attempt to find the table based on the `evarname`.
#'
#' @return A list containing the merged table, the exposure table, the phenotype table, and the series name.
#' @examples
#' \dontrun{
#' conn <- connect_pewas_data(...)
#' tables <- get_tables("LBXGLU", "LBXGTC", "C", conn, "L10AM_C", "L45VIT_C")
#' }
#' @export
get_tables <- function(pvarname, evarname, seriesName, con, pheno_table_name=NULL, expo_table_name=NULL) {
  variables <- dplyr::tbl(con, "variable_names_epcf")
  table_names <- dplyr::tbl(con, "table_names_epcf")
  series_begin_year <- seriesBeginYearMap(seriesName)
  demo_table_name <- table_names |> dplyr::filter(series == seriesName & component == 'DEMO') |> dplyr::collect()
  demo <- dplyr::tbl(con, demo_table_name$Data.File.Name)
  ## get table for the phenotype
  logger::log_info("Getting tables for {pvarname} ~ {evarname}")

  p_table <- tibble::tibble()
  if(!is.null(pheno_table_name)) {
    #p_table_name <- p_table_name |> filter(Data.File.Name == pheno_table_name)
    p_table <- tbl(con, pheno_table_name)
  } else{
    p_table_name <- variables |> dplyr::filter(Variable.Name == pvarname & Begin.Year == series_begin_year) |> dplyr::collect()
    pheno_table_name <- p_table_name$Data.File.Name
    p_table <- dplyr::tbl(con, pheno_table_name)
  }

  logger::log_info("Pheno table {pheno_table_name} has {p_table |> count() |> pull(n) } rows")
  ## get table for exposure
  e_table <- tibble::tibble()
  if(!is.null(expo_table_name)) {
    e_table <- dplyr::tbl(con, expo_table_name)
  } else {
    e_table_name <- variables |> dplyr::filter(Variable.Name == evarname & Begin.Year == series_begin_year) |> dplyr::collect()
    expo_table_name <- e_table_name$Data.File.Name
    e_table <- dplyr::tbl(con, expo_table_name)
  }

  if(expo_table_name == pheno_table_name) {
    p_table <- p_table |> dplyr::select(SEQN, pvarname)
  }

  etable_nr <- e_table |> dplyr::count() |> dplyr::pull(n)
  logger::log_info("Exposure table {expo_table_name} has {etable_nr} rows")
  small_tab <- demo |> dplyr::inner_join(p_table, by="SEQN") |> dplyr::inner_join(e_table, by="SEQN") |> dplyr::collect()
  logger::log_info("Merged table has {small_tab |> count() |> pull(n) } rows")
  return(list(merged_tab=small_tab, e_table=e_table, p_table=p_table, series=seriesName))
}

#' Get Table Names for a Given Variable Name
#'
#' This function fetches the table names from a database for a given variable name
#' and optionally for a specific series. It retrieves the `Data.File.Name` for each unique
#' year, and if there are multiple tables for a given year, it selects the table with the most rows.
#'
#' @param con A database connection object.
#' @param varname The name of the variable to search for in the `variable_names_epcf` table.
#' @param series Optional. The series to filter the tables by their `Begin.Year`. Default is NULL.
#'
#' @return A data frame with the selected `Data.File.Name` and `Begin.Year`.
#'
#' @details This function performs the following steps:
#' \enumerate{
#'   \item Filters the `variable_names_epcf` table to get the `Data.File.Name` and `Begin.Year` for the given `varname`.
#'   \item If `series` is provided, filters the tables by the year associated with the series.
#'   \item Groups the tables by `Begin.Year` and counts the number of tables for each year.
#'   \item For each year, if there is only one table, it selects that table. If there are multiple tables, it counts the number of rows in each table and selects the table with the most rows.
#' }
#'
#' @examples
#' \dontrun{
#' # Assuming you have a database connection `con`:
#' con <- dbConnect(RSQLite::SQLite(), ":memory:")
#'
#' # Get table names for a variable name "your_variable_name":
#' get_table_names_for_varname(con, "your_variable_name")
#'
#' # Get table names for a variable name "your_variable_name" and series "your_series":
#' get_table_names_for_varname(con, "your_variable_name", series = "your_series")
#' }
#'
#' @import dplyr DBI logger
#' @export
get_table_names_for_varname <- function(con, varname, series=NULL) {
  variables <- dplyr::tbl(con, "variable_names_epcf")
  all_tables <- variables |> dplyr::filter(Variable.Name == varname) |> select(Data.File.Name, Begin.Year) |> collect()
  if(!is.null(series)) {
    all_tables <- all_tables |> dplyr::filter(Begin.Year==seriesBeginYearMap(series))
  }
  table_names_per_year <-  all_tables |> group_by(Begin.Year) |> count() |> ungroup()
  ## if multiple tables are there for a series Begin.Year, choose the one with the largest N
  table_names <- vector(mode="list", length = nrow(table_names_per_year))
  for(table_name_row in 1:nrow(table_names_per_year)) {
    num_surveys <- (table_names_per_year |> dplyr::slice(table_name_row)) |> pull(n)
    yr <- table_names_per_year |> dplyr::slice(table_name_row) |> pull(Begin.Year)
    if(num_surveys == 1) {
      table_names[[table_name_row]] <- all_tables |>
        filter(Begin.Year == yr) |>
        select(Data.File.Name, Begin.Year)
    } else {
      # need to count the rows and get the table that has to the most data
      potential_tables <- all_tables |> filter(Begin.Year == yr) |> pull(Data.File.Name)
      nrow_per_table <- vector(mode="list", length = length(potential_tables))
      for(i in 1:length(potential_tables)) {
        tab_name <- potential_tables[[i]]
        if(!(tab_name %in% dbListTables(con))) {
          next;
        }
        nrows <- tbl(con, tab_name) |> filter(!is.na(varname)) |> tally() |> pull(n)
        logger::log_info("Number of rows in { tab_name }: { nrows } ")
        nrow_per_table[[i]] <- tibble(Data.File.Name=tab_name, nrows_in_table=nrows)
      }
      table_names[[table_name_row]] <- nrow_per_table |> dplyr::bind_rows() |> filter(nrows_in_table == max(nrows_in_table)) |>
        slice_head() |> select(Data.File.Name) |> mutate(Begin.Year = yr)
    }
  }
  return(table_names |> dplyr::bind_rows())
}


#' Retrieve and merge demographic, exposure, and phenotype tables from a database
#'
#' @param con The database connection.
#' @param pheno_table_name Character. The phenotype table name.
#' @param expo_table_name Character. The exposure table name.
#'
#' @return A list containing the merged table, the exposure table, the phenotype table, and the series name.
#' @examples
#' \dontrun{
#' conn <- connect_pewas_data()
#' results <- get_expo_pheno_tables(conn, pheno_table_name="L10AM_C", expo_table_name="L45VIT_C")
#' }
#' @export
get_expo_pheno_tables <- function(con, pheno_table_name, expo_table_name) {
  table_names <- dplyr::tbl(con, "table_names_epcf")
  seriesName <- table_names |> dplyr::filter(Data.File.Name == pheno_table_name) |> dplyr::pull(series)
  logger::log_info("Series of phenotype is { seriesName } ")
  demo_table_name <- table_names |> dplyr::filter(component == 'DEMO', series==seriesName) |> dplyr::collect() |> dplyr::pull(Data.File.Name)
  logger::log_info("Demographics table is { demo_table_name } ")
  demo <- dplyr::tbl(con, demo_table_name)
  e_table <- dplyr::tbl(con, expo_table_name)
  e_table_nrow <- e_table |> collect() |> nrow()
  logger::log_info("Exposure table {expo_table_name} has { e_table_nrow } rows")
  p_table <- dplyr::tbl(con, pheno_table_name)
  if(pheno_table_name == expo_table_name) { # hack to preserve data structure
    p_table <- p_table |> dplyr::select(SEQN)
  }
  p_table_nrow <- p_table |> collect() |> nrow()
  logger::log_info("Pheno table {pheno_table_name} has {p_table_nrow } rows")
  small_tab <- demo |> dplyr::inner_join(p_table, by="SEQN") |> dplyr::inner_join(e_table, by="SEQN") |> dplyr::collect()
  small_tab_nrow <- small_tab |> collect() |> nrow()
  logger::log_info("Merged table has { small_tab_nrow } rows")
  return(list(merged_tab=small_tab, e_table=e_table, p_table=p_table, series=seriesName))
}


#' Retrieve and merge demographic, and any two tables from the nhanes database
#'
#' @param con The database connection.
#' @param table_name1 Character. table name 1.
#' @param table_name2 Character. table name 2
#'
#' @return A list containing the merged table, the table 1, the table 2, and the series name.
#' @examples
#' \dontrun{
#' conn <- connect_pewas_data()
#' results <- get_expo_pheno_tables(conn, pheno_table_name="L10AM_C", expo_table_name="L45VIT_C")
#' }
#' @export
get_x_y_tables <- function(con, table_name1, table_name2) { ##
  table_names <- dplyr::tbl(con, "table_names_epcf")
  seriesName <- table_names |> dplyr::filter(Data.File.Name == table_name1) |> dplyr::pull(series)
  logger::log_info("Series is { seriesName } ")
  demo_table_name <- table_names |> dplyr::filter(component == 'DEMO', series==seriesName) |> dplyr::collect() |> dplyr::pull(Data.File.Name)
  logger::log_info("Demographics table is { demo_table_name } ")
  demo <- dplyr::tbl(con, demo_table_name)
  table1 <- dplyr::tbl(con, table_name1)
  table1_nrow <- table1 |> dplyr::collect() |> nrow()
  logger::log_info("Table 1 {table_name1} has { table1_nrow } rows")
  table2 <- dplyr::tbl(con, table_name2)
  if(table_name1 == table_name2) { # preserve the cols
    table2 <- table1 |> dplyr::select(SEQN)
  }
  table2_nrow <- table2 |> dplyr::collect() |> nrow()
  logger::log_info("Table 2 {table_name2} has {table2_nrow } rows")
  small_tab <- demo |> dplyr::inner_join(table1, by="SEQN") |> dplyr::inner_join(table2, by="SEQN") |> dplyr::collect()
  small_tab_nrow <- small_tab |> dplyr::collect() |> nrow()
  logger::log_info("Merged table has { small_tab_nrow } rows")
  return(list(merged_tab=small_tab, table1=table1, table2=table2, series=seriesName))
}



bind_tables <- function(con, table_name_list) {
  tables <- vector(mode="list", length=length(table_name_list))
  for(t in 1:length(table_name_list)) {
    table_name <- table_name_list[t]
    tables[[t]] <- dplyr::tbl(con, table_name) |> dplyr::collect()
  }
  deep_table <- bind_rows(tables)
}


#' Get Multiyear X and Y Tables -- possibly delete this, or use sparingly without weighting
#'
#' This function retrieves and merges demographic, table 1, and table 2 data for multiyear NHANES datasets from a specified database connection. It ensures the tables are appropriately joined and returns the merged result along with individual tables and series information.
#'
#' @param con A database connection object.
#' @param table_list_1 A character vector of table names for the first set of tables to be merged.
#' @param table_list_2 A character vector of table names for the second set of tables to be merged.
#' @return A list containing the following elements:
#' \itemize{
#'   \item `merged_tab`: A data frame of the merged tables.
#'   \item `table1`: A data frame of the first set of tables.
#'   \item `table2`: A data frame of the second set of tables.
#'   \item `series`: A data frame with series information.
#' }
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Retrieves table names from the `table_names_epcf` table.
#'   \item Obtains series information based on `table_list_1`.
#'   \item Retrieves demographic table names corresponding to the series.
#'   \item Binds demographic tables, `table_list_1`, and `table_list_2`.
#'   \item Merges the demographic data with the two sets of tables by the `SEQN` identifier.
#'   \item Collects the final merged table and returns it along with individual tables and series information.
#' }
#' @examples
#' \dontrun{
#' con <- DBI::dbConnect(RSQLite::SQLite(), dbname = "nhanes.db")
#' table_list_1 <- c("table1_1999", "table1_2001")
#' table_list_2 <- c("table2_1999", "table2_2001")
#' result <- get_multiyear_x_y_tables(con, table_list_1, table_list_2)
#' }
#' @export
get_multiple_x_y_tables <- function(con, table_list_1, table_list_2) {
  ## make sure table_list_1 and 2 are equal in size
  if(length(table_list_1) != length(table_list_2)) {
    stop("The lists of the tables to be binded must be of the same size!")
  }
  table_names <- dplyr::tbl(con, "table_names_epcf") |> dplyr::collect()
  ##get series
  table_series <- table_names |> right_join(tibble(Data.File.Name = table_list_1)) |> select(series)
  demo_table_names <- table_names |> dplyr::filter(component == "DEMO") |> right_join(table_series) |> pull(Data.File.Name)

  demo <- bind_tables(con, demo_table_names)
  table_1 <- bind_tables(con, table_list_1)
  table_2 <- bind_tables(con, table_list_2)
  small_tab <- demo |> dplyr::inner_join(table_1, by="SEQN") |> dplyr::inner_join(table_2, by="SEQN") |> dplyr::collect()
  return(list(merged_tab=small_tab, table1=table_1, table2=table_2, series=table_series))
}

get_x_y_tables_as_list <- function(con, table_list_1, table_list_2) {
  if(length(table_list_1) != length(table_list_2)) {
    stop("The lists of the tables to be binded must be of the same size!")
  }

  table_names_1 <- dplyr::tbl(con, "table_names_epcf") |> select(series, Data.File.Name) |>
    right_join(tibble(Data.File.Name = table_list_1), copy=T) |> dplyr::collect() |> rename(table_1=Data.File.Name)
  table_names_2 <- dplyr::tbl(con, "table_names_epcf") |>  select(series, Data.File.Name) |>
    right_join(tibble(Data.File.Name = table_list_2), copy=T) |> dplyr::collect() |> rename(table_2=Data.File.Name)
  table_names_series <- table_names_1 |> inner_join(table_names_2, by="series")
  #return(table_names_series)
  tab_obj_by_series <- vector(mode="list", length=nrow(table_names_series))
  for(rw in 1:nrow(table_names_series)) {
    tab_obj_by_series[[rw]] <- get_x_y_tables(con, table_names_series$table_1[rw], table_names_series$table_2[rw])
  }
  tab_obj_by_series
}

rbind_x_y_tables <- function(tab_obj_by_series, keep_xy_tables=F) {
  list_of_tibbles <- function(nested_list, name) {
    new_list <- map(nested_list, function(x) {x[[name]] |> dplyr::collect() })
  }
  if(keep_xy_tables) {
    tab_obj <- list(merged_tab=list_of_tibbles(tab_obj_by_series, "merged_tab") |> dplyr::bind_rows(),
                    table1=list_of_tibbles(tab_obj_by_series, "table1") |> dplyr::bind_rows(),
                    table2=list_of_tibbles(tab_obj_by_series, "table2") |> dplyr::bind_rows()
    )
  } else {
    tab_obj <- list(merged_tab=list_of_tibbles(tab_obj_by_series, "merged_tab") |> dplyr::bind_rows())
  }
  tab_obj$series <- map_chr(tab_obj_by_series, function(x) {x[["series"]]})
  tab_obj
}



#' Identify the appropriate weight column from the exposure and phenotype tables
#'
#' This function identifies the appropriate weight column from the exposure (e_table) and phenotype (p_table) tables, contained within the list object produced by the `get_tables()` function. The weight column is determined based on the following criteria:
#' 1. If a column named "WTDRD1" exists, it is selected.
#' 2. If not, and a column ending in "2YR" exists, that is selected.
#' 3. If neither of the above conditions are met, the first weight column is selected.
#' If there is no weight column in either table, a default weight column ("WTMEC2YR") is used.
#' Once the weight column is identified, a new column "wt" is added to the merged table (merged_tab), which is equivalent to the weight column.
#'
#' @param get_tables_obj A list object produced by the `get_tables()` function. Contains the exposure table (e_table), phenotype table (p_table), merged table (merged_tab), and series name.
#' @return The input list object with an additional column "wt" in the merged table (merged_tab) corresponding to the identified weight column.
#' @examples
#' \dontrun{
#' conn <- connect_pewas_data()
#' get_tables_results <- get_tables("LBXGLU", "LBXGTC", "C", conn, "L10AM_C", "L45VIT_C")
#' weighted_tables <- figure_out_weight(get_tables_results)
#' }
#' @export
figure_out_weight <- function(get_tables_obj) {
  ##  need to select weights that are not the bootstrap weights
  weight_name_demo <- 'WTMEC2YR'
  e_table <- get_tables_obj$e_table
  p_table <- get_tables_obj$p_table
  weight_name_e <- colnames(e_table)[grep("^WT", colnames(e_table))]
  weight_name_p <- colnames(p_table)[grep("^WT", colnames(p_table))]

  ## filter here for a single weight
  ## if e_weight is empty, move on
  ## if e_weight  has dietary weight (wtdrd1), select that (vs. the 2 day weight)
  ## else if e_weight has 2yr weight select that
  ## else select the first weight

  if(length(weight_name_e) > 1) {
    if('WTDRD1' %in% weight_name_e) { ## dietary variable
      weight_name_e <- 'WTDRD1'
    } else if(any(grepl('2YR$', weight_name_e))) {
      weight_name_e <- weight_name_e[grep('2YR$', weight_name_e)]
    } else {
      weight_name_e <- weight_name_e[1]
    }
  }

  if(length(weight_name_p) > 1) {
    if('WTDRD1' %in% weight_name_p) {
      weight_name_p <- 'WTDRD1'
    } else if(any(grepl('2YR$', weight_name_p))) {
      weight_name_p <- weight_name_p[grep('2YR$', weight_name_p)]
    } else {
      weight_name_p <- weight_name_p[1]
    }
  }

  logger::log_info("p weight name: { weight_name_p }")
  logger::log_info("e weight name: { weight_name_e }")
  etable_nrows <- e_table |> dplyr::count() |> pull(n)
  ptable_nrows <- p_table |> dplyr::count() |> pull(n)

  if(rlang::is_empty(weight_name_e) & rlang::is_empty(weight_name_p)) {
    logger::log_info("no weights in e or p table")
    weight_name <- weight_name_demo
  } else if(!rlang::is_empty(weight_name_e) & !rlang::is_empty(weight_name_p)) {
    logger::log_info("weights in both e and p table")

    if(etable_nrows < ptable_nrows) {
      weight_name <-  weight_name_e
    } else {
      weight_name <-  weight_name_p
    }

    if(weight_name_e == weight_name_p) {
      weight_name <- sprintf('%s.y', weight_name)
    }

  } else if(!rlang::is_empty(weight_name_p)) {
    logger::log_info("weight in p table")
    weight_name <- weight_name_p
  } else if (!rlang::is_empty(weight_name_e)) {
    logger::log_info("weight in e table")
    weight_name <- weight_name_e
  } else {
    weight_name <- weight_name_demo
  }
  logger::log_info("final weight name: { weight_name }")
  get_tables_obj$merged_tab <- get_tables_obj$merged_tab |> dplyr::mutate(wt = !!as.name(weight_name))
  get_tables_obj$weight_name <- weight_name
  get_tables_obj
}


#' Figure Out Multiyear Weight
#'
#' This function calculates appropriate weights for multiyear survey data, specifically handling special cases where both series A and B are present. It updates the `weight_name` and adjusts weights in the merged table accordingly.
#'
#' @param get_tables_obj A list object that contains the data tables and relevant metadata. This object should include:
#' \itemize{
#'   \item `series`: A data frame that contains the series information.
#'   \item `weight_name`: The name of the original weight variable.
#'   \item `merged_tab`: The merged data table that contains the survey data.
#' }
#' @return The updated `get_tables_obj` list with recalculated weights.
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Calls `figure_out_weight` to initialize weights.
#'   \item Checks if both series A and B are present. If so, it logs the use of 4-year weights and updates the weight variable name.
#'   \item If the original weight name is `WTDRD1`, it uses the special 4-year dietary weight `WTDR4YR`.
#'   \item Updates the weights in the `merged_tab`:
#'   \itemize{
#'     \item If only series A and B are present, sets the weight to the new 4-year weight.
#'     \item Otherwise, adjusts weights for series 1 and 2 and averages weights across all series.
#'   }
#'   \item If both series A and B are not present, logs that it is averaging 2-year weights across the surveys and divides the weight by the number of series.
#' }
#' @examples
#' \dontrun{
#' get_tables_obj <- figure_out_weight(get_tables_obj)
#' updated_obj <- figure_out_multiyear_weight(get_tables_obj)
#' }
#' @export
figure_out_multiyear_weight <- function(get_tables_as_list_obj) {
  series_weight <- tibble()
  for(ii in 1:length(get_tables_as_list_obj)) {
    get_tables_as_list_obj[[ii]]$e_table <- get_tables_as_list_obj[[ii]]$table2
    get_tables_as_list_obj[[ii]]$p_table <- get_tables_as_list_obj[[ii]]$table1
    get_tables_as_list_obj[[ii]] <- figure_out_weight(get_tables_as_list_obj[[ii]])
    series_weight <- rbind(series_weight, tibble(weight_name=get_tables_as_list_obj[[ii]]$weight_name, series=get_tables_as_list_obj[[ii]]$series))
  }
  logger::log_info("series: { series_weight$series } ")
  n <- length(series_weight$series)
  ab <- c("A", "B") %in% str_to_upper(series_weight$series)
  if(!any(ab) | sum(ab) == 1) {
    # easy scenario
    logger::log_info("averaging 2 year weights across { n } surveys")
    get_tables_obj <- rbind_x_y_tables(get_tables_as_list_obj)
    get_tables_obj$merged_tab <- get_tables_obj$merged_tab |> dplyr::mutate(orig_wt=wt) # save the old 2 year weight
    get_tables_obj$merged_tab <- get_tables_obj$merged_tab |> dplyr::mutate(wt=wt/n)
    return(get_tables_obj)
  } else {
    logger::log_info("Weights including A+B")
    for(table_index in 1:length(get_tables_as_list_obj)) {
        orig_weight_name <- get_tables_as_list_obj[[table_index]]$weight_name
        new_4yr_weight_name <- paste0(substr(orig_weight_name, 1, str_length(orig_weight_name) - 3), "4YR")
        ## special case
        if(orig_weight_name == 'WTDRD1') {
          new_4yr_weight_name <- "WTDR4YR"
        }
        get_tables_as_list_obj[[table_index]]$merged_tab <- get_tables_as_list_obj[[table_index]]$merged_tab |> dplyr::mutate(new_wt = ifelse(SDDSRVYR == 1 | SDDSRVYR == 2, !!sym(new_4yr_weight_name), wt))
    }
    get_tables_obj <- rbind_x_y_tables(get_tables_as_list_obj)
    get_tables_obj$merged_tab <- get_tables_obj$merged_tab |> dplyr::mutate(orig_wt=wt) |> mutate(wt=new_wt) |> select(-new_wt) # save the old 2 year weight
    if(length(get_tables_as_list_obj) > 2) {
      logger::log_info("averaging 2 year weights across { n } surveys including A+B")
      get_tables_obj$merged_tab <- get_tables_obj$merged_tab |> dplyr::mutate(wt=ifelse(SDDSRVYR == 1 | SDDSRVYR == 2,2*wt/n,wt/n))
    }
    return(get_tables_obj)
  }
}









#' Retrieve and merge demographic, multi-variable exposure, and phenotype tables from a database
#'
#' This function retrieves and merges tables from a database based on provided table names. Specifically, it merges demographic, multiple exposure variables and phenotype tables. The function logs the series of the phenotype, the name of the demographics table, the number of rows in each exposure and phenotype table, and finally the number of rows in the merged table.
#'
#' @param con The database connection.
#' @param pheno_table_name Character. The name of the phenotype table.
#' @param expo_table_names Character vector. The names of the exposure tables.
#'
#' @return A list containing the merged table, the final exposure table after joining all exposure tables, the phenotype table, and the series name.
#' @examples
#' \dontrun{
#' conn <- connect_pewas_data()
#' results <- get_mv_expo_pheno_tables(conn, "L10AM_C", c("L45VIT_C", "L06COT_C"))
#' }
get_mv_expo_pheno_tables <- function(con, pheno_table_name, expo_table_names) {
  table_names <- dplyr::tbl(con, "table_names_epcf")
  seriesName <- table_names |> dplyr::filter(Data.File.Name == pheno_table_name) |> dplyr::pull(series)
  logger::log_info("Series of phenotype is { seriesName } ")
  demo_table_name <- table_names |> dplyr::filter(component == 'DEMO', series==seriesName) |> dplyr::collect() |> dplyr::pull(Data.File.Name)
  logger::log_info("Demographics table is { demo_table_name } ")
  demo <- dplyr::tbl(con, demo_table_name) |> dplyr::collect()

  e_table <- NULL
  for(ii in 1:length(expo_table_names)) {
    e_table_lcl <- dplyr::tbl(con, expo_table_names[ii]) |> dplyr::collect()
    logger::log_info("Exposure table {expo_table_names[ii]} has {e_table_lcl |> count() |> pull(n) } rows")
    if(ii == 1) {
      e_table <- e_table_lcl
      next;
    }
    prev_colnames <- colnames(e_table)
    cols_to_keep <- c(setdiff(colnames(e_table_lcl), prev_colnames), "SEQN")
    e_table <- e_table |> dplyr::inner_join(e_table_lcl |> dplyr::select(tidyselect::all_of(cols_to_keep)), by="SEQN")
  }

  p_table <- dplyr::tbl(con, pheno_table_name) |> dplyr::collect()

  logger::log_info("Pheno table {pheno_table_name} has {p_table |> count() |> pull(n) } rows")
  small_tab <- demo |> dplyr::inner_join(p_table, by="SEQN") |> dplyr::inner_join(e_table, by="SEQN")
  logger::log_info("Merged table has { small_tab |> count() |> pull(n) } rows")

  return(list(merged_tab=small_tab, e_table=e_table, p_table=p_table, series=seriesName))
}




