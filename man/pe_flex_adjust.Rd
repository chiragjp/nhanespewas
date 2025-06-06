% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quantpe.R
\name{pe_flex_adjust}
\alias{pe_flex_adjust}
\title{Phenotype and Exposure Analysis under different modeling scenarios}
\usage{
pe_flex_adjust(
  pheno,
  exposure,
  adjustment_variables,
  con,
  series = NULL,
  logxform_p = T,
  logxform_e = T,
  scale_e = T,
  scale_p = F,
  pheno_table_name = NULL,
  expo_table_name = NULL,
  quantile_expo = NULL,
  exposure_levels = NULL,
  scale_type = 1
)
}
\arguments{
\item{pheno}{A character string specifying the phenotype variable.}

\item{exposure}{A character string specifying the exposure variable.}

\item{adjustment_variables}{A data frame containing adjustment variables and scenarios. It should have columns \code{variables} and \code{scenario}.}

\item{con}{A database connection object.}

\item{series}{A character vector specifying the series to consider (default is NULL).}

\item{logxform_p}{Logical, whether to log-transform the phenotype variable (default is TRUE).}

\item{logxform_e}{Logical, whether to log-transform the exposure variable (default is TRUE).}

\item{scale_e}{Logical, whether to scale the exposure variable (default is TRUE).}

\item{scale_p}{Logical, whether to scale the phenotype variable (default is FALSE).}

\item{pheno_table_name}{Optional character string specifying the phenotype table name (default is NULL).}

\item{expo_table_name}{Optional character string specifying the exposure table name (default is NULL).}

\item{quantile_expo}{Optional vector specifying quantiles for the exposure variable (default is NULL).}

\item{exposure_levels}{Optional vector specifying levels for the exposure variable (default is NULL).}

\item{scale_type}{Numeric: scale type 1 is mean and SD; scale type 2 == CLR scale type 3 is IVT}
}
\value{
A list containing the following elements:
\itemize{
\item \code{log_p}: Logical, whether the phenotype was log-transformed.
\item \code{log_e}: Logical, whether the exposure was log-transformed.
\item \code{scaled_p}: Logical, whether the phenotype was scaled.
\item \code{scaled_e}: Logical, whether the exposure was scaled.
\item \code{unweighted_n}: The number of unweighted observations.
\item \code{phenotype}: The phenotype variable.
\item \code{series}: The series information from the data.
\item \code{exposure}: The exposure variable.
\item \code{models}: A list of models run based on the scenarios.
\item \code{base_models}: A list of base models run for comparison.
\item \code{adjustment_variables}: The adjustment variables used in the models, a tibble: see adjustment_models
\item \code{demographic_breakdown}: A table with demographic breakdown of the survey design.
}
}
\description{
This function performs flexible adjustment for phenotype and exposure analysis using specified adjustment variables. It retrieves necessary data, handles weighting, and runs models based on the specified scenarios.
}
\details{
The function performs the following steps:
\itemize{
\item Retrieves table names for the phenotype and exposure variables.
\item Ensures the phenotype and exposure variables are collected in the same survey.
\item Gets table names for each series and binds the tables.
\item Handles weights by calling \code{figure_out_multiyear_weight}.
\item Transforms and names the phenotype and exposure variables if specified.
\item Filters data and creates a survey design object.
\item Performs demographic breakdown of the survey design.
\item Runs models based on the specified scenarios and adjustment variables.
}
}
\examples{
\dontrun{
con <- DBI::dbConnect(RSQLite::SQLite(), dbname = "nhanes.db")
adjustment_vars <- data.frame(variables = c("AGE", "SEX"), scenario = c("A", "A"))
result <- pe_flex_adjust(pheno = "BMXBMI", exposure = "LBXCOT", adjustment_variables = adjustment_vars, con = con)
}
}
