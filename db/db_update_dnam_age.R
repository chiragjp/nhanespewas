## Chirag
## 8/30/25
## prepare the DNAm (epigenetic clocks) data for cycles A (1999–2000) and B (2001–2002)

library(tidyverse)
library(DBI)
library(RSQLite)
library(glue)

# ---- CONFIG ----
path_to_db <- "../db/nhanes_031725.sqlite"

# RDS paths for the two cycles (adjust if needed)
path_a <- "../download/DNAm Age/DNAM_A.rds"  # 1999–2000
path_b <- "../download/DNAm Age/DNAM_B.rds"  # 2001–2002

# Expected columns
dnam_vars <- c(
  "SEQN", "XY_Estimation","HorvathAge","HannumAge","SkinBloodAge","PhenoAge",
  "GDF15Mort","B2MMort","CystatinCMort","TIMP1Mort","ADMMort","PAI1Mort",
  "LeptinMort","PACKYRSMort","CRPMort","logA1CMort","GrimAgeMort","GrimAge2Mort",
  "HorvathTelo","YangCell","ZhangAge","LinAge","WeidnerAge","VidalBraloAge",
  "DunedinPoAm","CD8TPP","CD4TPP","Nkcell","Bcell","MonoPP","NeuPP","WTDN4YR"
)

# Minimal descriptions
desc_lookup <- c(
  SEQN = "NHANES participant identifier",
  XY_Estimation = "Predicted biological sex from DNAm (XY estimation)",
  HorvathAge = "DNAm age (Horvath clock)",
  HannumAge = "DNAm age (Hannum clock)",
  SkinBloodAge = "DNAm age (Skin & Blood clock)",
  PhenoAge = "DNAm PhenoAge",
  GDF15Mort = "DNAm GrimAge surrogate: GDF15",
  B2MMort = "DNAm GrimAge surrogate: Beta-2 microglobulin",
  CystatinCMort = "DNAm GrimAge surrogate: Cystatin C",
  TIMP1Mort = "DNAm GrimAge surrogate: TIMP1",
  ADMMort = "DNAm GrimAge surrogate: Adrenomedullin",
  PAI1Mort = "DNAm GrimAge surrogate: PAI-1",
  LeptinMort = "DNAm GrimAge surrogate: Leptin",
  PACKYRSMort = "DNAm smoking pack-years surrogate",
  CRPMort = "DNAm GrimAge surrogate: CRP",
  logA1CMort = "DNAm GrimAge surrogate: log(HbA1c)",
  GrimAgeMort = "DNAm GrimAge",
  GrimAge2Mort = "DNAm GrimAge2",
  HorvathTelo = "DNAm telomere length (Horvath Telo)",
  YangCell = "DNAm mitotic age (PCHi-C/Yang replication clock)",
  ZhangAge = "DNAm age (Zhang clock)",
  LinAge = "DNAm age (Lin clock)",
  WeidnerAge = "DNAm age (Weidner clock)",
  VidalBraloAge = "DNAm age (Vidal-Bralo clock)",
  DunedinPoAm = "DNAm pace of aging (DunedinPoAm)",
  CD8TPP = "Estimated immune proportion: CD8 T cells",
  CD4TPP = "Estimated immune proportion: CD4 T cells",
  Nkcell = "Estimated immune proportion: NK cells",
  Bcell = "Estimated immune proportion: B cells",
  MonoPP = "Estimated immune proportion: Monocytes",
  NeuPP = "Estimated immune proportion: Neutrophils",
  WTDN4YR = "DNAm 4-year examination weight"
)

# ---- CONNECT DB ----
con <- DBI::dbConnect(RSQLite::SQLite(), dbname = path_to_db)
on.exit(DBI::dbDisconnect(con), add = TRUE)

table_description <- dplyr::tbl(con, "table_names_epcf")
variable_description <- dplyr::tbl(con, "variable_names_epcf")

# ---- READ & COERCE DATA (RDS) ----
read_dnam_rds <- function(path) {
  df <- readRDS(path)

  # Ensure required columns exist (create missing as NA)
  missing_cols <- setdiff(dnam_vars, names(df))
  if (length(missing_cols) > 0) df[missing_cols] <- NA

  # Keep only expected columns in defined order
  df <- df[, dnam_vars]

  # Coerce types (SEQN integer; others numeric)
  df <- df |>
    mutate(SEQN = as.integer(SEQN))

  num_cols <- setdiff(dnam_vars, "SEQN")
  df[num_cols] <- lapply(df[num_cols], function(x) suppressWarnings(as.numeric(x)))

  df
}

dnam_a <- read_dnam_rds(path_a)
dnam_b <- read_dnam_rds(path_b)

# ---- WRITE DNAM TABLES ----
DBI::dbWriteTable(con, "DNAM_A", dnam_a, append = FALSE, overwrite = TRUE)
DBI::dbWriteTable(con, "DNAM_B", dnam_b, append = FALSE, overwrite = TRUE)

# ---- UPDATE table_names_epcf ----
table_descs <- tibble::tibble(
  Data.File.Name = c("DNAM_A", "DNAM_B"),
  component = "LABORATORY",
  series = c("A", "B"),
  epcf = "p"
)

new_table_description <- table_description |> collect() |> bind_rows(table_descs)
DBI::dbWriteTable(con, "table_names_epcf", new_table_description, append = FALSE, overwrite = TRUE)

# ---- UPDATE variable_names_epcf ----
make_var_desc <- function(var, file, begin, end) {
  tibble::tibble(
    Variable.Name = var,
    Variable.Description = unname(ifelse(!is.na(desc_lookup[var]),
                                         desc_lookup[var],
                                         glue("DNAm-derived measure: {var}"))),
    Data.File.Name = file,
    Data.File.Description = "DNA methylation (epigenetic clocks and surrogates)",
    Begin.Year = begin,
    EndYear = end,
    Component = "Laboratory",
    Use.Constraints = "None"
  )
}

vars <- dnam_vars
variable_descs_a <- purrr::map_dfr(vars, make_var_desc,
                                   file = "DNAM_A", begin = 1999, end = 2000) |> filter(Variable.Name != "SEQN")
variable_descs_b <- purrr::map_dfr(vars, make_var_desc,
                                   file = "DNAM_B", begin = 2001, end = 2002) |> filter(Variable.Name != "SEQN")

new_variable_description <- variable_description |> collect() |> bind_rows(variable_descs_a, variable_descs_b)
DBI::dbWriteTable(con, "variable_names_epcf", new_variable_description, append = FALSE, overwrite = TRUE)
tbl(con, "variable_names_epcf") |> filter(Data. == 'HorvathAge')
#dbDisconnect(con)

tbl(con, "variable_names_epcf") |> filter(Data.File.Name == 'DNAM_A' | Data.File.Name == 'DNAM_B') |> collect() |> write_csv("../select/select_mage_pheno_variables.csv")


## append to domain
old_domain <-read_csv("../select/variable_domains_ep_2.csv")
new_variable_description <- bind_rows(variable_descs_a, variable_descs_b) |> filter(Variable.Name != 'WTDN4YR')
new_variable_description <- new_variable_description |> mutate(varname = Variable.Name, category = 'aging', subcategory = NA)
new_variable_description <- new_variable_description |> mutate(series = ifelse(Begin.Year == 1999, 'A', 'B'))
new_variable_description <- new_variable_description |> mutate(epcf = 'p')

new_domain <- old_domain |> rbind(new_variable_description)
write_csv(new_domain, "../select/variable_domains_ep_2.csv")
