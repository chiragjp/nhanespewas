library(dplyr)
library(readr)
library(nhanespewas)

ss_file <- system.file("extdata", "sample_size_pe_category_0824.csv", package = "nhanespewas")
ss <- read_csv(ss_file)

# 1) Build exposure list (same logic you already have)
exposures <- ss %>%
  select(evarname, e_table_name) %>%
  distinct() %>%
  count(evarname, name = "n_tables") %>%
  filter(n_tables >= 2) %>%
  transmute(evarname)

# 2) Split into ~10 equal chunks
K <- 10
exposures <- exposures %>% mutate(chunk = ntile(row_number(), K))

# 3) Write each chunk file (one column CSV with header evarname)
dir.create("exposure_chunks", showWarnings = FALSE)

exposures %>%
  group_by(chunk) %>%
  group_walk(~ write_csv(.x %>% select(evarname),
                         file = sprintf("exposure_chunks/exposures_%02d.csv", .y$chunk)))