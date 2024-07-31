
## Chirag J Patel
## perform meta analysis across survey years for output of pe_tables_adj.R
## 02/04/24

library(tidyverse)
library(broom)
library(metafor)
library(getopt)
library(logger)
library(this.path)


spec <- matrix(c(
  'model_number', 'm', 1, "integer"
), byrow=TRUE, ncol=4);
opt <- getopt(spec)

MODEL_NUMBER <- opt$model_number #'adjusted''

log_info("Model number: {MODEL_NUMBER}")
#load('pe_out_pscale.Rdata')
#load('pp_summary_032024.Rdata')

expos <- pe |> filter(grepl('expo', term)) |> filter(model_number == MODEL_NUMBER)
expos_nested <- NULL
expos_nested <- expos |> group_by(exposure, phenotype, term) |> nest()

log_info("Models to run: {nrow(expos_nested)}")

rma_safely <- safely(rma.uni)
pe_meta <- expos_nested |> mutate( ## change this to future_map?
  meta_model = map(data, ~rma_safely(yi=.x$estimate, sei=.x$std.error, method="REML"))
) |> select(-data)

pe_meta <- pe_meta |> mutate(
  mod = map(meta_model,  pluck, "result"),
  error = map(meta_model, pluck, "error"),
  tidied = map(mod, tidy),
  glanced = map(mod, glance)
) |> select(-mod) |> mutate(model_number = MODEL_NUMBER)


log_info("Done... saving file now.")
##
remove(expos_nested)
fileout <- str_glue('pp_meta_model_{MODEL_NUMBER}.rds')
saveRDS(pe_meta, file=fileout)

