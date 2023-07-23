
## Chirag J Patel
## perform meta analysis across survey years

MODEL_TYPE <- 'adjusted'
VARTYPE <- 'categorical'

library(tidyverse)
library(broom)
library(metafor)
library(getopt)
library(logger)
library(this.path)

spec <- matrix(c(
  'model_type', 'm', 1, "character",
  'vartype', 'v', 1, "character"
), byrow=TRUE, ncol=4);
opt <- getopt(spec)

MODEL_TYPE <- opt$model_type #'adjusted'
VARTYPE <- opt$vartype #'categorical'

log_info("{MODEL_TYPE} & {VARTYPE}")
load('pe_out_pscale.Rdata')

expos <- pe |> filter(grepl('expo', term)) |> filter(model_type == MODEL_TYPE, vartype == VARTYPE)
expos_nested <- NULL
if(VARTYPE == 'categorical') {
  expos_nested <- expos |> group_by(evarname, pvarname, term) |> nest()   
} else {
  expos_nested <- expos |> group_by(evarname, pvarname) |> nest()   
}



rma_safely <- safely(rma.uni)
pe_meta <- expos_nested |> mutate(
  meta_model = map(data, ~rma_safely(yi=.x$estimate, sei=.x$std.error, method="REML"))
) |> select(-data)

pe_meta <- pe_meta |> mutate(
  mod = map(meta_model,  pluck, "result"),
  error = map(meta_model, pluck, "error"),
  tidied = map(mod, tidy),
  glanced = map(mod, glance)
) |> select(-mod)


## 

remove(expos_nested)

fileout <- str_glue('pe_meta_{MODEL_TYPE}_{VARTYPE}.rds')
saveRDS(pe_meta, file=fileout) 

