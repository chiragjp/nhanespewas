
## Chirag J Patel
## perform meta analysis across survey years for output of pe_tables_adj.R
## 02/04/24

library(tidyverse)
library(getopt)
library(logger)
#devtools::load_all("..")


spec <- matrix(c(
  'model_number', 'm', 1, "integer"
), byrow=TRUE, ncol=4);
opt <- getopt(spec)

MODEL_NUMBER <- opt$model_number #'adjusted''


stanley_meta <- function(dat) {
  ## dat should have the statistic and std.error from the tidy() function for a pair of variables
  # Neither fixed nor random: weighted leastsquares meta-analysis T. D. Stanley and Hristos Doucouliagos, Stat Med 2013
  # Unrestricted weighted least squares represent medical research better than random effects in 67,308 Cochrane meta-analyses. TD Stanley, JCE 2023
  dat <- dat |> select(estimate, statistic, std.error, p.value) |> mutate(variance = std.error^2)
  k <- nrow(dat)
  if(k == 1) {
    return(list(estimate=dat$estimate, std.error=dat$std.error, zval = dat$statistic, pval=dat$p.value, pval.stouffer=dat$p.value , k=k, h.squared=NA, Q=NA, I.squared=NA))
  }
  mod <- lm(statistic ~ I(1/std.error)-1, dat)
  m <- nrow(dat) - 1
  mse <- 1/(nrow(dat)-1) * sum((fitted(mod) - dat$statistic)^2) # this is Q/k-1 and H^2
  sse <-  sum((fitted(mod) - dat$statistic)^2) # this is Q
  wt <- 1/(dat$variance)
  estimate <- sum(dat$estimate*wt)/sum(wt)
  var.stanley <- mse / (sum(1/dat$variance)) # same as a FE scaled by the H2
  se <- sqrt(var.stanley)
  zval <- estimate / se
  pval <- pt(abs(zval), df=m, lower.tail = F)*2
  pval.stouffer <- pnorm(abs(sum(dat$statistic)/(m+1)), lower.tail =F)*2
  i2 <- (mse-1)/mse
  i2 <- ifelse(i2 < 0, 0, i2)
  return(list(estimate=estimate, std.error=se, zval=zval, pval = pval, pval.stouffer=pval.stouffer, k=k, h.squared=mse, Q=sse, I.squared=i2))
}

log_info("Model number: {MODEL_NUMBER}")
#load('pe_out_pscale.Rdata')
#load('../pe_summary_020424/pe_summary_022524.Rdata')
load('../pe_summary_0824/gathered_by_series_0824.Rdata')
#load('pp_summary_032024.Rdata')

expos <- pe |> filter(grepl('expo', term)) |> filter(model_number == MODEL_NUMBER)
expos_nested <- NULL
expos_nested <- expos |> group_by(exposure, phenotype, term) |> nest()

log_info("Models to run: {nrow(expos_nested)}")

stanley_safety <- safely(stanley_meta)

pe_meta <- expos_nested |> mutate(
  meta_model = map(data, ~stanley_safety(.x))
) |> select(-data)


pe_meta <- pe_meta |> mutate(
  mod = map(meta_model,  pluck, "result"),
  error = map(meta_model, pluck, "error"),
  estimate.uwls = map_dbl(mod, pluck, "estimate"),
  std.error.uwls = map_dbl(mod, pluck, "std.error"),
  p.value.uwls = map_dbl(mod, pluck, "pval"),
  p.value.stouffer = map_dbl(mod, pluck, "pval.stouffer"),
  statistic.uwls = map_dbl(mod, pluck, "zval"),
  h.squared.uwls = map_dbl(mod, pluck, "h.squared"),
  Q.uwls = map_dbl(mod, pluck, "Q"),
  i.squared.uwls = map_dbl(mod, pluck, "I.squared"),
  k.uwls = map_dbl(mod, pluck, "k")
) |> select(-mod, -meta_model) |> mutate(model_number = MODEL_NUMBER)



log_info("Done... saving file now.")
##
remove(expos_nested)
fileout <- str_glue('pe_meta_model_uwls_{MODEL_NUMBER}.rds')
saveRDS(pe_meta, file=fileout)




