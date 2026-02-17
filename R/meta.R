
# meta-analysis using stanley et al meta analysis approach
# Unrestricted weighted least squares represent medical research better than random effects in 67,308 Cochrane meta-analyses
# T.D. Stanley John P.A. Ioannidis Maximilian Maier, Hristos Doucouliagos Willem M. Otte, Frantisek Bartos JCE 2023


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
  pval.stouffer <- pnorm(abs(sum(dat$statistic)/sqrt(m+1)), lower.tail =F)*2
  i2 <- (mse-1)/mse
  i2 <- ifelse(i2 < 0, 0, i2)
  return(list(estimate=estimate, std.error=se, zval=zval, pval = pval, pval.stouffer=pval.stouffer, k=k, h.squared=mse, Q=sse, I.squared=i2))
}
