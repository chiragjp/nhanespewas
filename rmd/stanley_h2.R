

## run the stanley method

dat <- dat |> filter(model_number == 2)



mod <- lm(statistic ~ I(1/std.error)-1, dat)
tidy(mod)
glance(mod)

mse <- 1/(nrow(dat)-1) * sum((fitted(mod) - dat$statistic)^2) # this is Q/k-1 or H^2
sse <-  sum((fitted(mod) - dat$statistic)^2) # this is Q


stanley_meta <- function(dat) {
  ## dat should have the statistic and std.error from the tidy() function
  mod <- lm(statistic ~ I(1/std.error)-1, dat)
  mse <- 1/(nrow(dat)-1) * sum((fitted(mod) - dat$statistic)^2) # this is Q/k-1 or H^2
  sse <-  sum((fitted(mod) - dat$statistic)^2) # this is Q
  return(list(h2=mse, Q=sse))
}
