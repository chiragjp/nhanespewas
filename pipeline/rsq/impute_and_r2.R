# Chirag 
# 4/14/25
# estimate R2
# see create_mv_data.R

library(tidyverse)
library(mice)
library(finalfit)
library(getopt)

MISSING_PCT_THRESHOLD <- 40 # default
spec <- matrix(c(
  'pheno_file', 'p', 1, "character",
  'missing_pct', 'm', 2, "double"
), byrow=TRUE, ncol=4);
opt <- getopt(spec)

pfile <- opt$pheno_file

big_data <- read_rds(pfile)

big_missing_gl <- missing_glimpse(big_data$big_data)
bd <- big_data$big_data
summ_stat <- big_data$summ_stats
phenotype <- big_data$phenotype
if(!is.na(opt$missing_pct)) {
  MISSING_PCT_THRESHOLD <- opt$missing_pct
}


exposures_selected <- big_missing_gl |> filter(missing_percent < MISSING_PCT_THRESHOLD) |> select(label) |> inner_join(summ_stat, by=c("label"="evarname"))

if(nrow(exposures_selected) == 0) {
  message("no exposures above the threshold of missingness to impute, qutting")
  quit("no", 0)
}

cols <- c("SDMVPSU", "WTMEC2YR", "INDFMPIR", "SDDSRVYR",
          "ETHNICITY_NONHISPANICWHITE", "ETHNICITY_MEXICAN", "ETHNICITY_OTHERHISPANIC",  "ETHNICITY_OTHER", "ETHNICITY_NONHISPANICBLACK", "AGE_SQUARED", "RIAGENDR", "RIDAGEYR",
          "EDUCATION_LESS9",  "EDUCATION_9_11", "EDUCATION_HSGRAD", "EDUCATION_AA", "EDUCATION_COLLEGEGRAD", "BORN_INUSA")


bd_exposure_tibble <- bd |> select(all_of(c("SEQN", phenotype, cols, exposures_selected$label)))
bd_exposure_tibble <- bd_exposure_tibble |> filter(!is.na(!!as.symbol(phenotype))) |> filter(!is.na(INDFMPIR))

bd_m <- mice(bd_exposure_tibble, maxit = 0, m = 10)
dat <- complete(bd_m)

dat <- dat |> rename(pheno=all_of(phenotype))

mod1 <- lm(scale(pheno) ~ RIAGENDR + RIDAGEYR + AGE_SQUARED + ETHNICITY_NONHISPANICBLACK + ETHNICITY_OTHER + ETHNICITY_MEXICAN + ETHNICITY_OTHERHISPANIC + EDUCATION_LESS9 + EDUCATION_9_11 + EDUCATION_HSGRAD, dat)
mod2 <- lm(scale(pheno) ~ ., dat |> select(-c(SEQN, SDMVPSU, BORN_INUSA, EDUCATION_AA, ETHNICITY_NONHISPANICWHITE)))


# Multivariate R2
mv_diff <- glance(mod2)$r.squared - glance(mod1)$r.squared

# sum of univariate
uni_r2 <- exposures_selected |> summarize(sum(rsq_adj)) 

to_save <- list(phenotype=phenotype, rsq_exposures=mv_diff, rsq_uni=uni_r2, number_exposures_in_model=nrow(exposures_selected), exposures_selected=exposures_selected, missing_threshold=MISSING_PCT_THRESHOLD)
print(mv_diff)
print(uni_r2)
print(nrow(exposures_selected))
print(MISSING_PCT_THRESHOLD)
file_out <- sprintf("%s_imp_R2.rds", phenotype)
save_rds(to_save, file=file_out)

