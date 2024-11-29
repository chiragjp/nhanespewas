# create  summary stats table for the shiny app

library(tidyverse)
library(DBI)
library(corrr)


ep_assoc_summary_across_models <- function(summary_stats, glanced_stats) {
  summary_stats_wide <- summary_stats |> pivot_wider(names_from = "model_type", values_from = c("estimate", "std.error", "statistic", "p.value"))
  summary_stats_wide <- summary_stats_wide |> mutate(estimate_diff = estimate_adjusted-estimate_unadjusted)
  adj_vs_base <- glanced |> select(-c(adj.r2, df.residual, null.deviance, df.null, deviance)) |> pivot_wider(names_from=model_type, values_from = c("rsq", "nobs", "AIC", "BIC"))
  adj_vs_base <- adj_vs_base |> mutate(rsq_adjusted_base_diff=rsq_adjusted-rsq_base, rsq_adjusted_diff = rsq_adjusted-rsq_unadjusted)
  summary_stats_wide |> left_join(adj_vs_base, by=c("evarname", "pvarname", "exposure_table_name", "phenotype_table_name"))
}

con <- DBI::dbConnect(RSQLite::SQLite(), dbname='./pe_summary_stats.sqlite')
varnames <- tbl(con, "variable_names_epcf")
adjusted_meta <- tbl(con, "adjusted_meta")
unadjusted_meta <- tbl(con, "unadjusted_meta")
adjusted_meta <- adjusted_meta |> left_join(unadjusted_meta |> select(evarname, pvarname, expo_name, vartype, estimate, p.value) |> rename(estimate_unadjusted=estimate, p.value_unadjusted=p.value), by=c("evarname", "pvarname", "expo_name", "vartype"))
mvr2 <- tbl(con, 'mvr2') |> mutate(mv = mve_rsq-base_rsq)
pe <- tbl(con, "pe")
glanced <- tbl(con, "glanced")
variable_domain <- tbl(con, "variable_domain")

expos <- pe |> filter(term %like% 'expo%')
expos_wide <- ep_assoc_summary_across_models(expos, glanced)
expos_wide <- expos_wide |> left_join(varnames, by=c("evarname"="Variable.Name", "exposure_table_name"="Data.File.Name"))
expos_wide <- expos_wide |> left_join(varnames |> select(Variable.Name, Data.File.Name, Variable.Description, Data.File.Description),
                                      by=c("pvarname"="Variable.Name", "phenotype_table_name"="Data.File.Name"))

expos_wide <- expos_wide |> select(-Use.Constraints) |> rename(e_data_file_desc=Data.File.Description.x, p_data_file_desc=Data.File.Description.y,
                                                               e_variable_description=Variable.Description.x,
                                                               p_variable_description=Variable.Description.y
)
expos_wide <- expos_wide |> collect()
expos_wide_summary <- expos_wide |> filter(term == 'expo' | term == 'expo1') |> group_by(evarname, pvarname) |> summarize(mean_adjusted_base_r2_diff = mean(rsq_adjusted_base_diff), mean_unadjusted_r2_diff=mean(rsq_adjusted_diff), total_n = sum(nobs_adjusted)) |> ungroup()
adjusted_meta <- adjusted_meta |> collect() |> left_join(expos_wide_summary, by=c("evarname", "pvarname"))


p_variable_domain <- variable_domain |> filter(epcf == 'p') |> collect() |> group_by(Variable.Name) |> summarise(pvardesc=first(Variable.Description),pcategory=first(category),psubcategory=first(subcategory))
e_variable_domain <- variable_domain |> filter(epcf == 'e') |> collect() |> group_by(Variable.Name) |> summarise(evardesc=first(Variable.Description),ecategory=first(category),esubcategory=first(subcategory))

adjusted_meta <- adjusted_meta |> left_join(p_variable_domain, by=c("pvarname"="Variable.Name"))
adjusted_meta <- adjusted_meta |> left_join(e_variable_domain, by=c("evarname"="Variable.Name"))

expos_wide <- expos_wide |> left_join(p_variable_domain, by=c("pvarname"="Variable.Name"))


#expos_wide <- expos_wide |> select(-c(Variable.Description.y, Data.File.Description.y)) |> rename(Data.File.Description=Data.File.Description.x, Variable.Description=Variable.Description.x)

adjusted_meta_2 <- adjusted_meta |> filter(nobs >= 2)
adjusted_meta_2 <- adjusted_meta_2 |> ungroup() |>  mutate(pval_BY=p.adjust(p.value, method="BY"), pvalue_bonferroni=p.adjust(p.value, method="bonferroni"))
adjusted_meta_2 <- adjusted_meta_2 |> mutate(sig_levels = case_when(
  pvalue_bonferroni < 0.05 ~ 'Bonf.<0.05',
  pval_BY < 0.05 ~ 'BY<0.05',
  TRUE ~ '> BY & Bonf.'
))


e_summary <- adjusted_meta_2 |> group_by(evarname) |> arrange(pvalue_bonferroni) |>
  summarize(mean_r2=mean(mean_adjusted_base_r2_diff, na.rm=T),  mean_estimate=mean(abs(estimate), na.rm=T),
            median_pvalue=median(p.value, na.rm=T), n_sig=sum(pvalue_bonferroni < 0.05, na.rm=T),
            n_tests=sum(!is.na(pvalue_bonferroni)),  median_i.squared=median(i.squared, na.rm=T),
            max_r2=first(mean_adjusted_base_r2_diff), max_pvarname=first(pvarname) , max_estimate=first(estimate), max_p.value=first(p.value)) |> mutate(n_sig_pct=n_sig/n_tests)


p_summary <- adjusted_meta_2 |> group_by(pvarname) |> arrange(pvalue_bonferroni) |>
  summarize(mean_r2=mean(mean_adjusted_base_r2_diff, na.rm=T), mean_estimate=mean(abs(estimate), na.rm=T),
            median_pvalue=median(p.value, na.rm=T), n_sig=sum(pvalue_bonferroni < 0.05, na.rm=T),
            n_tests=sum(!is.na(pvalue_bonferroni)),  median_i.squared=median(i.squared, na.rm=T),
            max_r2=first(mean_adjusted_base_r2_diff), max_evarname=first(evarname) , max_estimate=first(estimate), max_p.value=first(p.value)) |> mutate(n_sig_pct=n_sig/n_tests)

## deeper summary by group
p_group_summary <- adjusted_meta_2 |> unite(p_scategory, c(pcategory, psubcategory)) |> group_by(p_scategory) |> arrange(pvalue_bonferroni) |>
  summarize(mean_r2=mean(mean_adjusted_base_r2_diff, na.rm=T),  mean_estimate=mean(abs(estimate), na.rm=T),
            median_pvalue=median(p.value, na.rm=T), n_sig=sum(pvalue_bonferroni < 0.05, na.rm=T),
            n_tests=sum(!is.na(pvalue_bonferroni)),  median_i.squared=median(i.squared, na.rm=T),
            max_r2=first(mean_adjusted_base_r2_diff), max_evarname=first(evarname) , max_estimate=first(estimate), max_p.value=first(p.value)) |> mutate(n_sig_pct=n_sig/n_tests)


e_group_summary <- adjusted_meta_2 |> unite(e_scategory, c(ecategory, esubcategory)) |> group_by(e_scategory) |> arrange(pvalue_bonferroni) |>
  summarize(mean_r2=mean(mean_adjusted_base_r2_diff, na.rm=T),
            mean_abs_estimate=mean(abs(estimate), na.rm=T),
            mean_estimate=mean((estimate), na.rm=T),
            median_pvalue=median(p.value, na.rm=T),
            n_sig=sum(pvalue_bonferroni < 0.05, na.rm=T),
            n_tests=sum(!is.na(pvalue_bonferroni)),
            median_i.squared=median(i.squared, na.rm=T),
            max_r2=first(mean_adjusted_base_r2_diff),
            max_pvarname=first(pvarname),
            max_estimate=first(estimate), max_p.value=first(p.value)) |> mutate(n_sig_pct=n_sig/n_tests)



to_array <- adjusted_meta_2 |> filter(expo_name == 'expo', vartype == 'continuous') |>
  select(evarname, pvarname, estimate, p.value) |> mutate(estimate= ifelse(p.value >= 0.05, 0, estimate)) |>
  mutate(estimate = ifelse(is.na(estimate), 0, estimate))  |> select(-p.value) |> pivot_wider(names_from = pvarname, values_from = estimate)
exposure_correlation <- to_array |> select(-evarname) |> correlate(diagonal = 1) |> stretch(na.rm = TRUE, remove.dups = TRUE)

to_array <- adjusted_meta_2 |> filter(expo_name == 'expo', vartype == 'continuous') |> select(evarname, pvarname, estimate, p.value) |> mutate(estimate= ifelse(p.value >= 0.05, 0, estimate)) |> mutate(estimate = ifelse(is.na(estimate), 0, estimate))  |> select(-p.value) |> pivot_wider(names_from = evarname, values_from = estimate)
phenome_correlation <- to_array |> select(-pvarname) |> correlate(diagonal = 1) |> stretch(na.rm = TRUE, remove.dups = TRUE)



shiny_con <-  DBI::dbConnect(RSQLite::SQLite(), dbname='./pe_atlas_app/pe_summary_stats_shiny.sqlite')

dbWriteTable(shiny_con, 'p_variable_domain', p_variable_domain, temporary=F, overwrite=T)
dbWriteTable(shiny_con, 'e_variable_domain', e_variable_domain, temporary=F, overwrite=T)
dbWriteTable(shiny_con, 'exposure_correlation', exposure_correlation, temporary=F, overwrite=T)
dbWriteTable(shiny_con, 'phenotype_correlation', phenome_correlation, temporary=F, overwrite=T)
dbWriteTable(shiny_con, 'e_group_summary', e_group_summary, temporary=F, overwrite=T)
dbWriteTable(shiny_con, 'e_summary', e_summary, temporary=F, overwrite=T)
dbWriteTable(shiny_con, 'p_group_summary', p_group_summary, temporary=F, overwrite=T)
dbWriteTable(shiny_con, 'p_summary', p_summary, temporary=F, overwrite=T)
dbWriteTable(shiny_con, 'adjusted_meta_2', adjusted_meta_2, temporary=F, overwrite=T)
dbWriteTable(shiny_con, 'expos_wide', expos_wide, temporary=F, overwrite=T)
dbWriteTable(shiny_con, 'expos_wide_summary', expos_wide_summary, temporary=F, overwrite=T)

