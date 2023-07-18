
library(tidyverse)
library(DBI)
## extract the meta_analysis and the per-cohort PE tibbles

path_to_summary_stats <- './pe_summary_stats.sqlite'
summary_stats_con <- DBI::dbConnect(RSQLite::SQLite(), dbname=path_to_summary_stats)

ep_assoc_summary_across_models <- function(summary_stats, glanced_stats) {
  summary_stats_wide <- summary_stats |> pivot_wider(names_from = "model_type", values_from = c("estimate", "std.error", "statistic", "p.value")) 
  summary_stats_wide <- summary_stats_wide |> mutate(estimate_diff = estimate_adjusted-estimate_unadjusted)
  adj_vs_base <- glanced |> select(-c(adj.r2, df.residual, null.deviance, df.null, deviance)) |> pivot_wider(names_from=model_type, values_from = c("rsq", "nobs", "AIC", "BIC"))
  adj_vs_base <- adj_vs_base |> mutate(rsq_adjusted_base_diff=rsq_adjusted-rsq_base, rsq_adjusted_diff = rsq_adjusted-rsq_unadjusted)
  summary_stats_wide |> left_join(adj_vs_base, by=c("evarname", "pvarname", "exposure_table_name", "phenotype_table_name"))
}


varnames <- tbl(summary_stats_con, "variable_names_epcf")
adjusted_meta <- tbl(summary_stats_con, "adjusted_meta")
unadjusted_meta <- tbl(summary_stats_con, "unadjusted_meta")
adjusted_meta <- adjusted_meta |> left_join(unadjusted_meta |> select(evarname, pvarname, expo_name, vartype, estimate, p.value) |> rename(estimate_unadjusted=estimate, p.value_unadjusted=p.value), by=c("evarname", "pvarname", "expo_name", "vartype"))
pe <- tbl(summary_stats_con, "pe")
glanced <- tbl(summary_stats_con, "glanced")
variable_domain <- tbl(summary_stats_con, "variable_domain")


expos <- pe |> filter(term %like% 'expo%')
expos_wide <- ep_assoc_summary_across_models(expos, glanced)
expos_wide <- expos_wide |> left_join(varnames, by=c("evarname"="Variable.Name", "exposure_table_name"="Data.File.Name"))
expos_wide <- expos_wide |> left_join(varnames |> select(Variable.Name, Data.File.Name, Variable.Description, Data.File.Description), 
                                      by=c("pvarname"="Variable.Name", "phenotype_table_name"="Data.File.Name"))

expos_wide <- expos_wide |> select(-Use.Constraints) |> rename(e_data_file_desc=Data.File.Description.x, p_data_file_desc=Data.File.Description.y,
                                                               e_variable_description=Variable.Description.x, 
                                                               p_variable_description=Variable.Description.y
)

# expos_wide
expos_wide <- expos_wide |> collect()
expos_wide_summary <- expos_wide |> filter(term == 'expo' | term == 'expo1') |> group_by(evarname, pvarname) |> summarize(mean_adjusted_base_r2_diff = mean(rsq_adjusted_base_diff), mean_unadjusted_r2_diff=mean(rsq_adjusted_diff), total_n = sum(nobs_adjusted)) |> ungroup()
adjusted_meta <- adjusted_meta |> collect() |> left_join(expos_wide_summary, by=c("evarname", "pvarname"))


p_variable_domain <- variable_domain |> filter(epcf == 'p') |> collect() |> group_by(Variable.Name) |> summarise(pvardesc=first(Variable.Description),
                                                                                                                 pcategory=first(category),
                                                                                                                 psubcategory=first(subcategory))
e_variable_domain <- variable_domain |> filter(epcf == 'e') |> collect() |> group_by(Variable.Name) |> summarise(pvardesc=first(Variable.Description),
                                                                                                                 ecategory=first(category),
                                                                                                                 esubcategory=first(subcategory))
# adjusted_meta
adjusted_meta <- adjusted_meta |> left_join(p_variable_domain, by=c("pvarname"="Variable.Name"))
adjusted_meta <- adjusted_meta |> left_join(e_variable_domain, by=c("evarname"="Variable.Name"))

expos_wide <- expos_wide |> left_join(p_variable_domain, by=c("pvarname"="Variable.Name"))
expos_wide <- expos_wide |> left_join(e_variable_domain, by=c("evarname"="Variable.Name"))


adjusted_meta_2 <- adjusted_meta |> filter(nobs >= 2)
#adjusted_meta_2
adjusted_meta_2 <- adjusted_meta_2 |> ungroup() |>  mutate(pval_BY=p.adjust(p.value, method="BY"), pvalue_bonferroni=p.adjust(p.value, method="bonferroni"))
adjusted_meta_2 <- adjusted_meta_2 |> mutate(sig_levels = case_when(
  pvalue_bonferroni < 0.05 ~ 'Bonf.<0.05',
  pval_BY < 0.05 ~ 'BY<0.05',
  TRUE ~ '> BY & Bonf.'
))


dbWriteTable(summary_stats_con, 'adjusted_meta_2', adjusted_meta_2, overwrite=T, append=F)
dbWriteTable(summary_stats_con, 'expos_wide', expos_wide, overwrite=T, append=F)

dbDisconnect(summary_stats_con)
