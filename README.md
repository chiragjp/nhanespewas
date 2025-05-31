# Benchmarking exposome research: An atlas of exposome-phenome associations

04/07/25

## Summary

We build a comprehensive "atlas" of correlations between 278 phenotypes (e.g., body mass index, glucose, height, creatinine) and 651 exposures and behaviors (e.g., self-reported nutrient intake; blood lead levels; urinary phthalates; blood PFOA) across the National Health and Nutrition Examination Surveys (NHANES), estimating 127,817 associations in individuals from 1999-2017. For each survey, we associate the phenotype and exposure, and summarize the association size across the surveys. We assess the robustness of associations by estimating the false discovery rate, consistency with adjustments, and concordance across multiple waves of the surveys.

Specifically, the atlas provides a comprehensive list of exposure-phenome associations across the cohort, document statistical approaches used to associate exposures with phenotype, and systematically assesses replicability across independent cohorts.

<img src="img/pe_fig1.png" width="70%" height="70%"/>

## Phenome-Exposome Atlas and Browser

-   [PE Browsable Atlas](http://apps.chiragjpgroup.org/pe_atlas/)

## nhanespewas: run your own P-E Associations

### Installation

## Download the database
NHANES 1999-2017
https://doi.org/10.6084/m9.figshare.29182196.v1


## install the R package
library(devtools)
devtools::install_github('chiragjp/nhanespewas')

### Execute an E-P Association

```         
# 
library(nhanespewas)
# connect to nhanes data, built for PE analysis
pedata_con <- connect_pewas_data(PATH_TO_THE_DATABASE)

# conduct an association between c-reactive protein and cotinine
## -   See pe_quickstart.Rmd

crp_cot <- pe_flex_adjust("LBXCRP", "LBXCOT", pedata_con, scale_p=T)
# clean up
disconnect_pewas_data(pedata_con)
```

