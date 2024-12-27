## split the todo file into 10 equal lists of exposures

library(tidyverse)
n_groups <- 5
#ss_file <- '../select/sample_size_pe_category_060623.csv'
#ss_file <- '../select/sample_size_pe_category_0824.csv'
#ss_file <- '../select/sample_size_ubiome_e_112424.csv'
ss_file <- '../select/sample_size_ubiome_y_112424.csv'
to_do <- read_csv(ss_file) |> filter(grepl("genus", evarname)) |> filter(grepl("relative", evarname)) |> group_by(evarname) |> tally()

index <- rep(1:n_groups, floor(nrow(to_do)/n_groups))
index <- c(index, 1:( nrow(to_do)-length(index) ))
to_do <- to_do |> mutate(index=index)

for (i in 1:n_groups) {
  group_data <- to_do |> filter(index == i) |> select(evarname)
  file_name <- paste0(i, ".csv")
  write_csv(group_data, file_name, col_names=F)
}
