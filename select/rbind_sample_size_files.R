## combine sample_size files

#read_csv('./sample_size_pe.csv') |> rbind(read_csv('./sample_size_alphape.csv')) |> write_csv('./sample_size_pe_040523.csv')
library(tidyverse)
rbind(
  read_csv('sample_size_pe_category_A_0824.csv'),
  read_csv('sample_size_pe_category_B_0824.csv'),
  read_csv('sample_size_pe_category_C_0824.csv'),
  read_csv('sample_size_pe_category_D_0824.csv'),
  read_csv('sample_size_pe_category_E_0824.csv'),
  read_csv('sample_size_pe_category_F_0824.csv'),
  read_csv('sample_size_pe_category_G_0824.csv'),
  read_csv('sample_size_pe_category_H_0824.csv'),
  read_csv('sample_size_pe_category_I_0824.csv'),
  read_csv('sample_size_pe_category_J_0824.csv')
) |> write_csv("sample_size_pe_category_0824.csv")
