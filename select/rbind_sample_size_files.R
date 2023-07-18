## combine sample_size files

read_csv('./sample_size_pe.csv') |> rbind(read_csv('./sample_size_alphape.csv')) |> write_csv('./sample_size_pe_040523.csv')