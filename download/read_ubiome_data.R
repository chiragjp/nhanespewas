## microbiome data
## this was released on 11/27/24
## chirag 
## 11/28/24


directory_to_ubiome <- '~/Dropbox (HMS)/projects/nhanespewas/download/ubiome/'
count_files <- dir(directory_to_ubiome, pattern = "*count")
relative_files <- dir(directory_to_ubiome, pattern = "*relative")
annotate_files <- dir(directory_to_ubiome, pattern = "*annotate")
## rb and rsv files

files <- c(count_files, relative_files)

table_list <- vector("list", length=length(files))
for(ii in 1:length(files)) {
  print(files[[ii]])
  tab <- read_tsv(file.path(directory_to_ubiome, files[[ii]]))
  table_list[[ii]] <- tab
  table_name <- sub(".txt","" , gsub("-", "_", files[[ii]]))
  
}
annotation_rb <- read_tsv(file.path(directory_to_ubiome, "dada2rb-taxonomy-annotate.txt"))
annotation_rsv <- read_tsv(file.path(directory_to_ubiome, "dada2rsv-taxonomy-annotate.txt"))

annotation_rb <- annotation_rb |> rename(taxonomy_silva_123 = `Taxonomy in SILVA v123`)
annotation_rsv <- annotation_rsv |> rename(taxonomy_silva_123 = `Taxonomy in SILVA v123`)


