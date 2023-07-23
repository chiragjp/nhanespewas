## upload and download data from S3

download_pewas_nhanes_data <- function(path_to_database) {
  ## use S3?
  download.file("http://mywebsite-name.com/xxxx.sqlite", path_to_database)
  path_to_database
}

