## read log files from slurm and parse them for error string



# Function to read a log file and search for the word "error"
search_errors_in_log <- function(log_file_path) {
  # Open the log file
  con <- file(log_file_path, "r")
  
  # Initialize a counter for errors
  error_count <- 0
  
  # Read the log file line by line
  while(TRUE) {
    line <- readLines(con, n = 1, warn = FALSE)
    
    # Break the loop if the end of the file is reached
    if(length(line) == 0) {
      break
    }
    
    # Search for the word "error" (case-insensitive)
    if(grepl("error", line, ignore.case = TRUE)) {
      # Print the line containing the word "error"
      cat("Error found: ", line, "\n")
      error_count <- error_count + 1
    }
  }
  
  # Close the file connection
  close(con)
  
  # Print the total number of errors found
  cat("Total number of errors found: ", error_count, "\n")
}

## Example usage
#log_file_path <- "path/to/your/logfile.log"
#search_errors_in_log(log_file_path)


# Function to traverse all log files ending with ".out" in a directory
search_errors_in_directory <- function(directory_path) {
  # Get the list of all files in the directory ending with ".out"
  log_files <- list.files(path = directory_path, pattern = "\\.out$", full.names = TRUE)
  
  # Loop through each log file and call search_errors_in_log
  for (log_file in log_files) {
    cat("Log file", log_file, "\n")
    search_errors_in_log(log_file)
  }
}

# Example usage
directory_path <- "."
search_errors_in_directory(directory_path)

#Log file ./BPXPLS.out 
# URXUCR - outlier, need to check what is up
#Log file ./LBXSNASI.out 
#Log file ./LBXSKSI.out 
#Log file ./BMXWT.out 
#Log file ./BMXWAIST.out 
#Log file ./BMXHT.out 
#Log file ./BMXLEG.out 