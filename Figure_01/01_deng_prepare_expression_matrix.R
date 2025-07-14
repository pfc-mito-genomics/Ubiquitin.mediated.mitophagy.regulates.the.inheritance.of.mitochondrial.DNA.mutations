# Load necessary libraries
library(dplyr)
library(data.table)
library(purrr)
library(future)
library(future.apply)

# Set up parallel processing
plan(multisession)  # Use multicore processing

# Download Deng data
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE45719&format=file", "~/Mount/Suffolk/WorkGenomics/cdr42/MF_selection/deng/GSE45719_RAW/")
input_dir <- "~/Mount/Suffolk/WorkGenomics/cdr42/MF_selection/deng/GSE45719_RAW/"
output_file <- "~/Mount/Suffolk/WorkGenomics/cdr42/MF_selection/deng/combined_gene_counts.txt"

# List all .txt.gz files in the input directory
files <- list.files(input_dir, pattern = "_expression.txt.gz$", full.names = TRUE)

# Function to read and process a single file
process_file <- function(file) {
  # Extract base name for the column header
  base_name <- tools::file_path_sans_ext(basename(file))
  base_name <- sub("_expression", "", base_name)
  
  # Read and process the file
  df <- fread(file, header = TRUE, sep = "\t") %>%
    select(Gene_symbol = 1, reads = 4) %>% # Ensure correct column indices and names
    rename(!!base_name := reads)
  
  return(df)
}

# Read and process files in parallel
data_list <- future_lapply(files, process_file)

# Function to remove duplicates by summing reads
remove_duplicates <- function(df) {
  df %>%
    group_by(Gene_symbol) %>%
    summarize(across(everything(), sum, na.rm = TRUE))
}

# Apply the function to each data frame in the list
data_list <- lapply(data_list, remove_duplicates)

# Combine all data frames by Gene_symbol using a full join
combined_df <- reduce(data_list, function(x, y) {
  full_join(x, y, by = "Gene_symbol", relationship = "many-to-many")
})

# Clean up column names by removing .txt suffixes
colnames(combined_df) <- gsub("^.{11}", "", gsub("\\.txt$", "", colnames(combined_df)))

# Check the number of rows and columns in the combined data frame
print(dim(combined_df))
head(combined_df)
# Write to output file
fwrite(combined_df, output_file, sep = "\t")

cat("Processing complete. Results saved in", output_file, "\n")
