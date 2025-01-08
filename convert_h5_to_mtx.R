####################################
### Script Description
####################################
# This R script converts .h5 files containing single-cell RNA-seq data into 
# the 10X Genomics matrix format (matrix.mtx.gz, barcodes.tsv.gz and features.tsv.gz). It takes a directory of .h5 files as input
# and creates a structured output directory where each .h5 file is processed 
# into its corresponding 10X output format.
#
# Required libraries:
# - Seurat: For handling and subsetting single-cell data.
# - DropletUtils: For reading and writing 10X matrix files.
# - hdf5r: For reading .h5 files.
#
# To run this script, use the command line with the following syntax:
# Rscript convert_h5_to_10x.R <input_directory> <output_directory>
# - <input_directory>: Path to the directory containing .h5 files.
# - <output_directory>: Path to the directory where the converted files will be saved.
####################################

####################################
### Required libraries
####################################
library(Seurat) # To subset matrix
library(DropletUtils) # (Bioconductor) to handle reading and writing mtx files and downsample matrix
library(hdf5r) # To read .h5 files
####################################

# Get input and output directories from command line arguments
args <- commandArgs(trailingOnly = TRUE)
data_dir <- args[1]  # Directory containing .h5 files
output_dir <- args[2]  # Directory for outputs

# Create a main output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE)

# Create a subdirectory for 10X output
tenx_output_dir <- file.path(output_dir, "h5_to_10X")
dir.create(tenx_output_dir, showWarnings = FALSE)

# List .h5 files
h5_files <- list.files(data_dir, pattern = "\\.h5$", full.names = TRUE)

# Check if any .h5 files were found
if (length(h5_files) == 0) {
  stop("No .h5 files found in the specified directory.")
}

# Process each .h5 file
for (file in h5_files) {
  
  # Create sample name by removing .h5
  sample_name <- sub("\\.h5$", "", basename(file))
  sample_path <- file.path(tenx_output_dir, sample_name)  # Use the new 10X directory
  
  # Print the sample name to show progress
  print(sample_name)
  
  # Check if output path already exists
  if (dir.exists(sample_path)) {
    message("Removing existing directory: ", sample_path)
    unlink(sample_path, recursive = TRUE)  # Remove the existing directory
  }
  
  # Load the count matrices
  counts <- tryCatch({
    Seurat::Read10X_h5(file, use.names = TRUE)
  }, error = function(e) {
    message("Error reading ", file, ": ", e$message)
    return(NULL)
  })
  
  # If counts were loaded successfully
  if (!is.null(counts)) {
    gene_names <- rownames(counts)
    
    # Convert to 10X format
    DropletUtils::write10xCounts(sample_path, counts, version = "3", gene.symbol = gene_names)
  }
}

print("Processing completed.")
