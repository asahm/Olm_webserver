
create_folders_by_gene <- function(directory, dataframe, index_range, annotations_df) {
  # Create www directory if it doesn't exist
  www_dir <- "www"
  if (!dir.exists(www_dir)) {
    dir.create(www_dir, recursive = TRUE)
  }
  
  # Ensure the data directory exists under www
  full_directory <- file.path(www_dir, directory)
  if (!dir.exists(full_directory)) {
    dir.create(full_directory, recursive = TRUE)
  }
  
  # Check if the column "Gene ID" exists in the dataframe
  if (!"Gene ID" %in% colnames(dataframe)) {
    stop("The dataframe does not contain a column named 'Gene ID'.")
  }
  
  # Check if the required columns exist in the annotations dataframe
  required_columns <- c("gene_id", 
                        "prank_best_fas_translation", 
                        "prank_best_fas_translation_view_annotations", 
                        "prank_best_fas_renamed", 
                        "prank_best_fas_view_annotations")
  missing_columns <- setdiff(required_columns, colnames(annotations_df))
  if (length(missing_columns) > 0) {
    stop(paste("The annotations dataframe is missing the following columns:", paste(missing_columns, collapse = ", ")))
  }
  
  # Subset the dataframe based on the given index range
  if (any(index_range > nrow(dataframe)) || any(index_range < 1)) {
    stop("Index range is out of bounds.")
  }
  
  subset_df <- dataframe[index_range, ]
  
  # Create folders for each value in the "Gene ID" column
  for (gene_id in subset_df$`Gene ID`) {
    gene_folder <- file.path(full_directory, as.character(gene_id))
    if (!dir.exists(gene_folder)) {
      dir.create(gene_folder)
    }
    
    # Create subfolders "translation" and "nucleic_acid"
    translation_folder <- file.path(gene_folder, "translation")
    nucleic_acid_folder <- file.path(gene_folder, "nucleic_acid")
    
    if (!dir.exists(translation_folder)) {
      dir.create(translation_folder)
    }
    if (!dir.exists(nucleic_acid_folder)) {
      dir.create(nucleic_acid_folder)
    }
    
    # Find the corresponding row in the annotations dataframe
    annotation_row <- annotations_df[startsWith(annotations_df$gene_id, gene_id), ]
    if (nrow(annotation_row) == 0) {
      warning(paste("No annotation found for Gene ID:", gene_id))
      next
    }
    
    # Extract paths from the annotation row
    aa_alignment_path <- annotation_row$prank_best_fas_translation
    aa_alignment_annotation_path <- annotation_row$prank_best_fas_translation_view_annotations
    na_alignment_path <- annotation_row$prank_best_fas_renamed
    na_alignment_annotation_path <- annotation_row$prank_best_fas_view_annotations
    
    # Print the paths to the console
    cat("Gene ID:", gene_id, "\n")
    cat("  AA Alignment Path:", aa_alignment_path, "\n")
    cat("  AA Alignment Annotation Path:", aa_alignment_annotation_path, "\n")
    cat("  NA Alignment Path:", na_alignment_path, "\n")
    cat("  NA Alignment Annotation Path:", na_alignment_annotation_path, "\n")
    
    # Copy files if the paths are present
    if (!is.na(aa_alignment_path) && file.exists(aa_alignment_path)) {
      file.copy(aa_alignment_path, file.path(translation_folder, basename(aa_alignment_path)))
    }
    if (!is.na(aa_alignment_annotation_path) && file.exists(aa_alignment_annotation_path)) {
      file.copy(aa_alignment_annotation_path, file.path(translation_folder, basename(aa_alignment_annotation_path)))
    }
    if (!is.na(na_alignment_path) && file.exists(na_alignment_path)) {
      file.copy(na_alignment_path, file.path(nucleic_acid_folder, basename(na_alignment_path)))
    }
    if (!is.na(na_alignment_annotation_path) && file.exists(na_alignment_annotation_path)) {
      file.copy(na_alignment_annotation_path, file.path(nucleic_acid_folder, basename(na_alignment_annotation_path)))
    }
  }
  
  message("Folders and files created successfully under www/data/")
}


library(dplyr)
library(purrr)
library(readr)
library(seqinr)
library(tidyr)
library(stringr)
# MAIN SCRIPT
# Get the current working directory
current_directory <- getwd()
print(current_directory)
setwd(current_directory)

data_table_file <- "gene_data_table.csv"
annotation_paths_table_file <- "gene_annotation_paths.csv"

# Read CSV files into dataframe
gene_data_df <- read.csv(file.path(current_directory, "www", "data", "resources", data_table_file), header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
annotation_paths_df <- read.csv(file.path(current_directory, "www", "data","resources", annotation_paths_table_file), header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# Check the first few rows of the dataframe
head(gene_data_df)
head(annotation_paths_df)

# Call the function
create_folders_by_gene("data", dataframe=gene_data_df, 1:length(gene_data_df$`Gene ID`), annotations_df=annotation_paths_df)
