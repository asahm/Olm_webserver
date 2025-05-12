# script for visualizing one gene alignment only
# alignment for both amino acid and nucleic acid

library(cowplot)
library(grid)
library(seqinr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(gridExtra)
library(stringr)
library(plotly)
library(ggtext)
library(readr)
library(jsonlite)

create_annotation_information_df <- function(file_contents){
  
  # Split the text by newline
  lines <- strsplit(file_contents, "\n")[[1]]
  # Initialize variables for storing results
  positive_selection_probability <- NULL
  analyzed_positions <- NULL
  
  # Loop through the lines to find the relevant information
  for (i in 1:(length(lines) - 1)) {
    line <- lines[i]
    
    # Match "LINE_GRAPH" and "Positive selection probability" in the same line
    if (grepl("^LINE_GRAPH\\s.*Positive selection probability", line)) {
      positive_selection_probability <- line
    }
    
    # Match "NO_GRAPH" and "Analyzed positions" in the same line
    if (grepl("^NO_GRAPH\\s.*Analyzed positions", line)) {
      analyzed_positions <- line
    }
  }
  
  # Create a data frame with the extracted data
  result <- data.frame(
    positive_selection_probability = positive_selection_probability,
    analyzed_positions = analyzed_positions,
    stringsAsFactors = FALSE
  )
  #print("annotation_info_df")
  #print(head(result))
  return(result)
}

remove_first_pipe_after_E <- function(input_string) {
  # Remove the "NO_GRAPH\tAnalyzed positions\t" part
  without_prefix <- sub("^NO_GRAPH\\tAnalyzed positions\\t", "", input_string)
  # Use gsub with a regular expression to remove the first | after each E
  result <- gsub("E\\|", "E", without_prefix)
  return(result)
}

get_positive_selection_percentages <- function(input_string) {
  #print(input_string)
  
  # Remove the "LINE_GRAPH\tPositive selection probability\t" part
  without_prefix <- sub("LINE_GRAPH\\tPositive selection probability\\t", "", input_string)
  #print(without_prefix)
  
  # Step 1: Replace patterns like "0.229,.,0.23%" with "0.23"
  cleaned <- gsub("(\\d+\\.\\d+),\\.,(\\d+\\.\\d+)%", "\\2", without_prefix)
  #print(cleaned)
  
  # Step 2: Replace patterns like "90.737,*,90.74%" with "90.74"
  cleaned <- gsub("(\\d+\\.\\d+),\\*,(\\d+\\.\\d+)%", "\\2", cleaned)
  #print(cleaned)
  
  # Step 4: Extract numeric values and '|' characters into a list
  numeric_values <- unlist(strsplit(cleaned, "\\|"))
  #print(numeric_values)
  
  return(numeric_values)
}

# Function to ensure the last 3 rows are fixed
# Improved function to ensure the last 3 rows are fixed and in correct order
ensure_fixed_rows <- function(df) {
  # Find the Proteus row dynamically
  proteus_row <- rownames(df)[grepl("^Proteus anguinus", rownames(df))]
  
  if(length(proteus_row) == 0) {
    warning("No row starting with 'Proteus anguinus' found in the dataframe")
    # If no Proteus row, create a placeholder
    proteus_row <- "Proteus anguinus"
  }
  
  # Define the fixed row patterns with the actual Proteus row name
  fixed_rows <- c(proteus_row, "Analyzed Sites", "Positive Selection Probability")
  
  # Create missing fixed rows with empty values if not present
  for (row_name in fixed_rows) {
    if (!row_name %in% rownames(df)) {
      new_row <- setNames(
        data.frame(matrix(NA, nrow = 1, ncol = ncol(df))), 
        colnames(df)
      )
      # Set species_name for the new row
      new_row$species_name <- row_name
      
      # Bind to the original dataframe
      df <- rbind(df, new_row)
      rownames(df)[nrow(df)] <- row_name
    }
  }
  
  # Ensure the fixed rows are in the correct order at the end
  non_fixed_rows <- setdiff(rownames(df), fixed_rows)
  df <- df[c(non_fixed_rows, fixed_rows), , drop = FALSE]
  
  return(df)
}


generate_aa_alignment_df <- function(current_directory, gene_id){
  
  cat(gene_id, "\n")
  # STEP 1: PREPROCESSING - read in the data  
  aa_directory_path <- file.path(current_directory, "data", "genes",
                                 gene_id, "translation",
                                 "prank.best.fas.translation")
  cat(aa_directory_path, "\n")
  
  # Check if the file exists
  if (!file.exists(aa_directory_path)) {
    # Display a message if the file does not exist
    message("The file does not exist: ", aa_directory_path, "\n")
  } else {
    message("The file exists: ", aa_directory_path, "\n")
    
    # Load the alignment using msa package
    aas <- Biostrings::readAAStringSet(aa_directory_path, format = "fasta")
    
    # Get sequence names and lengths
    seq_names <- names(aas)
    seq_lengths <- sapply(aas, length)
    
    # Convert dataset to to matrix
    aa_alignment_matrix <- as.matrix(aas)
    
    rownames(aa_alignment_matrix) <- gsub("_", " ", rownames(aa_alignment_matrix))
    
    # Create the dataframe
    aa_alignment_df <- data.frame(
      species_name = rownames(aa_alignment_matrix),
      aa_alignment_matrix
    )
    
    # Read the annotation file content
    aa_annotations_directory <- file.path(current_directory,
                                         "data","genes",
                                          gene_id, "translation", 
                                          "prank.best.fas.translation.view.annotations")
    file_content <- readr::read_file(aa_annotations_directory)[1]
    aa_annotation_df <- create_annotation_information_df(file_content)
    
    # positive selection indicator
    # Extract and print the NO_GRAPH content along with its length
    aa_pos_selection_annot <- aa_annotation_df$analyzed_positions[1]
    
    # clean up the no_graph line
    aa_pos_selection_clean_info <- remove_first_pipe_after_E(aa_pos_selection_annot)
    aa_pos_selection_clean_info_vector <- strsplit(aa_pos_selection_clean_info, "")[[1]]
    
    # Add title to the beginning of vector
    aa_pos_selection_clean_info_vec <- c("Analyzed Sites", aa_pos_selection_clean_info_vector)
    
    # Calculate the difference in length
    length_diff <- ncol(aa_alignment_df) - length(aa_pos_selection_clean_info_vec)
    
    # Extend the vector with "|" if it's shorter
    if (length_diff > 0) {
      #print("analyzed positions vector is shorter")
      aa_pos_selection_clean_info_vec <- c(aa_pos_selection_clean_info_vec, rep("|", length_diff))
    }
    
    # Add the new row to the DataFrame using rbind()
    aa_combined_alignment_df <- rbind(aa_alignment_df, aa_pos_selection_clean_info_vec)
    
    # Identify the row index for "_positive_selection_indication"
    aa_pos_selection_row <- which(aa_combined_alignment_df$species_name == "Analyzed Sites")
    
    # Replace 'E' with '*' only in this specific row
    if (length(aa_pos_selection_row) > 0) {
      for (col in names(aa_combined_alignment_df)[names(aa_combined_alignment_df) != "species_name"]) {
        aa_combined_alignment_df[aa_pos_selection_row, col] <- ifelse(
          aa_combined_alignment_df[aa_pos_selection_row, col] == 'E', 
          '*', 
          aa_combined_alignment_df[aa_pos_selection_row, col]
        )
      }
    }
    
    # positive selection percentages
    
    aa_pos_selection_percentage_annot <- aa_annotation_df$positive_selection_probability[1]
    #print(aa_pos_selection_percentage_annot)
    
    # extract the percentage annotations
    aa_pos_selection_percentage_annot_vector <- get_positive_selection_percentages(aa_pos_selection_percentage_annot)
    aa_pos_selection_percentage_annot_vector <- c("Positive Selection Probability", aa_pos_selection_percentage_annot_vector)
    
    # Calculate the difference in length
    length_diff <- ncol(aa_combined_alignment_df) - length(aa_pos_selection_percentage_annot_vector)
    
    # Extend the vector with "|" if it's shorter
    if (length_diff > 0) {
      aa_pos_selection_percentage_annot_vector <- c(aa_pos_selection_percentage_annot_vector, rep("0.0", length_diff))
    }
    
    # add the percentage information to the df
    aa_combined_alignment_df <- rbind(aa_combined_alignment_df, aa_pos_selection_percentage_annot_vector)
    
    # fix the last 3 rows
    aa_combined_alignment_df <- ensure_fixed_rows(aa_combined_alignment_df)
    
    # Identify the row that starts with "Proteus anguinus"
    proteus_row <- grep("^Proteus anguinus", rownames(aa_combined_alignment_df), value = TRUE)
    
    # Update fixed_rows dynamically
    fixed_rows <- c(proteus_row, "Analyzed Sites", "Positive Selection Probability")
    
    # Get all row names in the desired order
    all_levels <- c(
      setdiff(rownames(aa_combined_alignment_df), fixed_rows),  # Non-fixed rows
      fixed_rows  # Fixed rows at the end
    )
    
    # Create a named vector for text formatting
    text_format <- setNames(
      rep("plain", length(all_levels)),
      all_levels
    )
    if (length(proteus_row) > 0) {
      text_format[proteus_row] <- "bold"
    }
    
    # Convert species_name to factor with these levels
    aa_combined_alignment_df$species_name <- factor(
      aa_combined_alignment_df$species_name,
      levels = rev(all_levels)  # Reverse to get the desired visual order
    )
    
    #output directory
    aa_output_dir = file.path(current_directory, "data", "genes",
                              gene_id, "translation")
    output_file <- file.path(aa_output_dir, "amino_acid_alignment.csv")
    
    aa_combined_alignment_df <- aa_combined_alignment_df[!rownames(aa_combined_alignment_df) %in% c("Analyzed Sites", "Positive Selection Probability"), ]
    
    # Save the dataframe as a CSV
    write.csv(aa_combined_alignment_df, file = output_file, row.names = FALSE)
    
    cat("AA alignment created and saved.", "\n")
  }
}


generate_alignment_csv <- function(
    current_directory, 
    gene_id
){
  
  #generate aa alignment
  generate_aa_alignment_df(current_directory = current_directory, 
                                         gene_id = gene_id)
  cat("AA alignment csv is created")
}

# MAIN SCRIPT

main <- function(){
  
  # Get the current working directory
  current_directory <- getwd()
  print(current_directory)
  setwd(current_directory)
  
  gene_id = "TRINITY_DN100027_c0_g2"
  
  generate_alignment_visuals(current_directory = current_directory,
                             gene_id = gene_id)
  
} 

main()

