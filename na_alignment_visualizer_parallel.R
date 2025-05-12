# Modified script for visualizing multiple gene alignments in parallel
# Handles both amino acid and nucleic acid alignments

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
library(parallel)
library(foreach)
library(doParallel)


plot_and_save_na_alignment <- function(
    alignment_df, 
    output_file = "alignment", 
    seq_length = NULL, 
    chunk_size = 120,
    export_format = c("pdf", "png"),
    output_dir,
    max_plot_height = 20,
    max_plot_width = 20
) {
  if (nrow(alignment_df) == 0) {
    stop("Alignment dataframe is empty")
  }
  
  if (is.null(seq_length)) {
    seq_length <- ncol(alignment_df) - 1
  }
  
  # Find Proteus anguinus row
  proteus_row <- grep("Proteus anguinus", alignment_df$species_name, value = TRUE)
  
  # Reorder the dataframe with the specified order for the last three rows
  other_species <- setdiff(
    unique(as.character(alignment_df$species_name)), 
    c(proteus_row, "Analyzed Sites", "Positive Selection Probability")
  )
  
  # Create ordered levels - this is crucial for consistent ordering
  ordered_levels <- c(
    other_species,
    proteus_row, 
    "Analyzed Sites", 
    "Positive Selection Probability"
  )
  
  # Set factor levels to ensure consistent order
  alignment_df$species_name <- factor(
    alignment_df$species_name,
    levels = ordered_levels
  )
  
  # The critical change: ensure alignment_df is actually sorted by the factor levels
  alignment_df <- alignment_df[order(alignment_df$species_name), ]
  
  #print(head(alignment_df, n=10))
  
  chunks <- seq(1, seq_length, by = chunk_size)
  num_chunks <- length(chunks)
  
  # Calculate the max height according to length of alignment
  if (seq_length <= 900) {
    max_plot_height = 10
  } else {
    max_plot_height = 30
  }
  
  # Color palette for annotations with Hex Codes
  aa_colors <- c(
    # Gap
    '-' = '#FFFFFF',   # White for gaps
    # Special symbol for * (navy blue background)
    '*' = '#000080'    # Navy blue
  )
  
  # Vectorized function to convert percentage strings to numeric values
  extract_percentage <- function(x) {
    as.numeric(gsub("%", "", x))
  }
  
  # Vectorized function to create white-to-dark-red gradient based on percentage
  get_red_gradient <- Vectorize(function(percentage) {
    if (is.na(percentage)) return("#FFFFFF")
    red_value <- 255
    green_value <- as.integer(255 * (1 - percentage/100))
    blue_value <- as.integer(255 * (1 - percentage/100))
    green_value <- max(0, min(255, green_value))
    blue_value <- max(0, min(255, blue_value))
    sprintf("#%02X%02X%02X", red_value, green_value, blue_value)
  })
  
  plot_list <- lapply(seq_along(chunks), function(chunk_idx) {
    start_index <- chunks[chunk_idx]
    end_index <- min(start_index + chunk_size - 1, seq_length)
    
    # Calculate the x-axis limits for consistent scaling
    x_min <- start_index - 0.5  # Adjust minimum to include first amino acid
    x_max <- start_index + chunk_size - 0.5  # Adjust maximum to include last amino acid
    
    plot_chunk <- alignment_df %>%
      mutate(across(
        .cols = -species_name, 
        .fns = as.character
      )) %>%
      select(species_name, all_of(paste0("X", start_index:end_index))) %>%
      pivot_longer(
        cols = -species_name, 
        names_to = "position", 
        values_to = "amino_acid"
      ) %>%
      mutate(
        position = as.numeric(sub("X", "", position)),
        amino_acid = ifelse(amino_acid %in% c('|'), NA_character_, amino_acid)
      )
    
    # CRITICAL: Make sure species_name factor levels are preserved
    plot_chunk$species_name <- factor(plot_chunk$species_name, levels = ordered_levels)
    
    # Handle colors and text display
    plot_chunk <- plot_chunk %>%
      mutate(
        percentage = case_when(
          species_name == "Positive Selection Probability" ~ extract_percentage(amino_acid),
          TRUE ~ NA_real_
        ),
        color = case_when(
          species_name == "Positive Selection Probability" ~ get_red_gradient(percentage),
          species_name == "Analyzed Sites" & amino_acid == '*' ~ '#000080',
          species_name == "Analyzed Sites" ~ '#FFFFFF',
          amino_acid == '-' ~ '#FFFFFF',
          amino_acid == '*' ~ '#000080',
          !is.na(amino_acid) ~ aa_colors[amino_acid],
          TRUE ~ 'white'
        ),
        text_color = case_when(
          grepl("_percentages", as.character(species_name)) & grepl("\\*", amino_acid) ~ 'black',
          grepl("_percentages", as.character(species_name)) & amino_acid == '*' ~ '#000080',
          grepl("_percentages", as.character(species_name)) ~ 'white',
          amino_acid == '*' ~ '#000080',
          TRUE ~ 'black'
        ),
        display_text = ifelse(species_name == "Positive Selection Probability", "", amino_acid)
      ) %>%
      filter(!is.na(amino_acid))
    
    # Create a lookup vector for face styles
    face_styles <- rep("plain", length(ordered_levels))
    names(face_styles) <- ordered_levels
    face_styles[grep("Proteus anguinus", names(face_styles))] <- "bold"
    
    # Create custom breaks for x-axis (every 25 positions)
    x_breaks <- seq(
      from = 0, #ceiling(x_min/20) * 20,  # Start from next multiple of 25
      to = x_max, #floor(x_max/25) * 25,      # End at last multiple of 25
      by = 20
    )
    
    # Create a lookup vector for face styles - more direct approach
    face_styles <- setNames(
      ifelse(ordered_levels == proteus_row, "bold", "plain"),
      ordered_levels
    )
    
    # Create the plot with fixed width elements
    p <- ggplot(plot_chunk, aes(x = position, y = species_name)) +
      geom_tile(
        aes(fill = I(color)), 
        color = "white", 
        width = 1, 
        height = 0.8
      ) +
      geom_text(
        aes(label = display_text, color = I(text_color)), 
        size = 2
      ) +
      scale_x_continuous(
        limits = c(x_min, x_max),
        breaks = x_breaks,  # Note: x_breaks seems undefined in your code
        expand = c(0, 0)
      ) +
      # The key fix: use rev(ordered_levels) to get bottom-to-top ordering
      scale_y_discrete(limits = rev(ordered_levels)) +
      labs(x = "", y = "") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(size = 6, angle = 0, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(5, 10, 5, 10),
        axis.text.y = element_text(
          size = 6,
          face = face_styles[as.character(rev(ordered_levels))]
        )
      )
    
    return(p)
  })
  
  # Rest of your function remains the same...
  # Create a custom heatmap legend
  legend_data <- data.frame(
    percentage = seq(0, 100, by = 10)
  )
  
  legend_plot <- ggplot(legend_data, aes(x = 1, y = percentage, fill = get_red_gradient(percentage))) +
    geom_tile() +
    scale_fill_identity() +
    labs(x = "", y = "Positive Selection Probability") +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(5, 5, 5, 5)
    ) +
    geom_text(aes(label = paste0(percentage, "%")), color = "black", size = 1.5, hjust = 0.5, vjust = -1)
  
  # Combine the plot and legend side by side using cowplot
  combined_plot <- plot_grid(
    do.call(grid.arrange, c(plot_list, ncol = 1)),
    legend_plot,
    ncol = 2,
    rel_widths = c(2, 0.2)  # Set the legend width to be thin
  )
  
  # Save plot
  export_format <- match.arg(export_format)
  plot_height <- max_plot_height
  plot_width <- min(max_plot_width, 10)
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  full_file_path <- file.path(output_dir, paste0(output_file, ".", export_format))
  
  ggsave(
    full_file_path, 
    plot = combined_plot, 
    width = plot_width, 
    height = plot_height, 
    units = "in", 
    dpi = 600,
    limitsize = FALSE
  )
  
  return(combined_plot)
}

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

generate_and_save_na_alignment_visuals <- function(current_directory, gene_id){
  tryCatch({
    cat(paste0("Processing gene: ", gene_id, "\n"))
    # STEP 1: PREPROCESSING - read in the data  
    na_directory_path <- file.path(current_directory, "www", "data", 
                                   gene_id, "nucleic_acid",
                                   "prank.best.fas.renamed")
    
    # Check if the file exists
    if (!file.exists(na_directory_path)) {
      message(paste0("The file does not exist: ", na_directory_path))
      return(FALSE)  # Return FALSE to indicate failure
    }
    
    message(paste0("The file exists: ", na_directory_path))
    
    # Load the alignment using Biostrings package
    aas <- Biostrings::readAAStringSet(na_directory_path, format = "fasta")
    
    # Get sequence names and lengths
    seq_names <- names(aas)
    seq_lengths <- sapply(aas, length)
    
    # Convert dataset to matrix
    na_alignment_matrix <- as.matrix(aas)
    
    rownames(na_alignment_matrix) <- gsub("_", " ", rownames(na_alignment_matrix))
    
    # Create the dataframe
    na_alignment_df <- data.frame(
      species_name = rownames(na_alignment_matrix),
      na_alignment_matrix
    )
    
    # Read the annotation file content
    na_annotations_directory <- file.path(current_directory,
                                          "www", "data",
                                          gene_id, "nucleic_acid", 
                                          "prank.best.fas.view.annotations")
    
    if (!file.exists(na_annotations_directory)) {
      message(paste0("Annotation file does not exist: ", na_annotations_directory))
      return(FALSE)
    }
    
    # Read file content safely
    file_content <- tryCatch({
      readr::read_file(na_annotations_directory)[1]
    }, error = function(e) {
      message(paste0("Error reading annotation file: ", e$message))
      return(NULL)
    })
    
    if (is.null(file_content)) {
      return(FALSE)
    }
    
    na_annotation_df <- create_annotation_information_df(file_content)
    
    # positive selection indicator
    # Extract and print the NO_GRAPH content along with its length
    na_pos_selection_annot <- na_annotation_df$analyzed_positions[1]
    
    # Check if annotation is valid
    if (is.null(na_pos_selection_annot) || is.na(na_pos_selection_annot)) {
      message("Invalid or missing analyzed positions annotation")
      return(FALSE)
    }
    
    # clean up the no_graph line
    na_pos_selection_clean_info <- remove_first_pipe_after_E(na_pos_selection_annot)
    na_pos_selection_clean_info_vector <- strsplit(na_pos_selection_clean_info, "")[[1]]
    
    # Add title to the beginning of vector
    na_pos_selection_clean_info_vec <- c("Analyzed Sites", na_pos_selection_clean_info_vector)
    
    # Calculate the difference in length
    length_diff <- ncol(na_alignment_df) - length(na_pos_selection_clean_info_vec)
    
    # Extend the vector with "|" if it's shorter
    if (length_diff > 0) {
      na_pos_selection_clean_info_vec <- c(na_pos_selection_clean_info_vec, rep("|", length_diff))
    } else if (length_diff < 0) {
      # If the vector is too long, truncate it
      na_pos_selection_clean_info_vec <- na_pos_selection_clean_info_vec[1:ncol(na_alignment_df)]
    }
    
    # Add the new row to the DataFrame using rbind()
    na_combined_alignment_df <- rbind(na_alignment_df, na_pos_selection_clean_info_vec)
    
    # Identify the row index for "_positive_selection_indication"
    na_pos_selection_row <- which(na_combined_alignment_df$species_name == "Analyzed Sites")
    
    # Replace 'E' with '*' only in this specific row
    if (length(na_pos_selection_row) > 0) {
      for (col in names(na_combined_alignment_df)[names(na_combined_alignment_df) != "species_name"]) {
        na_combined_alignment_df[na_pos_selection_row, col] <- ifelse(
          na_combined_alignment_df[na_pos_selection_row, col] == 'E', 
          '*', 
          na_combined_alignment_df[na_pos_selection_row, col]
        )
      }
    }
    
    # positive selection percentages
    na_pos_selection_percentage_annot <- na_annotation_df$positive_selection_probability[1]
    
    # Check if annotation is valid
    if (is.null(na_pos_selection_percentage_annot) || is.na(na_pos_selection_percentage_annot)) {
      message("Invalid or missing positive selection probability annotation")
      return(FALSE)
    }
    
    # extract the percentage annotations
    na_pos_selection_percentage_annot_vector <- get_positive_selection_percentages(na_pos_selection_percentage_annot)
    na_pos_selection_percentage_annot_vector <- c("Positive Selection Probability", na_pos_selection_percentage_annot_vector)
    
    # Calculate the difference in length
    length_diff <- ncol(na_combined_alignment_df) - length(na_pos_selection_percentage_annot_vector)
    
    # Extend the vector with "|" if it's shorter
    if (length_diff > 0) {
      na_pos_selection_percentage_annot_vector <- c(na_pos_selection_percentage_annot_vector, rep("0.0", length_diff))
    } else if (length_diff < 0) {
      # If the vector is too long, truncate it
      na_pos_selection_percentage_annot_vector <- na_pos_selection_percentage_annot_vector[1:ncol(na_combined_alignment_df)]
    }
    
    # add the percentage information to the df
    na_combined_alignment_df <- rbind(na_combined_alignment_df, na_pos_selection_percentage_annot_vector)
    
    # fix the last 3 rows
    na_combined_alignment_df <- ensure_fixed_rows(na_combined_alignment_df)
    
    # Identify the row that starts with "Proteus anguinus"
    proteus_row <- grep("^Proteus anguinus", rownames(na_combined_alignment_df), value = TRUE)
    
    # Update fixed_rows dynamically
    fixed_rows <- c(proteus_row, "Analyzed Sites", "Positive Selection Probability")
    
    # Get all row names in the desired order
    all_levels <- c(
      setdiff(rownames(na_combined_alignment_df), fixed_rows),  # Non-fixed rows
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
    na_combined_alignment_df$species_name <- factor(
      na_combined_alignment_df$species_name,
      levels = rev(all_levels)  # Reverse to get the desired visual order
    )
    
    # Output directory
    na_output_dir = file.path(current_directory, "www", "data", 
                              gene_id, "nucleic_acid")
    
    # Create output directory if it doesn't exist
    if (!dir.exists(na_output_dir)) {
      dir.create(na_output_dir, recursive = TRUE)
    }
    
    output_file <- file.path(na_output_dir, "amino_acid_alignment.csv")
    
    # Save the dataframe as a CSV
    write.csv(na_combined_alignment_df, file = output_file, row.names = FALSE)
    
    # Generate and save the plot
    tryCatch({
      plot_and_save_na_alignment(
        alignment_df = na_combined_alignment_df,
        output_dir = na_output_dir,
        output_file = "alignment", 
        export_format = "png"
      )
      
      cat(paste0("NA alignment created and saved for gene: ", gene_id, "\n"))
      return(TRUE)  # Return TRUE to indicate success
    }, error = function(e) {
      message(paste0("Error creating plot for gene: ", gene_id, ": ", e$message))
      return(FALSE)
    })
    
  }, error = function(e) {
    message(paste0("Error processing gene: ", gene_id, ": ", e$message))
    return(FALSE)  # Return FALSE to indicate failure
  })
}

# function to process a list of gene IDs in parallel with batch tracking
generate_alignment_visuals_parallel <- function(
    current_directory, 
    gene_ids,
    batch_size = 500,  # Process genes in batches of this size
    num_cores = NULL   # Number of cores to use (default: detect automatically)
){
  # Set up parallel backend
  if (is.null(num_cores)) {
    # Use 75% of available cores, but at least 1 and at most 8
    num_cores <- max(1, min(8, parallel::detectCores() * 0.75))
  }
  
  cat(paste0("Setting up parallel processing with ", num_cores, " cores\n"))
  
  # Create summary tracking variables
  total_genes <- length(gene_ids)
  all_successful_genes <- character(0)
  all_failed_genes <- character(0)
  all_error_messages <- list()
  
  # Split genes into batches
  num_batches <- ceiling(total_genes / batch_size)
  batches <- split(gene_ids, ceiling(seq_along(gene_ids) / batch_size))
  
  cat(paste0("Processing ", total_genes, " genes in ", num_batches, " batches of up to ", 
             batch_size, " genes each\n\n"))
  
  # Process each batch
  for (batch_idx in seq_along(batches)) {
    current_batch <- batches[[batch_idx]]
    batch_start_time <- Sys.time()
    
    cat(paste0("---- BATCH ", batch_idx, "/", num_batches, " ----\n"))
    cat(paste0("Starting batch at: ", format(batch_start_time), "\n"))
    cat(paste0("Processing ", length(current_batch), " genes in this batch\n"))
    
    # Set up cluster for this batch
    cl <- makeCluster(num_cores)
    registerDoParallel(cl)
    
    # Export necessary functions to the worker nodes
    clusterExport(cl, c(
      "create_annotation_information_df", 
      "remove_first_pipe_after_E", 
      "get_positive_selection_percentages", 
      "ensure_fixed_rows", 
      "plot_and_save_na_alignment",
      "generate_and_save_na_alignment_visuals"
    ))
    
    # Process batch in parallel with error handling
    results <- foreach(
      gene_id = current_batch,
      .packages = c("dplyr", "ggplot2", "tidyr", "readr", "Biostrings", "cowplot", "grid", "gridExtra", "stringr"),
      .errorhandling = "pass"
    ) %dopar% {
      # Wrap in tryCatch to ensure we always return a valid result structure
      tryCatch({
        # Process the gene
        success <- generate_and_save_na_alignment_visuals(current_directory, gene_id)
        # Return a structured result including the gene_id
        return(list(gene_id = gene_id, success = success, error = NULL))
      }, error = function(e) {
        # If an error occurs, capture it and return a structured result
        return(list(gene_id = gene_id, success = FALSE, error = e$message))
      })
    }
    
    # Stop the cluster for this batch
    stopCluster(cl)
    
    # Process batch results
    batch_successful_genes <- character(0)
    batch_failed_genes <- character(0)
    batch_error_messages <- list()
    
    for (i in seq_along(results)) {
      result <- results[[i]]
      # Check if the result is valid
      if (is.list(result) && !is.null(result$gene_id)) {
        gene_id <- result$gene_id
        if (isTRUE(result$success)) {
          batch_successful_genes <- c(batch_successful_genes, gene_id)
        } else {
          batch_failed_genes <- c(batch_failed_genes, gene_id)
          if (!is.null(result$error)) {
            batch_error_messages[[gene_id]] <- result$error
          }
        }
      } else {
        # Handle malformed result - record the index
        batch_failed_genes <- c(batch_failed_genes, paste0("unknown_batch", batch_idx, "_", i))
      }
    }
    
    # Add batch results to overall results
    all_successful_genes <- c(all_successful_genes, batch_successful_genes)
    all_failed_genes <- c(all_failed_genes, batch_failed_genes)
    all_error_messages <- c(all_error_messages, batch_error_messages)
    
    # Calculate batch processing time
    batch_end_time <- Sys.time()
    batch_duration <- difftime(batch_end_time, batch_start_time, units = "mins")
    
    # Print batch summary
    cat("\n---- BATCH SUMMARY ----\n")
    cat(paste0("Batch ", batch_idx, "/", num_batches, " completed in ", 
               round(as.numeric(batch_duration), 2), " minutes\n"))
    cat(paste0("Successfully processed: ", length(batch_successful_genes), "/", 
               length(current_batch), " genes\n"))
    cat(paste0("Failed to process: ", length(batch_failed_genes), "/", 
               length(current_batch), " genes\n"))
    
    # Save batch results to file
    batch_summary_file <- file.path(current_directory, 
                                    paste0("alignment_batch_", batch_idx, "_summary.txt"))
    sink(batch_summary_file)
    cat(paste0("---- ALIGNMENT BATCH ", batch_idx, "/", num_batches, " SUMMARY ----\n"))
    cat(paste0("Date: ", format(batch_end_time), "\n"))
    cat(paste0("Batch duration: ", round(as.numeric(batch_duration), 2), " minutes\n\n"))
    cat(paste0("Total genes in batch: ", length(current_batch), "\n"))
    cat(paste0("Successfully processed: ", length(batch_successful_genes), "\n"))
    cat(paste0("Failed to process: ", length(batch_failed_genes), "\n\n"))
    
    if (length(batch_successful_genes) > 0) {
      cat("Successfully processed genes in this batch:\n")
      cat(paste(batch_successful_genes, collapse = "\n"), "\n\n")
    }
    
    if (length(batch_failed_genes) > 0) {
      cat("Failed genes in this batch:\n")
      cat(paste(batch_failed_genes, collapse = "\n"), "\n\n")
      
      # Include error messages for this batch
      if (length(batch_error_messages) > 0) {
        cat("Error details:\n")
        for (gene_id in names(batch_error_messages)) {
          cat(paste0(gene_id, ": ", batch_error_messages[[gene_id]], "\n"))
        }
      }
    }
    sink()
    
    cat(paste0("\nBatch summary saved to: ", batch_summary_file, "\n\n"))
    
    # Estimate remaining time
    if (batch_idx < num_batches) {
      avg_time_per_batch <- as.numeric(batch_duration)
      est_remaining_time <- avg_time_per_batch * (num_batches - batch_idx)
      cat(paste0("Estimated remaining time: ~", round(est_remaining_time, 1), 
                 " minutes (", round(est_remaining_time/60, 1), " hours)\n\n"))
    }
  }
  
  # Create and save final summary
  overall_summary <- list(
    total_genes = total_genes,
    successful_genes = all_successful_genes,
    failed_genes = all_failed_genes,
    error_messages = all_error_messages,
    num_successful = length(all_successful_genes),
    num_failed = length(all_failed_genes),
    num_batches = num_batches
  )
  
  # Print final summary
  cat("\n---- FINAL PROCESSING SUMMARY ----\n")
  cat(paste0("Total genes processed: ", total_genes, "\n"))
  cat(paste0("Successfully processed: ", length(all_successful_genes), "\n"))
  cat(paste0("Failed to process: ", length(all_failed_genes), "\n"))
  
  if (length(all_failed_genes) > 0) {
    cat("\nFailed genes (first 10 shown):\n")
    cat(paste(head(all_failed_genes, 10), collapse = "\n"), "\n")
    if (length(all_failed_genes) > 10) {
      cat(paste0("... and ", length(all_failed_genes) - 10, " more\n"))
    }
  }
  
  # Save detailed final summary to file
  summary_file <- file.path(current_directory, "alignment_processing_summary.txt")
  sink(summary_file)
  cat("---- ALIGNMENT PROCESSING FINAL SUMMARY ----\n")
  cat(paste0("Date: ", Sys.time(), "\n\n"))
  cat(paste0("Total genes processed: ", total_genes, "\n"))
  cat(paste0("Total batches: ", num_batches, "\n"))
  cat(paste0("Successfully processed: ", length(all_successful_genes), "\n"))
  cat(paste0("Failed to process: ", length(all_failed_genes), "\n\n"))
  
  if (length(all_successful_genes) > 0) {
    cat("Successfully processed genes:\n")
    cat(paste(all_successful_genes, collapse = "\n"), "\n\n")
  }
  
  if (length(all_failed_genes) > 0) {
    cat("Failed genes:\n")
    cat(paste(all_failed_genes, collapse = "\n"), "\n\n")
    
    # Include error messages
    if (length(all_error_messages) > 0) {
      cat("Error details (first 30 shown):\n")
      count <- 0
      for (gene_id in names(all_error_messages)) {
        if (count >= 30) break
        cat(paste0(gene_id, ": ", all_error_messages[[gene_id]], "\n"))
        count <- count + 1
      }
    }
  }
  sink()
  
  cat(paste0("\nFinal summary saved to: ", summary_file, "\n"))
  
  return(overall_summary)
}

# Main function
main <- function(){
  # Get the current working directory
  current_directory <- getwd()
  cat(paste0("Current directory: ", current_directory, "\n"))
  
  # Or use a list of gene IDs directly
  data_table_file <- "gene_data_table.csv"
  gene_data_df <- read.csv(file.path(current_directory, "www", "data", "resources", data_table_file), header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  print(head(gene_data_df))
  #gene_ids <- gene_data_df$`Gene ID`
  gene_ids <- c('TRINITY_DN12397_c0_g1')
  
  # Process all genes in parallel
  summary <- generate_alignment_visuals_parallel(
    current_directory = current_directory,
    gene_ids = gene_ids,
    #num_cores = 4  # Adjust based on your machine's capabilities
  )
  
  cat("All alignment files have been processed\n")
}

# Run the main function
main()