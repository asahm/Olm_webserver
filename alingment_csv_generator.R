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


plot_and_save_aa_alignment <- function(
    alignment_df, 
    output_file = "alignment", 
    seq_length = NULL, 
    chunk_size = 75,
    export_format = c("pdf", "png", "svg"),
    output_dir,
    max_plot_height = 5,
    max_plot_width = 20,
    tile_aspect_ratio = 1,  # Parameter to control tile aspect ratio
    base_font_size = 2,      # parameter for base font size
    species_font_size = 5   # font size for species names
) {
  if (nrow(alignment_df) == 0) {
    stop("Alignment dataframe is empty")
  }
  
  if (is.null(seq_length)) {
    seq_length <- ncol(alignment_df) - 1
  }
  
  if (nrow(alignment_df)-3 > 5) {
    base_font_size = 1.5     # parameter for base font size
    species_font_size = 4
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
  
  # Calculate optimal chunk size based on sequence length
  if (seq_length > 500) {
    chunk_size <- min(chunk_size, 60)  # Reduce chunk size for long sequences
  }
  
  chunks <- seq(1, seq_length, by = chunk_size)
  num_chunks <- length(chunks)
  
  # Calculate the max height according to length of alignment
  num_species <- length(ordered_levels)
  
  if (seq_length <= 350) {
    max_plot_height <- min(5, num_species * 0.3)
  } else {
    max_plot_height <- min(20, num_species * 0.3)
  }
  
  # Adjust plot width based on chunk size
  chunk_plot_width <- min(max_plot_width, chunk_size * 0.12)
  
  # Calculate font scaling factor based on chunk size and tile aspect ratio
  # This will scale the font size proportionally with the tile size
  font_scaling_factor <- min(1, 90 / chunk_size) * sqrt(tile_aspect_ratio)  ## TODO CHECK THIS RATIO!!! - MIGHT BE TOO SMALL/BIG
  scaled_font_size <- base_font_size * font_scaling_factor
  scaled_font_size_y <- species_font_size * font_scaling_factor
  
  # Color palette for amino acids
  aa_colors <- c(
    'A' = '#32CD32', 'R' = '#696969', 'N' = '#FF6347',
    'D' = '#FFB6C1', 'C' = '#FFFFE0', 'Q' = '#DC143C',
    'G' = '#4B0082', 'E' = '#FF0000', 'H' = '#708090',
    'I' = '#228B22', 'L' = '#808080', 'K' = '#00709B',
    'M' = '#F4A460', 'F' = '#696969', 'P' = '#333333',
    'S' = '#FF8C00', 'T' = '#A52A2A', 'W' = '#98FB98',
    'Y' = '#F5DEB3', 'V' = '#FFB0CB', '-' = '#FFFFFF',
    '*' = '#000080'
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
  
  # Pre-calculate global xlim for ALL chunks to ensure consistent scaling
  global_xmin <- min(chunks) - 0.5
  global_xmax <- min(seq_length + 0.5, max(chunks) + chunk_size - 0.5)
  
  plot_list <- lapply(seq_along(chunks), function(chunk_idx) {
    start_index <- chunks[chunk_idx]
    end_index <- min(start_index + chunk_size - 1, seq_length)
    
    # Calculate the x-axis limits for this specific chunk
    # Important: limits will be different for each chunk but scales will be consistent
    x_min <- start_index - 0.5
    x_max <- start_index + chunk_size - 0.5
    
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
    face_styles <- setNames(
      ifelse(ordered_levels == proteus_row, "bold", "plain"),
      ordered_levels
    )
    
    # Create custom breaks for x-axis (every 20 positions)
    # Ensure break starting point is consistent between plots
    first_break <- ceiling(start_index / 20) * 20
    
    # Check if first_break is less than or equal to end_index before creating sequence
    if (first_break <= end_index) {
      x_breaks <- seq(
        from = first_break,
        to = end_index,
        by = 20
      )
    } else {
      # If first_break is greater than end_index, just use end_index
      x_breaks <- end_index
    }
    
    # For the final chunk that might be shorter, we still maintain the same scale
    # Importantly: use fixed limits for EACH chunk to ensure consistent scaling
    p <- ggplot(plot_chunk, aes(x = position, y = species_name)) +
      geom_tile(
        aes(fill = I(color)), 
        color = "white", 
        width = 0.9,
        height = 0.8
      ) +
      geom_text(
        aes(label = display_text, color = I(text_color)), 
        size = scaled_font_size  # Use the calculated scaled font size instead of fixed size 2
      ) +
      scale_x_continuous(
        limits = c(x_min, x_max),  # Fixed per-chunk limits for consistent scaling
        breaks = x_breaks,
        expand = c(0, 0)
      ) +
      scale_y_discrete(limits = rev(ordered_levels)) +
      labs(x = "", y = "") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(size = 6, angle = 0, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(5, 10, 5, 10),
        axis.text.y = element_text(
          size = scaled_font_size_y,,
          face = face_styles[as.character(rev(ordered_levels))]
        )
      ) +
      # Force square tiles with coord_fixed
      coord_fixed(ratio = tile_aspect_ratio)
    
    return(p)
  })
  
  # Create a custom heatmap legend
  legend_data <- data.frame(
    percentage = seq(0, 100, by = 10)
  )
  
  # Scale legend text size consistently with the main plot
  legend_text_size <- 1.5 * font_scaling_factor
  
  legend_plot <- ggplot(legend_data, aes(x = 1, y = percentage, fill = get_red_gradient(percentage))) +
    geom_tile() +
    scale_fill_identity() +
    labs(x = "", y = "Positive Selection Probability") +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_text(
        size = scaled_font_size_y,,
        face = ifelse(
          grepl("^Proteus anguinus", levels(alignment_df$species_name)), 
          "bold", 
          "plain"
        )
      ),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(5, 5, 5, 5)
    ) +
    geom_text(aes(label = paste0(percentage, "%")), 
              color = "black", 
              size = legend_text_size, 
              hjust = 0.5, 
              vjust = -1)
  
  # Use gtable to ensure all plots have identical widths
  # This is critical for consistent scaling across all plot rows
  plot_grobs <- lapply(plot_list, ggplotGrob)
  
  # Find the maximum width among all grobs
  max_width <- do.call(grid::unit.pmax, lapply(plot_grobs, function(g) g$widths))
  
  # Set all grobs to have the same width
  for (i in seq_along(plot_grobs)) {
    plot_grobs[[i]]$widths <- max_width
  }
  
  # Calculate overall dimensions based on number of chunks
  total_height <- max_plot_height * num_chunks
  
  # Use grid.arrange with the standardized grobs
  arranged_plot <- do.call(
    gridExtra::grid.arrange, 
    c(plot_grobs, ncol = 1)
  )
  
  # Combine the plot and legend side by side using cowplot
  combined_plot <- plot_grid(
    arranged_plot,
    legend_plot,
    ncol = 2,
    rel_widths = c(2, 0.2)  # Set the legend width to be thin
  )
  
  # Save plot with adjusted dimensions
  export_format <- match.arg(export_format)
  plot_height <- min(total_height, max_plot_height * 5)  # Limit maximum height
  plot_width <- chunk_plot_width + 2  # Add space for legend
  
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

generate_and_save_aa_alignment_visuals <- function(current_directory, gene_id){
  tryCatch({
    cat(paste0("Processing gene: ", gene_id, "\n"))
    # STEP 1: PREPROCESSING - read in the data  
    aa_directory_path <- file.path(current_directory, "www", "data", 
                                   gene_id, "translation",
                                   "prank.best.fas.translation")
    
    # Check if the file exists
    if (!file.exists(aa_directory_path)) {
      message(paste0("The file does not exist: ", aa_directory_path))
      return(FALSE)  # Return FALSE to indicate failure
    }
    
    message(paste0("The file exists: ", aa_directory_path))
    
    # Load the alignment using Biostrings package
    aas <- Biostrings::readAAStringSet(aa_directory_path, format = "fasta")
    
    # Get sequence names and lengths
    seq_names <- names(aas)
    seq_lengths <- sapply(aas, length)
    
    # Convert dataset to matrix
    aa_alignment_matrix <- as.matrix(aas)
    
    rownames(aa_alignment_matrix) <- gsub("_", " ", rownames(aa_alignment_matrix))
    
    # Create the dataframe
    aa_alignment_df <- data.frame(
      species_name = rownames(aa_alignment_matrix),
      aa_alignment_matrix
    )
    
    # Read the annotation file content
    aa_annotations_directory <- file.path(current_directory,
                                          "www", "data",
                                          gene_id, "translation", 
                                          "prank.best.fas.translation.view.annotations")
    
    if (!file.exists(aa_annotations_directory)) {
      message(paste0("Annotation file does not exist: ", aa_annotations_directory))
      return(FALSE)
    }
    
    # Read file content safely
    file_content <- tryCatch({
      readr::read_file(aa_annotations_directory)[1]
    }, error = function(e) {
      message(paste0("Error reading annotation file: ", e$message))
      return(NULL)
    })
    
    if (is.null(file_content)) {
      return(FALSE)
    }
    
    aa_annotation_df <- create_annotation_information_df(file_content)
    
    # positive selection indicator
    # Extract and print the NO_GRAPH content along with its length
    aa_pos_selection_annot <- aa_annotation_df$analyzed_positions[1]
    
    # Check if annotation is valid
    if (is.null(aa_pos_selection_annot) || is.na(aa_pos_selection_annot)) {
      message("Invalid or missing analyzed positions annotation")
      return(FALSE)
    }
    
    # clean up the no_graph line
    aa_pos_selection_clean_info <- remove_first_pipe_after_E(aa_pos_selection_annot)
    aa_pos_selection_clean_info_vector <- strsplit(aa_pos_selection_clean_info, "")[[1]]
    
    # Add title to the beginning of vector
    aa_pos_selection_clean_info_vec <- c("Analyzed Sites", aa_pos_selection_clean_info_vector)
    
    # Calculate the difference in length
    length_diff <- ncol(aa_alignment_df) - length(aa_pos_selection_clean_info_vec)
    
    # Extend the vector with "|" if it's shorter
    if (length_diff > 0) {
      aa_pos_selection_clean_info_vec <- c(aa_pos_selection_clean_info_vec, rep("|", length_diff))
    } else if (length_diff < 0) {
      # If the vector is too long, truncate it
      aa_pos_selection_clean_info_vec <- aa_pos_selection_clean_info_vec[1:ncol(aa_alignment_df)]
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
    
    # Check if annotation is valid
    if (is.null(aa_pos_selection_percentage_annot) || is.na(aa_pos_selection_percentage_annot)) {
      message("Invalid or missing positive selection probability annotation")
      return(FALSE)
    }
    
    # extract the percentage annotations
    aa_pos_selection_percentage_annot_vector <- get_positive_selection_percentages(aa_pos_selection_percentage_annot)
    aa_pos_selection_percentage_annot_vector <- c("Positive Selection Probability", aa_pos_selection_percentage_annot_vector)
    
    # Calculate the difference in length
    length_diff <- ncol(aa_combined_alignment_df) - length(aa_pos_selection_percentage_annot_vector)
    
    # Extend the vector with "|" if it's shorter
    if (length_diff > 0) {
      aa_pos_selection_percentage_annot_vector <- c(aa_pos_selection_percentage_annot_vector, rep("0.0", length_diff))
    } else if (length_diff < 0) {
      # If the vector is too long, truncate it
      aa_pos_selection_percentage_annot_vector <- aa_pos_selection_percentage_annot_vector[1:ncol(aa_combined_alignment_df)]
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
    aa_output_dir = file.path(current_directory, "www", "data",
                              gene_id, "translation")
    output_file <- file.path(aa_output_dir, "alignment.csv")
    
    aa_combined_alignment_df <- aa_combined_alignment_df[!rownames(aa_combined_alignment_df) %in% c("Analyzed Sites", "Positive Selection Probability"), ]
    
    
    # Create output directory if it doesn't exist
    if (!dir.exists(aa_output_dir)) {
      dir.create(aa_output_dir, recursive = TRUE)
    }
    
    # Save the dataframe as a CSV
    #write.csv(aa_combined_alignment_df, file = output_file, row.names = FALSE)
    
    # Generate and save the plot
    tryCatch({
      #plot_and_save_aa_alignment(
      #  alignment_df = aa_combined_alignment_df,
      #  output_dir = aa_output_dir,
      #  output_file = "alignment", 
      #  export_format = "png"
      #)
      
      # Save the dataframe as a CSV
      write.csv(aa_combined_alignment_df, file = output_file, row.names = FALSE) 
      
      cat(paste0("AA alignment csv created and saved for gene: ", gene_id, "\n"))
      return(TRUE)  # Return TRUE to indicate success
    }, error = function(e) {
      message(paste0("Error creating csv for gene: ", gene_id, ": ", e$message))
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
      "plot_and_save_aa_alignment",
      "generate_and_save_aa_alignment_visuals"
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
        success <- generate_and_save_aa_alignment_visuals(current_directory, gene_id)
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
  gene_ids <- gene_data_df$`Gene ID`
  #gene_ids <- c('TRINITY_DN12397_c0_g1')
  
  
  # Process all genes in parallel
  summary <- generate_alignment_visuals_parallel(
    current_directory = current_directory,
    gene_ids = gene_ids,
    #num_cores = 4  # Adjust based on your machine's capabilities
  )
  
  cat("All alignment csv files have been processed\n")
}

# Run the main function
main()