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


generate_and_save_aa_alignment_visuals <- function(current_directory, gene_id){
  
  cat(gene_id, "\n")
  # STEP 1: PREPROCESSING - read in the data  
  aa_directory_path <- file.path(current_directory, "www", "data", 
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
    
    #print(head(aa_alignment_df))
    
    # Read the annotation file content
    aa_annotations_directory <- file.path(current_directory,
                                          "www", "data",
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
    aa_output_dir = file.path(current_directory, "www", "data", 
                              gene_id, "translation")
    output_file <- file.path(aa_output_dir, "amino_acid_alignment.csv")
    
    #save the df
    # Save the dataframe as a CSV
    write.csv(aa_combined_alignment_df, file = output_file, row.names = FALSE)

    # STEP 2: VISUALIZATION - create the visuals for aa
    # Generate and save the plot
    plot_and_save_aa_alignment(
      alignment_df = aa_combined_alignment_df,
      output_dir = aa_output_dir,
      output_file = "alignment", 
      export_format = "png"
    )
    
    cat("AA alignment created and saved.", "\n")
  }
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

generate_alignment_visuals <- function(
    current_directory, 
    gene_id
){

  #generate aa alignment
  generate_and_save_aa_alignment_visuals(current_directory = current_directory, 
                                           gene_id = gene_id)
  #  
  #generate na alignment
  #generate_and_save_na_alignment_visuals(current_directory = current_directory,
  #                                         gene_id = gene_id)
  cat("Alignment files are created")
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

