# app.R
library(shiny)
library(htmlwidgets)
library(DT)
library(ggplot2)
library(plotly)
library(dplyr)
library(matrixStats)
library(purrr)
library(readr)
library(reshape2)
library(reactable)
library(grid)
library(seqinr)
library(tidyr)
library(gridExtra)
library(stringr)
library(scales)
library(htmltools)
library(shinyjs)
library(R.utils)
library(cowplot)
#library(ggtext)

# Get the current working directory
current_directory <- getwd()
setwd(current_directory)

# [Previous functions remain unchanged]
robust_numeric_conversion <- function(df, columns) {
  for (col in columns) {
    tryCatch({
      df[[col]] <- suppressWarnings(as.numeric(ifelse(df[[col]] == "-", NA, df[[col]])))
    }, error = function(e) {
      cat("Error in conversion for column", col, ":", e$message, "\n")
    })
  }
  return(df)
}

format_pvalues <- function(x) {
  result <- ifelse(
    is.na(x),
    NA,
    ifelse(
      x < 1e-2,
      sprintf("%.2e", x),
      sprintf("%.2f", x)
    )
  )
  return(result)
}

# Modify the plot_aa_alignment function to ensure readable text size and spacing
plot_aa_alignment <- function(alignment_df, seq_length = NULL, chunk_size = 50) {
  if (nrow(alignment_df) == 0) {
    stop("Alignment dataframe is empty")
  }
  
  if (is.null(seq_length)) {
    seq_length <- ncol(alignment_df) - 1
  }
  
  # Define color palette for amino acids
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
  
  # Extract Positive Selection Probability values
  positive_selection_values <- alignment_df %>%
    filter(species_name == "Positive Selection Probability") %>%
    select(-species_name) %>%
    pivot_longer(cols = everything(), names_to = "position", values_to = "positive_selection") %>%
    mutate(position = as.numeric(gsub("X", "", position)))
  
  # Reorder rows
  alignment_df <- alignment_df %>%
    mutate(
      row_order = case_when(
        grepl("Positive", species_name) ~ 3,
        grepl("Analyzed", species_name) ~ 2,
        grepl("Proteus", species_name) ~ 1,
        TRUE ~ 0
      )
    ) %>%
    arrange(row_order, species_name) %>%
    mutate(
      species_name = factor(
        species_name,
        levels = c(
          unique(species_name[row_order == 0]),
          unique(species_name[row_order == 1]),
          unique(species_name[row_order == 2]),
          unique(species_name[row_order == 3])
        )
      )
    ) %>%
    select(-row_order)
  
  # Split the data into chunks
  chunks <- seq(1, seq_length, by = chunk_size)
  num_chunks <- length(chunks)
  
  # Calculate optimal square size based on number of positions in chunk
  square_size <- 25  # Base size, adjust as needed for readability
  
  # Generate a list of Plotly plots, one for each chunk
plot_list <- lapply(seq_along(chunks), function(chunk_idx) {
  start_index <- chunks[chunk_idx]
  end_index <- min(start_index + chunk_size - 1, seq_length)
  
  # Calculate how many positions are in this chunk
  positions_in_chunk <- end_index - start_index + 1
  
  # Calculate width proportional to the number of positions in this chunk
  # This ensures the last plot will be shorter if it has fewer positions
  chunk_width <- positions_in_chunk * square_size
  
  # Prepare data for this chunk
  plot_data <- alignment_df %>%
    pivot_longer(
      cols = all_of(paste0("X", start_index:end_index)),
      names_to = "position",
      values_to = "amino_acid"
    ) %>%
    mutate(position = as.numeric(gsub("X", "", position))) %>%
    left_join(positive_selection_values, by = "position")
  
  # Add a helper column to safely convert positive_selection to numeric values for coloring
  plot_data <- plot_data %>%
    mutate(
      selection_value = as.numeric(positive_selection),
      # Now handle coloring for each row type
      color = case_when(
        # For Positive Selection row, apply gradient based on selection value
        grepl("Positive Selection Probability", species_name) & !is.na(selection_value) ~ 
          colorRampPalette(c("white", "red"))(100)[pmin(100, pmax(1, as.integer(selection_value)))],
        # For all other rows, use amino acid colors
        amino_acid %in% names(aa_colors) ~ aa_colors[amino_acid],
        TRUE ~ "rgb(255,255,255)"
      ),
      text_value = case_when(
        grepl("Positive Selection Probability|Analyzed Sites", species_name) ~ "",
        TRUE ~ amino_acid
      ),
      hover_text = case_when(
        grepl("Positive Selection Probability", species_name) ~ 
          paste0("Position: ", position, "<br>Selection Probability: ", positive_selection, "%"),
        grepl("Analyzed Sites", species_name) ~ 
          paste0(""),
        TRUE ~ paste0(
          "Position: ", position, "<br>",
          "Species: ", species_name, "<br>",
          "Amino Acid: ", amino_acid, "<br>",
          "Positive Selection: ", ifelse(is.na(positive_selection), "N/A", paste0(positive_selection, "%"))
        )
      )
    )
  
  # Create tickvals and ticktext that show only every 25th position
  index_positions <- seq(start_index, end_index)
  tickvals <- index_positions[index_positions %% 25 == 0 | index_positions == start_index]
  ticktext <- as.character(tickvals)
  
  # Generate the Plotly scatter plot with square markers
  p <- plot_ly(
    data = plot_data,
    x = ~position,
    y = ~species_name,
    type = "scatter",
    mode = "markers+text",
    marker = list(
      symbol = "square", 
      size = square_size,
      color = ~color,
      line = list(color = "black", width = 0.5)  # Thinner border
    ),
    text = ~text_value,
    textfont = list(
      size = 12,
      weight = 700,  # Bold text (700 = bold)
      color = "black"
    ),
    textposition = "middle",
    hoverinfo = "text",
    hovertext = ~hover_text,
    hoverlabel = list(
      bgcolor = "black", 
      font = list(color = "white")
    )
  ) %>%
    layout(
      width = chunk_width,  # Width is now proportional to positions
      height = 150 + nrow(alignment_df) * 30,  # Fixed height for consistent row spacing
      xaxis = list(
        title = NULL, 
        showgrid = FALSE, 
        zeroline = FALSE,
        tickvals = tickvals,
        ticktext = ticktext,
        tickfont = list(size = 10)
      ),
      yaxis = list(
        title = NULL, 
        showgrid = FALSE,
        zeroline = FALSE,
        autorange = "reversed",
        tickfont = list(size = 10)
      ),
      plot_bgcolor = "white",
      margin = list(l = 80, r = 40, b = 40, t = 20)  # Consistent margins
    )
  
  return(p)
})
}

# [Data loading code remains unchanged]
gene_expression_df <- read.csv(file.path(current_directory, "www", "data","resources", "gene_data_table.csv"),
                               header = TRUE,
                               na = c("", "NA", "N/A", "-"),
                               stringsAsFactors = FALSE)

raw_data_df <- read.csv(file.path(current_directory, "www", "data", "resources", "proteus_gene_exp_box_plot_data.csv"), 
                        header = TRUE, 
                        stringsAsFactors = FALSE, 
                        check.names = FALSE)

# Convert columns to numeric, replacing "-" with NA
numeric_cols <- c("Brain", "Gut", "Heart", "Liver", "Lung", "Skin", 
                  "X.Sites.under.selection", "dn.ds", "p.value", "FDR")

# Convert specified columns to numeric
for (col in numeric_cols) {
  # Replace '-' with NA explicitly
  gene_expression_df[[col]] <- ifelse(gene_expression_df[[col]] == "-", NA, gene_expression_df[[col]])
  gene_expression_df[[col]] <- as.numeric(gene_expression_df[[col]])
  
}

names(gene_expression_df) <- c("Gene ID", "Gene Symbol", "Brain", "Gut", "Heart", "Liver", "Lung", "Skin", 
                               "p-value", "FDR", "dn/ds", "#Sites under selection")

# Round numeric columns to 2 decimal places
numeric_cols_without_pvalue <- c("Brain", "Gut", "Heart", "Liver", "Lung", "Skin", 
                                 "FDR", "dn/ds", "#Sites under selection")
gene_expression_df[numeric_cols_without_pvalue] <- round(gene_expression_df[numeric_cols_without_pvalue], 2)

# Reorder columns
new_order <- c("Gene ID", "Gene Symbol", "Brain", "Gut", "Heart", "Liver", "Lung", "Skin", 
               "#Sites under selection", "dn/ds", "p-value", "FDR")
gene_expression_df <- gene_expression_df[, new_order]

gene_expression_df <- gene_expression_df[order(gene_expression_df$'dn/ds', decreasing = TRUE), ]



# UI definition
ui <- fluidPage(
  tags$head(
    # Load local CSS
    tags$link(rel = "stylesheet", type = "text/css", href = "css/styles.css"),
    # Load local JS
    tags$script(src = "js/formatters.js"),
    tags$script(HTML("
    $(document).ready(function() {
      // Add a global variable to track current view
      window.currentAlignmentView = 'png';
      
      // Function to handle view toggle
      window.toggleAlignmentView = function() {
        var currentView = window.currentAlignmentView;
        if (currentView === 'png') {
          Shiny.setInputValue('alignmentViewToggle', 'plotly');
        } else {
          Shiny.setInputValue('alignmentViewToggle', 'png');
        }
      };
    });
  ")),
    
    # Add overlay styles
    tags$style(HTML("
    .png-container {
      border: 1px solid #ddd;
      padding: 10px;
      margin-top: 20px;
      background-color: white;
      position: relative;
      overflow: auto;
    }
    .tooltip-custom {
      position: absolute;
      padding: 8px;
      background: rgba(0, 0, 0, 0.8);
      color: white;
      border-radius: 4px;
      font-size: 14px;
      z-index: 1000;
      pointer-events: none;
      transition: opacity 0.3s;
      max-width: 250px;
    }
    .overlay-grid {
      position: absolute;
      top: 0;
      left: 0;
      z-index: 10;
      pointer-events: none;
    }
    .cell-overlay {
      position: absolute;
      background-color: transparent;
      border: 1px solid transparent;
      cursor: pointer;
      z-index: 20;
      pointer-events: auto;
    }
    .cell-overlay:hover {
      border: 1px solid red;
      background-color: rgba(255, 0, 0, 0.1);
    }
  ")),
    
    # Add tooltip JavaScript
    tags$script(HTML("
    $(document).ready(function() {
      // Create tooltip div if it doesn't exist
      if ($('#aa-tooltip').length === 0) {
        $('body').append('<div id=\"aa-tooltip\" class=\"tooltip-custom\" style=\"display:none;\"></div>');
      }
      
      // Tooltip functionality
      $(document).on('mouseenter', '.cell-overlay', function(e) {
        var $this = $(this);
        var position = $this.attr('data-position');
        var row = $this.attr('data-row');
        var aminoAcid = $this.attr('data-amino-acid');
        var percentage = $this.attr('data-percentage');
        
        var tooltipContent = '';
        
        if (position && row) {
          tooltipContent = '<strong>Position:</strong> ' + position + '<br>' +
                          '<strong>Species:</strong> ' + row;
          
          if (aminoAcid && aminoAcid !== 'NA') {
            tooltipContent += '<br><strong>Amino Acid:</strong> ' + aminoAcid;
          }
          
          if (percentage && percentage !== 'NA') {
            tooltipContent += '<br><strong>Selection Probability:</strong> ' + percentage + '%';
          }
          
          $('#aa-tooltip').html(tooltipContent).show();
        }
      });
      
      $(document).on('mouseleave', '.cell-overlay', function() {
        $('#aa-tooltip').hide();
      });
      
      $(document).on('mousemove', function(e) {
        $('#aa-tooltip').css({
          left: e.pageX + 10,
          top: e.pageY + 10
        });
      });
    });
  "))
  ),
  tags$script(src = "js/imageZoom.js"),
  
  shinyjs::useShinyjs(),
  
  div(
    style = "position: relative; padding: 15px 0;",
    tags$a(
      href = "https://www.phenomics.ruhr-uni-bochum.de/",
      target = "_blank",
      tags$img(
        src = "phenocomp_logo.png",
        style = "position: absolute; top: 15px; right: 30px; height: 100px; z-index: 1000;"
      )
    ),
    titlePanel("Proteus anguinus Transcriptome Analysis")
  ),
  
  
  # Wrap existing content in tabsetPanel
  tabsetPanel(
    tabPanel("Analysis",
             fluidRow(
               column(12,
                      h3("Gene Expression Levels"),
                      DTOutput("geneTable"),
                      hr()
               )
             ),
             
             fluidRow(
               column(12,
                      h3("Expression Analysis"),
                      plotlyOutput("boxplot"),
                      hr()
               )
             ),
             
             fluidRow(
               column(12,
                      h3("Alignment Visualization"),
                      div(
                        style = "display: flex; justify-content: space-between; align-items: center; margin: 20px 0;",
                        div(
                          style = "display: flex; align-items: center;",
                          selectInput("seqType",
                                      "Select Alignment Type:",
                                      choices = c("Amino acid" = "aa",
                                                  "Coding sequence" = "cds"),
                                      selected = "aa",
                                      width = "300px"),
                          # New toggle button placed right next to the selector
                          actionButton("toggleAlignmentView", 
                                       "Analyze Alignment", 
                                       icon = icon("exchange-alt"),
                                       style = "margin-left: 10px; height: 34px;")
                        ),
                        div(
                          downloadButton("downloadAlignment", "Download Alignment FASTA")
                        ),
                      ),
                      # Use a conditional panel or uiOutput to switch between PNG and Plotly
                      uiOutput("alignmentDisplay"),
                      textOutput("errorMessage")
               )
             )
    ),
    
    tabPanel("Download",
             fluidRow(
               column(12,
                      h3("Download Options"),
                      div(
                        style = "display: flex; flex-direction: column;",
                        # download links
                        downloadLink("downloadCodingSequences", "Coding Sequences"),
                        downloadLink("downloadGeneDataTable", "Gene Expression Table"),
                        downloadLink("downloadTranscriptomeAssembly", "Transcriptome Assembly")
                      )
               )
             )
    )
  )
)

# Server logic
server <- function(input, output, session) {
  
  runjs <- shinyjs::runjs
  
  # Reactive value for selected gene
  selected_gene <- reactiveVal(gene_expression_df$'Gene ID'[1])
  
  
  # Reactive values to store cell information
  rv <- reactiveValues(
    grid_cells = NULL,
    image_dimensions = NULL,
    alignment_data = NULL
  )
  
  # Reactive value to track the current view mode
  view_mode <- reactiveVal("png")
  
  # Observe the toggle button
  observeEvent(input$toggleAlignmentView, {
    current_mode <- view_mode()
    new_mode <- if(current_mode == "png") "plotly" else "png"
    view_mode(new_mode)
  })
  
  # Observe changes in the data table selection to reset to PNG
  observeEvent(input$geneTable_rows_selected, {
    view_mode("png")
  })
  
  # Reactive to get the alignment file path
  alignment_file_path <- reactive({
    req(selected_gene(), input$seqType)
    
    # Construct the file path
    file.path(
      "www",
      "data", 
      selected_gene(),
      "translation",
      "alignment.csv"
    )
  })
  
  # Reactive to read alignment data
  alignment_data <- reactive({
    if(view_mode() == "plotly") {
      file_path <- alignment_file_path()
      
      # Check if file exists
      if(!file.exists(file_path)) {
        showNotification(
          paste("Alignment file not found:", file_path), 
          type = "error"
        )
        return(NULL)
      }
      
      # Read the CSV file
      alignment_df <- read.csv(file_path, stringsAsFactors = FALSE)
      
      return(alignment_df)
    }
    return(NULL)
  })
  
  # Conditional output for alignment display
  output$alignmentDisplay <- renderUI({
    req(view_mode())
    
    if(view_mode() == "png") {
      # Render PNG images as before
      uiOutput("imageOutput")
    } else {
      # Render Plotly plot with horizontal scrolling
      div(
        style = "width: 100%; overflow-x: auto; overflow-y: visible;", 
        plotlyOutput("alignmentPlotly", height = "auto")
      )
    }
  })
  
  # New Plotly plot rendering
  output$alignmentPlotly <- renderPlotly({
    # Get the alignment data
    alignment_df <- alignment_data()
    
    # Check if data is available
    req(alignment_df)
    
    # Dynamically calculate sequence length from the number of columns in the alignment
    # Subtract 1 to account for the species_name column
    seq_length <- ncol(alignment_df) - 1
    
    # Calculate the number of amino acids per chunk
    chunk_size <- 50  # You can adjust this value
    
    # Calculate square size for amino acids
    square_size <- 30  # Base size for each amino acid cell
    
    # Your custom function to generate Plotly visualization
    plot_result <- plot_aa_alignment(alignment_df, seq_length = seq_length, chunk_size = chunk_size)
    
    # If it's a list of plots, subplot them with dynamic height calculation
    if (is.list(plot_result)) {
      # Calculate the height based on the number of plots and species
      num_plots <- length(plot_result)
      num_species <- nrow(alignment_df)
      
      # Base height per plot plus additional height per species
      plot_base_height <- 150  # base height per plot in pixels
      height_per_species <- 30  # additional height per species in pixels
      
      total_height <- num_plots * (plot_base_height + num_species * height_per_species)
      
      # Calculate total width needed for each chunk
      positions_per_chunk <- min(chunk_size, seq_length)
      chunk_width <- max(1000, positions_per_chunk * square_size)
      
      # Apply config to each plot in the list before creating subplot
      plot_result <- lapply(plot_result, function(p) {
        p %>% layout(
          autosize = FALSE,
          width = chunk_width,
          height = plot_base_height + num_species * height_per_species,
          margin = list(l = 50, r = 50, b = 50, t = 50, pad = 4),
          dragmode = FALSE,
          showlegend = FALSE
        ) %>%
          config(
            displayModeBar = FALSE,
            displaylogo = FALSE,
            editable = FALSE
          )
      })
      
      # Create the subplot with the calculated height
      subplot(plot_result, nrows = length(plot_result), shareX = FALSE, shareY = TRUE) %>% 
        layout(
          autosize = FALSE,
          width = chunk_width + 100,  # Add some margin
          height = total_height,
          showlegend = FALSE
        )
    } else {
      # For a single plot, set height based on number of species
      num_species <- nrow(alignment_df)
      plot_height <- 150 + num_species * 30  # Adjust as needed
      
      # Calculate width based on sequence length
      positions_count <- min(seq_length, chunk_size)
      plot_width <- max(1000, positions_count * square_size)
      
      # Apply config to single plot
      plot_result %>% 
        layout(
          autosize = FALSE,
          width = plot_width,
          height = plot_height,
          margin = list(l = 50, r = 50, b = 50, t = 50, pad = 4),
          dragmode = FALSE,
          showlegend = FALSE
        ) %>%
        config(
          displayModeBar = FALSE,
          displaylogo = FALSE,
          editable = FALSE
        )
    }
  })
  
  
  # In your renderDataTable function, modify the options:
  output$geneTable <- DT::renderDataTable({
    
    sketch = htmltools::withTags(table(
      class = 'display',
      thead(
        tr(
          th(colspan = 2, 'Basic Information'),
          th(colspan = 6, 'Gene Expression'),
          th(colspan = 4, 'Positive Selection')
        ),
        tr(
          lapply(colnames(gene_expression_df), th)
        )
      )
    ))
    
    # Columns to handle
    numeric_cols <- c("Brain", "Gut", "Heart", "Liver", "Lung", "Skin", 
                      "#Sites under selection", "dn/ds", "p-value", "FDR")
    
    DT::datatable(gene_expression_df,
                  container = sketch,
                  rownames = FALSE,
                  escape = FALSE, 
                  selection = 'single',
                  options = list(
                    pageLength = 10,
                    processing = FALSE,      # Disable processing indicator
                    language = list(         # Remove loading text
                      processing = ''
                    ),
                    dom = 'lfrtip',         # Standard display without processing indicator
                    drawCallback = JS("function() { 
                      $('.dataTables_processing').hide(); 
                    }"),
                    columnDefs = list(
                      # Specific p-value formatting
                      list(
                        targets = 10,  # assuming p-value is column index 10
                        render = JS("formatPValue")
                      ),
                      list(  ## todo - this part is overridden by html of container/sketch
                        targets = 2,  # assuming GENE SYMBOL is column index 2
                        render = JS("formatGeneSymbol")
                      ),
                      list(
                        targets = 9,  # assuming dn/ds is column index 9
                        render = JS("formatNumber")
                      ),
                      list(
                        targets = 8,  # assuming sites under selection is column index 8
                        render = JS("formatNumber")
                      ),
                      list(
                        targets = 11,  # assuming fdr is column index 11
                        render = JS("formatNumber")
                      )
                    )
                  )
    )
  })
  
  # Update selected gene when table selection changes
  observeEvent(input$geneTable_rows_selected, {
    if (!is.null(input$geneTable_rows_selected)) {
      selected_gene(gene_expression_df$'Gene ID'[input$geneTable_rows_selected])
    }
  })
  
  # Render boxplot
  output$boxplot <- renderPlotly({
    
    # Find the matching row in proteus_gene_exp_df using Gene ID
    gene_row_index <- grep(paste0("^", selected_gene()), raw_data_df$`id`)

    # Get the selected gene's expression data
    selected_data <- raw_data_df[gene_row_index, ]
    
    # Reshape the data for boxplot
    tissue_names <- sub("\\..*", "", colnames(raw_data_df)[2:ncol(raw_data_df)])
    unique_tissues <- unique(tissue_names)
    
    # Prepare data for plotting
    plot_data <- data.frame(
      Tissue = character(),
      Expression = numeric()
    )
    
    for (tissue in unique_tissues) {
      tissue_cols <- grep(paste0("^", tissue, "\\."), colnames(raw_data_df), value = TRUE)
      tissue_expressions <- selected_data[, tissue_cols]
      
      plot_data <- rbind(plot_data, 
                         data.frame(
                           Tissue = tissue,
                           Expression = as.numeric(tissue_expressions)
                         )
      )
      
    }
    
    # Calculate medians for each tissue
    median_data <- plot_data %>%
      group_by(Tissue) %>%
      summarize(Median = median(Expression, na.rm = TRUE)) %>%
      ungroup()
    
    # Generate the scatter plot with median lines
    p <- plot_ly() %>%
      # Add scatter points
      add_trace(
        data = plot_data,
        x = ~Tissue,
        y = ~Expression,
        type = "scatter",
        mode = "markers",
        marker = list(size = 10),
        name = "Samples"
      ) %>%
      # Add median horizontal lines
      add_trace(
        data = median_data,
        x = ~Tissue,  # Use the tissue as x
        y = ~Median,  # Use the median as y line
        type = "scatter",
        mode = "markers",
        marker = list(
          symbol = "line-ew",  # Use dash symbol
          size = 15,  # Adjust size as needed
          line = list(color = "red", width = 2)
        ),
        name = "Median",
        showlegend = TRUE,
        hoverinfo = "y",
        error_x = list(type = "data", symmetric = TRUE, width = 0.2)
      ) %>%
      layout(
        title = paste("Expression of", selected_gene()),
        #xaxis = list(title = "Tissue Type"),
        xaxis = list(title = ""),
        yaxis = list(title = "Expression Level (TPM)"),
        showlegend = TRUE
      )
    
  })
  
  # Function to get the visualization file path 
  getVisualizationPath <- function(gene_id, type) {
    base_path <- file.path("www", "data", gene_id)
    
    if(type == "aa") {
      file_path <- file.path(base_path, "translation", "alignment.png")
    } else {
      file_path <- file.path(base_path, "nucleic_acid", "alignment.png")
    }
    
    return(file_path)
  }
  
  # Function to check if visualization file exists
  checkFile <- function(gene_id, type) {
    file_path <- getVisualizationPath(gene_id, type)
    return(file.exists(file_path))
  }
  
  # Function to get the web-accessible path for the image
  getWebPath <- function(gene_id, type) {
    # Convert file system path to web path
    if(type == "aa") {
      return(file.path("data", gene_id, "translation", "alignment.png"))
    } else {
      return(file.path("data", gene_id, "nucleic_acid", "alignment.png"))
    }
  }
  
  # Reactive expression for file checking
  fileStatus <- reactive({
    req(selected_gene())
    checkFile(selected_gene(), input$seqType)
  })
  
  # Output for image display
  output$imageOutput <- renderUI({
    req(selected_gene())
    if(fileStatus()) {
      web_path <- getWebPath(selected_gene(), input$seqType)
      div(
        style = "margin-top: 20px;",
        div(
          id = "imageContainer",
          class = "png-container",
          # Image display
          img(
            id = "alignmentImage",
            src = web_path,
            alt = paste(if(input$seqType == "aa") "Amino acid" else "Coding sequence",
                        "alignment for", selected_gene())
          ),
          # Overlay container
          tags$div(id = "overlayContainer", class = "overlay-grid")
        )
      )
    }
  })
  
  # Output for error message
  output$errorMessage <- renderText({
    req(selected_gene())
    if(!fileStatus()) {
      paste("Alignment visualization not found for", selected_gene(),".")
    }
  })
  
  # Set initial table selection
  observeEvent(1, {
    dataTableProxy("geneTable") %>%
      selectRows(1)
  }, once = TRUE)
  
  # download handler for coding sequences
  output$downloadCodingSequences <- downloadHandler(
    filename = function() { "proteus_anguinus_coding_sequences.tar.gz"},
    
    content = function(file) {
      file.copy(file.path(current_directory, 
                          "www", 
                          "download_data", 
                          "proteus_anguinus_coding_sequences.tar.gz"), 
                file) 
    },
    contentType = "tar.gz" 
  )
  
  # download handler for transcriptome assembly
  output$downloadTranscriptomeAssembly <- downloadHandler(
    filename = function() { "proteus_anguinus_transcriptome_assembly.tar.gz"},
    
    content = function(file) {
      file.copy(file.path(current_directory, 
                          "www", 
                          "download_data", 
                          "proteus_anguinus_transcriptome_assembly.tar.gz"), 
                file) 
    },
    contentType = "tar.gz" 
  )
  
  # download handler for complete gene data table
  output$downloadGeneDataTable <- downloadHandler(
    filename = function() { "proteus_anguinus_gene_expression_positive_selection_matrix.tar.gz"},
    
    content = function(file) {
      file.copy(file.path(current_directory, 
                          "www", 
                          "download_data", 
                          "proteus_anguinus_gene_expression_positive_selection_matrix.tar.gz"), 
                file)
    },
    contentType = "tar.gz" 
  )
  
  # First, create the hidden UI elements for downloads
#  observe({
#    # Create two hidden download buttons
#    insertUI(
#      selector = "#downloadAlignment",
#      where = "afterEnd",
#      ui = div(
#        downloadButton("downloadAA", ""), 
#        style = "visibility:hidden;",
#        class = "hidden-div"
#      )
#    )
#    insertUI(
#      selector = "#downloadAlignment",
#      where = "afterEnd",
#      ui = div(
#        downloadButton("downloadCDS", ""), 
#        style = "visibility:hidden;",
#       class = "hidden-div"
#      )
#    )
#  })
  
#  # When the main button is clicked, trigger both hidden downloads
#  observeEvent(input$downloadAlignment, {
#    shinyjs::click("downloadAA")
#    #shinyjs::click("downloadCDS")
#  })
  
  
  
  # Handler for translation fasta
  output$downloadAlignment <- downloadHandler(
    filename = function() {
      if (input$seqType=="aa") paste0(selected_gene(), "_translation.fasta") else paste0(selected_gene(), "_coding_sequence.fasta")
    },
    content = function(file) {
      # Get file from translation directory
         if (input$seqType=="aa") file.copy(file.path(current_directory, 
                          "www", 
                          "data",
                          selected_gene(),
                          "translation",
                          "prank.best.fas.translation"), 
                file) else file.copy(file.path(current_directory, 
                                               "www", 
                                               "data",
                                               selected_gene(),
                                               "nucleic_acid",
                                               "prank.best.fas.renamed"), 
                                     file)
    }
)
  
  # Handler for coding sequence fasta
  output$downloadCDS <- downloadHandler(
    filename = function() {
      paste0(selected_gene(), "_coding_sequence.fasta")
    },
    content = function(file) {
      # Get file from nucleic_acid directory
      file.copy(file.path(current_directory, 
                          "www", 
                          "data",
                          selected_gene(),
                          "nucleic_acid",
                          "prank.best.fas.renamed"), 
                file)
    }
  )
  
}

# Run the app
shinyApp(ui = ui, server = server)