# library(shiny)
# library(ggplot2)
# library(dplyr)
# library(shinyjs)
# library(shinyWidgets)
# library(shinybusy)
# library(data.table)

#' @title motifSimulationApp
#' @description Shows details of an object of class MyS4Class.
#' @param object An object of class MyS4Class.
#' @export
motifSimulationApp <- function(error_str,mot_details) {
  ui <- fluidPage(
    useShinyjs(),  # Initialize shinyjs
    useSweetAlert(),  # Initialize shinyWidgets
    
    # Custom CSS for styling
    tags$head(
      tags$style(HTML("
        /* General styling */
        body {
          background-color: #f8f9fa;
          color: #343a40;
          font-family: 'Arial', sans-serif;
        }
        .well {
          background-color: #ffffff;
          color: #343a40;
          border: 1px solid #ced4da;
          box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
          border-radius: 10px;
          padding: 20px;
        }
        .btn-primary {
          background-color: #007bff;
          border: none;
          padding: 10px 20px;
          font-size: 16px;
          cursor: pointer;
          border-radius: 5px;
          transition: background-color 0.3s, transform 0.2s;
          margin-right: 10px;
        }
        .btn-primary:hover {
          background-color: #0056b3;
          transform: scale(1.05);
        }
        .btn-primary:active {
          background-color: #004080;
          transform: scale(1.02);
        }
        #pagination_controls {
          display: flex;
          justify-content: center;
          margin: 20px 0;
        }
        .pagination-buttons {
          background-color: #007bff;
          border: none;
          color: white;
          padding: 10px 20px;
          font-size: 16px;
          border-radius: 5px;
          cursor: pointer;
          transition: background-color 0.3s, transform 0.2s;
          margin: 0 10px;
        }
        .pagination-buttons:disabled {
          background-color: #cccccc;
          cursor: not-allowed;
        }
        .pagination-buttons:hover:not(:disabled) {
          background-color: #0056b3;
          transform: scale(1.05);
        }
        .pagination-buttons:active:not(:disabled) {
          background-color: #004080;
          transform: scale(1.02);
        }
        .collapsible {
          cursor: pointer;
          padding: 10px;
          background-color: #007bff;
          border: none;
          color: white;
          text-align: left;
          font-size: 16px;
          width: 100%;
          border-radius: 5px;
          margin-bottom: 10px;
          transition: background-color 0.3s;
        }
        .collapsible:hover {
          background-color: #0056b3;
        }
        .content {
          display: none;
          padding: 10px;
          background-color: #f8f9fa;
          border-radius: 5px;
          margin-bottom: 20px;
        }
        .plot-container {
          margin-bottom: 30px;
        }
        .plot-title {
          text-align: center;
          margin-bottom: 10px;
          font-size: 20px;
          font-weight: bold;
        }
        .input-row {
          display: flex;
          flex-wrap: wrap;
          gap: 15px;
          margin-bottom: 15px;
        }
        .input-group {
          flex: 1;
          min-width: 220px;
        }
        #error_message {
          color: red;
          font-weight: bold;
          margin-top: 20px;
          text-align: center;
        }
      "))
    ),
    
    # Title and Main Panel
    titlePanel("Motif Simulation Plots"),
    
    mainPanel(
      # Dynamic UI for plots
      uiOutput("plots_ui"),
      
      # Plot and pagination buttons
      fluidRow(
        column(width = 12, 
               actionButton("plotBtn", "Generate Plots", class = "btn btn-primary"),
               tags$div(
                 id = "pagination_controls",
                 actionButton("prevPage", "Previous", class = "pagination-buttons", disabled = TRUE),
                 actionButton("nextPage", "Next", class = "pagination-buttons", disabled = TRUE)
               )
        )
      ),
      
      # Download buttons
      fluidRow(
        column(width = 12, 
               downloadButton("downloadPdf", "Download PDF", class = "btn btn-primary"),
               downloadButton("downloadRdata", "Download Rdata", class = "btn btn-primary"),
               div(id = "error_message", "")
        )
      ),
    ),
    
    # Input controls panel
    fluidRow(
      column(width = 12,
             wellPanel(
               textInput("path", "Directory Path to Save Plots", value = getwd()),
               div(class = "input-row",
                   div(class = "input-group", numericInput("N", "Number of Curves (N)", value = 20, min = 1)),
                   div(class = "input-group", numericInput("len", "Curve Length (len)", value = 300, min = 1))
               ),
               div(class = "input-row",
                   div(class = "input-group", numericInput("norder", "Order of Basis (norder)", value = 3, min = 1)),
                   div(class = "input-group", numericInput("coeff_min", "Minimum Coefficient (coeff_min)", value = -15))
               ),
               div(class = "input-row",
                   div(class = "input-group", numericInput("coeff_max", "Maximum Coefficient (coeff_max)", value = 15)),
                   div(class = "input-group", numericInput("dist_knots", "Distance Between Knots (dist_knots)", value = 10, min = 1))
               ),
               div(class = "input-row",
                   div(class = "input-group", numericInput("min_dist_motifs", "Minimum Distance Between Motifs (min_dist_motifs)", value = 30,min = 1)),
                   div(class = "input-group", selectInput("distribution", "Coefficient Distribution", choices = c("unif", "beta", "empirical"), selected = "beta"))
               ),
               
               # Conditional panel to show when "empirical" is selected
               conditionalPanel(
                 condition = "input.distribution == 'empirical'",
                 textInput("empirical_values", "Enter Empirical Coefficients (comma-separated)", value = "1, 2, 3, 4, 5")
               ),
               
               div(class = "input-row",
                   div(class = "input-group", selectInput("error_type", "Error Type", choices = c("pointwise", "coeff")))
               )
             )
      )
    )
  )
  
  server <- function(input, output, session) {
    
    observeEvent(input$plotBtn, {
      req(input$path)
      
      # Show progress bar
      show_modal_spinner(spin = "circle", text = "Generating plots, please wait...")
      
      tryCatch({
        N <- input$N
        len <- input$len
        norder <- input$norder
        coeff_min <- input$coeff_min
        coeff_max <- input$coeff_max
        dist_knots <- input$dist_knots
        min_dist_motifs <- input$min_dist_motifs
        error_type <- input$error_type
        
        # Handle empirical distribution
        if (input$distribution == "empirical") {
          distribution <- as.numeric(unlist(strsplit(input$empirical_values, ",")))
        } else {
          distribution <- input$distribution
        }
        
        if (!dir.exists(input$path)) {
          stop("Invalid directory path.")
        }
  
        builder <- mtfd::motifSimulationBuilder(N = N,len = len,
                                                mot_details = mot_details,
                                                norder = norder,
                                                coeff_min = coeff_min,
                                                coeff_max = coeff_max,
                                                dist_knots = dist_knots,
                                                min_dist_motifs = min_dist_motifs,
                                                distribution = distribution)
        curves <- mtfd::generateCurves(builder,error_type = error_type, error_str)
        
        output_file <- file.path(input$path, "plots.pdf")
        
        if (!dir.exists(input$path)) {
          dir.create(input$path)
        }
        n_error <- NULL
        plots <- lapply(seq_along(curves), function(k) {
          if(!is.null(curves[[k]]$with_error)) {
            n_error <<- length(curves[[k]]$with_error$error_y)
            curve_data_no_error <- data.frame(
              t = seq(0, curves[[k]]$basis$rangeval[2]-1),
              x = curves[[k]]$background$no_error_y,
              z = curves[[k]]$no_noise$motif_y)
            names(curve_data_no_error) <- c("t","x","z")
            
            curve_data_error <- NULL
            curve_data_error <- data.frame(
              t = seq(0, curves[[k]]$basis$rangeval[2]-1),
              x = curves[[k]]$with_error$error_y)
            names(curve_data_error) <- c("t",paste0("x",seq(length(curves[[k]]$with_error$error_y))))
            
            motif_lines <- mapply(function(id_motif, pos_motif, instance) {
              motif_t = seq((pos_motif - 1) * builder@dist_knots,
                            (pos_motif - 1) * builder@dist_knots + builder@mot_details[[id_motif]]$len)
              motif_x = lapply(curves[[k]]$with_error$error_y,function(curve){return(curve[motif_t + 1])})
              
              return(lapply(motif_x,function(motif){ data.frame(t = motif_t, x = motif, motif_id = factor(paste(id_motif, instance, sep = "_")),
                                                                initial_number = str_extract(as.character(id_motif), "^[^_]+"),
                                                                xmin = (pos_motif - 1) * builder@dist_knots,
                                                                xmax = (pos_motif - 1) * builder@dist_knots + builder@mot_details[[id_motif]]$len)}))
            }, builder@motifs_in_curves[[k]]$motif_id, builder@motifs_in_curves[[k]]$starting_coeff_pos, seq_along(builder@motifs_in_curves[[k]]$motif_id), SIMPLIFY = FALSE)
            
            motif_colors <- c( "1" = "red", "2" = "blue", "3" = "darkgreen", "4" = "orange",
                               "5" = "purple", "6" = "cyan", "7" = "magenta", "8" = "brown",
                               "9" = "pink", "10" = "grey")
            motif_colors <- rep(motif_colors,length.out = length(builder@mot_details))
            if(length(builder@mot_details) > 10 )
              attr(motif_colors,"names")[11:length(builder@mot_details)] <- as.character(as.integer(attr(motif_colors,"names")[11:length(builder@mot_details)]) + 10)
            
            max_dataframes <- max(sapply(motif_lines, function(sublist) length(sublist)))
            # Initialize a list to store results
            motif_data <- vector("list", max_dataframes)
            for (i in seq_len(max_dataframes)) {
              # Extract data frames at the i-th level
              dataframes_at_level_i <- lapply(motif_lines, function(sublist) sublist[[i]])
              # Combine the extracted data frames
              motif_data[[i]] <- bind_rows(dataframes_at_level_i)
              names(motif_data[[i]]) <- c("t", "x", "motif_id", "initial_number", "xmin", "xmax")
            }
            p <- lapply(1:length(motif_data), function(j) {
              # Create motif_id labels
              motif_labels <- paste("motif_id:", unique(motif_data[[j]]$initial_number))
              pic <- ggplot() +
                # Plot the main curve in gray30
                geom_line(data = curve_data_no_error, aes(x = t, y = x, color = 'background_curve'), linewidth = 0.5) +
                # Plot the curve with motifs in gold
                geom_line(data = curve_data_no_error, aes(x = t, y = z, color = 'with_motif'), linewidth = 0.5) +
                # Plot the error curve
                geom_line(data = curve_data_error, aes_string(x = "t", y = paste0("x", j)), color = "black", linewidth = 0.5) +
                # Add shaded rectangles for motif positions
                geom_rect(data = motif_data[[j]], 
                          aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = factor(initial_number)), 
                          alpha = 0.005) +
                # Add SNR text on top of each rectangle
                geom_text(data = curves[[k]]$SNR[[j]], 
                          aes(x = (xmin + xmax) / 2, y = Inf, 
                              label = paste("SNR:", round(SNR, 3))),
                          vjust = 1.5, color = "black", size = 3.5) +
                # Plot motifs with distinct colors
                geom_line(data = motif_data[[j]], aes(x = t, y = x, color = factor(initial_number), group = motif_id), linewidth = 1.0) + 
                # Add color and fill scales with custom labels for motif_id
                scale_color_manual(
                  values = c('background_curve' = scales::alpha('gray30', 0.15), 
                             'with_motif' = 'gold', 
                             motif_colors),
                  labels = c('background_curve' = 'background_curve', 
                             'with_motif' = 'with_motif', 
                             setNames(motif_labels, unique(motif_data[[j]]$initial_number)))
                ) +
                scale_fill_manual(
                  values = motif_colors, 
                  labels = setNames(motif_labels, unique(motif_data[[j]]$initial_number))
                ) +
                # Title
                labs(
                  title = paste0('<b><span style="color:#0073C2;">Random curve ', k, ' - type_error ', 
                                 ifelse(j %% length(curves[[k]]$with_error$error_y) == 0, 
                                        length(curves[[k]]$with_error$error_y), j %% length(curves[[k]]$with_error$error_y)), 
                                 '</span></b>'), 
                  x = "t", 
                  y = paste0("x", j)
                ) +
                # Clean theme with subtle grid
                theme_minimal(base_size = 15) +
                theme(
                  plot.title = element_markdown(size = 20, face = "bold", hjust = 0.5, margin = margin(b = 10)),
                  axis.title = element_text(size = 14, margin = margin(t = 10)),
                  axis.text = element_text(size = 12),
                  legend.text = element_text(size = 12),
                  legend.position = "right",
                  legend.box.margin = margin(10, 10, 10, 10),
                  plot.margin = margin(15, 15, 15, 15),
                  panel.grid.major = element_line(color = "gray90"),
                  panel.grid.minor = element_blank()
                ) +
                guides(
                  color = guide_legend(ncol = 1, byrow = TRUE, title = NULL),
                  fill = "none",
                )
              
              return(pic)
            })
          }
        })
     
        plots <- unlist(plots, recursive = FALSE)
        
        plots_per_page <- n_error
        num_pages <- ceiling(length(plots) / plots_per_page)
        current_page <- reactiveVal(1)
        
        output$plots_ui <- renderUI({
          page <- current_page()
          start <- (page - 1) * plots_per_page + 1
          end <- min(page * plots_per_page, length(plots))
          
          plot_output_list <- lapply(start:end, function(i) {
            plotname <- paste0("plot", i)
            fluidRow(
              column(
                width = 12,
                div(class = "plot-container",
                    div(class = "plot-title"),
                    plotOutput(outputId = plotname, height = "350px", width = "150%")
                )
              )
            )
          })
          
          do.call(tagList, plot_output_list)
        })
        
        lapply(seq_along(plots), function(i) {
          output[[paste0("plot", i)]] <- renderPlot({
            plots[[i]]
          })
        })
        
        observe({
          shinyjs::enable("prevPage")
          shinyjs::enable("nextPage")
          
          if (current_page() == 1) {
            shinyjs::disable("prevPage")
          }
          
          if (current_page() == num_pages) {
            shinyjs::disable("nextPage")
          }
        })
        
        observeEvent(input$prevPage, {
          current_page(max(1, current_page() - 1))
        })
        
        observeEvent(input$nextPage, {
          current_page(min(num_pages, current_page() + 1))
        })
        
        output$downloadPdf <- downloadHandler(
          filename = function() {
            "plots.pdf"
          },
          content = function(file) {
            pdf(file = file, width = 12, height = 14)
            for (plot in plots) {
              print(plot)
            }
            dev.off()
          }
        )
        
        output$downloadRdata <- downloadHandler(
          filename = function() {
            "simulation_results.Rdata"
          },
          content = function(file) {
            save(builder, curves, file = file)
          }
        )
        
        shinyjs::html("error_message", "")
        
      }, error = function(e) {
        shinyjs::html("error_message", paste("Error: ", e$message))
      }, finally = {
        remove_modal_spinner()
      })
    })
  }
  
  shinyApp(ui = ui, server = server)
}