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
    
    tags$head(
      tags$style(HTML("
    /* Your custom CSS */
    body {
      background-color: #001f3f;
      color: #ffffff;
      font-family: 'Arial', sans-serif;
    }
    .well {
      background-color: #003366;
      color: white;
      border: none;
      box-shadow: 0 4px 8px rgba(0, 0, 0, 0.2);
      border-radius: 10px;
      padding: 20px;
    }
    .btn-primary {
      background-color: #007bff;
      border: none;
      padding: 10px 20px;  /* Adjusted padding for larger clickable area */
      font-size: 18px;     /* Larger font size */
      cursor: pointer;
      border-radius: 8px;  /* Rounded corners */
      transition: background-color 0.3s, transform 0.2s; /* Smooth transitions */
      margin-right: 10px;  /* Space between buttons */
    }
    .btn-primary:hover {
      background-color: #0056b3;
      transform: scale(1.05);  /* Slightly increase button size on hover */
    }
    .btn-primary:active {
      background-color: #004080;
      transform: scale(1.02);  /* Slightly reduce button size on click */
    }
    #pagination_controls {
      display: flex;
      justify-content: center;
      margin-top: 20px;
      margin-bottom: 20px; /* Added bottom margin to separate from other elements */
    }
    .pagination-buttons {
      background-color: #007bff;
      border: none;
      color: white;
      padding: 10px 20px;  /* Adjusted padding for larger clickable area */
      font-size: 18px;     /* Larger font size */
      border-radius: 8px;  /* Rounded corners */
      cursor: pointer;
      transition: background-color 0.3s, transform 0.2s; /* Smooth transitions */
      margin: 0 10px;  /* Space between pagination buttons */
    }
    .pagination-buttons:disabled {
      background-color: #cccccc;
      cursor: not-allowed;
    }
    .pagination-buttons:hover:not(:disabled) {
      background-color: #0056b3;
      transform: scale(1.05);  /* Slightly increase button size on hover */
    }
    .pagination-buttons:active:not(:disabled) {
      background-color: #004080;
      transform: scale(1.02);  /* Slightly reduce button size on click */
    }
    .collapsible {
      cursor: pointer;
      padding: 12px;
      background-color: #003366;
      border: none;
      color: white;
      text-align: left;
      outline: none;
      font-size: 18px;  /* Larger font size */
      width: 100%;
      display: block;
      margin-bottom: 5px;
      border-radius: 8px;  /* Rounded corners */
      transition: background-color 0.3s; /* Smooth transitions */
    }
    .active, .collapsible:hover {
      background-color: #00509e;
    }
    .content {
      display: none;
      padding: 10px;
      background-color: #002244;
      border-radius: 5px;
      margin-bottom: 10px;
    }
    .plot-container {
      margin-bottom: 20px;
    }
    .plot-title {
      text-align: center;
      margin-bottom: 10px;
    }
    .input-row {
      display: flex;
      justify-content: space-between;
      flex-wrap: wrap;
    }
    .input-group {
      flex: 1;
      margin-right: 10px;
    }
    .input-group:last-child {
      margin-right: 0;
    }
    #error_message {
      color: red;
      font-weight: bold;
      margin-top: 20px;
    }
  "))
    ),
    
    titlePanel("Comprehensive Curve Plots for Each Error"),
    
    mainPanel(
      uiOutput("plots_ui"),
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
      downloadButton("downloadPdf", "Download PDF"),
      downloadButton("downloadRdata", "Download Rdata"),
      div(id = "error_message", "")
    ),
    
    fluidRow(
      column(width = 12,
             wellPanel(
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
                   div(class = "input-group", numericInput("min_dist_motifs", "Minimum Distance Between Motifs (min_dist_motifs)", value = 30, min = 1)),
                   div(class = "input-group", selectInput("distribution", "Coefficient Distribution", choices = c("unif", "beta"), selected = "unif"))
               ),
               textInput("path", "Output Directory", value = getwd())
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
        curve_details <- list(
          N = input$N,
          len = input$len,
          norder = input$norder,
          coeff_min = input$coeff_min,
          coeff_max = input$coeff_max,
          dist_knots = input$dist_knots,
          min_dist_motifs = input$min_dist_motifs
        )
        
        parse_weights <- function(weights_input) {
          if (weights_input == "") {
            return(NULL)
          }
          as.numeric(unlist(strsplit(weights_input, ",")))
        }
        
        distribution <- input$distribution
        if (!dir.exists(input$path)) {
          stop("Invalid directory path.")
        }
        builder <- mtfd::motifSimulationBuilder(curve_details, mot_details, distribution)
        curves <- mtfd::generateCurves(builder, error_str)
        
        output_file <- file.path(input$path, "plots.pdf")
        
        if (!dir.exists(input$path)) {
          dir.create(input$path)
        }
        
        plots <- lapply(seq_along(curves), function(k) {
          if (!is.null(curves[[k]]$with_error$error_y)) {
            curve_data_no_error <- data.frame(
              t = seq(0, curves[[k]]$no_error$basis$rangeval[2] - 1),
              x = curves[[k]]$no_error$no_error_y)
            names(curve_data_no_error) <- c("t", "x")
            
            curve_data_error <- NULL
            curve_data_error <- data.frame(
              t = seq(0, curves[[k]]$no_error$basis$rangeval[2] - 1),
              x = curves[[k]]$with_error$error_y,
              SNR = unlist(rep(curves[[k]]$SNR[1:length(curves[[k]]$SNR)], len = curves[[k]]$no_error$basis$rangeval[2])))
            names(curve_data_error) <- c("t", paste0("x", seq(length(curves[[k]]$with_error$error_y))), "SNR")
            
            motif_lines <- mapply(function(id_motif, pos_motif, instance) {
              motif_t = seq((pos_motif - 1) * builder@dist_knots,
                            (pos_motif - 1) * builder@dist_knots + builder@mot_details[[id_motif]]$len)
              motif_x = lapply(curves[[k]]$with_error$error_y, function(curve) { return(curve[motif_t + 1]) })
              
              return(lapply(motif_x, function(motif) { data.frame(t = motif_t, x = motif, motif_id = factor(paste(id_motif, instance, sep = "_")),
                                                                  initial_number = str_extract(as.character(id_motif), "^[^_]+"),
                                                                  xmin = (pos_motif - 1) * builder@dist_knots,
                                                                  xmax = (pos_motif - 1) * builder@dist_knots + builder@mot_details[[id_motif]]$len) }))
            }, builder@motifs_in_curves[[k]]$motif_id, builder@motifs_in_curves[[k]]$starting_coeff_pos, seq_along(builder@motifs_in_curves[[k]]$motif_id), SIMPLIFY = FALSE)
            
            motif_colors <- c("1" = "red", "2" = "blue", "3" = "darkgreen", "4" = "orange",
                              "5" = "purple", "6" = "cyan", "7" = "magenta", "8" = "brown",
                              "9" = "pink", "10" = "grey")
            motif_colors <- rep(motif_colors, length.out = length(builder@mot_details))
            if (length(builder@mot_details) > 10)
              attr(motif_colors, "names")[11:length(builder@mot_details)] <- as.character(as.integer(attr(motif_colors, "names")[11:length(builder@mot_details)]) + 10)
            
            max_dataframes <- max(sapply(motif_lines, function(sublist) length(sublist)))
            # Inizializza una lista per memorizzare i risultati
            motif_data <- vector("list", max_dataframes)
            # Cicla su ogni "livello" dei data frame
            for (i in seq_len(max_dataframes)) {
              # Estrai i data frame dal livello i-esimo
              dataframes_at_level_i <- lapply(motif_lines, function(sublist) {
                sublist[[i]]
              })
              # Fai il bind_rows sui data frame estratti
              motif_data[[i]] <- bind_rows(dataframes_at_level_i)
              names(motif_data[[i]]) <- c("t", "x", "motif_id", "initial_number", "xmin", "xmax")
            }
            return(lapply(1:length(motif_data), function(j) {
              pic <- ggplot() +
                # Plot the main curve in black
                geom_line(data = curve_data_no_error, aes(x = t, y = x), color = scales::alpha('gray30', 0.15), linewidth = 0.5) + 
                # Plot the error curve in black
                geom_line(data = curve_data_error, aes_string(x = "t", y = paste0("x", j)), color = "black", linewidth = 0.5) +
                # Add shaded rectangles for motif positions with transparency
                geom_rect(data = motif_data[[j]], aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = initial_number), alpha = 0.005) +
                # Plot motifs with distinct colors
                geom_line(data = motif_data[[j]], aes(x = t, y = x, color = factor(initial_number), group = motif_id), linewidth = 1.0) + 
                scale_color_manual(values = motif_colors) +
                scale_fill_manual(values = motif_colors) +
                # Add SNR to the title
                labs(
                  title = paste('Random curve', k, '- SNR:', round(curve_data_error$SNR[j], 3)), 
                  x = "t", 
                  y = paste0("x", j)
                ) +
                theme_minimal(base_size = 15) +
                guides(color = guide_legend(title = "Motif ID"), fill = guide_legend(title = "Motif ID"))
              
              return(pic)
            }))
          }
        })
        plots <- unlist(plots, recursive = FALSE)
        
        plots_per_page <- dim(error_str)[1]
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