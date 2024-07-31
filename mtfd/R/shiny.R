library(shiny)
library(ggplot2)
library(dplyr)
library(shinyjs)

ui <- fluidPage(
  useShinyjs(),  # Initialize shinyjs
  
  tags$head(
    tags$style(HTML("
      body {
        background-color: #001f3f;
        color: white;
      }
      .well {
        background-color: #003366;
        color: white;
      }
      #pagination_controls {
        display: flex;
        justify-content: center;
        margin-top: 20px;
      }
      .pagination-buttons {
        margin: 0 10px;
      }
      .collapsible {
        cursor: pointer;
        padding: 10px;
        background-color: #003366;
        border: none;
        color: white;
        text-align: left;
        outline: none;
        font-size: 16px;
        width: 100%;
        display: block;
        margin-bottom: 5px;
      }
      .active, .collapsible:hover {
        background-color: #00509e;
      }
      .content {
        display: none;
        padding: 10px;
        background-color: #002244;
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
      }
      .input-group {
        flex: 1;
        margin-right: 10px;
      }
      .input-group:last-child {
        margin-right: 0;
      }
    ")),
    tags$script(HTML("
      $(document).on('click', '.collapsible', function() {
        $(this).toggleClass('active');
        $(this).next('.content').slideToggle('fast');
      });
    "))
  ),
  
  titlePanel("Motif Simulation Plotter"),
  
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
    downloadButton("downloadPdf", "Download PDF")
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
                 div(class = "input-group", selectInput("distribution", "Coefficient Distribution", choices = c("unif", "norm"), selected = "unif"))
             ),
             numericInput("numMotifs", "Number of Motifs", value = 5, min = 1),
             uiOutput("motifInputs"),
             textInput("path", "Output Directory", value = tempdir())
           )
    )
  )
)

server <- function(input, output, session) {
  
  observeEvent(input$numMotifs, {
    output$motifInputs <- renderUI({
      numMotifs <- input$numMotifs
      lapply(1:numMotifs, function(i) {
        tagList(
          actionButton(paste0("toggleMotif", i), paste("Motif", i), class = "collapsible"),
          div(id = paste0("motif", i, "Content"), class = "content",
              numericInput(paste0("mot", i, "_len"), paste("Motif", i, "Length"), value = 100, min = 1),
              textInput(paste0("mot", i, "_weights"), paste("Motif", i, "Weights (comma-separated)"), value = ""),
              numericInput(paste0("mot", i, "_appearance"), paste("Motif", i, "Appearance"), value = 10, min = 0)
          )
        )
      })
    })
  })
  
  observeEvent(input$plotBtn, {
    req(input$path)
    
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
    
    mot_details <- lapply(1:input$numMotifs, function(i) {
      list(
        len = input[[paste0("mot", i, "_len")]],
        weights = parse_weights(input[[paste0("mot", i, "_weights")]]),
        appearance = input[[paste0("mot", i, "_appearance")]]
      )
    })
    
    distribution <- input$distribution
    path <- input$path
    
    # Run the motifSimulationBuilder algorithm
    builder <- motifSimulationBuilder(curve_details, mot_details, distribution)
    
    curves <- generateCurves(builder,0.1)
    
    output_file <- file.path(path, "plots.pdf")
    
    # Create the directory if it does not exist
    if (!dir.exists(path)) {
      dir.create(path)
    }
    
    # Generate the plots
    plots <- lapply(seq_along(curves), function(i) {
      curve_data <- data.frame(
        t = seq(0, curves[[i]]$basis$rangeval[2]),
        x = eval.fd(seq(0, curves[[i]]$basis$rangeval[2]), curves[[i]], Lfdobj = 0)
      )
      names(curve_data) <- c("t","x")
      
      motif_lines <- mapply(function(id_motif, pos_motif, instance) {
        motif_t = seq((pos_motif - 1) * builder@dist_knots,
                      (pos_motif - 1 + length(curves$coeff_motifs[[id_motif]]) - builder@norder + 1) * builder@dist_knots)
        motif_x = eval.fd(motif_t, curves[[i]], Lfdobj = 0)
        data.frame(t = motif_t, x = motif_x, motif_id = factor(paste(id_motif, instance, sep = "_")),
                   initial_number = str_extract(as.character(id_motif), "^[^_]+"),
                   xmin = (pos_motif - 1) * builder@dist_knots,
                   xmax = (pos_motif - 1 + length(curves$coeff_motifs[[id_motif]]) - builder@norder + 1) * builder@dist_knots)
      }, builder@motifs_in_curves[[i]]$motif_id, builder@motifs_in_curves[[i]]$starting_coeff_pos, seq_along(builder@motifs_in_curves[[i]]$motif_id), SIMPLIFY = FALSE)
      
      motif_data <- bind_rows(motif_lines)
      names(motif_data) <- c("t","x","motif_id","initial_number","xmin","xmax")
      
      motif_colors <- c( "1" = "red", "2" = "green", "3" = "blue", "4" = "orange",
                         "5" = "purple", "6" = "cyan", "7" = "magenta", "8" = "brown",
                         "9" = "pink", "10" = "grey")
      motif_colors <- rep(motif_colors,length.out = length(mot_details))
      if(length(mot_details) > 10 )
        attr(motif_colors,"names")[11:length(mot_details)] <- as.character(as.integer(attr(motif_colors,"names")[11:length(mot_details)]) + 10)
      
      if (nrow(motif_data) > 0) {
        p <- ggplot() +
          # Plot the main curve in black
          geom_line(data = curve_data, aes(x = t, y = x), color = "black", size = 0.5) +
          # Add shaded rectangles for motif positions with transparency
          geom_rect(data = motif_data, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = initial_number), alpha = 0.01) +
          # Plot motifs with distinct colors
          geom_line(data = motif_data, aes(x = t, y = x, color = factor(initial_number), group = motif_id), size = 1.5) + 
          scale_color_manual(values = motif_colors) +
          scale_fill_manual(values = motif_colors) +
          labs(title = paste('Random curve', i), x = 't', y = 'x(t)') +
          theme_minimal(base_size = 15) +
          ylim(-20, 20) +
          guides(color = guide_legend(title = "Motif ID"), fill = guide_legend(title = "Motif ID"))
      }
      p
    })
    
    # Pagination settings
    plots_per_page <- 2
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
  })
}

# Run the app
shinyApp(ui = ui, server = server)