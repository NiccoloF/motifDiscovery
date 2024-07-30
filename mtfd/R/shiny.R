library(shiny)
library(ggplot2)
library(dplyr)

ui <- fluidPage(
  useShinyjs(),
  titlePanel("Motif Simulation Plotter"),
  sidebarLayout(
    sidebarPanel(
      numericInput("N", "Number of Curves (N)", value = 20),
      numericInput("len", "Curve Length (len)", value = 300),
      numericInput("norder", "Order of Basis (norder)", value = 3),
      numericInput("coeff_min", "Minimum Coefficient (coeff_min)", value = -15),
      numericInput("coeff_max", "Maximum Coefficient (coeff_max)", value = 15),
      numericInput("dist_knots", "Distance Between Knots (dist_knots)", value = 10),
      numericInput("min_dist_motifs", "Minimum Distance Between Motifs (min_dist_motifs)", value = 30),
      selectInput("distribution", "Coefficient Distribution", choices = c("unif", "norm"), selected = "unif"),
      
      h4("Motif 1"),
      numericInput("mot1_len", "Motif 1 Length", value = 100),
      textInput("mot1_weights", "Motif 1 Weights (comma-separated)", value = ""),
      numericInput("mot1_appearance", "Motif 1 Appearance", value = 10),
      
      h4("Motif 2"),
      numericInput("mot2_len", "Motif 2 Length", value = 100),
      textInput("mot2_weights", "Motif 2 Weights (comma-separated)", value = ""),
      numericInput("mot2_appearance", "Motif 2 Appearance", value = 10),
      
      h4("Motif 3"),
      numericInput("mot3_len", "Motif 3 Length", value = 100),
      textInput("mot3_weights", "Motif 3 Weights (comma-separated)", value = ""),
      numericInput("mot3_appearance", "Motif 3 Appearance", value = 10),
      
      h4("Motif 4"),
      numericInput("mot4_len", "Motif 4 Length", value = 100),
      textInput("mot4_weights", "Motif 4 Weights (comma-separated)", value = ""),
      numericInput("mot4_appearance", "Motif 4 Appearance", value = 10),
      
      h4("Motif 5"),
      numericInput("mot5_len", "Motif 5 Length", value = 100),
      textInput("mot5_weights", "Motif 5 Weights (comma-separated)", value = ""),
      numericInput("mot5_appearance", "Motif 5 Appearance", value = 10),
      
      textInput("path", "Output Directory", value = tempdir()),
      actionButton("plotBtn", "Generate Plots"),
      div(
        id = "pagination_controls",
        actionButton("prevPage", "Previous", disabled = TRUE),
        actionButton("nextPage", "Next", disabled = TRUE)
      )
    ),
    mainPanel(
      uiOutput("plots_ui"),
      downloadButton("downloadPdf", "Download PDF")
    )
  )
)

server <- function(input, output, session) {
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
    
    mot1 <- list(
      len = input$mot1_len,
      weights = parse_weights(input$mot1_weights),
      appearance = input$mot1_appearance
    )
    
    mot2 <- list(
      len = input$mot2_len,
      weights = parse_weights(input$mot2_weights),
      appearance = input$mot2_appearance
    )
    
    mot3 <- list(
      len = input$mot3_len,
      weights = parse_weights(input$mot3_weights),
      appearance = input$mot3_appearance
    )
    
    mot4 <- list(
      len = input$mot4_len,
      weights = parse_weights(input$mot4_weights),
      appearance = input$mot4_appearance
    )
    
    mot5 <- list(
      len = input$mot5_len,
      weights = parse_weights(input$mot5_weights),
      appearance = input$mot5_appearance
    )
    mot_details <- list(mot1,mot2)
    #mot_details <- list(mot1, mot2, mot3, mot4, mot5)
    #mot_details <- mot_details[sapply(mot_details, function(mot) !is.null(mot$weights) || mot$appearance != 0)]
    
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
      
      motif_colors <- c("1" = "red", "2" = "green", "3" = "blue", "4" = "orange", "5" = "purple", "6" = "cyan", "7" = "magenta")
      
      if (nrow(motif_data) > 0) {
        p <- ggplot() +
          geom_line(data = curve_data, aes(x = t, y = x), color = "black", size = 0.5) +
          geom_rect(data = motif_data, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = initial_number), alpha = 0.01) +
          geom_line(data = motif_data, aes(x = t, y = x, color = factor(initial_number, levels = 1:7), group = motif_id), size = 1.5) + 
          scale_color_manual(values = motif_colors) +
          scale_fill_manual(values = motif_colors) +
          labs(title = paste('Random curve', i), x = 't', y = 'x(t)') +
          theme_minimal(base_size = 15) +
          ylim(-20, 20) +
          guides(color = guide_legend(title = "Motif ID"), fill = guide_legend(title = "Motif ID"))
        return(p)
      }
    })
    
    plots <- Filter(Negate(is.null), plots)
    
    # Pagination settings
    plots_per_page <- 3
    num_pages <- ceiling(length(plots) / plots_per_page)
    current_page <- reactiveVal(1)
    
    output$plots_ui <- renderUI({
      page <- current_page()
      start <- (page - 1) * plots_per_page + 1
      end <- min(page * plots_per_page, length(plots))
      
      plot_output_list <- lapply(start:end, function(i) {
        plotname <- paste0("plot", i)
        plotOutput(outputId = plotname, height = "600px")
      })
      
      tagList(plot_output_list)
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
        pdf(file = file, width = 8, height = 6)
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
