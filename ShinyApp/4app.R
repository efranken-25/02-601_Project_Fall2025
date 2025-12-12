# app.R

# Install packages if needed
if (!requireNamespace("shiny", quietly = TRUE)) {
  install.packages("shiny")
}
if (!requireNamespace("visNetwork", quietly = TRUE)) {
  install.packages("visNetwork")
}
if (!requireNamespace("colourpicker", quietly = TRUE)) {
  install.packages("colourpicker")
}
if (!requireNamespace("shinycssloaders", quietly = TRUE)) {
  install.packages("shinycssloaders")
}
if (!requireNamespace("here", quietly = TRUE)) {
  install.packages("here")
}

library(shiny)
library(visNetwork)
library(colourpicker)
library(shinycssloaders)
library(here)

# -------------------------------------------------------------------
# UI
# -------------------------------------------------------------------

# app.R

# -------------------------------------------------------------------
# Install packages if needed
# -------------------------------------------------------------------
packages <- c("shiny", "visNetwork", "colourpicker", "shinycssloaders", "here")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(shiny)
library(visNetwork)
library(colourpicker)
library(shinycssloaders)
library(here)

# -------------------------------------------------------------------
# UI
# -------------------------------------------------------------------
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      /* Style the visNetwork container + add PNG watermark */
      #network {
        position: relative;
        background-color: rgba(255,255,255,0.9);
        border-radius: 10px;
        box-shadow: 0 2px 6px rgba(0,0,0,0.12);
        background-image: url('watermark.png'); /* Make sure this is in www/ */
        background-repeat: no-repeat;
        background-position: bottom right;
        background-size: 180px auto;
      }
    "))
  ),
  
  # Title + subtitle + centered image
  div(
    style = "text-align: left; margin-left: 20px;",
    tags$h2("Gene Co-expression Network Explorer"),
    tags$h5(
      style = "color:#666; margin-top:-5px;",
      "Powered by Louvain Community Detection"
    ),
    tags$img(
      src = "watermark.png",  # Should be in www/
      style = "
        display: block;
        margin-top: 10px;
        margin-left: 90px;
        width: 150px;
      "
    )
  ),
  
  sidebarLayout(
    sidebarPanel(
      width = 4,
      
      # Run button
      actionButton(
        "run_go",
        "Run Full Analysis in Go",
        class = "btn btn-primary btn-lg btn-block"
      ),
      
      br(),
      h4("Status:"),
      textOutput("status_text"),
      tags$hr(),
      
      # Dataset selection
      h4("Dataset"),
      selectInput(
        "dataset",
        label   = NULL,
        choices = c("Dataset 1" = "dataset1", "Dataset 2" = "dataset2"),
        selected = "dataset1"
      ),
      tags$hr(),
      
      # Edge sign
      h4("Select edges (positive or negative)"),
      checkboxInput("edge_pos", "Positive (r > 0)", value = TRUE),
      checkboxInput("edge_neg", "Negative (r < 0)", value = TRUE),
      tags$hr(),
      
      # Edge colouring
      h4("Edge colouring"),
      checkboxInput("color_edges_by_sign", "Color edges by sign", TRUE),
      div(
        style = "margin-left: 10px;",
        colourInput("pos_color", "Positive edge color", "#1F77B4"),
        colourInput("neg_color", "Negative edge color", "#D62728")
      ),
      tags$hr(),
      
      # Node options
      h4("Node options"),
      checkboxInput("color_nodes_comm", "Color nodes by Louvain community", FALSE),
      tags$hr(),
      
      # Community selection
      h4("Community Selection"),
      tags$small(
        style = "color:#666;",
        "The dropdown lists communities in descending order of density (most dense → least dense)."
      ),
      selectInput(
        "community",
        label = NULL,
        choices = "All",
        selected = "All"
      ),
      
      checkboxInput("show_labels", "Show gene labels", TRUE)
    ),
    
    mainPanel(
      width = 8,
      withSpinner(
        visNetworkOutput("network", height = "700px"),
        type = 4
      )
    )
  )
)

# -------------------------------------------------------------------
# SERVER
# -------------------------------------------------------------------
server <- function(input, output, session) {
  
  status <- reactiveVal("Click 'Run Full Analysis in Go' to start.")
  output$status_text <- renderText(status())
  
  graph_data <- eventReactive(input$run_go, {
    
    status("Running Go pipeline (this may take ~10–30 seconds)...")
    
    go_exe   <- "/Users/emmafranken/02-601_Project_Fall2025/02-601_Project_Fall2025-exe"
    base_dir <- "/Users/emmafranken/02-601_Project_Fall2025//ShinyApp"
    
    if (!file.exists(go_exe)) {
      status("Error: Go executable not found.")
      validate(need(FALSE, "Go executable not found."))
    }
    
    cmd_out <- NULL
    status("Running Go pipeline: executing program...")
    run_err <- tryCatch({
      cmd_out <- system2(go_exe, stdout = TRUE, stderr = TRUE)
      NULL
    }, error = function(e) e)
    
    if (inherits(run_err, "error")) {
      status(paste("Error running Go pipeline:", run_err$message))
      validate(need(FALSE, "Go pipeline failed. See R console for details."))
    }
    
    message("Go output:")
    message(paste(cmd_out, collapse = "\n"))
    
    # CSV paths
    breast_nodes_path  <- file.path(base_dir, "dataset1_nodes_communities.csv")
    breast_edges_path  <- file.path(base_dir, "dataset1_edges.csv")
    ovarian_nodes_path <- file.path(base_dir, "dataset2_nodes_communities.csv")
    ovarian_edges_path <- file.path(base_dir, "dataset2_edges.csv")
    
    status("Loading CSV outputs...")
    
    if (!file.exists(breast_nodes_path) ||
        !file.exists(breast_edges_path) ||
        !file.exists(ovarian_nodes_path) ||
        !file.exists(ovarian_edges_path)) {
      status("Error: One or more CSV files not found.")
      validate(need(FALSE, "Required CSV files are missing."))
    }
    
    # Read CSVs
    breast_nodes  <- read.csv(breast_nodes_path, stringsAsFactors = FALSE)
    breast_edges  <- read.csv(breast_edges_path, stringsAsFactors = FALSE)
    ovarian_nodes <- read.csv(ovarian_nodes_path, stringsAsFactors = FALSE)
    ovarian_edges <- read.csv(ovarian_edges_path, stringsAsFactors = FALSE)
    
    # Ensure IDs/labels are character
    for (df in list(breast_nodes, ovarian_nodes)) {
      df$id <- as.character(df$id)
      df$label <- as.character(df$label)
    }
    
    for (df in list(breast_edges, ovarian_edges)) {
      df$from <- as.character(df$from)
      df$to   <- as.character(df$to)
    }
    
    # Ensure weight column exists
    if (!"weight" %in% names(breast_edges) ||
        !"weight" %in% names(ovarian_edges)) {
      status("Error: edges CSV must contain a 'weight' column.")
      validate(need(FALSE, "Missing 'weight' column in edges CSV."))
    }
    
    list(
      breast = list(nodes = breast_nodes, edges = breast_edges, communities = sort(unique(breast_nodes$community))),
      ovarian = list(nodes = ovarian_nodes, edges = ovarian_edges, communities = sort(unique(ovarian_nodes$community)))
    )
  })
  
  # Update community dropdown
  observeEvent(list(graph_data(), input$dataset), {
    dat_all <- graph_data()
    req(dat_all)
    
    current <- input$dataset
    comms   <- dat_all[[current]]$communities
    choices <- c("All", as.character(comms))
    
    updateSelectInput(
      session,
      "community",
      choices  = choices,
      selected = "All"
    )
  })
  
  # Render network
  output$network <- renderVisNetwork({
    dat_all <- graph_data()
    req(dat_all)
    
    dataset <- input$dataset
    dat     <- dat_all[[dataset]]
    req(dat)
    
    nodes <- dat$nodes
    edges <- dat$edges
    
    # Edge sign
    edges$sign <- ifelse(edges$weight >= 0, "pos", "neg")
    sel_sign <- c()
    if (isTRUE(input$edge_pos)) sel_sign <- c(sel_sign, "pos")
    if (isTRUE(input$edge_neg)) sel_sign <- c(sel_sign, "neg")
    edges <- edges[edges$sign %in% sel_sign, , drop = FALSE]
    
    # Node coloring
    if (isTRUE(input$color_nodes_comm)) nodes$group <- as.factor(nodes$community)
    
    # Edge coloring
    if (isTRUE(input$color_edges_by_sign)) {
      edges$color <- ifelse(edges$sign == "pos", input$pos_color, input$neg_color)
    }
    
    if (!isTRUE(input$show_labels)) nodes$label <- NA
    
    visNetwork(nodes, edges) %>%
      visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
      visPhysics(solver = "barnesHut", stabilization = list(enabled = TRUE, iterations = 1500)) %>%
      visEvents(stabilizationIterationsDone = "function () { this.setOptions({physics: false}); }") %>%
      visExport(type = "png", name = paste0(dataset, "_network")) %>%
      visLegend(enabled = TRUE, useGroups = isTRUE(input$color_nodes_comm))
  })
}

# -------------------------------------------------------------------
# Run the app
# -------------------------------------------------------------------
shinyApp(ui, server)
