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

ui <- fluidPage(
  titlePanel(
    div(
      style = "text-align: left;",
      tags$h2("Gene Coexpression Network Explorer"),
      tags$h5(
        style = "color:#666; margin-top:-5px;",
        "Powered by Louvain Community Detection"
      )
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
      
      # Status
      h4("Status:"),
      textOutput("status_text"),
      tags$hr(),
      
      # Dataset
      h4("Dataset"),
      selectInput(
        "dataset",
        label   = NULL,
        choices = c("Breast" = "breast", "Ovarian" = "ovarian"),
        selected = "breast"
      ),
      
      tags$hr(),
      
      # Edge sign section
      h4("Edge sign to include"),
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
      checkboxInput(
        "color_nodes_comm",
        "Color nodes by Louvain community",
        value = FALSE   # default OFF to avoid rainbow explosion
      ),
      
      selectInput(
        "community",
        "Focus on one community (optional):",
        choices = "All",
        selected = "All"
      ),
      
      sliderInput(
        "top_hubs",
        "Show only top N hub genes in current view (0 = show all):",
        min = 0, max = 2000, value = 0, step = 50
      ),
      
      checkboxInput("show_labels", "Show gene labels", TRUE),
      
      tags$hr(),
      
      # Edge limit
      h4("Edge limit"),
      sliderInput(
        "max_edges",
        "Maximum edges to plot:",
        min = 1000, max = 50000, value = 20000, step = 1000
      )
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
  
  # Status text reactive
  status <- reactiveVal("Click 'Run Full Analysis in Go' to start.")
  output$status_text <- renderText(status())
  
  # ---------------------- Run Go + load CSVs -----------------------
  
  graph_data <- eventReactive(input$run_go, {
    
    status("Running Go pipeline (this may take ~10–30 seconds)...")
    
    #directory changes per person (depends on where the github repo is downloaded)
    project_dir <- here()
    go_exe <- here("02-601_Project_Fall2025")
    base_dir <- here("ShinyApp")
    
    if (!file.exists(go_exe)) {
      status("Error: Go executable not found.")
      validate(need(FALSE, "Go executable not found."))
    }
    
    # Run Go
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
    breast_nodes_path  <- file.path(base_dir, "breast_nodes_communities.csv")
    breast_edges_path  <- file.path(base_dir, "breast_edges.csv")
    ovarian_nodes_path <- file.path(base_dir, "ovarian_nodes_communities.csv")
    ovarian_edges_path <- file.path(base_dir, "ovarian_edges.csv")
    
    status("Running Go pipeline: loading CSV outputs...")
    
    if (!file.exists(breast_nodes_path) ||
        !file.exists(breast_edges_path) ||
        !file.exists(ovarian_nodes_path) ||
        !file.exists(ovarian_edges_path)) {
      
      status("Error: One or more CSV files not found. Check Go output paths.")
      validate(need(FALSE, "Required CSV files are missing."))
    }
    
    # Read CSVs
    breast_nodes  <- read.csv(breast_nodes_path, stringsAsFactors = FALSE)
    breast_edges  <- read.csv(breast_edges_path, stringsAsFactors = FALSE)
    ovarian_nodes <- read.csv(ovarian_nodes_path, stringsAsFactors = FALSE)
    ovarian_edges <- read.csv(ovarian_edges_path, stringsAsFactors = FALSE)
    
    message("breast_nodes:  ", nrow(breast_nodes),  " rows")
    message("breast_edges:  ", nrow(breast_edges),  " rows")
    message("ovarian_nodes: ", nrow(ovarian_nodes), " rows")
    message("ovarian_edges: ", nrow(ovarian_edges), " rows")
    
    # Make IDs/labels character
    breast_nodes$id    <- as.character(breast_nodes$id)
    breast_nodes$label <- as.character(breast_nodes$label)
    ovarian_nodes$id    <- as.character(ovarian_nodes$id)
    ovarian_nodes$label <- as.character(ovarian_nodes$label)
    
    # Edge endpoints as character
    breast_edges$from <- as.character(breast_edges$from)
    breast_edges$to   <- as.character(breast_edges$to)
    ovarian_edges$from <- as.character(ovarian_edges$from)
    ovarian_edges$to   <- as.character(ovarian_edges$to)
    
    # Check for weight column
    if (!"weight" %in% names(breast_edges) ||
        !"weight" %in% names(ovarian_edges)) {
      status("Error: edges CSV must contain a 'weight' column.")
      validate(need(FALSE, "Missing 'weight' column in edges CSV."))
    }
    
    # Communities for dropdown
    breast_comm  <- sort(unique(breast_nodes$community))
    ovarian_comm <- sort(unique(ovarian_nodes$community))
    
    status("Finished Go pipeline. Ready to render network.")
    
    list(
      breast = list(
        nodes = breast_nodes,
        edges = breast_edges,
        communities = breast_comm
      ),
      ovarian = list(
        nodes = ovarian_nodes,
        edges = ovarian_edges,
        communities = ovarian_comm
      )
    )
  })
  
  # ------------------- Update community dropdown -------------------
  
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
  
  # ------------------------ Render network -------------------------
  
  output$network <- renderVisNetwork({
    dat_all <- graph_data()
    req(dat_all)
    
    status("Rendering network… applying filters and computing layout. This can take several seconds for large graphs.")
    
    dataset <- input$dataset
    dat     <- dat_all[[dataset]]
    req(dat)
    
    nodes <- dat$nodes
    edges <- dat$edges
    
    # 1) classify edge sign
    edges$sign <- ifelse(edges$weight >= 0, "pos", "neg")
    
    # 2) filter by sign
    sel_sign <- c()
    if (isTRUE(input$edge_pos)) sel_sign <- c(sel_sign, "pos")
    if (isTRUE(input$edge_neg)) sel_sign <- c(sel_sign, "neg")
    
    if (length(sel_sign) == 0) {
      status("No edge sign selected. Please choose positive and/or negative edges.")
      validate(need(FALSE, "Select at least one edge sign (positive / negative)."))
    }
    
    edges <- edges[edges$sign %in% sel_sign, , drop = FALSE]
    if (nrow(edges) == 0) {
      status("No edges remain after sign filtering. Try changing your options.")
      validate(need(FALSE, "No edges left after filtering by sign."))
    }
    
    # 3) cap number of edges by strongest |weight|
    if (!is.null(input$max_edges) && nrow(edges) > input$max_edges) {
      edges <- edges[order(-abs(edges$weight)), ][1:input$max_edges, , drop = FALSE]
    }
    
    # 4) community filter
    comm_choice <- input$community
    if (!is.null(comm_choice) && comm_choice != "All") {
      status(paste0("Rendering network for community ", comm_choice, "…"))
      comm_choice_num <- as.numeric(comm_choice)
      nodes <- nodes[nodes$community == comm_choice_num, , drop = FALSE]
      used_ids <- nodes$id
      edges <- edges[edges$from %in% used_ids & edges$to %in% used_ids, , drop = FALSE]
    }
    
    # 5) drop nodes not in any edge
    used_ids <- unique(c(edges$from, edges$to))
    nodes    <- nodes[nodes$id %in% used_ids, , drop = FALSE]
    if (nrow(nodes) == 0 || nrow(edges) == 0) {
      status("No nodes/edges remain after filters. Try relaxing your options.")
      validate(need(FALSE, "No nodes/edges left after filtering."))
    }
    
    # 6) hub filter
    top_n <- input$top_hubs
    if (!is.null(top_n) && top_n > 0) {
      status(paste0("Focusing on top ", top_n, " hubs in current view…"))
      deg      <- table(c(edges$from, edges$to))
      deg      <- sort(deg, decreasing = TRUE)
      keep_ids <- names(deg)[1:min(top_n, length(deg))]
      
      nodes <- nodes[nodes$id %in% keep_ids, , drop = FALSE]
      edges <- edges[edges$from %in% nodes$id & edges$to %in% nodes$id, , drop = FALSE]
      if (nrow(nodes) == 0 || nrow(edges) == 0) {
        status("Top hub filter left no nodes/edges. Try reducing N.")
        validate(need(FALSE, "No nodes/edges left after hub filtering."))
      }
    }
    
    # 7) node colouring
    if (isTRUE(input$color_nodes_comm)) {
      nodes$group <- as.factor(nodes$community)
    } else {
      nodes$group <- NA
    }
    
    # 8) edge colouring by sign
    if (isTRUE(input$color_edges_by_sign)) {
      pos_col <- input$pos_color
      neg_col <- input$neg_color
      edges$color <- ifelse(edges$sign == "pos", pos_col, neg_col)
    }
    
    # 9) hide labels if requested
    if (!isTRUE(input$show_labels)) {
      nodes$label <- NA
    }
    
    # 10) draw with browser physics first, then freeze
    net <- visNetwork(nodes, edges) %>%
      visOptions(
        highlightNearest = TRUE,
        nodesIdSelection = TRUE
      ) %>%
      visPhysics(
        solver = "barnesHut",
        stabilization = list(
          enabled    = TRUE,
          iterations = 1500  # how long it tries to find a nice layout
        )
      ) %>%
      visEvents(
        stabilizationIterationsDone = "function () {
          // stop physics once stabilized so it doesn't keep wiggling
          this.setOptions({physics: false});
        }"
      ) %>%
      visExport(
        type = "png",
        name = paste0(dataset, "_network")
      ) %>%
      visLegend(
        enabled   = TRUE,
        useGroups = isTRUE(input$color_nodes_comm)
      )
    
    status("Ready. Adjust filters to explore different network views.")
    net
  })
}

# -------------------------------------------------------------------
# Run the app
# -------------------------------------------------------------------

shinyApp(ui, server)

