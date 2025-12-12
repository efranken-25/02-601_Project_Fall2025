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
  tags$head(
    tags$style(HTML("
      /* Style the visNetwork container + add PNG watermark */
      #network {
        position: relative;
        background-color: rgba(255,255,255,0.9);
        border-radius: 10px;
        box-shadow: 0 2px 6px rgba(0,0,0,0.12);

        /* PNG watermark inside the panel */
        background-image: url('watermark.png');
        background-repeat: no-repeat;
        background-position: bottom right;   /* bottom-right corner */
        background-size: 180px auto;         /* tweak this size */
      }
    "))
  ),
  
  # Title + subtitle + centered-under-title image
  div(
    style = "text-align: left; margin-left: 20px;",
    
    tags$h2("Gene Co-expression Network Explorer"),
    
    tags$h5(
      style = "color:#666; margin-top:-5px;",
      "Powered by Louvain Community Detection"
    ),
    
    tags$img(
      src = "watermark.png",
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
      
      # Status
      h4("Status:"),
      textOutput("status_text"),
      tags$hr(),
      
      # Dataset
      h4("Dataset"),
      selectInput(
        "dataset",
        label   = NULL,
        choices = c("Dataset 1" = "dataset1", "Dataset 2" = "dataset2"),
        selected = "dataset1",
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
      checkboxInput(
        "color_nodes_comm",
        "Color nodes by Louvain community",
        value = FALSE
      ),
      
      tags$hr(),
      
      h4("Community Selection"),
      tags$small(
        style = "color:#666;",
        "The dropdown lists communities in descending order of density (most dense → least dense)."
      ),
      
      selectInput(
        "community",
        label = NULL,
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
  
  # Status text reactive
  status <- reactiveVal("Click 'Run Full Analysis in Go' to start.")
  output$status_text <- renderText(status())
  
  # ---------------------- Run Go + load CSVs -----------------------
  
  graph_data <- eventReactive(input$run_go, {
    
    status("Running Go pipeline (this may take ~10–30 seconds)...")
    
    # here() points to ShinyApp (where app.R lives)
    base_dir <- here("ShinyApp")  # /Users/bethany/ProgrammingProject/ShinyApp
    
    # Project root (contains ProgrammingProject binary + ShinyApp/)
    project_dir <- dirname(base_dir)  
    
    # Go executable in the project root
    go_exe <- file.path(project_dir, "02-601_Project_Fall2025")
    
    message("DEBUG: project_dir = ", project_dir)
    message("DEBUG: go_exe      = ", go_exe)
    message("DEBUG: base_dir    = ", base_dir)
    message("DEBUG: getwd() BEFORE running Go = ", getwd())
    
    if (!file.exists(go_exe) || file.info(go_exe)$isdir) {
      status(paste("Error: Go executable not found at:", go_exe))
      validate(need(FALSE, paste("Go executable not found at:", go_exe)))
    if (!file.exists(go_exe) || file.info(go_exe)$isdir) {
      status(paste("Error: Go executable not found at:", go_exe))
      validate(need(FALSE, paste("Go executable not found at:", go_exe)))
    }
    
    # Run Go from the project root so its relative paths work
    # Run Go from the project root so its relative paths work
    status("Running Go pipeline: executing program...")
    old_wd <- getwd()
    cmd_out <- NULL
    
    old_wd <- getwd()
    cmd_out <- NULL
    
    run_err <- tryCatch({
      setwd(project_dir)
      message("DEBUG: getwd() INSIDE Go run = ", getwd())
      cmd_out <- system2(
        go_exe,
        stdout = TRUE,
        stderr = TRUE
      )
      setwd(project_dir)
      message("DEBUG: getwd() INSIDE Go run = ", getwd())
      cmd_out <- system2(
        go_exe,
        stdout = TRUE,
        stderr = TRUE
      )
      NULL
    }, error = function(e) e)
    
    setwd(old_wd)
    message("DEBUG: getwd() AFTER running Go = ", getwd())
    
    setwd(old_wd)
    message("DEBUG: getwd() AFTER running Go = ", getwd())
    
    if (inherits(run_err, "error")) {
      status(paste("Error running Go pipeline:", run_err$message))
      validate(need(FALSE, "Go pipeline failed. See R console for details."))
    }
    
    message("===== RAW GO OUTPUT =====")
    if (length(cmd_out) > 0) {
      message(paste(cmd_out, collapse = "\n"))
    } else {
      message("(no stdout/stderr from Go)")
    }
    message("===== END GO OUTPUT =====")
    
    # Check Go exit status if present
    exit_status <- attr(cmd_out, "status")
    if (!is.null(exit_status) && exit_status != 0) {
      status(paste("Go pipeline failed with exit status", exit_status, "- see R console for details."))
      validate(need(FALSE, paste("Go pipeline failed (status", exit_status, "). Check R console for Go error output.")))
    }
    
    # Before checking files, show what's in ShinyApp
    message("===== FILES IN base_dir (ShinyApp) =====")
    message(paste(list.files(base_dir), collapse = ", "))
    message("===== END FILE LIST =====")
    
    # CSV paths – must match Go's dataset1/dataset2 filenames
    dataset1_nodes_path <- file.path(base_dir, "dataset1_nodes_communities.csv")
    dataset1_edges_path <- file.path(base_dir, "dataset1_edges.csv")
    dataset2_nodes_path <- file.path(base_dir, "dataset2_nodes_communities.csv")
    dataset2_edges_path <- file.path(base_dir, "dataset2_edges.csv")
    message("===== RAW GO OUTPUT =====")
    if (length(cmd_out) > 0) {
      message(paste(cmd_out, collapse = "\n"))
    } else {
      message("(no stdout/stderr from Go)")
    }
    message("===== END GO OUTPUT =====")
    
    # Check Go exit status if present
    exit_status <- attr(cmd_out, "status")
    if (!is.null(exit_status) && exit_status != 0) {
      status(paste("Go pipeline failed with exit status", exit_status, "- see R console for details."))
      validate(need(FALSE, paste("Go pipeline failed (status", exit_status, "). Check R console for Go error output.")))
    }
    
    # Before checking files, show what's in ShinyApp
    message("===== FILES IN base_dir (ShinyApp) =====")
    message(paste(list.files(base_dir), collapse = ", "))
    message("===== END FILE LIST =====")
    
    # CSV paths – must match Go's dataset1/dataset2 filenames
    dataset1_nodes_path <- file.path(base_dir, "dataset1_nodes_communities.csv")
    dataset1_edges_path <- file.path(base_dir, "dataset1_edges.csv")
    dataset2_nodes_path <- file.path(base_dir, "dataset2_nodes_communities.csv")
    dataset2_edges_path <- file.path(base_dir, "dataset2_edges.csv")
    
    status("Running Go pipeline: loading CSV outputs...")
    
    if (!file.exists(dataset1_nodes_path) ||
        !file.exists(dataset1_edges_path) ||
        !file.exists(dataset2_nodes_path) ||
        !file.exists(dataset2_edges_path)) {
    if (!file.exists(dataset1_nodes_path) ||
        !file.exists(dataset1_edges_path) ||
        !file.exists(dataset2_nodes_path) ||
        !file.exists(dataset2_edges_path)) {
      
      status("Error: One or more CSV files not found. Check Go output paths and filenames.")
      validate(need(FALSE, "Required CSV files are missing. See R console for file listing."))
      status("Error: One or more CSV files not found. Check Go output paths and filenames.")
      validate(need(FALSE, "Required CSV files are missing. See R console for file listing."))
    }
    
    # Read CSVs
    dataset1_nodes <- read.csv(dataset1_nodes_path, stringsAsFactors = FALSE)
    dataset1_edges <- read.csv(dataset1_edges_path, stringsAsFactors = FALSE)
    dataset2_nodes <- read.csv(dataset2_nodes_path, stringsAsFactors = FALSE)
    dataset2_edges <- read.csv(dataset2_edges_path, stringsAsFactors = FALSE)
    dataset1_nodes <- read.csv(dataset1_nodes_path, stringsAsFactors = FALSE)
    dataset1_edges <- read.csv(dataset1_edges_path, stringsAsFactors = FALSE)
    dataset2_nodes <- read.csv(dataset2_nodes_path, stringsAsFactors = FALSE)
    dataset2_edges <- read.csv(dataset2_edges_path, stringsAsFactors = FALSE)
    
    message("dataset1_nodes: ", nrow(dataset1_nodes), " rows")
    message("dataset1_edges: ", nrow(dataset1_edges), " rows")
    message("dataset2_nodes: ", nrow(dataset2_nodes), " rows")
    message("dataset2_edges: ", nrow(dataset2_edges), " rows")
    
    dataset1_stats_path <- file.path(base_dir, "dataset1_community_stats.csv")
    dataset2_stats_path <- file.path(base_dir, "dataset2_community_stats.csv")
    
    if (!file.exists(dataset1_stats_path) || !file.exists(dataset2_stats_path)) {
      status("Error: Community stats CSV files not found. Check Go output paths.")
      validate(need(FALSE, "Missing community stats CSVs. See R console for file listing."))
    }
    
    dataset1_stats <- read.csv(dataset1_stats_path, stringsAsFactors = FALSE)
    dataset2_stats <- read.csv(dataset2_stats_path, stringsAsFactors = FALSE)
    message("dataset1_nodes: ", nrow(dataset1_nodes), " rows")
    message("dataset1_edges: ", nrow(dataset1_edges), " rows")
    message("dataset2_nodes: ", nrow(dataset2_nodes), " rows")
    message("dataset2_edges: ", nrow(dataset2_edges), " rows")
    
    dataset1_stats_path <- file.path(base_dir, "dataset1_community_stats.csv")
    dataset2_stats_path <- file.path(base_dir, "dataset2_community_stats.csv")
    
    if (!file.exists(dataset1_stats_path) || !file.exists(dataset2_stats_path)) {
      status("Error: Community stats CSV files not found. Check Go output paths.")
      validate(need(FALSE, "Missing community stats CSVs. See R console for file listing."))
    }
    
    dataset1_stats <- read.csv(dataset1_stats_path, stringsAsFactors = FALSE)
    dataset2_stats <- read.csv(dataset2_stats_path, stringsAsFactors = FALSE)
    
    # Make IDs/labels character
    dataset1_nodes$id    <- as.character(dataset1_nodes$id)
    dataset1_nodes$label <- as.character(dataset1_nodes$label)
    dataset2_nodes$id    <- as.character(dataset2_nodes$id)
    dataset2_nodes$label <- as.character(dataset2_nodes$label)
    dataset1_nodes$id    <- as.character(dataset1_nodes$id)
    dataset1_nodes$label <- as.character(dataset1_nodes$label)
    dataset2_nodes$id    <- as.character(dataset2_nodes$id)
    dataset2_nodes$label <- as.character(dataset2_nodes$label)
    
    # Edge endpoints as character
    dataset1_edges$from <- as.character(dataset1_edges$from)
    dataset1_edges$to   <- as.character(dataset1_edges$to)
    dataset2_edges$from <- as.character(dataset2_edges$from)
    dataset2_edges$to   <- as.character(dataset2_edges$to)
    dataset1_edges$from <- as.character(dataset1_edges$from)
    dataset1_edges$to   <- as.character(dataset1_edges$to)
    dataset2_edges$from <- as.character(dataset2_edges$from)
    dataset2_edges$to   <- as.character(dataset2_edges$to)
    
    # Check for weight column
    if (!"weight" %in% names(dataset1_edges) ||
        !"weight" %in% names(dataset2_edges)) {
    if (!"weight" %in% names(dataset1_edges) ||
        !"weight" %in% names(dataset2_edges)) {
      status("Error: edges CSV must contain a 'weight' column.")
      validate(need(FALSE, "Missing 'weight' column in edges CSV."))
    }
    
    # Communities for dropdown
    dataset1_comm <- sort(unique(dataset1_nodes$community))
    dataset2_comm <- sort(unique(dataset2_nodes$community))
    dataset1_comm <- sort(unique(dataset1_nodes$community))
    dataset2_comm <- sort(unique(dataset2_nodes$community))
    
    status("Finished Go pipeline. Ready to render network.")
    
    list(
      dataset1 = list(
        nodes       = dataset1_nodes,
        edges       = dataset1_edges,
        communities = dataset1_comm,
        comm_stats  = dataset1_stats
      dataset1 = list(
        nodes       = dataset1_nodes,
        edges       = dataset1_edges,
        communities = dataset1_comm,
        comm_stats  = dataset1_stats
      ),
       dataset2 = list(
         nodes       = dataset2_nodes,
         edges       = dataset2_edges,
         communities = dataset2_comm,
         comm_stats  = dataset2_stats
       )
      )
  })
  
  # ------------------- Update community dropdown -------------------
  
  observeEvent(list(graph_data(), input$dataset), {
    dat_all <- graph_data()
    req(dat_all)
    
    current <- input$dataset   # "dataset1" or "dataset2"
    dat     <- dat_all[[current]]
    nodes   <- dat$nodes
    stats   <- dat$comm_stats
    
    # possible communities from the nodes
    comms <- sort(unique(nodes$community))
    
    # subset stats to only those communities present
    stats <- stats[stats$community_id %in% comms, , drop = FALSE]
    
    if (nrow(stats) == 0) {
      # fallback: no stats, just "All"
      choices <- c("All" = "All")
    } else {
      # sort communities by density (highest first) so densest are at the top
      stats <- stats[order(-stats$density), , drop = FALSE]
      
      # values are community IDs (as strings); labels show ID, size, density
      values <- as.character(stats$community_id)
      labels <- sprintf(
        "Community %d (n = %d, dens = %.3f)",
        stats$community_id,
        stats$num_nodes,
        stats$density
      )
      named_choices <- stats::setNames(values, labels)
      
      choices <- c("All" = "All", named_choices)
    }
    current <- input$dataset   # "dataset1" or "dataset2"
    dat     <- dat_all[[current]]
    nodes   <- dat$nodes
    stats   <- dat$comm_stats
    
    # possible communities from the nodes
    comms <- sort(unique(nodes$community))
    
    # subset stats to only those communities present
    stats <- stats[stats$community_id %in% comms, , drop = FALSE]
    
    if (nrow(stats) == 0) {
      # fallback: no stats, just "All"
      choices <- c("All" = "All")
    } else {
      # sort communities by density (highest first) so densest are at the top
      stats <- stats[order(-stats$density), , drop = FALSE]
      
      # values are community IDs (as strings); labels show ID, size, density
      values <- as.character(stats$community_id)
      labels <- sprintf(
        "Community %d (n = %d, dens = %.3f)",
        stats$community_id,
        stats$num_nodes,
        stats$density
      )
      named_choices <- stats::setNames(values, labels)
      
      choices <- c("All" = "All", named_choices)
    }
    
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
    
    # 3) community filter
    # 3) community filter
    comm_choice <- input$community
    if (!is.null(comm_choice) && comm_choice != "All") {
      status(paste0("Rendering network for community ", comm_choice, "…"))
      comm_choice_num <- as.numeric(comm_choice)
      nodes <- nodes[nodes$community == comm_choice_num, , drop = FALSE]
      used_ids <- nodes$id
      edges <- edges[edges$from %in% used_ids & edges$to %in% used_ids, , drop = FALSE]
    }
    
    # 4) drop nodes not in any edge
    # 4) drop nodes not in any edge
    used_ids <- unique(c(edges$from, edges$to))
    nodes    <- nodes[nodes$id %in% used_ids, , drop = FALSE]
    if (nrow(nodes) == 0 || nrow(edges) == 0) {
      status("No nodes/edges remain after filters. Try relaxing your options.")
      validate(need(FALSE, "No nodes/edges left after filtering."))
    }
    
    # 5) node colouring
    # 5) node colouring
    if (isTRUE(input$color_nodes_comm)) {
      nodes$group <- as.factor(nodes$community)
    } else {
      nodes$group <- NA
    }
    
    # 6) edge colouring by sign
    # 6) edge colouring by sign
    if (isTRUE(input$color_edges_by_sign)) {
      pos_col <- input$pos_color
      neg_col <- input$neg_color
      edges$color <- ifelse(edges$sign == "pos", pos_col, neg_col)
    }
    
    # 7) hide labels if requested
    # 7) hide labels if requested
    if (!isTRUE(input$show_labels)) {
      nodes$label <- NA
    }
    
    # 8) draw with browser physics first, then freeze
    # 8) draw with browser physics first, then freeze
    net <- visNetwork(nodes, edges) %>%
      visOptions(
        highlightNearest = TRUE,
        nodesIdSelection = TRUE
      ) %>%
      visPhysics(
        solver = "barnesHut",
        stabilization = list(
          enabled    = TRUE,
          iterations = 1500
          iterations = 1500
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