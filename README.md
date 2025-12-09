# Building and Analyzing RNA Co-Expression Networks 

## Table of Contents
* [Project Overview](#project-overview)
* [Why We Used These Technologies](#why-we-used-these-technologies)
* [Challenges We Faced](#challenges-we-faced)
* [Potential Future Features](#potential-future-features)
* [Installation](#installation)
* [Usage](#usage)
* [Features](#features)
* [Authors](#authors)


## Project Overview 
Our application takes Gene Expression Quantification RNA-Seq data and creates an interactive and customizable co-expression network that allows the user to explore modules, while adjusting degree density, visibility of positive or negative edges, color schemes, and exporting figures. Users can analyze a single network or directly compare two datasets (in our case, cancer types) side-by-side to asses differences in co-expression topology and structure. Alongside this, it will calculate the global network statistics, such as clustering coefficients and degree distribution, to assess overall co-expression structure between cancer types.

## Why We Used These Technologies
Go was chosen for the backend analysis because of its strong performance, concurrency model, and clean support for object-oriented principles such as modularization and method-based struct organization. These features make it easier to structure complex biological workflows while maintaining readable, maintainable code.

R Shiny was used to create the interactive interface due to its simplicity in linking statistical analysis with real-time visualizations. Various packages in R's library, specifically `visNetwork` and `colourpicker` for customization, allow for user-friendly visualizations and controls without requiring users to install additional frameworks. 

R Studio....

Integrating Go with R Shiny enables CPU-efficient backend processing while providing an accessible, interactive frontend tailored for exploratory biological workflows.

## Challenges We Faced
* Determining how to compare networks:
    - Comparing co-expression networks requires identifying which metrics are meaningful across datasets, since some statistics are designed to characterize a single network and are more biologically interpretable at that level. These include metrics such as module structure, network diameter, and centrality measures. In contrast, other metrics are more suitable for directly comparing two networks, such as differences in degree distributions, overall topology, density, or edge sign proportions. Selecting appropriate comparison strategies, and ensuring they were calculated consistently, was a major challenge.

* Handling large, high-dimensional data:
    - RNA-Seq correlation matrices include thousands of genes, resulting in millions of pairwise comparisons. Efficiently computing and filtering these correlations, managing memory usage, and ensuring rapid data transfer between Go and Shiny required careful tuning.

## Potential Future Features
* Module-level enrichment analysis
* Support for more than two datasets simultaneosuly, with comparitive visualization panels
* Implement soft modularity to allow calculation of the Jaccard similarity index within modules for a single network

## Installation

### Prerequisites
* **Go 1.24** or higher
    * Required Go packages (all listed in `go.mod` and `go.sum`; you might need to run `go get` to install them on your machine):
        * Go modules enabled (Go 1.11+)
        * `github.com/glycerine/golang-fisher-exact`
        * `golang.org/x/exp/stats`
        * `gonum.org/v1/gonum/stat`

* **R version 4.2.1** or higher
    * Required R packages: 
        * `shiny`
        * `visNetwork`
        * `colourpicker`
        * `shinycssloaders`
        * `here`


## Usage 

### Running the Go Analysis
1. Make sure you have Go 1.24 or higher installed and Go modules enabled. 
2. Clone this repository:
    ```bash
    git clone https://github.com/efranken-25/02-601_Project_Fall2025.git 
    cd 02-601_Project_Fall2025
3. Install any Go dependencies (if not already in `go.mod`):
    `go get ./...`
4. Run the main Go program:
    `./02-601_Project_Fall2025` for Mac or `02-601_Project_Fall2025.exe` for Windows
5. The output CSV files for gene networks and community/module assignments will be saved in the `ShinyApp/` folder. 

### Running the Shiny App (R)
1. Open **R** or **RStudio**
2. Make sure you have the required packages installed:
    `install.packages(c("shiny", "visNetwork", "colourpicker", "shinycssloaders", "here"))`
3. Set your working directory to the project folder (so that the Shiny app can access the CSVs):
    `setwd("/path/to/02-601_Project_Fall2025")`
4. Launch the Shiny app:
    `shiny::runApp()`
5. The app will open in your web browser and display the interactive gene network visualizations. 

## Features




## Authors
* **Beth Vazquez Smith** - [@bvazquezsmith](https://github.com/bvazquezsmith)
* **Noemi Banda** - [@b-noemi](https://github.com/b-noemi)
* **Emma Franken** - [@efranken-25](https://github.com/efranken-25)

