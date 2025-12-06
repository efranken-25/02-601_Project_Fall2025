# Building and Analyzing RNA Co-Expression Networks 

what application does, why we used the technologies, challenges we faced and features we hope to implement in the future

## Table of Contents
* [Installation](#installation)
* [Usage](#usage)
* [Features](#features)


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

# Include Credits
* Beth Vazquez Smith 
* Noemi Banda
* Emma Franken

