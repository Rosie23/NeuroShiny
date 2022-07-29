# NeuroShiny
 A Web app for looking at gene expression changes during human brain development

https://rosie-griffiths.shinyapps.io/NeuroShiny/

## Installation to run it locally

First run the following code to make sure you have all the packages to run the NeuroShiny app

```r
reqPkg = c("GeneOverlap", "visNetwork", "DT", "Matrix", "ggplot2", "VennDiagram","png","stringr",
           "shinydashboard", "plotly", "shiny", "RColorBrewer", "shinyalert", "shinyWidgets")
newPkg = reqPkg[!(reqPkg %in% installed.packages()[,"Package"])]
if(length(newPkg)){install.packages(newPkg)}
```

## Running the app

To run the app locally;
1. Download the repository to your machine and unzip it
2. In RStudio run the following code

```r
library(shiny)
setwd("~/NeuroShiny-main") #Change working directory to the NeuroShiny folder
runApp()
```


