# NeuroShiny
 A Web app for looking at gene expression changes during human brain development

[https://rosie-griffiths.shinyapps.io/NeuroShiny/]https://rosie-griffiths.shinyapps.io/NeuroShiny/

## Installation

First run the following code to make sure you have all the packages to run the NeuroShiny app

```r
reqPkg = c("GeneOverlap", "visNetwork", "DT", "Matrix", "ggplot2", "VennDiagram","png","stringr",
           "shinydashboard", "plotly", "shiny", "RColorBrewer", "shinyalert", "shinyWidgets")
newPkg = reqPkg[!(reqPkg %in% installed.packages()[,"Package"])]
if(length(newPkg)){install.packages(newPkg)}
```

## Running the app

To run the app locally download the reposity to your machine and unzip it. Then change directory to the folder
```r
library(shiny)
setwd("~/NeuroShiny-main")
runApp()
```


