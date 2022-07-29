options(rsconnect.http.trace = TRUE)
options(rsconnect.http.verbose = TRUE)
#library(BiocManager)
#options(repos = BiocManager::repositories())
#remotes::install_github("shenlab-sinai/GeneOverlap")

#Upload code and resources to github, allow users to host locally?

#library(directlabels)
library(GeneOverlap)
library(visNetwork)
library(DT)
library(Matrix)
library(shinydashboard)
library(plotly)
library(RColorBrewer)
library(shiny)
library(ggplot2)
library(shinyalert)
library(shinyWidgets)
library(dashboardthemes)
library(VennDiagram)
library(png)
library(stringr)

#Load Data
load("Supplementary Data.rda")
source("functions.R",local = TRUE)
stringDB_data <- read.table("stringDB_Data.txt.gz", sep="\t",header=T)
load("TTRUST_data.rda")
load("ASD DEGs.rda")
load("Gene lists.rda")
datasets <- readxl::read_xlsx("Datasets.xlsx",sheet=1)
dataset_table <- readxl::read_xlsx("Datasets.xlsx",sheet=2)
Allen_gene_code <- read.csv("Allen_gene_code.csv")
load("Gene_Info_Table.rda")
load("BrainSpan_regional_data.rda")
load("cell_type_plot_data.rda")
load("GO Info.rda")
load("Summary_scRNAseq.rda")

ui <- dashboardPage(
  
    dashboardHeader(title = "NeuroShiny"),

#list of all the items and subitems found in the side bar
    dashboardSidebar(
      tags$style(".skin-blue .sidebar .shiny-download-link { color: #444; }"),
      #sidebarMenu allows for tabs to restore and app to open directly on homepage 
      sidebarMenu(id = "menu",
     menuItem("Home", tabName = "Homepage", icon = icon("home")),
     menuItem("Explore Genes", selected = T, startExpanded = TRUE,
                 menuSubItem("Gene Table",tabName = "Genetab"),
                 menuSubItem("Gene Set Overlap (Two lists)", tabName="enrichment_test"),
                 menuSubItem("Gene Set Overlap (Multiple lists)", tabName="overlap_multiple"),
                 menuSubItem("Cluster Classification Tool", tabName = "cluster_overlap"),
                 menuSubItem("ASD DEGs (Velmeshev et al.,)",tabName = "ASD_DEGtab")
                 ),
     selectizeInput("Gene",
                    "Select Gene:",
                    choices= gene.list,
                    selected="PAX6"),
     menuItem("Regional Expression", icon = icon("chart-bar"),startExpanded = TRUE,
                 menuSubItem("Tissue", tabName = "tissue_region"),
                 menuSubItem("Regional DEGs (BrainSpan)",tabName = "Region_DEGs")),
     menuItem("Time-course", icon = icon("chart-bar"),startExpanded = TRUE,
                 menuSubItem("2D in vitro", tabName = "2D"),
                 menuSubItem("Organoid",tabName = "Organoid"),
                 menuSubItem("Primary", tabName = "tissue_timecourse"),
                 menuSubItem("Comparison", tabName = "Comparison")),
     menuItem("Cell Type", icon = icon("chart-bar"), startExpanded = TRUE,
                 menuSubItem("Combined Data", tabName = "combined_cell")),
     menuItem("Evolution", icon = icon("chart-bar"),startExpanded = TRUE,
                 menuSubItem("Kanton et al., 2019", tabName="Kanton"),
                 menuSubItem("Pollen et al., 2019", tabName = "Pollen"),
                 menuSubItem("Benito-Kwiencinski et al., 2021", tabName = "Benito")),
     menuItem("Gene Clusters",tabName = "modules"),
     menuItem("References", tabName = "refs"),
       
         #allows user to download report
        downloadButton("report", label="Generate report")
        
    )),
    
    body <- dashboardBody(
      tags$head(tags$style(HTML('
      .main-header .logo {
        font-family: "Georgia", Times, "Times New Roman", serif;
        font-weight: bold;
        font-size: 24px;
      }
    '))),
      
        tabItems(
          
          #homepage appearance
            tabItem(tabName = "Homepage", 
            tags$head(tags$style(HTML("
body {
font-size: 16px;
}

.box-header {
font-size: 17px;
font-weight: bold;
}
                    #image_intro{
background-image: url(background.png);
height: auto;
width: 100%;
position: absolute;
padding: 0;
left: 0;
margin-right: 0px;
margin-left: 0px;
margin-top: -100px;
box-shadow: 2px 2px 2px grey;
z-index: 0;
} 

h1 {
color: white;
  font-weight: bold;
position: relative;
font-size: 40px;
text-align: center;
z-index: 1;
text-shadow: 2px 2px #000000	;
margin-top: 100px;
}

#b_desc {
color: white;
position: relative;
text-align: center;
z-index: 1;
font-size: 25px;
text-shadow: 1px 1px #000000;
padding-bottom: 120px;

}

#home_box {
position: relative;
z-index: 0;
margin-bottom: 10px;
}
" ))),
                    imageOutput(outputId = "image_intro"),
                    fluidRow(h1("Welcome to NeuroShiny")),
                    fluidRow(textOutput(outputId = "b_desc")),
            fluidRow(box(id="home_box" ,"NeuroShiny allows you to explore gene expression changes during early brain development across multiple datasets.
                         Explore the Gene Table and Gene Overlap function to identify genes within multiple datasets of interest.
                         Once a gene of interest has been selected NeuroShiny generates multiple interactive plots of its expression across time, cell types and model systems.", status = "primary", width = 12)),
           
                    tags$br(),
                    fluidRow(
                      box(id="home_box", title = p("Getting started", 
                                    actionBttn(inputId = "help_info", label = NULL, style = "material-circle", color = "default", icon = icon("info"), size = "xs")), 
                          solidHeader = TRUE, status = "primary",
                          #useShinyalert(),
                          tags$br(),
                          actionButton(inputId = "home_gene_table2", "Explore gene table"),
                          tags$br(),
                          selectizeInput("Gene_h","Select Gene:",choices= gene.list, selected="PAX6"),
                          actionBttn(inputId = "gene_desc",label = "Gene Summary",style = "simple",color = "primary"),
                          #useShinyalert()
                          ), 
                      
                      box(id="home_box", title = "Feedback", status = "primary", solidHeader = TRUE, 
                          a(actionButton(inputId = "email1", label = "Contact us", icon = icon("envelope", lib = "font-awesome")), href="mailto:rosie.griffiths@ed.ac.uk"),
                          tags$br(),
                          "If you have any feedback such as: features or datasets to add, glitches identified, figures unclear etc. please contact us by clicking the button above.",
                          tags$br())
                    ),
                    
                    fluidRow(
                      box(id="home_box", title = "Generate and download report", width = 6, solidHeader = TRUE, status = "primary", 
                          "Generates a report containing all the plots of the selected gene of interest",
                          tags$br(),
                          downloadButton("report_h", label="Generate report")),
                      tabBox(width = 6,id="home_box",
                             #id="home_box", title = "Additional information", width = 6, solidHeader = TRUE, status = "primary", collapsible = TRUE,
                             tabPanel("Useful Links",
                                      tags$p(actionLink(inputId = "add_info_ref", "References of primary papers" )),
                                      tags$p(tags$a(href="https://gene.sfari.org", "SFARI (ASD)"))),
                             tabPanel("Human Links",
                                      tags$p(tags$a(href="https://www.brainspan.org", "BrainSpan")),
                                      tags$p(tags$a(href="https://cortecon.neuralsci.org", "CORTECON")),
                                      tags$p(tags$a(href="http://stemcell.libd.org/scb/", "LIBD Stem Cell Browser"))
                                      ),
                             tabPanel("Mouse Links",
                                      tags$p(tags$a(href="https://portal.brain-map.org", "Allen Brain Map")),
                                      tags$p(tags$a(href="https://gp3.mpg.de", "Gene Paint")),
                                      tags$p(tags$a(href="http://www.informatics.jax.org/expression.shtml", "Mouse Gene Expression Database (JAX)"))
                                      )
                      )
                    )
            ),
            
            #code for the gene table tab
            tabItem(tabName = "Genetab",
                    fluidRow(
                        tabBox(width=12,
                               #Choose Columns to display
                               tabPanel("Columns to show:",
                                        checkboxGroupInput("show_vars", "Columns to show:", inline=T,
                                                           choices = names(gene_info_column_lists), selected = c("Gene info","SFARI"))),
                               #Filter by Cell Type (Fan, Bhaduri) or Cell State (Bhaduri)
                               tabPanel("Filter by Cell Type", 
                                   checkboxGroupInput("select_celltype_fan", "Select Cell Type (Fan et al.):",
                                                      inline=T,choices = fan_celltypes),
                                   checkboxGroupInput("select_cellstate_bhaduri", "Select Cell State (Bhaduri et al.):",
                                                      inline=T,choices = bhaduri_cellstates),
                                   checkboxGroupInput("select_celltype_bhaduri", "Select Cell Type (Bhaduri et al.):",
                                                      inline=T,choices = bhaduri_celltypes),
                                   checkboxGroupInput("select_celltype_nowakowski", "Select Cell Type (Nowakowski et al.):",
                                                      inline=T,choices = nowakowski_celltypes),
                                   pickerInput(
                                     inputId = "select_cellsubtype_nowakowski",
                                     label = "Select Cell Subtype (Nowakowski et al.):", 
                                     choices = nowakowski_cellsubtypes,
                                     options = list(
                                       `actions-box` = TRUE,
                                       `live-search` = TRUE), 
                                     multiple = TRUE)
                                   ),

                        tabPanel("Filter by Region",
                                 pickerInput(
                                   inputId = "region_filter",
                                   label = "Brainspan:", 
                                   choices = brain.regions,
                                   options = list(
                                     title = "Nothing selected",
                                     `actions-box` = TRUE), 
                                   multiple = TRUE)
                                 ),

                    #Select GO term to filter gene list by
                               tabPanel("Other filters",
                                        pickerInput(
                                          inputId = "go_terms",
                                          choices = GO_term_list,
                                          label = "Filter by GO term",
                                          options = list(
                                            title = "No GO term selected",
                                            `actions-box` = TRUE,
                                            `live-search`=TRUE,
                                            size = 8),
                                          multiple = TRUE
                                        ),
                                        #selectizeInput("go_terms", "Filter by GO term (delete - to type)",choices= GO_term_list, multiple = FALSE),
            #if you set choices = c("select" = "", GO_term_list), you can type in the term you are looking for but the table doesn't show any gene
                                        checkboxInput("human_specific", "Human Specific", value = FALSE),
                                        checkboxInput("SFARI_gene", "SFARI Genes", value = FALSE),
                                        checkboxInput("user_gene_list_select", "Uploaded gene list", value = FALSE),
                                        fileInput("user_data", "Upload Gene List",multiple = FALSE,
                                                  accept = c("text/txt","text/tab-separated-values,text/plain",".txt"))
            ),
                        ),
                        
            box(DTOutput("mytable"), width=12),
            downloadButton("downloadGeneTable","Download Table")
                    )),
            
            #Gene Overlap Tabs
            tabItem(tabName = "enrichment_test",
                    tags$h2("Gene List Overlap Analysis",align = "center"),
                    box(width=6,title="Select Gene lists to overlap",
                        selectizeInput("enrichment_list1","List 1",choices=c(names(gene_lists)), selected=names(gene_lists)[4]),
                        selectizeInput("enrichment_list2","List 2",choices=names(gene_lists), selected=names(gene_lists)[6]),
                        selectizeInput("enrichment_background_list","Background List",choices=names(gene_lists)),
                        fileInput("user_gene_lists", "Upload Gene List",multiple = TRUE,
                                  accept = c("text/txt","text/tab-separated-values,text/plain",".txt"))),
                    box(width=6,title="Overlap Results",
                        htmlOutput("enrichment_text"),
                        downloadButton("downloadEnrichment", "Download Results")),
                    box(width=12,plotOutput("venn_diagram"))
            ),
            
            tabItem(tabName = "overlap_multiple",
                    tags$h2("Multiple Gene List Overlap Analysis",align = "center"),
                    box(width=3, title = "Select Gene lists",
                        selectInput("multiple_gene_list","Select Gene Lists",choices = names(gene_lists),
                                    multiple = T,selected=names(gene_lists)[3:5])),
                    box(width=9,title="Venn Diagram",
                        plotOutput("multipleVenn"),
                        downloadButton("downloadMultiVennOverlap","Download Overlap")),
                    box(width=12,title="UpSet Plot",
                        plotOutput("upsetPlot"))),
            
            tabItem(tabName = "cluster_overlap",
                    box(width=6,
                        actionButton(inputId = "help_cluster_overlap", "", icon = icon("info"), align = "right"),
                        fileInput("user_clusters_filename", "Upload Cluster Biomarkers (results of FindMarkers() from Seurat)",multiple = FALSE,
                                  accept = c("text/txt","text/tab-separated-values,text/plain",".txt")),
                        #actionButton("test_scRNAseq_cluster", "Run Cluster Test Data"),
                        #HTML(paste("",sep="<br/>")),
                        HTML(paste("File needs to be in the format .txt with column names;", 
                             	"p_val or power	| avg_log2FC | cluster |	gene",sep="<br/>"))),
                    box(width=3,
                        numericInput("cluster_FC_filter","Log2 Fold Change Filter", value=0),
                        radioButtons("stats_filter","Filter by P-value or Power",
                                     choices = c("P-value"="P-value","Power"="Power"))
                        ),
                    box(width = 3,
                        numericInput("cluster_pvalue_filter","P value Filter", value=0.0001),
                        numericInput("cluster_power_filter","Power Filter", value=0.4),
                        actionButton("run_cluster_analysis","Run Analysis")
                        ),
                    tabBox(width=12,
                           tabPanel("User Data", DTOutput("uploaded_cell_clusters")),
                           tabPanel("Results",DTOutput("fishers_cluster_overlap_table")),
                           tabPanel("Summary", DTOutput("summary_clusters_table")))
                    ),
            
            #code for the regional DEGs
            tabItem(tabName = "Region_DEGs", tags$h2("Regional DEGs (Differentially Expressed Genes)", align = "center"),
                    box(width=6, radioButtons("region_DEGs_age_range", "Select Age Range", choices = c("Early (8-13 p.c.w.)"=1, "Mid (16-21 p.c.w)"=2,"Late (24-37 p.c.w)"=3),selected=1), numericInput("region_DEGs_FDR","Select FDR cutoff", value=0.01)),
                    
                    box(width=6, selectInput("region_DEGs_xaxis","Select Region to plot (x-axis)", choices = brain.regions.DEGs.names, selected = "ACtx"),
                        selectInput("region_DEGs_yaxis","Select Region to plot (y-axis)", choices = brain.regions.DEGs.names, selected = "Cerebellum")),
                    
                    box(width=12, plotlyOutput("PlotDEGsRegion"))
                    #tabBox(width = 12, title = "Additional information", tabPanel("Study design", HTML(paste("method", sep="<br/>"))),
                    #       tabPanel(width=12, "Data processing", HTML(paste("analysis", sep="<br/>"))))
                    ),
            
            tabItem(tabName = "ASD_DEGtab",
                    h2("Cell type–specific gene expression changes in ASD",align = "center"),
                    box(plotlyOutput("PlotASD_DEG_neuronal")),
                    box(plotlyOutput("PlotASD_DEG_non_neuronal")),
                    tabBox(width=12, title = "Study design", tabPanel("Velmeshev", HTML(paste("Single-cell genomics identifies
cell type–specific molecular changes in autism",
                                                                                              "-------------------",
                                                                                              "DOI: 10.1126/science.aav8130",
                                                                                              "-------------------", "Protocol: Postmortem tissue samples (including the prefrontal cortex and anterior cingulate cortex) were gathered from 15 ASD patients and 16 controls. The experimental and control population were matched for age, sex, postmortem interval and RNA integrity number. Unbiased snRNA-seq and bulk RNA-seq were performed to identify ASD-associated genes and compare the results to the SFARI database.",  

"Additional snRN-seq data was performed to compare ASD and sporadic epilepsy using PFC samples from 8 patients with sporadic epilepsy and 7 matched controls. ", sep="<br/>"))))),
            
            tabItem(tabName = "2D",
                    h2("2D Stem Cell Neural Differentiation",align = "center"),
                    fluidRow(
                        box(plotlyOutput("PlotInoue", height = 250)),
                        box(plotlyOutput("PlotLIBD", height = 250)),
                        box(plotlyOutput("PlotCORTECON", height = 250)),
                        box(plotlyOutput("PlotMicali", height = 250)),
                        tabBox(title = "Study Design",width=12, id = "tabset1", height = "250px", # The id lets us use input$tabset1 on the server to find the current tab
                               tabPanel("Inoue", HTML(paste("Identification and Massively Parallel Characterization of Regulatory Elements Driving Neural Induction",
                                                            "-------------------",
                                                            "https://doi.org/10.1016/j.stem.2019.09.010",
                                                            "-------------------",
                                                            "Cell line: hESCs (H1)" ,

"Protocol: Cells were cultured on Matrigel and a Dual-Smad inhibition protocol was used to induce neural differentiation. Cells were harvested at 7 time points (0h, 3h, 6h, 12h, 24h, 48h, 72h) and RNA-seq, ChIP-seq and ATAC-seq were performed for each time point, followed by comprehensive functional assays (MPRAs) to understand the genes and regulatory elements involved in orchestrating neural induction. ", sep="<br/>"))),
                               
                               tabPanel("LIBD", HTML(paste("Dissecting transcriptomic signatures of neuronal differentiation and maturation using iPSCs",
                                                           "-------------------",
                                                           "https://doi.org/10.1038/s41467-019-14266-z",
                                                           "-------------------",
                                                           "Cell line: hiPSCs" ,

"Protocol: RNA-seq was performed at 9 different time points during differentiation corresponding to the transition stages between 5 different differentiation conditions (self-renewal, early neuronal differentiation, neuronal precursor cells, assembled rosettes, differentiated neuronal cells).", sep="<br/>"))),
                               
tabPanel("Micali", HTML(paste("Variation of Human Neural Stem Cells Generating Organizer States In Vitro before Committing to Cortical Excitatory or Inhibitory Neuronal Fates",
                                                             "-------------------",
                                                             "https://doi.org/10.1016/j.celrep.2020.107599",
                                                             "-------------------",
                                                             "Cell lines: hPSC (H9), hiPSC lines (2063-1, _2, 2053-2, _6, 2075-1, _3)", 

"Protocol: Dual-SMAD inhibition protocol was used to induce neural differentiation in hPSC. Forebrain NSCs were generated from hiPSCs using the protocol by Maroof et al. (2013). Bulk RNA-seq was conducted on cells collected at day 30 (neurons without astrocytes) and day 32 (neurons co-cultured with astrocytes). ", sep="<br/>"))),

                               tabPanel("CORTECON", HTML(paste("CORTECON: A Temporal Transcriptome Analysis of In Vitro Human Cerebral Cortex Development from Human Embryonic Stem Cells",
                                                               "-------------------",
                                                               "http://dx.doi.org/10.1016/j.neuron.2014.05.013",
                                                               "-------------------",
                                                               "Cell line: hESCs (WA-09, WiCell)", 

"Protocol: The dual SMAD-inhibition protocol was performed to initiate neural progenitors induction and dorsal telencephalic fate was ensured by adding cyclopamine from day 3 of induction. Cultures were maintained in N2B27 medium and supplemented with 10ng/ml of FGF2.", 

"RNA-seq was performed at several time points during differentiation to determine changes in gene expression during corticogenesis. ", sep="<br/>")))
                        )
                    )
            ),
            
            tabItem(tabName = "Organoid",
                    h2("Human Cortical Organoid",align = "center"),
                    box(plotlyOutput("PlotOrganoidBenito", height = 250)),
                    box(plotlyOutput("PlotOrganoidSloan", height = 250)),
                    box(plotlyOutput("PlotOrganoidPollen", height = 250)),
                    box(plotlyOutput("PlotOrganoidKanton", height = 250)),
                    box(plotlyOutput("PlotOrganoidCamp", height = 250)),
                    tabBox(title = "Study Design",width=12,
                           id = "tabset2", height = "250px",
                           tabPanel("Benito", HTML(paste(" 
                                                         An early cell shape transition drives evolutionary
expansion of the human forebrain",
                                                         "-------------------",
                                                         "https://doi.org/10.1016/j.cell.2021.02.050",
                                                         "-------------------",
                                                         "Cell lines: H9 (ESC, human), iPS(IMR90)-4 (iPSC, human), C3651 (iPSC, chimpanzee), goiPSC clone 1 (iPSC, gorilla), gorC1 (iPSC, gorilla)",

"Protocol: Cerebral organoids were generated using modified protocol described by Lancaster et al. (2017) to generate telencephalic identities. An identical method was used for all species to ensure comparability, but chimpanzee organoids were transferred into neural induction 2 days earlier than others. Bulk RNA-seq was performed at several time points from day 0 to 25 (fully committed neurogenesis). ",
"Note: while the study looked at multiple species, this graph only contains human data", sep="<br/>"))),

                           tabPanel("Sloan", HTML(paste("Human Astrocyte Maturation Captured in 3D Cerebral Cortical Spheroids Derived from
Pluripotent Stem Cells", 
                                                        "-------------------",
                                                        "http://dx.doi.org/10.1016/j.neuron.2017.07.035",
                                                        "-------------------",
                                                        "Cell lines: 4 iPSC lines (derived from fibroblasts)", 

"Protocol: 33 cortical spheroids (hCSs) containing astrocytes were generated from iPSCs following the protocol by Pasca et al. (2015), they were maintained in floating conditions on low attachment. Astrocyte lineage cells were purified through immunopanning using anti-HepaCAM and anti-Thy1 antibodies. The purified astrocytes were analysed using scRNA-seq at many time points between day 100 and 450 of culture. The gene expression profiles were compared to that of primary astrocytes isolated from foetal and adult CNS from 27 brain samples. "
                                                        , sep="<br/>"))),
                           tabPanel("Pollen", HTML(paste("Establishing Cerebral Organoids as Models of Human-Specific Brain Evolution", 
                                                         "-------------------",
                                                         "https://doi.org/10.1016/j.cell.2019.01.017",
                                                         "-------------------",
                                                         "Cell lines: iPSC (described by Bershteyn et al., 2017; Gallego Romero et al., 2015; Pavlovic et al., 2018), 4 additional iPSC lines (reprogrammed from fibroblasts using Okita et al., 2011 protocol) ", 
                                                         
                                                         "Protocol: scRNA-seq was performed on 48 human primary tissue samples, 6 macaque primary tissue samples and 56 organoids derived from 10 humans and 8 macaques. The samples were distributed across neurogenesis stages. The organoids were differentiated using protocols by Kadoshima et al. (2013) and Bersheteyn et al. (2017), the cortical differentiation medium was supplemented with Rho-Kinase inhibitor and TGFβ inhibitor.",
                                                         "Note: while the study looked at multiple species, this graph only contains human data", sep="<br/>"))),

                           tabPanel("Kanton", HTML(paste("Organoid single-cell genomic atlas uncovers human-specific features of brain development", 
                                                         "-------------------",
                                                         "https://doi.org/10.1038/s41586-019-1654-9",
                                                         "-------------------",
                                                         "Cell lines: ESC (H9, human), iPSCs (409b2, human), iPSCs (chimpanzee), iPSCs (bonobo, reprogrammed from fibroblasts with StemMACS mRNA transfection kit)", 
                                                         
                                                         "Protocol: Organoids were cultured using the protocol described by Lancaster et al. (2014). Every month, from 0 to 4 months organoids were dissociated to perform scRNA-seq with accessible chromatin profiling. Furthermore, snRNA-seq was performed on human, chimpanzee, and bonobo adult prefrontal cortex tissue to compare gene expression profiles between species and primary samples/organoids.",
                                                         "Note: while the study looked at multiple species, this graph only contains human data", sep="<br/>"))),

                           tabPanel("Camp", HTML(paste("Human cerebral organoids recapitulate gene expression programs of fetal neocortex development", 
                                                       "-------------------",
                                                       "https://doi.org/10.1073/pnas.1520760112",
                                                       "-------------------",
                                                       "Cell lines: PS409b2 (iPSC), ESCs ", 
                                                       
                                                       "Protocol: Primary foetal neocortex samples (12-13 wpc) and human organoids were analysed using scRNA-seq. 5 iPSC derived organoids were then generated following the protocol by Lancaster et al. (2013) and Lancaster & Knoblich (2014), additionally mTESR1 was added during EB formation. scRNA-seq was performed on whole organoids 35, 37, 41 and 65 days after EB formation, and on days 53 and 58 for 4 microdissected cortical regions (2 ESC derived and 2 iPSCs derived). Time lapse microscopy and immunostaining for bIP markers were also used on top of scRNA-seq for both organoids and primary samples.", sep="<br/>")))
                    )),
            
            tabItem(tabName = "tissue_region",
                    tags$h2("Regional Expression", align ="center"),
                    tabsetPanel(
                      tabPanel("Human (BrainSpan)", fluidRow(
                        box(width = 8, plotlyOutput("PlotBrainSpan", height = "500px", width = "500px")),
                        box(width = 4, title = "Acronyms", solidHeader = TRUE, status = "primary", collapsible = TRUE, collapsed = TRUE,
                            tags$ul(
                              tags$li("MFC: medial frontal cortex"),
                              tags$li("OFC: orbitofrontal cortex"),
                              tags$li("DFC: dorsofrontal cortex"),
                              tags$li("VFC: ventral frontal cortex"),
                              tags$li("M1C: primary motor cortex"),
                              tags$li("S1C: primary sensory cortex"),
                              tags$li("IPC: inferior parietal cortex"),
                              tags$li("TCx: temporal cortex"),
                              tags$li("ITC: caudal inferior temporal cortex"),
                              tags$li("STC: superior temporal cortex"),
                              tags$li("A1C: primary auditory cortex"),
                              tags$li("Ocx: occipital cortex"),
                              tags$li("V1C: primary visual cortex"),
                              tags$li("CGE: caudal ganglionic eminence"),
                              tags$li("LGE: lateral ganglionic eminence"),
                              tags$li("MGE: medial ganglionic eminence"),
                              tags$li("STR: striatum"),
                              tags$li("HIP: hippocampus"),
                              tags$li("AMY: amygdala"),
                              tags$li("DTH: dorsal thalamus"),
                              tags$li("MD: medial dorsal nucleus"),
                              tags$li("URL: upper rhombic lip"),
                              tags$li("CB: cerebellum"),
                              tags$li("CBC: cerebellar cortex")
                            )),
                        box(width = 4, offset = 8, title = "Link to original data", solidHeader = TRUE, status = "primary", collapsible = TRUE,
                            uiOutput(outputId = "bulk_links_bs")),
                        tabBox(width = 12, title = "", tabPanel("Study design", HTML(paste("Protocol: post-mortem brain tissue specimen we collected ranging from embryonic stages to late adulthood.
                                                                                           Samples were sectioned and different brain regions (8-16 brain structures) were separated. Bulk RNA-seq was conducted.",
                                                                                           "-------------------",
                                                                                           "Note: only this graph only contain pre-natal data points ", sep="<br/>")))
                               )
                      )),
                      
                      tabPanel("Human (HDBR)", fluidRow(
                        box(width = 8, plotlyOutput("PlotHDBR", height = "500px", width = "500px")),
                        box(width = 4, offset = 8, title = "Link to original data", solidHeader = TRUE, status = "primary", collapsible = TRUE,
                            tags$p(tags$a(href="https://www.hdbr.org/expression", "HDBR Homepage"))
                            ),
                        
                        tabBox(width = 12, title = "", tabPanel("Study Design",HTML(paste(
                        "Protocol:Bulk RNA-seq was performed on human fetal and embryonic tissue samples obtained from the HDBR tissue bank.
                        A total of 557 RNA-seq datasets were collected from different regions from 172 brains at ages between 4 and 17PCW.", sep="<br/>")))
                        )
                               
                      )),
                      
                      tabPanel("Mouse (Allen)", fluidRow(
                        box(width = 8, plotlyOutput("PlotAllen", height = "500px", width = "500px")),
                        box(width = 4, offset = 8, title = "Link to original data", solidHeader = TRUE, status = "primary", collapsible = TRUE,
                            uiOutput(outputId = "bulk_links_al"),
                            uiOutput(outputId = "bulk_links_paint")),
                        box(width = 4, title = "Acronyms", solidHeader = TRUE, status = "primary", collapsible = TRUE, collapsed = TRUE,
                            tags$ul(
                              tags$li("RSP: Rostral Secondary Prosencephalon"),
                              tags$li("Tel: Telencephalon"),
                              tags$li("PHy: Penduncular (caudal) Hypothalamus"),
                              tags$li("p3: Prosomere 3"),
                              tags$li("p2: Prosomere 2"),
                              tags$li("p1: Prosomere 1"),
                              tags$li("M: Midbrain"),
                              tags$li("PPH: Prepontine Hindbrain"),
                              tags$li("PH: Pontine Hindbrain"),
                              tags$li("PMH: Pontomedullary Hindbrain"),
                              tags$li("MH: Medullary Hindbrain (medulla)")
                            )),
                        tabBox(width = 12, title = "",tabPanel("Study Design",HTML(paste(
                                                                           "Protocol: mouse brain tissue was collected at 7 time points (4 pre-natal, 3 post-natal). The tissue was then frozen, sectioned and in situ hybridization was performed. The genes selected for analysis had strong established developmental roles.",
                                                                           "-------------------",
                                                                           "Note: this graph only contain pre-natal data", sep="<br/>"))))
                      ))
                      )),
            
            tabItem(tabName = "tissue_timecourse",
                    tags$h2("Primary tissue gene expression comparison", align = "center"),
                    tabBox(width=12,
                        tabPanel("Plot Method",radioButtons("plot_type_primary_timecourse","Plot Method",
                                                            choices = c("Smooth (loess)","Summary (mean)","Line","Line (smooth)"),selected="Summary (mean)")),
                        tabPanel("Region Filter",checkboxGroupInput("select_regions_brianspan_timecourse","Select Regions (BrainSpan):",
                                           choices = brain.regions,selected=brain.regions,inline=T),
                        checkboxGroupInput("select_regions_allen_timecourse","Select Regions (Allen):",
                                           choices = allen_regions, selected = allen_regions,inline=T))),
                    box(plotlyOutput("PlotBrainSpan_timecourse")),
                    box(plotlyOutput("PlotAllen_timecourse"))),
            

          tabItem(tabName = "combined_cell",
                  h2("Cell Type Gene Differential Expression", align="center"),
                  tabBox(width=8,
                         tabPanel("Scatter Plot",plotlyOutput("PlotCombinedCell")),
                         tabPanel("Dataset Summary",plotlyOutput("PlotSummaryCell"))),
                  imageOutput("CellColLedgend"),
                  box(width=12,tableOutput("DatasetSummaryTable")) #DTOutput
                      ),  

          tabItem(tabName="Comparison",
                  h2("Time-course Gene Expression in Different Model Systems", align = "center"),
                  box(width=6,
                      sliderInput("slider_comparison", "Day Range", min = 0, max = 100, value = c(0, 100)),
                      selectInput("Facet", "Sort by:", choices = c("Model","Species","none"),selected="Model")),
                  box(width=6,
                      selectInput("ColourPal", "Colour Palette", choices = rownames(brewer.pal.info), selected="Spectral"),
                      selectInput("ColourBy1", "Colour By",choices = c('Experiment',"Model","Species"), selected='Experiment')),
                  box(plotlyOutput("PlotComparison"),width=12)
                  #tabBox(width = 12, title = "Data processing", tabPanel("Method", HTML(paste("", sep="<br/>"))))
                  ),
            
            tabItem(tabName = "Kanton",
                    tags$h2("Human vs. Chimp Gene Expression Comparison across time course and cell types", align = "center"),
                    box(checkboxGroupInput("select_cells_kanton", "Select Cell Types:",
                                           inline=T,
                                           choices = kanton_cell_order, selected = kanton_cell_order), width=12),
                    fluidRow(
                        tabBox(width=12,
                               tabPanel("Time-course",plotlyOutput("PlotKanton_timecourse", height="800px")),
                               tabPanel("Cell Type Expression",plotlyOutput("PlotKanton_violin", height="750px", width = "950px")))),
                    tabBox(width = 12, title = "Study design", tabPanel("Kanton", HTML(paste("Organoid single-cell genomic atlas uncovers human-specific features of brain development", 
                                                                                             "-------------------",
                                                                                             "https://doi.org/10.1038/s41586-019-1654-9",
                                                                                             "-------------------",
                                                                                             "Cell lines: ESC (H9, human), iPSCs (409b2, human), iPSCs (chimpanzee), iPSCs (bonobo, reprogrammed from fibroblasts with StemMACS mRNA transfection kit)", 
                                                               "Protocol: Organoids were cultured using the protocol described by Lancaster et al. (2014). Every month, from 0 to 4 months organoids were dissociated to perform scRNA-seq with accessible chromatin profiling. Furthermore, snRNA-seq was performed on human, chimpanzee, and bonobo adult prefrontal cortex tissue to compare gene expression profiles between species and primary samples/organoids.",
                                                                                             "Note: while the study looked at multiple species, this graph only contains human data", sep="<br/>"))))),
            
            tabItem(tabName = "Pollen",
                    
                    fluidRow(
                      box(width=12,
                          checkboxGroupInput("select_cells_pollen", "Select Cell Types:",
                                             inline=T,choices = pollen_celltype_order, selected = pollen_celltype_order)), 
                        
                        tabBox(width=12,
                               tabPanel("Time-course",plotlyOutput("PollenPlot_timecourse", height="600px")),
                               tabPanel("Cell Type",plotlyOutput("PollenPlot_violin", height="600px")),
                               tabPanel("DEGs", plotlyOutput("DEGPlot")))
                    ),
                    tabBox(width = 12, height = "250px", title = "Study design", tabPanel("Pollen", HTML(paste("Establishing Cerebral Organoids as Models of Human-Specific Brain Evolution", 
                                                                                                               "-------------------",
                                                                                                               "https://doi.org/10.1016/j.cell.2019.01.017",
                                                                                                               "-------------------",
                                                                                                               "Cell lines: iPSC (described by Bershteyn et al., 2017; Gallego Romero et al., 2015; Pavlovic et al., 2018), 4 additional iPSC lines (reprogrammed from fibroblasts using Okita et al., 2011 protocol) ", 
            "Protocol: scRNA-seq was performed on 48 human primary tissue samples, 6 macaque primary tissue samples and 56 organoids derived from 10 humans and 8 macaques. The samples were distributed across neurogenesis stages. The organoids were differentiated using protocols by Kadoshima et al. (2013) and Bersheteyn et al. (2017), the cortical differentiation medium was supplemented with Rho-Kinase inhibitor and TGFβ inhibitor.", sep = "<br/>"))))
            ),
            
            tabItem(tabName = "Benito",
                    tags$h2("Human vs. Gorilla Gene Expression Comparison", align= "center"),
                    fluidRow(column(HTML(paste("","", sep="<br/>")), width=12),
                             
                             box(sliderInput("slider_benito", "Day Range",
                                             min = 5, max = 30, value = c(8, 15)), width=12)),
                    box(plotlyOutput("benito_plot"),width=12),
                    fluidRow(tabBox(width = 12, height = "250px", title = "Study design", tabPanel("Benito", HTML(paste(" 
                                                         An early cell shape transition drives evolutionary
expansion of the human forebrain",
                                                                                                               "-------------------",
                                                                                                               "https://doi.org/10.1016/j.cell.2021.02.050",
                                                                                                               "-------------------",
                                                                                                               "Cell lines: H9 (ESC, human), iPS(IMR90)-4 (iPSC, human), C3651 (iPSC, chimpanzee), goiPSC clone 1 (iPSC, gorilla), gorC1 (iPSC, gorilla)",
                        "Protocol: Cerebral organoids were generated using modified protocol described by Lancaster et al. (2017) to generate telencephalic identities. An identical method was used for all species to ensure comparability, but chimpanzee organoids were transferred into neural induction 2 days earlier than others. Bulk RNA-seq was performed at several time points from day 0 to 25 (fully committed neurogenesis). ", sep = "<br/>")))))),
            
            tabItem(tabName = "modules",
                    fluidRow(h3("Co-expression modules containing selected gene",align = "center")),
                    box( width=12,
                         actionButton(inputId = "help_cluster", "", icon = icon("info"), align = "right"),
                         DTOutput("module_table")),
                    selectizeInput(
                      inputId = "module", 
                      label = "Select Module",
                      multiple  = F,
                      choices = ""
                    ),
                    #uiOutput("select_module"),
                    tabBox(width=12,
                           tabPanel("PPI Network",visNetworkOutput("module_network_plot")),
                           tabPanel("Module Gene List",DTOutput("module_gene_table")),
                           tabPanel("TFs",DTOutput("module_TF_table")),
                           tabPanel("TFs Fishers",DTOutput("module_TF_table_fishers")),
                           tabPanel("Top TFs",DTOutput("top_TFs")),
                           tabPanel("Top TFs Plot", plotlyOutput("top_TFs_plot"))
                           
                    )
            ),
            
            tabItem(tabName = "refs",
                    tabItem(tabName = "refs",
                            HTML(paste("<b>",datasets$Name[1],"</b>",datasets$Paper[1])),
                            tags$a(href=datasets$doi[1], "doi"),
                            HTML(paste("","---------------------------------","", sep="<br/>")),
                            
                            HTML(paste("<b>",datasets$Name[2],"</b>",datasets$Paper[2])),
                            tags$a(href=datasets$doi[2], "doi"),
                            HTML(paste("","---------------------------------","", sep="<br/>")),
                            
                            HTML(paste("<b>",datasets$Name[3],"</b>",datasets$Paper[3])),
                            tags$a(href=datasets$doi[3], "doi"),
                            HTML(paste("","---------------------------------","", sep="<br/>")),
                            
                            HTML(paste("<b>",datasets$Name[5],"</b>",datasets$Paper[5])),
                            tags$a(href=datasets$doi[5], "doi"),
                            HTML(paste("","---------------------------------","", sep="<br/>")),
                            
                            HTML(paste("<b>",datasets$Name[4],"</b>",datasets$Paper[4])),
                            tags$a(href=datasets$doi[4], "doi"),
                            HTML(paste("","---------------------------------","", sep="<br/>")),
                            
                            HTML(paste("<b>",datasets$Name[6],"</b>",datasets$Paper[6])),
                            tags$a(href=datasets$doi[6], "doi"),
                            HTML(paste("","---------------------------------","", sep="<br/>")),
                            
                            HTML(paste("<b>",datasets$Name[7],"</b>",datasets$Paper[7])),
                            tags$a(href=datasets$doi[7], "doi"),
                            HTML(paste("","---------------------------------","", sep="<br/>")),
                            
                            HTML(paste("<b>",datasets$Name[8],"</b>",datasets$Paper[8])),
                            tags$a(href=datasets$doi[8], "doi"),
                            HTML(paste("","---------------------------------","", sep="<br/>")),
                            
                            HTML(paste("<b>",datasets$Name[9],"</b>",datasets$Paper[9])),
                            tags$a(href=datasets$doi[9], "doi"),
                            HTML(paste("","---------------------------------","", sep="<br/>")),
                            
                            HTML(paste("<b>",datasets$Name[10],"</b>",datasets$Paper[10])),
                            tags$a(href=datasets$doi[10], "doi"),
                            HTML(paste("","---------------------------------","", sep="<br/>")),
                            
                            HTML(paste("<b>",datasets$Name[11],"</b>",datasets$Paper[11])),
                            tags$a(href=datasets$doi[11], "doi"),
                            HTML(paste("","---------------------------------","", sep="<br/>")),
                            
                            HTML(paste("<b>",datasets$Name[12],"</b>",datasets$Paper[12])),
                            tags$a(href=datasets$doi[12], "doi"),
                            HTML(paste("","---------------------------------","", sep="<br/>")),
                            
                            HTML(paste("<b>",datasets$Name[13],"</b>",datasets$Paper[13])),
                            tags$a(href=datasets$doi[13], "doi"),
                            HTML(paste("","---------------------------------","", sep="<br/>")),
                            
                            HTML(paste("<b>",datasets$Name[14],"</b>",datasets$Paper[14])),
                            tags$a(href=datasets$doi[14], "doi"),
                            HTML(paste("","---------------------------------","", sep="<br/>")),
                            
                            HTML(paste("<b>",datasets$Name[15],"</b>",datasets$Paper[15])),
                            tags$a(href=datasets$doi[15], "doi"),
                            HTML(paste("","---------------------------------","", sep="<br/>")),
                            
                            HTML(paste("<b>",datasets$Name[16],"</b>",datasets$Paper[16])),
                            tags$a(href=datasets$doi[16], "doi"),
                            HTML(paste("","---------------------------------","", sep="<br/>")),
                            
                            HTML(paste("<b>",datasets$Name[17],"</b>",datasets$Paper[17])),
                            tags$a(href=datasets$doi[17], "doi"),
                            HTML(paste("","---------------------------------","", sep="<br/>")),
                            
                            HTML(paste("<b>",datasets$Name[18],"</b>",datasets$Paper[18])),
                            tags$a(href=datasets$doi[18], "doi"),
                            HTML(paste("","---------------------------------","", sep="<br/>")),
                            
                            HTML(paste("<b>",datasets$Name[19],"</b>",datasets$Paper[19])),
                            tags$a(href=datasets$doi[19], "doi"),
                            HTML(paste("","---------------------------------","", sep="<br/>")),
                            
                            HTML(paste("<b>",datasets$Name[20],"</b>",datasets$Paper[20])),
                            tags$a(href=datasets$doi[20], "doi"),
                            HTML(paste("","---------------------------------","", sep="<br/>"))
                    )     
) )))

server <- function(input, output, session) {
  
  output$image_intro <- renderImage(expr = list(src = "background.png"), deleteFile = FALSE)
  
  output$b_desc <- renderText("An interactive webapp to explore development transcriptomics")
    
  #links button in homepage to gene table
  observeEvent(input$home_gene_table2, {
    newtab2 <- switch(input$menu, "Homepage" = "Genetab","Genetab" = "Homepage")
    updateTabItems(session, "menu", newtab2)
  })

  #syncs the two select gene options
  observeEvent(input$Gene_h, {
    updateSelectizeInput(inputId = "Gene", selected = input$Gene_h)
  })
  observeEvent(input$Gene, {
    updateSelectizeInput(inputId = "Gene_h", selected = input$Gene)
  })
  
  #creates the info pop up for the getting started box
  observeEvent(input$help_info, {
    shinyalert(
      title = "Getting started",
      text = "If you have a gene in mind: Select a gene from the drop-down menu, an overview of the selected gene can be found at the 'Gene overview' button \n 
Explore the plots of the selected gene by navigating the sidebar menu. \n
If you do not have a gene in mind: click on the “Explore gene” button to view a list of interesting genes. \n
You can also upload your own gene list and identify genes within your list which overlap with other datasets within the app e.g. GO term genes, cell type genes, regional genes, ASD-genes, and/or human evolution genes.",
      size = "s", 
      closeOnEsc = TRUE,
      closeOnClickOutside = TRUE,
      html = FALSE,
      type = "info",
      showConfirmButton = TRUE,
      showCancelButton = FALSE,
      confirmButtonText = "Close",
      confirmButtonCol = "blue"
    )
  })
  
  #creates gene description pop up
  observeEvent(input$gene_desc, {
    gene_row <- which(gene_info$Gene == input$Gene)
    sfari <- if(is.na(gene_info[gene_row,5])) {"No"} else{"Yes"}
    #human <- if(is.na(gene_info[gene_row,18])) {"No"} else{"Yes"}
    shinyalert(
      title = paste0(input$Gene, " overview:"),
      text = paste0("Gene name: ", gene_info[gene_row,1], 
                    "\n Ensembl ID: ", gene_info[gene_row,2],
                    "\n Chromosome location: ", gene_info[gene_row,4],
                    "\n Is it a SFARI gene: ", sfari,
                    #"\n Class: ", gene_info[gene_row,11],
                    "\n Cell Type: ", gene_info[gene_row,16],
                    #"\n Cell subtype: ", gene_info[gene_row,17],
                    "\n Region enriched expression: ", gene_info[gene_row,21]
      ),
      size = "s", 
      closeOnEsc = TRUE,
      closeOnClickOutside = TRUE,
      html = FALSE,
      #type = "info",
      showConfirmButton = TRUE,
      showCancelButton = FALSE,
      confirmButtonText = "Close",
      confirmButtonCol = "blue"
    )
  })
  
  #help box in cluster page
  observeEvent(input$help_cluster, {
    shinyalert(
      title = "Co-expression modules",
      text = "Many studies have used Weighted correlation network analysis (WGCNA) to find clusters (modules)
      of highly correlated genes. Here we have collated modules identifed from multiple studies allowing users
      to easily find genes highly co-expressed with their gene of interest. Once a module is selected a Protein-Protein
      interaction (PPI) network is generated using interactions from the STRING database.",
      size = "s", 
      closeOnEsc = TRUE,
      closeOnClickOutside = TRUE,
      html = FALSE,
      type = "info",
      showConfirmButton = TRUE,
      showCancelButton = FALSE,
      confirmButtonText = "Close",
      confirmButtonCol = "blue"
    )
  })
  #links the links in additonal info box to the right pages
  observeEvent(input$add_info_ref, {
    home2ref <- switch(input$menu, "Homepage" = "refs","refs" = "Homepage")
    updateTabItems(session, "menu", home2ref)
    
  })
  
  #creates links text output in brainspan bulk RNA page
  output$bulk_links_bs <- renderUI(
    a("BrainSpan gene page", href = paste0(
        "https://www.brainspan.org/rnaseq/searches?exact_match=true&search_term=", input$Gene, "&search_type=gene"
      
    ))
    
  )
  
  output$bulk_links_al <- renderUI(
    
    a("Allen brain atlas gene page", href = paste0(
      if(!is.na("a")){
        num <- grep(input$Gene, Allen_gene_code$gene_symbol, ignore.case = TRUE, value = FALSE)
      paste0("http://developingmouse.brain-map.org/gene/show/", Allen_gene_code$id[num], "/"
        )}
      )
      )
  )
  
  output$bulk_links_paint <- renderUI(
    
    a("Genepaint gene page", href = paste0(
      "https://gp3.mpg.de/results/", input$Gene
    )
    )
  )
  
  
  #all of these link changes in the selected gene to changes in the graphs (I think)  
    NPC_df <- reactive({
      null_all_data()
      extract_data(gene=input$Gene, dataset_dir = "NPC")
      get(load(paste0("Data by Gene/NPC/files/",input$Gene,".RData")))
    })
    
    comparison_df <- reactive({
      null_all_data()
      extract_data(gene=input$Gene, dataset_dir = "Comparison")
        get(load(paste0("Data by Gene/Comparison/files/",input$Gene,".RData")))
    })
    
    brainspan_gene_df <- reactive({
      null_all_data()
      extract_data(gene=input$Gene, dataset_dir = "BrainSpan")
      get(load(paste0("Data by Gene/BrainSpan/files/",input$Gene,".RData")))
    })
    
    HDBR_gene_df <- reactive({
      null_all_data()
      extract_data(gene=input$Gene, dataset_dir = "HDBR")
      get(load(paste0("Data by Gene/HDBR/files/",input$Gene,".RData")))
    })
    
    kanton_gene <- reactive({
      null_all_data()
      extract_data(gene=input$Gene, dataset_dir = "Kanton")
        get(load(paste0("Data by Gene/Kanton/files/",input$Gene,".RData")))
    })
    
    benito_gene_df<- reactive({
      null_all_data()
      extract_data(gene=input$Gene, dataset_dir = "Benito")
        get(load(paste0("Data by Gene/Benito/files/",input$Gene,".RData")))
    })
    
    pollen_gene_df<- reactive({
      null_all_data()
      extract_data(gene=input$Gene, dataset_dir = "Pollen")
      get(load(paste0("Data by Gene/Pollen/files/",input$Gene,".RData")))
    pollen_gene_df <- subset(pollen_gene_df, ClusterName %in% input$select_cells_pollen)
        pollen_gene_df
    })
    
    
    clusters <- reactive({
        clusters <- as.character(unique(subset(gene_info, Gene == input$Gene)[,9]))
        clusters <- strsplit(clusters, ", ")
        clusters <- unlist(clusters)
    })
    
    groups <- reactive({
        groups <- as.character(unique(subset(gene_info, Gene == input$Gene)[,10]))
        groups <- strsplit(groups, ", ")
        groups <- unlist(groups)
    })
    
    #ASD DEG plots
    output$PlotASD_DEG_neuronal <- renderPlotly({
        plot = ggplot(ASD_DEGs_neuronal, aes(x=(FC), y=-log10(qvalue), colour=CellType, label=Gene))+
            geom_point()+scale_color_manual(values = neuronal_cols)+theme_bw()+
            theme(panel.border = element_blank(),panel.grid = element_blank(),
                  axis.line = element_line(colour = "black", size=2),
                  legend.title = element_blank())+
            geom_vline(xintercept = 0, colour="grey")+geom_hline(yintercept = 1.3, colour="grey")+
            xlab("log2(Fold Change)")+
            scale_x_continuous(limits = c(-1,1))+ggtitle("ASD-associated neuronal DEGs, \n n=458")
        ggplotly(plot, tooltip = c("Gene","CellType")) %>% layout(legend = list(x = 1, y = 1)) %>% 
            layout(legend=list(title=list(text='<b> Cell Type </b>')))
    })
    
    output$PlotASD_DEG_non_neuronal <- renderPlotly({
        plot = ggplot(ASD_DEGs_non_neuronal, aes(x=(FC), y=-log10(qvalue), colour=CellType, label=Gene))+
            geom_point()+scale_color_manual(values = non_neuronal_cols)+theme_bw()+
            theme(panel.border = element_blank(),panel.grid = element_blank(),
                  axis.line = element_line(colour = "black", size=2),
                  legend.title = element_blank())+
            geom_vline(xintercept = 0, colour="grey")+geom_hline(yintercept = 1.3, colour="grey")+
            xlab("log2(Fold Change)")+
            scale_x_continuous(limits = c(-1,1))+ggtitle("ASD-associated non-neuronal \n DEGs, n=234")
        ggplotly(plot, tooltip = c("Gene","CellType")) %>% layout(legend = list(x = 1, y = 1)) %>% 
            layout(legend=list(title=list(text='<b> Cell Type </b>')))
    })
    
    #Regional DEG plot BrainSpan
    
    index_x <- reactive({
        index_x <- as.numeric(brain.regions.DEGs[input$region_DEGs_xaxis])
    })
    
    index_y <- reactive({
        index_y <- as.numeric(brain.regions.DEGs[input$region_DEGs_yaxis])
    })
    
    x_name <- reactive({
        x_name = input$region_DEGs_xaxis
    })
    
    y_name <- reactive({
        y_name = input$region_DEGs_yaxis
    })
    
    region_DEGs_df <- reactive({
      null_all_data()
      load("BrainSpan DEGs.rda")
      
        if(input$region_DEGs_age_range == 1){
            region_DEGs_df <- early_brainspan_DEGs
        }
        
        if(input$region_DEGs_age_range == 2){
            region_DEGs_df <- mid_brainspan_DEGs
        }
        
        if(input$region_DEGs_age_range == 3){
            region_DEGs_df <- late_brainspan_DEGs
        }
        
        FDR <- input$region_DEGs_FDR
        
        region_DEGs_df <- region_DEGs_df[,c(1,index_x(),index_x()+1, index_y(), index_y()+1)]
        colnames(region_DEGs_df) <- c("Gene", "FC.x","FDR.x","FC.y","FDR.y")
        
        region_DEGs_df <- subset(region_DEGs_df, FDR.x<FDR | FDR.y<FDR)
        
        #region_DEGs_df$Signif = "No"
        region_DEGs_df$Signif[region_DEGs_df$FDR.x<FDR] = x_name()
        region_DEGs_df$Signif[region_DEGs_df$FDR.y<FDR] = y_name()
        region_DEGs_df$Signif[region_DEGs_df$FDR.y<FDR & region_DEGs_df$FDR.x<FDR] = "Both"
        
        region_DEGs_df$Signif <- factor(region_DEGs_df$Signif, levels = c(x_name(),y_name(),"Both"))
        region_DEGs_df
    })
    
    output$PlotDEGsRegion <- renderPlotly({
        plot=ggplot(region_DEGs_df(), aes(x=FC.x,y=FC.y,colour=Signif, label=Gene))+geom_point()+
            theme_black+
            geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
            scale_color_manual(values=c("steelblue","firebrick","forestgreen"))+
            xlab(paste0("log2FC (",x_name()," vs. Brain)"))+
            ylab(paste0("log2FC (",y_name()," vs. Brain)"))+
            theme(legend.title = element_blank(), legend.position = "right")
        
        ggplotly(plot, tooltip = c("Gene")) %>% layout(legend = list(x = 1, y = 1)) %>% 
            layout(legend=list(title=list(text='<b> Significant </b>')))
    })
    
    #Timecourse Plots - Primary
    
    allen_gene_df_timecourse <- reactive({
        
        allen_gene_df <- subset(allen_data, Symbol==input$Gene & name %in% input$select_regions_allen_timecourse)
        
    })
    
    output$PlotAllen_timecourse <- renderPlotly({
        validate(
            need((input$Gene %in% allen_data$Symbol), "Gene Expression data not avaliable"))
        if (input$plot_type_primary_timecourse == "Summary (mean)"){
            plot_type = stat_summary(fun="mean",geom="line", aes(group=Symbol),colour="black")
        } else if (input$plot_type_primary_timecourse == "Smooth (loess)"){
            plot_type = stat_smooth(method="loess", aes(group=Symbol))}
        else if (input$plot_type_primary_timecourse == "Line"){
            plot_type = geom_line(size=1)}
        else if (input$plot_type_primary_timecourse == "Line (smooth)"){
            plot_type = stat_smooth(method="loess")}
        
        plot <- ggplot(allen_gene_df_timecourse(), aes(x=days, y=expression_energy, colour=name, text=name))+
            theme_black+plot_type+
            geom_point(aes(fill=name), colour="black", size=2, shape=21)+
            scale_fill_brewer(palette = "Spectral")+scale_color_brewer(palette = "Spectral")+
            labs(y="Expression (ISH)", x="Age (days)")+#theme(legend.position = "bottom")+
            scale_x_continuous(breaks = c(11.5,13.5,15.5,18.5))+
            ggtitle("Allen Brain Atlas (mouse)")
        
        ggplotly(plot,tooltip = "text") #%>% layout(legend = list(orientation = "h", x = 0.4, y = -0.4))
    })
    
    brainspan_gene_df_timecourse <- reactive({
        brainspan_gene_df_timecourse <- subset(brainspan_gene_df(), Region %in% input$select_regions_brianspan_timecourse)
        brainspan_gene_df_timecourse$Group <- "A"
        brainspan_gene_df_timecourse
    })
    
    output$PlotBrainSpan_timecourse <- renderPlotly({
        
        if (input$plot_type_primary_timecourse == "Summary (mean)"){ #
            plot_type = stat_summary(fun="mean", geom="line", aes(group=Group),colour="black")
        } else if (input$plot_type_primary_timecourse == "Smooth (loess)"){
            plot_type = stat_smooth(method="loess", aes(group=Group))}
        else if (input$plot_type_primary_timecourse == "Line"){
            plot_type = geom_line(size=1)}
        else if (input$plot_type_primary_timecourse == "Line (smooth)"){ #
            plot_type = stat_smooth(method="loess")}
        
        plot <- ggplot(brainspan_gene_df_timecourse(), aes(x=age, y=value, fill=Region, colour=Region, text=Region))+
            geom_point(aes(fill=Region), colour="black", size=2, shape=21)+
            labs(y="Expression", x="Age (p.c.w)")+
            theme_black+plot_type+#theme(legend.position = "bottom")+
            scale_fill_brewer(palette="Spectral")+
            scale_colour_brewer(palette="Spectral")+
            ggtitle("BrainSpan (human)")
        
        ggplotly(plot,tooltip = "text")# %>% layout(legend = list(orientation = "h", x = 0.4, y = -0.4))
    })
    
    #Timecourse - Human Organoid
    pollen_function <- function(){
      pollen_org <- subset(pollen_gene_df(), CellSource == "in_vitro" & Organism == "Human" & Laboratory == "Pollen")
      y_limit_pollen = max(pollen_org$Expression)
      plot<- ggplot(pollen_org, aes(x=Age*7,y=Expression, colour=ClusterName, text=ClusterName, alpha=ClusterName))+
        stat_smooth(method = "loess",aes(group=Laboratory))+geom_point()+
        scale_color_manual(values = my_subcell_cols)+scale_alpha_manual(values = my_subcell_alpha)+
        xlab("Days")+ylab("Expression")+theme_black+
        ggtitle("Pollen et al.,")+ theme(plot.title = element_text(size = 14))+
        scale_y_continuous(limits=c(0,y_limit_pollen))
    }
    
    output$PlotOrganoidPollen <- renderPlotly({
      plot = pollen_function()
      ggplotly(plot, tooltip = "text")
    })
    
    benito_function <- function(){
      benito_org <- subset(benito_gene_df(), Experiment == "Benito-Kwiecinski-Human")
      benito_org$Batch <- c(rep("Batch 1",7),rep("Batch 2",7),rep("Batch 3",7))
      ggplot(benito_org, aes(x=Day,y=value, colour=Batch, text=Batch))+scale_color_brewer(palette = "Set1")+
        stat_smooth(method = "loess",aes(group=Model))+
        geom_line(size=2)+
        xlab("Days")+ylab("Expression")+theme_black + ggtitle("Benito-Kwiecinski et al.,")+ theme(plot.title = element_text(size = 14))
    }
    
    output$PlotOrganoidBenito <- renderPlotly({
      plot = benito_function()
      ggplotly(plot, tooltip = 'text')
    })
    
    kanton_function <- function(){
      kanton_org <- subset(kanton_gene(), Species == "Human")
      y_limit_kanton = max(kanton_org$value)
      
      plot <-ggplot(kanton_org, aes(x=Day,y=value, colour=CellType, text=CellType, alpha=CellType))+
        stat_smooth(method = "loess",aes(group=Species))+
        geom_point()+scale_alpha_manual(values = my_subcell_alpha)+
        xlab("Days")+ylab("Expression")+theme_black+
        ggtitle("Kanton et al.,")+ theme(plot.title = element_text(size = 14))+
        scale_y_continuous(limits=c(0,y_limit_kanton))+
        scale_color_manual(values = my_subcell_cols)
    }
    
    output$PlotOrganoidKanton <- renderPlotly({
      plot = kanton_function()
      ggplotly(plot, tooltip = 'text')
    })
    
    sloan_function <- function(){
      sloan_org <- subset(pollen_gene_df(), CellSource == "in_vitro" & Organism == "Human" & Laboratory == "Pasca")
      y_limit_sloan = max(sloan_org$Expression)
      
      plot <- ggplot(sloan_org, aes(x=Age*7,y=Expression, colour=ClusterName, text=ClusterName, alpha=ClusterName))+
        stat_smooth(method = "loess",aes(group=Laboratory))+geom_point()+
        scale_alpha_manual(values = my_subcell_alpha)+scale_color_manual(values = my_subcell_cols)+
        xlab("Days")+ylab("Expression")+theme_black+
        ggtitle("Sloan et al.,")+ theme(plot.title = element_text(size = 14))+
        scale_y_continuous(limits=c(0,y_limit_sloan))
    }
    
    output$PlotOrganoidSloan <- renderPlotly({
      plot = sloan_function()
      ggplotly(plot, tooltip = 'text')
    })
    
    camp_function <- function(){
      treutlein_org <- subset(pollen_gene_df(), CellSource == "in_vitro" & Organism == "Human" & Laboratory == "Treutlein")
      treutlein_org$Age <- round(treutlein_org$Age*7,0)
      y_limit_camp = max(treutlein_org$Expression)
      
      plot <-ggplot(treutlein_org, aes(x=Age,y=Expression))+
        geom_point()+
        stat_smooth(method = "loess",aes(group=Laboratory))+
        xlab("Days")+ylab("Expression")+theme_black+
        ggtitle("Camp et al.,")+ theme(plot.title = element_text(size = 14))+
        scale_y_continuous(limits=c(0,y_limit_camp))
    }
    
    output$PlotOrganoidCamp <- renderPlotly({
      plot = camp_function()
      ggplotly(plot)
    })
    
    #Comparison Plot
    comparison_function <- function(){
      if(input$Facet != "none"){
        facetwrap = facet_wrap(input$Facet, scales = "free_y")
      } else facetwrap = NULL
      
      plot <- ggplot(comparison_df(), aes(x=Day, y=log2(value+1), group=Experiment, Model=Model, Species=Species))+
        stat_smooth(se=F,aes_string(colour=input$ColourBy1))+theme_black+
        ylab("Expression")+#theme(legend.position="right")+
        scale_color_brewer(palette=input$ColourPal)+
        coord_cartesian(xlim = c(input$slider_comparison[1], input$slider_comparison[2]))+
        facetwrap
    }
    
    output$PlotComparison <- renderPlotly({
      plot = comparison_function()
        ggplotly(plot, tooltip = c("Experiment","Model","Organism","Species")) %>% layout(legend = list(orientation = "h", x = 0.4, y = -0.4))
    })
    
    brainspan_plot <- function(){
        ggplot(brainspan_gene_df(), aes(x=age, y=value, colour=structure_acronym, group=structure_acronym))+
            geom_point(aes(fill=structure_acronym), colour="black", size=2, shape=21)+
            geom_line(size=1)+
            facet_wrap(~Region, scales="free_x",  ncol = 3)+
            labs(y="Expression", x="Age (p.c.w)")+
            theme_black+
            scale_fill_manual(values=structure_cols)+scale_color_manual(values=structure_cols)+
            #directlabels::geom_dl(aes(label = structure_acronym), method ="smart.grid", cex = 0.5)+
            labs(title = paste0("Regional ", input$Gene, " expression in humans")) + theme(plot.title = element_text(size = 14))
    }
    
    HDBR_plot <- function(){
      ggplot(HDBR_gene_df(), aes(x=Age, y=value, fill=organism_part, colour=organism_part, group=organism_part))+
        geom_point(aes(fill=organism_part), colour="black", size=2, shape=21)+
        stat_smooth(fun="mean",geom = "line", colour="black",aes(group=Region))+
        coord_cartesian(ylim = c(0,max(HDBR_gene_df()$value)))+
        facet_wrap(~Region,  ncol = 3)+
        labs(y="Expression", x="Age (p.c.w)")+
        theme_black+
        theme(legend.position = "none")+
        scale_fill_manual(values=structure_cols_HDBR)+scale_color_manual(values=structure_cols_HDBR)+
        labs(title = paste0("Regional ",input$Gene, " expression in humans")) + theme(plot.title = element_text(size = 14))
    }
    
    plotInoue <- function(){
        title = "Inoue"
        data <- subset(NPC_df(), Experiment==title)
        timecourse_plot(data, title, inoue_cols)
    }
    
    plotLIBD <- function(){
        title = "LIBD"
        data <- subset(NPC_df(), Experiment==title)
        timecourse_plot(data, title, LIBD_cols)
    }
    
    plotCORTECON <- function(){
        title = "CORTECON"
        data <- subset(NPC_df(), Experiment==title)
        timecourse_plot(data, title, cortecon_cols)
    }
    
    plotMicali <- function(){
        title = "Micali"
        data <- subset(NPC_df(), Experiment==title)
        timecourse_plot(data, title, micali_cols)
    }
    
    plotAllen <- function(){
        validate(
            need((input$Gene %in% allen_data$Symbol), "Gene Expression data not avaliable")
        )
        
        allen_gene_df <- subset(allen_data, Symbol==input$Gene)
        
        ggplot(allen_gene_df, aes(x=days, y=expression_energy, colour=name, group=name))+
            geom_point(aes(fill=name), colour="black", size=2, shape=21)+
            geom_line(size=1)+
            facet_wrap(~acronym, scales="free_x",  ncol = 3)+
            labs(y="Expression (ISH)", x="Age (days)")+
            scale_x_continuous(breaks = c(11.5,13.5,15.5,18.5))+ 
          scale_y_continuous(breaks = c(0, 5, 10,15))+
            theme_black+
            scale_fill_brewer(palette = "Spectral")+scale_color_brewer(palette = "Spectral")+
            labs(title = paste0("Regional ", input$Gene, " expression in mice")) + theme(plot.title = element_text(size=14))
            
    }
    
    
    user_gene_list <- reactive({
       inFile <- input$user_data
        
        if (is.null(inFile))
            return(NULL)
        
       user_gene_list <- read.table(inFile$datapath, header = F, sep="\t")[,1]
    })
    
    #Make Gene Info Table
    output$mytable <- renderDataTable({
        
      #If no GO term selected skip GO filter
      if(!(is.null(input$go_terms))){
        for (i in 1:length(input$go_terms)) {
          select_go_id <- strsplit(input$go_terms[i]," - ")[[1]][1]
          go_genes <- subset(all_go_genes, go_id == select_go_id)[,1]
          gene_info = subset(gene_info, Gene %in% go_genes)
          }
        }
      if(!(is.null(input$region_filter))){
        region_genes = brainspan_regional_data[grepl(input$region_filter,brainspan_regional_data$`Regional Expression (BrainSpan)`),1]
        gene_info <- subset(gene_info, Gene %in% region_genes)
      }
      if(!(is.null(input$select_celltype_fan))){
        for (i in 1:length(input$select_celltype_fan)) {
          cluster_genes = gene_info[grepl(input$select_celltype_fan[i],gene_info$`Cell Type (Fan et al.,)`),1]
          gene_info = subset(gene_info, Gene %in% cluster_genes)
        }
      }
      if(!(is.null(input$select_cellstate_bhaduri))){
        for (i in 1:length(input$select_cellstate_bhaduri)) {
          cluster_genes = gene_info[grepl(input$select_cellstate_bhaduri[i],gene_info$`State (Bhaduri et al.,)`),1]
          gene_info = subset(gene_info, Gene %in% cluster_genes)
        }        }
      if(!(is.null(input$select_celltype_bhaduri))){
        for (i in 1:length(input$select_celltype_bhaduri)) {
          cluster_genes = gene_info[grepl(input$select_celltype_bhaduri[i],gene_info$`Cell Type (Bhaduri et al.,)`),1]
          gene_info = subset(gene_info, Gene %in% cluster_genes)
        }
      }
      if(!(is.null(input$select_celltype_nowakowski))){
        for (i in 1:length(input$select_celltype_nowakowski)) {
          cluster_genes = gene_info[grepl(input$select_celltype_nowakowski[i],gene_info$`Cell Type (Nowakowski et al.,)`),1]
          gene_info = subset(gene_info, Gene %in% cluster_genes)
        }
      }
      if(!(is.null(input$select_cellsubtype_nowakowski))){
        for (i in 1:length(input$select_cellsubtype_nowakowski)) {
          cluster_genes = gene_info[grepl(input$select_cellsubtype_nowakowski[i],gene_info$`Cell subtypes (Nowakowski et al.,)`),1]
          gene_info = subset(gene_info, Gene %in% cluster_genes)
        }
      }
      if(input$human_specific==TRUE){
        gene_info = subset(gene_info, Gene %in% unlist(gene_lists["Human specific regulation"]))
      }
      if(input$SFARI_gene==TRUE){
        gene_info = subset(gene_info, Gene %in% unlist(gene_lists["SFARI (ASD genes)"]))
      }
      if(input$user_gene_list_select==TRUE){
        gene_info = subset(gene_info, Gene %in% user_gene_list())
      }
      
      columns_to_show = gene_info_column_lists[input$show_vars]
      columns_to_show = c("Gene",as.character(unlist(columns_to_show)))
      gene_info <- unique(gene_info[,columns_to_show,drop = FALSE])
      datatable(gene_info,options = list(scrollX=TRUE))

    })
    
    #Download Gene Table
    output$downloadGeneTable <-downloadHandler(
      filename = function(){"Gene Table.xlsx"},
      content = function(file){
        writexl::write_xlsx(gene_info, path = file)
      })
    
    #Enrichment Analysis
    
    #Create excel workbook to download 
    output$downloadEnrichment <- downloadHandler(
      filename = function(){"Enrichment Results.xlsx"},
      content = function(file) {
        if(input$enrichment_list1 == "User list"){
          user_file_t <- input$user_gene_lists
          user_file_text <- read.csv(user_file_t$datapath)
          list1 = user_file_text
          list1_genes <- unique(as.character(unlist(list1)))
        } else {
          list1 = input$enrichment_list1
          list1_genes <- unique(as.character(unlist(gene_lists[list1])))
        }
        list2 = input$enrichment_list2
        list3 = input$enrichment_background_list
        
        list2_genes <- unique(as.character(unlist(gene_lists[list2])))
        list3_genes <- unique(as.character(unlist(gene_lists[list3])))
        
        go.obj <- newGeneOverlap(list1_genes,list2_genes,genome.size=length(list3_genes))
        go.obj <- testGeneOverlap(go.obj)
        
        str1 <- paste("List 1 size =", length(go.obj@listA))
        str2 <- paste("List 2 size =", length(go.obj@listB))
        str3 <- paste("Background gene list size =",go.obj@genome.size)
        str4 <- paste("Odds Ratio =", round(go.obj@odds.ratio,digits=2))
        str5 <- paste("P-value =", formatC(go.obj@pval, format = "e", digits = 2))
        str6 <- paste("Jaccard Index =", round(go.obj@Jaccard,digits=2))
        str7 <- paste("Background Gene list -", list3)
        str8 <- as.data.frame(getIntersection(go.obj))
        colnames(str8)=paste("n=",length(getIntersection(go.obj)),"genes")
        list1_df <- as.data.frame(list1_genes)
        colnames(list1_df)=paste("n=",length(list1_genes),"genes")
        list2_df <- as.data.frame(list2_genes)
        colnames(list2_df)=paste("n=",length(list2_genes),"genes")
        list3_df <- as.data.frame(list3_genes)
        colnames(list3_df)=paste("n=",length(list3_genes),"genes")
        
        tmp <- as.data.frame(c(str1,str2,str3,str4,str5,str6, str7))
        colnames(tmp) <- paste("Overlap between ",list1,"vs.",list2)
        
        sheets <- list(tmp, str8, list1_df, list2_df, list3_df)
        names(sheets) <- c("Summary", "Intersection Gene List", paste(list1,"genes"),paste(list2,"genes"),paste(list3,"genes"))
        writexl::write_xlsx(sheets, path = file)
      },
      contentType="application/xlsx" )
    
    output$downloadMultiVennOverlap <- downloadHandler(
      filename = function(){"Venn Overlap Results.xlsx"},
      content = function(file) {
        intersection_df = get.venn.partitions(x = gene_lists[input$multiple_gene_list])
        tmp <- data.frame(ID=as.numeric(),Genes=as.character())
        for(i in 1:length(intersection_df$..values..)){
          lst <- intersection_df$..values..[i]
          lst <- paste(lst[[1]], collapse = ", ")
          tmp2 <- data.frame(ID=i,Genes=lst)
          tmp <- rbind(tmp,tmp2)
        }
        
        venn_lists <- intersection_df$..values..
        venn_lists <- lapply(venn_lists, FUN=as.data.frame)
        intersection_df = intersection_df[,-5]
        intersection_df$ID = rownames(intersection_df)
        intersection_df = merge(intersection_df, tmp, by="ID", by.y="ID")
        sheets <- c(list(intersection_df),venn_lists)
        names(sheets) <- c("Summary", c(1:length(venn_lists)))
        writexl::write_xlsx(sheets, path = file)
      },
      contentType="application/xlsx")
    
    #Output text summary results
    output$enrichment_text <- renderText({
      if(input$enrichment_list1 == "User list"){
        user_file_t <- input$user_gene_lists
        user_file_text <- read.csv(user_file_t$datapath)
        go.obj <- newGeneOverlap(as.character(unlist(user_file_text)), 
                                 as.character(unlist(gene_lists[input$enrichment_list2])), 
                                 genome.size=length(as.character(unlist(gene_lists[input$enrichment_background_list]))))
      } else {
        go.obj <- newGeneOverlap(as.character(unlist(gene_lists[input$enrichment_list1])),
                               as.character(unlist(gene_lists[input$enrichment_list2])),
                               genome.size=length(as.character(unlist(gene_lists[input$enrichment_background_list]))))
      }
      
      go.obj <- testGeneOverlap(go.obj)
      str1 <- paste("List 1 size =", length(go.obj@listA))
      str2 <- paste("List 2 size =", length(go.obj@listB))
      str3 <- paste("Background gene list size =",go.obj@genome.size)
      str4 <- paste("Odds Ratio =", round(go.obj@odds.ratio,digits=2))
      str5 <- paste("P-value =", formatC(go.obj@pval, format = "e", digits = 2))
      str6 <- paste("Jaccard Index =", round(go.obj@Jaccard,digits=2))
      
      HTML(paste(str1, str2, str3, str4, str5, str6, sep = '<br/>'))
    })
    
    #Enrichment Venn Diagram - add function for multiple venn diagram overlap
    observeEvent(input$user_gene_lists, {
      updateSelectInput(
        session,
        "enrichment_list1",
        choices = c("User list", names(gene_lists)))
      #user_file <- input$user_gene_lists
    })
    
    output$venn_diagram <- renderPlot({
      if(input$enrichment_list1 == "User list"){
        user_file <- input$user_gene_lists
        venn_list_1 <- read.csv(user_file$datapath)
        venn_list_1 <- as.character(unlist(venn_list_1))
        venn_list_1 <-venn_list_1[!is.na(venn_list_1)]
      } else {
        venn_list_1 <-as.character(unlist(gene_lists[input$enrichment_list1]))
      venn_list_1 <-venn_list_1[!is.na(venn_list_1)]
      }
      venn_list_2 <-as.character(unlist(gene_lists[input$enrichment_list2]))
      venn_list_2 <-venn_list_2[!is.na(venn_list_2)]
      venn_list <- list(venn_list_1, venn_list_2)
      names(venn_list) <- c(input$enrichment_list1, input$enrichment_list2)
      
      display_venn(x = venn_list,
                   # Circles
                   fill = c("#009E73", "#56B4E9"),
                   lwd = 2,
                   lty = 'blank',
                   # Numbers
                   cex = .9,
                   fontface = "italic",
                   # Set names
                   cat.cex = 1,
                   cat.fontface = "bold",
                   cat.default.pos = "outer",
                   cat.dist = c(0, 0))
    })
    
    #Multiple gene lists overlap
    output$upsetPlot <- renderPlot({
      UpSetR::upset(UpSetR::fromList(gene_lists[input$multiple_gene_list]),order.by = "freq")
    })
    
    output$multipleVenn <- renderPlot({
      if(length(input$multiple_gene_list)==2){
        venn_cols <- c("#009E73", "#56B4E9")
      }else
        venn_cols <- brewer.pal(length(input$multiple_gene_list),"Spectral")
      display_venn(x = gene_lists[input$multiple_gene_list],
                   # Circles
                   fill = venn_cols,
                   lwd = 2,
                   lty = 'blank',
                   # Numbers
                   cex = .9,
                   fontface = "italic",
                   # Set names
                   cat.cex = 1,
                   cat.fontface = "bold",
                   cat.default.pos = "outer")
    })
    
    #Download scRNA-seq cluster Test Data
    output$downloadTest_scRNAmarkers <- downloadHandler(
      filename = function(){"scRNAseq_test_markers.txt"},
      content = function(file) {
        write.table(datasetInput(), file, row.names = FALSE)
      }
    )
    
    #help box in cluster overlap page
    observeEvent(input$help_cluster_overlap, {
      shinyalert(
        title = "Cluster Overlap Enrichement Tool",
        text = HTML("This tool helps to classify your cell clusters from a scRNA-seq experiment.
        Upload your table of differentially expressed genes within your clusters and set your filters
        for genes within your cluster (Log2 Fold Change and p value). The tool will overlap each cluster gene list with
        cell type gene lists from published scRNA-seq datasets of fetal human brain. 
        Results table gives the Fisher's Enrichment test results and Summary table gives the top cell type classification (by Enrichment) for each dataset."),
        size = "s", 
        closeOnEsc = TRUE,
        closeOnClickOutside = TRUE,
        html = TRUE,
        type = "info",
        showConfirmButton = TRUE,
        showCancelButton = FALSE,
        confirmButtonText = "Close",
        confirmButtonCol = "blue"
      )
    })
    
    #Cluster Identification Overlap Function
    user_cluster_data <- reactive({
      inFile <- input$user_clusters_filename
      
      if (is.null(inFile))
        return(NULL)
      
      user_clusters <- read.table(inFile$datapath, header = T, sep="\t")
      user_clusters <- subset(user_clusters, avg_log2FC>input$cluster_FC_filter)
      if(input$stats_filter=="Power"){
        user_clusters <- subset(user_clusters, power>input$cluster_power_filter)
      }else
      if(input$stats_filter=="P-value"){
        user_clusters <- subset(user_clusters, p_val<input$cluster_pvalue_filter)
      }
    })
    
   # user_cluster_data <- reactive({
  #    if(input$test_scRNAseq_cluster == TRUE)
  #    {
  #      user_clusters <- read.table("scRNAseq_test_markers.txt", header = T, sep="\t")
  #      user_clusters <- subset(user_clusters, p_val<input$cluster_pvalue_filter)
  #      user_clusters <- subset(user_clusters, avg_log2FC>input$cluster_FC_filter)
  #    }
  #  })
    
    observeEvent(input$run_cluster_analysis,{
      user_cluster_names <- as.character(unique(user_cluster_data()$cluster))
      #Make user cluster list
      user_cluster_list <- list()
      for(i in user_cluster_names){
        user_gene_list <- subset(user_cluster_data(), cluster== i)
        user_gene_list <- unique(user_gene_list$gene)
        user_gene_list <- list(user_gene_list)
        names(user_gene_list) <- i
        user_cluster_list <- c(user_cluster_list,user_gene_list)
      }
      
      
      #Comparing user lists to all loaded lists
      all_gene_lists <- c(gene_lists, user_cluster_list)

      background_list <- gene_lists[1][[1]]
      background_size=length(background_list)
      
      #Using only scRNA-seq fetal human datasets
      #Fan, Shi, Nowakowski, Zhong, Polioudakis, Brainspan 
      my_lists <- names(gene_lists)[16:68]
      
      fishers_cluster_overlap <- data.frame(Cluster=as.character(),Dataset=as.character(),pvalue=as.numeric(),
                                            Overlap=as.numeric(),Expected=as.numeric(),Enrichment=as.numeric())
      for(i in user_cluster_names){
        listA <- all_gene_lists[i][[1]]
        for(j in my_lists){
          if(i!=j){
            listB <- all_gene_lists[j][[1]]
            x=length(subset(listA,listA %in% listB))
            m=length(listA)
            k=length(listB)
            n=background_size-m
            expected <- round((m*k)/background_size,2)
            enrichment = round(x/expected,2)
            pval = dhyper(x, m, n, k)
            tmp <- data.frame(Cluster =i, Dataset=j, pvalue=pval,Overlap=x,Expected=expected,Enrichment = enrichment)
            fishers_cluster_overlap <- rbind(fishers_cluster_overlap, tmp) 
          }
        }
      }
      
      fishers_cluster_overlap$Enrichment <- round(fishers_cluster_overlap$Enrichment,2)
      fishers_cluster_overlap$pvalue <- signif(fishers_cluster_overlap$pvalue,digits=3)
      fishers_cluster_overlap <- fishers_cluster_overlap[order(fishers_cluster_overlap$Enrichment,decreasing = TRUE),]
      
      .render_fishers_cluster_overlap_table(output = output, data = fishers_cluster_overlap,filename="Cluster Results")
      
      output$uploaded_cell_clusters<- renderDataTable({
        datatable(user_cluster_data())
      })
      
      summary_clusters <- data.frame(Cluster=as.character(),Fan=as.character(),Shi=as.character(),
                                     Polioudakis=as.character(),Nowakowski=as.character(),BrainSpan=as.character())
      for (i in user_cluster_names) {
        tmp <- subset(fishers_cluster_overlap, Cluster==i & Overlap>3 & pvalue<0.0001)
        fan <- tmp[grepl("Fan et al",tmp$Dataset),][1,2]
        fan <- str_split(fan," genes")[[1]][1]
        shi <- tmp[grepl("Shi et al",tmp$Dataset),][1,2]
        shi <- str_split(shi," genes")[[1]][1]
        zhong <- tmp[grepl("Zhong et al",tmp$Dataset),][1,2]
        zhong <- str_split(zhong," genes")[[1]][1]
        polioudakis <- tmp[grepl("Polioudakis et al",tmp$Dataset),][1,2]
        polioudakis <- str_split(polioudakis," genes")[[1]][1]
        nowakowski <- tmp[grepl("Nowakowski et al",tmp$Dataset),][1,2]
        nowakowski <- str_split(nowakowski," genes")[[1]][1]
        brainspan <- tmp[grepl("BrainSpan",tmp$Dataset),][1,2]
        brainspan <- str_split(brainspan," genes")[[1]][1]
        
        tmp <- data.frame(Cluster=i,Fan=fan,Shi=shi,Zhong=zhong,Polioudakis=polioudakis,
                          Nowakowski=nowakowski,BrainSpan=brainspan)
        summary_clusters <- rbind(summary_clusters,tmp)
      }

      .render_summary_cluster_table(output = output, data = summary_clusters,filename="Summary Cluster")

    })
    
    
    output$downloadCellClusterResults <- downloadHandler(
      filename = function(){"Cell Cluster Classification Results.xlsx"},
      content = function(file) {
        writexl::write_xlsx(summary_clusters, path = file)
      },
      contentType="application/xlsx")
    

    #PSC Neural Differentiation Plot
    output$PlotInoue <- renderPlotly({
        plot = plotInoue()
        ggplotly(plot, tooltip = c("Day", "value", "group"))
    })
    
    output$PlotLIBD <- renderPlotly({
        plot = plotLIBD()
        ggplotly(plot, tooltip= c("Day", "value", "group"))
    })
    
    output$PlotCORTECON <- renderPlotly({
        plot= plotCORTECON()
      ggplotly(plot, tooltip= c("Day", "value", "group"))
    })
    
    output$PlotMicali <- renderPlotly({
        plot = plotMicali()
      ggplotly(plot, tooltip= c("Day", "value", "group"))
    })
    
    #BrainSpan Plot By Region
    output$PlotBrainSpan <- renderPlotly({
        plot = brainspan_plot()
      ggplotly(plot, tooltip = c("x","group"))
    })
    
    #HDBR Plot
    output$PlotHDBR <- renderPlotly({
      plot=HDBR_plot()
      ggplotly(plot, tooltip = c("x","group"))
    })
    
    #AllenBrain Atlas Plot By Region
    output$PlotAllen <- renderPlotly({
        plot=plotAllen()
        ggplotly(plot, tooltip = c("x","group"))
    })
    
    #Fan scRNAseq Data
    output$PlotFan_tsne <- renderPlot({
        plotFan_tsne()
    })
    
    output$PlotFan_violin <- renderPlot({
        plotFan_violin()
    })
    
  
    #Combined scRNAseq Differential Expression PLot
    output$CellColLedgend <- renderImage({
      return(list(src = "CellTypeLedgend.png",contentType = "image/png",alt = "Alignment",
                  width = "15%"))
    }, deleteFile = FALSE) #where the src is wherever you have the picture
    
    plotCombinedCell <- function(){
      gene_cell_plot_data <- subset(cell_plot_data, Gene==input$Gene)
      plot_limits=round(max(sqrt(gene_cell_plot_data$Fold.Change^2)),0)+1
      
      ggplot(gene_cell_plot_data, aes(x=Fold.Change,y=-(log10(Pval)), fill=`Cell Type`, colour=`Cell Type`, shape=Dataset, label=Cluster))+
        geom_point(size=2)+
        scale_fill_manual(values=my_cell_cols)+
        scale_shape_manual(values=dataset_shapes)+
        scale_color_manual(values=my_cell_cols)+theme_bw()+
        theme(panel.border = element_blank(),panel.grid = element_blank(),
              axis.line = element_line(colour = "black", size=1),
              legend.title = element_blank(),
              legend.position = "none")+
        geom_vline(xintercept = 0, colour="grey")+scale_x_continuous(limits = c(-plot_limits,plot_limits))+
        ylab("-Log10(P val)")+xlab("log2(Fold Change)")
    }
    
    output$PlotCombinedCell <- renderPlotly({
      plot = plotCombinedCell()
      plotly::ggplotly(plot,tooltip = c("Dataset","Cluster"))
    })
    
    output$PlotSummaryCell <- renderPlotly({
      plot=ggplot(summary_scRNAseq_studies, aes(x=Study, y=Percentage, fill=Cell_Type))+
        scale_fill_manual(values=my_cell_cols)+theme_black+
        coord_flip()+
        geom_col(width = 0.75, colour="white")+
        scale_y_continuous(labels=scales::percent)+
        labs(fill="Cell Type")
      plotly::ggplotly(plot,tooltip = c("Dataset","Cell_Type"))
      
    })
    
    output$DatasetSummaryTable <- renderTable({
        Key <- c(as.character(icon("circle")),  as.character("&#9650"),
                 as.character(icon("square")),  as.character(icon("plus")), as.character("&#9660"))
        dataset_table <- cbind(Key, dataset_table)
        dataset_col_names <- colnames(dataset_table)
        dataset_table <- as.data.frame(lapply(dataset_table, as.character)) #Convert to character to round numbers
        colnames(dataset_table) <- dataset_col_names
        xtable::xtable(dataset_table) 
      },sanitize.text.function = function(x)x)
    
    
    #Evolution Tab
    
    #Kanton Data
    kanton_gene_plot <- reactive({
      kanton_gene_plot <- subset(kanton_gene(), CellType %in% input$select_cells_kanton)
      kanton_gene_plot
    })
    
    output$PlotKanton_timecourse <- renderPlotly({
      plot = ggplot(kanton_gene_plot(), aes(x=Day, y=log2(value+1), colour=CellType, text=CellType,alpha=CellType))+
        scale_alpha_manual(values = my_subcell_alpha)+
        geom_point(position="jitter")+
        stat_smooth(aes(group=Species), method="loess", colour="black")+
        facet_wrap(~Species)+
        theme_black+#theme(legend.position = "bottom")+
        ylab("log2(Expression)")+
        scale_color_manual(values=my_subcell_cols)
      
      ggplotly(plot,tooltip = "text") %>% layout(legend = list(orientation = "h", x = 0.4, y = -0.4))
    })
    
    output$PlotKanton_violin <- renderPlotly({
      plot = ggplot(kanton_gene_plot(), aes(x=Species, y=value))+
        scale_color_manual(values=c("red","black"))+
        theme_bw()+geom_violin(aes(colour=Species))+
        facet_wrap(~CellType, labeller = labeller(CellType = label_wrap_gen(20)))+
        scale_fill_gradient2(low = "blue", mid="orange",high = "red",  midpoint = 60)+
        geom_point(aes(fill=Day), shape=21)+
        theme(text = element_text(size = 10))
      
      ggplotly(plot)
    })
    
    
    #Pollen et al 2019 data
    output$PollenPlot_timecourse <- renderPlotly({
        
        plot = ggplot(pollen_gene_df(), aes(x=Age, y=Expression, group=Organism, colour=ClusterName, text=ClusterName, alpha=ClusterName))+
          geom_point(position="jitter")+
          stat_smooth(aes(group=Organism), method="loess", colour="black")+
          facet_wrap(~Organism, scales = "free_x")+
          theme_black+#theme(legend.position = "bottom")+
          scale_alpha_manual(values = my_subcell_alpha)+
          ylab("log2(Expression)")+
          scale_y_continuous(limits = c(0,12))+
          scale_color_manual(values=my_subcell_cols)
        
        ggplotly(plot,tooltip = "text") %>% layout(legend = list(orientation = "h", x = 0.4, y = -0.4))
        
    })
    
    output$PollenPlot_violin <- renderPlotly({
        plot = ggplot(pollen_gene_df(), aes(x=Organism, y=Expression))+
            scale_color_brewer(palette="Set1")+
            theme_bw()+geom_violin(aes(colour=Organism))+
            facet_wrap(~ClusterName, labeller = labeller(ClusterName = label_wrap_gen(20)))+
            scale_fill_gradient2(low = "blue", mid="orange",high = "red",  midpoint = 30)+
            geom_point(aes(fill=Age), shape=21)+
            theme(text = element_text(size = 10))
        
        ggplotly(plot)
        
    })
    
    pollen_module_genes <- reactive({
        select_modules <- subset(pollen_modules, Gene == input$Gene)
        pollen_modules_genes <- subset(pollen_modules,pollen_modules$Module.ID %in% select_modules$Module.ID)[,1]
    })

    ##Gene Module Analysis
    module_table <- reactive({
        pollen_modules %>%
            filter(Gene == input$Gene) %>%
            merge(pollen_module_info, by.x="Module.ID", by.y="Network", all.x=T)
    })
    
    output$module_table <- renderDataTable({
        datatable(module_table())
    })
    
    observe({
      updateSelectizeInput(session, 'module', choices = unique(subset(pollen_modules, Gene == input$Gene))[,3])
    })
    
   # output$select_module <- renderUI({
  #      choice_module <- reactive({
  #          pollen_modules %>%
  #              filter(Gene == input$Gene) %>%
  #              merge(pollen_module_info, by.x="Module.ID", by.y="Network") %>%
  #              dplyr::pull(Module.ID) %>%
  #              as.character()
  #      })
  #      selectizeInput('module', 'Select module to plot PPI Network;', choices = c("select" = "", choice_module()))
  #  })
    
    module_network <- reactive({
        modules_genes <- subset(pollen_modules,pollen_modules$Module.ID == input$module)[,1]
        module_network <- subset(stringDB_data, Gene1 %in% modules_genes & Gene2 %in% modules_genes)
    })
    
    #my_module= "primary.macaque.ME.antiquewhite1"
    
    module_gene_table <- reactive({
        modules_genes <- subset(pollen_modules,pollen_modules$Module.ID == input$module)[,1]
        gene_info_slim <- unique(gene_info[,c("Gene","Gene Score","Cell Type (Nowakowski et al.,)","Human Specific Regulation Dataset","Notes")])
        module_gene_table <- subset(gene_info_slim, Gene %in% modules_genes)
        colnames(module_gene_table) <- c("Gene","SFARI.Gene.Score","Cell.Type","Human regulation (dataset)","Notes")
        #ene_info_modules <- dplyr::filter(genecard_description_summary, genecard_description_summary$gene %in% shQuote(modules_genes, type = c("cmd")))[,1:4]
        #gene_info_modules$gene <- noquote(gene_info_modules$gene)
        #module_gene_table <- merge(module_gene_table, gene_info_modules, by.x="Gene",by.y="gene")
        module_gene_table
    })
    
    output$module_gene_table <- renderDataTable({
        datatable(module_gene_table())
    })
    
    module_TF_table <- reactive({
        module_genes <- subset(pollen_modules, Module.ID==input$module)[,1]
        module_TF_table <- subset(ttrust, Target %in% module_genes)
        module_TF_table
    })
    
    output$module_TF_table <- renderDataTable({
        datatable(module_TF_table())
    })
    
    output$top_TFs <- renderDataTable({
        datatable(top_TFs)
    })
    
    output$top_TFs_plot <- renderPlotly({
        plot <- ggplot(TF_plot, aes(x=FDR,y=OddsRatio, label=Module.ID, fill=TF,size=nTargets,
                                    Interpretation_notes=Interpretation_notes, shape=Module.Type))+
            scale_x_log10()+ylab("Odds Ratio")+
            scale_shape_manual(values = c(21,24))+
            geom_point(aes(fill=TF), colour="black")+
            ggtitle("TFs regulating modules")+theme_black+scale_y_log10()
        
        ggplotly(plot, tooltip = c("Module.ID","Interpretation_notes","TF"))
        
    })
    
    module_TF_table_fishers <- reactive({
        module_TF_table_fishers <- subset(module_TF_fishers_df, Module.ID==input$module)
        module_TF_table_fishers <- module_TF_table_fishers[,-1]
        module_TF_table_fishers
    })
    
    output$module_TF_table_fishers <- renderDataTable({
        datatable(module_TF_table_fishers(),
                  caption="Transcription Factors whose targets are enriched for genes within the gene module.")
    })
    
    output$module_network_plot <- renderVisNetwork({
      validate(
        need((input$module != ""), "Select module to plot from dropdown list"))
        network = module_network()
        edges <- data.frame(from = network$Gene1, to = network$Gene2)
        nodes_genes <- unique(c(network$Gene1, network$Gene2))
        nodes <- data.frame(id = nodes_genes, label = nodes_genes)
        
        ##Set colour##
        nodes$group <- "none"
        nodes$group[nodes$id %in% Signalling]= "Signalling"
        nodes$group[nodes$id %in% Cellular_Respiration]= "Cellular Respiration"
        nodes$group[nodes$id %in% TF]= "TF"
        nodes$group[nodes$id %in% Receptor]= "Receptor"
        nodes$group[nodes$id %in% ECM]= "ECM"
        nodes$group[nodes$id %in% Glycolysis]= "Glycolysis"
        nodes$group[nodes$id %in% Synapse]= "Synapse"
        nodes$group[nodes$id %in% Neurotransmitter]= "Neurotransmitter"
        nodes$group[nodes$id %in% Voltage_Channel]= "Voltage Channel"
        
        visNetwork(nodes,edges) %>%
            
            visGroups(groupname = "Signalling", color = "violet") %>%
            visGroups(groupname = "Cellular Respiration", color = "cornflowerblue") %>%
            visGroups(groupname = "TF", color = "orange") %>%
            visGroups(groupname = "Receptor", color = "pink") %>%
            visGroups(groupname = "ECM", color = "brown") %>%
            visGroups(groupname = "Glycolysis", color = "cyan") %>%
            visGroups(groupname = "Synapse", color = "red") %>%
            visGroups(groupname = "Neurotransmitter", color = "darkolivegreen1") %>%
            visGroups(groupname = "Voltage Channel", color = "lightred") %>%
            visGroups(groupname = "none", color = "grey") %>%
            
            visLegend(width = 0.2, position = "right", main = "Group", ncol = 2) %>%
            
            visPhysics(
                forceAtlas2Based=list(
                    gravitationalConstant=-26,
                    centralGravity=0.005,
                    springLength=230,
                    springConstant=0.18,
                    avoidOverlap=1.5),
                maxVelocity=146,
                solver='forceAtlas2Based',
                timestep=0.35,
                enabled=FALSE,
                stabilization=list(
                    enabled=TRUE,
                    iterations=1000,
                    updateInterval=25))
        
        
    })
    
    output$DEGPlot <- renderPlotly({
        plot <- ggplot(pollen_DEGs, aes(x=log2FC_humanvsmac, y=log2FC_humanvschimp_org, text=Gene))+
            geom_point(alpha=0.5, colour="black")+
            geom_point(data = subset(pollen_DEGs, Gene==input$Gene), colour="red")+
            geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+
            scale_color_manual(values = c("red","black"))+
            scale_x_continuous(limits = c(-2,2))+scale_y_continuous(limits = c(-2,2))+
            xlab("Log2FC - Human vs. Macaque (tissue)")+
            ylab("Log2FC - Human vs. Chimp (organoid)")+
            ggtitle("Pollen et al., 2019")+
            theme_black+theme(plot.title = element_text(size = 15))
        
        ggplotly(plot, tooltip ="text")
        
    })
    
    output$benito_plot <- renderPlotly({
        ggplot(benito_gene_df(), aes(x=Day,y=log2(value), colour=Species, group=Species))+
            stat_summary(fun.data = mean_se, geom = "ribbon", fill="grey",alpha=0.5,colour="white", size=0.01)+
            stat_summary(geom = "line", aes(colour=Species, group=Species))+
            scale_colour_manual(values = c("black","red"))+
            scale_x_continuous(limits = c(input$slider_benito[1],input$slider_benito[2]))+
          coord_cartesian(xlim = NULL, ylim = NULL, expand = TRUE, default = FALSE, clip = "on")+
            theme_black+theme(legend.position = "right")+
            ggtitle("Benito-Kwiencinski et al., 2021")+xlab("Day")+ylab("Log2(CPM)")+ 
            theme(plot.title = element_text(size = 14))
    })
    
    #Make Rmarkdown report (sidebar)
    output$report <- downloadHandler(
      # For PDF output, change this to "report.pdf"
      filename = "report.html",
      content = function(file) {
        
        tempReport <- file.path(tempdir(), "report.Rmd")
        file.copy("report.Rmd", tempReport, overwrite = TRUE)
        tempReport2 <- file.path(tempdir(), "CellTypeLedgend_horizontal.png")
        file.copy("CellTypeLedgend_horizontal.png", tempReport2, overwrite = TRUE)
        
        # Set up parameters to pass to Rmd document
        params <- list(set_title = input$Gene,
                       plotInoue = plotInoue(),
                       plotCORTECON = plotCORTECON(),
                       plotLIBD = plotLIBD(),
                       plotMicali = plotMicali(),
                       plotBrainSpan = brainspan_plot(),
                       plotComparison = comparison_function(),
                       plotAllen = plotAllen(),
                       plotBenitoOrganoid = benito_function(),
                       plotCampOrganoid = camp_function(),
                       plotPollenOrganoid = pollen_function(),
                       plotSloanOrganoid = sloan_function(),
                       plotKantonOrganoid = kanton_function(),
                       plotCombinedCell = plotCombinedCell())
        
        rmarkdown::render(tempReport, output_file = file,
                          params = params,
                          envir = new.env(parent = globalenv())
        )
      })
      
      #Make Rmarkdown report (homepage)
      output$report_h <- downloadHandler(
        # For PDF output, change this to "report.pdf"
        filename = "report.html",
        content = function(file) {
          
          tempReport <- file.path(tempdir(), "report.Rmd")
          file.copy("report.Rmd", tempReport, overwrite = TRUE)
          tempReport2 <- file.path(tempdir(), "CellTypeLedgend_horizontal.png")
          file.copy("CellTypeLedgend_horizontal.png", tempReport2, overwrite = TRUE)

          # Set up parameters to pass to Rmd document
          params <- list(set_title = input$Gene,
                         plotInoue = plotInoue(),
                         plotCORTECON = plotCORTECON(),
                         plotLIBD = plotLIBD(),
                         plotMicali = plotMicali(),
                         plotBrainSpan = brainspan_plot(),
                         plotComparison = comparison_function(),
                         plotAllen = plotAllen(),
                         plotBenitoOrganoid = benito_function(),
                         plotCampOrganoid = camp_function(),
                         plotPollenOrganoid = pollen_function(),
                         plotSloanOrganoid = sloan_function(),
                         plotKantonOrganoid = kanton_function(),
                         plotCombinedCell = plotCombinedCell())
          
          rmarkdown::render(tempReport, output_file = file,
                            params = params,
                            envir = new.env(parent = globalenv())
          )
        })
}

shinyApp(ui, server)
