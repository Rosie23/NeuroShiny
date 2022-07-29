###My aesthetics
null_all_data <- function(){
  NPC_df = NULL
  comparison_df = NULL
  brainspan_gene_df = NULL
  HDBR_gene_df = NULL
  kanton_gene = NULL
  pollen_gene_df = NULL
  benito_gene_df = NULL
  region_DEGs_df = NULL
}

theme_black <- theme(panel.background=element_blank(),
                     panel.border = element_blank(),
                     axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
                     axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
                     plot.title = element_text(hjust = 0.5, size = 20),
                     legend.position = "none")

theme_axisblank = theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  axis.text.x = element_blank(), 
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  axis.line = element_blank()
)

#Make strings
module_list <- unique(pollen_modules$Module.ID)
brain.regions.DEGs <- c("ACtx"=2,"Parietal Cortex"=4,"Temporal Cortex"=6, "Occipital Cortex"=8,"Ventral Forebrain"=10, "Hippocampus"=12,"Amygdala"=14,"Thalamus"=16,"Cerebellum"=18)
brain.regions.DEGs.names <- c("ACtx","Parietal Cortex","Temporal Cortex", "Occipital Cortex","Ventral Forebrain", "Hippocampus","Amygdala","Thalamus","Cerebellum")
allen_regions <- c("rostral secondary prosencephalon","telencephalic vesicle","peduncular (caudal) hypothalamus", "prosomere 3","prosomere 2","prosomere 1","midbrain","prepontine hindbrain","pontine hindbrain", "pontomedullary hindbrain","medullary hindbrain (medulla)")

info_inoue <- ("Identification and Massively Parallel Characterization of Regulatory Elements Driving Neural Induction")

inoue_cols = c("#E41A1C", "#377EB8", "#4DAF4A")
micali_cols = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C") #brewer.pal(n=6, "Paired")
cortecon_cols = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")#brewer.pal(n=4, "Set1")
LIBD_cols = c(brewer.pal(n=9, "Reds")[c(3:7)], #line 165 n=5
              brewer.pal(n=9, "Greens")[c(4:6)], #line 21 n=3
              brewer.pal(n=9, "Blues")[c(4,6)], #line 3 n=2
              brewer.pal(n=9, "Oranges")[c(4,6)], #line 66 n=2
              brewer.pal(n=9, "Purples")[c(4,6)]) #line90 n=2

structure_cols <- c(brewer.pal("YlOrRd",n=9)[5:9], #anterior_ctx
                    brewer.pal("YlGnBu",n=9)[6:9], #parietal_ctx
                    brewer.pal("Greens",n=9)[6:9],#temporal_ctx
                    brewer.pal("Purples",n=9)[6:9],#vental_forebrain
                    brewer.pal("Oranges",n=9)[5:6],#occipital
                    brewer.pal("Set1",n=9)[9], #HIP
                    brewer.pal("Set1",n=9)[7],#AMY
                    brewer.pal("BrBG",n=11)[c(8,11)],#thalamus
                    brewer.pal("PiYG",n=11)[1:3]) #cerebellum

#HDBR Scheme
structure_cols_HDBR <- c(brewer.pal("YlOrRd",n=9)[5:8], #forebrain
                         brewer.pal("YlGnBu",n=9)[6:7], #cerebral cortex
                         brewer.pal("Purples",n=9)[6],#basal ganglion
                         brewer.pal("Set1",n=9)[4],#choroid plexus
                         brewer.pal("Set1",n=9)[5],#hippocampus
                         brewer.pal("BrBG",n=11)[c(8,11)], #diencephalon
                         brewer.pal("PiYG",n=11)[1:3],#midbrain
                         brewer.pal("PiYG",n=11)[c(4:7)],#hindbrain
                         brewer.pal("PiYG",n=11)[9:11]) #spinal

names(structure_cols_HDBR) <- c("forebrain","temporal lobe","telencephalon","forebrain fragment",
                                "cerebral cortex","brain fragment",
                                "basal ganglion",
                                "choroid plexus",
                                "hippocampus",
                                "diencephalon","pituitary and diencephalon",
                                "midbrain","forebrain and midbrain","diencephalon and midbrain",
                                "hindbrain","hindbrain fragment","hindbrain without cerebellum","cerebellum",
                                "spinal cord","pons","medulla oblongata")


###My Functions
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

extract_data <- function(gene, dataset_dir){
  files = list.files(paste0("Data by Gene/",dataset_dir,"/files"))
  file_name = paste0(gene,".RData")
  my_letter = substr(gene, 1, 1)
  if(!(file_name %in% files)){
    untar(paste0("Data by Gene/",dataset_dir,"/",my_letter,".tar.gz"),exdir=paste0("Data by Gene/",dataset_dir,"/files"))
  }
}


timecourse_plot <- function(data,title,col_scheme){
  ggplot(data, aes(x=Day, y=value, colour=Batch, group=Batch))+
    geom_point(aes(fill=Batch), colour="black", size=2, shape=21)+
    geom_line(size=1)+
    labs(y="Expression ( log2(RPKM+1))", x="Days post differentiation")+
    theme_black+ggtitle(title)+
    scale_fill_manual(values= col_scheme)+scale_color_manual(values= col_scheme)
}


.render_fishers_cluster_overlap_table <- function(output, data, filename)
{
  output$fishers_cluster_overlap_table <- DT::renderDataTable(
    server = FALSE,
    expr = {
      dt_options <- list(
        dom = "lftipB",
        searchHighlight = TRUE,
        scrollX = TRUE,
        search = list(regex = TRUE, caseInsensitive = TRUE, smart = FALSE),
        language = list(searchPlaceholder = "case insensitive search"),
        buttons = list(
          list(extend = "csv", filename = filename),
          list(extend = "excel", filename = filename)
        )
      )
      
      out <- DT::datatable(
        data = data,
        rownames = FALSE,
        filter = list(position = "top", clear = FALSE),
        selection = list(mode = "single"),
        options = dt_options,
        class = "compact stripe nowrap",
        extensions = "Buttons"
      ) 
      
      out <- DT::formatStyle(
        table = out,
        columns = colnames(data),
        fontSize = "90%"
      )
      
      return(out)
    }
  )
}

.render_summary_cluster_table <- function(output, data, filename)
{
  output$summary_clusters_table <- DT::renderDataTable(
    server = FALSE,
    expr = {
      dt_options <- list(
        dom = "lftipB",
        searchHighlight = TRUE,
        scrollX = TRUE,
        search = list(regex = TRUE, caseInsensitive = TRUE, smart = FALSE),
        language = list(searchPlaceholder = "case insensitive search"),
        buttons = list(
          list(extend = "csv", filename = filename),
          list(extend = "excel", filename = filename)
        )
      )
      
      out <- DT::datatable(
        data = data,
        rownames = FALSE,
        filter = list(position = "top", clear = FALSE),
        selection = list(mode = "single"),
        options = dt_options,
        class = "compact stripe nowrap",
        extensions = "Buttons"
      ) 
      
      out <- DT::formatStyle(
        table = out,
        columns = colnames(data),
        fontSize = "90%"
      )
      
      return(out)
    }
  )
}
                      