#Split Data up by Gene
setwd("~/Desktop/NeuroShiny_Complete/Data Files")
#####Step 1 Load Data and Packages#####
library(reshape2)
library(tidyverse)
load("Fan sample info.rda")
load("Pollen info.rda")
load("all_data.rda") #all data for 2D plots
load("brainspan_data.rda")
load("allen_data.rda")
load("FanData_v2.rda")
load("Pollen2019Data_updated.rda")
load("BenitoData.rda")
load("Pasca_Data.rda")
load("Luo_Data.rda")
load("Kanton_data.rda")
load("Gene lists.rda")

#####Step 2 Get a list of genes you want to make files for######
library(biomaRt)
ensembl=useMart("ensembl")
ensemblHuman = useDataset("hsapiens_gene_ensembl",mart=ensembl)
#Set filters biotype=protein_coding to select just protein-coding genes
#gene.list = getBM(attributes=c("external_gene_name"), filters='biotype', values=c('protein_coding'), mart=ensemblHuman)
gene.list = getBM(attributes=c("external_gene_name"), mart=ensemblHuman)
gene.list <- unique(gene.list$external_gene_name)
length(gene.list)

#####Split gene list by letter
library(dplyr)

gene_az_list <- list()
for (letter in LETTERS) {
tmp <- gene.list[(substr(gene.list, 1, 1) == letter)]
tmp = list(tmp)
names(tmp) = letter
gene_az_list <- c(gene_az_list, tmp)
}
head(gene_az_list)

#Test gene list 
#gene.list <- c("PAX6","FOXG1","SOX2","OTX2","SOX1","CDH2")

#Neural Gene List
#load("Gene_Info_Table.rda")
#sfari_genes <- unique(subset(gene_info, gene_info$SAFARI.Gene.Score>0)[,1])
#length(sfari_genes)
#Neural genes
#go_genes <- subset(all_go_genes, go_id == "GO:0007420")[,1] #brain development go term
#length(go_genes)
#gene.list <- unique(c(sfari_genes, go_genes, gene.list))

#####Step 3 Save files within one Supplementary Data file#####
setwd("~/Desktop/NeuroShiny_Complete")
save(file="Supplementary Data.rda",fan_sample_info,gene.list,
     pollen_modules,pollen_module_info,kanton_cell_order,pollen_celltype_order,
     cluster_col_fan,allen_data,pollen_DEGs)
#####Step 4 Create Folders for data files#########
dir.create("Data by Gene AZ")
setwd("~/Desktop/NeuroShiny_Complete/Data by Gene AZ")
dir.create("Kanton")
dir.create("Fan")
dir.create("Pollen")
dir.create("BrainSpan")
dir.create("NPC")
dir.create("Comparison")
dir.create("Benito")

#####Step 5 Make files#####
brainspan_counts$Gene <- as.character(brainspan_counts$Gene)

#Function for making a file for all datasets for one gene
make_gene_file_all <- function(my_gene){
  #in vitro 2D plot
  if(my_gene %in% rownames(counts_matrix)){
    NPC_df  <- counts_matrix[my_gene,]
    NPC_df <- reshape2::melt(NPC_df)
    colnames(NPC_df) <- c("Sample","value")
    NPC_df <- inner_join(NPC_df, sample_info, by="Sample")
  } else NPC_df=NULL
  
  #BrainSpan
  if(my_gene %in% brainspan_counts$Gene){
    brainspan_gene_df <- filter(brainspan_counts, Gene == my_gene)[,-1]
    brainspan_gene_df <- reshape2::melt(brainspan_gene_df)
    colnames(brainspan_gene_df) <- c("column_num","value")
    brainspan_gene_df$column_num <- as.integer(brainspan_gene_df$column_num)
    brainspan_gene_df <- inner_join(brainspan_gene_df, brainspan_sample, by="column_num")
  }else brainspan_gene_df=NULL
  
  #Kanton
  if(my_gene %in% rownames(kanton_chimp_average$RNA)){
    kanton_chimp <- kanton_chimp_average$RNA[my_gene,]
    kanton_chimp <- data.frame(Sample_Cluster = names(kanton_chimp), value=kanton_chimp,Species="Chimp",stringsAsFactors = F)
    kanton_chimp <- inner_join(kanton_chimp, kanton_chimp_info, by="Sample_Cluster")
    kanton_chimp <- kanton_chimp[,c("value","Species","CellType","Day")]
  } else kanton_chimp <- data.frame(value=as.numeric(), Species=as.character(), CellType=as.character(),Day=as.numeric())
  
  if(my_gene %in% rownames(kanton_human_average$RNA)){
    kanton_human <- kanton_human_average$RNA[my_gene,]
    kanton_human <- data.frame(Sample_Cluster = names(kanton_human), value=kanton_human, Species="Human",stringsAsFactors = F)
    kanton_human <- inner_join(kanton_human, kanton_human_info, by="Sample_Cluster")
    kanton_human <- kanton_human[,c("value","Species","CellType","Day")]
  } else kanton_human <- data.frame(value=as.numeric(), Species=as.character(), CellType=as.character(),Day=as.numeric())
  
  kanton_gene <- rbind(kanton_chimp, kanton_human)
  kanton_gene$CellType <- factor(kanton_gene$CellType, levels = kanton_cell_order)
  
  #Pollen
  if(my_gene %in% rownames(pollen_df)){
    pollen_gene_df <- pollen_df[rownames(pollen_df) == my_gene,]
    pollen_gene_df <- data.frame(Cell.ID =names(pollen_gene_df), Expression=pollen_gene_df)
    pollen_gene_df <- inner_join(pollen_gene_df, pollen_sample_info, by="Cell.ID")
    pollen_gene_df$ClusterName <- factor(pollen_gene_df$ClusterName, levels = pollen_celltype_order)
    pollen_gene_df <- pollen_gene_df[,c("Expression","CellSource","Organism","Age","Laboratory","Dataset","ClusterName")]
    pollen_gene_df$Organism <- factor(pollen_gene_df$Organism, levels = c("Human","Chimp","Macaque"))
    pollen_gene_df$Laboratory <- factor(pollen_gene_df$Laboratory, levels = c("Kriegstein","Pollen","Pasca","Treutlein"))
    pollen_gene_df$Dataset <- factor(pollen_gene_df$Dataset, levels = c("1.Hp","2.Mp","3.Ho","4.Co","5.Ho","6.Co", "7.Hs"))
    pollen_gene_df$ClusterName <- factor(pollen_gene_df$ClusterName, levels = pollen_celltype_order)
    
  } else pollen_gene_df = NULL
  
  #Fan
  if(my_gene %in% rownames(fan_counts)){
    fan_gene_df <- fan_counts[rownames(fan_counts) == my_gene,]
    fan_gene_df <- data.frame(Expression=fan_gene_df, Cell = names(fan_gene_df))
    #fan_gene_df$Expression[fan_gene_df$Expression==0] <- NA
    #fan_gene_df <- fan_gene_df[complete.cases(fan_gene_df),]
    fan_gene_df <- inner_join(fan_gene_df, fan_sample_info, by="Cell")
    fan_gene_df <- fan_gene_df[,c("Expression","tSNE1","tSNE2", "Cluster","Group","subRegion","Region","Sample")]
  }else fan_gene_df = NULL
  
  #Comparison
  if(my_gene %in% rownames(counts_matrix)){
    NPC_df2 <- NPC_df[,c("Experiment","value","Day")]
    NPC_df2$Model <- "2D in-vitro"
    NPC_df2$Species <- "Human"
  }else {NPC_df2 <- data.frame(Experiment=as.character(),value=as.numeric(),
                               Day=as.numeric(),Model=as.character(),Species=as.character())}
  
  if(my_gene %in% brainspan_counts$Gene){
    brainspan <- subset(brainspan_gene_df, brainspan_gene_df$Region=="Anterior Cortex")
    brainspan$Experiment <- "BrainSpan"
    brainspan$Day <- (brainspan$age*7)
    brainspan <- brainspan[,c("Experiment","value","Day")]
    brainspan$Model <- "Tissue"
    brainspan$Species <- "Human"
  } else {brainspan <- data.frame(Experiment=as.character(),value=as.numeric(),
                                  Day=as.numeric(),Model=as.character(),Species=as.character())}
  
  if(my_gene %in% allen_data$Symbol){
    allen_gene_df <- filter(allen_data, Symbol==my_gene)
    allen_gene_df$Day <- 0.0004*(allen_gene_df$days)^3+0.1787*(allen_gene_df$days)^2+0.7331*(allen_gene_df$days)+0.4441
    allen_gene_df <- subset(allen_gene_df, allen_gene_df$name == "telencephalic vesicle")
    allen_gene_df <- allen_gene_df[,c("expression_energy","Day")]
    colnames(allen_gene_df) <- c("value","Day")
    allen_gene_df$Experiment <- "Allen"
    allen_gene_df <- allen_gene_df[,c("Experiment","value","Day")]
    allen_gene_df$Model <- "Tissue"
    allen_gene_df$Species <- "Mouse"
  } else {allen_gene_df <- data.frame(Experiment=as.character(),value=as.numeric(),
                                      Day=as.numeric(),Model=as.character(),Species=as.character())}
  
  if(my_gene %in% row.names(benito_counts)){
    benito_gene_df <- benito_counts[my_gene,]
    benito_gene_df <- reshape2::melt(benito_gene_df)
    colnames(benito_gene_df) <- c("Sample","value")
    benito_gene_df <- dplyr::inner_join(benito_gene_df, benito_info, by="Sample")
    benito_gene_df$Model <- "Organoid"
    benito_gene_df <- benito_gene_df[,c("Experiment", "value", "Day","Model", "Species")]
  } else {benito_gene_df <- data.frame(Experiment=as.character(),value=as.numeric(),
                                       Day=as.numeric(),Model=as.character(),Species=as.character())}
  
  if(my_gene %in% luo_counts$Gene){
    luo_gene_df  <- subset(luo_counts, Gene==my_gene)
    luo_gene_df <- reshape2::melt(luo_gene_df)
    luo_gene_df <- merge(luo_gene_df, luo_sample_info, by.x="variable", by.y="Sample")
    luo_gene_df$Model <- "Organoid"
    luo_gene_df$Experiment <- "Luo"
    luo_gene_df$Species <- "Human"
    luo_gene_df <- luo_gene_df[,c("Experiment", "value", "Day","Model", "Species")]
  } else {luo_gene_df <- data.frame(Experiment=as.character(),value=as.numeric(),
                                    Day=as.numeric(),Model=as.character(),Species=as.character())}
  
  if(my_gene %in% row.names(pasca_data)){
    pasca_gene_df <- pasca_data[my_gene,]
    pasca_gene_df <- data.frame(Cell=names(pasca_gene_df),value=pasca_gene_df)
    row.names(pasca_gene_df) = NULL
    pasca_gene_df <- dplyr::inner_join(pasca_gene_df, pasca_info, by="Cell")
    pasca_gene_df$Model <- "Organoid"
    pasca_gene_df$Experiment <- "Sloan"
    pasca_gene_df$Species <- "Human"
    pasca_gene_df <- pasca_gene_df[,c("Experiment", "value", "Day","Model", "Species")]
  } else {pasca_gene_df <- data.frame(Experiment=as.character(),value=as.numeric(),
                                      Day=as.numeric(),Model=as.character(),Species=as.character())}
  
  if(my_gene %in% row.names(kanton_chimp_average_byStage$RNA)){
    kanton_chimp_mean <- kanton_chimp_average_byStage$RNA[my_gene,]
    kanton_chimp_mean <- reshape2::melt(kanton_chimp_mean)
    kanton_chimp_mean$Stage <- rownames(kanton_chimp_mean)
    kanton_chimp_mean <- inner_join(kanton_chimp_mean, kanton_chimp_info, by="Stage")
    kanton_chimp_mean$Model <- "Organoid"
    kanton_chimp_mean$Experiment <- "Kanton-Chimp"
    kanton_chimp_mean$Species <- "Chimp"
    kanton_chimp_mean <- kanton_chimp_mean[,c("Experiment", "value", "Day","Model", "Species")]
  } else {kanton_chimp_mean <- data.frame(Experiment=as.character(),value=as.numeric(),
                                          Day=as.numeric(),Model=as.character(),Species=as.character())}
  
  if(my_gene %in% row.names(kanton_human_average_byStage)){
    kanton_human_mean <- kanton_human_average_byStage$RNA[my_gene,]
    kanton_human_mean <- reshape2::melt(kanton_human_mean)
    kanton_human_mean$Stage <- rownames(kanton_human_mean)
    kanton_human_mean <- dplyr::inner_join(kanton_human_mean, kanton_human_info, by="Stage")
    kanton_human_mean$Model <- "Organoid"
    kanton_human_mean$Experiment <- "Kanton-Human"
    kanton_human_mean$Species <- "Human"
    kanton_human_mean <- kanton_human_mean[,c("Experiment", "value", "Day","Model", "Species")]
  } else {kanton_human_mean <- data.frame(Experiment=as.character(),value=as.numeric(),
                                          Day=as.numeric(),Model=as.character(),Species=as.character())}
  
  comparison_df <- bind_rows(NPC_df2, brainspan, allen_gene_df, benito_gene_df,kanton_chimp_mean,
                             luo_gene_df,pasca_gene_df,kanton_human_mean)
  
  save(file=paste0("NPC/",my_gene,".RData"),NPC_df)
  
  save(file=paste0("BrainSpan/",my_gene,".RData"),brainspan_gene_df)
  
  save(file=paste0("Kanton/",my_gene,".RData"),kanton_gene)
  
  save(file=paste0("Pollen/",my_gene,".RData"),pollen_gene_df)
  
  save(file=paste0("Comparison/",my_gene,".RData"),comparison_df)
  
  save(file=paste0("Fan/",my_gene,".RData"),fan_gene_df)
  
  save(file=paste0("Benito/",my_gene,".RData"),benito_gene_df)
}


#Loop through your list of genes to make a file for each gene
for(i in 1:length(gene.list)){
  my_gene = gene.list[i]
  make_gene_file_all(my_gene)
  percent_complete = (i/length(gene.list))*100
  print(paste(round(percent_complete,2),"% complete"))
}

#Compress genes of same letter into one zip
#setwd("~/Desktop/NeuroShiny_Complete")
#for (letter in LETTERS) {
#  genes <- gene_az_list[letter][[1]]
#  for(gene in genes){
#    
#    gene <- get(load(paste0("Data by Gene/Kanton/",gene,".RData")))
#    gene <- get(load(paste0("Data by Gene/Kanton/",gene,".RData")))
#    save(gene, file = paste0(letter,".rda")
#  }
#
#}

#Open compress files
#In terminal compress each file into a folder
# cd ~/Desktop/NeuroShiny_Complete/Data_by_Gene/Pollen
#for i in {A..Z}; do echo "$i"; tar -zcf ${i}.tar.gz ${i}*; done

#for i in "${A..Z}"; do echo i; done


#tar -zcf A.tar.gz A*
my_letter="A"
my_gene = "PAX6"
my_letter = substr(my_gene, 1, 1)
setwd("~/Desktop/NeuroShiny_AZ")

files = list.files("Data by Gene/Benito/files")
file_name = paste0(my_gene,".RData")
if(!(file_name %in% files)){
  untar(paste0("Data by Gene/Benito/",my_letter,".tar.gz"),exdir="Data by Gene/Benito/files")
}


my_gene = "A1BG"
setwd("~/Desktop/NeuroShiny_AZ/Data by Gene/Benito")
my_letter = substr(my_gene, 1, 1)
untar(paste0(my_letter,".tar.gz")) # exdir="files"
kanton_gene <- get(load(paste0("files/",my_gene,".RData")))
head(kanton_gene)

#####Optimizing - to see what bit of code is running slow#####
library(profvis)
system.time(make_gene_file("GBX2"))
profvis(make_gene_file("GBX2"))



