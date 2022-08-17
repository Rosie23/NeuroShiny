####Making data files for NeuroShiny App
setwd("~/Desktop/NeuroShiny_Complete/Data Files")
library(readxl)
### Gene Info Table####
library(biomaRt)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
gene_info_data <- getBM(attributes=c('external_gene_name','ensembl_gene_id','entrezgene_id', 'chromosome_name','transcript_biotype'), mart = ensembl, useCache = FALSE)
#To View other attributes to add use listAttributes(ensembl) function
gene_info_colnames <- c("Ensembl ID","Entrez ID","Chromosome","BioType")
colnames(gene_info_data) <- c("Gene",gene_info_colnames)
head(gene_info_data)
dim(gene_info_data)
gene_info_data = unique(gene_info_data)
#Remove empty gene cells
gene_info_data <- gene_info_data[gene_info_data$Gene!="",]

#Collapse gene info table
gene_info_data_collapsed <- data.frame(Gene=as.character(), `Ensembl ID`=as.character(),
                                       `Entrez ID`=as.character(),Chromosome=as.character(),BioType=as.character())
for (i in unique(gene_info_data$Gene)) {
  tmp <- unique(subset(gene_info_data, Gene==i))
  df <- data.frame(Gene=i, `Ensembl ID`=paste(unique(tmp$`Ensembl ID`),collapse=", "),`Entrez ID`=paste(unique(tmp$`Entrez ID`),collapse=", "),
                   Chromosome=paste(unique(tmp$Chromosome),collapse=", "),BioType=paste(unique(tmp$BioType),collapse=", "))
  gene_info_data_collapsed <- rbind(gene_info_data_collapsed,df)
}

gene_info_data <- gene_info_data_collapsed
colnames(gene_info_data) <- c("Gene",gene_info_colnames)
save(gene_info_data, file="gene_info_data.rda")

####Make SFARI Data####
#Current release 09-02-2021
SFARI_data <- read.table("SFARI-Genes.csv", sep=",",header = T) 
SFARI_data <- SFARI_data[,-c(1,3,4,5)]
SFARI_colnames <- c("Genetic Category","Gene Score","Syndromic","Number of Reports")
colnames(SFARI_data) <- c("Gene",SFARI_colnames)
head(SFARI_data)
sfari_genes <- unique(SFARI_data$Gene)
save(SFARI_data, file="SFARI_data.rda")

###Make Bhaduri Data#####
primary_gene_cluster = read_xlsx("41586_2020_1962_MOESM3_ESM.xlsx", sheet=2)[,1:7]
primary_gene_cluster <- subset(primary_gene_cluster, `adjusted p_val`<0.05)
primary_cluster_info = read_xlsx("41586_2020_1962_MOESM3_ESM.xlsx", sheet=2)[,10:15]
primary_cluster_info <- primary_cluster_info[complete.cases(primary_cluster_info),]
tmp <- colnames(primary_cluster_info[2:6]) 
colnames(primary_cluster_info) <- c("Cluster", tmp)
primary_gene_cluster <- merge(primary_gene_cluster, primary_cluster_info[,1:5],
                              by.x="cluster", by.y="Cluster")
bhaduri_celltype_data <- primary_gene_cluster[,c(2,8:11)]
bhaduri_celltypes = unique(bhaduri_celltype_data$Type)
bhaduri_cellstates = unique(bhaduri_celltype_data$State)
bhaduri_cellclass = unique(bhaduri_celltype_data$Class)

#Collapse Data
genes = unique(bhaduri_celltype_data$Gene)
df_1 <- data.frame(Gene=as.character(), Cell.Type =as.character())
for(i in genes){
  cell_types <- subset(bhaduri_celltype_data, Gene==i)
  cell_types <- unique(cell_types$Type)
  tmp <- data.frame(Gene=i, Cell.Type = paste(cell_types, collapse = ", "))
  df_1 <- rbind(df_1, tmp)
}

df_2 <- data.frame(Gene=as.character(), State =as.character())
for(i in genes){
  cell_state <- subset(bhaduri_celltype_data, Gene==i)
  cell_state <- unique(cell_state$State)
  tmp <- data.frame(Gene=i, State = paste(cell_state, collapse = ", "))
  df_2 <- rbind(df_2, tmp)
}

df_3 <- data.frame(Gene=as.character(), Subtype =as.character())
for(i in genes){
  cell_Subtype <- subset(bhaduri_celltype_data, Gene==i)
  cell_Subtype <- unique(cell_Subtype$Subtype)
  tmp <- data.frame(Gene=i, Subtype = paste(cell_Subtype, collapse = ", "))
  df_3 <- rbind(df_3, tmp)
}
df_4 <- data.frame(Gene=as.character(), Class =as.character())
for(i in genes){
  cell_class <- subset(bhaduri_celltype_data, Gene==i)
  cell_class <- unique(cell_class$Class)
  tmp <- data.frame(Gene=i, Class = paste(cell_class, collapse = ", "))
  df_4 <- rbind(df_4, tmp)
}

df <- merge(df_1,df_2, by="Gene")
df <- merge(df,df_3, by="Gene")
df <- merge(df,df_4, by="Gene")
bhaduri_celltype_data_collapsed <- df[,c(1,5,3,2,4)]
bhaduri_colnames <- c("Class (Bhaduri et al.,)","State (Bhaduri et al.,)","Cell Type (Bhaduri et al.,)","Subtype (Bhaduri et al.,)")
colnames(bhaduri_celltype_data_collapsed) <- c("Gene",bhaduri_colnames)
colnames(bhaduri_celltype_data) <- c("Gene",bhaduri_colnames)
head(bhaduri_celltype_data_collapsed)
bhaduri_celltype_data_collapsed <- unique(bhaduri_celltype_data_collapsed)


save(bhaduri_celltype_data_collapsed, file="Bhaduri_CellType_data.rda")

###Fan Data####
fan_cell_DEG <- as.data.frame(read_xlsx("41422_2018_53_MOESM9_ESM.xlsx", sheet=1)) #8 Groups
fan_cell_DEG2 <- as.data.frame(read_xlsx("41422_2018_53_MOESM9_ESM.xlsx", sheet=2)) #12 Groups
clusters <- unique(fan_cell_DEG$cluster)
fan_celltypes = clusters

gene_marker_list = fan_cell_DEG$gene
fan_cell_markers <- data.frame(Gene=as.character(), Cell.Type =as.character())
for(i in 1:length(gene_marker_list)){
  marker_gene = gene_marker_list[i]
  cell_type = fan_cell_DEG[fan_cell_DEG$gene == marker_gene,][,6]
  tmp <- data.frame(Gene=marker_gene, Cell.Type = paste(cell_type, collapse = ", "))
  fan_cell_markers <- rbind(fan_cell_markers, tmp)
}

#Subcluster markers
fan_cell_DEG2 <- as.data.frame(read_xlsx("41422_2018_53_MOESM9_ESM.xlsx", sheet=2)) #12 Groups
fan_cell_DEG3 <- as.data.frame(read_xlsx("41422_2018_53_MOESM9_ESM.xlsx", sheet=3)) #Glial
fan_cell_DEG4 <- as.data.frame(read_xlsx("41422_2018_53_MOESM9_ESM.xlsx", sheet=4)) #non-neural

fan_subcell_DEG <- rbind(fan_cell_DEG2, fan_cell_DEG3)
fan_subcell_DEG <- rbind(fan_subcell_DEG, fan_cell_DEG4)

gene_marker_list = fan_subcell_DEG$gene
fan_subcell_markers <- data.frame(Gene=as.character(), Cell.Subtype =as.character())
for(i in 1:length(gene_marker_list)){
  marker_gene = gene_marker_list[i]
  cell_type = fan_subcell_DEG[fan_subcell_DEG$gene == marker_gene,][,6]
  tmp <- data.frame(Gene=marker_gene, Cell.Subtype = paste(cell_type, collapse = ", "))
  fan_subcell_markers <- rbind(fan_subcell_markers, tmp)
}

fan_cell_marker_df <- unique(merge(fan_subcell_markers, fan_cell_markers, by.x="Gene",by.y="Gene", all=T))
fan_cell_marker_df <- fan_cell_marker_df[,c(1,3,2)]
fan_colnames <- c("Cell Type (Fan et al.,)","Cell Subtype (Fan et al.,)")
colnames(fan_cell_marker_df) <- c("Gene",fan_colnames)
fan_cell_marker_df <- unique(fan_cell_marker_df)

#head(fan_cell_marker_df)
save(fan_cell_marker_df, file="Fan_CellType_data.rda")

###Zhong Data####
zhong <- as.data.frame(readxl::read_xlsx("41586_2018_BFnature25980_MOESM2_ESM.xlsx", sheet=1))
zhong <- unique(zhong[,c(1,5)])
unique(zhong$cluster)

zhong_npc <- as.data.frame(readxl::read_xlsx("41586_2018_BFnature25980_MOESM2_ESM.xlsx", sheet=2))
zhong_npc$Cell_Type <- "NPCs"
zhong_ex_n <- as.data.frame(readxl::read_xlsx("41586_2018_BFnature25980_MOESM2_ESM.xlsx", sheet=4))
zhong_ex_n$Cell_Type <- "Excitatory neurons"
zhong_interneurons <- as.data.frame(readxl::read_xlsx("41586_2018_BFnature25980_MOESM2_ESM.xlsx", sheet=6))
zhong_interneurons$Cell_Type <- "Interneurons"
zhong_OPC <- as.data.frame(readxl::read_xlsx("41586_2018_BFnature25980_MOESM2_ESM.xlsx", sheet=7))
zhong_OPC$Cell_Type <- "OPC"
zhong_astro <- as.data.frame(readxl::read_xlsx("41586_2018_BFnature25980_MOESM2_ESM.xlsx", sheet=8))
zhong_astro$Cell_Type <- "Astrocytes"
zhong_microglia <- as.data.frame(readxl::read_xlsx("41586_2018_BFnature25980_MOESM2_ESM.xlsx", sheet=9))
zhong_microglia$Cell_Type <- "Microglia"
zhong_all <- as.data.frame(readxl::read_xlsx("41586_2018_BFnature25980_MOESM2_ESM.xlsx", sheet=1))
zhong_all$Cell_Type <- zhong_all$cluster
zhong_all <- zhong_all[,c(1,3,7)]
colnames(zhong_all) <- colnames(zhong_microglia[,c(1,3,7)])

zhong <- dplyr::bind_rows(zhong_microglia[,c(1,3,7)],zhong_astro[,c(1,3,7)],zhong_interneurons[,c(1,3,7)],
                          zhong_ex_n[,c(1,3,7)],zhong_npc[,c(1,3,7)], zhong_OPC[,c(1,3,7)],zhong_all)

zhong <- unique(subset(zhong, avg_diff>0)[,-2])

unique(zhong$Cell_Type)

####Nowakowski Data####
load("Nowakowski Cell markers.rda")
nowakowski_cell_markers$`Cluster Interpretation`[nowakowski_cell_markers$`Cluster Interpretation` =="Micrgolia"] = "Microglia"
positive_nowakowski_cell_markers <- subset(nowakowski_cell_markers, avg_diff>0 & p_val<0.001)
nowakowski_cell_data <- unique(positive_nowakowski_cell_markers[,c(6,1)])
nowakowski_genes <- unique(nowakowski_cell_data$gene)
nowakowski_celltypes <- unique(nowakowski_cell_markers$`Cluster Interpretation`)

nowakowski_cell_data$Cell_Type <- nowakowski_cell_data$`Cluster Interpretation`
nowakowski_cell_data$Cell_Type[nowakowski_cell_data$`Cluster Interpretation` %in% nowakowski_celltypes[c(3:7,22,23)]] = "Excitatory Neuron"
nowakowski_cell_data$Cell_Type[nowakowski_cell_data$`Cluster Interpretation` %in% nowakowski_celltypes[c(10:13)]] = "Inhibitory Neuron"
nowakowski_cell_data$Cell_Type[nowakowski_cell_data$`Cluster Interpretation` %in% nowakowski_celltypes[c(26:30,32)]] = "Radial Glia"
nowakowski_cell_data$Cell_Type[nowakowski_cell_data$`Cluster Interpretation` %in% nowakowski_celltypes[c(14:16)]] = "Intermediate Progenitor Cell"
nowakowski_cell_data$Cell_Type[nowakowski_cell_data$`Cluster Interpretation` %in% nowakowski_celltypes[c(17:19,24)]] = "MGE progenitors/neurons"

nowakowski_celltypes <- unique(nowakowski_cell_data$Cell_Type)
nowakowski_celltypes_summary <- unique(nowakowski_cell_data[,c(2,3)])
unique(nowakowski_cell_data$Cell_Type)
nowakowski_celltypes_summary

#Generate cell type gene lists
unique(nowakowski_celltypes_summary$Cell_Type)
nowakowski_astocyte_genes = unique(subset(nowakowski_cell_data, Cell_Type== "Astocyte")[,1])
nowakowski_choroid_genes = unique(subset(nowakowski_cell_data, Cell_Type== "Choroid")[,1])
nowakowski_ex_neuron_genes = unique(subset(nowakowski_cell_data, Cell_Type== "Excitatory Neuron")[,1])
nowakowski_endo_genes = unique(subset(nowakowski_cell_data, Cell_Type== "Endothelial")[,1])
nowakowski_inhibitory_neuron_genes = unique(subset(nowakowski_cell_data, Cell_Type== "Inhibitory Neuron")[,1])
nowakowski_IPC_genes = unique(subset(nowakowski_cell_data, Cell_Type== "Intermediate Progenitor Cell")[,1])
nowakowski_microglia_genes = unique(subset(nowakowski_cell_data, Cell_Type== "Micrgolia")[,1])
nowakowski_mural_genes = unique(subset(nowakowski_cell_data, Cell_Type== "Mural/Pericyte")[,1])
nowakowski_OPC_genes = unique(subset(nowakowski_cell_data, Cell_Type== "Oligodendrocyte progenitor cell")[,1])
nowakowski_radial_glia_genes = unique(subset(nowakowski_cell_data, Cell_Type== "Radial Glia")[,1])
nowakowski_glycolysis_genes = unique(subset(nowakowski_cell_data, Cell_Type== "Glycolysis")[,1])

nowakowski_celltypes <- unique(nowakowski_cell_data$Cell_Type)
nowakowski_cellsubtypes <- unique(nowakowski_cell_data$`Cluster Interpretation`)

#Collapse data, multiple cell types in one field
nowakowski_cell_data_collpased_1 <- data.frame(Gene="",Cell.Type="")
for (i in nowakowski_genes) {
  tmp <- unique(subset(nowakowski_cell_data[,c(1,3)], gene==i))
  df <- data.frame(Gene=i, Cell.Type=paste(tmp$Cell_Type, collapse = ", "))
  nowakowski_cell_data_collpased_1 <- rbind(nowakowski_cell_data_collpased_1,df)
}

nowakowski_cell_data_collpased_2 <- data.frame(Gene="",Cell.Subtype="")
for (i in nowakowski_genes) {
  tmp <- unique(subset(nowakowski_cell_data[,1:2], gene==i))
  df <- data.frame(Gene=i, Cell.Subtype=paste(tmp$`Cluster Interpretation`, collapse = ", "))
  nowakowski_cell_data_collpased_2 <- rbind(nowakowski_cell_data_collpased_2,df)
}
nowakowski_cell_data_collpased <- merge(nowakowski_cell_data_collpased_1, nowakowski_cell_data_collpased_2, by.x="Gene",by.y="Gene")
nowakowski_cell_data = nowakowski_cell_data_collpased[-1,]
nowakowski_colnames <- c("Cell Type (Nowakowski et al.,)","Cell subtypes (Nowakowski et al.,)")
colnames(nowakowski_cell_data) <- c("Gene",nowakowski_colnames)
nowakowski_cell_data <- unique(nowakowski_cell_data)
head(nowakowski_cell_data)

save(nowakowski_cell_data, file="Nowakowski_CellType_data.rda")

#####Human Specific Regulation #####
# load("Pollen2019Data.rda")
# load("BenitoData.rda")
# 
# kanton_human_genes <- read_xlsx("~/Downloads/41586_2019_1654_MOESM3_ESM/Supplementary_Table_8.xlsx", skip=1)
# kanton_other_studies <- read_xlsx("~/Downloads/41586_2019_1654_MOESM3_ESM/Supplementary_Table_7.xlsx")[1:24,]
# benito_genes <- read_xlsx("~/Downloads/1-s2.0-S0092867421002397-mmc2.xlsx", sheet=7,skip=3)
# 
# kanton_human_genes <- kanton_human_genes[1:98,]
# 
# pollen_human_genes_up <- subset(pollen_DEGs, pollen_DEGs$log2FC_humanvsmac>0 &
#                                   pollen_DEGs$log2FC_humanvschimp_org>0)[,1]
# 
# pollen_human_genes_down <- subset(pollen_DEGs, pollen_DEGs$log2FC_humanvsmac<0 &
#                                     pollen_DEGs$log2FC_humanvschimp_org<0)[,1]
# 
# human_specific_genes1 <- data.frame(Gene=pollen_human_genes_up, Dataset="Pollen", Notes="Up in human org vs. chimp org & Up in human tissue vs. macaque tissue")
# human_specific_genes2 <- data.frame(Gene=pollen_human_genes_down, Dataset="Pollen", Notes="Down in human org vs. chimp org & Down in human tissue vs. macaque tissue")
# human_specific_genes3 <- data.frame(Gene=kanton_human_genes$Symbol, Dataset="Kanton", Notes="Human-specific DE genes compared to chimp")
# human_specific_genes4 <- data.frame(Gene=kanton_other_studies$symbol, Dataset=kanton_other_studies$source, Notes=kanton_other_studies$`expression pattern`)
# human_specific_genes5 <- data.frame(Gene=benito_genes$Gene, Dataset="Benito-Kwiecinski", Notes="Expression pattern different to Gorilla (orgs)")
# 
# human_specific_genes <- human_specific_genes[complete.cases(human_specific_genes),]
# 
# human_gene_list <-c(kanton_other_studies$symbol, kanton_human_genes$Symbol,
#                     pollen_human_genes_up,pollen_human_genes_down, benito_genes$Gene) 
# 
# 
# 
# human_specific_genes <- dplyr::bind_rows(human_specific_genes1,human_specific_genes2,human_specific_genes3,
#                                          human_specific_genes4,human_specific_genes5)
# human_genes_summary <- as.data.frame(table(human_gene_list))
# 
# save(human_specific_genes, human_genes_summary, file="Human_specific_genes.rda")
# 
load("Human_specific_genes.rda")

#Condense Data
df_1 <- data.frame(Gene=as.character(),Dataset=as.character())
for (i in unique(human_specific_genes$Gene)) {
  tmp_df <- subset(human_specific_genes, Gene==i)[,c(1,2)]
  tmp_1 <- data.frame(Gene=i, Dataset = paste(tmp_df$Dataset, collapse = ", "))
  df_1 <- rbind(df_1, tmp_1)
}

df_2 <- data.frame(Gene=as.character(),Notes=as.character())
for (i in unique(human_specific_genes$Gene)) {
  tmp_df <- subset(human_specific_genes, Gene==i)[,c(1,3)]
  tmp_2 <- data.frame(Gene=i, Notes = paste(tmp_df$Notes, collapse = ", "))
  df_2 <- rbind(df_2, tmp_2)
}

human_specific_genes <- merge(df_1, df_2, by.x = "Gene", by.y="Gene")
human_specific_genes <- human_specific_genes[complete.cases(human_specific_genes),]

library(stringr)
human_specific_genes$Freq <- str_count(human_specific_genes$Dataset, ',')+1

human_colnames <- c("Human Specific Regulation Dataset","Notes","Number of Datasets")
colnames(human_specific_genes) <- c("Gene",human_colnames)

save(human_specific_genes,human_colnames, file="Human_specific_genes_collapsed.rda")

#####BrainSpan Regional Expression#####
load("BrainSpan DEGs.rda")

FDR_cutoff = 0.05
log2FC_cutoff = 2

#Filter for protein coding genes
#library(biomaRt)
#ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
#gene_ids <- getBM(attributes=c('external_gene_name','ensembl_gene_id','description','gene_biotype'), mart = ensembl, useCache = FALSE)
#protein_coding_genes <- unique(subset(gene_ids, gene_biotype=="protein_coding")[,1])
#length(protein_coding_genes)

#Top genes for each region
###Top Anterior Cortex Genes####
region = "Anterior Cortex"
genes_1 <- subset(early_brainspan_DEGs, FDR.ACtx<FDR_cutoff & FC.ACtx>log2FC_cutoff)[1]
genes_1$Time <- "Early"
genes_2 <- subset(late_brainspan_DEGs, FDR.ACtx<FDR_cutoff & FC.ACtx>log2FC_cutoff)[1]
genes_2$Time <- "Mid"
genes_3 <- subset(mid_brainspan_DEGs, FDR.ACtx<FDR_cutoff & FC.ACtx>log2FC_cutoff)[1]
genes_3$Time <- "Late"
genes <- rbind(genes_1,genes_2)
genes <- rbind(genes,genes_3)
anterior_ctx_genes <- unique(genes$Gene)
genes_df <- data.frame(Gene=as.character(), Time =as.character())
for(i in genes$Gene){
  genes_tmp <- subset(genes, Gene == i)
  tmp <- data.frame(Gene=i, Time = paste(genes_tmp$Time, collapse = ", "))
  genes_df <- rbind(genes_df, tmp)
}
genes_df$Region <- region
genes_df$Info <- paste0(genes_df$Region," (",genes_df$Time,")")
regional_genes <- genes_df[,c(1,4)]
head(regional_genes)

###Top Parietal Cortex Genes####
region = "Parietal Cortex"
genes_1 <- subset(early_brainspan_DEGs, FDR.Parietal<FDR_cutoff & FC.Parietal>log2FC_cutoff)[1]
genes_1$Time <- "Early"
genes_2 <- subset(late_brainspan_DEGs, FDR.Parietal<FDR_cutoff & FC.Parietal>log2FC_cutoff)[1]
genes_2$Time <- "Mid"
genes_3 <- subset(mid_brainspan_DEGs, FDR.Parietal<FDR_cutoff & FC.Parietal>log2FC_cutoff)[1]
genes_3$Time <- "Late"
genes <- rbind(genes_1,genes_2)
genes <- rbind(genes,genes_3)
parietal_ctx_genes <- unique(genes$Gene)

genes_df <- data.frame(Gene=as.character(), Time =as.character())
for(i in genes$Gene){
  genes_tmp <- subset(genes, Gene == i)
  tmp <- data.frame(Gene=i, Time = paste(genes_tmp$Time, collapse = ", "))
  genes_df <- rbind(genes_df, tmp)
}
genes_df$Region <- region
genes_df$Info <- paste0(genes_df$Region," (",genes_df$Time,")")
regional_genes <- rbind(regional_genes, genes_df[,c(1,4)])
head(regional_genes)

#Top Temporal Cortex Genes
region = "Temporal Cortex"
genes_1 <- subset(early_brainspan_DEGs, FDR.Temporal<FDR_cutoff & FC.Temporal>log2FC_cutoff)[1]
genes_1$Time <- "Early"
genes_2 <- subset(late_brainspan_DEGs, FDR.Temporal<FDR_cutoff & FC.Temporal>log2FC_cutoff)[1]
genes_2$Time <- "Mid"
genes_3 <- subset(mid_brainspan_DEGs, FDR.Temporal<FDR_cutoff & FC.Temporal>log2FC_cutoff)[1]
genes_3$Time <- "Late"
genes <- rbind(genes_1,genes_2)
genes <- rbind(genes,genes_3)
genes_df <- data.frame(Gene=as.character(), Time =as.character())
for(i in genes$Gene){
  genes_tmp <- subset(genes, Gene == i)
  tmp <- data.frame(Gene=i, Time = paste(genes_tmp$Time, collapse = ", "))
  genes_df <- rbind(genes_df, tmp)
}
genes_df$Region <- region
genes_df$Info <- paste0(genes_df$Region," (",genes_df$Time,")")
regional_genes <- rbind(regional_genes, genes_df[,c(1,4)])
head(regional_genes)









###Top Temporal Cortex Genes####
region = "Temporal Cortex"
genes_1 <- subset(early_brainspan_DEGs, FDR.Temporal<FDR_cutoff & FC.Temporal>log2FC_cutoff)[1]
genes_1$Time <- "Early"
genes_2 <- subset(late_brainspan_DEGs, FDR.Temporal<FDR_cutoff & FC.Temporal>log2FC_cutoff)[1]
genes_2$Time <- "Mid"
genes_3 <- subset(mid_brainspan_DEGs, FDR.Temporal<FDR_cutoff & FC.Temporal>log2FC_cutoff)[1]
genes_3$Time <- "Late"
genes <- rbind(genes_1,genes_2)
genes <- rbind(genes,genes_3)
temporal_ctx_genes <- unique(genes$Gene)

genes_df <- data.frame(Gene=as.character(), Time =as.character())
for(i in genes$Gene){
  genes_tmp <- subset(genes, Gene == i)
  tmp <- data.frame(Gene=i, Time = paste(genes_tmp$Time, collapse = ", "))
  genes_df <- rbind(genes_df, tmp)
}
genes_df$Region <- region
genes_df$Info <- paste0(genes_df$Region," (",genes_df$Time,")")
regional_genes <- rbind(regional_genes, genes_df[,c(1,4)])
head(regional_genes)

#####Top Occipital####
region = "Occipital Cortex"
genes_1 <- subset(early_brainspan_DEGs, FDR.Occipital<FDR_cutoff & FC.Occipital>log2FC_cutoff)[1]
genes_1$Time <- "Early"
genes_2 <- subset(late_brainspan_DEGs, FDR.Occipital<FDR_cutoff & FC.Occipital>log2FC_cutoff)[1]
genes_2$Time <- "Mid"
genes_3 <- subset(mid_brainspan_DEGs, FDR.Occipital<FDR_cutoff & FC.Occipital>log2FC_cutoff)[1]
genes_3$Time <- "Late"
genes <- rbind(genes_1,genes_2)
genes <- rbind(genes,genes_3)
occipital_ctx_genes <- unique(genes$Gene)
genes_df <- data.frame(Gene=as.character(), Time =as.character())
for(i in genes$Gene){
  genes_tmp <- subset(genes, Gene == i)
  tmp <- data.frame(Gene=i, Time = paste(genes_tmp$Time, collapse = ", "))
  genes_df <- rbind(genes_df, tmp)
}
genes_df$Region <- region
genes_df$Info <- paste0(genes_df$Region," (",genes_df$Time,")")
regional_genes <- rbind(regional_genes, genes_df[,c(1,4)])
head(regional_genes)
#####Top Ventral####
region = "Ventral Forebrain"
genes_1 <- subset(early_brainspan_DEGs, FDR.Ventral<FDR_cutoff & FC.Ventral>log2FC_cutoff)[1]
genes_1$Time <- "Early"
genes_2 <- subset(late_brainspan_DEGs, FDR.Ventral<FDR_cutoff & FC.Ventral>log2FC_cutoff)[1]
genes_2$Time <- "Mid"
genes_3 <- subset(mid_brainspan_DEGs, FDR.Ventral<FDR_cutoff & FC.Ventral>log2FC_cutoff)[1]
genes_3$Time <- "Late"
genes <- rbind(genes_1,genes_2)
genes <- rbind(genes,genes_3)
ventral_genes <- unique(genes$Gene)
genes_df <- data.frame(Gene=as.character(), Time =as.character())
for(i in genes$Gene){
  genes_tmp <- subset(genes, Gene == i)
  tmp <- data.frame(Gene=i, Time = paste(genes_tmp$Time, collapse = ", "))
  genes_df <- rbind(genes_df, tmp)
}
genes_df$Region <- region
genes_df$Info <- paste0(genes_df$Region," (",genes_df$Time,")")
regional_genes <- rbind(regional_genes, genes_df[,c(1,4)])
head(regional_genes)


#####Top Hippocampus####
region = "Hippocampus"
genes_1 <- subset(early_brainspan_DEGs, FDR.Hippocampus<FDR_cutoff & FC.Hippocampus>log2FC_cutoff)[1]
genes_1$Time <- "Early"
genes_2 <- subset(late_brainspan_DEGs, FDR.Hippocampus<FDR_cutoff & FC.Hippocampus>log2FC_cutoff)[1]
genes_2$Time <- "Mid"
genes_3 <- subset(mid_brainspan_DEGs, FDR.Hippocampus<FDR_cutoff & FC.Hippocampus>log2FC_cutoff)[1]
genes_3$Time <- "Late"
genes <- rbind(genes_1,genes_2)
genes <- rbind(genes,genes_3)
hippocampus_genes <- unique(genes$Gene)
genes_df <- data.frame(Gene=as.character(), Time =as.character())
for(i in genes$Gene){
  genes_tmp <- subset(genes, Gene == i)
  tmp <- data.frame(Gene=i, Time = paste(genes_tmp$Time, collapse = ", "))
  genes_df <- rbind(genes_df, tmp)
}
genes_df$Region <- region
genes_df$Info <- paste0(genes_df$Region," (",genes_df$Time,")")
regional_genes <- rbind(regional_genes, genes_df[,c(1,4)])
head(regional_genes)


#####Top Amygdala####
region = "Amygdala"
genes_1 <- subset(early_brainspan_DEGs, FDR.Amygdala<FDR_cutoff & FC.Amygdala>log2FC_cutoff)[1]
genes_1$Time <- "Early"
genes_2 <- subset(late_brainspan_DEGs, FDR.Amygdala<FDR_cutoff & FC.Amygdala>log2FC_cutoff)[1]
genes_2$Time <- "Mid"
genes_3 <- subset(mid_brainspan_DEGs, FDR.Amygdala<FDR_cutoff & FC.Amygdala>log2FC_cutoff)[1]
genes_3$Time <- "Late"
genes <- rbind(genes_1,genes_2)
genes <- rbind(genes,genes_3)
amygdala_genes <- unique(genes$Gene)
genes_df <- data.frame(Gene=as.character(), Time =as.character())
for(i in genes$Gene){
  genes_tmp <- subset(genes, Gene == i)
  tmp <- data.frame(Gene=i, Time = paste(genes_tmp$Time, collapse = ", "))
  genes_df <- rbind(genes_df, tmp)
}
genes_df$Region <- region
genes_df$Info <- paste0(genes_df$Region," (",genes_df$Time,")")
regional_genes <- rbind(regional_genes, genes_df[,c(1,4)])
head(regional_genes)


#####Top Thalamus####
region = "Thalamus"
genes_1 <- subset(early_brainspan_DEGs, FDR.Thalamus<FDR_cutoff & FC.Thalamus>log2FC_cutoff)[1]
genes_1$Time <- "Early"
genes_2 <- subset(late_brainspan_DEGs, FDR.Thalamus<FDR_cutoff & FC.Thalamus>log2FC_cutoff)[1]
genes_2$Time <- "Mid"
genes_3 <- subset(mid_brainspan_DEGs, FDR.Thalamus<FDR_cutoff & FC.Thalamus>log2FC_cutoff)[1]
genes_3$Time <- "Late"
genes <- rbind(genes_1,genes_2)
genes <- rbind(genes,genes_3)
thalamus_genes <- unique(genes$Gene)
genes_df <- data.frame(Gene=as.character(), Time =as.character())
for(i in genes$Gene){
  genes_tmp <- subset(genes, Gene == i)
  tmp <- data.frame(Gene=i, Time = paste(genes_tmp$Time, collapse = ", "))
  genes_df <- rbind(genes_df, tmp)
}
genes_df$Region <- region
genes_df$Info <- paste0(genes_df$Region," (",genes_df$Time,")")
regional_genes <- rbind(regional_genes, genes_df[,c(1,4)])
head(regional_genes)


#####Top Cerebellum####
region = "Cerebellum"
genes_1 <- subset(early_brainspan_DEGs, FDR.Cerebellum<FDR_cutoff & FC.Cerebellum>log2FC_cutoff)[1]
genes_1$Time <- "Early"
genes_2 <- subset(late_brainspan_DEGs, FDR.Cerebellum<FDR_cutoff & FC.Cerebellum>log2FC_cutoff)[1]
genes_2$Time <- "Mid"
genes_3 <- subset(mid_brainspan_DEGs, FDR.Cerebellum<FDR_cutoff & FC.Cerebellum>log2FC_cutoff)[1]
genes_3$Time <- "Late"
genes <- rbind(genes_1,genes_2)
genes <- rbind(genes,genes_3)
cerebellum_genes <- unique(genes$Gene)
genes_df <- data.frame(Gene=as.character(), Time =as.character())
for(i in genes$Gene){
  genes_tmp <- subset(genes, Gene == i)
  tmp <- data.frame(Gene=i, Time = paste(genes_tmp$Time, collapse = ", "))
  genes_df <- rbind(genes_df, tmp)
}
genes_df$Region <- region
genes_df$Info <- paste0(genes_df$Region," (",genes_df$Time,")")
regional_genes <- rbind(regional_genes, genes_df[,c(1,4)])
head(regional_genes)

#####Tidy BrainSpan Regional Data####
regional_genes <- unique(regional_genes)
genes_df <- data.frame(Gene=as.character(), `Regional Expression` =as.character())
for(i in regional_genes$Gene){
  genes_tmp <- subset(regional_genes, Gene == i)
  tmp <- data.frame(Gene=i, `Regional Expression` = paste(genes_tmp$Info, collapse = ", "))
  genes_df <- rbind(genes_df, tmp)
}

brainspan_regional_data <- unique(genes_df)
colnames(brainspan_regional_data) <- c("Gene","Regional Expression (BrainSpan)")

brain.regions <- c("Anterior Cortex", "Parietal Cortex", "Temporal Cortex","Ventral","Hippocampus","Amygdala","Thalamus","Cerebellum")

save(brainspan_regional_data,brain.regions, file="BrainSpan_regional_data.rda")

###Make Gene Table####
#Show columns for certain datasets, user chooses dataset and then all the columns from that dataset are shown
gene_info_column_lists <- list(`Gene info`=gene_info_colnames, SFARI=SFARI_colnames, Fan=fan_colnames,
                               Bhaduri = bhaduri_colnames, Nowakowski=nowakowski_colnames, `Human Regulation`=human_colnames,
                               BrainSpan= "Regional Expression (BrainSpan)")

#save(gene_info_column_lists,file="Gene_info_column_lists.rda")

load("gene_info_data.rda")
load("SFARI_data.rda")
load("Bhaduri_CellType_data.rda")
load("Human_specific_genes_collapsed.rda")
load("Fan_CellType_data.rda")
load("Nowakowski_CellType_data.rda")
load("BrainSpan_regional_data.rda")
load("Gene_info_column_lists.rda")

#load("gene_info.rda")
gene_info <- merge(gene_info_data, SFARI_data, by.x = "Gene", by.y="Gene", all.x = T)
gene_info <- merge(gene_info, bhaduri_celltype_data_collapsed, by.x = "Gene", by.y="Gene", all.x = T)
gene_info <- merge(gene_info, fan_cell_marker_df, by.x = "Gene", by.y="Gene", all.x = T)
gene_info <- merge(gene_info, nowakowski_cell_data, by.x = "Gene", by.y="Gene", all.x = T)
gene_info <- merge(gene_info, human_specific_genes, by.x = "Gene", by.y="Gene", all.x = T)
gene_info <- merge(gene_info, brainspan_regional_data, by.x = "Gene", by.y="Gene", all.x = T)

#Order Gene Info Table by SFARI Score
gene_info <- gene_info[order(gene_info$`Number of Reports`,decreasing = TRUE),]

gene_info <- unique(gene_info)
head(gene_info)

######Make Gene Lists#####
all_genes <- unique(gene_info_data$Gene)
length(all_genes)
protein_coding_genes <- unique(subset(gene_info_data, BioType=="protein_coding")[,1])
length(protein_coding_genes)

#Polioudakis Gene Lists
polioudakis_data <- as.data.frame(readxl::read_xlsx("Data Files/1-s2.0-S0896627319305616-mmc5.xlsx", sheet=3))
polioudakis_ex_neurons <- unique(polioudakis_data[polioudakis_data$Cluster %in% c("ExN","ExDp1","ExM-U","ExDp2","ExM"),2])
polioudakis_interneuron <- unique(polioudakis_data[polioudakis_data$Cluster %in% c("InMGE","InCGE"),2])
polioudakis_RGC <- unique(polioudakis_data[polioudakis_data$Cluster %in% c("oRG","vRG"),2])
polioudakis_progenitors <- unique(polioudakis_data[polioudakis_data$Cluster %in% c("PgG2M","PgS"),2])
polioudakis_microglial <- unique(polioudakis_data[polioudakis_data$Cluster %in% c("Mic"),2])
polioudakis_OPC <- unique(polioudakis_data[polioudakis_data$Cluster %in% c("OPC"),2])
polioudakis_endo <- unique(polioudakis_data[polioudakis_data$Cluster %in% c("End"),2])
polioudakis_IPC <- unique(polioudakis_data[polioudakis_data$Cluster %in% c("IP"),2])

#Zhong Gene Lists
unique(zhong$Cell_Type)
zhong_NPCs <- unique(subset(zhong, Cell_Type== "NPCs")[,1])
zhong_ExN <- unique(subset(zhong, Cell_Type== "Excitatory neurons")[,1])
zhong_interneuron <- unique(subset(zhong, Cell_Type== "Interneurons")[,1])
zhong_OPC <- unique(subset(zhong, Cell_Type== "OPC")[,1])
zhong_astro <- unique(subset(zhong, Cell_Type== "Astrocytes")[,1])
zhong_microglia <- unique(subset(zhong, Cell_Type== "Microglia")[,1])

#Make Bhaduri Gene Lists
neuron_genes = unique(subset(bhaduri_celltype_data, `Class (Bhaduri et al.,)`== "Neuron")[,1])
non_neuron_genes = unique(subset(bhaduri_celltype_data, `Class (Bhaduri et al.,)`== "Non-neuronal")[,1])
cell_types <- unique(bhaduri_celltype_data$`Cell Type (Bhaduri et al.,)`)
bhaduri_radial_glia_genes = unique(subset(bhaduri_celltype_data, `Cell Type (Bhaduri et al.,)`== cell_types[1])[,1])
bhaduri_ex_neuron_genes = unique(subset(bhaduri_celltype_data, `Cell Type (Bhaduri et al.,)`== cell_types[2])[,1])
bhaduri_inhib_neuron_genes = unique(subset(bhaduri_celltype_data, `Cell Type (Bhaduri et al.,)`== cell_types[3])[,1])
bhaduri_IPC_genes = unique(subset(bhaduri_celltype_data, `Cell Type (Bhaduri et al.,)`== cell_types[4])[,1])
bhaduri_microglia_genes = unique(subset(bhaduri_celltype_data, `Cell Type (Bhaduri et al.,)`== cell_types[5])[,1])
bhaduri_RBC_genes = unique(subset(bhaduri_celltype_data, `Cell Type (Bhaduri et al.,)`== cell_types[6])[,1])
bhaduri_mural_genes = unique(subset(bhaduri_celltype_data, `Cell Type (Bhaduri et al.,)`== cell_types[7])[,1])
bhaduri_endo_genes = unique(subset(bhaduri_celltype_data, `Cell Type (Bhaduri et al.,)`== cell_types[8])[,1])
bhaduri_OPC_genes = unique(subset(bhaduri_celltype_data, `Cell Type (Bhaduri et al.,)`== cell_types[9])[,1])

#Fan Cell Type List
fan_all_DEG <- rbind(fan_subcell_DEG, fan_cell_DEG)

fan_inhib_neuron_genes = unique(subset(fan_all_DEG, cluster %in% c('Inhibitory_Neuron',"In_1","In_2",
                                                                  "In_3","In_4","In_5","In_6","In_7","In_8"))[,7])
fan_ex_neuron_genes = unique(subset(fan_all_DEG, cluster %in% c('Excitatory_Neuron',"Ex_1","Ex_2","Ex_4"))[,7])
fan_cajal_neuron_genes = unique(subset(fan_all_DEG, cluster== 'Cajal_Retzius')[,7])
fan_glia_genes = unique(subset(fan_all_DEG, cluster %in% c('Glia_Cell',"Astro_1","Astro_2"))[,7])
fan_microglia_genes = unique(subset(fan_all_DEG, cluster %in% c('Microglia',"Micro_1","Micro_2","Micro_3","Immune"))[,7])
fan_endo_genes = unique(subset(fan_all_DEG, cluster %in% c('Endothelial_Cell',"Endo_1","Endo_2"))[,7])
fan_blood_genes = unique(subset(fan_all_DEG, cluster %in% c('Blood',"B cell","Myeloiid cell",
                                                            "Naïve-like T cell","Effector T cell"))[,7])
fan_stem_cell_genes = unique(subset(fan_all_DEG, cluster %in% c('NSC_1',"NSC_2"))[,7])
fan_OPC_genes = unique(subset(fan_all_DEG, cluster %in% c('OPC',"Olig"))[,7])

#Remove NAs
human_genes = unique(human_specific_genes$Gene)
human_genes <- human_genes[!is.na(human_genes)]
human_genes <- as.character(human_genes)

gene_lists <- list(`All Genes`=all_genes,`Protein coding genes`=protein_coding_genes,
                   `SFARI (ASD genes)`= sfari_genes, `Human specific regulation`=human_genes,
                   `Neuronal (Bhaduri et al.,)`=neuron_genes,`Non-neuronal (Bhaduri et al.,)`=non_neuron_genes,
                   `Radial Glial genes (Bhaduri et al.,)`=bhaduri_radial_glia_genes,
                   `Excitatory Neuron genes (Bhaduri et al.,)`=bhaduri_ex_neuron_genes,
                   `Inhibitory Neuron genes (Bhaduri et al.,)`=bhaduri_inhib_neuron_genes,
                   `IPC genes (Bhaduri et al.,)`=bhaduri_IPC_genes,
                   `Microglial genes (Bhaduri et al.,)`=bhaduri_microglia_genes,
                   `RBC genes (Bhaduri et al.,)`=bhaduri_RBC_genes,
                   `OPC genes (Bhaduri et al.,)`=bhaduri_OPC_genes,
                   `Mural genes (Bhaduri et al.,)`=bhaduri_mural_genes,
                   `Endothelial genes (Bhaduri et al.,)`=bhaduri_endo_genes,
                   `Excitatory Neuron genes (Fan et al.,)`=fan_ex_neuron_genes,
                   `Inhibitor Neuron genes (Fan et al.,)`=fan_inhib_neuron_genes,
                   `Cajal Retzius genes (Fan et al.,)`=fan_cajal_neuron_genes,
                   `Astrocyte genes (Fan et al.,)`=fan_glia_genes,
                   `Microglial genes (Fan et al.,)`=fan_microglia_genes,
                   `Endothelial genes (Fan et al.,)`=fan_endo_genes,
                   `NPC genes (Fan et al.,)`=fan_stem_cell_genes,
                   `OPC genes (Fan et al.,)`=fan_OPC_genes,
                   `Blood genes (Fan et al.,)`=fan_blood_genes,
                   `Astrocyte genes (Nowakowski et al.,)`=nowakowski_astocyte_genes,
                   `Choroid genes (Nowakowski et al.,)`=nowakowski_choroid_genes,
                   `Excitatory Neuron genes (Nowakowski et al.,)`=nowakowski_ex_neuron_genes,
                   `Endothelial genes (Nowakowski et al.,)`=nowakowski_endo_genes,
                   `Inhibitory Neuron genes (Nowakowski et al.,)`=nowakowski_inhibitory_neuron_genes,
                   `IPC genes (Nowakowski et al.,)`=nowakowski_IPC_genes,
                   `Micrgolia genes (Nowakowski et al.,)`=nowakowski_microglia_genes,
                   `Mural/Pericyte genes (Nowakowski et al.,)`=nowakowski_mural_genes,
                   `OPC genes (Nowakowski et al.,)`=nowakowski_OPC_genes,
                   `Radial Glial genes (Nowakowski et al.,)`=nowakowski_radial_glia_genes,
                   `Glycolysis genes (Nowakowski et al.,)`=nowakowski_glycolysis_genes,
                   `NPCs genes (Zhong et al.,)`=zhong_NPCs,
                   `Excitatory Neuron genes (Zhong et al.,)`=zhong_ExN,
                   `Interneuron genes (Zhong et al.,)`=zhong_interneuron,
                   `OPC genes (Zhong et al.,)`=zhong_OPC,
                   `Astrocyte genes (Zhong et al.,)`=zhong_astro,
                   `Microglia genes (Zhong et al.,)`=zhong_microglia,
                   `Excitatory Neuron genes (Polioudakis et al.,)`=polioudakis_ex_neurons,
                   `Interneuron genes (Polioudakis et al.,)`=polioudakis_interneuron,
                   `Radial Glial genes (Polioudakis et al.,)`=polioudakis_RGC,
                   `Neural Progenitor genes (Polioudakis et al.,)`=polioudakis_progenitors,
                   `Microglia genes (Polioudakis et al.,)`=polioudakis_microglial,
                   `OPC genes (Polioudakis et al.,)`=polioudakis_OPC,
                   `Endothelial genes (Polioudakis et al.,)`=polioudakis_endo,
                   `IPC genes (Polioudakis et al.,)`=polioudakis_IPC,
                   `Anterior Ctx genes (BrainSpan)`=anterior_ctx_genes,
                   `Parietal Ctx genes (BrainSpan)`=parietal_ctx_genes,
                   `Occipital Ctx genes (BrainSpan)`=occipital_ctx_genes,
                   `Temporal Ctx genes (BrainSpan)`=temporal_ctx_genes,
                   `Ventral Forebrain genes (BrainSpan)`=ventral_genes,
                   `Amygdala genes (BrainSpan)`=amygdala_genes,
                   `Hippocampus genes (BrainSpan)`=hippocampus_genes,
                   `Thalamus genes (BrainSpan)`=thalamus_genes,
                   `Cerebellum genes (BrainSpan)`=cerebellum_genes)

names(gene_lists)

#####Save Data for app####
setwd("~/Desktop/NeuroShiny_Complete")
save(gene_info,gene_info_column_lists,gene_lists,
     nowakowski_celltypes, nowakowski_cellsubtypes,fan_celltypes,
     bhaduri_celltypes,bhaduri_cellstates, bhaduri_cellclass,
     file="Gene_Info_Table.rda")
