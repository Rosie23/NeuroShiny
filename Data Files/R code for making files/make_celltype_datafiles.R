#Make Collective Cell Type Summary
setwd("~/Desktop/NeuroShiny")
library(RColorBrewer)
library(ggplot2)
library(readxl)
########Make Cell Type Colour Scheme#####
#NSC --> RGC --> IPC --> Neuron Progenitor --> Ex/Inhib Neuron
#             --> Astrocyte/OPC
#1. Stem Cell
c1 = "Stem Cell"
col1 = "lightgreen"
#2. Radial Glia
c2="Radial Glia"
col2="springgreen4"
#11. Intermediate Progenitor Cell
c11="Intermediate Progenitor Cell"
col11="lightslateblue"
#14.Progenitor
c14 = "Neural Progenitor"
col14 = "lightseagreen"
#3. Excitatory Neuron
c3="Excitatory Neuron"
col3="darkred"
#4. Inhibitory Neuron
c4="Inhibitory Neuron"
col4="steelblue4"
#6. Interneuron
c6="Interneuron"
col6="steelblue1"
#5. Microglial
c5="Microglia"
col5="salmon"
#7. OPC
c7="OPC"
col7="purple"
#8. Astrocyte
c8="Astrocyte"
col8="orange"
#9. Endothelial
c9="Endothelial"
col9="wheat"
#10. Choroid
c10="Choroid"
col10="yellow3"
#12. Cajal Retzuis
c12="Cajal Retzuis"
col12="yellowgreen"
#13. Other - Blood, Glycolysis, Unknown, Mural, Pericyte
c13="Other"
col13="grey50"
c15 = "Thalamic neurons"
col15 = "pink"

cell_type_color <- data.frame(Cell_Type = c(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14,c15),
                                 Colour =c(col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13, col14,col15))

cell_types <- cell_type_color$Cell_Type
my_cell_cols <- cell_type_color$Colour
names(my_cell_cols) <- cell_type_color$Cell_Type
#Set order of ledgend
my_cell_cols <- my_cell_cols[c("Stem Cell","Radial Glia","Intermediate Progenitor Cell","Neural Progenitor",
                               "Interneuron","Inhibitory Neuron","Excitatory Neuron","Thalamic neurons","Cajal Retzuis","Choroid",
                               "OPC","Astrocyte","Endothelial","Microglia","Other")]

#View colour pallette
ggplot(cell_type_color, aes(x=Cell_Type, y=Cell_Type, colour=Cell_Type))+geom_point(size=5)+theme_classic()+
  scale_colour_manual(values=my_cell_cols)+theme(legend.position = "right", axis.text.x = element_text(angle=90))

ggplot(cell_type_color, aes(x=Cell_Type, y=Cell_Type, colour=Cell_Type))+geom_point(size=5)+theme_classic()+
  scale_colour_manual(values=my_cell_cols)+theme(legend.position = "bottom", axis.text.x = element_text(angle=90))

#####Nowakowski######
load("Data Files/Nowakowski Cell markers.rda")
head(nowakowski_cell_markers)
unique(nowakowski_cell_markers$`Cluster Interpretation`)
nowakowski_celltypes <- unique(nowakowski_cell_markers$`Cluster Interpretation`)
nowakowski_cell_data<- nowakowski_cell_markers
nowakowski_cell_data$Cell_Type <- "Other"
nowakowski_cell_data$Cell_Type[nowakowski_cell_data$`Cluster Interpretation` %in% nowakowski_celltypes[c(3:7,22,23)]] = "Excitatory Neuron"
nowakowski_cell_data$Cell_Type[nowakowski_cell_data$`Cluster Interpretation` %in% nowakowski_celltypes[c(10:13,24)]] = "Inhibitory Neuron"
nowakowski_cell_data$Cell_Type[nowakowski_cell_data$`Cluster Interpretation` %in% nowakowski_celltypes[c(26:30,32,19)]] = "Radial Glia"
nowakowski_cell_data$Cell_Type[nowakowski_cell_data$`Cluster Interpretation` %in% nowakowski_celltypes[c(14:16)]] = "Intermediate Progenitor Cell"
nowakowski_cell_data$Cell_Type[nowakowski_cell_data$`Cluster Interpretation` %in% nowakowski_celltypes[c(17,18)]] = c14
nowakowski_cell_data$Cell_Type[nowakowski_cell_data$`Cluster Interpretation` == "Micrgolia"] = "Microglia"
nowakowski_cell_data$Cell_Type[nowakowski_cell_data$`Cluster Interpretation` == "Astocyte"] = "Astrocyte"
nowakowski_cell_data$Cell_Type[nowakowski_cell_data$`Cluster Interpretation` == "Endothelial"] = "Endothelial"
nowakowski_cell_data$Cell_Type[nowakowski_cell_data$`Cluster Interpretation` == "Choroid"] = "Choroid"
nowakowski_cell_data$Cell_Type[nowakowski_cell_data$`Cluster Interpretation` == "Oligodendrocyte progenitor cell"] = "OPC"

nowakowski_cell_subtype_df <- unique(nowakowski_cell_data[,c(1,7)])
colnames(nowakowski_cell_subtype_df) <- c("Cell Subtype","Cell Type")
save(nowakowski_cell_subtype_df, file="Data Files/Nowakowski_cell_types.rda")

#####Fan ####
load("~/Desktop/NeuroShiny Archive/NeuroShiny/Data by Gene/Fan/PAX6.RData")
head(fan_gene_df)
fan_cluster_group <- unique(fan_gene_df[,c(4,5)])
fan_cluster_group$Group <- as.character(fan_cluster_group$Group)

unique(fan_cluster_group$Group)
fan_cluster_group$Group[fan_cluster_group$Group == "Excitatory"] = "Excitatory Neuron"
fan_cluster_group$Group[fan_cluster_group$Group == "Inhibitory"] = "Inhibitory Neuron"
fan_cluster_group$Group[fan_cluster_group$Group == "Micro"] = c5
fan_cluster_group$Group[fan_cluster_group$Group == "Endo"] = c9
fan_cluster_group$Group[fan_cluster_group$Group == "Cajal"] = c12
fan_cluster_group$Group[fan_cluster_group$Group == "NSC"] = c1
fan_cluster_group$Group[fan_cluster_group$Group == "Cajal"] = c12
fan_cluster_group$Group[fan_cluster_group$Group == "Astro"] = c8
fan_cluster_group$Group[fan_cluster_group$Group == "Oligo"] = c7
fan_cluster_group$Group[fan_cluster_group$Group %in% c("Immune","Blood")] = "Other"
unique(fan_cluster_group$Group)
colnames(fan_cluster_group) <- c("Cell Subtype","Cell Type")

#####Polioudakis#####
polioudakis_data <- as.data.frame(readxl::read_xlsx("Data Files/1-s2.0-S0896627319305616-mmc5.xlsx", sheet=3))
polioudakis_data$FDR <- as.numeric(polioudakis_data$FDR)
polioudakis_data$`P-value` <- as.numeric(polioudakis_data$`P-value`)
polioudakis_data$Log2_fold_change <- as.numeric(polioudakis_data$Log2_fold_change)
polioudakis_data$FC <- 2^(polioudakis_data$Log2_fold_change)

pvals <- as.numeric(polioudakis_data$`P-value`)
pvals <- pvals[pvals!= 0]
min_pval <- min(pvals)
min_pval

unique(polioudakis_data$Cluster)
polioudakis_data$Cell_Type = "Other"
polioudakis_data$Cell_Type[polioudakis_data$Cluster %in% c("ExN","ExDp1","ExM-U","ExDp2","ExM")] = "Excitatory Neuron"
polioudakis_data$Cell_Type[polioudakis_data$Cluster %in% c("InMGE","InCGE")] = "Interneuron"
polioudakis_data$Cell_Type[polioudakis_data$Cluster %in% c("oRG","vRG")] = c2
polioudakis_data$Cell_Type[polioudakis_data$Cluster %in% c("PgG2M","PgS")] = "Neural Progenitor"
polioudakis_data$Cell_Type[polioudakis_data$Cluster %in% c("Mic")] = "Microglia"
polioudakis_data$Cell_Type[polioudakis_data$Cluster %in% c("OPC")] = "OPC"
polioudakis_data$Cell_Type[polioudakis_data$Cluster %in% c("End")] = "Endothelial"
polioudakis_data$Cell_Type[polioudakis_data$Cluster %in% c("IP")] = "Intermediate Progenitor Cell"

cell_types
polioudakis_cell_types <- unique(polioudakis_data[,c(3,10)])

colnames(polioudakis_cell_types) <- c("Cell Subtype","Cell Type")

######Zhong Data####
zhong_npc <- as.data.frame(readxl::read_xlsx("Data Files/41586_2018_BFnature25980_MOESM2_ESM.xlsx", sheet=2))
zhong_npc$Cluster <- paste0("NPC-",zhong_npc$cluster)
zhong_npc$Cell_Type <- "Stem Cell"
zhong_ex_n <- as.data.frame(readxl::read_xlsx("Data Files/41586_2018_BFnature25980_MOESM2_ESM.xlsx", sheet=4))
zhong_ex_n$Cluster <- paste0("ExNeun-",zhong_ex_n$cluster)
zhong_ex_n$Cell_Type <- c3
zhong_interneurons <- as.data.frame(readxl::read_xlsx("Data Files/41586_2018_BFnature25980_MOESM2_ESM.xlsx", sheet=6))
zhong_interneurons$Cluster <- paste0("Interneurons-",zhong_interneurons$cluster)
zhong_interneurons$Cell_Type <- c6
zhong_OPC <- as.data.frame(readxl::read_xlsx("Data Files/41586_2018_BFnature25980_MOESM2_ESM.xlsx", sheet=7))
zhong_OPC$Cluster <- paste0("OPC-",zhong_OPC$cluster)
zhong_OPC$Cell_Type <- c7
zhong_astro <- as.data.frame(readxl::read_xlsx("Data Files/41586_2018_BFnature25980_MOESM2_ESM.xlsx", sheet=8))
zhong_astro$Cluster <- paste0("Astrocyte-",zhong_astro$cluster)
zhong_astro$Cell_Type <- c8
zhong_microglia <- as.data.frame(readxl::read_xlsx("Data Files/41586_2018_BFnature25980_MOESM2_ESM.xlsx", sheet=9))
zhong_microglia$Cluster <- paste0("Microglia-",zhong_microglia$cluster)
zhong_microglia$Cell_Type <- c5

zhong_data <- dplyr::bind_rows(zhong_microglia[,-c(1,5)],zhong_astro[,-c(1,5)],zhong_interneurons[,-c(1,5)],
                               zhong_ex_n[,-c(1,5)],zhong_npc[,-c(1,5)], zhong_OPC[,-c(1,5)])

zhong_cell_types <- unique(zhong_data[,c(5,6)])
colnames(zhong_cell_types) <- c("Cell Subtype", "Cell Type")
#####Shi Data####
shi_df <- read_excel("Data Files/shi_science.abj6641_table_s3.xlsx", skip=1)

shi_cells <- shi_df %>%
  select("Cell Subtype"= cluster) %>%
  mutate("Cell Type"="Other") %>% unique()

shi_cells$'Cell Type'[shi_cells$'Cell Subtype' %in% c("MGE","CGE","LGE")] = "Inhibitory Neuron"
shi_cells$'Cell Type'[shi_cells$'Cell Subtype' == "progenitor"] = "Neural Progenitor"
shi_cells$'Cell Type'[shi_cells$'Cell Subtype' == "OPC"] = "OPC"
shi_cells$'Cell Type'[shi_cells$'Cell Subtype' == "Excitatory neuron"] = "Excitatory Neuron"
shi_cells$'Cell Type'[shi_cells$'Cell Subtype' == "Excitatory IPC"] = "Intermediate Progenitor Cell"
shi_cells$'Cell Type'[shi_cells$'Cell Subtype' == "Microglia"] = "Microglia"
shi_cells$'Cell Type'[shi_cells$'Cell Subtype' == "Endothelial"] = "Endothelial"
shi_cells$'Cell Type'[shi_cells$'Cell Subtype' == "Thalamic neurons"] = "Thalamic neurons"




######Merge Dataset Cell Types######
cell_subtype_df <- dplyr::bind_rows(nowakowski_cell_subtype_df,shi_cells,
                                    polioudakis_cell_types,zhong_cell_types)

#####Make Combined Cell Data to plot####
###Get data to plot - merge datasets 

psuedo_min_pval = 1e-16

#Fan Data
fan_cell_DEG2 <- as.data.frame(read_xlsx("Data Files/41422_2018_53_MOESM9_ESM.xlsx", sheet=2)) #12 Groups
fan_cell_DEG3 <- as.data.frame(read_xlsx("Data Files/41422_2018_53_MOESM9_ESM.xlsx", sheet=3)) #Glial
fan_cell_DEG4 <- as.data.frame(read_xlsx("Data Files/41422_2018_53_MOESM9_ESM.xlsx", sheet=4)) #non-neural
fan_cell_DEG <- rbind(fan_cell_DEG2,fan_cell_DEG3)
fan_cell_DEG <- rbind(fan_cell_DEG,fan_cell_DEG4)
unique(fan_cell_DEG$cluster)
#Convert Power to p-value?
#Differentially expressed genes of each cluster were identified by Seurat function ‘find_all_markers’ tested by ‘roc’ .
#The ROC test returns the 'power' for any individual genes (ranging from 0 - random, to 1 - perfect),
#and positive ‘avg_diff’ means upregulation in the corresponding clusters while negative ‘avg_diff’ means downregulation,
#'myAUC' is the area under the ROC curve (ranging from 0 to 1, 1-perfect),
#'#the ROC curve shows the true positive rate (TPR) against the false positive rate (FPR).

fan_cell_DEG$Pval <- 100^-(1/(1-fan_cell_DEG$power))+psuedo_min_pval #rough conversion of power value to p-value
cell_plot_data3 <- fan_cell_DEG[,c(7,8,2,6)]
colnames(cell_plot_data3) <- c("Gene","Pval","Fold.Change","Cluster")
cell_plot_data3$Dataset <- "Fan"

#Nowakowski
nowakowski_cell_markers$`Cluster Interpretation`[nowakowski_cell_markers$`Cluster Interpretation` =="Micrgolia"] = "Microglia"
cell_plot_data1 <- nowakowski_cell_markers[,c(6,2,3,1)]
colnames(cell_plot_data1) <- c("Gene","Pval","Fold.Change","Cluster")
cell_plot_data1$Dataset <- "Nowakowski"

#Polioudakis
cell_plot_data2 <- polioudakis_data[,c(2,5,4,3)]
colnames(cell_plot_data2) <- c("Gene","Pval","Fold.Change","Cluster")
#Add a psuedo minium pval so there is data to plot
cell_plot_data2$Pval <- as.numeric(cell_plot_data2$Pval)+min(cell_plot_data1$Pval)
cell_plot_data2$Dataset <- "Polioudakis"

#Zhong
cell_plot_data4 <- zhong_data
cell_plot_data4$Pval <- 100^-(1/(1-cell_plot_data4$power))+psuedo_min_pval
cell_plot_data4 <- cell_plot_data4[,c(4,7,2,5)]
colnames(cell_plot_data4) <- c("Gene","Pval","Fold.Change","Cluster")
cell_plot_data4$Dataset <- "Zhong"

c(min(cell_plot_data1$Pval), min(cell_plot_data2$Pval), min(cell_plot_data3$Pval), min(cell_plot_data4$Pval))

#shi
cell_plot_data5 <- shi_df %>%
  select(Gene = gene, Pval= p_val, Fold.Change = avg_logFC, Cluster=cluster) %>%
  mutate(Dataset = "Shi", Pval= Pval+psuedo_min_pval)

cell_plot_data <- dplyr::bind_rows(cell_plot_data1, cell_plot_data2, cell_plot_data3,cell_plot_data4, cell_plot_data5)



#Add Cell Type Group
cell_plot_data <- merge(cell_plot_data, cell_subtype_df, by.x="Cluster", by.y="Cell Subtype")
cell_plot_data <- unique(cell_plot_data)

#Save data
save(cell_plot_data, my_cell_cols, file="Data Files/cell_type_plot_data.rda")

my_gene = "DLX2"
gene_cell_plot_data <- subset(cell_plot_data, Gene==my_gene)
plot_limits=round(max(sqrt(gene_cell_plot_data$Fold.Change^2)),0)+1

plot=ggplot(gene_cell_plot_data, aes(x=Fold.Change,y=-(log10(Pval)), colour=`Cell Type`, shape=Dataset, label=Cluster))+
  geom_point(size=3)+
  scale_color_manual(values=my_cell_cols)+theme_bw()+
  theme(panel.border = element_blank(),panel.grid = element_blank(),
        axis.line = element_line(colour = "black", size=1),
        legend.title = element_blank(),
        legend.position = "none")+
  geom_vline(xintercept = 0, colour="grey")+scale_x_continuous(limits = c(-plot_limits,plot_limits))+
  ylab("-Log10(P val)")+xlab("Log2(Average Difference)")

  
plotly::ggplotly(plot,tooltip = c("Dataset","Cluster"))



######Add all cell subtypes for consistent colour scheme####
#####Kanton####
load("~/Desktop/NeuroShiny Archive/NeuroShiny_Complete/Data by Gene/Kanton/PAX6.RData")
kanton_cell_types <- data.frame(Kanton=unique(kanton_gene$CellType), Cell_Type="Other")
kanton_cell_types$Cell_Type[kanton_cell_types$Kanton %in% c("cortical neurons","deep layer neurons", "upper layer neurons")] = "Excitatory Neuron"
kanton_cell_types$Cell_Type[kanton_cell_types$Kanton %in% c("RGCs","RGCs early")] = "Radial Glia"
kanton_cell_types$Cell_Type[kanton_cell_types$Kanton %in% c("stem cells","NSCs","neuroepithelial-like cells")] = "Stem Cell"
kanton_cell_types$Cell_Type[kanton_cell_types$Kanton %in% c("IPs and early cortical neurons")] = "Intermediate Progenitor Cell"
kanton_cell_types$Cell_Type[kanton_cell_types$Kanton %in% c("choroid plexus/mesenchymal-like cells")] = "Choroid"
kanton_cell_types$Cell_Type[kanton_cell_types$Kanton %in% c("ventral progenitors and neurons")] = c4
kanton_cell_types$Cell_Type[kanton_cell_types$Kanton %in% c("dorsal progenitors","cycling ventral progenitors")] = c14
kanton_cell_types$Cell_Type[kanton_cell_types$Kanton %in% c("gliogenic/outer RGCs and astrocytes")] = "Astrocyte"

kanton_cell_types <- unique(kanton_cell_types)
colnames(kanton_cell_types) <- c("Cell Subtype", "Cell Type")

####Pollen####
load("~/Desktop/NeuroShiny Archive/NeuroShiny/Data by Gene/Pollen/PAX6.RData")
pollen_cell_types <- data.frame(`Cell Subtype`=unique(pollen_gene_df$ClusterName), `Cell Type`="Other")
colnames(pollen_cell_types) <- c("Cell Subtype", "Cell Type")
pollen_cell_types$`Cell Type`[pollen_cell_types$`Cell Subtype`%in% c("Newborn Excitatory Neuron","Maturing Excitatory Neuron",
"upper layer neurons","Excitatory neuron - deep Layer","Upper Layers","Hindbrain")] = "Excitatory Neuron"
pollen_cell_types$`Cell Type`[pollen_cell_types$`Cell Subtype`%in% c("Interneuron","Interneurons (MGE derived)",
                                                                     "Interneurons (CGE derived)")] = c6
pollen_cell_types$`Cell Type`[pollen_cell_types$`Cell Subtype`%in% c("Microglia")] = c5
pollen_cell_types$`Cell Type`[pollen_cell_types$`Cell Subtype`%in% c("Choroid/Ependymal","Choroid")] = "Choroid"
pollen_cell_types$`Cell Type`[pollen_cell_types$`Cell Subtype`%in% c("Astrocytes")] = c8
pollen_cell_types$`Cell Type`[pollen_cell_types$`Cell Subtype`%in% c("Endothelial")] = c9
pollen_cell_types$`Cell Type`[pollen_cell_types$`Cell Subtype`%in% c("Early RG G1/S","RG")] = c2
pollen_cell_types$`Cell Type`[pollen_cell_types$`Cell Subtype`%in% c("Oligodendrocyte precursor cells")] = c7
pollen_cell_types$`Cell Type`[pollen_cell_types$`Cell Subtype`%in% c("Intermediate Progenitor Cells")] = c11
pollen_cell_types$`Cell Type`[pollen_cell_types$`Cell Subtype`%in% c("RELN Excitatory Neuron")] = c12

cell_subtype_df <- dplyr::bind_rows(nowakowski_cell_subtype_df, fan_cluster_group,polioudakis_cell_types,
                                    zhong_cell_types,pollen_cell_types,kanton_cell_types)

#subcell col scheme
#cell_type_color <- data.frame(Cell_Type = c(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14),
#                              Colour =c(col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13, col14))

subcell_type_color <- merge(cell_subtype_df, cell_type_color, by.x="Cell Type", by.y="Cell_Type")
subcell_type_color$alpha <- 1
subcell_type_color$alpha[subcell_type_color$`Cell Type`=="Other"] = 0.5

my_subcell_alpha <- subcell_type_color$alpha
names(my_subcell_alpha) <- subcell_type_color$`Cell Subtype`


my_subcell_cols <- subcell_type_color$Colour
names(my_subcell_cols) <- subcell_type_color$`Cell Subtype`

  
cell_plot_data$Dataset <- factor(cell_plot_data$Dataset, levels = c("Fan","Nowakowski","Polioudakis","Zhong","Shi"))
dataset_shapes <- c(21,24,22,3,25)
names(dataset_shapes) <-levels(cell_plot_data$Dataset)

######Save all data#####
save(cell_plot_data, my_cell_cols,my_subcell_cols,my_subcell_alpha, dataset_shapes, file="cell_type_plot_data.rda")


