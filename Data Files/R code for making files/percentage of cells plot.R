#Load Libraries
setwd("~/Desktop/NeuroShiny")
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(readxl)
load("cell_type_plot_data.rda")
theme_black <- theme(panel.background=element_blank(),
                     panel.border = element_blank(),
                     axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
                     axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
                     plot.title = element_text(hjust = 0.5, size = 20),
                     legend.position = "none")

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
#14. Neuron Progenitor
c14 = "Neuron Progenitor"
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

#####Nowakowski######
nowakowski_cells1 <- read_excel("Data Files/aap8809_nowakowski_sm-tables-s1-s11.xlsx",sheet=3)
nowakowski_cells2 <- read_excel("Data Files/aap8809_nowakowski_sm-tables-s1-s11.xlsx",sheet=4)
load("Data Files/Nowakowski_cell_types.rda")

nowakowski_cells <- nowakowski_cells1 %>%
  inner_join(nowakowski_cells2, by=c("WGCNAcluster"="Cluster Name")) %>%
  inner_join(nowakowski_cell_subtype_df, by=c("Cluster Interpretation"="Cell Subtype")) %>%
  select(Cell_Type='Cell Type') %>%
  group_by(Cell_Type) %>%
  summarise(Cells = n()) %>%
  mutate(Study = "Nowakowski")

####Fan Data######
load("Data Files/Fan sample info.rda")

fan_cells <- fan_sample_info %>%
  select(Cell, Group)

fan_cells$Cell_Type <- "Other"
fan_cells$Cell_Type[fan_cells$Group == "Excitatory"] = "Excitatory Neuron"
fan_cells$Cell_Type[fan_cells$Group == "Inhibitory"] = "Inhibitory Neuron"
fan_cells$Cell_Type[fan_cells$Group == "Micro"] = c5
fan_cells$Cell_Type[fan_cells$Group == "Endo"] = c9
fan_cells$Cell_Type[fan_cells$Group == "Cajal"] = c12
fan_cells$Cell_Type[fan_cells$Group == "NSC"] = c1
fan_cells$Cell_Type[fan_cells$Group == "Cajal"] = c12
fan_cells$Cell_Type[fan_cells$Group == "Astro"] = c8
fan_cells$Cell_Type[fan_cells$Group== "Oligo"] = c7


fan_cells <-  fan_cells %>%
  group_by(Cell_Type) %>%
  summarise(Cells = n()) %>%
  mutate(Study = "Fan")

#####Polioudakis#####
polioudakis_data <- read_excel("Data Files/1-s2.0-S0896627319305616-mmc5.xlsx", sheet=1)

polioudakis_data$Cell_Type = "Other"
polioudakis_data$Cell_Type[polioudakis_data$Cluster %in% c("ExN","ExDp1","ExM-U","ExDp2","ExM")] = "Excitatory Neuron"
polioudakis_data$Cell_Type[polioudakis_data$Cluster %in% c("InMGE","InCGE")] = "Interneuron"
polioudakis_data$Cell_Type[polioudakis_data$Cluster %in% c("oRG","vRG")] = c2
polioudakis_data$Cell_Type[polioudakis_data$Cluster %in% c("PgG2M","PgS")] = "Neural Progenitor"
polioudakis_data$Cell_Type[polioudakis_data$Cluster %in% c("Mic")] = "Microglia"
polioudakis_data$Cell_Type[polioudakis_data$Cluster %in% c("OPC")] = "OPC"
polioudakis_data$Cell_Type[polioudakis_data$Cluster %in% c("End")] = "Endothelial"
polioudakis_data$Cell_Type[polioudakis_data$Cluster %in% c("IP")] = "Intermediate Progenitor Cell"

polioudakis_cells <-  polioudakis_data %>%
  select(Cell_Type, Cells = 'Number of cells') %>%
  group_by(Cell_Type) %>%
  summarise(Cells = sum(Cells)) %>%
  mutate(Study = "Polioudakis")

####Shi Data ####
shi_cells <- read_excel("Data Files/shi_science.abj6641_table_s2.xlsx", sheet = 1, skip=1) %>%
  select("Cell Subtype"= 'Major types') %>%
  mutate(Cell_Type="Other")
shi_cells$Cell_Type[shi_cells$'Cell Subtype' %in% c("MGE","CGE","LGE")] = "Inhibitory Neuron"
shi_cells$Cell_Type[shi_cells$'Cell Subtype' == "progenitor"] = "Neural Progenitor"
shi_cells$Cell_Type[shi_cells$'Cell Subtype' == "OPC"] = "OPC"
shi_cells$Cell_Type[shi_cells$'Cell Subtype' == "Excitatory neuron"] = "Excitatory Neuron"
shi_cells$Cell_Type[shi_cells$'Cell Subtype' == "Excitatory IPC"] = "Intermediate Progenitor Cell"
shi_cells$Cell_Type[shi_cells$'Cell Subtype' == "Microglia"] = "Microglia"
shi_cells$Cell_Type[shi_cells$'Cell Subtype' == "Endothelial"] = "Endothelial"
shi_cells$Cell_Type[shi_cells$'Cell Subtype' == "Thalamic neurons"] = "Thalamic neurons"

shi_cells  <- shi_cells %>%
  group_by(Cell_Type) %>%
  summarise(Cells = n())%>%
  mutate(Study = "Shi")



######Zhong Data####
zhong_cells <- read_excel("Data Files/zhong_GSE104276_readme_sample_barcode.xlsx", sheet=5) %>%
  select('Cell Subtype' = cell_types) %>%
  mutate(Cell_Type="Other")

zhong_cells$Cell_Type[zhong_cells$'Cell Subtype' =="Neurons"] = "Excitatory Neuron"
zhong_cells$Cell_Type[zhong_cells$'Cell Subtype' =="Stem cells"] = "Stem Cell"
zhong_cells$Cell_Type[zhong_cells$'Cell Subtype' =="GABAergic neurons"] = "Inhibitory Neuron"
zhong_cells$Cell_Type[zhong_cells$'Cell Subtype' =="Microglia"] = "Microglia"
zhong_cells$Cell_Type[zhong_cells$'Cell Subtype' =="OPC"] = "OPC"
zhong_cells$Cell_Type[zhong_cells$'Cell Subtype' =="Astrocytes"] = "Astrocyte"

zhong_cells  <- zhong_cells %>%
  group_by(Cell_Type) %>%
  summarise(Cells = n())%>%
  mutate(Study = "Zhong")



#### Plot Data ####
all_studies <-
  bind_rows(polioudakis_cells, fan_cells, nowakowski_cells, shi_cells,zhong_cells) %>%
  group_by(Study) %>%
  mutate(Cell_Type = factor(Cell_Type, levels = names(my_cell_cols)), Percentage = Cells/sum(Cells))

ggplot(all_studies, aes(x=Study, y=Percentage, fill=Cell_Type))+
  scale_fill_manual(values=my_cell_cols)+theme_black+
  coord_flip()+
  ggchicklet::geom_chicklet(width = 0.75)+
  scale_y_continuous(labels=scales::percent)+
  labs(fill="Cell Type")


ggplot(all_studies, aes(x=Study, y=Percentage, fill=Cell_Type))+
  scale_fill_manual(values=my_cell_cols)+theme_black+
  coord_flip()+
  geom_col(width = 0.75, colour="white")+
  scale_y_continuous(labels=scales::percent)+
  labs(fill="Cell Type")


summary_scRNAseq_studies <- all_studies

save(summary_scRNAseq_studies, file="Summary_scRNAseq.rda")

