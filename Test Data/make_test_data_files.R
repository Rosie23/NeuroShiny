
clusters <- readxl::read_excel("~/Downloads/41556_2022_929_MOESM3_ESM.xlsx", sheet =7, skip = 1) %>%
  select(p_val = pval_adj, avg_log2FC=logfoldchange,cluster=cell_type, gene) %>%
  filter(avg_log2FC>0 & p_val<0.001) %>%
  write.table("~/Downloads/test_clusters_rosebrock.txt", quote=FALSE, row.names = FALSE, col.names = TRUE, sep="\t")
