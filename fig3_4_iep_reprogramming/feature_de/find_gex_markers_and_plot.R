#RNA markers

library(Seurat)
library(stringr)

b4rna <- readRDS("./inter_files_rev1/all_rna_new_rev1.Rds")
b4rna_sub <- b4rna[,b4rna$md_fate_rev1!='na']

Idents(b4rna_sub) <- b4rna_sub$md_fate_rev1
markers_mdfate <- FindAllMarkers(object = b4rna_sub, max.cells.per.ident=1000)
saveRDS(markers_mdfate, "inter_files_rev1/mdfate_markers_new_rev1.Rds")

Idents(b4rna_sub) <- b4rna_sub$md_fate_coarse_rev1
markers_mdfate_c <- FindAllMarkers(object = b4rna_sub, max.cells.per.ident=1000)
saveRDS(markers_mdfate_c, "inter_files_rev1/mdfate_coarse_markers_new_rev1.Rds")

#da3 vs mefs
b4rna_sub <- b4rna_sub[,b4rna_sub$sample %in% c("D0","B4D3")]
Idents(b4rna_sub) <- b4rna_sub$md_fate_rev1
markers_mdfate_d3 <- FindAllMarkers(object = b4rna_sub, max.cells.per.ident=1000)
saveRDS(markers_mdfate_d3, "inter_files_rev1/mdfate_d3_markers_new_rev1.Rds")

#plot uniquely enriched markers on a heatmap
md_fate_new <- readRDS("./inter_files_rev1//mdfate_d3_markers_new_rev1.Rds")

md_fate_new <- md_fate_new[md_fate_new$p_val_adj < 0.05,]
md_fate_new <- md_fate_new[md_fate_new$avg_log2FC > 0,]
length(unique(md_fate_new$gene))
write.csv(md_fate_new,"./final_files_rev1//supp_table_d3_de_genes_rev1.csv")

b4rna_sub <- b4rna[,b4rna$md_fate_rev1!='na']
b4rna_sub2 <- b4rna_sub[,b4rna_sub$sample %in% c("D0","B4D3")]
Idents(b4rna_sub2) <- b4rna_sub2$md_fate_rev1

md_fate_new <- md_fate_new[!(md_fate_new$gene %in% c("CellTag.UTR",'GFP.CDS','FoxA1.HNF4a')),]

gene_freq <- table(md_fate_new$gene)
d3_unique <- md_fate_new[md_fate_new$gene %in% names(gene_freq[gene_freq==1]),]
write.csv(d3_unique,"./final_files_rev1/supp_table_d3_de_unique_genes_rev1.csv")
d3_unique <- read.csv("./final_files_rev1/supp_table_d3_de_unique_genes_rev1.csv")

View(d3_unique[d3_unique$gene %in% up_genes_rep,])
View(d3_unique[d3_unique$gene %in% up_genes,])


#plot heatmap on a cluster level
d3_unique$cluster <- factor(d3_unique$cluster, levels = c("MEF","Day3 reprogramming","Day3 dead-end"))
mat<- b4rna_sub2[["RNA"]]@data[d3_unique$gene, ] %>% as.matrix()
rownames(mat) <- rownames(d3_unique)
mat <- t(mat)
mat <- as.data.frame(mat)

#sanity check
all(rownames(mat) == rownames(b4rna_sub2@meta.data))
mat$idents <- Idents(b4rna_sub2)

x <- mat %>% group_by(idents) %>% dplyr::summarise_all(mean)
temp<- x$idents
x$idents <- NULL
x <- scale(x)
rownames(x) <- temp


col_fun <- get_circ_palette(viridis(10), -1.2,1.2)
x <- t(x)
x <- as.data.frame(x)
colnames(x) <- c('D3-on','D3-off','MEF')
write.csv(x,"./plot_data_rev1/unique_de_genes_heatmap_rev1.csv")
x <- read.csv("./plot_data_rev1/unique_de_genes_heatmap_rev1.csv")

annot_column <- as.data.frame(rownames(x))
colnames(annot_column) <- c('rowname')
annot_column$idx <- rownames(annot_column)
annot_column$gene <- d3_unique[annot_column$rowname,]$gene
annot_column$de_annot <- 'none'
annot_column[annot_column$gene %in% up_genes,]$de_annot <- 'up'
annot_column[annot_column$gene %in% down_genes,]$de_annot <- 'down'

select_genes <- c("Krt19","Wnt4","Ezr","Acta2",'Anxa8',"Tagln","Ncam1",'Lgals1','S100a4','Vim')
hb = rowAnnotation(link=anno_mark(at=as.numeric(annot_column[annot_column$gene %in% select_genes,]$idx),
                                  labels = annot_column[annot_column$gene %in% select_genes,]$gene,
                                  labels_gp = gpar(fontsize = 13, fontface='italic'), padding = unit(1, "mm")))

ComplexHeatmap::Heatmap(x[rownames(d3_unique[d3_unique$gene %in% down_genes,]),], show_row_names = FALSE, cluster_columns = FALSE, cluster_rows = FALSE, col = col_fun)

#  + Heatmap(as.matrix(annot_column$de_annot), width = unit(5,'mm'), col = c('none' = 'white','up' = 'green','down' = 'red'))
pdf("./plots_rev1/unique_de_genes_heatmap_rev1.pdf", bg = "transparent", width = 3, height = 4)
ComplexHeatmap::Heatmap(x, show_row_names = FALSE, cluster_columns = TRUE, column_dend_reorder = c(1,2,3),column_names_rot = 45, column_names_centered = FALSE,
                        cluster_rows = FALSE, col = col_fun,border = TRUE,border_gp = gpar(lwd=1)) + hb 
dev.ofF()



#Find markers at each day
Idents(b4rna) <- b4rna$md_fate_rev1
m1 <- FindMarkers(b4rna, ident.1 = "Day12 dead-end", ident.2 = "Day3 dead-end",
                  max.cells.per.ident = 1000, only.pos = TRUE, logfc.threshold = 0.1)
m2 <- FindMarkers(b4rna, ident.1 = "Day21 dead-end", ident.2 = "Day12 dead-end",
                  max.cells.per.ident = 1000, only.pos = TRUE, logfc.threshold = 0.1)

m1_d <- FindMarkers(b4rna, ident.1 = "Day3 dead-end", ident.2 = "Day12 dead-end",
                    max.cells.per.ident = 1000, only.pos = TRUE, logfc.threshold = 0.1)
m2_d <- FindMarkers(b4rna, ident.1 = "Day12 dead-end", ident.2 = "Day21 dead-end",
                    max.cells.per.ident = 1000, only.pos = TRUE, logfc.threshold = 0.1)


m3 <- FindMarkers(b4rna, ident.1 = "Day12 reprogramming", ident.2 = "Day3 reprogramming",
                  max.cells.per.ident = 1000, only.pos = TRUE, logfc.threshold = 0.1)
m4 <- FindMarkers(b4rna, ident.1 = "Day21 reprogramming", ident.2 = "Day12 reprogramming",
                  max.cells.per.ident = 1000, only.pos = TRUE, logfc.threshold = 0.1)

m3_d <- FindMarkers(b4rna, ident.1 = "Day3 reprogramming", ident.2 = "Day12 reprogramming",
                    max.cells.per.ident = 1000, only.pos = TRUE, logfc.threshold = 0.1)
m4_d <- FindMarkers(b4rna, ident.1 = "Day12 reprogramming", ident.2 = "Day21 reprogramming",
                    max.cells.per.ident = 1000, only.pos = TRUE, logfc.threshold = 0.1)