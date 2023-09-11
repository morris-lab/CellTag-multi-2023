library(ArchR)
addArchRGenome("mm10")
b4.atac <- loadArchRProject("../archr_objs/b4_atac_txg/", showLogo = FALSE)

#find marker peaks for each Day 3 imputed grp and MEFs
markersMEF <- getMarkerFeatures(
  ArchRProj = b4.atac, 
  useMatrix = "MotifMatrix",
  groupBy = "md_fate_impute_rev1",
  useGroups=c("MEFs"),
  bgdGroups=c('Day3 reprogramming','Day3 dead-end'), 
  useSeqnames='z',
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  binarize=TRUE
)

markersRep <- getMarkerFeatures(
  ArchRProj = b4.atac, 
  useMatrix = "MotifMatrix",
  groupBy = "md_fate_impute_rev1",
  useGroups=c('Day3 reprogramming'),
  bgdGroups=c('Day3 dead-end','MEFs'), 
  useSeqnames='z',
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  binarize=TRUE
)

markersDead <- getMarkerFeatures(
  ArchRProj = b4.atac, 
  useMatrix = "MotifMatrix",
  groupBy = "md_fate_impute_rev1",
  useGroups=c('Day3 dead-end'),
  bgdGroups=c('Day3 reprogramming','MEFs'), 
  useSeqnames='z',
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  binarize=TRUE
)
marker_list <- list(markersMEF, markersDead, markersRep)

saveRDS(marker_list,"../proc_files/Day3_tf_markers_md_fate_impute_rev1.Rds")

#plot heatmap for TF activity markers
tf_markers <- readRDS("../proc_files/Day3_tf_markers_md_fate_impute_rev1.Rds")
mef_markers <- as.data.frame(getMarkers(tf_markers[[1]], cutOff = "FDR <= 0.05 & MeanDiff > 0.5")$MEFs)
de_markers <- as.data.frame(getMarkers(tf_markers[[2]], cutOff = "FDR <= 0.05 & MeanDiff > 0.5")$`Day3 dead-end`)
rep_markers <- as.data.frame(getMarkers(tf_markers[[3]], cutOff = "FDR <= 0.05 & MeanDiff > 0.5")$`Day3 reprogramming`)

mef_markers$cluster <- 'MEF'
rep_markers$cluster <- 'Day3 Reprogramming'
de_markers$cluster <- 'Day3 deadend'

all_markers <- rbind(rep_markers, de_markers,mef_markers)
corr_vals <- read.csv("../proc_files/GSM_MM_correlation_b4_full.csv", row.names = 1)
corr_vals <- corr_vals[(corr_vals$cor > 0.4) | (corr_vals$MotifMatrix_name=='Foxa1_306'),]

tf_freq <- table(all_markers$name)
all_markers <- all_markers[all_markers$name %in% corr_vals$MotifMatrix_name,]
all_markers_unique <- all_markers[all_markers$name %in% names(tf_freq[tf_freq ==1]),]

#load TF activity Z scores for marker TFs
b4.atac.d3 <- loadArchRProject("../archr_objs/b4_atac_d3/")
tf_mtx <- LoadH5Seurat("../proc_files/b4atac_data_mtx/all.cv.h5Seurat", assays='Z_SCORE')
cells_keep <- b4.atac@cellColData$alt_cellNames[b4.atac@cellColData$md_fate_impute_rev1 %in% c("Day3 reprogramming", "Day3 dead-end","MEFs")]
tf_mtx_sub <- tf_mtx[,as.character(cells_keep)]
tf_mtx_sub <- tf_mtx_sub[gsub(pattern = "_", replacement = "-", unique(all_markers$name)),]
tf_mtx_sub <- t(tf_mtx_sub@assays$Z_SCORE@counts)
metadata <- as.data.frame(b4.atac@cellColData$md_fate_impute_rev1, col.names = c("annot"))
rownames(metadata) <- b4.atac@cellColData$alt_cellNames
tf_mtx_sub <- as.data.frame(tf_mtx_sub)
tf_mtx_sub$annot <- metadata[rownames(tf_mtx_sub),]

#create matrix for heatmap
mat <- tf_mtx_sub %>% group_by(annot) %>% dplyr::summarise_all(mean)
temp<- mat$annot
mat$annot <- NULL
mat <- scale(mat)
rownames(mat) <- temp

#create heatmap
col_fun <- get_circ_palette(viridis(10,option="A"), -1.2,1.2)
mat <- t(mat)
mat <- as.data.frame(mat)

rownames(mat) <- gsub("-","_",rownames(mat))
colnames(mat) <- c("D3-Off",'D3-On','MEF')

#save data
write.csv(mat, "../proc_files/plot_data_rev1/unique_de_genes_heatmap_rev1.csv")
mat <- read.csv("../proc_files/plot_data_rev1/unique_de_genes_heatmap_rev1.csv", row.names = 1)
colnames(mat) <- c("D3-Off",'D3-On','MEF')
rownames(mat) <- unlist(lapply(rownames(mat), function(x) {str_split(x,"_")[[1]][1]}))

#save data
saveRDS(all_markers_unique,"../proc_files/plot_data_rev1/unique_de_genes_heatmap_rowdata_rev1.Rds")
all_markers_unique <- readRDS("../proc_files/plot_data_rev1/unique_de_genes_heatmap_rowdata_rev1.Rds")
write.csv(all_markers_unique,"../proc_files/plot_data_rev1/unique_de_genes_heatmap_rowdata_rev1.csv")

all_markers_unique2 <- all_markers_unique
all_markers_unique2$name <- unlist(lapply(all_markers_unique2$name, function(x) {str_split(x,"_")[[1]][1]}))


#row annots
hb = rowAnnotation(link=anno_mark(at=as.numeric(which(all_markers_unique2$name %in% c("Foxa1","Foxd2",'Hnf4a','Foxn2','Nr1h3','Fosl1'))),
                                  labels = c("FOXA1","FOXD2",'HNF4A','FOXN2','NR1H3','FOSL1'),
                                  labels_gp = gpar(fontsize = 11, fontface='italic'), padding = unit(1, "mm")))


write.csv(mat[all_markers_unique2$name,], "../proc_files/plot_data_rev1/unique_de_genes_heatmap_final_plot_rev1.csv")

pdf("../proc_files/plots_rev1/unique_de_TF_heatmap_rev1.pdf", bg = "transparent", width = 3.5, height = 4)
ComplexHeatmap::Heatmap(mat[all_markers_unique2$name,], show_row_names = TRUE,column_names_centered = TRUE,
                        cluster_columns = TRUE, cluster_rows = FALSE, col = col_fun, column_names_rot = 0,
                        border = TRUE,border_gp = gpar(lwd=1)) + hb
dev.off()


#plot another heatmap for only dead-end core TFs, but across all days
tfs_to_plot <- all_markers_unique[all_markers_unique$cluster=='Day3 deadend',]$name
cells_keep <- b4.atac@cellColData$alt_cellNames[!b4.atac@cellColData$md_fate_impute_rev1 %in% c("na")]
tf_mtx_sub <- tf_mtx[,as.character(cells_keep)]
tf_mtx_sub <- tf_mtx_sub[gsub(pattern = "_", replacement = "-", tfs_to_plot),]
tf_mtx_sub <- t(tf_mtx_sub@assays$Z_SCORE@counts)
metadata <- as.data.frame(b4.atac@cellColData$md_fate_impute_rev1, col.names = c("annot"))
rownames(metadata) <- b4.atac@cellColData$alt_cellNames
tf_mtx_sub <- as.data.frame(tf_mtx_sub)
tf_mtx_sub$annot <- metadata[rownames(tf_mtx_sub),]

mat <- tf_mtx_sub %>% group_by(annot) %>% dplyr::summarise_all(mean)
temp<- mat$annot
mat$annot <- NULL
mat <- scale(mat)
rownames(mat) <- temp

#create heatmap
col_fun <- get_circ_palette(viridis(10,option="A"), -1.2,1.2)
mat <- t(mat)
mat <- as.data.frame(mat)

rownames(mat) <- gsub("-","_",rownames(mat))
colnames(mat) <- c("D12-Off", "D12-On","D21-Off","D21-On","D3-Off","D3-On","MEF")
rownames(mat) <- unlist(lapply(rownames(mat), function(x) {str_split(string = x, "_")[[1]][1]}))
rownames(mat) <- toupper(rownames(mat))

pdf("../proc_files/plots_rev1/deadend_core_tf_heatmap_rev1.pdf", bg = "transparent", width = 3, height = 3)
ComplexHeatmap::Heatmap(mat[,c("MEF","D3-On","D12-On","D21-On","D3-Off","D12-Off","D21-Off")],
                        show_row_names = TRUE,column_names_centered = FALSE, row_names_side = 'left',
                        cluster_columns = FALSE, cluster_rows = FALSE, col = col_fun, column_names_rot = 90,
                        border = TRUE,border_gp = gpar(lwd=1),column_split = c(1,2,2,2,3,3,3))
dev.off()