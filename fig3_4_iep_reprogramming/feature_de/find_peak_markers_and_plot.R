

b4.atac <- loadArchRProject("../archr_objs/b4_atac_txg/", showLogo = FALSE)


#find marker peaks for each Day 3 imputed grp and MEFs
markersMEF <- getMarkerFeatures(
  ArchRProj = b4.atac, 
  useMatrix = "PeakMatrix",
  groupBy = "md_fate_impute_rev1",
  useGroups=c("MEFs"),
  bgdGroups=c('Day3 reprogramming','Day3 dead-end'),  
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  binarize=TRUE
)

markersRep <- getMarkerFeatures(
  ArchRProj = b4.atac, 
  useMatrix = "PeakMatrix",
  groupBy = "md_fate_impute_rev1",
  useGroups=c('Day3 reprogramming'),
  bgdGroups=c('Day3 dead-end','MEFs'),  
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  binarize=TRUE
)

markersDead <- getMarkerFeatures(
  ArchRProj = b4.atac, 
  useMatrix = "PeakMatrix",
  groupBy = "md_fate_impute_rev1",
  useGroups=c('Day3 dead-end'),
  bgdGroups=c('Day3 reprogramming','MEFs'),  
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  binarize=TRUE
)

saveRDS(markersMEF, "../proc_files/DARs/MEF_DARs_rev1.Rds")
saveRDS(markersDead, "../proc_files/DARs/D3_DE_DARs_rev1.Rds")
saveRDS(markersRep, "../proc_files/DARs/D3_Rep_DARs_rev1.Rds")

#save as bed files
markersMEF <- readRDS("../proc_files/DARs/MEF_DARs_rev1.Rds")
markersDead <- readRDS("../proc_files/DARs/D3_DE_DARs_rev1.Rds")
markersRep <- readRDS("../proc_files/DARs/D3_Rep_DARs_rev1.Rds")

markersMEF.gr <- getMarkers(markersMEF, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)
markersRep.gr <- getMarkers(markersRep, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)
markersDead.gr <- getMarkers(markersDead, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)

markersMEF.gr$MEFs <- markersMEF.gr$MEFs[order(markersMEF.gr$MEFs$Log2FC, decreasing = TRUE)]
markersRep.gr$`Day3 reprogramming` <- markersRep.gr$`Day3 reprogramming`[order(markersRep.gr$`Day3 reprogramming`$Log2FC, decreasing = TRUE)]
markersDead.gr$`Day3 dead-end` <- markersDead.gr$`Day3 dead-end`[order(markersDead.gr$`Day3 dead-end`$Log2FC, decreasing = TRUE)]



export.bed(markersMEF.gr$MEFs,"../proc_files/DARs/MEF_DARs_rev1.bed")
export.bed(markersRep.gr$`Day3 reprogramming`,"../proc_files/DARs/D3_Rep_DARs_rev1.bed")
export.bed(markersDead.gr$`Day3 dead-end`,"../proc_files/DARs/D3_DE_DARs_rev1.bed")

#plot heatmap
source("./highlight_archr.R")

markersMEF = readRDS("../proc_files/DARs/MEF_DARs_rev1.Rds")
markersDead = readRDS("../proc_files/DARs/D3_DE_DARs_rev1.Rds")
markersRep = readRDS("../proc_files/DARs/D3_Rep_DARs_rev1.Rds")

#create DAR column annot
markersMEF.gr <- import.bed("../proc_files/DARs/MEF_DARs_rev1.bed")
markersRep.gr <- import.bed("../proc_files/DARs/D3_Rep_DARs_rev1.bed")
markersDead.gr <- import.bed("../proc_files/DARs/D3_DE_DARs_rev1.bed")


df1 <- as.data.frame(cbind(unlist(paste0(seqnames(markersMEF.gr),":",start(markersMEF.gr),"-",end(markersMEF.gr))),rep('MEF',length(markersMEF.gr))))
df2 <- as.data.frame(cbind(unlist(paste0(seqnames(markersRep.gr),":",start(markersRep.gr),"-",end(markersRep.gr))),rep('D3-Rep',length(markersRep.gr))))
df3 <- as.data.frame(cbind(unlist(paste0(seqnames(markersDead.gr),":",start(markersDead.gr),"-",end(markersDead.gr))),rep('D3-DE',length(markersDead.gr))))

df_all <- rbind(df1,df2,df3)

dars <- cbind(markersRep,markersDead,markersMEF)
rownames(dars) <- paste0(rowData(dars)[,1],":" ,rowData(dars)[,3],"-" ,rowData(dars)[,4])
dorc_peaks <- read.csv("../../multimodal_analysis/proc_files/dorc_calling/dorc_filt_p2g_list.csv")

x <- getPeak2GeneLinks(b4.atac)

x$Peak2GeneLinks


df <- read.csv("../proc_files/dorc_p2g_loops.bed",sep="\t", header = FALSE)
colnames(df) <- c("chr","start","end","rObs")
df_bed <- makeGRangesFromDataFrame(df)


df_bed$value <- df$rObs
df_bed <- df_bed[df_bed$value > 0.03]
highlight_gr <- GRanges(c("chr6","chrX","chr12"),
                        IRanges(c(91259711-200,164330441-200,45070133-200),
                                c(91260211+200,164330941+200,45070633+200)))

custom_pal <- c("#636363","#54278f","#756bb1","#9e9ac8","#bd0026","#f03b20","#fd8d3c")
names(custom_pal) <- as.character(c(1:7))
custom_pal <- c("Day3 reprogramming"= "#54278f", "Day3 dead-end"= "#bd0026",
                "Day12 reprogramming"= "#756bb1", "Day12 dead-end"= "#f03b20",
                "Day21 reprogramming"= "#9e9ac8", "Day21 dead-end"= "#fd8d3c","MEFs"= "#636363")

#
p1 <- plotBrowserTrack(
  ArchRProj = b4.atac, 
  groupBy = "md_fate_rev1", pal = custom_pal,
  useGroups = c("MEFs","Day3 reprogramming","Day12 reprogramming", "Day21 reprogramming",
                "Day3 dead-end","Day12 dead-end", "Day21 dead-end"),
  geneSymbol = c("Vegfd"),downstream = 7000,
  plotSummary = c("bulkTrack",'loopTrack','geneTrack','featureTrack'), threads = 2,
  features = highlight_gr,
  loops = df_bed, sizes = c(5,0.5,0.5,0.5))
grid::grid.newpage()
grid::grid.draw(p1$Vegfd)
plotPDF(p1, name = "vegfd_encode_enh_rev1", width = 4.5, height = 3, ArchRProj = b4.atac, addDOC = FALSE)

p2 <- plotBrowserTrack(
  ArchRProj = b4.atac, 
  groupBy = "md_fate_rev1", pal = custom_pal,
  useGroups = c("MEFs","Day3 reprogramming","Day12 reprogramming", "Day21 reprogramming",
                "Day3 dead-end","Day12 dead-end", "Day21 dead-end"),
  geneSymbol = c("Col28a1"),downstream = 7000,
  plotSummary = c("bulkTrack",'loopTrack','geneTrack','featureTrack'), threads = 2,
  features = highlight_gr,
  loops = df_bed, sizes = c(5,0.5,0.5,0.5))
grid::grid.newpage()
grid::grid.draw(p2$Col28a1)
plotPDF(p2, name = "col28a1_encode_enh_rev1", width = 4, height = 3, ArchRProj = b4.atac, addDOC = FALSE)

p3 <- plotBrowserTrack(
  ArchRProj = b4.atac, 
  groupBy = "md_fate_rev1", pal = custom_pal,
  useGroups = c("MEFs","Day3 reprogramming","Day12 reprogramming", "Day21 reprogramming",
                "Day3 dead-end","Day12 dead-end", "Day21 dead-end"),
  geneSymbol = c("Aox3"),downstream = 40000, upstream = 1000,
  plotSummary = c("bulkTrack",'loopTrack','geneTrack','featureTrack'), threads = 2,
  features = highlight_gr,
  loops = df_bed, sizes = c(5,0.5,0.5,0.5))
grid::grid.newpage()
grid::grid.draw(p3$Aox3)
plotPDF(p3, name = "Aox3_encode_enh_rev1", width = 4, height = 3, ArchRProj = b4.atac, addDOC = FALSE)

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = dars, 
  cutOff = "FDR <= 0.05 & Log2FC >= 1", returnMatrix = TRUE,
  transpose = TRUE, labelMarkers = c("chr6:91259711-91260211", "chr12:45070133-45070633"), clusterCols = FALSE
)

heatmapPeaks_2 <- plotMarkerHeatmap(
  seMarker = dars, 
  cutOff = "FDR <= 0.05 & Log2FC >= 1", returnMatrix = FALSE,
  transpose = FALSE, clusterCols = FALSE,
)

heatmapPeaks <- t(heatmapPeaks)

#get DORC annotation for all DAR peaks
dar_peaks <- as.data.frame(gsub(":","-",rownames(heatmapPeaks)))
colnames(dar_peaks) <- c("dar")

#add DORC info to DARs
dar_peaks <- merge(dar_peaks, dorc_peaks, by.x='dar', by.y='Peak', all.x=TRUE)
rownames(dar_peaks) <- dar_peaks$dar
dar_peaks$bool <- 'NO'
dar_peaks$bool[!is.na(dar_peaks$Gene)] <- 'YES'
dar_peaks <- dar_peaks[gsub(":","-",rownames(heatmapPeaks)),]

#reorder to ensure the correct order is preserved
all(dar_peaks$dar == gsub(":","-",rownames(heatmapPeaks)))

dar_peaks$idx <- c(1:dim(dar_peaks)[[1]])

dar_peaks_sub <- dar_peaks[!is.na(dar_peaks$Gene),]
# dar_peaks_sub <- dar_peaks_sub[(!rev(duplicated(rev(dar_peaks_sub$Gene)))) | (dar_peaks_sub$idx %in% c(16009,2757)),]

#select which genes-peaks to plot
#rep_genes: 'Abhd2','Alb','Creb3l2','Stxbp6','Myo5c','Foxq1'
#de_genes: 'Fbln2','Ncam1','Vegfd','Cd34','Acvrl1', 'Wnt5a'
#mef_genes: 'Col1a1','Col1a2','Serpine1','Thbs1'

plotgenes <- c('Abhd2','Alb','Creb3l2','Myo5c','Foxq1',
               'Ncam1','Vegfd','Cd34','Fbln2',
               'Col1a2','Serpine1','Thbs1', 'Stxbp6')

# plotgenes <- c("Stxbp6")
# c("chr6:91259711-91260211", "chr12:45070133-45070633")
# dar_peaks_sub <- dar_peaks_sub[(dar_peaks_sub$Gene %in% plotgenes) | (dar_peaks_sub$idx %in% c(16009,2757)),]

#subset to peaks linked to genes of interest
dar_peaks_sub <- dar_peaks_sub[(dar_peaks_sub$Gene %in% plotgenes),]

dar_peaks_sub_unique <- list()
for (gene_curr in plotgenes){
  mtx_curr <- dar_peaks_sub[dar_peaks_sub$Gene==gene_curr,]
  highest_val <- mtx_curr[order(mtx_curr$rObs, decreasing = FALSE)[1],]
  dar_peaks_sub_unique[[gene_curr]] <- highest_val
}
dar_peaks_sub_unique <- do.call(rbind, dar_peaks_sub_unique)


#create matrix for DAR label
# rownames(df_all) <- df_all$V1

# df_all <- df_all[colnames(heatmapPeaks),]
# heatmap_data <- 1*(gsub(":","-",colnames(heatmapPeaks)) %in% dorc_peaks$Peak)
# ha = rowAnnotation(bar =as.matrix(df_all$V2))
hb = rowAnnotation(link=anno_mark(at=dar_peaks_sub_unique$idx,
                                  labels = dar_peaks_sub_unique$Gene,
                                  labels_gp = gpar(fontsize = 13, fontface='italic'), padding = unit(1, "mm")))


col_fun <- get_circ_palette(jdb_palettes$solar_extra, -1.15,1.15)
colnames(heatmapPeaks) <- c('D3-On','D3-Off','MEF')
write.csv(heatmapPeaks,"../proc_files/plot_data_rev1/dar_heatmap_genenames_rev1.csv")
heatmapPeaks <- read.csv("../proc_files/plot_data_rev1/dar_heatmap_genenames_rev1.csv", row.names = 1)

#removed
#+Heatmap(as.matrix(df_all$V2), width = unit(5,'mm'), col = c('D3-Rep' = '#bdd7f0','De-DE' = '#f9c6cf','MEF' = '#b9dc98'))


pdf("../proc_files/plots_rev1/dar_heatmap_genenames_rev1.pdf", bg = "transparent", width = 4, height = 6)
Heatmap(heatmapPeaks, cluster_rows = FALSE, cluster_columns = TRUE, column_dend_reorder = c(1,2,3),
        column_names_rot = 0, column_names_centered = TRUE,
        show_row_names = FALSE, height = 3, col = col_fun, border = TRUE,border_gp = gpar(lwd=1)) + Heatmap(as.matrix(dar_peaks$bool), width = unit(5,'mm'), col = c('YES' = 'red','NO' = 'white'),border = TRUE,border_gp = gpar(lwd=1)) + hb
dev.ofF()