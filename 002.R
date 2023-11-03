library(EnrichedHeatmap)
library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(circlize)
library(Guitar)
library(ChIPseeker)
setwd("e:/K9homer/")
bedfile<-readPeakFile("old/mergeLK9KO4_LWT15_new.bed")

bw.files <- c("d:/LXHjiace/Bw_R_R/L_K9_KO4_new3.bw",
              "e:/XGY0918/Bw/L_WT15_new3.bw")

sample <- c("L_K9KO4_new3","L_WT15_new3")

bg_col <- c('#660303','#031d66')

ht_list <- lapply(seq_along(bw.files), function(x){
  bw <- import.bw(bw.files[x])
  
  mat_trim = normalizeToMatrix(bw, bedfile, value_column = "score",
                               extend = 5000, mean_mode = "w0", w = 20,
                               # keep = c(0, 0.99),
                               background = 0, smooth = TRUE)
  
  col_fun = colorRamp2(quantile(mat_trim, c(0, 0.99)), c("grey95", bg_col[x]))
  # col_fun = colorRamp2(c(0,30), c("grey95", bg_col[x]))
  
  # plot
  ht <-
    EnrichedHeatmap(mat_trim, col = col_fun, name = sample[x],
                    column_title = sample[x],
                    column_title_gp = gpar(fontsize = 10,
                                           fill = bg_col[x]),
                    # top_annotation = HeatmapAnnotation(
                    #   enriched = anno_enriched(
                    #     ylim = c(0, 30))),
                    use_raster = F,
                    # show_heatmap_legend = lgd[x],
                    show_heatmap_legend = T
    )
  
  return(ht)
})
##### unfinished