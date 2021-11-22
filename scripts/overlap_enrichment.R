library(tidyverse)
library(GeneOverlap)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
library(STRINGdb)
library(ggvenn)
library(ggpubr)

en <- read_tsv("./data/overlap_enrichment.tsv")
en <- en %>% filter(`background gene count`<1000)
en$genes <- sapply(en$`matching proteins in your network (labels)`, function(x) strsplit(x, split = ","))
en <- en %>% dplyr::select(`term description`,genes) %>% unnest(genes)
en$on <- 1

mat <- reshape2::acast(en, `term description`~genes, value.var = "on")
mat[is.na(mat)] <- 0

col_fun <- colorRamp2(c(0,1),c("white","grey50"))

pdf(file = paste0("./output/overlap_enrichment.pdf"), height = 4, width = 7) 
Heatmap(mat,
        col = col_fun,
        rect_gp = gpar(type = "none"),
        row_dend_side = "right",
        column_dend_side = "bottom",
        row_names_side = "left",
        column_names_side = "top",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "black", fill = NA, lty = 3, lwd = 1))
          grid.circle(x = x, y = y, r = mat[i, j]*0.4 * min(unit.c(width, height)),gp = gpar(fill = col_fun(mat[i, j]), col = NA))
        },
        border = TRUE)
dev.off()
