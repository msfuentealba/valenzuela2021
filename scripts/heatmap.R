library(tidyverse)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)

data <- read_rds("./data/differential_abundance.rds")
data <- data[!grepl("Krt",data$symbol),]
mat <- reshape2::dcast(data, symbol~group, value.var = "logfc")
rownames(mat) <- mat$symbol
mat$symbol <- NULL
mat <- mat %>% dplyr::select(g3_g1,g2_g1,g4_g2)

r <- cor(mat, method = "s")
cor.test(mat$g3_g1,mat$g2_g1, method = "s")$p.value # r = 0.3
cor.test(mat$g2_g1,mat$g4_g2, method = "s")$p.value # r = 0.3
cor.test(mat$g3_g1,mat$g4_g2, method = "s")$p.value # r = 0.3


col_fun <- colorRamp2(seq(-1,1,2/10),rev(brewer.pal(11, "RdBu")))
min(mat)
colnames(mat) <- c("Non-Tg XBP1s / Non-Tg Mock","Tg Mock / Non-Tg Mock", "Tg XBP1s / Tg Mock")
pdf(file = paste0("./output/heatmap.pdf"), height = 8, width = 2) 
Heatmap(as.matrix(mat),
        border = TRUE,
        cluster_columns = FALSE,
        cluster_rows = TRUE,
        show_row_names = FALSE,
        show_column_names = TRUE,
        heatmap_legend_param = list(direction = "vertical", title = "Log2FC", title_position = "topcenter"),
        col = col_fun)
dev.off()
