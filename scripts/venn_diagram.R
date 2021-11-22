library(tidyverse)
library(GeneOverlap)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
library(STRINGdb)
library(ggvenn)
library(ggpubr)

data <- read_rds("./data/differential_abundance.rds")
data <- data[!grepl("Krt",data$symbol),]
g2g1_up <- data %>% filter(group=="g2_g1"&fdr<0.05&logfc>0) %>% ungroup %>% dplyr::select(symbol) %>% unlist %>% as.character()
g2g1_dn <- data %>% filter(group=="g2_g1"&fdr<0.05&logfc<0) %>% ungroup %>% dplyr::select(symbol) %>% unlist %>% as.character()
g4g2_up <- data %>% filter(group=="g4_g2"&fdr<0.05&logfc>0) %>% ungroup %>% dplyr::select(symbol) %>% unlist %>% as.character()
g4g2_dn <- data %>% filter(group=="g4_g2"&fdr<0.05&logfc<0) %>% ungroup %>% dplyr::select(symbol) %>% unlist %>% as.character()

testGeneOverlap(newGeneOverlap(g2g1_up,g4g2_dn,genome.size = length(unique(data$symbol))))
testGeneOverlap(newGeneOverlap(g2g1_dn,g4g2_up,genome.size = length(unique(data$symbol))))
intersect(g2g1_up,g4g2_dn)
#write.table(intersect(g2g1_dn,g4g2_up), file = "./output/overlap_dn_up.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)

rb_pal <- brewer.pal(11, "RdBu")

p1 <- ggvenn(list(`Tg Mock /\nNon-Tg Mock` = g2g1_up, `Tg XBP1s /\nTg Mock` = g4g2_dn),
             show_percentage = FALSE,
             stroke_color = c(rep(rb_pal[2],100),rep(rb_pal[10],100)),
             fill_color = c("white","white"),
             stroke_size = 1, set_name_size = 4) 

p2 <- ggvenn(list(`Tg Mock /\nNon-Tg Mock` = g2g1_dn, `Tg XBP1s /\nTg Mock` = g4g2_up),
             show_percentage = FALSE,
            
             stroke_color = c(rep(rb_pal[10],100),rep(rb_pal[2],100)),
             fill_color = c("white","white"),
             stroke_size = 1, set_name_size = 4) 

pdf(file = paste0("./output/venn_diagrams.pdf"), height = 3.5, width = 5.5) 
ggarrange(p1,p2,ncol = 2)
dev.off()

g3g1_up <- data %>% filter(group=="g3_g1"&fdr<0.05&logfc>0) %>% ungroup %>% dplyr::select(symbol) %>% unlist %>% as.character()
g3g1_dn <- data %>% filter(group=="g3_g1"&fdr<0.05&logfc<0) %>% ungroup %>% dplyr::select(symbol) %>% unlist %>% as.character()
intersect(g3g1_dn,g4g2_dn)
intersect(g3g1_up,g4g2_up)
#write.table(g3g1_up, file = "./output/top_nontg_xbp1.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
#write.table(g4g2_up, file = "./output/top_tg_xbp1.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
#write.table(union(g3g1_up,g4g2_up), file = "./output/top_both_xbp1.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
testGeneOverlap(newGeneOverlap(g3g1_up,g4g2_up,genome.size = length(unique(data$symbol))))

pdf(file = paste0("./output/venn_xbp1_up.pdf"), height = 3, width = 5) 
ggvenn(list(`Non-Tg XBP1s /\nNon-Tg Mock` = g3g1_up, `Tg XBP1s /\nTg Mock` = g4g2_up),
       show_percentage = FALSE,
       #stroke_color = c(rep(brewer.pal(3, "Set3")[2],100),rep(rb_pal[2],100)),
       fill_color = c(brewer.pal(3, "Set3")[3],brewer.pal(3, "Set3")[1]),
       stroke_size = 1, set_name_size = 4) 
dev.off()
