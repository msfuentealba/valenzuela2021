library(tidyverse)
orto <- read_rds("./data/ortologues_human_mouse.rds")[,3:4] %>% na.omit %>% unique %>% set_names(c("symbol","mouse"))

hump <- readxl::read_xlsx("./data/humphrey_2021.xlsx", sheet = 3, skip = 1)[,c(1,3,7,11)]
hump$cervical.lfc <- as.numeric(hump$cervical.lfc)
hump$thoracic.lfc <- as.numeric(hump$thoracic.lfc)
hump$lumbar.lfc <- as.numeric(hump$lumbar.lfc)
#hump <- hump %>% group_by(genename) %>% mutate_all(as.numeric)
colnames(hump)[1] <- "symbol"
hump <- hump %>% left_join(orto)
hump <- hump[,2:5] %>% group_by(mouse) %>% summarise_all(mean)
colnames(hump) <- c("mouse","Cervical","Thoracic","Lumbar")

data <- read_rds("./data/differential_abundance.rds")
data <- data[!grepl("Krt",data$symbol),]
g4g2_up <- data %>% filter(group=="g4_g2"&fdr<0.05&logfc>0) %>% ungroup %>% dplyr::select(symbol) %>% unlist %>% as.character()
g3g1_up <- data %>% filter(group=="g3_g1"&fdr<0.05&logfc>0) %>% ungroup %>% dplyr::select(symbol) %>% unlist %>% as.character()
inter <- intersect(g4g2_up,g3g1_up)

rehump <- reshape2::melt(hump) %>% na.omit
rehump <- rehump %>% filter(mouse%in%unique(c(g4g2_up,g3g1_up)))
rehump$group <- ifelse(rehump$mouse%in%inter,"inter",ifelse(rehump$mouse%in%g4g2_up, "g4g2_up", "g3g1_up"))
rehump$group <- factor(rehump$group, levels = c("g4g2_up","inter","g3g1_up"))

library(RColorBrewer)
pdf(file = paste0("./output/human_comparison.pdf"), height = 4, width = 3) 
ggplot(rehump, aes(x=variable, y=value, color = group)) + 
  ylab("Log 2 fold change (ALS patients vs controls)")+
  xlab("Spinal cord region")+
  ylim(-1.2,1.35)+
  geom_hline(yintercept = 0, linetype="dotted", color = "black", size = 0.5) +
  #stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="crossbar", width=0.5) +
  geom_jitter(position=position_jitter(0.25), shape = 16, size = 1.5) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="pointrange", color="black")+
  theme_bw() +
  #ylim(-0.7,0.7)+
  scale_color_manual(values = brewer.pal(3, "Set3")[1:3]) +
  #coord_capped_cart(bottom='both', left='both') +
  theme(legend.position = "none")
dev.off()
    
t.test(rehump$value[rehump$variable=="Cervical"])$p.value
t.test(rehump$value[rehump$variable=="Thoracic"])$p.value
t.test(rehump$value[rehump$variable=="Lumbar"])$p.value


