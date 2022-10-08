library(tidyverse)
library(scatterpie)
source("/data/scripts/utility.R")
genus = read.table("../taxa/infant_stool_e10k_L6.txt", skip = 1, check.names = F, comment.char = "", sep = "\t", header = T) %>% 
  column_to_rownames(var = "#OTU ID") %>% 
  sjmisc::rotate_df()

# 各物种丰度排序
top <- sort(colSums(genus), decreasing=TRUE)

# 选择丰度前n的物种进行绘制图形
n = 15
need <- names(top)[1:n]
L6_RA <- genus %>% 
  select(all_of(need))

colnames(L6_RA) = taxa_shortname(colnames(L6_RA))


Others <- 1-rowSums(L6_RA)

# topN物种+others丰度表
genus_abdc = cbind(L6_RA, Others) %>% 
  rownames_to_column(var = "#SampleID")


meta = read_tsv("../../../Metadata_infant_stool_microbiome.txt") %>% 
  select(`#SampleID`, ID, Group, Day) %>% 
  arrange(Group) %>% 
  mutate(Day_p = as.numeric(factor(Day, levels=c("D3","D7","D30","D42"))),
         Patient_p = as.numeric(factor(ID, levels=unique(ID)))) %>% 
  filter(Group == "VD")



plotdata = meta %>% inner_join(genus_abdc, by = "#SampleID") 

plotdata %>% distinct(ID, .keep_all = T) %>% group_by(Group) %>% 
  summarise(n = n())

colors <-  c("#FFFFB3","#80B1D3","#B3DE69","#8DD3C7","#4daf4a","#377eb8","#BEBADA","#FB8072","#FDB462","#FCCDE5","#BC80BD","#CCEBC5","#FFED6F","#CD4F39","#BC41A4","#4F94CD","#E41A1C","#00CD66","#CD3278","#CD8A96","#00C5CD","#CDCD00","#CD85CD","#CD853F","#8B5A2B","#5CACEE","#EE5C42","#00EE76","#EE4A8C","#EED8AE","#00E5EE","#EEEE00","#EED2EE","#EE9A49","#E41A1C","#377EB8","#FF6A6A","#87CEFA","#6E8B3D","#FFEBCD","#B2DFEE")
color <- c(colors[1:15], "#D9D9D9")

p3 <- ggplot() + geom_scatterpie(aes(x=Day_p, y=Patient_p), data = plotdata, color=NA, pie_scale = 4.5, 
                                 cols = colnames(plotdata)[7:ncol(plotdata)]) + 
  scale_y_continuous(breaks = plotdata$Patient_p, labels = plotdata$ID, expand = expansion(mult = 0.01)) + 
  scale_x_continuous(breaks = plotdata$Day_p, labels=plotdata$Day, expand = expansion(mult = 0.1)) + 
  theme_bw() + 
  scale_fill_manual(values = color) +
  #geom_hline(yintercept = 36+0.5, linetype=5) + 
  #geom_hline(yintercept = 68+0.5, linetype=5) + 
  theme(panel.grid.minor = element_blank(), legend.title = element_blank(), 
        #legend.position = c(0.9, 0.9), 
        legend.key.size = unit(1.5,"cm"), 
        legend.text = element_text(face = "italic", size = 30), 
        axis.text = element_text(size = 18),) +
  coord_equal() +
  labs(x="Day", y="")
ggsave("VD_pie.pdf", width=10, height = 10, limitsize = FALSE)

library(patchwork)
p1 + p2 + p3 + plot_layout(guides = "collect")
  
ggsave("Pie.pdf", width = 40, height = 10)
