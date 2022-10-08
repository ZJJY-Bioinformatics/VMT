library(tidyverse)
source("/data/scripts/utility.R")

meta = read_tsv("../../VMT_Con_metadata_infant_stool_microbiome.txt") %>% 
  filter(Group != "VD")
pc = read.pc("bdiv/bray_curtis_pc.txt")
PC = pc$pc_data %>% as.data.frame() %>% select(PC1:PC5) %>% 
  rownames_to_column(var = "#SampleID")
proportion = pc$pc_percent

PC_plot = meta %>% select(`#SampleID`, Day, Time_point, Group) %>% 
  inner_join(PC, by = "#SampleID") %>% 
  mutate(Day = fct_relevel(Day, c("D3", "D7", "D30", "D42")))

## 三组样本JSD分成两型
ET = read_tsv("../enterotype/2_cluster.sample.txt") %>% 
  select(`#SampleID`, ET)
ET_plot = ET %>% inner_join(PC_plot, by = "#SampleID") %>% 
  mutate(ET = as.factor(ET))

pal = RColorBrewer::brewer.pal(12, "Paired")
ggplot(ET_plot, aes(PC1, PC2)) + 
  geom_point(aes(color = ET, shape = Group)) +
  theme_minimal() + 
  stat_ellipse(aes(PC1, PC2, fill = ET), geom = "polygon", show.legend = F, alpha = 0.1) +
  scale_fill_manual(values = c("#7bc2c4", "#8562cc")) + 
  scale_color_manual(values = c("#7bc2c4", "#8562cc")) + 
  scale_shape_manual(values = c(15, 16)) + 
  labs(
       x = paste0("PC1 (", proportion[1], "%)"), 
       y = paste0("PC2 (", proportion[2], "%)")) +
  theme(aspect.ratio = 1, 
        panel.border = element_rect(color = "black", fill = "transparent")) 

ggsave("fig1b.pdf", width =4, height = 3)



# adonis by time
library(vegan)
bc = read.table("bdiv/bray_curtis_dm.txt", row.names = 1) %>% as.dist()
group = ET_plot %>% column_to_rownames(var = "#SampleID")
adonis2(bc ~ ET, group, permutations = 999)


