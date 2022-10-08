library(tidyverse)
source("/data/scripts/utility.R")

meta = read_tsv("../../../VMT_Con_metadata_infant_stool_metabolome.txt") 

pc = read.pc("bray_curtis_pc.txt")
PC = pc$pc_data %>% as.data.frame() %>% select(PC1:PC5) %>% 
  rownames_to_column(var = "#SampleID")
proportion = pc$pc_percent

PC_plot = meta %>% select(`#SampleID`, Day, Time_point, Group) %>% 
  inner_join(PC, by = "#SampleID") %>% 
  mutate(Day = fct_relevel(Day, c("D3", "D7", "D30", "D42")))


ggplot(PC_plot, aes(PC1, PC2)) + 
  geom_point(aes(color = Time_point, shape = Group)) +
  theme_minimal()+
  scale_color_distiller(palette = "Spectral", direction = 1, guide = guide_legend(), label = c("D3", "D7", "D30", "D42")) + 
  scale_shape_manual(values = c(15, 16)) + 
  labs(color = "Day",
       x = paste0("PC1 (", proportion[1], "%)"), 
       y = paste0("PC2 (", proportion[2], "%)")) +
  theme(aspect.ratio = 1)
ggsave("PCoA.pdf", width = 4, height = 4)




##############


coor_summary = PC_plot %>% group_by(Group, Day) %>% 
  summarise(PC1.se = sd(PC1)/sqrt(n()), 
            PC2.se = sd(PC2)/sqrt(n()), 
            PC1_centroid = mean(PC1), 
            PC2_centroid = mean(PC2), 
            PC1.sd = sd(PC1), 
            PC2.sd = sd(PC2)) 
ggplot() + 
  geom_errorbar(data = coor_summary, mapping = aes(PC1_centroid, PC2_centroid, ymin = PC2_centroid - PC2.se, ymax = PC2_centroid + PC2.se, color = Day),width=.01, show.legend = F) + 
  geom_errorbarh(data = coor_summary, mapping = aes(PC1_centroid, PC2_centroid, xmin = PC1_centroid - PC1.se, xmax = PC1_centroid + PC1.se, color = Day),height=.005, show.legend = F) + 
  geom_point(data = coor_summary, mapping = aes(PC1_centroid, PC2_centroid, color = Day, shape = Group), size = 3)+ 
  scale_color_brewer(palette = "Spectral", direction = 1, guide = guide_legend(), label = c("D3", "D7", "D30", "D42")) +
  scale_shape_manual(values = c(21,22)) +
  theme_minimal()+
  theme(aspect.ratio = 1) + 
  labs(x = paste0("PC1 (", proportion[1], "%)"), 
       y = paste0("PC2 (", proportion[2], "%)"))

ggsave("PCoA_errorbar.pdf", width = 4, height = 4)



pal = RColorBrewer::brewer.pal(12, "Paired")
data = read.pc("bray_curtis_pc.txt") 
coor = data$pc_data %>% 
  select(1:3) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "#SampleID") %>% 
  inner_join(meta, by = "#SampleID")
centroid = coor %>% group_by(Group, Day) %>% 
  summarise(PC1_centroid = mean(PC1), PC2_centroid = mean(PC2)) 
coor = coor %>% inner_join(centroid, by = c("Group", "Day")) %>% 
  mutate(Day = fct_relevel(Day, c("D3", "D7", "D30", "D42")), 
         Time_point = as.numeric(Day), 
         color = factor(str_c(Group, Time_point)))

proportion = data$pc_percent
coor_summary = coor %>% group_by(Group, Day) %>% 
  summarise(PC1.se = sd(PC1)/sqrt(n()), 
            PC2.se = sd(PC2)/sqrt(n()), 
            PC1_centroid = mean(PC1), 
            PC2_centroid = mean(PC2), 
            PC1.sd = sd(PC1), 
            PC2.sd = sd(PC2)) %>% 
  mutate(Day = fct_relevel(Day, c("D1", "D2", "D3", "D7", "D30", "D42")), 
         Time_point = as.numeric(Day), 
         color = factor(str_c(Group, Time_point)))
ggplot() + 
  stat_ellipse(coor, mapping = aes(PC1, PC2, fill=color), type = "norm", geom ="polygon", alpha=0.1,show.legend = F, color = NA, level = 0.6) + 
  geom_segment(coor, mapping = aes(x = PC1, y = PC2, xend = PC1_centroid, yend = PC2_centroid, color = color), alpha = 0.3, show.legend = F) + 
  geom_point(coor, mapping = aes(PC1, PC2, color = color), show.legend = F) + 
  geom_point(data = coor_summary, mapping = aes(PC1_centroid, PC2_centroid, fill = color), size = 4, shape = 22, color = "black") + 
  scale_color_manual(values = c(pal[c(1,3,5,7)], pal[c(2,4,6,8)])) + 
  scale_fill_manual(values = c(pal[c(1,3,5,7)], pal[c(2,4,6,8)]),
                    label = c(str_c(c("Con"), c("D3", "D7", "D30", "D42"), sep = "_"), str_c(c("VMT"), c("D3", "D7", "D30", "D42"), sep = "_"))
  )+
  theme_minimal() +
  labs(
    fill = "",
    x = paste0("PC1 (", proportion[1], "%)"),
    y = paste0("PC2 (", proportion[2], "%)")) +
  theme(aspect.ratio = 1, 
        axis.text = element_text(color = "black"))
ggsave("所有样本PCoA.pdf", width = 5, height = 5)


library(scatterplot3d)
library(RColorBrewer)
colormap <- brewer.pal(6, "Spectral")

PC_plot_3d = PC_plot %>% mutate(shape = if_else(Group == "Con", 21, 22),
                                Group = fct_relevel(Group, c("Con", "VMT"))) 
colors = colormap[PC_plot_3d$Time_point]
pdf("3D_PCoA.pdf", width = 6, height = 6)
p = scatterplot3d(PC_plot_3d$PC1, PC_plot_3d$PC2, PC_plot_3d$PC3, 
                  bg = colors, pch = PC_plot_3d$shape, 
                  xlab = paste0("PC1 (", proportion[1], "%)"),
                  ylab = paste0("PC2 (", proportion[2], "%)"), 
                  zlab = paste0("PC3 (", proportion[3], "%)"),
                  grid = T, box = T, 
                  col.axis = "black", 
                  angle = 24, cex.symbols = 1.5)
legend("topright", legend = levels(PC_plot_3d$Group), 
       pch = c(21, 22))
legend("right", legend = levels(PC_plot_3d$Day), col = colormap, pch = 15)
dev.off()

# 每天的
pcoa = function(pc_file, r2, p_value){
  data = read.pc(pc_file) 
  coor = data$pc_data %>% 
    select(1:3) %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "#SampleID") %>% 
    inner_join(meta, by = "#SampleID")
  centroid = coor %>% group_by(Group) %>% 
    summarise(PC1_centroid = mean(PC1), PC2_centroid = mean(PC2)) 
  coor = coor %>% inner_join(centroid, by = "Group")
  proportion = data$pc_percent
  coor_summary = coor %>% group_by(Group) %>% 
    summarise(PC1.se = sd(PC1)/sqrt(n()), 
              PC2.se = sd(PC2)/sqrt(n()), 
              PC1_centroid = mean(PC1), 
              PC2_centroid = mean(PC2), 
              PC1.sd = sd(PC1), 
              PC2.sd = sd(PC2))
  
  ggplot() + 
    stat_ellipse(coor, mapping = aes(PC1, PC2, fill=Group), type = "norm", geom ="polygon", alpha=0.1,show.legend = F, color = NA) + 
    geom_segment(coor, mapping = aes(x = PC1, y = PC2, xend = PC1_centroid, yend = PC2_centroid, color = Group), alpha = 0.6, show.legend = F) + 
    geom_point(coor, mapping = aes(PC1, PC2, fill = Group), color = "black", shape = 21, show.legend = T) + 
    geom_point(data = coor_summary, mapping = aes(PC1_centroid, PC2_centroid, fill = Group), size = 3, shape = 22, color = "black") + 
    scale_fill_manual(values = c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467BC'))+
    scale_color_manual(values = c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467BC'))+
    annotate(geom = "text", x = Inf, y = Inf, label = paste0("R2=", r2, "\n", "P = ", p_value), vjust = 1.1, hjust = 1.1) + 
    theme_bw() +
    labs(
      x = paste0("PC1 (", proportion[1], "%)"),
      y = paste0("PC2 (", proportion[2], "%)")) +
    theme(aspect.ratio = 1, 
          axis.text = element_text(color = "black"))
  
}

D3 = pcoa("D3_bdiv/bray_curtis_pc.txt", 0.007, 0.86)
D7 = pcoa("D7_bdiv/bray_curtis_pc.txt", 0.006, 0.96)
D30 = pcoa("D30_bdiv/bray_curtis_pc.txt", 0.017, 0.43)
D42 = pcoa("D42_bdiv/bray_curtis_pc.txt", 0.013, 0.70)

library(patchwork)
P = D1 + D2 + D3 +D7 + D30 + D42 + plot_layout(guides = 'collect')
P

ggsave("各时间点PCoA_bray_curtis.pdf", width = 9, height = 6)



library(ggridges)
library(ggpubr)
pvalue1 = PC_plot %>% group_by(Day) %>% 
  rstatix::wilcox_test(PC1 ~ Group)
pvalue2 = PC_plot %>% group_by(Day) %>% 
  rstatix::wilcox_test(PC2 ~ Group)

p1 = ggplot(PC_plot, aes(PC1, y = Day)) + 
  geom_density_ridges2(aes(fill=Group), alpha = 0.6) +
  theme_ridges() + 
  theme(axis.title.x = element_text(hjust = 0.5), 
        axis.title.y = element_text(hjust = 0.5)) + 
  geom_text(data = pvalue1, mapping = aes(y = Day, label = paste0("P=", p)), x = 0.4, hjust = .5, vjust = -1.1) + 
  #scale_fill_manual(values = c('#1f77b4', '#ff7f0e')) + 
  scale_fill_manual(values = c("#7BC2CE", "#8562CC")) + 
  xlab("Bray-Curtis (PCoA1)")
p2 = ggplot(PC_plot, aes(PC2, y = Day)) + 
  geom_density_ridges2(aes(fill=Group), alpha = 0.6) +
  theme_ridges() + 
  theme(axis.title.x = element_text(hjust = 0.5), 
        axis.title.y = element_text(hjust = 0.5)) + 
  geom_text(data = pvalue2, mapping = aes(y = Day, label = paste0("P=", p)), x = 0.5, hjust = .5, vjust = -1.1) + 
  #scale_fill_manual(values = c('#1f77b4', '#ff7f0e')) + 
  scale_fill_manual(values = c("#7BC2CE", "#8562CC")) + 
  xlab("Bray-Curtis (PCoA2)")
ggarrange(p1, p2, common.legend = T)
ggsave("PC1-2峰峦图-color2.pdf", width = 8,height = 4)











library(vegan)
bc_matrix = read.table("bray_curtis_dm.txt", row.names = 1, header = T)
meta = read_tsv("../../VMT_Con_metadata_infant_stool_metabolome.txt")

sample = intersect(colnames(bc_matrix), meta$`#SampleID`)
meta = meta %>% column_to_rownames(var = "#SampleID") %>% 
  .[sample,]

bc_matrix = bc_matrix[sample, sample]

permanova = adonis(bc_matrix ~ Group, data = meta)
permanova$aov.tab
