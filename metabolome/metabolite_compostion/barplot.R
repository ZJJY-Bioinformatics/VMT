library(tidyverse)
library(RColorBrewer)
meta = read_tsv("../../Metadata_infant_stool_metabolome.txt") %>% 
  select(`#SampleID`, Group, Day, ID) %>% 
  mutate(Day = fct_relevel(Day, c("D3", "D7", "D30", "D42")), 
         Group = fct_relevel(Group, c("VD", "Con", "VMT"))) 

# 20大类代谢物浓度均值
metabolites_category = read_csv("03.2_AllMet_with_Raw_Metabolite.csv") %>% 
  select(Class, Raw_Metabolite)
metabolome_category = read_tsv("../../../metabolome/metabolome-raw.txt") %>% 
  inner_join(metabolites_category, ., by = c("Raw_Metabolite" = "#OTU ID")) %>% 
  group_by(Class) %>% 
  summarise(across(where(is.numeric), mean)) %>% 
  ungroup() %>% 
  sjmisc::rotate_df(cn = T, rn = "#SampleID")
fct_lv = names(sort(colSums(metabolome_category[,-1]), decreasing = F))

metabolites_category_long = metabolome_category %>% 
  inner_join(meta, ., by = "#SampleID") %>% 
  pivot_longer(-1:-4, names_to = "Class", values_to = "Concentration") %>% 
  mutate(Class = fct_relevel(Class, fct_lv))



cols = rev(colorRampPalette(brewer.pal(9,"Paired"))(20))
#cols = c("#FF0101","#FFAC01","#0193C6","#5CD8FD","#01B958","#96D531","#BB63B2","#FD84DE","#A7A7A7","#E1E1E1","#B35401","#F09146","#019393","#01D8DB","#FFDA0D","#FFFDBA","#FF90C8","#FFE6FF" ,"#FF8D7E","#FFBB6B","#6FBEDC","#A0E6EB","#57D38B","#ACFFC8","#FFFFFF")
cols = c("#FDDC7B","#4169B2","#B1A4C0","#479E9B","#ACCDDC","#DD5F60","#F2B379","#7CB6DB","#EE4A25","#BCEE68","#B1A4C0","#479E9B","#FDDC7B","#4169B2","#ACCDDC","#DD5F60","#F2B379","#7CB6DB","#FA3105","#CD3700", "#BC8F8F","#1C86EE","#A4D3EE","blue")

cols = c("#80B1D3","#B3DE69","#FFFFB3","#8DD3C7","#4daf4a","#377eb8","#BEBADA","#FB8072","#FDB462","#FCCDE5","#BC80BD","#CCEBC5","#FFED6F","#CD4F39","#BC41A4","#4F94CD","#E41A1C","#00CD66","#CD3278","#CD8A96","#00C5CD","#CDCD00","#CD85CD","#CD853F","#8B5A2B","#5CACEE","#EE5C42","#00EE76","#EE4A8C","#EED8AE","#00E5EE","#EEEE00","#EED2EE","#EE9A49","#E41A1C","#377EB8","#FF6A6A","#87CEFA","#6E8B3D","#FFEBCD","#B2DFEE")
cols = cols[1:20]

ggplot(metabolites_category_long, aes(ID, Concentration)) + 
  geom_col(aes(fill = Class))+ 
  facet_grid(Day ~ Group, scales = "free") +
  theme_bw() + 
  labs(y = "Concentration (μmol/g)") + 
  scale_fill_manual(values = cols) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 5), aspect.ratio = 1)

ggsave("metabolite_category_barplot.pdf", width = 12, height = 12)




# 单个代谢物浓度
metabolome = read_tsv("../../../metabolome/metabolome-raw.txt") %>% 
  sjmisc::rotate_df(cn = T, rn = "#SampleID") %>% 
  filter(`#SampleID` %in% meta$`#SampleID`)

top20 = names(sort(colSums(metabolome[,-1]), decreasing = T))[1:20]

top20_metabolome = metabolome %>% select(`#SampleID`, all_of(top20))
others = metabolome %>% select(-all_of(top20))
Others = rowSums(others[, -1])


data = top20_metabolome %>% mutate(Others = Others)

data = meta %>% inner_join(data, by = "#SampleID")

mtbl_lv = c(top20, "Others")
data_long = data %>% 
  pivot_longer(-1:-4, names_to = "metabolites", values_to = "Concentration") %>% 
  mutate(metabolites = fct_relevel(metabolites, rev(mtbl_lv)))

cols = c("#D9D9D9", rev(cols[1:20]))
ggplot(data_long, aes(ID, Concentration)) + 
  geom_col(aes(fill = metabolites))+ 
  facet_grid(Day ~ Group, scales = "free") +
  theme_bw() + 
  labs(y = "Concentration (μmol/g)") + 
  scale_fill_manual(values = cols) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 5), aspect.ratio = 1)
ggsave("metabolite_barplot.pdf", width = 12, height = 12)

group_data = data %>% group_by(Group, Day) %>% 
  summarise(across(where(is.numeric), mean)) %>% 
  pivot_longer(-1:-2, names_to = "metabolites", values_to = "Concentration") %>% 
  mutate(metabolites = fct_relevel(metabolites, rev(mtbl_lv)))
ggplot(group_data, aes(Group, Concentration)) + 
  geom_col(aes(fill = metabolites))+ 
  facet_grid(.~Day, scales = "free") +
  theme_bw() + 
  labs(y = "Concentration (μmol/g)") + 
  scale_fill_manual(values = cols) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5), 
        legend.title = element_blank())
ggsave("supp_fig7_各时间每组代谢物组成柱状图.pdf", width = 10, height = 3)
