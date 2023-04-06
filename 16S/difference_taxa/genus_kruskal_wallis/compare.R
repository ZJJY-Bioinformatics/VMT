library(tidyverse)
library(ggpubr)




map = read_tsv("../../../Metadata_infant_stool_metabolome.txt") %>% 
  select(`#SampleID`, Group, Day) %>% 
  filter(Group == "Con" | Group == "VMT") %>% 
  mutate(Day = fct_relevel(Day, c("D3", "D7", "D30", "D42")), 
         Group = fct_relevel(Group, c("Con", "VMT", "VD"))) 



# ============菌群===============
source("/data/scripts/utility.R")

map = read_tsv("../../../Metadata_infant_stool_microbiome.txt") %>% 
  select(`#SampleID`, Group, Day) %>% 
  filter(Group == "Con" | Group == "VMT") %>% 
  mutate(Day = fct_relevel(Day, c("D3", "D7", "D30", "D42")), 
         Group = fct_relevel(Group, c("Con", "VMT", "VD")))

genus = read_tsv("../Con_VMT/差异菌/infant_stool_e10k_L6_s46.biom.txt", skip = 1) %>% 
  filter(!str_detect(`#OTU ID`, ".*Other$")) %>% 
  filter(!str_detect(`#OTU ID`, ".*g__$"))
genus$`#OTU ID` = taxa_shortname(genus$`#OTU ID`)
genus = genus %>% sjmisc::rotate_df(cn = T, rn = "#SampleID")




# Day3
Day3_kw_res = read_tsv("genus_kruskal_wallis/Day3_kruskal_walls.txt") %>% 
  filter(!str_detect(OTU, ".*Other$")) %>% 
  filter(!str_detect(OTU, ".*g__$")) %>% 
  filter(P < 0.05)
Day3_kw_res$OTU = taxa_shortname(Day3_kw_res$OTU)
Day3_sig = Day3_kw_res$OTU

wilcox = function(microbiota){
  data = genus %>% select(`#SampleID`, all_of(microbiota)) %>% 
    inner_join(map, by = "#SampleID") %>% 
    mutate(Group = fct_relevel(Group, c("Con", "VMT", "VD")))
  res = data %>% filter(Day == "D3") %>% 
    rename("microbiota" = all_of(microbiota)) %>% 
    compare_means(microbiota ~ Group, data = ., p.adjust.method = "fdr") %>% 
    mutate(.y. = microbiota)
  return(res)
}

wilcox_res3 = map_dfr(Day3_sig, ~wilcox(.x)) 
wilcox_res3$Day = "Day3"


# Day7
Day7_kw_res = read_tsv("genus_kruskal_wallis/Day7_kruskal_walls.txt") %>% 
  filter(!str_detect(OTU, ".*Other$")) %>% 
  filter(!str_detect(OTU, ".*g__$")) %>% 
  filter(P < 0.05)
Day7_kw_res$OTU = taxa_shortname(Day7_kw_res$OTU)
Day7_sig = Day7_kw_res$OTU



wilcox = function(microbiota){
  data = genus %>% select(`#SampleID`, all_of(microbiota)) %>% 
    inner_join(map, by = "#SampleID") %>% 
    mutate(Group = fct_relevel(Group, c("Con", "VMT", "VD")))
  res = data %>% filter(Day == "D7") %>% 
    rename("microbiota" = all_of(microbiota)) %>% 
    compare_means(microbiota ~ Group, data = ., p.adjust.method = "fdr") %>% 
    mutate(.y. = microbiota)
  return(res)
}

wilcox_res7 = map_dfr(Day7_sig, ~wilcox(.x)) 
wilcox_res7$Day = "Day7"


# Day30
Day30_kw_res = read_tsv("genus_kruskal_wallis/Day30_kruskal_walls.txt") %>% 
  filter(!str_detect(OTU, ".*Other$")) %>% 
  filter(!str_detect(OTU, ".*g__$")) %>% 
  filter(P < 0.05)
Day30_kw_res$OTU = taxa_shortname(Day30_kw_res$OTU)
Day30_sig = Day30_kw_res$OTU

wilcox = function(microbiota){
  data = genus %>% select(`#SampleID`, all_of(microbiota)) %>% 
    inner_join(map, by = "#SampleID") %>% 
    mutate(Group = fct_relevel(Group, c("Con", "VMT", "VD")))
  res = data %>% filter(Day == "D30") %>% 
    rename("microbiota" = all_of(microbiota)) %>% 
    compare_means(microbiota ~ Group, data = ., p.adjust.method = "fdr") %>% 
    mutate(.y. = microbiota)
  return(res)
}

wilcox_res30 = map_dfr(Day30_sig, ~wilcox(.x)) 
wilcox_res30$Day = "Day30"


# Day42
Day42_kw_res = read_tsv("genus_kruskal_wallis/Day42_kruskal_walls.txt") %>% 
  filter(!str_detect(OTU, ".*Other$")) %>% 
  filter(!str_detect(OTU, ".*g__$")) %>% 
  filter(P < 0.05)
Day42_kw_res$OTU = taxa_shortname(Day42_kw_res$OTU)
Day42_sig = Day42_kw_res$OTU

wilcox = function(microbiota){
  data = genus %>% select(`#SampleID`, all_of(microbiota)) %>% 
    inner_join(map, by = "#SampleID") %>% 
    mutate(Group = fct_relevel(Group, c("Con", "VMT", "VD")))
  res = data %>% filter(Day == "D42") %>% 
    rename("microbiota" = all_of(microbiota)) %>% 
    compare_means(microbiota ~ Group, data = ., p.adjust.method = "fdr") %>% 
    mutate(.y. = microbiota)
  return(res)
}

wilcox_res42 = map_dfr(Day42_sig, ~wilcox(.x)) 
wilcox_res42$Day = "Day42"

total_res = bind_rows(wilcox_res3,wilcox_res7, wilcox_res30, wilcox_res42)
write_tsv(total_res, "菌在各时间点wilcox检验.tsv")



data = genus %>% 
  inner_join(map, by = "#SampleID")

need = colnames(genus)[-1]


plot_func = function(metabolite){
  ggline(data, x = "Day", y = metabolite, add = "mean_se", color = "Group", palette = c('#1f77b4', '#ff7f0e', '#2ca02c')) + 
    stat_compare_means(aes(group=Group), label = "p.signif", method = "wilcox.test", label.y.npc = 0.05)
}

plot_func(need[2])
ps = map(need, plot_func)

pdf("各菌两组间errorbar图.pdf", width = 4, height = 3)
walk(ps, print)
dev.off()


boxplot_func = function(metabolite){
  ggboxplot(data, x = "Day", y = metabolite, fill = "Group", palette = c('#1f77b4', '#ff7f0e', '#2ca02c')) + 
    stat_compare_means(aes(group=Group), 
                       label = "p", 
                       method = "wilcox.test")
}

boxplot_func(need[2])
ps = map(need, boxplot_func)

pdf("各菌两组间boxplot_P.pdf", width = 6, height = 5)
walk(ps, print)
dev.off()


##bifido
library(introdataviz)
ggplot(data, aes(Day, asin(sqrt(Bifidobacterium)), fill = Group)) + 
  introdataviz::geom_split_violin(alpha = 0.6, scale = "width", trim = T) + 
  #geom_boxplot(width = 0.2, alpha = 0.5) + 
  scale_fill_manual(values = c("#7BC2CE", "#8562CC")) + 
  #scale_fill_manual(values = c('#1f77b4', '#ff7f0e')) + 
  stat_summary(fun.data = "median_q1q3", geom = "pointrange", show.legend = F, position = position_dodge(0.35), size = 0.5) + 
  theme_minimal() +
  stat_compare_means(aes(group=Group),  method = "wilcox.test") + 
  labs(y = "Asin(sqrt(Bifidobacterium))")

ggsave("Bifido_color2.pdf", width = 5, height = 3.5)

ggplot(data, aes(Day, asin(sqrt(Lactobacillus)), fill = Group)) + 
  introdataviz::geom_split_violin(alpha = 0.6, scale = "width", trim = T) + 
  #geom_boxplot(width = 0.2, alpha = 0.5) + 
  scale_fill_manual(values = c("#7BC2CE", "#8562CC")) + 
  #scale_fill_manual(values = c('#1f77b4', '#ff7f0e')) + 
  stat_summary(fun.data = "median_q1q3", geom = "pointrange", show.legend = F, position = position_dodge(0.35), size = 0.5) + 
  theme_minimal() +
  stat_compare_means(aes(group=Group), label = "p.signif", method = "wilcox.test") + 
  labs(y = "Asin(sqrt(Lactobacillus))") +
  ylim(quantile(asin(sqrt(data$Lactobacillus)), c(0, 0.9)))


ggsave("Lacto_color2.pdf", width = 5, height = 3.5)





