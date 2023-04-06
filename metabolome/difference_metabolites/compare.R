library(tidyverse)
library(ggpubr)




map = read_tsv("../../Metadata_infant_stool_metabolome.txt") %>% 
  select(`#SampleID`, Group, Day) %>% 
  filter(Group == "Con" | Group == "VMT") %>% 
  filter(Day != "D1" & Day != "D2") %>% 
  mutate(Day = fct_relevel(Day, c("D1", "D2", "D3", "D7", "D30", "D42")), 
         Group = fct_relevel(Group, c("Con", "VMT", "VD"))) 
metabolites = read_tsv("../../../metabolome/metabolome-tsfm.txt") %>% filter(Sample %in% map$`#SampleID`)

mtblt_abs_prs = metabolites %>% column_to_rownames(var = "Sample") %>% 
  mutate(across(everything(), ~if_else(.x > 0, 1, 0))) %>% 
  colSums()

idx_10pc =  mtblt_abs_prs > 44


metabolites_s56 = metabolites %>% column_to_rownames(var = "Sample") %>% 
  .[, idx_10pc]

# 在少于10%样本中出现的代谢物
metabolites_none_s44 = metabolites %>% column_to_rownames(var = "Sample") %>% 
  .[, !idx_10pc] %>% colnames()






# Day1
Day1_kw_res = read_tsv("day1_Group_kruskal.txt") %>% 
  filter(P < 0.05)
Day1_sig = Day1_kw_res$OTU

wilcox = function(metabolite){
  data = metabolites %>% select(Sample, all_of(metabolite)) %>% 
    inner_join(map, by = c("Sample" = "#SampleID")) %>% 
    mutate(Group = fct_relevel(Group, c("Con", "VMT", "VD")))
  res = data %>% filter(Day == "D1") %>% 
    rename("metabolite" = all_of(metabolite)) %>% 
    compare_means(metabolite ~ Group, data = ., p.adjust.method = "fdr") %>% 
    mutate(.y. = metabolite)
  return(res)
}


wilcox_res1 = map_dfr(Day1_sig, ~wilcox(.x)) 
wilcox_res1$Day = "Day1"

# Day2
Day2_kw_res = read_tsv("day2_Group_kruskal.txt") %>% 
  filter(P < 0.05)
Day2_sig = Day2_kw_res$OTU

wilcox = function(metabolite){
  data = metabolites %>% select(Sample, all_of(metabolite)) %>% 
    inner_join(map, by = c("Sample" = "#SampleID")) %>% 
    mutate(Group = fct_relevel(Group, c("Con", "VMT", "VD")))
  res = data %>% filter(Day == "D2") %>% 
    rename("metabolite" = all_of(metabolite)) %>% 
    compare_means(metabolite ~ Group, data = ., p.adjust.method = "fdr") %>% 
    mutate(.y. = metabolite)
  return(res)
}


wilcox_res2 = map_dfr(Day2_sig, ~wilcox(.x)) 
wilcox_res2$Day = "Day2"


# Day3
Day3_kw_res = read_tsv("day3_Group_kruskal.txt") %>% 
  filter(P < 0.05)
Day3_sig = Day3_kw_res$OTU

wilcox = function(metabolite){
  data = metabolites %>% select(Sample, all_of(metabolite)) %>% 
    inner_join(map, by = c("Sample" = "#SampleID")) %>% 
    mutate(Group = fct_relevel(Group, c("Con", "VMT", "VD")))
  res = data %>% filter(Day == "D3") %>% 
    rename("metabolite" = all_of(metabolite)) %>% 
    compare_means(metabolite ~ Group, data = ., p.adjust.method = "fdr") %>% 
    mutate(.y. = metabolite)
  return(res)
}


wilcox_res3 = map_dfr(Day3_sig, ~wilcox(.x)) 
wilcox_res3$Day = "Day3"

# Day7
Day7_kw_res = read_tsv("day7_Group_kruskal.txt") %>% 
  filter(P < 0.05)
Day7_sig = Day7_kw_res$OTU

wilcox = function(metabolite){
  data = metabolites %>% select(Sample, all_of(metabolite)) %>% 
    inner_join(map, by = c("Sample" = "#SampleID")) %>% 
    mutate(Group = fct_relevel(Group, c("Con", "VMT", "VD")))
  res = data %>% filter(Day == "D7") %>% 
    rename("metabolite" = all_of(metabolite)) %>% 
    compare_means(metabolite ~ Group, data = ., p.adjust.method = "fdr") %>% 
    mutate(.y. = metabolite)
  return(res)
}


wilcox_res7 = map_dfr(Day7_sig, ~wilcox(.x)) 
wilcox_res7$Day = "Day7"

# Day30
Day30_kw_res = read_tsv("day30_Group_kruskal.txt") %>% 
  filter(P < 0.05)
Day30_sig = Day30_kw_res$OTU

wilcox = function(metabolite){
  data = metabolites %>% select(Sample, all_of(metabolite)) %>% 
    inner_join(map, by = c("Sample" = "#SampleID")) %>% 
    mutate(Group = fct_relevel(Group, c("Con", "VMT", "VD")))
  res = data %>% filter(Day == "D30") %>% 
    rename("metabolite" = all_of(metabolite)) %>% 
    compare_means(metabolite ~ Group, data = ., p.adjust.method = "fdr") %>% 
    mutate(.y. = metabolite)
  return(res)
}


wilcox_res30 = map_dfr(Day30_sig, ~wilcox(.x)) 
wilcox_res30$Day = "Day30"


# Day42
Day42_kw_res = read_tsv("day42_Group_kruskal.txt") %>% 
  filter(P < 0.05)
Day42_sig = Day42_kw_res$OTU

wilcox = function(metabolite){
  data = metabolites %>% select(Sample, all_of(metabolite)) %>% 
    inner_join(map, by = c("Sample" = "#SampleID")) %>% 
    mutate(Group = fct_relevel(Group, c("Con", "VMT", "VD")))
  res = data %>% filter(Day == "D42") %>% 
    rename("metabolite" = all_of(metabolite)) %>% 
    compare_means(metabolite ~ Group, data = ., p.adjust.method = "fdr") %>% 
    mutate(.y. = metabolite)
  return(res)
}


wilcox_res42 = map_dfr(Day42_sig, ~wilcox(.x)) 
wilcox_res42$Day = "Day42"


total_res = bind_rows(wilcox_res1, wilcox_res2, wilcox_res3,wilcox_res7, wilcox_res30, wilcox_res42)
write_tsv(total_res, "各代谢物在各时间点wilcox检验.tsv")


data = metabolites_s56 %>% rownames_to_column(var = "Sample") %>% 
  inner_join(map, by = c("Sample" = "#SampleID"))

need = colnames(metabolites_s56)

plot_func = function(metabolite){
  ggline(data, x = "Day", y = metabolite, add = "mean_se", color = "Group", palette = c('#1f77b4', '#ff7f0e', '#2ca02c')) + 
    stat_compare_means(aes(group=Group), label = "p.signif", method = "kruskal.test", label.y.npc = 0.2)
}


plot_func(need[2])
ps = map(need, plot_func)

pdf("各代谢物errorbar图.pdf", width = 4, height = 3)
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

pdf("各代谢物boxplot.pdf", width = 4, height = 3)
walk(ps, print)
dev.off()

# ============菌群===============
source("/data/scripts/utility.R")

map = read_tsv("../../Metadata_infant_stool_microbiome.txt") %>% 
  select(`#SampleID`, Group, Day) %>% 
  filter(Day != "D1" & Day != "D2") %>% 
  filter(Group == "Con" | Group == "VMT") %>% 
  mutate(Day = fct_relevel(Day, c("D1", "D2", "D3", "D7", "D30", "D42")), 
         Group = fct_relevel(Group, c("Con", "VMT", "VD")))

genus = read_tsv("infant_stool_e30k_L6_s46.biom.txt", skip = 1) %>% 
  filter(!str_detect(`#OTU ID`, ".*Other$")) %>% 
  filter(!str_detect(`#OTU ID`, ".*g__$"))
genus$`#OTU ID` = taxa_shortname(genus$`#OTU ID`)
genus = genus %>% sjmisc::rotate_df(cn = T, rn = "#SampleID")


# Day1
Day1_kw_res = read_tsv("genus_kruskal_wallis/Day1_kruskal_wallis.txt") %>% 
  filter(!str_detect(OTU, ".*Other$")) %>% 
  filter(!str_detect(OTU, ".*g__$")) %>% 
  filter(P < 0.05)
Day1_kw_res$OTU = taxa_shortname(Day1_kw_res$OTU)
Day1_sig = Day1_kw_res$OTU

wilcox = function(microbiota){
  data = genus %>% select(`#SampleID`, all_of(microbiota)) %>% 
    inner_join(map, by = "#SampleID") %>% 
    mutate(Group = fct_relevel(Group, c("Con", "VMT", "VD")))
  res = data %>% filter(Day == "D1") %>% 
    rename("microbiota" = all_of(microbiota)) %>% 
    compare_means(microbiota ~ Group, data = ., p.adjust.method = "fdr") %>% 
    mutate(.y. = microbiota)
  return(res)
}

wilcox_res1 = map_dfr(Day1_sig, ~wilcox(.x)) 
wilcox_res1$Day = "Day1"




# Day2
Day2_kw_res = read_tsv("genus_kruskal_wallis/Day2_kruskal_wallis.txt") %>% 
  filter(!str_detect(OTU, ".*Other$")) %>% 
  filter(!str_detect(OTU, ".*g__$")) %>% 
  filter(P < 0.05)
Day2_kw_res$OTU = taxa_shortname(Day2_kw_res$OTU)
Day2_sig = Day2_kw_res$OTU

wilcox = function(microbiota){
  data = genus %>% select(`#SampleID`, all_of(microbiota)) %>% 
    inner_join(map, by = "#SampleID") %>% 
    mutate(Group = fct_relevel(Group, c("Con", "VMT", "VD")))
  res = data %>% filter(Day == "D2") %>% 
    rename("microbiota" = all_of(microbiota)) %>% 
    compare_means(microbiota ~ Group, data = ., p.adjust.method = "fdr") %>% 
    mutate(.y. = microbiota)
  return(res)
}

wilcox_res2 = map_dfr(Day2_sig, ~wilcox(.x)) 
wilcox_res2$Day = "Day2"


# Day3
Day3_kw_res = read_tsv("genus_kruskal_wallis/Day3_kruskal_wallis.txt") %>% 
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
Day7_kw_res = read_tsv("genus_kruskal_wallis/Day7_kruskal_wallis.txt") %>% 
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
Day30_kw_res = read_tsv("genus_kruskal_wallis/Day30_kruskal_wallis.txt") %>% 
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
Day42_kw_res = read_tsv("genus_kruskal_wallis/Day42_kruskal_wallis.txt") %>% 
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

total_res = bind_rows(wilcox_res1, wilcox_res2, wilcox_res3,wilcox_res7, wilcox_res30, wilcox_res42)
write_tsv(total_res, "菌在各时间点wilcox检验.tsv")



data = genus %>% 
  inner_join(map, by = "#SampleID")

need = colnames(genus)[-1]


plot_func = function(metabolite){
  ggline(data, x = "Day", y = metabolite, add = "mean_se", color = "Group", palette = c('#1f77b4', '#ff7f0e', '#2ca02c')) + 
    stat_compare_means(aes(group=Group), label = "p", method = "wilcox.test", label.y.npc = 0.05)
}

plot_func(need[2])
ps = map(need, plot_func)

pdf("各菌errorbar图.pdf", width = 4, height = 3)
walk(ps, print)
dev.off()

boxplot_func(need[2])
ps = map(need, boxplot_func)

pdf("各菌boxplot——kruskalP.pdf", width = 6, height = 5)
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
  stat_compare_means(aes(group=Group), label = "p.signif", method = "wilcox.test") + 
  #scale_y_continuous(limits = quantile(data$Bifidobacterium, c(0.1, 0.9))) +
  labs(y = "Asin(sqrt(Bifidobacterium))")

ggsave("Bifido_color2.pdf", width = 4, height = 3.5)

ggplot(data, aes(Day, asin(sqrt(Lactobacillus)), fill = Group)) + 
  introdataviz::geom_split_violin(alpha = 0.6, scale = "width", trim = T) + 
  #geom_boxplot(width = 0.2, alpha = 0.5) + 
  scale_fill_manual(values = c("#7BC2CE", "#8562CC")) + 
  #scale_fill_manual(values = c('#1f77b4', '#ff7f0e')) + 
  stat_summary(fun.data = "median_q1q3", geom = "pointrange", show.legend = F, position = position_dodge(0.35), size = 0.5) + 
  theme_minimal() +
  stat_compare_means(aes(group=Group), label = "p.signif", method = "wilcox.test") + 
  labs(y = "Asin(sqrt(Lactobacillus))")

ggsave("Lacto_color2.pdf", width = 5, height = 3.5)





