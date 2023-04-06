library(tidyverse)
library(ggpubr)
source("/data/scripts/utility.R")
# ============菌群===============


map = read_tsv("../../../../Metadata_infant_stool_microbiome.txt") %>% 
  select(`#SampleID`, Group, Day) %>% 
  filter(Group == "VMT" | Group == "Con") %>% 
  rename(Day = Day) %>% 
  mutate(Day = fct_relevel(Day, c("D3", "D7", "D30", "D42")), 
         Group = fct_relevel(Group, c("Con", "VMT", "VD")))

family = read_tsv("../../Con_VMT/差异菌/infant_stool_e10k_L5_s46.biom.txt", skip = 1) %>% 
  filter(!str_detect(`#OTU ID`, ".*Other$")) %>% 
  filter(!str_detect(`#OTU ID`, ".*f__$"))
family$`#OTU ID` = taxa_shortname(family$`#OTU ID`)
family = family %>% sjmisc::rotate_df(cn = T, rn = "#SampleID")




wilcox = function(microbiota){
  data = family %>% select(`#SampleID`, all_of(microbiota)) %>% 
    inner_join(map, by = "#SampleID") %>% 
    mutate(Group = fct_relevel(Group, c("Con", "VMT", "VD")))
  res = data %>% filter(Day == "D1") %>% 
    rename("microbiota" = all_of(microbiota)) %>% 
    compare_means(microbiota ~ Group, data = ., p.adjust.method = "fdr") %>% 
    mutate(.y. = microbiota)
  return(res)
}



# Day3
Day3_kw_res = read_tsv("Day3_kruskal_walls.txt") %>% 
  filter(!str_detect(OTU, ".*Other$")) %>% 
  filter(!str_detect(OTU, ".*f__$")) %>% 
  filter(P < 0.05)
Day3_kw_res$OTU = taxa_shortname(Day3_kw_res$OTU)
Day3_sig = Day3_kw_res$OTU

wilcox = function(microbiota){
  data = family %>% select(`#SampleID`, all_of(microbiota)) %>% 
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
Day7_kw_res = read_tsv("Day7_kruskal_walls.txt") %>% 
  filter(!str_detect(OTU, ".*Other$")) %>% 
  filter(!str_detect(OTU, ".*f__$")) %>% 
  filter(P < 0.05)
Day7_kw_res$OTU = taxa_shortname(Day7_kw_res$OTU)
Day7_sig = Day7_kw_res$OTU

wilcox = function(microbiota){
  data = family %>% select(`#SampleID`, all_of(microbiota)) %>% 
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
Day30_kw_res = read_tsv("Day30_kruskal_walls.txt") %>% 
  filter(!str_detect(OTU, ".*Other$")) %>% 
  filter(!str_detect(OTU, ".*f__$")) %>% 
  filter(P < 0.05)
Day30_kw_res$OTU = taxa_shortname(Day30_kw_res$OTU)
Day30_sig = Day30_kw_res$OTU

wilcox = function(microbiota){
  data = family %>% select(`#SampleID`, all_of(microbiota)) %>% 
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
Day42_kw_res = read_tsv("Day42_kruskal_walls.txt") %>% 
  filter(!str_detect(OTU, ".*Other$")) %>% 
  filter(!str_detect(OTU, ".*f__$")) %>% 
  filter(P < 0.05)
Day42_kw_res$OTU = taxa_shortname(Day42_kw_res$OTU)
Day42_sig = Day42_kw_res$OTU

wilcox = function(microbiota){
  data = family %>% select(`#SampleID`, all_of(microbiota)) %>% 
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
write_tsv(total_res, "family在各时间点wilcox检验.tsv")



data = family %>% 
  inner_join(map, by = "#SampleID")

need = colnames(family)[-1]


plot_func = function(metabolite){
  ggline(data, x = "Day", y = metabolite, add = "mean_se", color = "Group", palette = c('#1f77b4', '#ff7f0e', '#2ca02c')) + 
    stat_compare_means(aes(group=Group), label = "p.signif", method = "wilcox.test", label.y.npc = 0.05)
}

plot_func(need[2])
ps = map(need, plot_func)

pdf("两组各菌errorbar图.pdf", width = 4, height = 3)
walk(ps, print)
dev.off()


boxplot_func = function(metabolite){
  ggboxplot(data, x = "Day", y = metabolite, fill = "Group", palette = c('#1f77b4', '#ff7f0e', '#2ca02c')) + 
    stat_compare_means(aes(group=Group), 
                       label = "p.signif", 
                       method = "wilcox.test")
}
boxplot_func(need[2])
ps = map(need, boxplot_func)

pdf("两组各菌boxplot.pdf", width = 4, height = 3)
walk(ps, print)
dev.off()










