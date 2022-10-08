library(tidyverse)
source("/data/scripts/utility.R")

# ===================菌群====================

mcb_res = read_tsv("family在各时间点wilcox检验.tsv")
family = read_tsv("../../Con_VMT/差异菌/infant_stool_e10k_L5_s46.biom.txt", skip = 1) %>% 
  filter(!str_detect(`#OTU ID`, ".*Other$")) %>% 
  filter(!str_detect(`#OTU ID`, ".*f__$"))
family$`#OTU ID` = taxa_shortname(family$`#OTU ID`)
family = family %>% sjmisc::rotate_df(cn = T, rn = "#SampleID")

map = read_tsv("../../../../Metadata_infant_stool_microbiome.txt") %>% 
  select(`#SampleID`, Group, Day) %>% 
  mutate(Day = fct_relevel(Day, c("D1", "D2", "D3", "D7", "D30", "D42")), 
         Group = fct_relevel(Group, c("Con", "VMT", "VD")))

pvalue = function(day){
  pvalue = mcb_res %>% filter(Day == day) %>% 
    mutate(compare = str_c(group1, group2, "p", sep = "_")) %>% 
    select(1, compare, p) %>% 
    pivot_wider(everything(), names_from = compare, values_from = p)
  
}



Day1_p = pvalue("Day1")

Day2_p = pvalue("Day2")
Day3_p = pvalue("Day3")
Day7_p = pvalue("Day7")
Day30_p = pvalue("Day30")
Day42_p = pvalue("Day42")
padjvalue = function(day){
  pvalue = mcb_res %>% filter(Day == day) %>% 
    mutate(compare = str_c(group1, group2, "p.adj", sep = "_")) %>% 
    select(1, compare, p.adj) %>% 
    pivot_wider(everything(), names_from = compare, values_from = p.adj)
  
}
Day1_padj = padjvalue("Day1")
Day1_padj$Day = "Day1"
Day2_padj = padjvalue("Day2")
Day2_padj$Day = "Day2"
Day3_padj = padjvalue("Day3")
Day3_padj$Day = "Day3"
Day7_padj = padjvalue("Day7")
Day7_padj$Day = "Day7"
Day30_padj = padjvalue("Day30")
Day30_padj$Day = "Day30"
Day42_padj = padjvalue("Day42")
Day42_padj$Day = "Day42"

d1 = Day1_p %>% inner_join(Day1_padj, by = ".y.")
d2 = Day2_p %>% inner_join(Day2_padj, by = ".y.")
d3 = Day3_p %>% inner_join(Day3_padj, by = ".y.")
d7 = Day7_p %>% inner_join(Day7_padj, by = ".y.")
d30 = Day30_p %>% inner_join(Day30_padj, by = ".y.")
d42 = Day42_p %>% inner_join(Day42_padj, by = ".y.")

median_value = map %>% inner_join(family, by = "#SampleID") %>% 
  select(-`#SampleID`) %>% 
  group_by(Group, Day) %>% 
  summarise(across(where(is.numeric), median, na.rm = T))


d1_kruskal = read_tsv("Day1_kruskal_walls.txt") %>% 
  filter(!str_detect(OTU, ".*Other$")) %>% 
  filter(!str_detect(OTU, ".*f__$")) %>% 
  select(OTU, P, VMT_mean, VD_mean, Con_mean) %>% 
  rename(kruskal_P = P) 
d1_kruskal$OTU = taxa_shortname(d1_kruskal$OTU)
d2_kruskal = read_tsv("Day2_kruskal_walls.txt") %>% 
  filter(!str_detect(OTU, ".*Other$")) %>% 
  filter(!str_detect(OTU, ".*f__$")) %>% 
  select(OTU, P, VMT_mean, VD_mean, Con_mean) %>% 
  rename(kruskal_P = P) 
d2_kruskal$OTU = taxa_shortname(d2_kruskal$OTU)

d3_kruskal = read_tsv("Day3_kruskal_walls.txt") %>% 
  filter(!str_detect(OTU, ".*Other$")) %>% 
  filter(!str_detect(OTU, ".*f__$")) %>% 
  select(OTU, P, VMT_mean, VD_mean, Con_mean) %>% 
  rename(kruskal_P = P) 
d3_kruskal$OTU = taxa_shortname(d3_kruskal$OTU)

d7_kruskal = read_tsv("Day7_kruskal_walls.txt")  %>% 
  filter(!str_detect(OTU, ".*Other$")) %>% 
  filter(!str_detect(OTU, ".*f__$")) %>% 
  select(OTU, P, VMT_mean, VD_mean, Con_mean) %>% 
  rename(kruskal_P = P) 
d7_kruskal$OTU = taxa_shortname(d7_kruskal$OTU)

d30_kruskal = read_tsv("Day30_kruskal_walls.txt") %>% 
  filter(!str_detect(OTU, ".*Other$")) %>% 
  filter(!str_detect(OTU, ".*f__$")) %>% 
  select(OTU, P, VMT_mean, VD_mean, Con_mean) %>% 
  rename(kruskal_P = P) 
d30_kruskal$OTU = taxa_shortname(d30_kruskal$OTU)
d42_kruskal = read_tsv("Day42_kruskal_walls.txt")  %>% 
  filter(!str_detect(OTU, ".*Other$")) %>% 
  filter(!str_detect(OTU, ".*f__$")) %>% 
  select(OTU, P, VMT_mean, VD_mean, Con_mean) %>% 
  rename(kruskal_P = P) 
d42_kruskal$OTU = taxa_shortname(d42_kruskal$OTU)

D1_median = median_value %>% filter(Day == "D1") %>% 
  select(-Day) %>% 
  ungroup() %>% 
  sjmisc::rotate_df(cn = T, rn = ".y.") %>% 
  rename(Con_median = Con, VMT_median = VMT, VD_median = VD) %>% 
  inner_join(d1, ., by = ".y.") %>% 
  inner_join(d1_kruskal, by = c(".y." = "OTU"))

D2_median = median_value %>% filter(Day == "D2") %>% 
  select(-Day) %>% 
  ungroup() %>% 
  sjmisc::rotate_df(cn = T, rn = ".y.") %>% 
  rename(Con_median = Con, VMT_median = VMT, VD_median = VD) %>% 
  inner_join(d2, ., by = ".y.")%>% 
  inner_join(d2_kruskal, by = c(".y." = "OTU"))

D3_median = median_value %>% filter(Day == "D3") %>% 
  select(-Day) %>% 
  ungroup() %>% 
  sjmisc::rotate_df(cn = T, rn = ".y.") %>% 
  rename(Con_median = Con, VMT_median = VMT, VD_median = VD) %>% 
  inner_join(d3, ., by = ".y.") %>% 
  inner_join(d3_kruskal, by = c(".y." = "OTU"))

D7_median = median_value %>% filter(Day == "D7") %>% 
  select(-Day) %>% 
  ungroup() %>% 
  sjmisc::rotate_df(cn = T, rn = ".y.") %>% 
  rename(Con_median = Con, VMT_median = VMT, VD_median = VD) %>% 
  inner_join(d7, ., by = ".y.") %>% 
  inner_join(d7_kruskal, by = c(".y." = "OTU"))

D30_median = median_value %>% filter(Day == "D30") %>% 
  select(-Day) %>% 
  ungroup() %>% 
  sjmisc::rotate_df(cn = T, rn = ".y.") %>% 
  rename(Con_median = Con, VMT_median = VMT, VD_median = VD) %>% 
  inner_join(d30, ., by = ".y.") %>% 
  inner_join(d30_kruskal, by = c(".y." = "OTU"))
D42_median = median_value %>% filter(Day == "D42") %>% 
  select(-Day) %>% 
  ungroup() %>% 
  sjmisc::rotate_df(cn = T, rn = ".y.") %>% 
  rename(Con_median = Con, VMT_median = VMT, VD_median = VD) %>% 
  inner_join(d42, ., by = ".y.") %>% 
  inner_join(d42_kruskal, by = c(".y." = "OTU"))
out = bind_rows(D1_median, D2_median, D3_median, D7_median, D30_median, D42_median)
write_tsv(out, "各family在各时间点wilcox检验2.tsv")

