library(tidyverse)


Day3_p = pvalue("Day3")
Day7_p = pvalue("Day7")
Day30_p = pvalue("Day30")
Day42_p = pvalue("Day42")
padjvalue = function(day){
  pvalue = mtblt_res %>% filter(Day == day) %>% 
    mutate(compare = str_c(group1, group2, "p.adj", sep = "_")) %>% 
    select(1, compare, p.adj) %>% 
    pivot_wider(everything(), names_from = compare, values_from = p.adj)
  
}

Day3_padj = padjvalue("Day3")
Day3_padj$Day = "Day3"
Day7_padj = padjvalue("Day7")
Day7_padj$Day = "Day7"
Day30_padj = padjvalue("Day30")
Day30_padj$Day = "Day30"
Day42_padj = padjvalue("Day42")
Day42_padj$Day = "Day42"

d3 = Day3_p %>% inner_join(Day3_padj, by = ".y.")
d7 = Day7_p %>% inner_join(Day7_padj, by = ".y.")
d30 = Day30_p %>% inner_join(Day30_padj, by = ".y.")
d42 = Day42_p %>% inner_join(Day42_padj, by = ".y.")





median_value = map %>% inner_join(metabolites, by = c("#SampleID" = "Sample")) %>% 
  group_by(Group, Day) %>% 
  summarise(across(where(is.numeric), median, na.rm = T))

d3_kruskal = read_tsv("day3_Group_kruskal.txt") %>% select(OTU, P, VMT_mean, VD_mean, Con_mean)  %>% rename(kruskal_P = P, VMT_mean = VMT_mean, Con_mean = Con_mean)
d7_kruskal = read_tsv("day7_Group_kruskal.txt") %>% select(OTU, P, VMT_mean, VD_mean, Con_mean)  %>% rename(kruskal_P = P, VMT_mean = VMT_mean, Con_mean = Con_mean)
d30_kruskal = read_tsv("day30_Group_kruskal.txt") %>% select(OTU, P, VMT_mean, VD_mean, Con_mean)  %>% rename(kruskal_P = P, VMT_mean = VMT_mean, Con_mean = Con_mean)
d42_kruskal = read_tsv("day42_Group_kruskal.txt") %>% select(OTU, P, VMT_mean, VD_mean, Con_mean)  %>% rename(kruskal_P = P, VMT_mean = VMT_mean, Con_mean = Con_mean)




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


out = bind_rows(D3_median, D7_median, D30_median, D42_median)
write_tsv(out, "各代谢物在各时间点wilcox检验3.tsv")


# ===================菌群====================
source("/data/scripts/utility.R")
mcb_res = read_tsv("菌在各时间点wilcox检验.tsv")
genus = read_tsv("../Con_VMT/差异菌/infant_stool_e10k_L6_s46.biom.txt", skip = 1) %>% 
  filter(!str_detect(`#OTU ID`, ".*Other$")) %>% 
  filter(!str_detect(`#OTU ID`, ".*g__$"))
genus$`#OTU ID` = taxa_shortname(genus$`#OTU ID`)
genus = genus %>% sjmisc::rotate_df(cn = T, rn = "#SampleID")

map = read_tsv("../../../Metadata_infant_stool_microbiome.txt") %>% 
  select(`#SampleID`, Group, Day) %>% 
  mutate(Day = fct_relevel(Day, c("D3", "D7", "D30", "D42")), 
         Group = fct_relevel(Group, c("Con", "VMT", "VD")))

pvalue = function(day){
  pvalue = mcb_res %>% filter(Day == day) %>% 
    mutate(compare = str_c(group1, group2, "p", sep = "_")) %>% 
    select(1, compare, p) %>% 
    pivot_wider(everything(), names_from = compare, values_from = p)
  
}




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

Day3_padj = padjvalue("Day3")
Day3_padj$Day = "Day3"
Day7_padj = padjvalue("Day7")
Day7_padj$Day = "Day7"
Day30_padj = padjvalue("Day30")
Day30_padj$Day = "Day30"
Day42_padj = padjvalue("Day42")
Day42_padj$Day = "Day42"


d3 = Day3_p %>% inner_join(Day3_padj, by = ".y.")
d7 = Day7_p %>% inner_join(Day7_padj, by = ".y.")
d30 = Day30_p %>% inner_join(Day30_padj, by = ".y.")
d42 = Day42_p %>% inner_join(Day42_padj, by = ".y.")

median_value = map %>% inner_join(genus, by = "#SampleID") %>% 
  select(-`#SampleID`) %>% 
  group_by(Group, Day) %>% 
  summarise(across(where(is.numeric), median, na.rm = T))



d3_kruskal = read_tsv("genus_kruskal_wallis/Day3_kruskal_walls.txt") %>% 
  filter(!str_detect(OTU, ".*Other$")) %>% 
  filter(!str_detect(OTU, ".*g__$")) %>% 
  select(OTU, P, VMT_mean, VD_mean, Con_mean) %>% 
  rename(kruskal_P = P) 
d3_kruskal$OTU = taxa_shortname(d3_kruskal$OTU)

d7_kruskal = read_tsv("genus_kruskal_wallis/Day7_kruskal_walls.txt")  %>% 
  filter(!str_detect(OTU, ".*Other$")) %>% 
  filter(!str_detect(OTU, ".*g__$")) %>% 
  select(OTU, P, VMT_mean, VD_mean, Con_mean) %>% 
  rename(kruskal_P = P) 
d7_kruskal$OTU = taxa_shortname(d7_kruskal$OTU)

d30_kruskal = read_tsv("genus_kruskal_wallis/Day30_kruskal_walls.txt") %>% 
  filter(!str_detect(OTU, ".*Other$")) %>% 
  filter(!str_detect(OTU, ".*g__$")) %>% 
  select(OTU, P, VMT_mean, VD_mean, Con_mean) %>% 
  rename(kruskal_P = P) 
d30_kruskal$OTU = taxa_shortname(d30_kruskal$OTU)
d42_kruskal = read_tsv("genus_kruskal_wallis/Day42_kruskal_walls.txt")  %>% 
  filter(!str_detect(OTU, ".*Other$")) %>% 
  filter(!str_detect(OTU, ".*g__$")) %>% 
  select(OTU, P, VMT_mean, VD_mean, Con_mean) %>% 
  rename(kruskal_P = P) 
d42_kruskal$OTU = taxa_shortname(d42_kruskal$OTU)


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
out = bind_rows(D3_median, D7_median, D30_median, D42_median)
write_tsv(out, "各菌在各时间点wilcox检验3.tsv")

