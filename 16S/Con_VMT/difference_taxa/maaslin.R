library(tidyverse)
library(Maaslin2)
source("/data/scripts/utility.R")

map_file = "../../VMT_Con_metadata_infant_stool_microbiome.txt"
other_confounder = readxl::read_xlsx("confounders.xlsx")

L6 = read_tsv("../infant_stool_e10k_L6_s46.biom.txt", skip = 1) %>% 
  filter(!str_detect(`#OTU ID`, ".*Other")) %>% 
  filter(!str_detect(`#OTU ID`, ".*g__$"))
L6$`#OTU ID` = taxa_shortname(L6$`#OTU ID`)
L6 = L6 %>% sjmisc::rotate_df(cn = T)

L5 = read_tsv("../infant_stool_e10k_L5_s46.biom.txt", skip = 1) %>% 
  filter(!str_detect(`#OTU ID`, ".*Other")) %>% 
  filter(!str_detect(`#OTU ID`, ".*f__$"))
L5$`#OTU ID` = taxa_shortname(L5$`#OTU ID`)
L5 = L5 %>% sjmisc::rotate_df(cn = T)





time = "D3"


map = read_tsv(map_file) %>% 
  filter(Group == "VMT" | Group == "Con") %>% 
  left_join(other_confounder, by = c("ID" = "PatientID")) %>% 
  filter(Day == time) %>% 
  column_to_rownames(var = "#SampleID") %>% 
  mutate(Gender = fct_relevel(Gender, c("Boy", "Girl")),
         Group = fct_relevel(Group, c("Con", "VMT")))
sample = intersect(rownames(map), rownames(L6))

map = map[sample,]
L6_subset = L6[sample,]

Maaslin2(L6_subset, map, output = paste0(time, "_L6"), 
         fixed_effects = c("Group", "Gestational_weeks", "Birth_weight", "Gender", "Feeding_3"),
         #fixed_effects = c("Group"),
         reference = c("Feeding_3,BF"),
         normalization = "TSS", 
         transform = "LOG")


# ========L5=============

time = "D42"
map = read_tsv(map_file) %>% 
  filter(Group == "VMT" | Group == "Con") %>% 
  left_join(other_confounder, by = c("ID" = "PatientID")) %>% 
  filter(Day == time) %>% 
  column_to_rownames(var = "#SampleID") %>% 
  mutate(Gender = fct_relevel(Gender, c("Boy", "Girl")),
         Group = fct_relevel(Group, c("Con", "VMT")))
sample = intersect(rownames(map), rownames(L5))

map = map[sample,]
L5_subset = L5[sample,]

Maaslin2(L5_subset, map, output = paste0(time, "_L5"), 
         fixed_effects = c("Group", "Gestational_weeks", "Birth_weight", "Gender", "Feeding_42"),
         #fixed_effects = c("Group"),
         reference = c("Feeding_42,BF"),
         normalization = "TSS", 
         transform = "LOG")




# 提取六月龄总分的结果
time = "3"
"D3_L6/"
func1 = function(time){
  read_tsv(paste0("D", time, "_L6/all_results.tsv")) %>%
    filter(metadata == "Group") %>% 
    mutate(Day = paste0("Day", {{time}}))
}

total = map_dfr(c(3,7,30,42), ~func1(.x))
write_tsv(total, "Group_with_L6_Maaslin_multivariate.tsv")



# 重新校正
read_tsv("Group_with_L6_Maaslin_multivariate.tsv") %>% 
  group_by(Day) %>% 
  mutate(p_adjust = p.adjust(pval, method = "BH")) %>% 
  write_tsv("Group_with_L6_Maaslin_multivariate.tsv")
read_tsv("Group_with_L5_Maaslin_multivariate.tsv") %>% 
  group_by(Day) %>% 
  mutate(p_adjust = p.adjust(pval, method = "BH")) %>% 
  write_tsv("Group_with_L5_Maaslin_multivariate.tsv")


