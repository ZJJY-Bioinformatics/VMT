library(tidyverse)
library(broom)
lm_feature = function(i){
  data_fit = data %>% select(`#SampleID`, ID, M6_total, all_of(i)) %>% 
    rename("x" = all_of(i))
  fit = lm(M6_total ~ x, data = data_fit)
  result2 = tidy(fit) %>% 
    filter(term != "(Intercept)") %>% 
    mutate(term = {{i}})
  #smry = summary(fit)

  return(result2)
  
}

metabo288 = read_rds("../interaction/LMM/metabolites288.RDS")
metabolite = read.table("../两组差异代谢物/VMT_Con_metabolome__Day_D42__.txt", sep = "\t", comment.char = "", check.name = F, skip = 1, header = T) %>% 
  filter(`#OTU ID` %in% metabo288)

meta = read_tsv("../../../Metadata_infant_stool_metabolome(ASQ).txt")
confounders = readxl::read_xlsx("拟调整混杂因素(1).xlsx")

#feature_ab = metabolite %>% mutate(across(where(is.numeric), ~if_else(.x>0,1,0)))
#idx_0.1 = rowSums(feature_ab[, -1]) > ceiling((ncol(metabolite)-1)*0.1)
#prevalent_metabolites = metabolite$`#OTU ID`[idx_0.1]

sample <- intersect(meta$`#SampleID`, colnames(metabolite))

feature2 <- metabolite %>% 
  #filter(`#OTU ID` %in% prevalent_metabolites) %>% 
  select(`#OTU ID`, all_of(sample)) %>%
  sjmisc::rotate_df(cn = TRUE, rn = "#SampleID")

data = meta %>% 
  inner_join(feature2, by = "#SampleID") %>% 
  inner_join(confounders, by = c("ID" = "PatientID"))


metabolite_name = colnames(feature2)[-1]
#i = "5_Hydroxylysine"



#D3
total = map_dfr(metabolite_name, ~lm_feature(.x))
fdr = p.adjust(total$p.value, method = "fdr")
total$FDR = fdr
total$Day = "D3"
total_D3 = total

#D7
total = map_dfr(metabolite_name, ~lm_feature(.x))
fdr = p.adjust(total$p.value, method = "fdr")
total$FDR = fdr
total$Day = "D7"
total_D7 = total


#D30
total = map_dfr(metabolite_name, ~lm_feature(.x))
fdr = p.adjust(total$p.value, method = "fdr")
total$FDR = fdr
total$Day = "D30"
total_D30 = total

#D42
total = map_dfr(metabolite_name, ~lm_feature(.x))
fdr = p.adjust(total$p.value, method = "fdr")
total$FDR = fdr
total$Day = "D42"
total_D42 = total


all = bind_rows(total_D3, total_D7, total_D30, total_D42)
write_tsv(all, "各代谢物与M6_total线性回归(univariate)0917.tsv")


#========multivariate==============

lm_feature_multi = function(i){
  data_fit = data %>% select(`#SampleID`, ID, M6_total, all_of(i), Gestational_weeks, Birth_weight, Gender, starts_with("Feeding")) %>% 
    rename("x" = all_of(i))
  fit = lm(M6_total ~ x + Gestational_weeks + Birth_weight + Gender + Feeding_42, data = data_fit)
  result2 = tidy(fit) %>% 
    filter(term == "x") %>% 
    mutate(term = {{i}})
  return(result2)
}


metabolite = read.table("../两组差异代谢物/VMT_Con_metabolome__Day_D42__.txt", sep = "\t", comment.char = "", check.name = F, skip = 1, header = T) %>% 
  filter(`#OTU ID` %in% metabo288)


meta = read_tsv("../../../Metadata_infant_stool_metabolome(ASQ).txt")
confounders = readxl::read_xlsx("C:\\data\\Documents\\BaiduNetdiskWorkspace\\newborn_NH\\202207revise\\16S\\原来\\maaslin\\拟调整混杂因素(1).xlsx")

# feature_ab = metabolite %>% mutate(across(where(is.numeric), ~if_else(.x>0,1,0)))
# idx_0.1 = rowSums(feature_ab[, -1]) > ceiling((ncol(metabolite)-1)*0.1)
# prevalent_metabolites = metabolite$`#OTU ID`[idx_0.1]

sample <- intersect(meta$`#SampleID`, colnames(metabolite))

feature2 <- metabolite %>% 
  #filter(`#OTU ID` %in% prevalent_metabolites) %>% 
  select(`#OTU ID`, all_of(sample)) %>%
  sjmisc::rotate_df(cn = TRUE, rn = "#SampleID")

data = meta %>% 
  inner_join(feature2, by = "#SampleID") %>% 
  inner_join(confounders, by = c("ID" = "PatientID"))


metabolite_name = colnames(feature2)[-1]
#i = "5_Hydroxylysine"



#D3
total = map_dfr(metabolite_name, ~lm_feature_multi(.x))
fdr = p.adjust(total$p.value, method = "fdr")
total$FDR = fdr
total$Day = "D3"
total_D3 = total

#D7
total = map_dfr(metabolite_name, ~lm_feature_multi(.x))
fdr = p.adjust(total$p.value, method = "fdr")
total$FDR = fdr
total$Day = "D7"
total_D7 = total


#D30
total = map_dfr(metabolite_name, ~lm_feature_multi(.x))
fdr = p.adjust(total$p.value, method = "fdr")
total$FDR = fdr
total$Day = "D30"
total_D30 = total

#D42
total = map_dfr(metabolite_name, ~lm_feature_multi(.x))
fdr = p.adjust(total$p.value, method = "fdr")
total$FDR = fdr
total$Day = "D42"
total_D42 = total


all = bind_rows(total_D3, total_D7, total_D30, total_D42)
write_tsv(all, "288代谢物与M6_total线性回归(multivariate)0917.tsv")



res = read_tsv("各代谢物与M6_total线性回归(multivariate).tsv", skip = 1) %>% 
  filter(Day != "D1" & Day != "D2") %>% 
  mutate(all_FDR = p.adjust(`Pr(>|t|)`, method = "fdr"))

write_tsv(res, "D3-D42代谢物与M6_total线性回归(multivariate).tsv")


res = read_tsv("各代谢物与M6_total线性回归(univariate).tsv") %>% 
  filter(Day != "D1" & Day != "D2") %>% 
  mutate(all_FDR = p.adjust(`Pr(>|t|)`, method = "fdr"))

write_tsv(res, "D3-D42代谢物与M6_total线性回归(univariate).tsv")
