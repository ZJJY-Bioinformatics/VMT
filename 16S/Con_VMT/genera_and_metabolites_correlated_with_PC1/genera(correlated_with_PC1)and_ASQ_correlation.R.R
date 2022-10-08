library(tidyverse)
library(rstatix)
library(broom)
source("/data/scripts/utility.R")
confounders = readxl::read_xlsx("C:\\data\\Documents\\BaiduNetdiskWorkspace\\newborn_NH\\202207revise\\16S\\原来\\maaslin\\拟调整混杂因素(1).xlsx")


taxa = read_tsv("../../taxa/infant_stool_e10k_L6.txt", skip = 1) %>% 
  filter(!str_detect(`#OTU ID`, ".*Other")) %>% 
  filter(!str_detect(`#OTU ID`, ".*g__$")) 

taxa_ap = taxa %>% mutate(across(where(is.numeric), ~if_else(.x>0,1,0)))
idx_0.1 = rowSums(taxa_ap[, -1]) > 23
taxa = taxa[idx_0.1, ]
taxa$`#OTU ID` = taxa_shortname(taxa$`#OTU ID`)



meta = read_tsv("../../../../Metadata_infant_stool_microbiome(ASQ).txt") %>% 
  filter(Day != "D1" & Day != "D2") %>% 
  filter(Group != "VD")

corr = read_tsv("与PC1相关的菌—L6.tsv") %>% 
  filter(FDR < 0.25)

#与菌群PC1有相关性的代谢物与ASQ的相关性  
cor_sig = taxa %>% filter(`#OTU ID` %in% corr$feature) %>% 
  sjmisc::rotate_df(cn = T, rn = "#SampleID") 

data = cor_sig %>% inner_join(meta, by = "#SampleID")

feature = colnames(cor_sig)[-1]

lm_feature = function(day, i){
  fit_data = data %>% filter(Day == paste0("D", day)) %>% 
    rename("x" = all_of(i))
  
  fit = lm(M6_total ~ x, data = fit_data)
  res = tidy(fit) %>% 
    filter(term == "x") %>% 
    mutate(y = "M6_total", x = {{i}}, Day = paste0("D", day)) %>% 
    select(y, x, estimate, std.error, statistic, p.value, Day)
  return(res)
  
}

lm_feature(3, feature[2])

D3 = map_dfr(feature, ~lm_feature(3, .x))
D3 = D3 %>% mutate(FDR = p.adjust(p.value, method = "fdr"))
D7 = map_dfr(feature, ~lm_feature(7, .x))
D7 = D7 %>% mutate(FDR = p.adjust(p.value, method = "fdr"))
D30 = map_dfr(feature, ~lm_feature(30, .x))
D30 = D30 %>% mutate(FDR = p.adjust(p.value, method = "fdr"))
D42 = map_dfr(feature, ~lm_feature(42, .x))
D42 = D42 %>% mutate(FDR = p.adjust(p.value, method = "fdr"))

all = bind_rows(D3, D7, D30, D42) %>% 
  mutate(all_FDR = p.adjust(p.value, method = "fdr"))
write_tsv(all, "与PC1相关的genus和M6_total线性回归(univariate).tsv")

## multivariable
data2 = cor_sig %>% inner_join(meta, by = "#SampleID") %>% 
  inner_join(confounders, by = c("ID" = "PatientID"))

multi_lm_feature = function(day, i){
  fit_data = data2 %>% filter(Day == paste0("D", day)) %>% 
    rename("x" = all_of(i))
  fml = as.formula(paste0("M6_total ~ x + Gestational_weeks + Birth_weight + Gender + Feeding_", day))
  fit = lm(fml, data = fit_data)
  res = tidy(fit) %>% 
    filter(term == "x") %>% 
    mutate(y = "M6_total", x = {{i}}, Day = paste0("D", day)) %>% 
    select(y, x, estimate, std.error, statistic, p.value, Day)
  return(res)
  
}

multi_lm_feature(3, feature[2])

D3 = map_dfr(feature, ~multi_lm_feature(3, .x))
D3 = D3 %>% mutate(FDR = p.adjust(p.value, method = "fdr"))
D7 = map_dfr(feature, ~multi_lm_feature(7, .x))
D7 = D7 %>% mutate(FDR = p.adjust(p.value, method = "fdr"))
D30 = map_dfr(feature, ~multi_lm_feature(30, .x))
D30 = D30 %>% mutate(FDR = p.adjust(p.value, method = "fdr"))
D42 = map_dfr(feature, ~multi_lm_feature(42, .x))
D42 = D42 %>% mutate(FDR = p.adjust(p.value, method = "fdr"))

all = bind_rows(D3, D7, D30, D42) %>% 
  mutate(all_FDR = p.adjust(p.value, method = "fdr"))
#write_tsv(all, "与PC1相关的genus和M6_total线性回归(multivariate).tsv")




## 与PC1相关的菌和ASQ的关系
data = read_tsv("与PC1相关的genus和M6_total线性回归(univariate).tsv")

lm_res = data %>% 
  mutate(fulog10fdr = -log10(FDR), 
         direction = sign(estimate), 
         color_deep = fulog10fdr*direction, 
         Day = fct_relevel(Day, c("D3", "D7", "D30", "D42")), 
         zscore = -qnorm(FDR/2), 
         color_deep2 = zscore*direction, 
         color_deep3 = cut(color_deep2, c(-4, -3, -2, -1, 0, 1)),
         color_deep4 = cut(color_deep, c(-3,  -1, 0, 1)))

label = lm_res %>% filter(FDR<0.25) %>% 
  mutate(labels = if_else(estimate > 0, "+", "-"))

ggplot(lm_res, aes(Day, x)) + 
  geom_point(aes(fill = color_deep4, size = fulog10fdr), shape = 21, key_glyph = draw_key_rect) + 
  geom_text(data=label, mapping = aes(Day, x, label = labels), size = 5) +
  scale_fill_manual(values = rev(c("#d6604d",  "#4393c3", "#2166ac"))) + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 10, face = "italic"), 
        axis.text = element_text(color = "black"), 
        aspect.ratio = 10/4, 
        panel.border = element_rect(fill = "transparent", color = "black")) +
  scale_size_continuous(range = c(5,10)) + 
  labs(fill = "-log(FDR) * effect direction", 
       size = "-log(FDR)")+ 
  guides(fill = guide_bins(axis = F, show.limits = T, reverse = T), 
         size = guide_none())

ggsave("与PC1相关的genus和M6_total的相关性.pdf", width = 4.5, height = 5)  



