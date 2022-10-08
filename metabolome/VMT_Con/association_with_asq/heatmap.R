library(tidyverse)
library(RColorBrewer)
lm_res = read_tsv("288代谢物与M6_total线性回归(multivariate)0917.tsv") 


category = readxl::read_xlsx("代谢物类别.xlsx") 

# FDR小于0.25的feature
fdr0.25 = lm_res %>% filter(FDR < 0.25) %>% 
  pull(term)

lm_res2 = lm_res %>% filter(term %in% fdr0.25) %>% 
  mutate(fulog10fdr = -log10(FDR), 
         direction = sign(estimate), 
         color_deep = fulog10fdr*direction, 
         Day = fct_relevel(Day, c("D3", "D7", "D30", "D42")), 
         zscore = -qnorm(FDR/2), 
         color_deep2 = zscore*direction, 
         color_deep3 = cut(color_deep2, c(-4, -3, -2, -1, 0, 1)),
         color_deep4 = cut(color_deep, c(-3.5, -2, -1, 0,1))) %>% 
  inner_join(category, by = c("term" = "metabolites")) %>% 
  arrange(category) %>%
  mutate(term = fct_relevel(term, unique(term)))

label = lm_res %>% filter(FDR<0.25) %>% 
  mutate(labels = if_else(estimate > 0, "+", "-"))


ggplot(lm_res2, aes(Day, term)) + 
  geom_point(aes(fill = color_deep4, size = fulog10fdr), shape = 21, key_glyph = draw_key_rect) + 
  geom_text(data=label, mapping = aes(Day, term, label = labels), size = 5) +
  scale_fill_manual(values = rev(c("#e08214", "#8073ac", "#542788", "#2d004b"))) + 
  theme_bw() +
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 10), 
        axis.text = element_text(color = "black")) +
  scale_size_continuous(range = c(5,10)) + 
  labs(fill = "-log(FDR) * effect direction", 
       size = "-log(FDR)") +
  guides(fill = guide_bins(axis = F, show.limits = T, reverse = T), 
         size = guide_none())
#ggsave("多元线性回归与ASQ相关性.pdf", width = 6, height = 5)  
ggsave("多元线性回归与ASQ相关性(排序).pdf", width = 6, height = 5)  
