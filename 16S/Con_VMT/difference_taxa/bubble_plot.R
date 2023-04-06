library(tidyverse)
##########属水平###################
#univariate
uni = read_tsv("Group_with_L6_Maaslin_univariate.tsv")

fdr0.1 = uni %>% filter(qval < 0.1) %>% 
  pull(feature)

uni_res2 = uni %>% filter(feature %in% fdr0.1) %>% 
  mutate(fulog10fdr = -log10(qval),
         direction = sign(coef), 
         color_deep = fulog10fdr * direction, 
         Day = fct_relevel(Day, c("Day3", "Day7", "Day30", "Day42")),
         color_deep4 = cut(color_deep, c(-1.5, -1, -0.5, 0, 0.5, 1, 3)))
label = uni %>% filter(qval<0.1) %>% 
  mutate(labels = if_else(coef > 0, "*", "*"))

pal = colorRampPalette(c("#2d004b", "white","#7f3b08" ))(6)
ggplot(uni_res2, aes(Day, feature)) + 
  geom_point(aes(fill = color_deep4, size = fulog10fdr), shape = 21, key_glyph = draw_key_rect) + 
  geom_text(data=label, mapping = aes(Day, feature, label = labels), size = 5, color = "white") +
  scale_fill_manual(values = rev(c("#7f3b08", "#b35806", "#e08214", "#8073ac", "#542788", "#2d004b"))) + 
  #scale_fill_gradient2(mid = "white", low = "#542788", high = "#b35806", midpoint = 0) +
  theme_bw() +
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 10), 
        axis.text = element_text(color = "black"), aspect.ratio = 24/6) +
  scale_size_continuous(range = c(5,10)) + 
  labs(fill = "-log(FDR) * effect direction", 
       size = "-log(FDR)") + 
  scale_x_discrete(labels = c("D3", "D7", "D30", "D42")) + 
  guides(fill = guide_bins(axis = F, show.limits = F, reverse = T), 
         size = guide_none())
ggsave("maaslin2_两组间差异菌_univarite.pdf", width = 6, height = 7)


##multi
multi = read_tsv("Group_with_L6_Maaslin_multivariate.tsv")

fdr0.1 = multi %>% filter(p_adjust < 0.1) %>% 
  pull(feature)

multi_res2 = multi %>% filter(feature %in% fdr0.1) %>% 
  mutate(fulog10fdr = -log10(p_adjust),
         direction = sign(coef), 
         color_deep = fulog10fdr * direction, 
         Day = fct_relevel(Day, c("Day3", "Day7", "Day30", "Day42")),
         color_deep4 = cut(color_deep, c(-1.5, -1, -0.5, 0, 0.5, 1, 2))
         )
label = multi %>% filter(p_adjust<0.1) %>% 
  mutate(labels = if_else(coef > 0, "*", "*"))


ggplot(multi_res2, aes(Day, feature)) + 
  geom_point(aes(fill = color_deep4, size = fulog10fdr), shape = 21, key_glyph = draw_key_rect) + 
  geom_text(data=label, mapping = aes(Day, feature, label = labels), size = 5, color = "white") +
  scale_fill_manual(values = rev(c("#7f3b08", "#b35806", "#e08214", "#8073ac", "#542788", "#2d004b"))) + 
  #scale_fill_gradient2(mid = "white", low = "#542788", high = "#b35806", midpoint = 0) +
  theme_bw() +
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 10), 
        axis.text = element_text(color = "black"), aspect.ratio = 28/6) +
  scale_size_continuous(range = c(5,10)) + 
  labs(fill = "-log(FDR) * effect direction", 
       size = "-log(FDR)") + 
  scale_x_discrete(labels = c("D3", "D7", "D30", "D42")) + 
  guides(fill = guide_bins(axis = F, show.limits = F, reverse = T), 
         size = guide_none())
ggsave("maaslin2_两组间差异菌_mulivarite.pdf", width = 7, height = 8)


############科水平###############################

#univariate
uni = read_tsv("Group_with_L5_Maaslin_univariate.tsv")

fdr0.1 = uni %>% filter(qval < 0.1) %>% 
  pull(feature)

uni_res2 = uni %>% filter(feature %in% fdr0.1) %>% 
  mutate(fulog10fdr = -log10(qval),
         direction = sign(coef), 
         color_deep = fulog10fdr * direction, 
         Day = fct_relevel(Day, c("Day3", "Day7", "Day30", "Day42")),
         color_deep4 = cut(color_deep, c(-1.5, -1, -0.5, 0, 0.5, 1, 3)))
label = uni %>% filter(qval<0.1) %>% 
  mutate(labels = if_else(coef > 0, "*", "*"))


ggplot(uni_res2, aes(Day, feature)) + 
  geom_point(aes(fill = color_deep4, size = fulog10fdr), shape = 21, key_glyph = draw_key_rect) + 
  geom_text(data=label, mapping = aes(Day, feature, label = labels), size = 5, color = "white") +
  scale_fill_manual(values = rev(c("#7f3b08", "#b35806", "#e08214", "#8073ac", "#542788", "#2d004b"))) + 
  #scale_fill_gradient2(mid = "white", low = "#542788", high = "#b35806", midpoint = 0) +
  theme_bw() +
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 10), 
        axis.text = element_text(color = "black"), aspect.ratio = 17/6) +
  scale_size_continuous(range = c(5,10)) + 
  labs(fill = "-log(FDR) * effect direction", 
       size = "-log(FDR)") + 
  scale_x_discrete(labels = c("D3", "D7", "D30", "D42")) + 
  guides(fill = guide_bins(axis = F, show.limits = F, reverse = T), 
         size = guide_none())
ggsave("maaslin2_两组间差异菌_univarite_L5.pdf", width = 5.5, height = 6)


##multi
multi = read_tsv("Group_with_L5_Maaslin_multivariate.tsv")

fdr0.1 = multi %>% filter(p_adjust < 0.1) %>% 
  pull(feature)

multi_res2 = multi %>% filter(feature %in% fdr0.1) %>% 
  mutate(fulog10fdr = -log10(p_adjust),
         direction = sign(coef), 
         color_deep = fulog10fdr * direction, 
         Day = fct_relevel(Day, c("Day3", "Day7", "Day30", "Day42")),
         color_deep4 = cut(color_deep, c(-1.8, -1, -0.5, 0, 0.5, 1, 2.5))
  )
label = multi %>% filter(p_adjust<0.1) %>% 
  mutate(labels = if_else(coef > 0, "*", "*"))


ggplot(multi_res2, aes(Day, feature)) + 
  geom_point(aes(fill = color_deep4, size = fulog10fdr), shape = 21, key_glyph = draw_key_rect) + 
  geom_text(data=label, mapping = aes(Day, feature, label = labels), size = 5, color = "white") +
  scale_fill_manual(values = rev(c("#7f3b08", "#b35806", "#e08214", "#8073ac", "#542788", "#2d004b"))) + 
  #scale_fill_gradient2(mid = "white", low = "#542788", high = "#b35806", midpoint = 0) +
  theme_bw() +
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 10), 
        axis.text = element_text(color = "black"), aspect.ratio = 22/6) +
  scale_size_continuous(range = c(5,10)) + 
  labs(fill = "-log(FDR) * effect direction", 
       size = "-log(FDR)") + 
  scale_x_discrete(labels = c("D3", "D7", "D30", "D42")) + 
  guides(fill = guide_bins(axis = F, show.limits = F, reverse = T), 
         size = guide_none())
ggsave("maaslin2_两组间差异菌_mulivarite_L5.pdf", width = 6, height = 6)

