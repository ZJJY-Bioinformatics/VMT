library(vegan)
library(tidyverse)
library(ggsci)
dist1 = read.table("bray_curtis_dm.txt", row.names = 1) %>% 
  select(ends_with("BSD42"))
dist2 = read.table("bray_curtis_dm_metabolome_absolute.txt", row.names = 1) %>% 
  select(ends_with("BSD42"))

sample = intersect(colnames(dist1), colnames(dist2))

dist1 = dist1[sample, sample]
dist2 = dist2[sample, sample]


bc1 = as.dist(dist1)
bc2 = as.dist(dist2)



pcoa1 = cmdscale(bc1,  k = (nrow(dist1) - 1))
pcoa2 = cmdscale(bc2,  k = (nrow(dist2) - 1))


pro = procrustes(pcoa1, pcoa2, symmetric = T)

summary(pro)
set.seed(2022)
prot = protest(pcoa1, pcoa2, permutations = 999)
prot$signif
prot$ss


Pro_Y <- cbind(data.frame(pro$Yrot), data.frame(pro$X))
Pro_X <- data.frame(pro$rotation)

D42 = ggplot(Pro_Y) + 
  geom_vline(xintercept = 0, color = 'gray', linetype = 2) +
  geom_hline(yintercept = 0, color = 'gray', linetype = 2) +
  geom_abline(intercept = 0, slope = Pro_X[1,2]/Pro_X[1,1]) +
  geom_abline(intercept = 0, slope = Pro_X[2,2]/Pro_X[2,1]) +
  geom_segment(aes(x = X1, y = X2,
                   xend = (X1 + Dim1)/2, yend = (X2 + Dim2)/2),
               # geom_segment 绘制两点间的直线
              # arrow = arrow(length = unit(0, 'cm')),
               color = "#BC3C29", alpha = 0.5) +
  geom_segment(aes(x = (X1 + Dim1)/2, y = (X2 + Dim2)/2,
                   xend = Dim1, yend = Dim2),
               #arrow = arrow(length = unit(0.1, 'cm')),
               color = "#0072B5", alpha = 0.5) +
  geom_point(aes(X1, X2), color = "#BC3C29") + 
  geom_point(aes(Dim1, Dim2), color = "#0072B5") + 
  annotate('text', label = paste0("M^2=", round(prot$ss, 3), "\nP=", prot$signif), x = -Inf, y = Inf, vjust = 1.1, hjust = -0.1) + 
  theme_minimal() +
  labs(x = "Dimension1", y = "Dimension2") +
  theme(axis.text = element_text(color = "black"))
D42
ggsave("supp_fig8_所有样本菌群和代谢物(absolute)的procrustes analysis.pdf", width = 4, height = 4)



library(patchwork)
all + D3 + D7 + D30 + D42 + plot_annotation(tag_levels = "a")
ggsave("supp_fig_菌群和代谢物(absolute)的procrustes analysis.pdf", width = 9, height = 6)

