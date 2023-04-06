library(tidyverse)
library(broom)
library(ggpubr)
library(ggdist)
source("/data/scripts/utility.R")

confounders = readxl::read_xlsx("confounders.xlsx")


pc = read.pc("bray_curtis_pc.txt")
PC = pc$pc_data %>% as.data.frame() %>% select(PC1:PC5) %>% 
  rownames_to_column(var = "#SampleID")

meta = read_tsv("Metadata_infant_stool_microbiome(ASQ).txt") %>% 
  filter(Group != "VD")

D3 = meta %>% filter(Day == "D3") %>% 
  inner_join(PC, by= "#SampleID") %>% 
  column_to_rownames(var = "ID") %>% 
  select(starts_with("PC"))
D7 = meta %>% filter(Day == "D7") %>% 
  inner_join(PC, by= "#SampleID") %>% 
  column_to_rownames(var = "ID")%>% 
  select(starts_with("PC"))
D30 = meta %>% filter(Day == "D30") %>% 
  inner_join(PC, by= "#SampleID") %>% 
  column_to_rownames(var = "ID") %>% 
  select(starts_with("PC"))
D42 = meta %>% filter(Day == "D42") %>% 
  inner_join(PC, by= "#SampleID") %>% 
  column_to_rownames(var = "ID") %>% 
  select(starts_with("PC"))


sample_front = intersect(rownames(D7), rownames(D3))
D7 = D7[sample_front,]
D3 = D3[sample_front,]
front_mean = (D7+D3)/2

sample_later = intersect(rownames(D30), rownames(D42))
D30 = D30[sample_later,]
D42 = D42[sample_later,]
later_mean = (D30+D42)/2


sample = intersect(rownames(front_mean), rownames(later_mean))
front_mean = front_mean[sample,]
later_mean = later_mean[sample,]

chazhi = later_mean - front_mean



patient = meta %>% distinct(ID, .keep_all = T)

data = chazhi %>% rownames_to_column(var="ID") %>% 
  inner_join(patient, by = "ID")


fit = lm(M6_total ~ PC1, data = data)
summary(fit)
# tidy(fit) %>% filter(term != "(Intercept)") %>% 
#   write_tsv("前后期均值的差值_braycurtis_PC1差值与M6_total线性回归.tsv")
p1 = ggplot(data, aes(M6_total, PC1)) + 
  geom_point(shape = 21, fill = "#7bc2ce", size = 2) + 
  geom_smooth(method = "lm", color ="#7bc2ce", size = 2) + 
  annotate(geom="text", x = -Inf, y = Inf, label = "beta = 75.26\nP = 0.008", hjust =-0.1, vjust = 1.1) + 
  theme_minimal() +
  labs(y = "The change of PC1",
       x = "The ASQ-3 total score at six months") + 
  theme(aspect.ratio = 1, 
        axis.text = element_text(color = "black")) 
ggsave("PC1差值与六月龄总分.pdf", width = 3, height = 3)




data_confounder = chazhi %>% rownames_to_column(var="ID") %>% 
  inner_join(patient, by = "ID") %>% 
  inner_join(confounders, by = c("ID" = "PatientID"))

fit = lm(M6_total ~ PC1 + Gestational_weeks + Birth_weight + Gender, data = data_confounder)

tidy(fit)%>% filter(term != "(Intercept)") %>% 
  write_tsv("前后期均值的差值_braycurtis_PC1差值与M6_total线性回归_adjust_confounders.tsv")


# 组间比较
fit1 = lm(PC1 ~ Group, data = data)
summary(fit1)
tidy(fit)

fit2 = lm(PC1 ~ Group + Gestational_weeks + Birth_weight + Gender, data = data_confounder)
summary(fit2)
tidy(fit)

ggplot(data, aes(Group, PC1)) + 
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot() + 
  stat_compare_means(comparisons = list(c(1,2)), label = "p.signif") + 
  theme_minimal()


library(ggbeeswarm)
ggplot(data, aes(Group, PC1)) + 
  geom_beeswarm(cex = 4) + 
  stat_compare_means(comparisons = list(c(1,2)), label = "p.signif") + 
  theme_minimal()

library(gghalves)
p2 = ggplot(data, aes(Group, PC1, fill = Group)) + 
  geom_half_violin(position = position_nudge(x = 0.1),side=2,alpha = 0.8)+
  geom_point(aes(fill = Group), 
             position = position_jitter(width = 0.05),
             size = 1.5,alpha = 0.8, shape = 21)+
  labs(x=NULL, y = "The change of PC1")+
  geom_boxplot(width = 0.1, outlier.shape = NA,alpha = 0.8) +
  stat_compare_means(comparisons = list(c(1,2)), label = "p.signif") + 
  scale_fill_manual(values = c("#7bc2ce", "#8562cc")) + 
  theme_minimal() + 
  theme(legend.position = "none", 
        axis.text = element_text(color = "black")) 
p2  
wilcox.test(PC1 ~ Group, data = data)
ggsave("PC1差值组间对比图.pdf", width = 3, height = 3)

library(patchwork)
p1+p2
ggsave("plot.pdf", width = 6, height = 3)
