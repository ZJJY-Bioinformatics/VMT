library(tidyverse)
source("/data/scripts/utility.R")

meta = read_tsv("../../VMT_Con_metadata_infant_stool_microbiome.txt")

pc = read.pc("bdiv/bray_curtis_pc.txt")
PC = pc$pc_data %>% as.data.frame() %>% select(PC1:PC5) %>% 
  rownames_to_column(var = "#SampleID")
proportion = pc$pc_percent

PC_plot = meta %>% select(`#SampleID`, Day, Time_point, Group) %>% 
  inner_join(PC, by = "#SampleID") %>% 
  mutate(Day = fct_relevel(Day, c("D3", "D7", "D30", "D42")))


## 求各时间点双峰分布中间的最小值的坐标，以此来对样本进行分组
each_day_group_min_x = function(day){
  Con = PC_plot %>% filter(Day == day & Group == "Con") %>% pull(PC1)
  VMT = PC_plot %>% filter(Day == day & Group == "VMT") %>% pull(PC1)
  
  ####Con
  Con_dst = density(Con) 
  Con_dst_df = data.frame(x = Con_dst$x, y = Con_dst$y)
  #左峰最大值
  Con_dst_left = Con_dst_df %>% filter(x<0.2) 
  Con_dst_left_max = Con_dst_left %>% pull(y) %>% max()
  Con_dst_left_max_x = Con_dst_left$x[which(Con_dst_left$y == Con_dst_left_max)]
  
  #右峰最大值
  Con_dst_right = Con_dst_df %>% filter(x >= 0.2)
  Con_dst_right_max = Con_dst_right %>% pull(y) %>% max()
  Con_dst_right_max_x = Con_dst_right$x[which(Con_dst_right$y == Con_dst_right_max)]
  
  # 两峰之间的最小值
  middle = Con_dst_df %>% filter(x > Con_dst_left_max_x & x < Con_dst_right_max_x)
  middle_min = middle %>% pull(y) %>% min()
  Con_middle_min_x = middle$x[which(middle$y == middle_min)]
  
  #####VMT
  VMT_dst = density(VMT)
  VMT_dst_df = data.frame(x = VMT_dst$x, y = VMT_dst$y)
  #左峰最大值
  VMT_dst_left = VMT_dst_df %>% filter(x<0.2) 
  VMT_dst_left_max = VMT_dst_left %>% pull(y) %>% max()
  VMT_dst_left_max_x = VMT_dst_left$x[which(VMT_dst_left$y == VMT_dst_left_max)]
  
  #右峰最大值
  VMT_dst_right = VMT_dst_df %>% filter(x >= 0.2)
  VMT_dst_right_max = VMT_dst_right %>% pull(y) %>% max()
  VMT_dst_right_max_x = VMT_dst_right$x[which(VMT_dst_right$y == VMT_dst_right_max)]
  
  # 两峰之间的最小值
  middle = VMT_dst_df %>% filter(x > VMT_dst_left_max_x & x < VMT_dst_right_max_x)
  middle_min = middle %>% pull(y) %>% min()
  VMT_middle_min_x = middle$x[which(middle$y == middle_min)]
  
  return(list(Con_middle_min_x = Con_middle_min_x, VMT_middle_min_x = VMT_middle_min_x))
  
}
D3_x = each_day_group_min_x("D3")
D7_x = each_day_group_min_x("D7")
D30_x = each_day_group_min_x("D30")
D42_x = each_day_group_min_x("D42")

day ="D3"
x = D3_x

D3_Con = PC_plot %>% filter(Day == day & Group == "Con") %>% 
  mutate(peak = if_else(PC1 < median(x$Con_middle_min_x), "Left", "Right"))
D3_VMT = PC_plot %>% filter(Day == day & Group == "VMT") %>% 
  mutate(peak = if_else(PC1 < median(x$VMT_middle_min_x), "Left", "Right"))

each_day_group = function(day, x){
  Day_Con = PC_plot %>% filter(Day == day & Group == "Con") %>% 
    mutate(peak = if_else(PC1 < median(x$Con_middle_min_x), "Left", "Right"))
  Day_VMT = PC_plot %>% filter(Day == day & Group == "VMT") %>% 
    mutate(peak = if_else(PC1 < median(x$VMT_middle_min_x), "Left", "Right"))
  
  res = bind_rows(Day_VMT, Day_Con)
  return(res)
}
D3 = each_day_group("D3", D3_x)
D7 = each_day_group("D7", D7_x)
D30 = each_day_group("D30", D30_x)
D42 = each_day_group("D42", D42_x)

all = bind_rows(D3,D7,D30,D42)




##############双峰和菌型的分布################

ET = read_delim("../enterotype/all.2_clusters.data.cluster.txt", delim = " ")

genus = read.table("infant_stool_e10k_L6.txt", skip = 1, comment.char = "", sep = "\t", header = T, check.names = F) %>% 
  sjmisc::rotate_df(cn = T) %>% 
  rownames_to_column(var = "#SampleID") 


ET$`#SampleID` = genus$`#SampleID`
data = ET %>% select(`#SampleID`, ET) %>% 
  inner_join(all, by = "#SampleID")

data %>% filter(Day == "D42")

ggplot(data, aes(factor(ET), fill = peak)) + 
  geom_bar(stat = "count", position = "fill") + 
  facet_grid(.~Day) + 
  scale_fill_manual(values = c("#7bc2ce", "#8562cc")) + 
  theme_minimal() + 
  labs(x = "ET", y = "Proportion")
library(ggstatsplot)
library(patchwork)
p1 = ggbarstats(data = filter(data, Day == "D3"), x = peak, y = ET) 
p2 = ggbarstats(data = filter(data, Day == "D7"), x = peak, y = ET) 
p3 = ggbarstats(data = filter(data, Day == "D30"), x = peak, y = ET) 
p4 = ggbarstats(data = filter(data, Day == "D42"), x = peak, y = ET) 
p1+p2+p3+p4
ggsave("两种菌型中双峰各占比例.pdf", width = 8, height = 8)

p1 = ggbarstats(data = filter(data, Day == "D3"), x = peak, y = Group) 
p2 = ggbarstats(data = filter(data, Day == "D7"), x = peak, y = Group) 
p3 = ggbarstats(data = filter(data, Day == "D30"), x = peak, y = Group) 
p4 = ggbarstats(data = filter(data, Day == "D42"), x = peak, y = Group) 
p1+p2+p3+p4
ggsave("两组中双峰各占比例.pdf", width = 8, height = 8)


p1 = ggbarstats(data = filter(data, Day == "D3"), x = Group, y = peak) 
p2 = ggbarstats(data = filter(data, Day == "D7"), x = Group, y = peak) 
p3 = ggbarstats(data = filter(data, Day == "D30"), x = Group, y = peak) 
p4 = ggbarstats(data = filter(data, Day == "D42"), x = Group, y = peak) 
p1+p2+p3+p4
ggsave("双峰中两组的比例各占比例.pdf", width = 8, height = 8)


p1 = ggbarstats(data = filter(data, Day == "D3"), x = Group, y = peak) 
p2 = ggbarstats(data = filter(data, Day == "D7"), x = Group, y = peak) 
p3 = ggbarstats(data = filter(data, Day == "D30"), x = Group, y = peak) 
p4 = ggbarstats(data = filter(data, Day == "D42"), x = Group, y = peak) 
p1+p2+p3+p4
ggsave("双峰中两组的比例各占比例.pdf", width = 8, height = 8)



p1 = ggbarstats(data = filter(data, Day == "D3"), x = ET, y = Group) 
p2 = ggbarstats(data = filter(data, Day == "D7"), x = ET, y = Group) 
p3 = ggbarstats(data = filter(data, Day == "D30"), x = ET, y = Group) 
p4 = ggbarstats(data = filter(data, Day == "D42"), x = ET, y = Group) 
p1+p2+p3+p4
ggsave("两组中两个肠型的比例.pdf", width = 8, height = 8)

ET = read_tsv("../enterotype/ET_meta.tsv")
ggbarstats(data = filter(ET, Group != "VD"), x = ET, y = Group, ggtheme = theme_minimal()) +
  scale_fill_manual(values = c('#1f77b4', '#ff7f0e')) + 
  theme(aspect.ratio = 1)
  
ggsave("fig1d-两组中所有样本两个肠型的比例.pdf", width = 4, height = 4)


ET = ET %>% mutate(delivery = if_else(Group == "VD", "VD", "CS"))
ggbarstats(data = ET, x = ET, y = delivery, ggtheme = theme_minimal()) +
  scale_fill_manual(values = c('#1f77b4', '#ff7f0e')) + 
  theme(aspect.ratio = 1)
ggsave("supp_fig2-不同分娩方式所有样本两个肠型的比例.pdf", width = 4, height = 4)



#JSD分型，两组分型

ET = read_delim("../enterotype/all.2_clusters.data.cluster.txt")

data = ET %>% select(X.SampleID, ET) %>% 
  inner_join(all, by = c("X.SampleID" = "#SampleID"))

data %>% filter(Day == "D42")

ggplot(data, aes(factor(ET), fill = peak)) + 
  geom_bar(stat = "count", position = "fill") + 
  facet_grid(.~Day) + 
  scale_fill_manual(values = c("#7bc2ce", "#8562cc")) + 
  theme_minimal() + 
  labs(x = "ET", y = "Proportion")
library(ggstatsplot)
library(patchwork)

p1 = ggbarstats(data = filter(data, Day == "D3"), x = peak, y = ET) 
p2 = ggbarstats(data = filter(data, Day == "D7"), x = peak, y = ET) 
p3 = ggbarstats(data = filter(data, Day == "D30"), x = peak, y = ET) 
p4 = ggbarstats(data = filter(data, Day == "D42"), x = peak, y = ET) 
p1+p2+p3+p4
ggsave("两种菌型中双峰各占比例.pdf", width = 8, height = 8)

p1 = ggbarstats(data = filter(data, Day == "D3"), x = peak, y = Group) 
p2 = ggbarstats(data = filter(data, Day == "D7"), x = peak, y = Group) 
p3 = ggbarstats(data = filter(data, Day == "D30"), x = peak, y = Group) 
p4 = ggbarstats(data = filter(data, Day == "D42"), x = peak, y = Group) 
p1+p2+p3+p4
ggsave("两组中双峰各占比例.pdf", width = 8, height = 8)


p1 = ggbarstats(data = filter(data, Day == "D3"), x = Group, y = peak) 
p2 = ggbarstats(data = filter(data, Day == "D7"), x = Group, y = peak) 
p3 = ggbarstats(data = filter(data, Day == "D30"), x = Group, y = peak) 
p4 = ggbarstats(data = filter(data, Day == "D42"), x = Group, y = peak) 
p1+p2+p3+p4
ggsave("双峰中两组的比例各占比例.pdf", width = 8, height = 8)





################
# Gaussian mix model

library(mixtools)
m <- normalmixEM(
     PC_plot$PC1)
summary(m)

m_diff_sq <- function (x) {
  m_diff <- dnorm(x, mean = m$mu[1], sd = m$sigma[1]) - 
    dnorm(x, mean = m$mu[2], sd = m$sigma[2])
  m_diff ^ 2
}
optim(list(x=0.5), m_diff_sq)



gmm_dist <- function(xs, m, k) {
  m$lambda[k] * dnorm(xs, mean = m$mu[k], m$sigma[k])
}
s0_clinical_gmm_pred <- data_frame(PC1=seq(min(PC_plot$PC1), max(PC_plot$PC1), 0.01)) %>%
  mutate(comp.1 = gmm_dist(PC1, m, 1)) %>%
  mutate(comp.2 = gmm_dist(PC1, m, 2)) %>%
  gather(Component, Probability, starts_with("comp")) %>%
  mutate(PredictedDensity = nrow(PC_plot) * Probability * 0.1)


PC_plot %>%
  cbind(m$posterior) %>%
  mutate(Component = if_else(comp.1 > 0.5, "comp.1", "comp.2")) %>%
  ggplot(aes(x=PC1)) +
  geom_histogram(aes(fill=Component), binwidth = .01, boundary=0) +
  geom_line(
    aes(y=PredictedDensity, color=Component), 
    data=s0_clinical_gmm_pred) +
  geom_vline(xintercept=0.55, linetype="dashed", color="#666666") +
  # scale_color_manual(values=host_palette) +
  # scale_fill_manual(values=host_palette) +
  labs(x="Fraction of human reads", y="Number of samples") +
  theme_bw()



mboot <- boot.comp(
  PC_plot$PC1, B=1000,
  max.comp=1, mix.type = "normalmix")



mboot$obs.log.lik
mboot$p.values

data_frame(GmmLogLik = mboot$log.lik[[1]]) %>%
  ggplot() +
  geom_histogram(aes(x=GmmLogLik), binwidth = 2, boundary=0) +
  geom_vline(xintercept = mboot$obs.log.lik, linetype="dashed") +
  annotate(
    "text", x=15, y=100, 
    label="Bootstrap\nsimulation\nvalues", hjust=0) +
  annotate(
    "text", x=118, y=200, 
    label="Observed value: 133.284\nP < 0.001", hjust=1) +
  theme_bw()

