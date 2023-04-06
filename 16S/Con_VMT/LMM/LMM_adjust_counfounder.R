library(tidyverse)
library(lmerTest)
library(broom.mixed)
source("/data/scripts/utility.R")
confounder = readxl::read_xlsx("C:\\data\\Documents\\BaiduNetdiskWorkspace\\newborn_NH\\202207revise\\16S\\原来\\maaslin\\拟调整混杂因素(1).xlsx")

taxa = read_tsv("../../../taxa/infant_stool_e10k_L6.txt", skip = 1) %>% 
  filter(!str_detect(`#OTU ID`, ".*Other")) %>% 
  filter(!str_detect(`#OTU ID`, ".*g__$")) %>% 
  select(!(ends_with("BSD1") | ends_with("BSD2")))
enterb = read_tsv("../../../taxa/infant_stool_e10k_L6.txt", skip = 1) %>% 
  filter(str_detect(`#OTU ID`, ".*f__Enterobacteriaceae;Other"))%>% 
  select(!(ends_with("BSD1") | ends_with("BSD2")))
taxa = taxa %>% bind_rows(enterb)

taxa_ap = taxa %>% mutate(across(where(is.numeric), ~if_else(.x>0,1,0)))
idx_0.1 = rowSums(taxa_ap[, -1]) > 23
taxa = taxa[idx_0.1, ]

meta = read_tsv("../../VMT_Con_metadata_infant_stool_microbiome.txt") %>% 
  #select(`#SampleID`, ID, Group, Day) %>% 
  filter(!Day %in% c("D1", "D2")) %>% 
  mutate(Day = fct_relevel(Day, c("D1", "D2", "D3", "D7", "D30", "D42")), 
         Time = as.numeric(Day),
         Day_time = as.numeric(str_extract(Day, "\\d+")))

feeding3 = meta %>% 
  filter(Day == "D3") %>% 
  select(`#SampleID`, Feeding_3) %>% 
  rename(Feeding = Feeding_3)
feeding7 = meta %>% 
  filter(Day == "D7") %>% 
  select(`#SampleID`, Feeding_7) %>% 
  rename(Feeding = Feeding_7)
feeding30 = meta %>% 
  filter(Day == "D30") %>% 
  select(`#SampleID`, Feeding_30) %>% 
  rename(Feeding = Feeding_30)
feeding42 = meta %>% 
  filter(Day == "D42") %>% 
  select(`#SampleID`, Feeding_42) %>% 
  rename(Feeding = Feeding_42)
feeding = bind_rows(feeding3, feeding7, feeding30, feeding42) 



meta = meta %>% inner_join(feeding, by = "#SampleID") %>% 
  inner_join(confounder, by = c("ID" = "PatientID"))



taxa$`#OTU ID` = taxa_shortname(taxa$`#OTU ID`)

taxa = taxa %>% sjmisc::rotate_df(cn = T, rn = "#SampleID")

data = meta %>% inner_join(taxa, by = "#SampleID") %>% 
  mutate(ratio = Klebsiella/Escherichia)





fit = lmer(Bifidobacterium ~ Group * Time + Gestational_weeks + Birth_weight + Gender + Feeding + (1|ID), data = data)
res = summary(fit)
res = tidy(fit)

label = paste0("Group: P=", round(res$p.value[2], 3), "\nTime: P=", round(res$p.value[3], 3), "\nGroup * Time: P=", round(res$p.value[9], 3))
ggplot(data, aes(Time, Bifidobacterium, color = Group, fill = Group)) + 
  geom_smooth(alpha = 0.2) + 
  annotate(geom = "text", x = -Inf, y = Inf, vjust = 1.1, hjust = -0.1, label = label) + 
  theme_minimal() + 
  scale_color_manual(values = c("#7bc2ce", "#8562cc")) + 
  scale_fill_manual(values = c("#7bc2ce", "#8562cc")) + 
  scale_x_continuous(label = c("D3", "D7", "D30", "D42")) + 
  xlab("Day")


lmm_table = function(taxa){
  data2 = data %>% rename(feature = all_of(taxa))
  fit = lmer(feature ~ Group + Time + Gestational_weeks + Birth_weight + Gender + Feeding + (1|ID), data = data2)
  res = tidy(fit) %>% 
    filter(term %in% c("GroupVMT", "Time", "GroupVMT:Time")) %>% 
    mutate(y = {{taxa}}) %>% 
    select(y, everything())
  return(res)
  
}
taxa_plot = colnames(taxa)[-1]
all = map_dfr(taxa_plot, ~lmm_table(.x))
write_tsv(all, "LMM-Group + Time_adjust_confounder_result.tsv")






#taxa = "[Clostridium]"
lmm_func = function(taxa){
  data2 = data %>% rename(feature = all_of(taxa))
  fit = lmer(feature ~ Group * Time + Gestational_weeks + Birth_weight + Gender + Feeding + (1|ID), data = data2)
  res = tidy(fit)
  
  label = paste0("Group: P=", round(res$p.value[2], 3), "\nTime: P=", round(res$p.value[3], 3), "\nGroup * Time: P=", round(res$p.value[9], 3))
  p = ggplot(data2, aes(Time, feature, color = Group, fill = Group)) + 
    geom_point() +
    geom_smooth(alpha = 0.2) + 
    annotate(geom = "text", x = -Inf, y = Inf, vjust = 1.1, hjust = -0.1, label = label) + 
    theme_minimal() + 
    scale_color_manual(values = c("#7bc2ce", "#8562cc")) + 
    scale_fill_manual(values = c("#7bc2ce", "#8562cc")) + 
    scale_x_continuous(label = c("D3", "D7", "D30", "D42")) + 
    theme(axis.text = element_text(color = "black")) +
    labs(x = "Day", y = taxa)
  return(p)
}
a = lmm_func("Lactococcus")
a
taxa_plot = colnames(taxa)[-1]

Ps = map(taxa_plot, ~lmm_func(.x))
 pdf("interaction_scatterplot_Time_point_confounder.pdf", width = 5, height = 4)
 walk(Ps, print)
 dev.off()


ratio = lmm_func("ratio")
ratio
ratio + ylab("Klebsiella/Escherichia")
ggsave("KlebsiellaEscherichia_ratio_Timepoint_confounder.pdf", width = 5, height = 4)



lmm_func = function(taxa){
  data2 = data %>% rename(feature = all_of(taxa))
  fit = lmer(feature ~ Group * Time + Gestational_weeks + Birth_weight + Gender + Feeding + (1|ID), data = data2)
  res = tidy(fit)
  
  label = paste0("Group: P=", round(res$p.value[2], 3), "\nTime: P=", round(res$p.value[3], 3), "\nGroup * Time: P=", round(res$p.value[9], 3))
  p = ggplot(data2, aes(Time, feature, color = Group, fill = Group)) + 
    geom_smooth(alpha = 0.2) + 
    annotate(geom = "text", x = -Inf, y = Inf, vjust = 1.1, hjust = -0.1, label = label) + 
    theme_minimal() + 
    scale_color_manual(values = c("#7bc2ce", "#8562cc")) + 
    scale_fill_manual(values = c("#7bc2ce", "#8562cc")) + 
    scale_x_continuous(label = c("D3", "D7", "D30", "D42")) + 
    theme(axis.text = element_text(color = "black")) +
    labs(x = "Day", y = taxa)
  return(p)
}

taxa = "Bifidobacterium"

lmm_group_func = function(taxa){
  data2 = data %>% rename(feature = all_of(taxa))
  fit = lmer(feature ~ Group + Gestational_weeks + Birth_weight + Gender + Feeding + (1|ID), data = data2)
  res = tidy(fit)
  
  label = paste0("Group: P=", round(res$p.value[2], 3))
  p = ggplot(data2, aes(Time, feature, color = Group, fill = Group)) + 
    geom_smooth(alpha = 0.2) + 
    annotate(geom = "text", x = -Inf, y = Inf, vjust = 1.1, hjust = -0.1, label = label) + 
    theme_minimal() + 
    scale_color_manual(values = c("#7bc2ce", "#8562cc")) + 
    scale_fill_manual(values = c("#7bc2ce", "#8562cc")) + 
    scale_x_continuous(label = c("D3", "D7", "D30", "D42")) + 
    theme(axis.text = element_text(color = "black")) +
    labs(x = "Day", y = taxa)
  return(p)
}


lmm_group_func("Bifidobacterium")

Ps = map(taxa_plot, ~lmm_group_func(.x))
pdf("Group_lineplot_Time_point_confounder.pdf", width = 5, height = 4)
walk(Ps, print)
dev.off()


lmm_time_func = function(taxa){
  data2 = data %>% rename(feature = all_of(taxa))
  fit = lmer(feature ~ Time + Gestational_weeks + Birth_weight + Gender + Feeding + (1|ID), data = data2)
  res = tidy(fit)
  
  label = paste0("Time: P=", round(res$p.value[2], 3))
  p = ggplot(data2, aes(Time, feature, color = Group, fill = Group)) + 
    geom_smooth(alpha = 0.2) + 
    annotate(geom = "text", x = -Inf, y = Inf, vjust = 1.1, hjust = -0.1, label = label) + 
    theme_minimal() + 
    scale_color_manual(values = c("#7bc2ce", "#8562cc")) + 
    scale_fill_manual(values = c("#7bc2ce", "#8562cc")) + 
    scale_x_continuous(label = c("D3", "D7", "D30", "D42")) + 
    theme(axis.text = element_text(color = "black")) +
    labs(x = "Day", y = taxa)
  return(p)
}
Ps = map(taxa_plot, ~lmm_time_func(.x))
pdf("Time_lineplot__confounder.pdf", width = 5, height = 4)
walk(Ps, print)
dev.off()

