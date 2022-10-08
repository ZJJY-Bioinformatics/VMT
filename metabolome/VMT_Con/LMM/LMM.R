library(tidyverse)
library(lmerTest)
library(broom.mixed)
source("/data/scripts/utility.R")

meta = read_tsv("../../VMT_Con_metadata_infant_stool_metabolome.txt")%>% 
  select(`#SampleID`, ID, Group, Day) %>% 
  mutate(Day = fct_relevel(Day, c("D3", "D7", "D30", "D42")), 
         Time = as.numeric(Day),
         Day_time = as.numeric(str_extract(Day, "\\d+")))


taxa = read_tsv("../../../../../metabolome/metabolome-tsfm.txt") %>% 
  column_to_rownames(var = "Sample") %>% 
  sjmisc::rotate_df() 





taxa_ap = taxa %>% mutate(across(where(is.numeric), ~if_else(.x>0,1,0)))
idx_0.1 = rowSums(taxa_ap[, -1]) > 37
taxa = taxa[idx_0.1, ]



taxa = taxa %>% sjmisc::rotate_df(rn = "#SampleID")

data = meta %>% inner_join(taxa, by = "#SampleID")

fit = lmer(Bifidobacterium ~ Group*Time + (1|ID), data = data)
res = summary(fit)
res = tidy(fit)


lmm_table = function(taxa){
  data2 = data %>% rename(feature = all_of(taxa))
  fit = lmer(feature ~ Group*Time + (1|ID), data = data2)
  res = tidy(fit) %>% 
    filter(term %in% c("GroupVMT", "Time", "GroupVMT:Time")) %>% 
    mutate(y = {{taxa}}) %>% 
    select(y, everything())
  return(res)

}
taxa_plot = colnames(taxa)[-1]
all = map_dfr(taxa_plot, ~lmm_table(.x))
write_tsv(all, "LMM-GroupTime_interaction_result.tsv")




label = paste0("Group: P=", round(res$p.value[2], 3), "\nTime: P=", round(res$p.value[3], 3), "\nGroup * Time: P=", round(res$p.value[4], 3))
ggplot(data, aes(Time, Bifidobacterium, color = Group, fill = Group)) + 
  geom_smooth(alpha = 0.2) + 
  annotate(geom = "text", x = -Inf, y = Inf, vjust = 1.1, hjust = -0.1, label = label) + 
  theme_minimal() + 
  scale_color_manual(values = c("#7bc2ce", "#8562cc")) + 
  scale_fill_manual(values = c("#7bc2ce", "#8562cc")) + 
  scale_x_continuous(label = c("D3", "D7", "D30", "D42")) + 
  xlab("Day")

#taxa = "[Clostridium]"
lmm_func = function(taxa){
  data2 = data %>% rename(feature = all_of(taxa))
  fit = lmer(feature ~ Group * Time + (1|ID), data = data2)
  res = tidy(fit)
  
  label = paste0("Group: P=", round(res$p.value[2], 3), "\nTime: P=", round(res$p.value[3], 3), "\nGroup * Time: P=", round(res$p.value[4], 3))
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
a = lmm_func(taxa_plot[1])
a
taxa_plot = colnames(taxa)[-1]

Ps = map(taxa_plot, ~lmm_func(.x))
pdf("interaction_lineplot_Day.pdf", width = 5, height = 4)
walk(Ps, print)
dev.off()


lmm_func("ratio") + 
  ylab("Klebsiella/Escherichia")
ggsave("KlebsiellaEscherichia_ratio_Timepoint.pdf", width = 5, height = 4)





lmm_group_func = function(taxa){
  data2 = data %>% rename(feature = all_of(taxa))
  fit = lmer(feature ~ Group + (1|ID), data = data2)
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
pdf("Group_lineplot_timepoint.pdf", width = 5, height = 4)
walk(Ps, print)
dev.off()



lmm_time_func = function(taxa){
  data2 = data %>% rename(feature = all_of(taxa))
  fit = lmer(feature ~ Time + (1|ID), data = data2)
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
lmm_time_func("Bifidobacterium")


Ps = map(taxa_plot, ~lmm_time_func(.x))
pdf("time_lineplot.pdf", width = 5, height = 4)
walk(Ps, print)
dev.off()


res = read_tsv("LMM-GroupTime_interaction_adjust_confounder_result.tsv") %>% 
  filter(term == "GroupVMT:Time") %>% 
  mutate(FDR = p.adjust(p.value, method = "fdr"))
write_tsv(res, "LMM-GroupTime_interaction_term_adjust_confounder_FDR.tsv")

res = read_tsv("LMM-GroupTime_interaction_result.tsv") %>% 
  filter(term == "GroupVMT:Time") %>% 
  mutate(FDR = p.adjust(p.value, method = "fdr"))
write_tsv(res, "LMM-GroupTime_interaction_term_FDR.tsv")
