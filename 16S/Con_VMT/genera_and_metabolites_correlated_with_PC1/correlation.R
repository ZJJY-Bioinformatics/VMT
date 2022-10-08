library(tidyverse)
library(psych)
library(ggrepel)
source("/data/scripts/utility.R")

pc = read.pc("bray_curtis_pc.txt")
PC = pc$pc_data %>% as.data.frame() %>% select(PC1:PC5) %>% 
  rownames_to_column(var = "#SampleID")
PC = PC %>% select(`#SampleID`, PC1) %>% 
  column_to_rownames(var = "#SampleID")

taxa = read_tsv("../../taxa/infant_stool_e10k_L6.txt", skip = 1) %>% 
  filter(!str_detect(`#OTU ID`, ".*Other")) %>% 
  filter(!str_detect(`#OTU ID`, ".*g__$")) 

taxa_ap = taxa %>% mutate(across(where(is.numeric), ~if_else(.x>0,1,0)))
idx_0.1 = rowSums(taxa_ap[, -1]) > 23
taxa = taxa[idx_0.1, ]

taxa$`#OTU ID` = taxa_shortname(taxa$`#OTU ID`)


taxa = taxa %>% sjmisc::rotate_df(cn = T)

sample = intersect(rownames(taxa), rownames(PC))

taxa = taxa[sample, ]
PC = PC[sample,]

res = corr.test(taxa, PC, method = "spearman", adjust = "fdr")
r = res$r %>% as.data.frame() %>% rownames_to_column(var = "feature")
p = res$p.adj
fdr = res$p

genera = bind_cols(r, p, fdr) %>% 
  rename(spearman_r = V1,
         p.value = `...3`,
         FDR = `...4`)
#write_tsv(genera, "与PC1相关的菌—L6.tsv")
data = read_tsv("与PC1相关的菌—L6.tsv") %>% 
  mutate(pn = case_when(spearman_r>0 & FDR < 0.05 ~ "Positive",
                        spearman_r < 0 & FDR <0.05 ~ "Negative", 
                        TRUE ~ "ns"),
         label = ifelse(pn != "ns", feature, NA))
color = c("#4194C4", "#d9d9d9", "#FC4D25")
ggplot(data, aes(spearman_r, -log10(FDR), color=pn)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") + 
  geom_point() + 
  geom_text_repel(aes(label = label), max.overlaps = Inf, max.iter = 20000, show.legend = F, box.padding = 0.5, ylim = c(-Inf, Inf)) +
  theme_minimal()+
  scale_color_manual(values = color) +
  scale_y_continuous(expand = expansion(mult = 0.1))
ggsave("与PC1相关的genusFDR0.05.pdf", width = 3.8, height = 3)





#=========代谢物===============
taxa = read_tsv("../../../../../metabolome/metabolome-tsfm.txt") %>% 
  column_to_rownames(var = "Sample") %>% 
  sjmisc::rotate_df() %>% 
  select(!(ends_with("BSD1") | ends_with("BSD2")))


taxa_ap = taxa %>% mutate(across(where(is.numeric), ~if_else(.x>0,1,0)))
idx_0.1 = rowSums(taxa_ap[, -1]) > 37
taxa = taxa[idx_0.1, ]

taxa = taxa %>% sjmisc::rotate_df()

PC = pc$pc_data %>% as.data.frame() %>% select(PC1:PC5) %>% 
  rownames_to_column(var = "#SampleID")
PC = PC %>% select(`#SampleID`, PC1) %>% 
  column_to_rownames(var = "#SampleID")

sample = intersect(rownames(taxa), rownames(PC))

taxa = taxa[sample, ]
PC = PC[sample,]


res = corr.test(taxa, PC, method = "spearman", adjust = "fdr")
r = res$r %>% as.data.frame() %>% rownames_to_column(var = "feature")
p = res$p.adj
fdr = res$p

metabolites = bind_cols(r, p, fdr) %>% 
  rename(spearman_r = V1,
         p.value = `...3`,
         FDR = `...4`)
#write_tsv(metabolites, "与PC1相关的代谢物.tsv")
library(ggrepel)
data = read_tsv("与PC1相关的代谢物.tsv") %>% 
  mutate(pn = case_when(spearman_r>0 & FDR < 0.25 ~ "Positive",
                        spearman_r < 0 & FDR <0.25 ~ "Negative", 
                        TRUE ~ "ns"))
label = data %>% filter(spearman_r>0.4|spearman_r < -0.3)
color = c("#4194C4", "#d9d9d9", "#FC4D25")
ggplot(data, aes(spearman_r, -log10(FDR), color=pn)) +
  geom_hline(yintercept = -log10(0.25), linetype = "dashed") + 
  geom_point() + 
  geom_text_repel(data = label, mapping = aes(spearman_r, -log10(FDR), label = feature), show.legend = F) +
  theme_minimal()+
  scale_color_manual(values = color) + 
  theme(legend.title = element_blank())
ggsave("与PC1相关的代谢物FDR0.25.pdf", width = 3.8, height = 3)

