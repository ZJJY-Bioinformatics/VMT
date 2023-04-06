library(tidyverse)
library(psych)
source("/data/scripts/utility.R")


meta = read_tsv("../../VMT_Con_metadata_infant_stool_microbiome.txt") 

pc = read.pc("bdiv/bray_curtis_pc.txt")
PC = pc$pc_data %>% as.data.frame() %>% select(PC1:PC5) %>% 
  rownames_to_column(var = "#SampleID")
proportion = pc$pc_percent

PC_plot = meta %>% select(`#SampleID`, Day, Time_point, Group) %>% 
  inner_join(PC, by = "#SampleID") %>% 
  mutate(Day = fct_relevel(Day, c("D3", "D7", "D30", "D42")))


genus = read_tsv("infant_stool_e10k_L6_s46.biom.txt", skip = 1) %>% 
  sjmisc::rotate_df(cn = T) 

PC = PC %>% select(`#SampleID`, PC1) %>% 
  column_to_rownames(var = "#SampleID")

sample = intersect(rownames(PC), rownames(genus))

PC = PC[sample, ]
genus = genus[sample, ]

res = corr.test(genus, PC, method = "spearman", adjust = "fdr")

r  = res$r %>% as.data.frame() %>% rownames_to_column(var = "genus")
p = res$p.adj
fdr = res$p

all = bind_cols(r, p, fdr) %>% 
  rename(spearman_r = V1, p.value = `...3`, FDR = `...4`)
write_tsv(all, "bray_curtis_PC1_and_genus_correlation.tsv")
