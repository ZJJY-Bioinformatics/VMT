source("C:/data/scripts/utility.R")

library(tidyverse)

L6 = read.table("infant_stool_e10k_L6.txt", comment.char = "", skip = 1, header = T, sep = "\t", check.names = F) %>% 
  column_to_rownames(var = "#OTU ID") %>% 
  sjmisc::rotate_df(rn = "#SampleID")
meta = read_tsv("../../../Metadata_infant_stool_microbiome.txt") %>% 
  select(`#SampleID`, Group, Day) %>% 
  #filter(Group != "VD")

L6 = L6 %>% filter(`#SampleID` %in% meta$`#SampleID`) %>% 
  column_to_rownames(var = "#SampleID")

# 各物种丰度排序
top <- sort(colSums(L6), decreasing=TRUE)

# 选择丰度前n的物种进行绘制图形
n = 13
need <- names(top)[1:n]
L6_RA <- L6 %>% 
  select(all_of(need))
colnames(L6_RA) = taxa_shortname(colnames(L6_RA))
Others <- 1-rowSums(L6_RA)


# topN物种+others丰度表
genus_abdc = cbind(L6_RA, Others) %>% 
  rownames_to_column(var = "#SampleID")


data = meta %>% 
  inner_join(genus_abdc, by = "#SampleID") %>% 
  group_by(Day, Group) %>% 
  summarise(across(where(is.numeric), mean)) 






data = data %>% 
  ungroup() %>% 
  pivot_longer(-1:-2, names_to = "otus", values_to = "value")%>% 
  mutate(Day = fct_relevel(Day, c("D3", "D7", "D30", "D42")),
         otus = fct_relevel(otus, rev(c(colnames(L6_RA), "Others"))), 
         Group = fct_relevel(Group, c("VD","Con", "VMT")))


colors = c("#80B1D3","#B3DE69","#FFFFB3","#8DD3C7","#4daf4a",
           "#377eb8","#BEBADA","#FB8072","#FDB462","#FCCDE5",
           "#BC80BD","#CCEBC5","#FFED6F","#CD4F39","#BC41A4",
           "#4F94CD","#E41A1C","#00CD66","#CD3278","#CD8A96",
           "#00C5CD","#CDCD00","#CD85CD","#CD853F","#8B5A2B",
           "#5CACEE","#EE5C42","#00EE76","#EE4A8C","#EED8AE",
           "#00E5EE","#EEEE00","#EED2EE","#EE9A49","#E41A1C",
           "#377EB8","#FF6A6A","#87CEFA","#6E8B3D","#FFEBCD",
           "#B2DFEE")

mycolors <- c("#D9D9D9",rev(colors[1:length(unique(data$otus))-1]))


p1 <- ggplot(data,aes(x=Day,y=value,fill=otus))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = mycolors)+
  theme(legend.background = element_rect(size = 0.1),
        legend.key.size= unit(5,"mm"),
        legend.text = element_text(size = 8))+
  theme_classic()+
  ylab(NULL)+
  xlab(NULL)+
  theme(legend.title=element_blank(),
        axis.title.x=element_text(vjust=3.5),
        axis.text.x=element_text(size=8,color="black"),
        axis.line=element_line(colour ="grey"),
        panel.background = element_rect(fill = "transparent",colour = "grey", size = 1, linetype = 1),
        panel.border=element_rect(fill="transparent",color="grey"), 
        aspect.ratio = 1.2)+ 
  facet_wrap(.~Group)
p1



ggsave("各组随时间物种组成.pdf", width = 8, height = 6)




data = meta %>% 
  inner_join(genus_abdc, by = "#SampleID") %>% 
  group_by(Group) %>% 
  summarise(across(where(is.numeric), mean)) %>% 
  ungroup() %>% 
  pivot_longer(-1, names_to = "otus", values_to = "value")%>% 
  filter(Group != "VD") %>% 
  mutate(otus = fct_relevel(otus, rev(c(colnames(L6_RA), "Others"))))

mycolors <- c("#D9D9D9",rev(colors[1:length(unique(data$otus))-1]))
ggplot(data,aes(x=Group,y=value,fill=otus))+
  geom_bar(stat = "identity", width = 0.6)+
  scale_fill_manual(values = mycolors)+

  theme_minimal()+
  ylab("Relative abundance")+
  xlab(NULL)+
  theme(legend.title=element_blank(),
        axis.text=element_text(color="black"),
        panel.border=element_rect(fill="transparent",color="grey"), 
        aspect.ratio = 1.2,
        legend.text = element_text(face = "italic")) + 
  scale_y_continuous(breaks = seq(0,1, 0.2))
ggsave("fig1g-两组所有样本物种组成图.pdf", width = 5, height = 5)

#################L2

L6 = read.table("infant_stool_e10k_L2.txt", comment.char = "", skip = 1, header = T, sep = "\t", check.names = F) %>% 
  column_to_rownames(var = "#OTU ID") %>% 
  sjmisc::rotate_df()
meta = read_tsv("../../../Metadata_infant_stool_microbiome.txt") %>% 
  select(`#SampleID`, Group, Day)

# 各物种丰度排序
top <- sort(colSums(L6), decreasing=TRUE)

# 选择丰度前n的物种进行绘制图形
n = 4
need <- names(top)[1:n]

L6_RA <- L6 %>% 
  select(all_of(need))
colnames(L6_RA) = taxa_shortname(colnames(L6_RA))
Others <- 1-rowSums(L6_RA)


# topN物种+others丰度表
genus_abdc = cbind(L6_RA, Others) %>% 
  rownames_to_column(var = "#SampleID")


data = meta %>% 
  inner_join(genus_abdc, by = "#SampleID") %>% 
  group_by(Day, Group) %>% 
  summarise(across(where(is.numeric), mean)) 




data = data %>% 
  ungroup() %>% 
  pivot_longer(-1:-2, names_to = "otus", values_to = "value")%>% 
  mutate(Day = fct_relevel(Day, c("D3", "D7", "D30", "D42")),
         otus = fct_relevel(otus, rev(c(colnames(L6_RA), "Others"))), 
         Group = fct_relevel(Group, c("VD","Con", "VMT")))


colors = c("#80B1D3","#B3DE69","#FFFFB3","#8DD3C7","#4daf4a",
           "#377eb8","#BEBADA","#FB8072","#FDB462","#FCCDE5",
           "#BC80BD","#CCEBC5","#FFED6F","#CD4F39","#BC41A4",
           "#4F94CD","#E41A1C","#00CD66","#CD3278","#CD8A96",
           "#00C5CD","#CDCD00","#CD85CD","#CD853F","#8B5A2B",
           "#5CACEE","#EE5C42","#00EE76","#EE4A8C","#EED8AE",
           "#00E5EE","#EEEE00","#EED2EE","#EE9A49","#E41A1C",
           "#377EB8","#FF6A6A","#87CEFA","#6E8B3D","#FFEBCD",
           "#B2DFEE")

mycolors <- c("#D9D9D9",rev(colors[1:length(unique(data$otus))-1]))


p1 <- ggplot(data,aes(x=Day,y=value,fill=otus))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = mycolors)+
  theme(legend.background = element_rect(size = 0.1),
        legend.key.size= unit(5,"mm"),
        legend.text = element_text(size = 8))+
  theme_classic()+
  ylab(NULL)+
  xlab(NULL)+
  theme(legend.title=element_blank(),
        axis.title.x=element_text(vjust=3.5),
        axis.text.x=element_text(size=8,color="black"),
        axis.line=element_line(colour ="grey"),
        panel.background = element_rect(fill = "transparent",colour = "grey", size = 1, linetype = 1),
        panel.border=element_rect(fill="transparent",color="grey"), 
        aspect.ratio = 1.2)+ 
  facet_wrap(.~Group)
p1



ggsave("各组随时间物种组成_L2.pdf", width = 8, height = 6)




