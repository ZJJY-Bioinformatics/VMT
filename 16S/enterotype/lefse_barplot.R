library(tidyverse)

lefse = read_tsv("tmp/stat.res") %>% 
  filter(!is.na(g)) %>% 
  filter(LDA > 4) %>% 
  mutate(LDA = case_when(enriched_in == "ET1" ~ LDA, 
                         enriched_in == "ET2" ~ -LDA)) %>% 
  arrange(LDA) %>% 
  mutate(g = fct_relevel(g, g))

ggplot(lefse, aes(g, LDA, fill = enriched_in)) +
  geom_hline(yintercept =  c(-3, -1.5, 0, 1.5, 3), linetype = 2, color = "grey60") + 
  geom_bar(stat = "identity", color = "black", alpha = 0.8) +
  coord_flip() +
  scale_fill_manual(values = c('#ff7f0e', '#1f77b4')) + 
  theme_test() + 
  theme(axis.title.y = element_blank(),
        legend.title = element_blank(), 
        axis.text = element_text(color="black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        panel.border = element_blank(), 
        legend.position = "top") + 
  geom_text(data = lefse[1:8, ], aes(x=g, y = 0.1, label =g), hjust =0, fontface = "italic") +
  geom_text(data = lefse[9:17, ], aes(x=g, y = -0.1, label =g), hjust =1, fontface = "italic") + 
  labs(y = "LDA Score (log 10)", fill = "enriched_in") + 
  scale_y_continuous(breaks = c(-3, -1.5, 0, 1.5, 3))

ggsave("ET_lefse.pdf", width = 4, height = 3.5)


