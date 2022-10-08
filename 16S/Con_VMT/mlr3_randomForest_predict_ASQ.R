library(tidyverse)
library(mlr3)
library(mlr3verse)
library(magrittr)
confounders = readxl::read_xlsx("拟调整混杂因素(1).xlsx") %>% 
  select(PatientID, Gestational_weeks, Birth_weight)

map_metabolome = read_tsv("../../../../Metadata_infant_stool_metabolome(ASQ).txt") 
map_microbiome = read_tsv("../../../../Metadata_infant_stool_microbiome(ASQ).txt") 
sample = intersect(map_metabolome$`#SampleID`, map_microbiome$`#SampleID`)
map_metabolome %<>% filter(`#SampleID` %in% sample)
map_microbiome %<>% filter(`#SampleID` %in% sample)

map = map_microbiome %>% select(`#SampleID`, M6_total, ID) %>% 
  filter(!is.na(M6_total))

metabolite = read_tsv("metabolome-tsfm.txt") %>% 
  filter(Sample %in% sample) %>% 
  sjmisc::rotate_df(cn = T, rn = "metabolites") %>% 
  mutate(alias = paste0("metabolites", row_number()))%>% 
  select(-metabolites) %>% 
  select(alias, everything()) %>% 
  sjmisc::rotate_df(rn = "#SampleID", cn = T)


L6 = read_tsv("infant_stool_e10k_L6.txt", skip = 1)

taxa_ap = L6 %>% mutate(across(where(is.numeric), ~if_else(.x>0,1,0)))
idx_0.1 = rowSums(taxa_ap[, -1]) > 23
taxa = L6[idx_0.1, ]


L6 = L6 %>% 
   select(`#OTU ID`, all_of(sample)) %>% 
   mutate(alias = paste0("genus", row_number())) %>% 
   select(-`#OTU ID`) %>% 
   select(alias, everything()) %>% 
   sjmisc::rotate_df(rn = "#SampleID", cn = T)

data = map %>% inner_join(metabolite, by = "#SampleID") %>% 
  inner_join(L6, by = "#SampleID") %>% 
  inner_join(confounders, by = c("ID" = "PatientID")) %>% 
  rename(SampleID = "#SampleID")

VD_task = data %>% filter(Group == "VD") %>% 
  select(-ID ) %>% 
  as_task_regr(target = "M6_total")

task = data %>% 
  filter(Group != "VD") %>% 
  select(-ID ) %>% 
  as_task_regr(target = "M6_total")
task$set_col_roles("SampleID", roles = "name")
learner <- lrn("regr.ranger", num.threads = 10, importance = "impurity")
search_space <- ps(
  mtry = p_int(lower = 10, upper = length(task$feature_names))
)

resampling <- rsmp("cv", folds = 10)
measure <- msr("regr.rsq")
none <- trm("none")
tuner <- tnr("grid_search", resolution=10)
at <- AutoTuner$new(
  learner = learner,
  resampling = resampling,
  measure = measure,
  search_space = search_space,
  terminator = none,
  tuner = tuner
)

lgr::get_logger("bbotk")$set_threshold("warn")
lgr::get_logger("mlr3")$set_threshold("warn")

set.seed(2022)
at$train(task)
best_nr <- at$archive$best()$batch_nr
best_pred <- at$archive$predictions(best_nr)

best_tuned_pred <- map_dfr(best_pred, ~ bind_rows(as.data.table(.x)))
write_tsv(best_tuned_pred, "代谢物和菌预测六月龄总分10-fold_CV.tsv")
spearman = cor.test(best_tuned_pred$truth, best_tuned_pred$response, method = "spearman")
spearman$estimate
spearman$p.value
fit = lm(truth ~ response, data = best_tuned_pred)
summary(fit)
ggplot(best_tuned_pred, aes(truth, response)) + 
  geom_point(shape = 21, fill = "#7bc2ce", size = 2) + 
  geom_smooth(method = lm, color ="#7bc2ce", size = 2) + 
  theme_minimal()+ 
  theme(aspect.ratio = 1) + 
  labs(x = "Actual values", y = "Predicted values") +
  annotate(geom = "text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = "Spearman r = 0.35 \n P<0.001")
ggsave("代谢物和菌预测六月龄总分_cv.pdf", width = 3, height = 3)
res = c(Day, score, spearman$estimate[[1]], spearman$p.value)


# 使用调优后的超参数建模
learner$param_set$values <- at$tuning_result$learner_param_vals[[1]]

learner$train(task)

# 特征重要性
importance <-  as.data.table(learner$importance(), keep.rownames = TRUE)

metabolite = read_tsv("../../../../../metabolome/metabolome-tsfm.txt") %>% 
  filter(Sample %in% sample) %>% 
  sjmisc::rotate_df(cn = T, rn = "metabolites") %>% 
  mutate(alias = paste0("metabolites", row_number())) %>% 
  select(metabolites, alias) %>% 
  rename(`#OTU ID` = metabolites)


# 
# L6_alias = taxa %>% 
#   filter(!str_detect(`#OTU ID`, ".*g__$")) %>% 
#   filter(!str_detect(`#OTU ID`, ".*Other$")) %>% 
#   mutate(alias = paste0("genus", row_number())) %>% 
#   select(`#OTU ID`, alias) 
L6 = read_tsv("infant_stool_e10k_L6.txt", skip = 1)
taxa_ap = L6 %>% mutate(across(where(is.numeric), ~if_else(.x>0,1,0)))
idx_0.1 = rowSums(taxa_ap[, -1]) > 23
taxa = L6[idx_0.1, ]
L6_alias = L6 %>% 
   select(`#OTU ID`, all_of(sample)) %>% 
   mutate(alias = paste0("genus", row_number())) %>% 
   select(`#OTU ID`, alias) 

alias = metabolite %>% bind_rows(L6_alias)

imp = importance %>% as.data.frame() %>% left_join(alias, by = c("V1" = "alias")) %>% 
  #select(`#OTU ID`, V2) %>% 
  rename(importance = V2) 
write.table(imp, "importance.txt", sep = "\t", quote = F, row.names = F)

#colnames(importance) = c("Feature", "Importance")
imp_plt <- ggplot(head(imp, 10), aes(x = reorder(`#OTU ID`, importance), y = importance)) +
  geom_col(fill = "#1f77b4", width = 0.7, color = "black") + coord_flip() + xlab("") + 
  theme_light()
imp_plt
ggsave("top10_importance.pdf", width = 10, height = 6)



