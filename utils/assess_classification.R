## Assess model classification ##

library(tidyverse)
library(magrittr)
library(broom)
library(emmeans)
library(tidymodels)
library(shadowtext)
library(patchwork)
library(gt)
library(vegan)
library(vip)

plot_confMat <- function(x, font_size = 8){
  x %>%
    tidy %>%
    separate(name, into = c('tmp', 'row', 'col')) %>%
    mutate(across(c(row, col), as.integer),
           value2 = formatC(value, format="f", big.mark=",", digits=0),
           value2 = str_c('Prediction = ', c('Sand', 'Reef', 'Sand', 'Reef'), '\n',
                          'Truth = ', c('Sand', 'Sand', 'Reef', 'Reef'), '\n', value2)) %>%
    ggplot(aes(x = col, y = -row, fill = value)) +
    geom_tile(show.legend = FALSE) +
    geom_shadowtext(aes(label = value2), colour = 'white', fontface = 'bold', size = font_size) +
    theme_void()
}

veganCovEllipse<-function (x, se = TRUE, conf = 0.95, npoints = 100) 
{
  #X is a dataframe of 2 coordinates
  
  covariance_mat <- cov.wt(x, wt=rep(1/nrow(x), nrow(x)))
  
  cov <- covariance_mat$cov
  
  if(se){cov <- cov * sum(covariance_mat$wt^2)}
  
  center <- covariance_mat$center
  
  scale <- sqrt(qchisq(conf, 2))
  
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov))) %>% as_tibble()
}


#### Read in Relevant Files ####
metrics_to_keep <- c('accuracy', 'roc_auc', 'j_index', 'sens', 'spec', 'mcnemar')
models_to_remove <- c('rbf_svm', 'cart')

overall_stats <- read_rds('../../Results/Habitat_Classification/All_topParam_model_assessment_small.rds') %>%
  filter(!model_name %in% models_to_remove) %>%
  select(model_name, best_params, confusion, contains(metrics_to_keep), -bal_accuracy)

site_stats <- read_rds('../../Results/Habitat_Classification/All_topParam_model_bySite.rds') %>%
  filter(!model_name %in% models_to_remove) %>%
  select(model_name, Site, confusion, contains(metrics_to_keep), -bal_accuracy)

#### Summarize training data ####


#### Plot Overall Stats highlighting best model for each metric ####
overall_stats %>%
  mutate(log_mcnemar_statistic = log(mcnemar_statistic, base = 10)) %>%
  select(-confusion, -best_params, -starts_with('mcnemar')) %>%
  pivot_longer(cols = -model_name) %>%
  group_by(name) %>%
  mutate(top_model = case_when(str_detect(name, 'mcnemar_statistic', negate = TRUE) ~ value == max(value, na.rm = TRUE),
                               TRUE ~ value == min(value, na.rm = TRUE))) %>%
  
  ungroup %>%
  ggplot(aes(x = model_name, y = value, colour = top_model)) +
  geom_text(aes(label = model_name), position = position_dodge(width = 0.2), show.legend = FALSE) +
  geom_hline(data = tibble(name = 'log_mcnemar_statistic', value = log(qchisq(0.05, 1, lower.tail = FALSE), base = 10)),
             aes(yintercept = value), linetype = 'dashed') +
  facet_wrap(~ name, scales = 'free_y')

##
overall_stats %>%
  select(model_name, confusion) %>%
  mutate(confusion = map(confusion, plot_confMat, font_size = 3),
         confusion = map2(confusion, model_name, ~.x + labs(title = .y))) %>%
  pull(confusion) %>%
  wrap_plots()

overall_stats %>%
  filter(model_name == 'c5') %>%
  select(model_name, confusion) %>%
  mutate(confusion = map(confusion, plot_confMat, font_size = 8),
         confusion = map2(confusion, model_name, ~.x + labs(title = .y))) %>%
  pull(confusion) %>%
  wrap_plots()

## Cluster models to see which overall looks 'best' ##
overall_nmds_mat <- overall_stats %>%
  # mutate(log_mcnemar_statistic = log(mcnemar_statistic, base = 10)) %>%
  select(-confusion, -best_params, -starts_with('mcnemar_p')) %>%
  recipe(~ ., data = .) %>%
  step_YeoJohnson(all_numeric()) %>%
  step_normalize(all_numeric()) %>%
  prep %>%
  juice %>%
  select(-model_name) %>%
  as.matrix() %>%
  set_rownames(overall_stats$model_name) 

overall_nmds <- vegan::metaMDS(overall_nmds_mat, distance = 'euclidian', k = 2, autotransform = FALSE) 

overall_arrows <- envfit(overall_nmds, overall_nmds_mat) %$%
  vectors$arrows %>%
  as_tibble(rownames = 'metric') %>%
  rename_with(~str_remove(., '^N'))
  
overall_stats <-overall_nmds %$%
  points %>%
  as_tibble(rownames = 'model_name') %>%
  right_join(overall_stats)

overall_stats %>%
  ggplot(aes(x = MDS1, y = MDS2)) +
  geom_segment(data = overall_arrows, aes(xend = 0, yend = 0), arrow = arrow(ends = 'first')) +
  geom_text(data = overall_arrows, aes(label = metric)) +
  
  geom_text(aes(label = model_name), show.legend = FALSE)

the_pca <- overall_stats %>%
  # mutate(log_mcnemar_statistic = log(mcnemar_statistic, base = 10)) %>%
  select(-confusion, -best_params, -starts_with('mcnemar_p'), -starts_with('MDS')) %>%
  recipe(~ ., data = .) %>%
  update_role(model_name, new_role = "id") %>%
  step_YeoJohnson(all_numeric()) %>%
  step_normalize(all_numeric()) %>%
  step_pca(all_predictors(), num_comp = 6) %>%
  prep

library(tidytext)

the_pca %>%
  tidy(3) %>%
  group_by(component) %>%
  ungroup() %>%
  mutate(terms = reorder_within(terms, abs(value), component)) %>%
  ggplot(aes(abs(value), terms, fill = value > 0)) +
  geom_col() +
  facet_wrap(~component, scales = "free_y") +
  scale_y_reordered() +
  labs(
    x = "Absolute value of contribution",
    y = NULL, fill = "Positive?"
  )

the_pca$steps[[3]]$res %>%
  tidy('eigenvalues') %>%
  ggplot(aes(PC, percent)) +
  geom_col(fill = "#56B4E9", alpha = 0.8)

the_segments <- the_pca$steps[[3]]$res %>%
  tidy(matrix = 'rotation') %>%
  pivot_wider(names_from = "PC", 
              names_prefix = "PC",
              values_from = "value") 

the_pca %>%
  juice %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_segment(data = the_segments, xend = 0, yend = 0) +
  geom_text(data = the_segments, aes(label = column)) +
  geom_text(aes(label = model_name))


## Make table of overall metrics
overall_stats %>%
  arrange(MDS1) %>%
  select(-best_params, -confusion, -starts_with('MDS')) %>%
  mutate(mcnemar_p.value = p.adjust(mcnemar_p.value, 'holm')) %>% 
  select(model_name, j_index, everything()) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  write_csv('tmp_table.csv')

#### Use sites as replicates to find best model for each metric ####
source('~/R/R Functions/tukeyNonAdditivityTest.R')

model_metric_stats <- site_stats %>%
  select(-confusion, -mcnemar_p.value) %>%
  pivot_longer(cols = -c(model_name, Site)) %>%
  nest(data = -c(name)) %>%
  mutate(test = map(data, ~aov(value ~ model_name + Site, data = .x))) %>%

  rowwise() %>%
  mutate(non_additivity = list(tukey.nonadditivity.test(test) %>% tidy %>% select(p.value) %>% rename(interaction_p = p.value)),
         sw_test = list(shapiro.test(resid(test)) %>% tidy %>% select(p.value) %>% rename(sw_p = p.value))) %>%
  unnest(c(sw_test, non_additivity)) %>%
  
  #add something here to log transform??
  #add something here to deal with site interaction??
  rowwise %>%
  mutate(test_post_model = list(emmeans(test, pairwise ~ model_name)),
         test_post_site = list(emmeans(test, pairwise ~ Site)),
         test = list(tidy(test)))

top_metric_models <- model_metric_stats %>%
  unnest(test) %>%
  filter(term == 'model_name') %>%
  select(name, test_post_model) %>%
  mutate(means = map(test_post_model, ~.x$emmeans %>% as_tibble),
         pairwise = map(test_post_model, ~.x$contrasts %>% as_tibble)) %>%
  select(-test_post_model) %>%
  mutate(best_model = map2_chr(means, name, ~.x %>%
                                 {if(str_detect(.y, 'mcnemar')) filter(., emmean == min(emmean)) else filter(., emmean == max(emmean))} %>%
                                pull(model_name) %>%
                                as.character()),
         top_models = map2(pairwise, best_model, ~ .x %>%
                             filter(str_detect(contrast, .y)) %>%
                             filter(p.value >= 0.1) %>%
                             pull(contrast) %>%
                             str_split(' - ') %>%
                             unlist %>%
                             unique)) %>%
  select(name, best_model, top_models) %>%
  unnest(top_models)


site_stats %>%
  select(-confusion, -mcnemar_p.value) %>%
  rename_with(~str_replace(., '_', '.')) %>%
  group_by(model.name) %>%
  summarise(across(-Site, list(mean = ~mean(., na.rm = TRUE), 
                               sd = ~sd(., na.rm = TRUE), 
                               n = ~sum(!is.na(.)))),
            .groups = 'drop') %>%
  pivot_longer(cols = -model.name,
               names_to = c('name', '.value'),
               names_pattern = '(.*)_(.*)') %>%
  rename(model_name = model.name) %>%
  mutate(name = str_replace(name, '\\.', '_')) %>%
  mutate(se = sd/sqrt(n),
         lwr = mean - se,
         upr = mean + se) %>%
  mutate(across(c(mean, lwr, upr), ~case_when(str_detect(name, 'mcnemar') ~ log(., base = 10),
                                              TRUE ~ .))) %>%
  left_join(select(top_metric_models, name, best_model) %>%
              distinct %>%
              mutate(is_best = TRUE),
            by = c('model_name' = 'best_model', 'name')) %>%
  mutate(is_best = !is.na(is_best)) %>%
  
  left_join(select(top_metric_models, name, top_models) %>%
              distinct %>%
              mutate(is_top = TRUE),
            by = c('model_name' = 'top_models', 'name')) %>%

  mutate(is_top = !is.na(is_top)) %>%
  mutate(name = str_replace(name, 'mcnemar_statistic', 'log_mcnemar_statistic')) %>%
  mutate(name2 = case_when(name == 'accuracy' ~ 'Accuracy',
                           name == 'j_index' ~ "Youden's J",
                           name == 'log_mcnemar_statistic' ~ 'log10(McNemar X2)',
                           name == 'roc_auc' ~ 'ROC/AUC',
                           name == 'sens' ~ 'Sensitivity',
                           name == 'spec' ~ 'Specificity')) %>% 
  
  # mutate(model_name = case_when(model_name == 'bag_mars' ~ 'Accuracy',
  #                               model_name == 'bag_tree' ~ 'Accuracy',
  #                               model_name == 'boostRF' ~ 'Accuracy',
  #                               model_name == 'c5' ~ 'Accuracy',
  #                               model_name == 'fda' ~ 'Accuracy',
  #                               model_name == 'glm' ~ 'Accuracy',
  #                               model_name == 'knn' ~ 'Accuracy',
  #                               model_name == 'lda' ~ 'Accuracy',
  #                               model_name == 'mars' ~ 'Accuracy',
  #                               model_name == 'nb' ~ 'Accuracy',
  #                               model_name == 'rda' ~ 'Accuracy',
  #                               model_name == 'rf' ~ 'Accuracy'
  #                               ))
  
  
  ggplot(aes(x = model_name, y = mean, ymin = lwr, ymax = upr)) +
  geom_linerange(show.legend = FALSE) +
  geom_point() +
  geom_hline(data = tibble(name2 = 'log10(McNemar X2)', mean = log(qchisq(0.05, 1, lower.tail = FALSE), base = 10)),
             aes(yintercept = mean), linetype = 'dashed') +
  facet_wrap(~ name2, scales = 'free_y') +
  labs(x = NULL,
       y = 'Value') +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  theme_classic() +
  theme(axis.text = element_text(colour = 'black'))


site_stats %>%
  filter(model_name == 'c5') %>%
  pivot_longer(cols = accuracy:spec) %>%
  ggplot(aes(y = Site, x = value, colour = name)) +
  geom_point()


#### Cluster models/sites ####
site_mat <- site_stats %>%
  filter(!model_name %in% c('rbf_svm')) %>%
  select(-confusion, -mcnemar_p.value) %>%
  mutate(model_site = str_c(model_name, Site, sep = '.')) %>%
  select(model_name, Site, model_site, everything()) %>%
  recipe(~ ., data = .) %>%
  step_YeoJohnson(all_numeric()) %>%
  step_normalize(all_numeric()) %>%
  prep %>%
  juice

vegan::adonis(select(site_mat, -model_name, -Site, -model_site) %>%
                as.matrix() ~ Site + model_name, data = site_mat, method = 'euclidean')


nmds_results <- site_mat %>%
  select(-model_name, -Site, -model_site) %>%
  as.matrix() %>%
  set_rownames(site_mat$model_site) %>%
  vegan::metaMDS(distance = 'euclidian', k = 2, autotransform = FALSE, trymax = 25) %$%
  points %>%
  as_tibble(rownames = 'model_site') %>%
  separate(model_site, into = c('model_name', 'site'), sep = '\\.')


nmds_results %>%
  nest(data = c(-model_name)) %>%
  mutate(ellipse = map(data, function(x) x %>%
                       dplyr::select(contains('MDS')) %>%
                       dplyr::select(1:2) %>%
                       veganCovEllipse())) %>%
  unnest(ellipse) %>%
  ggplot(aes(x = MDS1, y = MDS2, colour = model_name)) +
  geom_path() +
  geom_point(data = nmds_results)


site_mat %>%
  select(-model_name, -Site, -model_site) %>%
  as.matrix() %>%
  set_rownames(site_mat$model_site) %>%
  Rtsne::Rtsne(perplexity = 10) %$%
  Y %>%
  as_tibble() %>%
  bind_cols(site_mat) %>%
  ggplot(aes(x = V1, y = V2)) +
  geom_text(aes(label = model_name), show.legend = FALSE)


#### Site plots/tables for top model (c5) ####
site_stats %>%
  filter(model_name == 'c5') %>%
  select(-model_name, -confusion) %>%
  mutate(mcnemar_p.value = p.adjust(mcnemar_p.value, 'holm')) %>%
  select(Site, j_index, everything()) %>%
  mutate(Site = factor(Site, levels = str_c('BZ17', c('10KN', '5KN', '1KN', '500N', '100N', 
                                                      '0A', '0B', '60S', '100S', '500S', '1KS', '5KS'),
                                            sep = '-'))) %>%
  arrange(Site) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  write_csv('tmp_table.csv')

#### Parameter selection ####
overall_stats %>%
  select(model_name, MDS1, best_params) %>%
  mutate(best_params = map(best_params, ~.x %>%
                             mutate(across(everything(), as.character)) %>%
                             pivot_longer(cols = everything(),
                                          names_to = 'parameter',
                                          values_to = 'value'))) %>%
  unnest(best_params) %>%
  filter(parameter != '.config') %>%
  mutate(model_name = case_when(model_name == 'rf' ~ 'Random Forest',
                                model_name == 'knn' ~ 'k-Nearest Neighbors',
                                model_name == 'boostRF' ~ 'Boosted Random Forest',
                                model_name == 'glm' ~ 'Logistic Regression',
                                model_name == 'mars' ~ 'Multivariate Adaptive Regression Splines (MARS)',
                                model_name == 'fda' ~ 'Flexible Discriminant Analysis',
                                model_name == 'lda' ~ 'Linear Discriminant Analysis',
                                model_name == 'rda' ~ 'Regularized Discriminant Analysis',
                                model_name == 'nb' ~ 'Naive Bayes',
                                model_name == 'c5' ~ 'Decision Tree',
                                model_name == 'bag_mars' ~ 'Bagged MARS',
                                model_name == 'bag_tree' ~ 'Bagged Decision Tree',
                                )) %>%
  arrange(MDS1) %>%
  select(-MDS1) %>%
  write_csv('tmp_table.csv')

#### Classifier interpretation ####
c5_model <- read_rds('../Results/Habitat_Classification/All_c5_completeModel.rds')

c5_model %>%
  pull_workflow_preprocessor()

tmp <- pull_workflow_prepped_recipe(c5_model)
tmp$steps[[3]]$res$sdev / sum(tmp$steps[[3]]$res$sdev)
tmp$steps[[3]]$res$rotation



c5_model %>%
  pull_workflow_fit() %$%
  vi(fit, "model") %>%
  mutate(Variable = str_remove(Variable, 'Site_'),
         # Importance = Importance - 100,
         Variable = fct_reorder(Variable, Importance)) %>%
  
  ggplot(aes(x = Importance, y = Variable)) +
  geom_col()

library(C50)
rule_set <- c5_model %>%
  pull_workflow_fit() %$%
  summary(fit) %$%
  output %>%
  extract2(1) %>% 
  str_replace_all('\n\t', ';;') %>% 
  str_split('\n') %>%
  unlist %>%
  str_subset('^Rule [0-9].*') %>%
  tibble(lines = .) %>%
  
  separate(lines, sep = ';;->', into = c('lines', 'Result')) %>%
  mutate(Result = str_trim(Result)) %>%
  
  # slice(c(274, 1535)) %>%
  separate(lines, sep = ';;', into = LETTERS[1:15], fill = 'right') %>%
  separate(A, into = c('Trial', 'Rule'), sep = '/', extra = 'merge') %>%
  separate(Rule, into = c('Rule', 'Stat'), sep = ': ') %>%
  mutate(Trial = str_remove(Trial, 'Rule ')) %>%
  
  # filter(Trial == '0', Rule == '1') %>%
  
  pivot_longer(cols = -c(Trial:Stat, Result)) %>%
  filter(!is.na(value)) %>%
  mutate(parameter = str_extract(value, '[a-zA-Z0-9\\._]+'),
         value = str_remove(value, parameter) %>% str_trim,
         parameter = str_remove(parameter, 'Site_'),
         number = parse_number(value),
         sign = str_extract(value, '[=<>]+'),
         value = str_remove(value, sign) %>% str_trim) %>%
  
  group_by(Trial, Rule, Stat, Result, parameter) %>%
  mutate(number = round(number, 4),
         min_val = min(number),
         max_val = max(number)) %>%
  mutate(out = case_when(min_val == max_val ~ str_c(sign, value, sep = ' '),
                         TRUE ~ str_c(value[number == min(number)], ' < x <= ', value[number == max(number)], sep = ''))) %>%
  summarise(out = unique(out), .groups = 'drop') %>%
  mutate(out = str_c("'", out)) %>%
  
  pivot_wider(names_from = 'parameter',
              values_from = 'out', 
              values_fill = NULL) %>%
  mutate(across(Trial:Rule, as.integer)) %>%
  arrange(Trial, Rule) %>%
  select(Trial:Result, starts_with('PC'), everything()) %>%
  separate(Result, into = c('Result', 'Probability'), sep = '  ') %T>%
  write_csv('tmp_table.csv', na = "")


rule_set %>%
  mutate(no_site = rowSums(across(starts_with('BZ17'), is.na)) == 11) %>%
  
  mutate(keeps = rowSums(across(starts_with('BZ17'), ~. == "'> 0" & !is.na(.))),
         cancels = rowSums(across(starts_with('BZ17'), ~. == "'<= 0" & !is.na(.)))) %>%
  
  mutate(across(starts_with('BZ17'), ~(. == "'> 0")),
         across(starts_with('BZ17'), ~(. | no_site)),
         across(starts_with('BZ17'), ~if_else(is.na(.), (keeps < 1 & cancels > 0), .)),
         BZ17.0A = no_site | (keeps < 1 & cancels > 0)) %>%
  select(-no_site:-cancels) %>%
  pivot_longer(cols = starts_with('BZ17'), 
               names_to = 'Site',
               values_to = 'present') %>%
  filter(present) %>%
  select(-present) %>%
  mutate(Site = factor(Site, levels = str_c('BZ17', c('10KN', '5KN', '1KN', '500N', '100N', 
                                                      '0A', '0B', '60S', '100S', '500S', '1KS', '5KS'),
                                            sep = '.')),
         Site = LETTERS[as.integer(Site)]) %>%
  group_by(Trial, Rule, Result, Probability, PC1, PC2, PC3) %>%
  arrange(Trial, Rule, Site) %>%
  summarise(Site = str_c(Site, collapse = ', '), .groups = 'drop') %>%
  mutate(across(starts_with('PC'), ~replace_na(., '')),
         Site = str_replace_all(Site, ' +', ' ')) %>%
  rowwise() %>%
  mutate(Excluded = LETTERS[which(str_detect(Site, LETTERS[1:12], negate = TRUE))] %>%
           str_c(collapse = ', ') %>%
           str_replace_all(' +', ' ')) %>%
  ungroup %>%
  mutate(across(c(Site, Excluded), ~case_when(. == str_c(LETTERS[1:12], collapse = ', ') ~ 'Any',
                                              . == '' ~ 'None',
                                              TRUE ~ .)),
         Result = if_else(Result == 'class C1', 'Sand', 'Reef'),
         Probability = parse_number(Probability)) %>%
  write_csv('tmp_table.csv', na = "")


library(magrittr)

LETTERS[Site %>% str_detect(LETTERS[1:12]) %>% not %>% which] %>%
  str_c(collapse = ', ') %>%
  str_replace_all(' +', ' ')


LETTERS[which(str_detect(Site, LETTERS[1:12], negate = TRUE))] %>%
  str_c(collapse = ', ') %>%
  str_replace_all(' +', ' ')



Site = if_else(Site == str_c(LETTERS[1:12], collapse = ', '), 'Any', Site),


tmp <- rule_set %>%
  mutate(no_site = rowSums(across(starts_with('BZ17'), is.na)) == 11) %>%
  
  mutate(keeps = rowSums(across(starts_with('BZ17'), ~. == "'> 0" & !is.na(.))),
         cancels = rowSums(across(starts_with('BZ17'), ~. == "'<= 0" & !is.na(.)))) %>%
  
  mutate(across(starts_with('BZ17'), ~(. == "'> 0")),
         across(starts_with('BZ17'), ~(. | no_site)),
         across(starts_with('BZ17'), ~if_else(is.na(.), (keeps < 1 & cancels > 0), .)),
         BZ17.0A = no_site | (keeps < 1 & cancels > 0)) %>%
  select(-no_site:-cancels) %>%
  
  pivot_longer(cols = starts_with('BZ17'), 
               names_to = 'Site',
               values_to = 'present') %>%
  filter(present) %>%
  select(-present) %>%
  
  nest(rules = -c(Site))


tmp$rules[[1]] %>%
  filter(Trial == 0)


tmp %>%
  mutate(Site = factor(Site, levels = str_c('BZ17', c('10KN', '5KN', '1KN', '500N', '100N', 
                                                      '0A', '0B', '60S', '100S', '500S', '1KS', '5KS'),
                                            sep = '.')),
         Site = LETTERS[as.integer(Site)]) %>%
  arrange(Site) %>%
  unnest(rules) %>%
  group_by(Site, Trial) %>%
  mutate(Rule = 1:n()) %>%
  ungroup %>%
  select(Trial:Probability, Site, starts_with("PC")) %>%
  arrange(Trial, Site, Rule)



tmp$rules[[1]] %>%
  filter(Trial == 0) %>%
  filter(Result == 'class C1') %>%
  pivot_longer(cols = starts_with('PC')) %>%
  filter(!is.na(value)) %>%
  mutate(direction = str_extract(value, '[><=]+'),
         value = parse_number(value)) %>%
  mutate(xend = case_when(direction == '>' ~ Inf,
                          direction == '<=' ~ -Inf)) %>%
  
  ggplot(aes(x = value, y = Rule, colour = direction, xend = xend, yend = Rule)) +
  geom_segment() +
  facet_wrap(~name)



rule_set %>%
  filter(Trial == 0,
         Rule %in% tmp$rules[[1]]$Rule[tmp$rules[[1]]$Trial == 0]) %>% View


tmp$rules[[1]] %>%
  filter(Trial == 1)
