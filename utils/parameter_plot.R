library(tidyverse)
library(units)

#### Fish diagnostic ####
read_csv()


#### Parameters ####
## Distribution Model 
posterior_hpd_distribution <- read_rds('../../overall_parameters.rds') %>%
  select(-sd:-quant0.975, -distributions) %>%
  rowwise %>%
  mutate(hpd = list(pivot_wider(hpd, 
                                names_from = 'level',
                                values_from = c('low', 'high')))) %>%
  unnest(hpd) %>%
  select(type, parameter, index, mean, 
         contains(str_c(c('low'), c(0.99, 0.95, 0.9, 0.75, 0.5), sep = '_')),
         contains(str_c(c('high'), rev(c(0.99, 0.95, 0.9, 0.75, 0.5)), sep = '_')),
         significant,
         everything())

posterior_hpd_distribution %>%
  select(parameter, mean, ends_with('95')) %>%
  filter(parameter == 'Intercept') %>%
  mutate(across(c(mean, ends_with('95')), ~exp(.))) %>%
  mutate(across(c(mean, ends_with('95')), ~set_units(., '1/m2'))) %>%
  mutate(across(c(mean, ends_with('95')), ~set_units(., '1/km2'))) 


posterior_hpd_distribution %>%
  filter(parameter != 'Intercept', parameter != 'Site', type != 'hyper_par') %>%
  mutate(across(c(mean, ends_with('95')), ~exp(.) * 100 - 100)) %>%
  left_join(read_csv('../Results/Topography/overall_metrics_c5.csv') %>%
              select(metric, sd) %>%
              mutate(metric = case_when(metric == 'moran' ~ 'Coarse Complexity',
                                        metric == 'complexity' ~ 'Fine Complexity',
                                        metric == 'boundaryDistance' ~ 'Distance to Sand/Reef Margin',
                                        metric == 'habitat' ~ 'Sand',
                                        TRUE ~ str_to_sentence(metric))),
            by = c('parameter' = 'metric')) %>%
  mutate(across(c(mean, ends_with('95')), ~if_else(is.na(sd), ., .*sd))) %>%
  select(parameter, mean, ends_with('95'), starts_with('prob'), starts_with('evidence'))

posterior_hpd_distribution %>%
  filter(parameter == 'Site')


posterior_hpd_distribution %>%
  filter(type == 'hyper_par') %>%
  select(parameter, mean, ends_with('95'))

## Composition Model
composition_model <- read_rds('../Intermediate_Files/BRMS_models/shoal_composition_model.rds')


posterior_hpd_compostion <- posterior_interval(composition_model, pars = '^b_') %>%
  as_tibble(rownames = 'parameter') %>%
  rename(quant0.025 = '2.5%',
         quant0.975 = '97.5%') %>%
  mutate(significant = (quant0.025 < 0 & quant0.975 < 0) | (quant0.025 > 0 & quant0.975 > 0)) %>%
  mutate(parameter = str_remove(parameter, '^b_')) %>%
  select(-starts_with('quant')) %>%
  full_join(bayestestR::hdi(composition_model, ci = c(0.99, 0.95, 0.9, 0.75, 0.5)) %>%
              as_tibble() %>%
              select(-Component) %>%
              rename(type = Effects,
                     parameter = Parameter,
                     low = CI_low,
                     high = CI_high) %>%
              pivot_wider(names_from = 'CI', values_from = c('low', 'high')) %>%
              mutate(parameter = str_remove(parameter, '^b_')),
            by = 'parameter') %>%
  full_join(hypothesis(composition_model, c('Shoal.Size > 0', 'depth > 0', 
                                            'boundaryDistance > 0', 'viewshed > 0', 
                                            'coarse_complexity > 0', 'fine_complexity > 0',
                                            'habitatSand > 0'), alpha = 0.05) %$%
              hypothesis %>%
              as_tibble %>%
              mutate(parameter = str_remove(Hypothesis, '\\).*0') %>% str_remove('\\(')) %>%
              select(parameter, Estimate, Evid.Ratio, Post.Prob) %>%
              rename(mean = Estimate,
                     evidence_positive = Evid.Ratio, 
                     prob_positive = Post.Prob),
            by = 'parameter') %>%
  full_join(hypothesis(composition_model, c('Shoal.Size < 0', 'depth < 0', 
                                            'boundaryDistance < 0', 'viewshed < 0', 
                                            'coarse_complexity < 0', 'fine_complexity < 0',
                                            'habitatSand < 0'), alpha = 0.05) %$%
              hypothesis %>%
              as_tibble %>%
              mutate(parameter = str_remove(Hypothesis, '\\).*0') %>% str_remove('\\(')) %>%
              select(parameter, Estimate, Evid.Ratio, Post.Prob) %>%
              rename(mean = Estimate,
                     evidence_negative = Evid.Ratio, 
                     prob_negative = Post.Prob) ,
            by = c('parameter', 'mean')) %>%
  mutate(parameter = case_when(parameter == 'coarse_complexity' ~ 'Coarse Complexity',
                               parameter == 'fine_complexity' ~ 'Fine Complexity',
                               parameter == 'boundaryDistance' ~ 'Distance to Sand/Reef Margin',
                               parameter == 'habitatSand' ~ 'Sand',
                               parameter == 'Shoal.Size' ~ 'Shoal Size',
                               TRUE ~ str_to_sentence(parameter)))


#Both

bind_rows(Distribution = posterior_hpd_distribution,
          Composition = posterior_hpd_compostion,
          .id = 'model')


bind_rows(Distribution = posterior_hpd_distribution,
          Composition = posterior_hpd_compostion,
          .id = 'model') %>%
  filter(parameter != 'Intercept', parameter != 'Site') %>%
  mutate(model = factor(model, levels = c('Distribution', 'Composition')),
         parameter = fct_reorder(parameter, mean)) %>%
  ggplot(aes(y = parameter, x = mean, xmin = low_0.95, xmax = high_0.95)) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_pointrange(linetype = 'dashed') +
  geom_linerange(aes(xmin = low_0.5, xmax = high_0.5)) +
  facet_wrap(~model) +
  labs(y = NULL,
       x = 'Standardized Linear Effect on Density') +
  theme_classic() +
  theme(axis.text = element_text(size = 14, colour = 'black'),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))
ggsave('../Manuscript/Figures/parameter_estimates.png')

