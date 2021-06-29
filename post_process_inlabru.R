#### Post-Process many inlabru models ####
model_prefix <- '17.5.21'

##TODO Figure out how to merge hyper params
##TODO posterior predictive envelopes

library(sf)
library(terra)
library(tidyverse)
library(magrittr)
library(inlabru)
library(INLA)
library(patchwork)
library(lubridate)
library(units)
library(furrr)

plan('multicore')

#### Get Results ####
RESULTS_DIR <- '/work/hobi/jselwyn/Habitat/Results/INLA_models'
OUT_DIR <- str_c(RESULTS_DIR, '/models_Overall.', model_prefix)
dir.create(OUT_DIR, showWarnings = FALSE)

#### Runs which didn't manage to Run ####

model_fitting <- list.dirs(RESULTS_DIR, full.names = TRUE) %>%
  tibble(model_folders = .) %>%
  filter(str_detect(model_folders, model_prefix),
         str_detect(model_folders, 'Overall', negate = TRUE)) %>%
  mutate(inner_files = future_map(model_folders, ~list.files(.x, pattern = 'rds$', recursive = FALSE, full.names = FALSE)),
         model_finished = map_lgl(inner_files, ~length(.x) > 0)) %>%
  select(-inner_files)

model_fitting %>%
  filter(!model_finished)



#### Number of Fish Diagnostic ####
site_areas <- list.files(RESULTS_DIR,
                         pattern = 'shp$', recursive = TRUE, full.names = TRUE) %>%
  str_subset(model_prefix) %>%
  str_subset('siteOutline') %>%
  str_subset('Overall', negate = TRUE) %>%
  tibble(file = .) %>%
  mutate(version = str_extract(file, '[0-9\\.]+')) %>%
  mutate(date = str_extract(version, '[0-9]+\\.[0-9]+\\.[0-9]+') %>% dmy,
         version = str_remove(version, str_c(model_prefix, '.')) %>% as.integer) %>%
  select(-date) %>%
  rowwise %>%
  summarise(version = version, st_read(file, quiet = TRUE), .groups = 'drop') %>%
  mutate(area = st_area(geometry)) %>%
  select(version, Site, area) %>%
  mutate(Site = as.character(Site))

fish_site <- list.files(RESULTS_DIR,
                        pattern = 'fishCount', 
                        full.names = TRUE, recursive = TRUE) %>%
  str_subset(model_prefix) %>%
  str_subset('Overall', negate = TRUE) %>%
  tibble(file = .) %>%
  mutate(version = str_extract(file, '[0-9\\.]+\\.csv')) %>%
  mutate(date = str_extract(version, '[0-9]+\\.[0-9]+\\.[0-9]+') %>% dmy,
         version = str_extract(version, '[0-9]+\\.csv') %>% str_remove('\\.csv') %>% as.integer) %>%
  select(-date) %>%
  rowwise %>%
  summarise(version = version,
            read_csv(file,
                     col_types = cols(
                       .default = col_double(),
                       Site = col_character()
                     )),
            .groups = 'drop') %>%
  left_join(site_areas, by = c('version', 'Site')) %>%
  mutate(across(c(obs, mean, q0.025, median, q0.975, smin, smax), ~./area %>% as.numeric)) 

fish_site_plot <- fish_site %>%
  ggplot(aes(x = mean, y = Site, group = as.character(version))) +
  geom_linerange(aes(xmin = q0.025, xmax = q0.975), 
                  linetype = 'dashed', position = position_dodge(0.1),
                  colour = 'grey50', size = 0.1) +
  geom_point(aes(x = obs), colour = 'black', size = 5) +
  # scale_x_continuous(trans = scales::log10_trans()) +
  labs(x = '# Fish / m2',
       y = NULL,
       colour = 'Model Run') +
  theme_classic() +
  theme(axis.text = element_text(size = 14, colour = 'black'),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))

#### PPP Diagnostic ####
ppp_results <- list.files(RESULTS_DIR,
                          pattern = 'Envelopes',
                          full.names = TRUE, recursive = TRUE) %>%
  str_subset(model_prefix) %>%
  str_subset('Overall', negate = TRUE) %>%
  tibble(file = .) %>%
  mutate(version = str_extract(file, '[0-9\\.]+'),
         type = if_else(str_detect(file, 'Site'), 'Site', 'Global')) %>%
  mutate(date = str_extract(version, '[0-9]+\\.[0-9]+\\.[0-9]+') %>% dmy,
         version = str_remove(version, str_c(model_prefix, '.')) %>% as.integer) %>%
  select(-date) %>%
  mutate(data = map(file, read_rds)) %>%
  dplyr::select(-file) %>%
  pivot_wider(names_from = 'type',
              values_from = 'data') %>%
  rowwise %>%
  filter(!(is.null(Global) | is.null(Site))) %>%
  mutate(Global = list(mutate(Global, Site = 'all') %>%
                         pivot_wider(names_from = 'metric', values_from = 'envelope'))) %>%
  pivot_longer(cols = -version,
               names_to = 'scale',
               values_to = 'ppp') %>%
  unnest(ppp) %>%
  pivot_longer(cols = Kest:Jest,
               names_to = 'metric',
               values_to = 'envelope') %>%
  mutate(extras = map(envelope, ~attr(., 'einfo')),
         numerator = map_dbl(extras, ~2 * .x$nr),
         denominator = map_dbl(extras, ~(.x$nsim + 1)),
         p.value = numerator / denominator) %>%
  select(-extras)

ppp_results$envelope[[1]]

ppp_results %>%
  mutate(envelope = map(envelope, as_tibble)) %>%
  unnest(envelope) %>%
  group_by(scale, Site, metric, r) %>%
  summarise(across(obs:p.value, mean))

%>%
  
  mutate(ind_plot = pmap(list(envelope, metric, version, numerator, denominator, p.value),
                         ~..1 %>%
                           ggplot(aes(x = r, y = obs, ymin = lo, ymax = hi)) +
                           geom_ribbon(alpha = 0.4) +
                           geom_line() +
                           labs(title = str_c(..2, ': ', ..3),
                                subtitle = str_c('MC Sig: ', ..4, '/', ..5, ' = ', ..6)) +
                           theme_classic()))


#### Parameter Estimates ####
parameter_margins <- list.files(RESULTS_DIR,
           pattern = 'paramMargin', 
           full.names = TRUE, recursive = TRUE) %>%
  str_subset(model_prefix) %>%
  str_subset('Overall', negate = TRUE) %>%
  tibble(file = .) %>%
  # sample_n(4) %>% 
  mutate(version = str_extract(file, '[0-9\\.]+\\.csv')) %>%
  mutate(date = str_extract(version, '[0-9]+\\.[0-9]+\\.[0-9]+') %>% dmy,
         version = str_extract(version, '[0-9]+\\.csv') %>% str_remove('\\.csv') %>% as.integer) %>%
  select(-date) %>%
  mutate(data = future_map(file, ~read_csv(.x, 
                                           col_types = cols(
                                             type = col_character(),
                                             parameter = col_character(),
                                             index = col_character(),
                                             x = col_double(),
                                             y = col_double()
                                           )) %>% 
                             filter(parameter != 'spatialSmooth'),
                           .progress = TRUE, .options = furrr_options(seed = TRUE))) %>%
  select(-file) %>%
  unnest(data)

#Get parameter estimates for each model run
parameter_distributions_version <- parameter_margins %>%
  filter(parameter != 'spatialSmooth') %>%
  filter(type != 'hyper_par') %>% #Figure out how to keep this in!
  mutate(x = case_when(parameter %in% c('beta.moran', 'beta.depth') ~ x*-1,
                       TRUE ~ x)) %>%
  group_by(version, type, parameter, index) %>%
  summarise(inla.zmarginal(cbind(x, y), silent = TRUE) %>%
              bind_rows,
            prob_positive = 1 - inla.pmarginal(0, cbind(x, y)),
            prob_negative = inla.pmarginal(0, cbind(x, y)),
            
            posterior = list(inla.smarginal(cbind(x, y))), 
            hpd = list(inla.hpdmarginal(p = c(0.5, 0.75, 0.9, 0.95, 0.99), marginal = cbind(x, y)) %>%
                         as_tibble(rownames = 'level') %>%
                         mutate(level = str_remove(level, 'level:') %>% as.numeric)),
            .groups = 'drop') %>%
  rowwise() %>%
  summarise(version = version,
            type = type,
            parameter = str_remove(parameter, 'beta\\.'),
            index = index,
            across(c(mean, sd,
                     starts_with('quant'), 
                     starts_with('prob')), 
                   identity),
            evidence_positive = prob_positive / prob_negative,
            evidence_negative = prob_negative / prob_positive, 
            hpd = list(hpd),
            distributions = list(bind_rows(posterior) %>%
                                   rename(posterior = y) %>%
                                   mutate(posterior = posterior/sum(posterior))),
            .groups = 'drop') %>%
  mutate(significant = (quant0.025 < 0 & quant0.975 < 0) | (quant0.025 > 0 & quant0.975 > 0)) %>%
  select(version, type, parameter, index, where(is.numeric), where(is.logical), where(is.list)) %>%
  mutate(parameter = case_when(parameter == 'moran' ~ 'Coarse Complexity',
                               parameter == 'complexity' ~ 'Fine Complexity',
                               parameter == 'boundaryDistance' ~ 'Distance to Sand/Reef Margin',
                               parameter == 'habitat' ~ 'Sand',
                               TRUE ~ str_to_sentence(parameter)))

#Get overall parameter estimates
parameter_distributions <- parameter_margins %>%
  filter(parameter != 'spatialSmooth') %>%
  mutate(x = case_when(parameter %in% c('beta.moran', 'beta.depth') ~ x*-1,
                       TRUE ~ x)) %>%
  mutate(x = if_else(str_detect(parameter, 'Precision'), log(x), x)) %>%
  group_by(version, type, parameter, index) %>%
  mutate(distribution = cbind(x, y)) %>%
  select(-x, -y) %>%
  summarise(inla.smarginal(distribution) %>% bind_cols(), .groups = 'keep') %>%
  mutate(y = y/sum(y)) %>%
  
  mutate(marginal = cbind(x = x, y = y)) %>%
  group_by(type, parameter, index) %>%
  mutate(min_x = min(x),
         max_x = max(x)) %>%
  select(-x, -y) %>%
  ungroup %>%
  nest_by(version, type, parameter, index, min_x, max_x) %>%
  summarise(posterior = bind_cols(x = seq(min_x, max_x, length.out = 1000),
                                  y = inla.dmarginal(seq(min_x, max_x, length.out = 1000), data$marginal)),
            .groups = 'drop') %>%
  mutate(x = posterior$x,
         y = posterior$y) %>%
  select(-posterior) %>%
  group_by(type, parameter, index, x) %>%
  summarise(y = sum(y), .groups = 'drop_last') %>%
  mutate(x = if_else(str_detect(parameter, 'Precision'), sqrt(1/exp(x)), x),
         y = y/sum(y),
         posterior = cbind(x, y)) %>%
  select(-x, -y) %>%
  
  summarise(inla.zmarginal(posterior, silent = TRUE) %>%
              bind_rows,
            prob_positive = 1 - inla.pmarginal(0, posterior),
            prob_negative = inla.pmarginal(0, posterior),
            
            hpd = list(inla.hpdmarginal(p = c(0.5, 0.75, 0.9, 0.95, 0.99), marginal = posterior) %>%
                         as_tibble(rownames = 'level') %>%
                         mutate(level = str_remove(level, 'level:') %>% as.numeric)),
            
            posterior = list(inla.smarginal(posterior)), 
            .groups = 'drop') %>%
  mutate(parameter = str_replace(parameter, 'Precision', 'Stdev')) %>%
  
  rowwise() %>%
  summarise(type = type,
            parameter = str_remove(parameter, 'beta\\.'),
            index = index,
            across(c(mean, sd,
                     starts_with('quant'), 
                     starts_with('prob')), 
                   identity),
            evidence_positive = prob_positive / prob_negative,
            evidence_negative = prob_negative / prob_positive, 
            hpd = list(hpd),
            distributions = list(bind_rows(posterior) %>%
                                   rename(posterior = y) %>%
                                   mutate(posterior = posterior/sum(posterior))),
            .groups = 'drop') %>%
  mutate(significant = (quant0.025 < 0 & quant0.975 < 0) | (quant0.025 > 0 & quant0.975 > 0)) %>%
  select(type, parameter, index, where(is.numeric), where(is.logical), where(is.list)) %>%
  mutate(parameter = case_when(parameter == 'moran' ~ 'Coarse Complexity',
                               parameter == 'complexity' ~ 'Fine Complexity',
                               parameter == 'boundaryDistance' ~ 'Distance to Sand/Reef Margin',
                               parameter == 'habitat' ~ 'Sand',
                               TRUE ~ str_to_sentence(parameter)))

#Put together for plot to show versions well represented by summary
parameter_version_plot <- filter(parameter_distributions, type == 'fixed' | parameter == 'Sand', parameter != 'Intercept') %>%
  mutate(parameter = fct_reorder(parameter, mean)) %>%
  
  ggplot(aes(y = parameter, x = mean, xmin = quant0.025, xmax = quant0.975)) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  
  geom_linerange(data = filter(parameter_distributions_version, type == 'fixed' | parameter == 'Sand', parameter != 'Intercept') %>%
                   mutate(parameter = fct_reorder(parameter, mean)), 
                  aes(group = as.character(version)),
                 linetype = 'dashed', colour = 'grey50', size = 0.1,
                 position = position_dodge(0.1)) +
  
  geom_pointrange(linetype = 'dashed', colour = 'black') +
  
  labs(y = NULL,
       x = 'Standardized Linear Effect on Density') +
  theme_classic() +
  theme(axis.text = element_text(size = 14, colour = 'black'),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))


#### Output Summary Files ####
write_csv(fish_site, str_c(OUT_DIR, '/fish_density_diagnostic.csv'))
ggsave(str_c(OUT_DIR, '/version_posteriorFish.png'), plot = fish_site_plot, height = 15, width = 7)

write_rds(parameter_distributions, str_c(OUT_DIR, '/overall_parameters.rds'))
ggsave(str_c(OUT_DIR, '/version_parameters.png'), plot = parameter_version_plot, height = 15, width = 7)

#### Merge Prediction Rasters ####
prediction_out <- list.files(RESULTS_DIR,
           pattern = 'tif$', recursive = TRUE, full.names = TRUE) %>%
  str_subset(model_prefix) %>%
  str_subset('unsmoothed') %>%
  str_subset('noSpatial', negate = TRUE) %>%
  
  tibble(file = .) %>%
  mutate(version = str_extract(file, '[0-9\\.]+'),
         site = str_extract(file, 'BZ17-[0-9ABKNS]+')) %>%
  mutate(date = str_extract(version, '[0-9]+\\.[0-9]+\\.[0-9]+') %>% dmy,
         version = str_remove(version, str_c(model_prefix, '.')) %>% as.integer) %>%
  select(-date) %>%
  
  # sample_n(15) %>%
  
  rowwise %>%
  mutate(preds = list(rast(file)$mean)) %>%
  group_by(site) %>%
  summarise(preds = list(do.call(c, preds) %>% mean),
            .groups = 'drop') %>%
  rename(unsmoothed = preds) %>%
  mutate(circle_mat = future_map(unsmoothed, 
                                 ~raster::focalWeight(raster::raster(.x$mean), sqrt(1/pi), 'circle'))) %>%
  rowwise %>%
  mutate(smoothed = list(map(names(unsmoothed), ~terra::focal(unsmoothed[[.x]], w = circle_mat, fun = 'sum') %>%
                               cover(y = unsmoothed[[.x]])) %>%
                           rast
  )) %>%
  ungroup %>%
  select(-circle_mat) %>%
  pivot_longer(cols = c(unsmoothed, smoothed)) %>%
  mutate(out_name = str_c(OUT_DIR, '/', site, '_', name, '.tif',
                          sep = '')) %>%
  mutate(write = future_map2(value, out_name, ~writeRaster(.x, .y, overwrite = TRUE))) %>%
  select(-write)
