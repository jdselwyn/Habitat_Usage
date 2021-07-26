#### Post-Process many inlabru models ####


if(interactive()){
  model_prefix <- '9.7.21'
} else {
  args <- commandArgs(trailingOnly=TRUE)
  model_prefix <- args[1]
}

message(model_prefix)

suppressMessages(library(sf))
suppressMessages(library(terra))
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(inlabru))
suppressMessages(library(INLA))
suppressMessages(library(patchwork))
suppressMessages(library(lubridate))
suppressMessages(library(units))
suppressMessages(library(furrr))
suppressMessages(library(maptools))
suppressMessages(library(spatstat))

plan('multicore')

#TODO - figure out way to not have to randomly pick one observation for last bit

#### Get Results ####
RESULTS_DIR <- '/work/hobi/jselwyn/Habitat/Results/INLA_models'
OUT_DIR <- str_c(RESULTS_DIR, '/models_Overall.', model_prefix)
dir.create(OUT_DIR, showWarnings = FALSE)


site_rename <- tibble(Site = str_c('BZ17-', c('0A', '0B', '100N', '100S', '10KN', '1KN',
                                              '1KS', '500N', '500S', '5KN', '5KS', '60S')),
                      site = c('F', 'G', 'E', 'I', 'A', 'C', 'K', 'D', 'J', 'B', 'L', 'H'))

#### Runs which didn't manage to Run ####
model_fitting <- list.dirs(RESULTS_DIR, full.names = TRUE) %>%
  tibble(model_folders = .) %>%
  filter(str_detect(model_folders, model_prefix),
         str_detect(model_folders, 'Overall', negate = TRUE)) %>%
  mutate(inner_files = future_map(model_folders, ~list.files(.x, 
                                                             pattern = 'pointProcessSample.*shp$',
                                                             recursive = FALSE, 
                                                             full.names = FALSE)),
         model_finished = map_lgl(inner_files, ~length(.x) > 0)) %>%
  select(-inner_files) %>%
  mutate(model_number = str_extract(model_folders, '[0-9]+$') %>% as.integer)

message('Total number of attempted model fits: ', length(model_fitting$model_number))
message('Total number of successful model fits: ', sum(model_fitting$model_finished))

model_folders <- model_fitting %>%
  filter(model_finished) %>%
  rename(version = model_number) %>%
  select(-model_finished) %>%
  # sample_n(25) %>% #for testing
  identity()

#### Number of Fish Diagnostic ####
message('Starting fish diagnostics: ', Sys.time())
site_areas <- model_folders %>%
  mutate(file = map_chr(model_folders, ~list.files(.x, pattern = 'shp$', recursive = FALSE, full.names = TRUE) %>%
                          str_subset(model_prefix) %>%
                          str_subset('siteOutline') %>%
                          str_subset('Overall', negate = TRUE))) %>%
  mutate(sites = future_map(file, st_read, quiet = TRUE, .progress = TRUE, .options = furrr_options(seed = TRUE))) %>%
  unnest(sites) %>%
  mutate(area = st_area(geometry)) %>%
  select(version, Site, area) %>%
  mutate(Site = as.character(Site)) %>%
  group_by(Site) %>%
  summarise(area = mean(area), .groups = 'drop')


fish_site <- model_folders %>%
  mutate(file = map_chr(model_folders, ~list.files(.x, pattern = 'fishCount', recursive = FALSE, full.names = TRUE) %>%
                          str_subset(model_prefix) %>%
                          str_subset('Overall', negate = TRUE))) %>%
  mutate(sites = future_map(file, ~read_csv(.x,
                                            col_types = cols(
                                              .default = col_double(),
                                              Site = col_character()
                                            )),
                            .progress = FALSE, .options = furrr_options(seed = TRUE))) %>%
  unnest(sites) %>%
  group_by(Site, obs) %>%
  summarise(mean = mean(predicted_number),
            sd = sd(predicted_number),
            q0.025 = quantile(predicted_number, 0.025),
            q0.25 = quantile(predicted_number, 0.25),
            median = median(predicted_number),
            q0.75 = quantile(predicted_number, 0.75),
            q0.975 = quantile(predicted_number, 0.975),
            smin = min(predicted_number),
            smax = max(predicted_number), 
            cv = sd / mean, 
            .groups = 'drop') %>%
  left_join(site_areas, by = c('Site')) %>%
  mutate(across(c(obs, mean, q0.025, q0.25, median, q0.75, q0.975, smin, smax), ~./area %>% as.numeric)) %>%
  left_join(site_rename, by = c('Site')) %>%
  select(Site, site, everything())

fish_site_plot <- fish_site %>%
  mutate(site = fct_rev(site)) %>%
  ggplot(aes(x = mean, y = site)) +
  geom_linerange(aes(xmin = q0.025, xmax = q0.975), 
                  linetype = 'dashed') +
  geom_linerange(aes(xmin = q0.25, xmax = q0.75), 
                 linetype = 'solid') +
  geom_point(aes(x = obs), colour = 'black') +
  # scale_x_continuous(trans = scales::log10_trans()) +
  labs(#x = '# Fish / m2',
       x = expression(Density~(fish/m^2)),
       y = NULL,
       colour = 'Model Run') +
  theme_classic() +
  theme(axis.text = element_text(size = 14, colour = 'black'),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))

write_csv(fish_site, str_c(OUT_DIR, '/fish_density_diagnostic.csv'))
ggsave(str_c(OUT_DIR, '/version_posteriorFish.png'), plot = fish_site_plot, height = 7, width = 4)
message('Finished fish diagnostics: ', Sys.time())

#### Parameter Estimates ####
message('Starting parameter estimates: ', Sys.time())

parameter_margins <- model_folders %>%
  mutate(file = map_chr(model_folders, ~list.files(.x, pattern = 'paramMargin', recursive = FALSE, full.names = TRUE) %>%
                          str_subset(model_prefix) %>%
                          str_subset('Overall', negate = TRUE))) %>%
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
  select(-model_folders, -file) %>%
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

write_rds(parameter_distributions, str_c(OUT_DIR, '/overall_parameters.rds'))
ggsave(str_c(OUT_DIR, '/version_parameters.png'), plot = parameter_version_plot, height = 15, width = 7)

message('Finish parameter estimates: ', Sys.time())

#### PPP Diagnostic ####
message('Starting Secondary PPP diagnostics: ', Sys.time())

subsampling_design <- sample(nrow(model_folders) * 100, 100) %>%
  tibble(posterior_sample = .) %>%
  mutate(version = as.integer((posterior_sample - 1) %/% 100 + 1),
         posterior_sample = as.integer((posterior_sample - 1) %% 100 + 1))

obs_version <- sample(subsampling_design$version, 1)

site_outlines <- model_folders %>%
  arrange(version) %>%
  mutate(version = 1:n()) %>%
  inner_join(select(subsampling_design, version), by = 'version') %>%
  filter(version == obs_version) %>%
  mutate(file = map_chr(model_folders, ~list.files(.x, pattern = 'siteOutline.*shp$', 
                                                   recursive = FALSE, full.names = TRUE) %>%
                          str_subset(model_prefix) %>%
                          str_subset('Overall', negate = TRUE))) %>%
  mutate(sites = future_map(file, st_read, quiet = TRUE, .options = furrr_options(seed = TRUE))) %>%
  unnest(sites) %>%
  select(-model_folders, -file) %>%
  rename(obs_version = version)


observed_fish <- model_folders %>%
  arrange(version) %>%
  mutate(version = 1:n()) %>%
  inner_join(select(subsampling_design, version), by = 'version') %>%
  filter(version == obs_version) %>%
  mutate(file = map_chr(model_folders, ~list.files(.x, pattern = 'observed.*shp$', 
                                                   recursive = FALSE, full.names = TRUE) %>%
                          str_subset(model_prefix) %>%
                          str_subset('Overall', negate = TRUE))) %>%
  mutate(observed = future_map(file, st_read, quiet = TRUE, .options = furrr_options(seed = TRUE))) %>%
  unnest(observed) %>%
  select(-model_folders, -file) %>%
  rename(obs_version = version)


posterior_samples <- model_folders %>%
  arrange(version) %>%
  mutate(version = 1:n()) %>%
  inner_join(subsampling_design, by = 'version') %>%
  mutate(file = map_chr(model_folders, ~list.files(.x, pattern = 'pointProcessSample.*shp$', 
                                                   recursive = FALSE, full.names = TRUE) %>%
                          str_subset(model_prefix) %>%
                          str_subset('Overall', negate = TRUE))) %>%
  mutate(sim = future_map(file, st_read, quiet = TRUE, .options = furrr_options(seed = TRUE))) %>%
  mutate(sim = map2(sim, posterior_sample, ~filter(.x, sample == .y))) %>%
  select(-model_folders, -file, -posterior_sample) %>%
  mutate(sim = future_map(sim, ~nest(.x, sim = geometry), .options = furrr_options(seed = TRUE))) %>%
  unnest(sim) %>%
  nest(simulations = c(version, sample, sim))

message('Finished reading in PPP diagnostic files: ', Sys.time())

## Calculate point process diagnostics
diag_metrics <- c('Kest','Lest', 'pcf', 'Fest', 'Gest', 'Jest')

site_point_processes <- observed_fish %>%
  nest(obs = geometry) %>%
  left_join(site_outlines, by = c('obs_version', 'site' = 'Site')) %>%
  
  left_join(posterior_samples, 
            by = c('site' = 'Site')) %>%
  rename(Site = site) %>%
  
  mutate(obs = future_map(obs, ~st_as_sf(.x) %>% as('Spatial'), .options = furrr_options(seed = TRUE)),
         site = future_map(geometry, ~as(.x, 'Spatial'), .options = furrr_options(seed = TRUE)),
         simulations = future_map(simulations, ~pull(.x, sim) %>% map(~st_as_sf(.x) %>% as('Spatial')),
                                  .options = furrr_options(seed = TRUE))) %>%
  select(-geometry)

message('Finished sampling posterior for site PPP metrics: ', Sys.time())

ppp_metrics <- site_point_processes %>%
  
  rowwise %>%
  mutate(site = list(as(site, 'owin')),
         obs = list(ppp(x = coordinates(obs)[,1],
                        y = coordinates(obs)[,2],
                        window = site)),
         simulations = list(map(simulations, ~ppp(x = coordinates(.x)[,1],
                                                  y = coordinates(.x)[,2],
                                                  window = site))
         )) %>%
  ungroup %>%
  full_join(expand_grid(Site = .$Site,
                        metric = diag_metrics),
            by = 'Site') %>%
  select(Site, metric, everything()) %>%
  mutate(nsim = map_int(simulations, length)) %>%
  mutate(model = future_pmap(list(obs, simulations, nsim, metric),
                             ~envelope(..1, ..4, simulate = ..2, nsim = (..3 - 1)),
                             .progress = TRUE, .options = furrr_options(seed = TRUE))) %>%
  pivot_wider(names_from = 'metric',
              values_from = 'model') %>%
  select(-obs, -simulations, -site)

message(str_c('Finished measuring site level posterior point processes: ', Sys.time()))



#Plot ppp results 
ppp_results <- ppp_metrics %>%
  pivot_longer(cols = Kest:Jest,
               names_to = 'metric',
               values_to = 'envelope') %>%
  mutate(extras = map(envelope, ~attr(., 'einfo')),
         numerator = map_dbl(extras, ~2 * .x$nr),
         denominator = map_dbl(extras, ~(.x$nsim + 1)),
         p.value = numerator / denominator) %>%
  select(-extras) %>%
  
  mutate(ind_plot = pmap(list(envelope, Site, numerator, denominator, p.value),
                       ~..1 %>%
                         ggplot(aes(x = r, y = obs, ymin = lo, ymax = hi)) +
                         geom_ribbon(alpha = 0.4) +
                         geom_line() +
                         labs(title = str_c(..2),
                              subtitle = str_c('MC Sig: ', ..3, '/', ..4, ' = ', ..5)) +
                         theme_classic())) %>%
  group_by(metric) %>%
  summarise(metric_plots = list(wrap_plots(ind_plot) + plot_annotation(title = metric)), .groups = 'drop') %>%
  
  mutate(out_name = str_c(OUT_DIR, '/', metric, '.png',sep = '')) %>%
  mutate(write = future_map2(metric_plots, out_name, ~ggsave(.y, plot = .x, height = 15, width = 15),
                             .options = furrr_options(seed = TRUE))) %>%
  select(-write)


message(str_c('Finished All Post-Processing: ', Sys.time()))
