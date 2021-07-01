## Calculate site stats ###
library(terra)
library(sf)
library(tidyverse)
library(magrittr)
library(yardstick)
library(car)
library(performance)

select <- dplyr::select


#### Shoal size differences ####
shoals <- list.files("~/Coryphopterus/Maps/COPE_Sites/Shoals", pattern = '.shp$', recursive = TRUE, full.names = TRUE) %>%
  tibble(file = .) %>%
  mutate(Site = str_extract(file, 'BZ17-[0-9ABNSK]+'),
         Site = factor(Site, levels = str_c('BZ17', c('10KN', '5KN', '1KN', '500N', '100N', 
                                                      '0A', '0B', '60S', '100S', '500S', '1KS', '5KS'),
                                            sep = '-'))) %>%
  mutate(shoals = map(file, st_read, quiet = TRUE)) %>%
  select(-file) %>%
  unnest(shoals) %>%
  select(Site, Shoal.Size) %>%
  filter(Shoal.Size > 0)

which.max(table(shoals$Shoal.Size))

shoals %>%
  summarise(n_shoal = n(),
            max_shoal = max(Shoal.Size),
            mean_shoal = mean(Shoal.Size),
            sd_shoal = sd(Shoal.Size),
            se_shoal = sd_shoal / sqrt(n_shoal),
            median_shoal = median(Shoal.Size),
            modal_shoal = names(which.max(table(Shoal.Size))),
            .groups = 'drop')


site_shoals <- shoals %>%
  group_by(Site) %>%
  summarise(n_shoal = n(),
            cope = sum(Shoal.Size),
            max_shoal = max(Shoal.Size),
            mean_shoal = mean(Shoal.Size),
            sd_shoal = sd(Shoal.Size),
            se_shoal = sd_shoal / sqrt(n_shoal),
            median_shoal = median(Shoal.Size),
            modal_shoal = names(which.max(table(Shoal.Size))),
            .groups = 'drop') %>%
  arrange(Site) %>%
  mutate(across(c(mean_shoal, sd_shoal), ~round(., 2))) %T>%
  write_csv('tmp_table.csv')


aov(log(Shoal.Size) ~ Site, data = shoals) %>% check_model()
aov(log(Shoal.Size) ~ Site, data = shoals) %>% anova


shoal_size <- MASS::glm.nb(Shoal.Size ~ Site, data = shoals)

Anova(shoal_size, test="LR")
check_model(shoal_size)

library(emmeans)
library(multcomp)
emmeans(shoal_size, ~ Site, type = 'link') %>%
  contrast('pairwise')

emmeans(shoal_size, ~ Site, type = 'response') %>%
  as_tibble() %>%
  full_join(emmeans(shoal_size, ~ Site, type = 'link') %>%
              cld(alpha = 0.1, Letters = LETTERS) %>%
              as_tibble %>%
              dplyr::select(Site, .group),
            by = 'Site') %>%
  mutate(.group = str_trim(.group)) %>%
  ggplot(aes(x = response, xmin = asymp.LCL, xmax = asymp.UCL, y = Site)) +
  
  geom_rect(data = as_tibble(emmeans(shoal_size, ~ 1, type = 'response')),
            aes(ymin = -Inf, ymax = Inf, xmin = asymp.LCL, xmax = asymp.UCL),
            inherit.aes = FALSE) +
  geom_vline(data = as_tibble(emmeans(shoal_size, ~ 1, type = 'response')),
             aes(xintercept = response), linetype = 'dashed') +
  
  geom_pointrange() +
  geom_text(aes(label = .group), vjust = -0.5, hjust = 0.5) +
  scale_y_discrete(limits = rev) +
  labs(x = 'Average Shoal Size',
       y = NULL) +
  theme_classic()

#### DEM RMSE ####
error_measure <- read_csv('~/Coryphopterus/Maps/COPE_Sites/DEM/Error meassurements.csv') %>%
  mutate(Site = str_c('BZ17', Site, sep = '-'),
         Site = factor(Site, levels = str_c('BZ17', c('10KN', '5KN', '1KN', '500N', '100N', 
                                                      '0A', '0B', '60S', '100S', '500S', '1KS', '5KS'),
                                            sep = '-'))) %>%
  mutate(Slate = str_replace(Slate, '10', 'x'))


slates <- read_csv('~/Coryphopterus/Maps/COPE_Sites/DEM/Slate Measurements.csv', skip = 1, col_names = c('Slate', 'L', 'W')) %>%
  mutate(Slate = str_to_lower(Slate),
         Slate = str_remove(Slate, ' pipe')) %>%
  pivot_longer(cols = -Slate,
               names_to = 'type',
               values_to = 'true_length') %>%
  # filter(str_detect(Slate, 'pvc')) %>%
  mutate(type = if_else(Slate == 'pvc middle section' & type == 'L', 'M', type)) %>%
  filter(!(str_detect(Slate, 'pvc') & type == 'W')) %>%
  mutate(Slate = str_remove(Slate, ' middle section'))

error_measure %>%
  left_join(slates) %$%
  rmse_vec(true_length, length)

error_measure %>%
  left_join(slates) %>%
  group_by(Site) %>%
  summarise(rmse = rmse_vec(true_length, length)) %T>%
  write_csv('tmp_table.csv') %>%
  summarise(rmse_mean = mean(rmse),
            rmse_sd = sd(rmse),
            rmse_se = rmse_sd/sqrt(n()))


error_measure %>%
  left_join(slates) %>%
  group_by(Site) %>%
  summarise(rmse = rmse_vec(true_length, length)) %>%
  lm(rmse ~ 1, data = .) %>%
  summary

error_measure %>%
  left_join(slates) %>%
  mutate(difference = sqrt((length - true_length)^2)) %>%
  lm(difference ~ Site, data = .) %>%
  # anova
  emmeans(~Site) %>%
  contrast('pairwise') %>%
  as_tibble() %>%
  filter(p.value < 0.05)

#### Overall Site ####
#To convert radians to grade - grade = rad * 63.661977 - that gives # cm / m

site_metrics <- list.files("~/Coryphopterus/Maps/COPE_Sites/DEM", pattern = '.tif$', full.names = TRUE) %>%
  tibble(file = .) %>%
  mutate(Site = str_extract(file, 'bz17-[0-9ABNSK]+') %>% str_to_upper,
         Site = factor(Site, levels = str_c('BZ17', c('10KN', '5KN', '1KN', '500N', '100N', 
                                                      '0A', '0B', '60S', '100S', '500S', '1KS', '5KS'),
                                            sep = '-'))) %>%
  mutate(dem = map(file, rast)) %>%
  select(-file) %>%
  mutate(min_depth = map(dem, global, fun = 'min', na.rm = TRUE),
         max_depth = map(dem, global, fun = 'max', na.rm = TRUE)) %>%
  mutate(area = map_dbl(dem, expanse)) %>%
  
  mutate(across(c(min_depth, max_depth), ~map_dbl(., as.numeric)),
         across(c(min_depth, max_depth), ~-1 * abs(.))) %>%
  mutate(relief = abs(max_depth - min_depth)) %>%
  select(-dem) %>%
  pivot_longer(cols = ends_with('depth'),
               values_to = 'depth') %>%
  group_by(Site, relief, area) %>%
  summarise(depth = min(depth), .groups = 'drop') %>%
  ungroup %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  select(Site, depth, relief, area) %T>%
  write_csv('tmp_table.csv')

site_metrics %>%
  summarise(across(depth:area, list(mean = mean, se = ~sd(.)/sqrt(12), min = min, max = max))) %>%
  pivot_longer(cols = everything(),
               names_to = c('metric', '.value'),
               names_pattern = '(.*)_(.*)')


#### Site habitat metrics - for model ####
site_level_data <- read_csv('~/Coryphopterus/Habitat Association (Paper Y - PhD Ch. Y)/Results/Topography/site_metrics_c5.csv') %>%
  mutate(Site = factor(Site, levels = str_c('BZ17', c('10KN', '5KN', '1KN', '500N', '100N', 
                                                      '0A', '0B', '60S', '100S', '500S', '1KS', '5KS'),
                                            sep = '-'))) %>%
  select(-contains('_smoothed')) 


site_level_data %>%
  select(Site, classified_habitat_unsmoothed.prop_sand, classified_habitat_unsmoothed.prop_reef) %>%
  mutate(across(contains('unsmoothed'), ~.*100)) %>%
  summarise(across(contains('unsmoothed'), list(mean = mean, se = ~sd(.)/sqrt(12)))) %>%
  pivot_longer(cols = everything(),
               names_to = c('metric', '.value'),
               names_pattern = '(.*)_(mean|se)') %>%
  mutate(metric = str_extract(metric, 'sand|reef'))


site_level_data %>%
  pivot_longer(cols = -c(Site, contains('prop'), contains('globalMoran')), 
               names_to = c('metric', '.value'), 
               names_pattern = '(.*)\\.(.*)') %>%
  filter(metric %in% c('depth', 'slope', 'moran', 'relief', 'viewshed', 
                       'boundaryDistance', 'rugosity', 'vectorDispersion',
                       'rayShade', 'PC1')) %>%
  select(Site, metric, mean, sd) %>%
  mutate(across(c(mean, sd), ~case_when(metric == 'relief' ~ .*100,
                                        metric == 'depth' ~ abs(.),
                                        metric == 'rayShade' ~.*100,
                                        metric == 'viewshed' ~.*100,
                                        metric == 'slope' ~ . * 180/pi,
                                        TRUE ~ .))) %>%
  group_by(metric) %>%
  summarise(mean_metric = mean(mean),
            se_metric = sd(mean)/sqrt(12),
            mean_sd = mean(sd),
            se_sd = sd(sd)/sqrt(12),
            .groups = 'drop')


site_level_data %>%
  pivot_longer(cols = -c(Site, contains('prop'), contains('globalMoran')), 
               names_to = c('metric', '.value'), 
               names_pattern = '(.*)\\.(.*)') %>%
  filter(metric %in% c('depth', 'slope', 'moran', 'relief', 'viewshed', 
                       'boundaryDistance', 'rugosity', 'vectorDispersion',
                       'rayShade', 'PC1')) %>%
  mutate(across(c(mean, sd), ~case_when(metric == 'relief' ~ .*100,
                                        metric == 'depth' ~ abs(.),
                                        metric == 'rayShade' ~.*100,
                                        metric == 'viewshed' ~.*100,
                                        metric == 'slope' ~ . * 180/pi,
                                        TRUE ~ .))) %>%
  mutate(contents = str_c(round(mean, 3), ' +- ', round(sd, 3))) %>%
  select(-mean:-sd) %>%
  pivot_wider(names_from = 'metric',
              values_from = 'contents') %>%
  mutate(sand = as.character(round(100 * classified_habitat_unsmoothed.prop_sand, 1)),
         reef = as.character(round(100 * classified_habitat_unsmoothed.prop_reef, 1)),
         # globalMoran = round(globalMoran, 3),
         ) %>%
  select(-starts_with('classified')) %>%
  arrange(Site) %>%
  write_csv('tmp_table.csv')
