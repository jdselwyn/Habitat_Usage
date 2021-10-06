library(terra)
library(sf)
library(tidyverse)
library(emmeans)
library(broom)
library(patchwork)
library(brms)
library(metR)
library(lme4)
library(lmerTest)
library(RColorBrewer)


#### Read in data ####
prediction <- list.files('~/Coryphopterus/Habitat Association (Paper Y - PhD Ch. Y)/Results/INLA_models/model_Overall.9.7.21', 
                         pattern = 'tif$', full.names = TRUE) %>%
  str_subset('unsmoothed') %>%
  tibble(predictions = .) %>%
  mutate(site = str_extract(predictions, 'BZ17-[0-9KNSAB]+')) %>%
  group_by(site) %>%
  summarise(predictions = list(rast(predictions)))

quality <- list.files('~/Coryphopterus/Dispersal/Data/', pattern = 'tif', full.names = TRUE) %>%
  tibble(quality = .) %>%
  mutate(site = str_extract(quality, 'BZ17-[0-9KNSAB]+')) %>%
  group_by(site) %>%
  summarise(quality = list(rast(quality)))

habitat_type <- list.files('~/Coryphopterus/Habitat Association (Paper Y - PhD Ch. Y)/Intermediate_Files/Topography/',
                           pattern = 'tif', full.names = TRUE) %>%
  tibble(habitat = .) %>%
  mutate(site = str_extract(habitat, 'BZ17-[0-9KNSAB]+')) %>%
  group_by(site) %>%
  summarise(habitat = list(rast(habitat)$classified_habitat_unsmoothed))


shoals <- list.files('~/Coryphopterus/Maps/COPE_Sites/Shoals', pattern = 'shp$', recursive = TRUE, full.names = TRUE) %>%
  tibble(shoals = .) %>%
  mutate(site = str_extract(shoals, 'BZ17-[0-9ABSNK]+')) %>%
  rowwise(site) %>%
  summarise(st_read(shoals, quiet = TRUE), .groups = 'drop') %>%
  select(site, Shoal.Size, geometry) %>%
  filter(Shoal.Size > 0) %>%
  nest(data = -site) %>%
  rowwise(site) %>%
  summarise(shoals = list(st_as_sf(data) %>% as('Spatial') %>% vect),
            .groups = 'drop')

#### Make Data agree on resolution and all be rasters ####
full_data <- full_join(quality, prediction, by = 'site') %>%
  full_join(habitat_type, by = 'site') %>%
  full_join(shoals, by = 'site') %>%
  mutate(site = factor(site, levels = str_c('BZ17', c('10KN', '5KN', '1KN',
                                                      '500N', '100N','0A',
                                                      '0B', '60S', '100S',
                                                      '500S', '1KS', '5KS'),
                                            sep = '-')),
         site = LETTERS[as.integer(site)]) %>%
  rowwise %>%
  mutate(quality = list(trim(quality)),
         areas = list(cellSize(quality)),
         
         predictions = list(trim(predictions)),
         predictions = list(predictions * areas), #Make # fish not fish/area
         shoals = list(project(shoals, quality)),
         shoals = list(rasterize(shoals, quality, field = 'Shoal.Size')),
         site_area = global(cellSize(quality), "sum", na.rm = TRUE)$sum,
         habitat = list(resample(habitat, quality) %>% round))

#### Create Zones - split up each raster into grid ####
grid_zoner <- function(r, x, y){
  z <- r
  z <- raster::raster(z)
  
  Y <- nrow(r)
  X <- ncol(r)
  
  Y_cuts <- c(1, round((1:y) * (Y / y)))
  X_cuts <- c(1, round((1:x) * (X / x)))
  
  #Almost certainly a faster way to do this - either recursion or premake the grid and use apply maybe?
  counter <- 0
  for(i in 1:(length(Y_cuts) - 1)){
    for(j in 1:(length(X_cuts) - 1)){
      counter <- counter + 1
      z[Y_cuts[i]:Y_cuts[i + 1], X_cuts[j]:X_cuts[j + 1]] <- counter
    }
  }
  
  names(z) <- "zone"
  z <- rast(z)
  
  
  mask(z, r)
}

zones <- full_data %>%
  rowwise %>%
  mutate(zones = list(grid_zoner(quality, 5, 5)))

zone_stats <- zones %>%
  rowwise(site) %>%
  summarise(zonal(areas, zones, sum), # %>% rename(area = median)
            zonal(quality, zones, 'mean', na.rm = TRUE) %>%
              rename(quality = median),
            zonal(predictions, zones, 'sum', na.rm = TRUE),
            zonal(shoals, zones, 'sum', na.rm = TRUE) %>% 
              rename(fish = median),
            .groups = 'drop') 

zone_stats %>%
  # mutate(across(c(median, fish, starts_with('q')), ~./area)) %>%
  ggplot(aes(x = median, y = fish)) +
  geom_abline(intercept = 0, slope = 1, 
              colour = 'black', linetype = 'dashed') +
  geom_smooth(method = 'lm', formula = y ~ 0 + x) +
  geom_point(aes(colour = site))


anova(lmer(mean ~ fish + (0 + fish | site), data = zone_stats),
      lm(mean ~ fish, data = zone_stats))

summary(lmer(mean ~ fish + (0 + fish | site), data = zone_stats))

lmer(mean ~ fish + (0 + fish | site), data = zone_stats) %>%
  emtrends(~1, var = 'fish', infer = TRUE, null = 1)

red_blue_colours <- rev(colorRampPalette(brewer.pal(9, 'RdBu'))(12))

lmer(mean ~ fish + (0 + fish | site), data = zone_stats) %>%
  emmeans(~fish, at=list(fish = modelr::seq_range(zone_stats$fish, n = 100))) %>%
  tidy(conf.int = TRUE) %>%
  
  ggplot(aes(x = fish, y = estimate)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.5) +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, 
              colour = 'black', linetype = 'dashed') +
  geom_linerange(data = zone_stats, aes(ymin = q0.25, ymax = q0.75, 
                                         x = fish, y = mean)) +
  geom_point(data = zone_stats, aes(x = fish, y = mean,
                                    fill = site),
             pch = 21, size = 4) +
  # scale_x_continuous(trans = scales::pseudo_log_trans()) +
  # scale_y_continuous(trans = scales::pseudo_log_trans()) +
  scale_fill_manual(values = red_blue_colours) +
  guides(fill = guide_legend(override.aes = list(size = 4))) +
  labs(y = "Modelled # Fish",
       x = "Observed # Fish",
       fill = 'Site') +
  theme_classic() +
  theme(axis.text = element_text(size = 14, colour = 'black'),
        axis.title = element_text(size = 16),
        legend.position = c(0.1, 0.68),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        strip.background = element_blank(),
        strip.text = element_text(size = 20, hjust = 0),
        panel.border = element_rect(colour = 'black', size = 1, fill = NA))
ggsave('../../Manuscript/Figures/Figure 4.svg', height = 6, width = 6)

chisq.test(zone_stats$fish, 
           p = zone_stats$mean, 
           rescale.p = TRUE, 
           simulate.p.value = TRUE,
           B = 10000)

chisq.test(zone_stats$fish, 
           p = zone_stats$mean, 
           rescale.p = TRUE, 
           simulate.p.value = FALSE)

#### Site Stats ####
site_data <- full_data %>%
  rowwise(site, site_area) %>%
  mutate(reef = list(classify(habitat == 2, cbind(0, NA))),
         quality_reef = list(mask(quality, reef))) %>%
  summarise(global(predictions, 'sum', na.rm = TRUE) %>%
              as_tibble(rownames = 'parameter') %>%
              pivot_wider(names_from = 'parameter',
                          values_from = 'sum', 
                          names_prefix = 'pred_'),
            quality_overall_mean = global(quality, 'mean', na.rm = TRUE)$mean,
            quality_overall_sd = global(quality, 'sdpop', na.rm = TRUE)$sd,
            quality_reef_mean = global(quality_reef, 'mean', na.rm = TRUE)$mean,
            quality_reef_sd = global(quality_reef, 'sdpop', na.rm = TRUE)$sd,
            fish = global(shoals, 'sum', na.rm = TRUE)$sum,
            .groups = 'drop') %>%
  mutate(quality_overall_cv = quality_overall_sd / quality_overall_mean,
         quality_reef_cv = quality_reef_sd / quality_reef_mean)

chisq.test(site_data$fish, 
           p = site_data$pred_mean, 
           rescale.p = TRUE, 
           simulate.p.value = FALSE)


quality_glm <- brm(fish ~ I(quality_overall_mean - mean(quality_overall_mean)) *
                     I(quality_overall_sd - mean(quality_overall_sd)) + 
                     (1 | site) +
                     offset(log(site_area)),
                   family = 'poisson',
                   data = site_data,
                   backend = 'cmdstanr',
                   chains = 4,
                   cores = 2,
                   seed = 1234)
plot(quality_glm)
pp_check(quality_glm)
bayes_R2(quality_glm)
summary(quality_glm)
hypothesis(quality_glm, 
           c('mean' = "Iquality_overall_meanMmeanquality_overall_mean > 0",
             'sd' = "Iquality_overall_sdMmeanquality_overall_sd > 0",
             'interaction' = 'Iquality_overall_meanMmeanquality_overall_mean:Iquality_overall_sdMmeanquality_overall_sd < 0'))


mean_plot <- emmeans(quality_glm, ~quality_overall_mean | quality_overall_sd, type = 'response', offset = log(1),
        at=list(quality_overall_mean = modelr::seq_range(site_data$quality_reef_mean, n = 100))) %>%
  tidy(conf.int = TRUE) %>%
  
  ggplot(aes(x = quality_overall_mean, y = rate)) +
  geom_ribbon(aes(ymin = lower.HPD, ymax = upper.HPD), alpha = 0.2) +
  geom_line() +
  geom_point(data = site_data, aes(y = fish / site_area), colour = 'black', fill = 'black') +
  labs(x = 'Mean Reef Quality',
       y = expression(paste("# Fish (", m^-2, ')'))) +
  theme_classic() +
  theme(axis.text = element_text(size = 14, colour = 'black'),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        strip.background = element_blank(),
        strip.text = element_text(size = 20, hjust = 0),
        panel.border = element_rect(colour = 'black', size = 1, fill = NA))


sd_plot <- emmeans(quality_glm, ~quality_overall_sd | quality_overall_mean, type = 'response', offset = log(1),
        at=list(quality_overall_sd = modelr::seq_range(site_data$quality_overall_sd, n = 100))) %>%
  tidy(conf.int = TRUE) %>%
  
  ggplot(aes(x = quality_overall_sd, y = rate)) +
  geom_ribbon(aes(ymin = lower.HPD, ymax = upper.HPD), alpha = 0.2) +
  geom_line() +
  geom_point(data = site_data, aes(y = fish / site_area), colour = 'black', fill = 'black') +
  labs(x = 'Standard Deviation Reef\nQuality',
       y = expression(paste("# Fish (", m^-2, ')'))) +
  theme_classic() +
  theme(axis.text = element_text(size = 14, colour = 'black'),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        strip.background = element_blank(),
        strip.text = element_text(size = 20, hjust = 0),
        panel.border = element_rect(colour = 'black', size = 1, fill = NA))


contour_plot <- emmeans(quality_glm, ~quality_overall_sd | quality_overall_mean, type = 'response', offset = log(1),
                        at=list(quality_overall_sd = modelr::seq_range(site_data$quality_overall_sd, n = 50),
                                quality_overall_mean = modelr::seq_range(site_data$quality_reef_mean, n = 50))) %>%
  tidy(conf.int = FALSE) %>%
  ggplot(aes(x = quality_overall_mean, y = quality_overall_sd, z = rate)) +
  geom_contour_fill() +
  scale_fill_binned(type = "viridis", alpha = 0.75, breaks = seq(0, 5, by = 1)) +
  geom_contour(colour = 'black') +
  geom_label_contour(label.placer = label_placer_random()) +
  geom_point(data = site_data,
             aes(x = quality_overall_mean, y = quality_overall_sd,
                 fill = fish / site_area),
             pch = 21, size = 5, alpha = 1,
             inherit.aes = FALSE) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  guides(fill = guide_colourbar(override.aes = list(alpha = 1, shape = NA))) +
  labs(x = 'Mean Reef Quality',
       y = 'Standard Deviation\nReef Quality',
       fill = expression(paste("# Fish (", m^-2, ')'))) +
  theme_classic() +
  theme(axis.text = element_text(size = 14, colour = 'black'),
        axis.title = element_text(size = 16),
        legend.position = 'bottom',
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        strip.background = element_blank(),
        strip.text = element_text(size = 20, hjust = 0),
        panel.border = element_rect(colour = 'black', size = 1, fill = NA))

(mean_plot + sd_plot) / contour_plot + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 20))
ggsave('../../Manuscript/Figures/Figure 6.svg', height = 6, width = 7)
