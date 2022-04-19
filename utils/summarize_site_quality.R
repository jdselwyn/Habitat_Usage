library(terra)
library(sf)
library(tidyverse)

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


fish_count <- list.files('~/Coryphopterus/Maps/COPE_Sites/Shoals', pattern = 'shp$', recursive = TRUE, full.names = TRUE) %>%
  tibble(shoals = .) %>%
  mutate(site = str_extract(shoals, 'BZ17-[0-9ABSNK]+')) %>%
  rowwise(site) %>%
  summarise(st_read(shoals, quiet = TRUE), .groups = 'drop') %>%
  filter(Shoal.Size > 0) %>%
  group_by(site) %>%
  summarise(fish = sum(Shoal.Size))


site_quality <- full_join(quality, habitat_type, by = 'site') %>%
  rowwise %>%
  mutate(site_area = global(cellSize(habitat), "sum"),
         simple_mean = global(quality, 'mean', na.rm = TRUE)$mean,
         habitat = list(resample(habitat, quality) %>% round),
         reef = list(classify(habitat == 2, cbind(0, NA))),
         quality = list(mask(quality, reef)),
         upper_tier = list(classify(quality > 0.75, cbind(0, NA)))) %>%
  mutate(just_reef_mean = global(quality, 'mean', na.rm = TRUE)$mean, 
            good_area = global(cellSize(upper_tier), "sum")) %>%
  mutate(percent_good = good_area / site_area * 100) %>%
  ungroup %>%
  arrange(-percent_good)

site_quality$upper_tier[[11]] %>% plot

site_quality %>%
  select(-where(is.list)) %>%
  arrange(just_reef_mean)

library(broom)
full_join(site_quality, fish_count, by = 'site') %>%
  select(-where(is.list)) %>%
  mutate(density = fish / site_area)  %>%
  select(-good_area) %>%
  pivot_longer(cols = -c(site, fish, density, site_area)) %>%
  group_by(name) %>%
  summarise(tidy(MASS::glm.nb(fish ~ value + offset(log(site_area)))),
            .groups = 'drop') %>%
  filter(term == 'value')


tmp_glm <- MASS::glm.nb(fish ~ just_reef_mean + offset(log(site_area)),
    data = full_join(site_quality, fish_count, by = 'site'))

glm_pred <- predict(tmp_glm, newdata = expand_grid(just_reef_mean = seq(0.35, 0.6, length.out = 100), site_area = 1), 
        type = 'link', se.fit = TRUE) %>%
  as_tibble() %>%
  select(-residual.scale) %>%
  mutate(lwr = exp(fit - 1.96 * se.fit),
         upr = exp(fit + 1.96 * se.fit),
         fit = exp(fit)) %>%
  bind_cols(just_reef_mean = seq(0.35, 0.6, length.out = 100))


full_join(site_quality, fish_count, by = 'site') %>%
  select(-where(is.list)) %>%
  mutate(density = fish / site_area) %>%
  
  ggplot(aes(x = just_reef_mean, y = density)) +
  
  geom_ribbon(data = glm_pred, aes(y = fit, ymin = lwr, ymax = upr)) +
  geom_line(data = glm_pred, aes(y = fit)) +
  
  geom_point()
