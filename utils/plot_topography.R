library(tidyverse)
library(sf)
library(stars)
library(viridis)
library(patchwork)

vars_to_keep <- c('depth', 'moran', 
                  'viewshed', 'boundaryDistance',
                  'rugosity', 'vectorDispersion',
                  'relief', 'slope',
                  'classified_habitat_unsmoothed')

dir_to_look_in <- '/work/hobi/jselwyn/Habitat/Intermediate_Files/Topography'

site_rename <- tibble(Site = str_c('BZ17-', c('0A', '0B', '100N', '100S', '10KN', '1KN',
                                              '1KS', '500N', '500S', '5KN', '5KS', '60S')),
                      site = c('F', 'G', 'E', 'I', 'A', 'C', 'K', 'D', 'J', 'B', 'L', 'H'))



all_topography <- list.files(dir_to_look_in, pattern = 'stack.*tif$', full.names = TRUE) %>%
  tibble(file = .) %>%
  mutate(Site = str_extract(file, 'BZ17-[0-9ABNSK]+')) %>%
  rowwise %>%
  mutate(names = list(terra::rast(file) %>% names)) %>%
  mutate(topography = list(terra::rast(file) %>% 
                             terra::trim() %>% 
                             raster::stack() %>%
                             st_as_stars)) %>%
  mutate(topography = list(topography %>%
                             split('band') %>%
                             set_names(names) %>%
                             select(one_of(vars_to_keep)) %>%
                             mutate(depth = depth * -1,
                                    moran = moran * -1,
                                    classified_habitat_unsmoothed = if_else(classified_habitat_unsmoothed == 1, 
                                                                            'Sand', 'Reef')))) %>%
  select(-file, -names) %>%
  left_join(site_rename, by = 'Site') %>%
  arrange(site)



all_plots <- all_topography %>%
  mutate(plot_names = list(vars_to_keep)) %>%
  rowwise %>%
  mutate(all_plots = list(map(plot_names, ~ggplot() +
                                geom_stars(data = select(topography, one_of(.x))) +
                                coord_sf() +
                                theme_void() +
                                scale_x_discrete(expand = c(0, 0)) +
                                scale_y_discrete(expand = c(0, 0)) +
                                labs(fill = NULL) +
                                theme(plot.background = element_blank(),
                                      panel.background = element_blank())
                              
                              ))) %>%
  ungroup %>%
  select(-topography) %>%
  unnest(c(plot_names, all_plots)) %>%
  mutate(plot_titles = case_when(plot_names == 'moran' ~ 'Coarse Complexity',
                                plot_names == 'boundaryDistance' ~ 'Distance to Sand/Reef Margin',
                                plot_names == 'classified_habitat_unsmoothed' ~ 'Habitat Type',
                                plot_names == 'vectorDispersion' ~ 'Vector Dispersion',
                                TRUE ~ str_to_sentence(plot_names))) %>%
  rowwise %>%
  mutate(all_plots = list(all_plots + labs(title = plot_titles)),
         all_plots = if_else(plot_names == 'classified_habitat_unsmoothed',
                             list(all_plots + scale_fill_viridis(na.value = "white", discrete = TRUE, 
                                                                 na.translate = FALSE)),
                             list(all_plots + scale_fill_viridis(na.value = "white", discrete = FALSE))))



wrapped_plots <- all_plots %>%
  group_by(Site, site) %>%
  summarise(plot = list(wrap_plots(all_plots) + plot_annotation(title = site)),
            .groups = 'drop') %>%
  mutate(file_out = str_c('/work/hobi/jselwyn/Habitat/Results/Topography/', Site, '.png')) %>%
  ungroup %>%
  mutate(write_out = map2(plot, file_out, ~ggsave(filename = .y, plot = .x, height = 15, width = 15))) %>%
  select(-write_out)

