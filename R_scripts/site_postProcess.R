#### Post-Process many inlabru models ####


if(interactive()){
  model_prefix <- '8.7.21.test'
  site_choice <- 'BZ17-5KN'
} else {
  args <- commandArgs(trailingOnly=TRUE)
  model_prefix <- args[1]
  site_choice <- args[2]
}

message(model_prefix)
message(site_choice)

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
plan('multicore')

#### Get Results ####
RESULTS_DIR <- '/work/hobi/jselwyn/Habitat/Results/INLA_models'
OUT_DIR <- str_c(RESULTS_DIR, '/models_Overall.', model_prefix)

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

message('Total number of attempted model fits: ', max(model_fitting$model_number) + 1)
message('Total number of successful model fits: ', sum(model_fitting$model_finished))

model_folders <- model_fitting %>%
  filter(model_finished) %>%
  rename(version = model_number) %>%
  select(-model_finished) %>%
  # sample_n(25) %>% #for testing
  identity()

#### Merge Prediction Rasters ####
message('Start reading files: ', Sys.time())
prediction_out <- model_folders %>%
  mutate(file = map(model_folders, ~list.files(.x, pattern = 'tif$', recursive = FALSE, full.names = TRUE) %>%
                      str_subset(model_prefix) %>%
                      str_subset('unsmoothed') %>%
                      str_subset('noSpatial', negate = TRUE))) %>%
  unnest(file) %>%
  mutate(Site = str_extract(file, 'BZ17-[0-9ABKNS]+')) %>%
  left_join(site_rename, by = c('Site')) %>%
  filter(Site == site_choice) %>%
  rowwise() %>%
  mutate(preds = list(rast(file))) 


message('Start Merging files: ', Sys.time())
prediction_out2 <- prediction_out %>%
  group_by(Site, site) %>%
  summarise(all_pred = list(rast(preds)), .groups = 'rowwise') 


message('Start Summarizing files: ', Sys.time())
prediction_out3 <- prediction_out2 %>%
  mutate(mean = list(app(all_pred, 'mean', nodes = 20)),
         sd = list(app(all_pred, sd, na.rm = TRUE, nodes = 20)),
         smin = list(app(all_pred, 'min', nodes = 20)),
         smax = list(app(all_pred, 'max', nodes = 20)),
         median = list(app(all_pred, median, na.rm = TRUE, nodes = 20)),
         q0.025 = list(app(all_pred, quantile, probs = c(0.025), na.rm = TRUE, nodes = 20)),
         q0.25 = list(app(all_pred, quantile, probs = c(0.25), na.rm = TRUE, nodes = 20)),
         q0.75 = list(app(all_pred, quantile, probs = c(0.75), na.rm = TRUE, nodes = 20)),
         q0.975 = list(app(all_pred, quantile, probs = c(0.975), na.rm = TRUE, nodes = 20)),
         cv = list(sd / mean)) %>%
  ungroup %>%
  select(-all_pred) %>%
  pivot_longer(cols = -c(Site, site),
               names_to = 'metric',
               values_to = 'unsmoothed') %>%
  rowwise %>%
  mutate(unsmoothed = list(set_names(unsmoothed, metric))) %>%
  group_by(Site, site) %>%
  summarise(unsmoothed = list(rast(unsmoothed)), .groups = 'drop')


message('Start Smoothing files: ', Sys.time())
prediction_out4 <- prediction_out3 %>%
  
  mutate(circle_mat = future_map(unsmoothed, 
                                 ~raster::focalWeight(raster::raster(.x$mean), sqrt(1/pi), 'circle'),
                                 .options = furrr_options(seed = TRUE))) %>%
  rowwise %>%
  mutate(smoothed = list(map(names(unsmoothed), ~terra::focal(unsmoothed[[.x]], w = circle_mat, fun = 'sum') %>%
                               cover(y = unsmoothed[[.x]])) %>%
                           rast
  )) %>%
  ungroup %>%
  select(-circle_mat) %>%
  pivot_longer(cols = c(unsmoothed, smoothed)) 



message('Start Writing files: ', Sys.time())
prediction_out5 <- prediction_out4 %>%
  mutate(out_name = str_c(OUT_DIR, '/', Site, '_', name, '.tif',
                          sep = '')) %>%
  mutate(write = future_map2(value, out_name, ~writeRaster(.x, .y, overwrite = TRUE),
                             .options = furrr_options(seed = TRUE))) %>%
  select(-write)

message('Finished Post-Processing: ', Sys.time())
