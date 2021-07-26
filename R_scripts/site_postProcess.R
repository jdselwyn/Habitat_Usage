#### Post-Process many inlabru models ####

if(interactive()){
  model_prefix <- '9.7.21'
  site_choice <- 'BZ17-5KS'
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

message('Total number of attempted model fits: ', nrow(model_fitting))
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

#### Sample from posterior and summarize ####
message('Start Subsample Posterior: ', Sys.time())
samples_to_use <- sample(nlyr(prediction_out2$all_pred[[1]]), 1000)
posteriorSample <- subset(prediction_out2$all_pred[[1]], samples_to_use, 
              filename = '/tmp/posteriorSamples.tif', overwrite = TRUE,
              wopt = list(verbose = TRUE, progress = 1, memfrac = 0.9))

message('Start Summarizing Posterior: ', Sys.time())
r_mean <- app(posteriorSample, 'mean', nodes = 20)
message('Finished Mean ', Sys.time())

r_smin <- app(posteriorSample, 'min', nodes = 20)
message('Finished Min ', Sys.time())

r_smax <- app(posteriorSample, 'max', nodes = 20)
message('Finished Max ', Sys.time())

r_sd <- app(posteriorSample, sd, na.rm = TRUE, nodes = 1)
message('Finished SD ', Sys.time())

r_cv <- r_sd / r_mean
message('Finished CV ', Sys.time())

r_q0.025 <- app(posteriorSample, quantile, probs = c(0.025), na.rm = TRUE, nodes = 1)
message('Finished Qualtile 2.5% ', Sys.time())

r_q0.25 <- app(posteriorSample, quantile, probs = c(0.25), na.rm = TRUE, nodes = 1)
message('Finished Qualtile 25% ', Sys.time())

r_median <- app(posteriorSample, median, na.rm = TRUE, nodes = 1)
message('Finished Median ', Sys.time())

r_q0.75 <- app(posteriorSample, quantile, probs = c(0.75), na.rm = TRUE, nodes = 1)
message('Finished Qualtile 75% ', Sys.time())

r_q0.975 <- app(posteriorSample, quantile, probs = c(0.975), na.rm = TRUE, nodes = 1)
message('Finished Qualtile 97.5% ', Sys.time())


summary_stats <- list(mean = r_mean, sd = r_sd,
                      smin = r_smin, smax = r_smax,
                      median = r_median, q0.025 = r_q0.025,
                      q0.25 = r_q0.25, q0.75 = r_q0.75,
                      q0.975 = r_q0.975, cv = r_cv) %>%
  rast %>%
  set_names(c('mean', 'sd', 'smin', 'smax', 'median', 'q0.025', 'q0.25', 'q0.75', 'q0.975', 'cv'))



message('Start Smoothing files: ', Sys.time())

circle_mat <- raster::focalWeight(raster::raster(summary_stats$mean), sqrt(1/pi), 'circle')
smoothed_summary_stats <- map(names(summary_stats), ~terra::focal(summary_stats[[.x]], w = circle_mat, fun = 'sum') %>%
                                cover(y = summary_stats[[.x]])) %>%
  rast


message('Start Writing files: ', Sys.time())
prediction_out3 <- prediction_out2 %>%
  select(-all_pred) %>%
  mutate(unsmoothed = list(summary_stats),
         smoothed = list(smoothed_summary_stats)) %>%
  pivot_longer(cols = c(unsmoothed, smoothed)) %>%
  mutate(out_name = str_c(OUT_DIR, '/', Site, '_', name, '.tif',
                          sep = '')) %>%
  mutate(write = future_map2(value, out_name, ~writeRaster(.x, .y, overwrite = TRUE),
                             .options = furrr_options(seed = TRUE))) %>%
  select(-write)

message('Finished Post-Processing: ', Sys.time())
file.remove('/tmp/posteriorSamples.tif')

