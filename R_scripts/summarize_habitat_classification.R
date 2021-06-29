## Summarize Habitat Classification ##
args <- commandArgs(trailingOnly=TRUE)
model_choice <- args[1]
# model_choice <- 'c5'

COMPUTER <- 'HPC'
save_suffix <- 'All'

#### Libraries ####
library(sf)
library(raster)
library(terra)
library(rgdal)
# library(gdalUtils)
library(tidyverse)
library(magrittr)
library(tidymodels)
library(workflows)
library(butcher)
library(discrim)
library(rules)
library(baguette)
library(tidyposterior)
library(emmeans)
library(furrr)
library(doFuture)
library(stars)
library(exactextractr)
library(tidytext)
library(whitebox)

#### Set up Computer ####
if(COMPUTER == 'Gawain'){
  # all_cores <- parallel::detectCores(logical = FALSE)
  registerDoFuture()
  cl <- makeClusterPSOCK(10)
  plan('cluster', workers = cl)
  # plan('sequential')
  
  DATA_folder<-'~/Documents/Coryphopterus/Maps' #Gawain
  INTERMEDIATE_FILES<-'~/Documents/Coryphopterus/Habitat Association (Paper Y - PhD Ch. Y)/Intermediate_Files'
  SAVE_LOCATION<-'~/Documents/Coryphopterus/Habitat Association (Paper Y - PhD Ch. Y)/Results'
  
  # options(future.globals.maxSize= 5000*1024^2)
  
} else if(COMPUTER == 'Desktop' | COMPUTER == 'Laptop'){
  # plan(list(multisession, sequential))
  plan(multiprocess)
  
  DATA_folder<-'~/Coryphopterus/Maps' #Gawain
  INTERMEDIATE_FILES<-'~/Coryphopterus/Habitat Association (Paper Y - PhD Ch. Y)/Intermediate_Files'
  SAVE_LOCATION<-'~/Coryphopterus/Habitat Association (Paper Y - PhD Ch. Y)/Results'
} else if(COMPUTER == 'HPC'){
  nodes <- str_extract(Sys.getenv()["SLURM_NODELIST"], '[0-9-,]+') %>%
    str_split(',', simplify = TRUE) %>%
    t %>%
    set_colnames('node_range') %>%
    as_tibble %>%
    mutate(node1 = str_extract(node_range, '^[0-9]+') %>% parse_integer,
           node2 = str_extract(node_range, '[0-9]+$') %>% parse_integer) %>%
    mutate(nodes = map2(node1, node2, ~seq(.x, .y) %>% 
                          str_pad(2, side = 'left', pad = '0') %>%
                          str_c('hpcc', .))) %>%
    pull(nodes) %>%
    unlist
  
  # plan(list(tweak(cluster, workers = nodes), multisession))
  registerDoFuture()
  cl <- makeClusterPSOCK(detectCores())
  plan('cluster', workers = cl)
  
  DATA_folder<-'/work/hobi/jselwyn/Habitat'
  INTERMEDIATE_FILES<-'/work/hobi/jselwyn/Habitat/Intermediate_Files'
  SAVE_LOCATION<-'/work/hobi/jselwyn/Habitat/Results'
}

#### Data Sources ####
DATA_folder<-path.expand(DATA_folder)
INTERMEDIATE_FILES<-path.expand(INTERMEDIATE_FILES)
SAVE_LOCATION<-path.expand(SAVE_LOCATION)
classified_habitat <- if_else(save_suffix == 'justJason', 'Classification/Training Data/Just_Jason', 'Classification/Training Data')

newproj<-"+proj=utm +zone=16Q +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
oldproj <- "+proj=longlat +datum=WGS84 +no_defs"

all_site_locations <- tibble(training_file = list.files(path = str_c(DATA_folder, classified_habitat, sep='/'), pattern = 'shp$', full.names = TRUE),
                             shoals_file = list.files(str_c(DATA_folder,'COPE_Sites/Shoals',sep='/'), recursive = TRUE, pattern = 'shp$', full.names = TRUE)) %>%
  mutate(Site = str_extract(shoals_file, 'BZ17-[0-9]*.[ABKNS]')) %>%
  mutate(unsmoothed = str_c(INTERMEDIATE_FILES, '/Habitat_Classification/', Site, "_", save_suffix, '_', model_choice, '_New_habitat_class.tif',sep=''),
         smoothed = str_c(INTERMEDIATE_FILES, '/Habitat_Classification/', Site, "_", save_suffix, '_', model_choice, '_New_habitat_class_smoothed.tif',sep='')) %>%
  dplyr::select(Site, everything()) %>%
  filter(file.exists(unsmoothed),
         file.exists(smoothed))

#### Plot classifications ####

pdf(str_c(SAVE_LOCATION, '/Habitat_Classification/', save_suffix, '_', model_choice, '_classified_reefs.pdf',sep=''),
    onefile = TRUE)
for(i in 1:nrow(all_site_locations)){
  raster(all_site_locations$unsmoothed[i]) %>%
    plot(main = all_site_locations$Site[i])
}
dev.off()

pdf(str_c(SAVE_LOCATION, '/Habitat_Classification/', save_suffix, '_', model_choice, 
          '_classified_reefs_smoothed.pdf',sep=''),
    onefile = TRUE)
for(i in 1:nrow(all_site_locations)){
  raster(all_site_locations$smoothed[i])  %>%
    plot(main = all_site_locations$Site[i])
}
dev.off()


#### Evaluate training shapes ####

eval_met <- metric_set(sens, spec, accuracy, bal_accuracy, f_meas, j_index, kap, detection_prevalence, mcc, npv, ppv, precision, recall)

train_data <- all_site_locations %>%
  select(Site, training_file, unsmoothed, smoothed) %>%
  
  mutate(training = map(training_file, ~st_read(.x, quiet = TRUE) %>%
                          st_transform(crs = newproj) %>%
                          filter(classStr!='flag',classStr!='GCP') %>%
                          mutate(class=if_else(classStr=='sand','C1','C2')))) %>%
  mutate(across(c(smoothed, unsmoothed), ~map2(.x, training, ~raster(.x) %>%
                                                 exact_extract(.y, 'majority', progress = FALSE)))) %>%
  unnest(c(smoothed, unsmoothed, training)) %>%
  select(-training_file) %>%
  mutate(across(c(smoothed, unsmoothed), ~str_c('C', .)),
         across(c(class, smoothed, unsmoothed), ~factor(., levels = c('C1', 'C2')))) 

train_eval <- train_data %>%
  select(-Site, -geometry, -classStr) %>%
  nest(smoothed = c(smoothed, class),
       unsmoothed = c(unsmoothed, class)) %>%
  mutate(Site = 'overall') %>%
  mutate(across(c(smoothed, unsmoothed), ~map(.x, ~rename(.x, est = 1))),
         across(c(smoothed, unsmoothed), ~map(.x, ~eval_met(.x, truth = class, estimate = est)))) %>%
  mutate(combined = map2(smoothed, unsmoothed, ~full_join(.x, .y, 
                                                          by = c('.estimator', '.metric'),
                                                          suffix = c('.smoothed', '.unsmoothed')))) %>%
  select(Site, combined) %>%
  bind_rows(train_data %>%
              group_by(Site) %>%
              select(-geometry, -classStr) %>%
              nest(smoothed = c(smoothed, class),
                   unsmoothed = c(unsmoothed, class)) %>%
              mutate(across(c(smoothed, unsmoothed), ~map(.x, ~rename(.x, est = 1))),
                     across(c(smoothed, unsmoothed), ~map(.x, ~eval_met(.x, truth = class, estimate = est)))) %>%
              mutate(combined = map2(smoothed, unsmoothed, ~full_join(.x, .y, 
                                                                      by = c('.estimator', '.metric'),
                                                                      suffix = c('.smoothed', '.unsmoothed')))) %>%
              select(Site, combined) 
  ) %>%
  unnest(combined) %T>%
  write_csv(str_c(SAVE_LOCATION, '/Habitat_Classification/', save_suffix, '_', model_choice, 
                  '_eval_metrics_shapes.csv',sep=''))
