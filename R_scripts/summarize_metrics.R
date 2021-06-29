#### Get mean/sd of parameters for all sites together to z-score transform sites for inlabru model ####

args <- commandArgs(trailingOnly=TRUE)
model_choice <- args[1]
# model_choice <- 'c5'
##TODO why aren't there any negative z-scores for viewshed - check the rest?

#### Libraries ####
library(raster)
library(terra)
library(tidyverse)
library(magrittr)
library(furrr)

#### Set up Computer ####
save_suffix <- 'All'
COMPUTER <- if_else(Sys.info()['nodename'] == 'JASONDELL', 'laptop', 'HPC')
PROG <- if_else(COMPUTER == 'HPC', FALSE, TRUE)

if(COMPUTER == 'Gawain'){
  DATA_folder<-'~/Documents/Coryphopterus/Maps' #Gawain
  INTERMEDIATE_FILES<-'~/Documents/Coryphopterus/Habitat Association (Paper Y - PhD Ch. Y)/Intermediate_Files'
  SAVE_LOCATION<-'~/Documents/Coryphopterus/Habitat Association (Paper Y - PhD Ch. Y)/Results'
  
} else if(COMPUTER == 'laptop'){
  DATA_folder<-'~/Coryphopterus/Maps' #Gawain
  INTERMEDIATE_FILES<-'~/Coryphopterus/Habitat Association (Paper Y - PhD Ch. Y)/Intermediate_Files'
  SAVE_LOCATION<-'~/Coryphopterus/Habitat Association (Paper Y - PhD Ch. Y)/Results'
  
} else if (COMPUTER == 'HPC'){
  plan('multicore')
  
  options(future.globals.maxSize= 5000*1024^2) #number of mb
  
  DATA_folder<-'/work/hobi/jselwyn/Habitat'
  INTERMEDIATE_FILES<-'/work/hobi/jselwyn/Habitat/Intermediate_Files'
  SAVE_LOCATION<-'/work/hobi/jselwyn/Habitat/Results'
}

#### Get Data ####
full_data <- tibble(shoals = list.files(str_c(DATA_folder,'COPE_Sites/Shoals',sep='/'), recursive = TRUE, 
                                             pattern = 'shp$', full.names = TRUE)) %>%
  mutate(Site = str_extract(shoals, 'BZ17-[0-9]+[A-Z]+'),
         topography = str_c(INTERMEDIATE_FILES, '/Topography/', Site, '_', model_choice, '_stack.cm.tif')) %>%
  select(-shoals) %>%
  mutate(topography = map(topography, rast),
         topography = future_map(topography, as.matrix),
         topography = future_map(topography, ~as_tibble(.x))) %>%
  unnest(topography) 

overall_stats <- full_data %>%
  select(-starts_with('classified')) %>%
  summarise(across(-Site, .fns = list(mean = ~mean(.x, na.rm = TRUE), sd = ~sd(.x, na.rm = TRUE),
                                      n = ~sum(!is.na(.))), 
                   .names = "{.col}.{.fn}")) %>%
  pivot_longer(cols = everything(), 
               names_to = c('metric', '.value'), 
               names_pattern = '(.*)\\.(.*)') 

write_csv(select(overall_stats, -n), str_c(INTERMEDIATE_FILES, '/Topography/overall_metrics_', model_choice, '.csv'))
write_csv(overall_stats, str_c(SAVE_LOCATION, '/Topography/overall_metrics_', model_choice, '.csv'))



#### Make site metrics

site_stats <- full_data %>%
  group_by(Site) %>%
  summarise(across(-starts_with('classified'), .fns = list(mean = ~mean(.x, na.rm = TRUE), sd = ~sd(.x, na.rm = TRUE),
                                                           n = ~sum(!is.na(.))), 
                   .names = "{.col}.{.fn}"),
            across(starts_with('classified'), .fns = list(prop_reef = ~sum(. == 2, na.rm = TRUE)/sum(!is.nan(.)),
                                                          prop_sand = ~sum(. == 1, na.rm = TRUE)/sum(!is.nan(.))),
                   .names = "{.col}.{.fn}"),
            .groups = 'drop')

write_csv(select(site_stats, -ends_with('.n')), str_c(INTERMEDIATE_FILES, '/Topography/site_metrics_', model_choice, '.csv'))



global_moran <- tibble(shoals = list.files(str_c(DATA_folder,'COPE_Sites/Shoals',sep='/'), recursive = TRUE, 
                           pattern = 'shp$', full.names = TRUE)) %>%
  mutate(Site = str_extract(shoals, 'BZ17-[0-9]+[A-Z]+'),
         topography = str_c(INTERMEDIATE_FILES, '/Topography/', Site, '_', model_choice, '_stack.cm.tif')) %>%
  dplyr::select(-shoals) %>%
  mutate(topography = map(topography, rast),
         globalMoran = future_map_dbl(topography, ~Moran(raster(.x$depth)), .options = future_options(seed = TRUE))) %>%
  dplyr::select(-topography)

site_stats %>%
  left_join(global_moran, by = 'Site') %>%
  write_csv(str_c(SAVE_LOCATION, '/Topography/site_metrics_', model_choice, '.csv'))
