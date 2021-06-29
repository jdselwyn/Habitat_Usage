#### Run after running "classify_habitat.R" and choose prefered model to extract final top model and fit to site maps ####

args <- commandArgs(trailingOnly=TRUE)
model_choice <- args[1]
site_choice <- args[2]
# model_choice <- 'c5'
# site_choice <- 'BZ17-0A'

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

#### Functions ####
zero_to_NA<-function(x){
  x$Red[x$Red == 0] <- NA
  x$Green[x$Green == 0] <- NA
  x$Blue[x$Blue == 0] <- NA
  x
}

ortho_to_dem_res <- function(x, y){
  out_name <- str_replace(x, '\\.tif', '_demRez.tif')
  
  if(file.exists(out_name)){
    out_name
  } else {
    dem <- raster(y) %>%
      projectRaster(crs = newproj) %>% 
      rast
    
    ortho <- rast(x) %>%
      terra::subset(1:3) %>%
      set_names(c('Red', 'Green', 'Blue')) %>%
      zero_to_NA %>%
      terra::project(y = dem)
    
    for(color in c('Red', 'Green', 'Blue')){
      start_color <- ortho[[color]]
      start_color <- mask(start_color, dem)
      
      z <- 0
      gap <- 10
      while(gap != 0){
        #Fill in missing values with nearest neighbor mean
        if(z == 0){
          tmp <- start_color
          na_check <- is.na(tmp)
          previous_count_na <- global(na_check, "sum")
        } else {
          previous_count_na <- current_count_na
        }
        
        tmp <- focal(tmp, w = 9, na.only = TRUE, fun = 'mean')
        tmp <- mask(tmp, dem)
        na_check <- is.na(tmp)
        current_count_na <- global(na_check, "sum")
        gap <- previous_count_na - current_count_na
        z <- z + 1
      }
      ortho[[color]] <- tmp
    }
    
    out <- writeRaster(ortho, out_name, overwrite = TRUE)
    
    out_name
  }
  
}


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

all_site_locations <- tibble(dem_file = list.files(path = str_c(DATA_folder,'COPE_Sites/DEM',sep='/'), pattern = 'tif$', full.names = TRUE),
                             ortho_file = list.files(path = str_c(DATA_folder,'COPE_Sites/Orthomosaics',sep='/'), pattern = 'orthomosaic.tif$', 
                                                     full.names = TRUE),
                             training_file = list.files(path = str_c(DATA_folder, classified_habitat, sep='/'), pattern = 'shp$', full.names = TRUE),
                             shoals_file = list.files(str_c(DATA_folder,'COPE_Sites/Shoals',sep='/'), recursive = TRUE, pattern = 'shp$', full.names = TRUE)) %>%
  mutate(Site = str_extract(shoals_file, 'BZ17-[0-9]*.[ABKNS]')) %>%
  dplyr::select(Site, everything()) %>%
  
  filter(Site == site_choice) %>%
   
  mutate(ortho_file = map2_chr(ortho_file, dem_file, ortho_to_dem_res)) #Change resolution of orthomosaic to match the DEM

top_model <- read_rds(str_c(SAVE_LOCATION, '/Habitat_Classification/', save_suffix, "_", model_choice, '_completeModel.rds', sep = ''))

#### Predict onto sites ####
site_predictions <- all_site_locations %>%
  # slice(1) %>% #for testing
  select(Site, ortho_file) %>%
  mutate(rast = map(ortho_file, ~stack(.x) %>%
                      dropLayer(4) %>%
                      set_names(c('Red', 'Green', 'Blue'))),
         rows = map_int(rast, nrow),
         cols = map_int(rast, ncol)) %>%
  mutate(xBuf = cols/1, #Amount of resolution reduction - save memory but loose res. Bigger constant is lower res
         yBuf = rows/1) %>%
  mutate(rastIO = pmap(list(cols, rows, xBuf, yBuf), ~list(nXOff = 1, nYOff = 1, 
                                                           nXSize = ..1, nYSize = ..2, 
                                                           nBufXSize = ..3, nBufYSize = ..4, 
                                                           bands = c(1, 2, 3)))) %>%
  
  mutate(the_star = pmap(list(ortho_file, rastIO, Site), ~read_stars(..1, RasterIO = ..2) %>%
                           split('band') %>%
                           mutate(Red = X1, 
                                  Green = X2, 
                                  Blue = X3) %>%
                           select(-X1:-X3) %>%
                           mutate(Site = ..3))) %>%
  mutate(the_predictions = map(the_star, ~as_tibble(.x) %>%
                                 mutate(Red = if_else(is.na(Red), 0, Red),
                                        Green = if_else(is.na(Green), 0, Green),
                                        Blue = if_else(is.na(Blue), 0, Blue)) %>%
                                 mutate(class = NA) %>%
                                 predict(top_model, new_data = .))) %>%
  mutate(the_star = map2(the_star, the_predictions, ~mutate(.x, class = .y$.pred_class,
                                                            class = ifelse(Red == 0 | Green == 0 | Blue == 0, NA, class)) %>%
                           select(class))) %T>%
  mutate(file_out = map2(Site, the_star, ~write_stars(.y, str_c(INTERMEDIATE_FILES, '/Habitat_Classification/', .x, "_", 
                                                                save_suffix, '_', model_choice, '_New_habitat_class.tif',sep=''))))



#### Make smoothed classification version ####
input <- str_c(INTERMEDIATE_FILES, '/Habitat_Classification/', site_choice, "_", save_suffix, '_', model_choice, '_New_habitat_class.tif',sep='')
output <- str_c(INTERMEDIATE_FILES, '/Habitat_Classification/', site_choice, "_", save_suffix, '_', model_choice, '_New_habitat_class_smoothed.tif',sep='')

wbt_majority_filter(input = input,
                    output = output,
                    filterx = 101, filtery = 101,
                    # filterx = 11, filtery = 11, 
                    verbose_mode = FALSE)

out_rast <- mask(x = raster(output), 
                 mask = raster(input),
                 filename = output,
                 overwrite = TRUE)


parallel::stopCluster(cl)
