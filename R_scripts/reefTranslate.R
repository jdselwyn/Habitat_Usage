## Code to take all site raster stacks and translate them into a single grid of raster stacks


if(!interactive()){
  args <- commandArgs(trailingOnly=TRUE)
  model_choice <- args[1]
} else {
  model_choice <- 'c5'
}


#### Libraries ####
suppressWarnings(suppressMessages(library(raster)))
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(sf)))
suppressWarnings(suppressMessages(library(terra)))


#### Set up Computer ####
spacer <- if_else(Sys.info()['nodename'] == 'JASONDELL', 0, 1000/4) #divide by 4 is because you add in cells not cm's
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
  DATA_folder<-'/work/hobi/jselwyn/Habitat'
  INTERMEDIATE_FILES<-'/work/hobi/jselwyn/Habitat/Intermediate_Files'
  SAVE_LOCATION<-'/work/hobi/jselwyn/Habitat/Results'
}

#### Read in Data ####
DATA_folder<-path.expand(DATA_folder)
INTERMEDIATE_FILES<-path.expand(INTERMEDIATE_FILES)
SAVE_LOCATION<-path.expand(SAVE_LOCATION)

newproj<-"+proj=utm +zone=16Q +ellps=WGS84 +datum=WGS84 +units=m +no_defs" #PROBABLY NEEDS TO CHANGE

add_site <- function(brick, site){
  brick$site <- brick$depth
  values(brick$site) <- site
  brick$site <- mask(brick$site, brick$depth)
  
  brick
}

all_site_locations <- tibble(topography_stack = list.files(path = str_c(INTERMEDIATE_FILES,'Topography/',sep='/'), 
                                                           pattern = str_c(model_choice, 'stack.cm.tif$', sep = '_'), 
                                                           full.names = TRUE)) %>%
  mutate(Site = str_extract(topography_stack, 'BZ17-[0-9]*.[ABKNS]'),
         shoals = str_c(DATA_folder,'/COPE_Sites/Shoals/', Site, '.shp', sep='')) %>%
  dplyr::select(Site, everything()) %>%
  mutate(site_number = case_when(str_detect(Site, '0A') ~ 0,
                                 str_detect(Site, '0B') ~ 1,
                                 str_detect(Site, '100N') ~ 100,
                                 str_detect(Site, '100S') ~ -100,
                                 str_detect(Site, '10KN') ~ 10000,
                                 str_detect(Site, '1KN') ~ 1000,
                                 str_detect(Site, '1KS') ~ -1000,
                                 str_detect(Site, '500N') ~ 500,
                                 str_detect(Site, '500S') ~ -500,
                                 str_detect(Site, '5KN') ~ 5000,
                                 str_detect(Site, '5KS') ~ -5000,
                                 str_detect(Site, '60S') ~ -60)) %>%
  rowwise %>%
  mutate(topography_stack = list(rast(topography_stack)),
         # topography_stack = list(topography_stack$depth), #for testing purposes
         topography_stack = list(add_site(topography_stack, site_number)),
         shoals = list(st_read(shoals, quiet = TRUE)),
         shoals = list(st_transform(shoals, crs(topography_stack)))) %>%
  ungroup 


# plot(all_site_locations$topography_stack[[2]])
# plot(all_site_locations$shoals[[2]], col="red", pch=20, add = TRUE)

#### Center Stacks and Shoals ####

centering_step <- all_site_locations %>%
  # sample_frac(size = 1) %>% #random site shuffle
  rowwise %>%
  mutate(topography_stack = list(trim(topography_stack, padding = 0L)),
         topography_stack = list(expand(topography_stack, c(spacer, spacer))),
         min_coords = list(apply(bbox(topography_stack), 1, min)),
         site_dim = list(apply(bbox(topography_stack), 1, diff)),
         # center_translation = list(-1 * min_coords - site_dim/2),
         center_translation = list(-1 * min_coords),
         topography_stack = list(shift(topography_stack, dx = center_translation[1], dy = center_translation[2]))) %>%
  
  #Add in shoal translation to center 
  mutate(shoals = list(mutate(shoals, geometry = geometry + center_translation))) %>%
  
  ungroup

translation_instructions <- centering_step %>%
  dplyr::select(Site, min_coords, center_translation, site_dim)

# plot(centering_step$topography_stack[[2]])
# plot(centering_step$shoals[[2]], add = TRUE, pch = 20, col = 'red')


#### Fill in spacer gaps with nearest Neighbor ####
fill_in_buffer_factors <- function(stack){
  stack$site[is.na(stack$site)] <- as.numeric(unique(stack$site))
  
  
  prob_reef <- global(stack$classified_habitat_unsmoothed - 1, 'mean', na.rm = TRUE)$mean
  number_missing <- global(is.na(stack$classified_habitat_unsmoothed), 'sum')$sum
  values(stack$classified_habitat_unsmoothed)[is.na(values(stack$classified_habitat_unsmoothed))] <- 
    sample(c(1,2), size = number_missing, replace = TRUE, prob = c(1 - prob_reef, prob_reef))
  
  # prob_reef <- global(stack$classified_habitat_smoothed - 1, 'mean', na.rm = TRUE)$mean
  # number_missing <- global(is.na(stack$classified_habitat_smoothed), 'sum')$sum
  # values(stack$classified_habitat_smoothed)[is.na(values(stack$classified_habitat_smoothed))] <- 
  #   sample(c(1,2), size = number_missing, replace = TRUE, prob = c(1 - prob_reef, prob_reef))
 
  stack 
}

centering_step_2 <- centering_step %>%
  rowwise %>%
  mutate(topography_stack = list(fill_in_buffer_factors(topography_stack))) %>%
  ungroup


#### Translate Stacks and Shoals ####

smart_merger <- function(x){
  if(length(x) == 1){
    out <- x[[1]]
  } else {
    out <- do.call(merge, x)
  }
  out
}

#Create Rows of 3
horizontal_joining <- centering_step_2 %>%
  dplyr::select(-min_coords, -center_translation, -site_dim) %>%
  mutate(row_group = (row_number() - 1) %/% 3 + 1) %>%
  
  #Shift across to make rows of sites
  rowwise %>%
  mutate(original_bbox = list(apply(bbox(topography_stack), 1, diff))) %>%
  mutate(width = original_bbox[1]) %>%
  ungroup %>%
  
  mutate(width = max(width)) %>%
  
  group_by(row_group) %>%
  mutate(horizontal_shifter = cumsum(width),
         horizontal_shifter = lag(horizontal_shifter, 1, default = 0)) %>%
  ungroup %>%
  
  rowwise %>%
  mutate(topography_stack = list(shift(topography_stack , dx = horizontal_shifter, dy = 0)),
         ori = list(origin(topography_stack)),
         topography_stack = list(shift(topography_stack, dx = -1 * ori[1], dy = -1 * ori[2]))) %>%
  
  mutate(shoals = list(mutate(shoals, geometry = geometry + c(horizontal_shifter, 0))),
         shoals = list(mutate(shoals, geometry = geometry + -1 * ori))) %>%
  ungroup 

translation_instructions <- translation_instructions %>%
  full_join(horizontal_joining %>%
              dplyr::select(Site, row_group, original_bbox, width, horizontal_shifter, ori),
            by = 'Site')

horizontally_joined <- horizontal_joining %>%
  group_by(row_group) %>%
  summarise(reef_row = list(smart_merger(topography_stack)), 
            shoals = list(bind_rows(shoals)),
            .groups = 'drop')


# plot(horizontally_joined$reef_row[[2]])
# plot(horizontally_joined$shoals[[2]], add = TRUE, pch = 20, col = 'red')

#Join rows into single raster
reef_mosaic_precursor <- horizontally_joined %>%
  rowwise %>%
  mutate(row_bbox = list(apply(bbox(reef_row), 1, diff)),
         height = row_bbox[2]) %>%
  ungroup %>%
  mutate(height = max(height)) %>%
  mutate(vertical_shifter = cumsum(height),
         vertical_shifter = lag(vertical_shifter, 1, default = 0)) %>%
  
  rowwise %>%
  mutate(reef_row = list(shift(reef_row , dx = 0, dy = vertical_shifter)),
         ori = list(origin(reef_row)),
         reef_row = list(shift(reef_row, dx = -1 * ori[1], dy = -1 * ori[2]))) %>%
  
  mutate(shoals = list(mutate(shoals, geometry = geometry + c(0, vertical_shifter))),
         shoals = list(mutate(shoals, geometry = geometry + -1 * ori))) %>%
  ungroup 

translation_instructions <- translation_instructions %>%
  rename(horizontal_ori = ori) %>%
  full_join(reef_mosaic_precursor %>%
              dplyr::select(row_group, row_bbox, height, vertical_shifter, ori) %>%
              rename(vertical_ori = ori),
            by = 'row_group')

reef_mosaic <- reef_mosaic_precursor %>%
  summarise(reef_mosaic = list(smart_merger(reef_row)), 
            shoals = list(bind_rows(shoals)),
            .groups = 'drop') 

# plot(reef_mosaic$reef_mosaic[[1]])
# plot(reef_mosaic$shoals[[1]], add = TRUE, pch = 20, col = 'red')


#### Write Files ####
out_raster <- writeRaster(reef_mosaic$reef_mosaic[[1]], str_c(INTERMEDIATE_FILES, "/Topography/reef_mosaic_", model_choice, ".tif"), overwrite=TRUE)
st_write(reef_mosaic$shoals[[1]], str_c(INTERMEDIATE_FILES, "/Topography/shoal_mosaic_", model_choice,".shp"), 
         delete_dsn = file.exists(str_c(INTERMEDIATE_FILES, "/Topography/shoal_mosaic_", model_choice,".shp")), 
         quiet = TRUE)

pdf(str_c(SAVE_LOCATION, "/Topography/reef_mosaic_", model_choice,".pdf"), height = 10, width = 10, onefile = TRUE)
for(i in 1:length(names(reef_mosaic$reef_mosaic[[1]]))){
  plot(reef_mosaic$reef_mosaic[[1]][[i]], main = names(reef_mosaic$reef_mosaic[[1]])[i])
  reef_mosaic$shoals[[1]] %>%
    dplyr::select(contains("Shoal")) %>%
    vect %>%
    points()
}
dev.off()

write_rds(translation_instructions, str_c(INTERMEDIATE_FILES, "/Topography/translation_instructions_", model_choice, ".rds"))
