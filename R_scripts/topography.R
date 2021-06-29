## Create SpatialPixelsDataFrames for each habitat metric of interest to be used in INLA model ####

##TODO Evaluate viewshed - is it representing what I want? ex 5kn isn't showing greater degree of visibility at the top or consistently on the edges. that may be fine. Also inconsistent visibility of sand flats which should at least be able to see the rest of the sandy area.
##TODO Fix boundary distance to not count distance from the edge of raster!
##TODO Output pdf with all habitat metric plots
##TODO Boundary distance still not completely removing site boundary


args <- commandArgs(trailingOnly=TRUE)
model_choice <- args[1]
site_choice <- args[2]
# model_choice <- 'c5'
# site_choice <- 'BZ17-10KN'
print(paste('Start full topography', Sys.time(), sample(1000, 1), sep = ': '))

#### Libraries ####
library(sf)
library(raster)
library(terra)
# library(exactextractr)
library(tidyverse)
library(magrittr)
library(furrr)
library(rayshader)
# library(fractaldim)
library(whitebox)
library(spatialEco)
library(smoothr)

#Variable summary


#### Set up Computer ####
save_suffix <- 'All'
COMPUTER <- if_else(Sys.info()['nodename'] == 'JASONDELL', 'laptop', 'HPC')
PROG <- if_else(COMPUTER == 'HPC', FALSE, TRUE)

if(COMPUTER == 'Gawain'){
  cl <- makeClusterPSOCK(12)
  plan('cluster', workers = cl)
  
  DATA_folder<-'~/Documents/Coryphopterus/Maps' #Gawain
  INTERMEDIATE_FILES<-'~/Documents/Coryphopterus/Habitat Association (Paper Y - PhD Ch. Y)/Intermediate_Files'
  SAVE_LOCATION<-'~/Documents/Coryphopterus/Habitat Association (Paper Y - PhD Ch. Y)/Results'
  options(future.globals.maxSize = 1500*1024^2)
  
} else if(COMPUTER == 'laptop'){
  plan('multiprocess')
  # plan('sequential')
  
  DATA_folder<-'~/Coryphopterus/Maps' #Gawain
  INTERMEDIATE_FILES<-'~/Coryphopterus/Habitat Association (Paper Y - PhD Ch. Y)/Intermediate_Files'
  SAVE_LOCATION<-'~/Coryphopterus/Habitat Association (Paper Y - PhD Ch. Y)/Results'
  
  options(future.globals.maxSize = 1500*1024^2)
  
} else if (COMPUTER == 'HPC'){
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
  plan(multicore)
  
  DATA_folder<-'/work/hobi/jselwyn/Habitat'
  INTERMEDIATE_FILES<-'/work/hobi/jselwyn/Habitat/Intermediate_Files'
  SAVE_LOCATION<-'/work/hobi/jselwyn/Habitat/Results'
  
  options(future.globals.maxSize = 5000*1024^2)
}

#### Data Sources ####
DATA_folder<-path.expand(DATA_folder)
INTERMEDIATE_FILES<-path.expand(INTERMEDIATE_FILES)
SAVE_LOCATION<-path.expand(SAVE_LOCATION)

newproj<-"+proj=utm +zone=16Q +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

all_site_locations<-tibble(dem_file = list.files(path = str_c(DATA_folder,'COPE_Sites/DEM',sep='/'), pattern = 'tif$', full.names = TRUE),
                           # ortho_file = list.files(path = str_c(DATA_folder,'COPE_Sites/Orthomosaics',sep='/'), pattern = 'orthomosaic_demRez.tif$',
                           #                         full.names = TRUE),
                           # training_file = list.files(path = str_c(DATA_folder,'Classification/Training Data',sep='/'), pattern = 'shp$', full.names = TRUE),
                           
                           classified_habitat_smoothed = list.files(str_c(INTERMEDIATE_FILES, 'Habitat_Classification', sep = '/'), 
                                                                    recursive = TRUE, pattern = 'tif$', full.names = TRUE) %>%
                             str_subset(model_choice) %>%
                             str_subset('smoothed'),
                           classified_habitat_unsmoothed = list.files(str_c(INTERMEDIATE_FILES, 'Habitat_Classification', sep = '/'), 
                                                                      recursive = TRUE, pattern = 'tif$', full.names = TRUE) %>%
                             str_subset(model_choice) %>%
                             str_subset('smoothed', negate = TRUE),
                           
                           shoals_file = list.files(str_c(DATA_folder,'COPE_Sites/Shoals',sep='/'), recursive = TRUE, pattern = 'shp$', full.names = TRUE)) %>%
  mutate(Site = str_extract(shoals_file, 'BZ17-[0-9]*.[ABKNS]')) %>%
  dplyr::select(Site, everything(), -shoals_file) %>%
  inner_join(st_read(str_c(DATA_folder,'/COPE_Sites/Deepwater_Points.shp'), quiet = TRUE) %>%
               st_transform(crs=newproj) %>%
               dplyr::select(-id) %>%
               mutate(Site=as.character(Site)) %>%
               rename(deepwater = geometry), 
             by='Site')

#### Functions ####
Scale <- 1/100
Sun_Angle <- c(5) #This need to change to not have 3 specifically coded in but rather change the number dynamically
View_Height <- 1
View_Resolution <- 1e0
raster_window <- matrix(1, nrow = 11, ncol = 11) #queens case 3x3 matrix
#gaussian.kernel(sigma = 5, n = 15)
min_sand_area <- 0.25 #the minimum contiguous sand area for boundary distance calculations. If no sand area larger than this will default to max sand area size

## Functions
make_negative<-function(x){
  if(global(x, 'max', na.rm = TRUE) > 0){x <- -1 * x}
  x
}

fix_depth<-function(x, site, denoise = FALSE, ...){
  message(str_c('Start reprojecting Depth', site, Sys.time(), sep = ': '))
  the_file<-str_c(INTERMEDIATE_FILES,'/Topography/',site,if_else(denoise, '_denoised', ''),'_depth.cm.tif',sep='')
  
  if(file.exists(the_file)){  #If csv does exist then read it in.
    
    out <- rast(the_file)
    #the_file
    
  } else {#If it doesn't exist then process site and write it
    if(!denoise){
      out <- rast(x) %>%
        set_names('depth') %>%
        make_negative %>%
        # trim %>% 
        raster %>%
        projectRaster(crs = newproj) %>%
        rast %T>%
        writeRaster(filename = the_file, overwrite=T)
    } else {
      #NOT FUNCTIONAL
      tmp_out <- tempfile(pattern = 'out', fileext = '.tif') 
      
      wbt_feature_preserving_smoothing(dem = x, output = tmp_out, filter = 3)
      
      out <- rast(tmp_out) %>%
        set_names('depth') %>%
        make_negative %>%
        trim %>%
        projectRaster(crs=newproj) %T>%
        writeRaster(filename=the_file, overwrite=T)
    }
  }
  message(str_c('Finished reprojecting Depth', site, Sys.time(), sep = ': '))
  out
}

get_basic_terrain<-function(x, site, version, ...){
  message(str_c('Start calculating ', version, ' ', site, ': ', Sys.time()))
  the_file<-str_c(INTERMEDIATE_FILES,'/Topography/',site,'_',version,'.cm.tif',sep='')
  
  #Check if csv of site already exists
  if(file.exists(the_file)){  #If csv does exist then read it in.
    
    #raster(the_file)
    # the_file
    1+1
    
  } else {#If it doesn't exist then process site and write it
    
    if(version=='moran'){
      tmp <- MoranLocal(raster(x)) %T>%
        writeRaster(filename = the_file,overwrite=TRUE)
    } else {
      tmp <- terrain(raster(x), version, ..., filename = the_file, overwrite=TRUE)
    }
    
    # the_file
    
  }
  message(str_c('Finished calculating ', version, ' ', site, ': ', Sys.time()))
  the_file
}

get_Directionality <- function(aspect, site, version){
  message(str_c('Start calculating ', version, ' ', site, ': ', Sys.time()))
  the_file<-str_c(INTERMEDIATE_FILES,'/Topography/',site,'_',version,'.cm.tif',sep='')
  
  #Check if csv of site already exists
  if(file.exists(the_file)){  #If csv does exist then read it in.
    
    #raster(the_file)
    # the_file
    1+1
    
  } else {#If it doesn't exist then process site and write it
    
    tmp <- (rast(aspect) * pi / 180)
    
    if(version=='northness'){
      tmp2 <- cos(tmp) %T>%
        writeRaster(filename = the_file, overwrite=TRUE)
      
    } else {
      tmp2 <- sin(tmp) %T>%
        writeRaster(filename = the_file, overwrite=TRUE)
    }
    
    # the_file
    
  }
  message(str_c('Finished calculating ', version, ' ', site, ': ', Sys.time()))
  the_file
  
}

get_viewshed <- function(site, H, res){
  message(str_c('Start calculating ', 'Viewshed', ' ', site, ': ', Sys.time()))
  the_file<-str_c(INTERMEDIATE_FILES,'/Topography/',site,'_viewshed_',H,'.cm.tif',sep='')
  
  if(file.exists(the_file)){
    1+1
  } else {
    setwd(INTERMEDIATE_FILES)
    tmp_out <- tempfile(fileext = '.tif')
    
    wbt_visibility_index(dem = str_c('./Topography/',site, '_depth.cm.tif', sep=''), 
                         output = tmp_out, 
                         height = H, res_factor = res, verbose_mode = FALSE)
    
    views <- raster(tmp_out)
    depth <- raster(str_c('./Topography/',site, '_depth.cm.tif', sep=''))
    mask(views, depth) %>%
      writeRaster(the_file, overwrite=TRUE)
  }
  message(str_c('Finished calculating ', 'Viewshed', ' ', site, ': ', Sys.time()))
  the_file
}

calc_sandDist <- function(r, site, min_area, smoothed = FALSE){
  message(str_c('Start calculating ', 'Sand Distance', ' ', site, ': ', Sys.time()))
  the_file<-str_c(INTERMEDIATE_FILES,'/Topography/',site,'_','sandDistance_', min_area, '_', 
                  model_choice,'.', if_else(smoothed, 'smoothed.', ''),
                  'cm.tif',sep='')
  
  if(file.exists(the_file)){  #If csv does exist then read it in.
    
    1+1
    
  } else {#If it doesn't exist then process site and write it
    
    # plot(r)
    ## Create external border
    border <- r %>%
      is.na %>% 
      equals(0) %>% 
      classify(matrix(c(0, NA), ncol = 2)) %>% 
      as.polygons()
    
    
    tmpf4 <- tempfile()
    writeVector(border, tmpf4, overwrite=TRUE)
    border2 <- st_read(tmpf4, crs = crs(border), quiet = TRUE) %>%
      select(geometry) 
    
    if(any(str_detect(class(border2$geometry), 'MULTIPOLYGON'))){
      border2 <- border2 %>%
        st_cast("MULTILINESTRING")
    }
    
    border2 <- border2 %>%
      st_cast('LINESTRING') %>%
      st_cast('POLYGON') %>%
      mutate(area = st_area(geometry)) %>%
      filter(area == max(area)) %>%
      st_cast('LINESTRING') %>%
      select(geometry)
    
    
    ## Create internal borders between sand/reef
    r2 <- as.polygons(r)
    tmpf2 <- tempfile()
    writeVector(r2, tmpf2, overwrite=TRUE)
    
    boundary <- st_read(tmpf2, crs = crs(r), quiet = TRUE) %>%
      rename(habitat = 1) %>%
      mutate(habitat = if_else(habitat == 1, 'Sand', 'Reef')) %>%
      st_cast('POLYGON') %>%
      mutate(area = map(geometry, st_area),
             area = as.numeric(area)) %>% 
      filter(habitat == 'Sand')
    
    min_area <- if_else(min_area >= max(boundary$area, na.rm = TRUE), max(boundary$area, na.rm = TRUE), min_area)
    
    boundary2 <- boundary %>%
      filter(area >= min_area) %>% 
      fill_holes(threshold = min_area) %>%
      select(-habitat) 
    
    boundary3 <- boundary2 %>%
      # select(-area) %>%
      st_cast('LINESTRING') %>%
      vect 
    
    
    tmpf3 <- tempfile()
    writeVector(boundary3, tmpf3, overwrite=TRUE)
    boundary4 <- st_read(tmpf3, crs = crs(boundary3), quiet = TRUE) %>%
      select(-area)
    
    tmpf5 <- str_replace(the_file, '\\.tif$', '.shp')
    tmp <- st_difference(boundary4, border2) %>%
      st_write(tmpf5, quiet = TRUE, delete_dsn = file.exists(tmpf5))
    
    
    tmp_bound <- vect(tmpf5)
    
    out <- distance(r, tmp_bound)
    out2 <- terra::mask(out, r)
    
    writeRaster(out2, the_file, overwrite=TRUE)
    
    
    #out
  }
  message(str_c('Finished calculating ', 'Sand Distance', ' ', site, ': ', Sys.time()))
  the_file
}




calc_rugosity <- function(r, site, resolution){
  message(str_c('Start calculating ', 'rugosity', ' ', site, ': ', Sys.time()))
  the_file<-str_c(INTERMEDIATE_FILES,'/Topography/',site,'_','rugosity','.',
                  nrow(resolution),'x', ncol(resolution),'.cm.tif',sep='')
  
  if(file.exists(the_file)){  #If csv does exist then read it in.
    
    1+1
    
  } else {#If it doesn't exist then process site and write it
    
    r_sp <- as(raster(r), 'SpatialPixelsDataFrame')
    r_sa <- raster(surfaceArea(r_sp, byCell = TRUE))
    r_rug <- focal(r_sa, resolution) /(prod(res(r)) * length(resolution))
    
    writeRaster(r_rug, the_file)
    
    
    #out
  }
  
  message(str_c('Finished calculating ', 'rugosity', ' ', site, ': ', Sys.time()))
  the_file
}

build_ShadowRealm <- function(dem, site, distant_point, sun_angle, scale, internal_parallel){
  message(str_c('Finished calculating ', 'Shadow Realm', ' ', site, ': ', Sys.time()))
  the_file<-str_c(INTERMEDIATE_FILES,'/Topography/',site,'_','rayShade','.',
                  sun_angle,'.cm.tif',sep='')
  
  if(file.exists(the_file)){  #If csv does exist then read it in.
    
    1+1
    
  } else {#If it doesn't exist then process site and write it
    distant_point <- distant_point %>%
      st_sfc(crs = newproj) %>%
      st_sf %>%
      st_transform(crs = as.character(crs(dem))) %>% 
      st_coordinates()
    
    sun_heading <- (c(xmax(dem) + xmin(dem), ymax(dem) + ymin(dem))/2) %>%
      st_point %>%
      st_sfc(crs = as.character(crs(dem))) %>%
      st_sf %>%
      st_coordinates() %>%
      as_tibble %>%
      mutate(bearing = atan2(distant_point[1] - X, distant_point[2] - Y) * 180 / pi) %>%
      pull(bearing)
    
    out <- raster(dem) %>%
      as.matrix() %>%
      ray_shade(heightmap = .,
                lambert = FALSE, 
                anglebreaks = seq(0, sun_angle, length.out = 25), #Determine best/most meaningful setting for this (angle of sun over the horizon)
                sunangle = sun_heading, 
                zscale = scale,
                multicore = internal_parallel) %>% #Set up to work on multiple cores at once
      # .[,ncol(.):1] %>% 
      # .[nrow(.):1,ncol(.):1] %>%
      .[nrow(.):1,] %>%
      rast(crs = crs(dem)) %>%
      raster %>%
      setExtent(ext = extent(raster(dem))) %>% 
      mask(mask = raster(dem))
    
    writeRaster(out, the_file, overwrite = TRUE)
  }
  
  message(str_c('Finished calculating ', 'Shadow Realm', ' ', site, ': ', Sys.time()))
  the_file
}

calc_normal <- function(x, site, scale){
  
  remove_padding <- function(x){
    x[2:(nrow(x)-1), 2:(ncol(x)-1)]
  }
  
  message(str_c('Start calculating ', 'Normal Vectors', ' ', site, ': ', Sys.time()))
  the_file<-str_c(INTERMEDIATE_FILES,'/Topography/',site,'_normalVectors.tif',sep='')
  
  if(file.exists(the_file)){  #If csv does exist then read it in.
    
    out <- stack(the_file)
    #the_file
    
  } else {#If it doesn't exist then process site and write it
    x <- raster(x)
    
    out <- as.matrix(x) %>%
      calculate_normal(zscale = scale) %>% 
      map(t) %>%
      map(remove_padding) %>%
      map(raster, crs=crs(x)) %>% 
      stack %>%
      setExtent(ext=extent(x)) %>%
      mask(x) %T>%
      writeRaster(filename=the_file, overwrite=T)
  }
  
  message(str_c('Finished calculating ', 'Normal Vectors', ' ', site, ': ', Sys.time()))
  out
}

split_normal <- function(x, site, direction){
  
  the_file<-str_c(INTERMEDIATE_FILES,'/Topography/',site,'_normalVector',str_to_upper(direction),'.tif',sep='')
  
  if(file.exists(the_file)){  #If csv does exist then read it in.
    
    #stack(the_file)
    the_file
    
  } else {#If it doesn't exist then process site and write it
    tmp <- x %>%
      subset(which(str_detect(c('x','y','z'), direction))) %>%
      writeRaster(the_file)
    
    the_file
    
  }
}

calc_VectorDispersion <- function(x, y, z, site, resolution){
  message(str_c('Start calculating ', 'Vector Dispersion', ' ', site, ': ', Sys.time()))
  the_file<-str_c(INTERMEDIATE_FILES,'/Topography/',site,'_','vectorDispersion','.',
                  nrow(resolution),'x', ncol(resolution),'.cm.tif',sep='')
  
  if(file.exists(the_file)){  #If csv does exist then read it in.
    
    1+1
    
  } else {#If it doesn't exist then process site and write it
    
    focal_stack <- list(x, y, z) %>%
      map(rast) %>%
      map(~focal(.x, w = resolution, fun = 'sum')) %>%
      map(~.x^2)
    
    counter <- focal((!is.na(rast(x))), w = raster_window, fun = sum)
    
    vectDisp <- (counter - sqrt(focal_stack[[1]] + focal_stack[[2]] + focal_stack[[3]])) / (counter - 1)
    
    writeRaster(vectDisp, the_file)
    #out
  }
  message(str_c('Finished calculating ', 'Vector Dispersion', ' ', site, ': ', Sys.time()))
  the_file
}

calc_Curvature <- function(x, site, version, resolution){
  message(str_c('Start calculating ', version, ' curvature ', site, ': ', Sys.time()))
  
  the_file<-str_c(INTERMEDIATE_FILES,'/Topography/',site,'_', version,'_','curvature','.',
                  nrow(resolution),'x', ncol(resolution),'.cm.tif',sep='')
  
  if(file.exists(the_file)){  #If csv does exist then read it in.
    
    1+1
    
  } else {#If it doesn't exist then process site and write it
    
    out <- raster(x) %>%
      curvature(type = version) %>%
      rast %>%
      focal(w = resolution/length(resolution), fun = 'sum')
    
    writeRaster(out, the_file)
    
    #out
  }
  message(str_c('Finished calculating ', version, ' curvature ', site, ': ', Sys.time()))
  the_file
}


#### Process DEM ####
habitat_metrics <- all_site_locations %>%
  filter(Site == site_choice) %>%
  # slice(1) %>%
  
  mutate(depth = map2(dem_file, Site, fix_depth),
         across(starts_with('classified'), ~map(., ~rast(.x)))) %>%
  select(-dem_file) %>%
  mutate(aspect = map2_chr(depth, Site, get_basic_terrain, version='aspect', unit = 'degrees'),
         slope = map2_chr(depth, Site, get_basic_terrain, version='slope'),
         moran = map2_chr(depth, Site, get_basic_terrain, version='moran'),
         relief = map2_chr(depth, Site, get_basic_terrain, version='roughness')) %>%
  mutate(northness = map2_chr(aspect, Site, get_Directionality, version = 'northness'),
         eastness = map2_chr(aspect, Site, get_Directionality, version = 'eastness')) %>%
  
  mutate(viewshed = map_chr(Site, get_viewshed, H = View_Height, res = View_Resolution)) %>%
  mutate(boundaryDistance = map2_chr(classified_habitat_unsmoothed, Site, calc_sandDist, min_area = min_sand_area)) %>%
  mutate(boundaryDistanceSmoothed = map2_chr(classified_habitat_smoothed, Site, calc_sandDist,
                                             min_area = min_sand_area, smoothed = TRUE)) %>%
  mutate(rugosity = map2_chr(depth, Site, calc_rugosity, resolution = raster_window)) %>%
  mutate(normalVectors = future_map2(depth, Site, calc_normal, 
                                     scale = Scale), #Check scale - should be right
         normVectX = future_map2_chr(normalVectors, Site, split_normal, direction='x'),
         normVectY = future_map2_chr(normalVectors, Site, split_normal, direction='y'),
         normVectZ = future_map2_chr(normalVectors, Site, split_normal, direction='z'),
         vectorDispersion = pmap_chr(list(normVectX, normVectY, normVectZ, Site), 
                                     calc_VectorDispersion, resolution = raster_window)) %>%
  mutate(rayShade = pmap_chr(list(depth, Site, deepwater), build_ShadowRealm,
                             sun_angle = Sun_Angle, scale = Scale, internal_parallel = TRUE)) %>%
  mutate(planCurvature = map2_chr(depth, Site, calc_Curvature, 
                                  version = 'planform', resolution = raster_window),
         profCurvature = map2_chr(depth, Site, calc_Curvature, 
                                  version = 'profile', resolution = raster_window),
         totalCurvature = map2_chr(depth, Site, calc_Curvature, 
                                   version = 'total', resolution = raster_window)) %>%
  
  #Put new metrics here!
  
  select(-starts_with('norm'), -deepwater) %>%
  mutate(across(c(where(is.character), -Site), ~map(., rast))) %>%
  pivot_longer(cols = -Site,
               names_to = 'metric',
               values_to = 'raster') %>%
  mutate(raster = map2(raster, metric, ~set_names(.x, .y))) %>%
  group_by(Site) %>%
  summarise(raster_stack = list(do.call(c, raster)), .groups = 'drop') %>%
  mutate(file_name = str_c(INTERMEDIATE_FILES,'/Topography/', Site, '_', model_choice, 
                           '_','stack','.cm.tif',sep='')) %>%
  mutate(file_out = map2(raster_stack, file_name, ~writeRaster(.x, .y, overwrite = TRUE)))

#Try spatialEco metrics
#### Make plots ####

topography <- rast(str_c(INTERMEDIATE_FILES,'/Topography/', site_choice, '_', model_choice, 
                         '_','stack','.cm.tif',sep=''))


message(str_c('Start making data plots: ', Sys.time()))

pdf(str_c(SAVE_LOCATION,'/Topography/', site_choice, '_', model_choice, '_topography.pdf',sep=''),
    height = 10, width = 10, onefile = TRUE)

for(i in 1:length(names(topography))){
  
  plot(topography[[i]], main = names(topography)[i])
  
}
dev.off()
message(str_c('Finish making expanded data plots: ', Sys.time()))
