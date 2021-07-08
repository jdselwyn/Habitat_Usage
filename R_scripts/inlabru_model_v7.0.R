#### Code to fit and evaluate inlabru model ####

##TODO include mesh settings in sourced parameter file
##TODO make it so if there's already a param file it takes that one (right now can just specify the one that's already there)

if(interactive()){
  param_file <- '/home/jselwyn/Habitat/params/inlabru_params.R'
  run_name <- '7.7.21_test.8'
} else {
  args <- commandArgs(trailingOnly=TRUE)
  param_file <- args[1]
  run_name <- args[2]
}

message(run_name)
message(param_file)

newproj<-"+proj=utm +zone=16Q +ellps=WGS84 +datum=WGS84 +units=m +no_defs"


#### Set up Computer ####
COMPUTER <- ifelse(Sys.info()['nodename'] == 'JASONDELL', 'laptop', 
                    ifelse(Sys.info()['nodename'] == 'jdselwyn', 'Gawain', 'HPC'))

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

#### Set up Directory Structure ####
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
dir.create(path = str_c(INTERMEDIATE_FILES, '/INLA_models/', 
                        str_c('models_', run_name, '/', sep = ''), 
                        sep = ''), 
           showWarnings = FALSE)

dir.create(path = str_c(SAVE_LOCATION, '/INLA_models/', 
                        str_c('models_', run_name, '/', sep = ''), 
                        sep = ''), 
           showWarnings = FALSE)

#### Raw copy files from topography to inlabru to keep everything together ####
#shoals
#site polygon
#topography
#habitat
#sites
#translation instructions

source(param_file)

file_copying <- c(
  param_file, 
  
  str_c(INTERMEDIATE_FILES, '/INLA_models/', 
        c('site', 'topography', 'habitat_type', 'site_selection'),
        '.rds', sep = ''),
  
  str_c(INTERMEDIATE_FILES, '/Topography/', 
        c(
          str_c('translation_instructions_', model_choice, '.rds'),
          str_c('shoal_mosaic_', model_choice, c('.shp', '.dbf', '.shx'))
        ), 
        sep = '')
) %>%
  file.copy(to = str_c(INTERMEDIATE_FILES, '/INLA_models/models_', run_name, sep = ''))


#### Libraries ####
suppressMessages(library(inlabru))
suppressMessages(library(INLA))
suppressMessages(library(sf))
suppressMessages(library(future))
suppressMessages(library(furrr))
suppressMessages(library(doFuture))
suppressMessages(library(raster))
suppressMessages(library(terra))
suppressMessages(library(doRNG))
suppressMessages(library(progressr))
suppressMessages(library(maptools))
suppressMessages(library(spatstat))

select <- dplyr::select


registerDoFuture()
plan('multicore', workers = number_workers)

#### Read in data ####
message(str_c('Start Reading in data: ', Sys.time()))
site_poly <- read_rds(str_c(INTERMEDIATE_FILES, '/INLA_models/models_', run_name, 
                            '/site.rds',
                            sep = ''))
shoals <- st_read(str_c(INTERMEDIATE_FILES, '/INLA_models/models_', run_name, 
                        '/shoal_mosaic_', model_choice, '.shp',
                        sep = ''),
                  quiet = TRUE)


topography <- read_rds(str_c(INTERMEDIATE_FILES, '/INLA_models/models_', run_name, 
                             '/topography.rds',
                             sep = ''))
habitat_type <- read_rds(str_c(INTERMEDIATE_FILES, '/INLA_models/models_', run_name, 
                               '/habitat_type.rds',
                               sep = ''))
site_track <- read_rds(str_c(INTERMEDIATE_FILES, '/INLA_models/models_', run_name, 
                             '/site_selection.rds',
                             sep = ''))

translation_instructions <- read_rds(str_c(INTERMEDIATE_FILES, '/INLA_models/models_', 
                                           run_name, 
                                           '/translation_instructions_', model_choice, '.rds',
                                           sep = '')) %>%
  rowwise %>%
  mutate(movement = list(-1 * ((-1 * vertical_ori) + c(0, vertical_shifter) + (-1 * horizontal_ori) + c(horizontal_shifter, 0) + center_translation))) %>%
  ungroup %>%
  select(Site, movement)

message(str_c('Finished Reading in data: ', Sys.time()))

#### Distribute Fish around shoals ####
message(str_c('Start Distributing Fish: ', Sys.time()))
distribute_fish <- function(number_fish, centroid, spread_m = 0.25){
  #Distribute fish around centroid within a shoal. Assumes crs is in "m" and will distribute with stdev of spread_m (in m)
  if(number_fish == 1){
    out <- as.numeric(centroid) %>%
      as_tibble() %>%
      mutate(name = c('x', 'y')) %>%
      pivot_wider(names_from = 'name', values_from = 'value')
  } else {
    out <- MASS::mvrnorm(n = number_fish, mu = as.numeric(centroid), 
                         Sigma = diag(spread_m, nrow = 2, ncol = 2)) %>%
      set_colnames(c('x', 'y')) %>%
      as_tibble()
  }
  out
}

fish_raw <- shoals %>%
  select(Shoal_Size) %>%
  rename(fish = Shoal_Size) %>%
  filter(fish > 0) %>%
  st_join(st_as_sf(site_poly), st_within) %>%
  select(-weight) %>%
  rename(site = names) %>%
  mutate(points = map2(fish, geometry, distribute_fish, spread_m = 0.06)) %>%
  as_tibble %>%
  dplyr::select(-fish) %>%
  dplyr::select(-geometry) %>%
  unnest(points) %>%
  st_as_sf(coords = c('x', 'y'),
           crs = NA)

fish_sf <- fish_raw %>%
  left_join(translation_instructions, by = c('site' = 'Site')) %>%
  rowwise %>%
  mutate(geometry = geometry + movement) %>%
  ungroup %>%
  select(-movement) %>%
  st_sf(crs = newproj) %T>%
  st_write(str_c(SAVE_LOCATION, '/INLA_models/models_',
                 run_name, '/observed_fish.shp',
                 sep = ''),
           quiet = TRUE,
           delete_dsn = file.exists(str_c(SAVE_LOCATION, '/INLA_models/models_',
                                          run_name, '/observed_fish.shp',
                                          sep = '')))

fish <- fish_raw %>%
  select(geometry) %>%
  as_Spatial()
colnames(fish@coords) <- c('x', 'y')

message(str_c('Finished Distributing Fish: ', Sys.time()))

#### Set up mesh ####
message(str_c('Start Building Mesh: ', Sys.time()))
boundary.loc <- rbind(SpatialPoints(INLA::inla.sp2segment(site_poly)$loc,
                                    proj4string = CRS(SRS_string = wkt(site_poly))),
                      as(fish, "SpatialPoints"))

# boundary <- list(inla.nonconvex.hull(coordinates(boundary.loc), 2.5))
# 
# mesh <- inla.mesh.2d(#boundary = inla.sp2segment(sites),
#   boundary = boundary,
#   max.edge=c(0.5),
#   cutoff=0.25,
#   offset = c(2.5),
#   crs = CRS(NA_character_))
# 
# matern <- inla.spde2.pcmatern(mesh,
#                               prior.range = c(5, 0.5),
#                               prior.sigma = c(2, 0.01))

boundary <- list(inla.nonconvex.hull(coordinates(boundary.loc), mesh_setting$offset[1]),
                 inla.nonconvex.hull(coordinates(boundary.loc), mesh_setting$offset[2]))

mesh <- inla.mesh.2d(boundary=boundary,
                     max.edge = mesh_setting$max.edge, 
                     min.angle = mesh_setting$min.angle,
                     max.n = mesh_setting$max.n,
                     max.n.strict = mesh_setting$max.n.strict, 
                     cutoff = mesh_setting$cutoff,
                     offset = mesh_setting$offset, ## Offset for extra boundaries, if needed.
                     crs = CRS(mesh_setting$crs)) 

matern <- inla.spde2.pcmatern(mesh, 
                              prior.range = matern_priors$prior.range,
                              prior.sigma = matern_priors$prior.sigma)

mesh_plot <- ggplot() + 
  gg(mesh) +
  gg(site_poly) +
  gg(fish) + 
  coord_fixed() +
  theme_void()
ggsave(str_c(SAVE_LOCATION, '/INLA_models/models_',
             run_name, '/mesh_', run_name, '.png',
             sep = ''), 
       plot = mesh_plot, height = 30, width = 30)

for(i in 1:nrow(site_poly)){
  site_box <- bbox(site_poly[i,])
  
  mesh_plot_site <- mesh_plot + 
    coord_fixed(xlim = site_box['x',], ylim = site_box['y',]) +
    theme_void()
  
  ggsave(str_c(SAVE_LOCATION, '/INLA_models/models_',
               run_name, '/mesh_', site_poly$names[i], '_', 
               run_name, '.png',
               sep = ''), 
         plot = mesh_plot, height = 30, width = 30)
}

message(str_c('Finished Building Mesh: ', Sys.time()))

#### Fit Specified Model ####
message(str_c('Start Fitting Model: ', Sys.time()))
continuous_covars <- function(x, y, metric){
  # turn coordinates into SpatialPoints object:
  spp <- SpatialPoints(data.frame(x = x, y = y)) 
  # attach the appropriate coordinate reference system (CRS)
  proj4string(spp) <- CRS(proj4string(topography))
  # Extract metric values at spp coords, from SpatialGridDataFrame
  v <- over(spp, topography) 
  out <- v[[metric]]
  out[is.na(out)] <- median(out, na.rm = TRUE) # NAs are a problem! Remove them
  
  return(out)
}

model_fit <- lgcp(model_formula, 
                  fish, 
                  samplers = site_poly,
                  E = 1,
                  domain = list(coordinates = mesh),
                  options = list(verbose = TRUE, 
                                 quantiles = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975),
                                 control.compute = list(config = TRUE, 
                                                        dic = TRUE, 
                                                        waic = TRUE,
                                                        cpo = TRUE),
                                 control.inla = list(strategy = "simplified.laplace",
                                                     int.strategy = integration_strategy)))

write_rds(model_fit,
          str_c(INTERMEDIATE_FILES, 
                '/INLA_models/models_', run_name, 
                '/model_fit_', run_name, '_', model_choice, '.rds',
                sep = ''),
          compress = 'xz')

# model_fit <- read_rds(str_c(INTERMEDIATE_FILES,
#                '/INLA_models/models_', run_name,
#                '/model_fit_', run_name, '_', model_choice, '.rds',
#                sep = ''))

message(str_c('Finished Fitting Model: ', Sys.time()))

#### Output parameter estimates etc. ####
message(str_c('Start Parameter Extraction: ', Sys.time()))
parameter_estimates <- bind_rows(
  
  fixed = model_fit$marginals.fixed %>%
    map(as_tibble) %>%
    bind_rows(.id = 'parameter'),
  
  hyper_par = model_fit$marginals.hyperpar %>%
    map(as_tibble) %>%
    bind_rows(.id = 'parameter'),
  
  random = model_fit$marginals.random %>%
    map(~map(.x, as_tibble) %>%
          bind_rows(.id = 'index')) %>%
    bind_rows(.id = 'parameter'),
  
  .id = 'type'
) %>%
  
  select(type, parameter, index, everything()) %T>%
  write_csv(str_c(SAVE_LOCATION, '/INLA_models/models_',
                  run_name, '/paramMargin_', 
                  run_name, '.csv',
                  sep = ''))


bind_rows(
  INLA:::summary.inla(model_fit)$fixed %>%
    as_tibble(rownames = 'parameter'),
  
  INLA:::summary.inla(model_fit)$hyperpar %>%
    as_tibble(rownames = 'parameter'),
  
  model_fit$summary.random$habitat %>%
    as_tibble %>%
    mutate(parameter = 'habitat') %>%
    mutate(ID = as.character(ID)),
  
  model_fit$summary.random$site %>%
    as_tibble %>%
    mutate(parameter = 'site')  %>%
    mutate(ID = as.character(ID))
  
) %>%
  rename_with(~str_c('q_', str_remove(., 'quant')), .cols = ends_with('quant'))  %>%
  select(parameter, ID, everything()) %T>%
  write_csv(str_c(SAVE_LOCATION, '/INLA_models/models_',
                  run_name, '/summarizedParams', 
                  run_name, '.csv',
                  sep = ''))


message(str_c('Finished Parameter Extraction: ', Sys.time()))

#### Posterior Predictive Check - site & overall # of fish ####
message(str_c('Start Posterior Predictive Counts: ', Sys.time()))
waic_dic <- tibble(cpo = -1 * sum(log(model_fit$cpo$cpo)),
              cpo2 = 2 * cpo,
              waic = model_fit$waic$waic,
              waic_p.eff = model_fit$waic$p.eff,
              waic_sd = sd(model_fit$waic$local.waic),
              model_dic = model_fit$dic$dic,
              model_dic_p.eff = model_fit$dic$p.eff,
              model_deviance = model_fit$dic$mean.deviance,
              
              # null_dic = null_fit$dic$dic,
              # null_deviance = null_fit$dic$mean.deviance,
              # null_dic_p.eff = null_fit$dic$p.eff,
              # 
              # spatial_dic = spatial_fit$dic$dic,
              # spatial_deviance = spatial_fit$dic$mean.deviance,
              # spatial_dic_p.eff = spatial_fit$dic$p.eff
              ) %>%
  
  # mutate(pseudoR2 = 1 - model_deviance / null_deviance) %>%
  
  pivot_longer(cols = everything()) %T>%
  write_csv(str_c(SAVE_LOCATION, '/INLA_models/models_',
                  run_name, '/waicDic_', 
                  run_name, '.csv',
                  sep = ''))


site_preds <- site_poly %>%
  st_as_sf %>%
  select(-weight) %>%
  rename(Site = names) %>%
  # add_row(Site = 'all', geometry = st_as_sf(site_poly) %>% summarise %>% pull(geometry),
  #         .before = 1) %>%
  as_tibble() %>%
  rowwise %>%
  mutate(poly = list(as_Spatial(geometry))) %>%
  ungroup %>%
  left_join(fish_sf %>%
              as_tibble() %>%
              count(site) %>%
              rename(obs = n),
            by = c('Site' = 'site')) %>%
  mutate(spots = future_map(poly, ~ipoints(.x, mesh))) %>%
  mutate(predicted_number = future_map(spots, ~generate(model_fit, .x, pred_form, 
                                             n.samples = number_samples) %>%
                              as.numeric,
                            .options = furrr_options(seed = TRUE))) %>%
  select(Site, obs, predicted_number) %>%
  unnest(predicted_number) %T>%
  write_csv(str_c(SAVE_LOCATION, '/INLA_models/models_',
                  run_name, '/fishCountPostPred_', 
                  run_name, '.csv',
                  sep = ''))
message(str_c('Finished Posterior Predictive Counts: ', Sys.time()))

#### Predict model onto sites ####
message(str_c('Start Prediction onto Sites: ', Sys.time()))
set_crs_terra <- function(x, crs_proj){
  crs(x) <- crs_proj
  x
}

predictions_to_coords <- function(preds, coords){
  coords <- as(coords, 'SpatialPixelsDataFrame')
  coords@data <- as.data.frame(preds)
  coords
}

crop_pixels <- function(pix, pol){
  as(crop(as(pix, 'SpatialPoints'), pol), 'SpatialPixels')
}



dem_info <- list.files(path = str_c(DATA_folder,'COPE_Sites/DEM',sep='/'), pattern = 'tif$',
                            full.names = TRUE) %>%
  tibble(dem = .) %>%
  mutate(Site = str_extract(dem, 'bz17-[0-9KNSAB]+') %>% str_to_upper) %>%
  mutate(dem = map(dem, ~raster(.x)),
         dem = future_map(dem, ~projectRaster(.x, crs = newproj), .options = furrr_options(seed = TRUE))) %>%
  mutate(dem = map(dem, rast),
         dem = map(dem, ~trim(.x))) %>%
  rowwise %>%
  summarise(Site = Site, 
            rows = terra::nrow(dem),
            columns = terra::ncol(dem),
            .groups = 'drop') %>%
  left_join(translation_instructions, by = 'Site') %>%
  left_join(st_as_sf(site_poly) %>%
              as_tibble %>%
              select(-weight),
            by = c('Site' = 'names')) %>%
  rowwise %>%
  mutate(poly = list(as_Spatial(geometry))) %>%
  select(-geometry) %>%
  ungroup

#block overloads memory with 1k samples - 11:43 - 3:33
plan('sequential')
plan('multicore', workers = 3)

spatial_predictions <- dem_info %>%
  mutate(preds_pixels = future_pmap(list(poly, rows, columns),
                                    ~pixels(mesh, mask = ..1, 
                                            nx = ..3, ny = ..2) %>%
                                      crop_pixels(pol = ..1)
  )) %>%
  # mutate(model_pred = future_map(preds_pixels, ~predict(model_fit, #THIS NEEDS TO CHANGE TO GENERATE
  #                                                       .x,
  #                                                       spat_pred_form,
  #                                                       n.samples = number_samples),
  #                                .options = furrr_options(seed = TRUE))) %>%
  
  mutate(model_pred = future_map2(preds_pixels, movement, ~generate(model_fit,
                                                                    .x,
                                                                    spat_pred_form,
                                                                    n.samples = number_samples) %>%
                                    predictions_to_coords(coords = .x) %>%
                                    stack %>%
                                    set_names(str_c('V', 1:number_samples, sep = '')) %>%
                                    raster::shift(., dx = .y[1], dy = .y[2]) %>%
                                    trim,
                                  .options = furrr_options(seed = TRUE)
                                  )) %>%
  mutate(model_pred = map(model_pred, ~rast(.x) %>%
                            set_crs_terra(., newproj)))%>%
  
  select(Site, model_pred) %>%
  # mutate(circle_mat = future_map(unsmoothed, 
  #                                ~raster::focalWeight(raster(.x$V1), sqrt(1/pi), 'circle'),
  #                                .options = furrr_options(seed = TRUE))) %>%
  # rowwise %>%
  # mutate(smoothed = list(map(names(unsmoothed), ~terra::focal(unsmoothed[[.x]], w = circle_mat, fun = 'sum') %>%
  #                              cover(y = unsmoothed[[.x]])) %>%
  #                          rast
  #                        )) %>%
  # ungroup %>%
  # select(-circle_mat) %>%
  # pivot_longer(cols = c(unsmoothed)) %>% #smoothed
  mutate(out_name = str_c(SAVE_LOCATION, '/INLA_models/models_',
                          run_name, '/', Site, '_unsmoothed_', 
                          run_name, '.tif',
                          sep = '')) %>%
  mutate(write = future_map2(model_pred, out_name, ~writeRaster(.x, .y, overwrite = TRUE))) %>%
  select(-write)
message(str_c('Finished Prediction onto Sites: ', Sys.time()))

#### Posterior Samples ####
message(str_c('Start sampling posterior point processes: ', Sys.time()))
sample.lgcp_jds <- function(mesh, predictions, samplers){
  big_matrix_parallel <- function(x_mat, y_vec, shards = detectCores()){
    #Not used but may be helpful
    indicies <- sort(1:nrow(x_mat) %% shards + 1)
    
    y <- foreach(i = 1:max(indicies), .combine = c) %dopar% {
      as.vector(x_mat[indicies == i, ] %*% y_vec)
    }
    
    unlist(y)
  }
  
  predictions <- predictions %>%
    st_as_sf %>%
    st_join(st_as_sf(samplers), st_within) %>%
    mutate(inside_site = !is.na(names)) %>%
    select(-weight, -names) %>%
    as('Spatial')
  
  n.samples <- 1
  loglambda <- as.vector(log(predictions$mean))
  in_site <- as.vector(predictions$inside_site)
  mesh$crs <- NULL
  
  
  strategy <- "rectangle"
  internal.crs <- CRS(as.character(NA))
  mesh$crs <- NULL
  
  #Make the area just the area inside sites not outside too
  loc <- mesh$loc
  xmin <- min(loc[, 1])
  xmax <- max(loc[, 1])
  ymin <- min(loc[, 2])
  ymax <- max(loc[, 2])
  area <- (xmax - xmin) * (ymax - ymin)
  
  # loc <- mesh$loc[,1:2] %>% 
  #   set_colnames(c('x', 'y')) %>% 
  #   as_tibble %>%
  #   st_as_sf(coords = c('x', 'y')) %>%
  #   st_join(st_as_sf(site), st_within) %>%
  #   dplyr::select(-weight) %>%
  #   rename(Site = names) %>%
  #   filter(!is.na(Site)) %>%
  #   as_tibble %>%
  #   mutate(x = map_dbl(geometry, ~.x[1]),
  #          y = map_dbl(geometry, ~.x[2])) %>%
  #   group_by(Site) %>%
  #   summarise(across(c(x, y), list(min = min, max = max)),
  #             .groups = 'drop') %>%
  #   mutate(area = (x_max - x_min) * (y_max - y_min))
  #  
  # xmin <- loc$x_min
  # xmax <- loc$x_max
  # ymin <- loc$y_min
  # ymax <- loc$y_max
  # area <- loc$area
  
  lambda_max <- max(loglambda[in_site])
  
  # Npoints <- rpois(1, lambda = area * exp(lambda_max)) #PROBLEMS HERE MAYBE
  Npoints <- ifelse(0.999*.Machine$integer.max < area * exp(lambda_max), #If the mean is too big use a normal approximation
                    round(rnorm(1, mean = area * exp(lambda_max), sd = sqrt(area * exp(lambda_max)))),
                    rpois(1, lambda = area * exp(lambda_max)))
  
  shards <- Npoints %/% 1e8 + 1
  n_vec <- rep(Npoints %/% shards, shards)
  if(Npoints %% shards != 0){n_vec[1:(Npoints %% shards)] <- n_vec[1:(Npoints %% shards)] + 1}
  
  
  # ret_store <- vector('list', shards)
  # z <- 0
  # for(N in n_vec){
  #   if(z %% floor(0.01 * shards) == 0 & verbose){time_in <- Sys.time()}
  #   z <- z + 1
  #   
  #   points <- sp::SpatialPoints(cbind(x = runif(n = N, min = xmin, max = xmax), 
  #                                     y = runif(n = N, min = ymin, max = ymax)),
  #                               proj4string = internal.crs)
  #   proj <- INLA::inla.mesh.project(mesh, points)
  #   
  #   lambda_ratio <- exp(as.vector(proj$A %*% loglambda) - lambda_max)
  #   
  #   keep <- proj$ok & (runif(N) <= lambda_ratio)
  #   ret <- points[keep]
  #   waste_ratio <- sum(keep)/length(keep)
  #   
  #   
  #   proj4string(ret) <- CRS(as.character(NA))
  #   ret_store[[z]] <- ret[unique(which(!is.na(over(ret, samplers)), arr.ind = TRUE)[,1])]
  #   
  #   if(z %% floor(0.01 * shards) == 0 & verbose){
  #     time_spent <- Sys.time() - time_in
  #     message(str_c('Finished ', 'shards ', z - floor(0.01 * shards), ' - ', z, 
  #                   ' Time spent: ', round(time_spent, 3), ' ', attr(time_spent, 'units')))
  #   }
  # }
  # 
  # ret <- do.call(rbind, ret_store)
  
  p <- progressor(along = n_vec)
  
  ret <- foreach(N = n_vec, .combine = rbind) %dorng% {
    points <- sp::SpatialPoints(cbind(x = runif(n = N, min = xmin, max = xmax),
                                      y = runif(n = N, min = ymin, max = ymax)),
                                proj4string = internal.crs)
    proj <- INLA::inla.mesh.project(mesh, points)
    
    lambda_ratio <- exp(as.vector(proj$A[,in_site] %*% loglambda[in_site]) - lambda_max)
    
    keep <- proj$ok & (runif(N) <= lambda_ratio)
    rm(proj)
    
    
    out <- points[keep]
    rm(points)
    
    waste_ratio <- sum(keep)/length(keep)
    
    proj4string(out) <- CRS(as.character(NA))
    
    p()
    
    
    if(sum(keep) > 0){
      out <- out[unique(which(!is.na(over(out, samplers)), arr.ind = TRUE)[,1])]
    }
    
    out
  }
  
  
  ret
}

site_sf <- site_poly %>%
  st_as_sf %>%
  select(-weight) %>%
  rename(Site = names) %>%
  left_join(translation_instructions, by = 'Site') %>%
  rowwise %>%
  mutate(geometry = geometry + movement) %>%
  ungroup %>%
  select(-movement) %>%
  st_sf(crs = newproj) %T>%
  st_write(str_c(SAVE_LOCATION, '/INLA_models/models_',
                 run_name, '/siteOutline_', run_name, '.shp',
                 sep = ''),
           quiet = TRUE,
           delete_dsn = file.exists(str_c(SAVE_LOCATION, '/INLA_models/models_',
                                          run_name, '/siteOutline_', run_name, '.shp',
                                          sep = '')))

mesh_prediction <- predict(model_fit,
                           vertices(mesh),
                           spat_pred_form,
                           n.samples = number_samples)

# posterior_sample_test <- sample.lgcp_jds(mesh = mesh,
#                                     predictions = mesh_prediction,
#                                     samplers = site_poly)

plan('sequential')
plan(list(tweak(multicore, workers = number_workers), sequential))
# plan(list(tweak(multicore, workers = number_workers), 
#           tweak(multicore, workers = detectCores() - number_workers - 1)))

posterior_sample <- tibble(sim = 1:pp_samples) %>%
  mutate(posterior_sample = future_map(sim, function(x) sample.lgcp_jds(mesh = mesh,
                                                             predictions = mesh_prediction,
                                                             samplers = site_poly),
                                       .options = furrr_options(seed = TRUE))) %>%
  rowwise %>%
  summarise(sample = sim,
            as_tibble(posterior_sample),
            .groups = 'drop') %>%
  st_as_sf(coords = c('x', 'y')) %>%
  st_join(st_as_sf(site_poly), st_within) %>%
  select(-weight) %>%
  rename(Site = names) %>%
  left_join(translation_instructions, by = 'Site') %>%
  rowwise %>%
  mutate(geometry = geometry + movement) %>%
  ungroup %>%
  select(-movement) %>%
  st_sf(crs = newproj) %T>%
  st_write(str_c(SAVE_LOCATION, '/INLA_models/models_',
                 run_name, '/pointProcessSample_', run_name, '.shp',
                 sep = ''),
           quiet = TRUE,
           delete_dsn = file.exists(str_c(SAVE_LOCATION, '/INLA_models/models_',
                                          run_name, '/pointProcessSample_', run_name, '.shp',
                                          sep = '')))

message(str_c('Finished sampling posterior point processes: ', Sys.time()))


# Below transfer to post-processing script to analyze all together
# #### Posterior Predictive - secondary measures ####
# message(str_c('Start measuring site level posterior point processes: ', Sys.time()))
# 
# site_point_processes <- full_join(
#   fish_sf %>%
#     nest(obs = geometry),
#   
#   posterior_sample %>%
#     nest(sim = geometry),
#   
#   by = c('site' = 'Site')  
# ) %>%
#   nest(simulations = c(sample, sim)) %>%
#   full_join(site_sf, by = c('site' = 'Site')) %>%
#   rename(Site = site) %>%
#   
#   mutate(nsim = map_int(simulations, nrow),
#          obs = map(obs, ~as(.x, 'Spatial')),
#          site = map(geometry, ~as(.x, 'Spatial')),
#          simulations = map(simulations, ~pull(.x, sim))) %>%
#   select(-geometry) %>%
#   rowwise %>%
#   mutate(simulations = list(map(simulations, ~as(.x, 'Spatial')))) %>%
#   
#   mutate(site = list(as(site, 'owin')),
#          obs = list(ppp(x = coordinates(obs)[,1],
#                         y = coordinates(obs)[,2],
#                         window = site)),
#          simulations = list(map(simulations, ~ppp(x = coordinates(.x)[,1],
#                                                   y = coordinates(.x)[,2],
#                                                   window = site))
#          )) %>%
#   
#   ungroup %>%
#   full_join(expand_grid(Site = .$Site, 
#                         metric = c('Kest','Lest', 'pcf', 'Fest', 'Gest', 'Jest')),
#             by = 'Site') %>%
#   select(Site, metric, everything()) %>%
#   mutate(model = future_pmap(list(obs, simulations, nsim, metric), 
#                              ~envelope(..1, ..4, simulate = ..2, nsim = (..3 - 1)),
#                              .progress = TRUE, .options = furrr_options(seed = TRUE))) %>%
#   pivot_wider(names_from = 'metric', 
#               values_from = 'model') %>%
#   select(-obs, -simulations, -site) %T>%
#   write_rds(str_c(SAVE_LOCATION, '/INLA_models/models_',
#                   run_name, '/pointProcessSiteEnvelopes_', run_name, '.rds',
#                   sep = ''), 
#             compress = 'xz')
# message(str_c('Finished measuring site level posterior point processes: ', Sys.time()))  
# 
# ## All Sites together 
# 
# message(str_c('Start measuring overall posterior point processes: ', Sys.time()))
# 
# site_win <- site_sf %>%
#   select(geometry) %>%
#   summarise() %>%
#   as_Spatial() %>%
#   as('owin')
# 
# observed_ppp <- fish_sf %>%
#   select(geometry) %>%
#   summarise() %>%
#   as_Spatial() %>%
#   coordinates() %>%
#   ppp(x = .[,1],
#       y = .[,2],
#       window = site_win)
# 
# 
# posterior_sample_ppp <- posterior_sample %>%
#   select(sample, geometry) %>%
#   group_by(sample) %>%
#   group_split() %>%
#   future_map(~as_Spatial(.) %>% coordinates(), .options = furrr_options(seed = TRUE)) %>%
#   future_map(~ppp(x = .x[,1],
#                   y = .x[,2],
#                   window = site_win),
#              .options = furrr_options(seed = TRUE))
#  
# 
# overall_point_processes <- tibble(metric = c('Kest','Lest', 'pcf', 'Fest', 'Gest', 'Jest')) %>%
#   mutate(envelope = future_map(metric, ~envelope(observed_ppp, .x, 
#                                                  simulate = posterior_sample_ppp, 
#                                                  nsim = length(posterior_sample_ppp) - 1))) %T>%
#   write_rds(str_c(SAVE_LOCATION, '/INLA_models/models_',
#                   run_name, '/pointProcessEnvelopes_', run_name, '.rds',
#                   sep = ''), 
#             compress = 'xz')
# 
# message(str_c('Finished measuring overall posterior point processes: ', Sys.time()))
