#### Code to set up the inlabru model as rds files which can then be run as multiple slurm scripts simultaniously ####

#### Technical TODO List ####


#### Model TODO List ####
##TODO is mesh too fine??
##TODO set up site CV


if(Sys.info()['sysname'] == 'Windows'){
  model_choice <- 'c5'
} else {
  args <- commandArgs(trailingOnly=TRUE)
  model_choice <- args[1]
}

print(paste('Start LGCP Model for everything', Sys.time(), sample(1000, 1), sep = ': '))

#### Set up computer ####
# rm(list = ls())

#### Libraries ####
suppressMessages(library(raster))
suppressMessages(library(terra))
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(ggforce))
suppressMessages(library(corrr))
suppressMessages(library(sf))
suppressMessages(library(inlabru))
suppressMessages(library(INLA))
suppressMessages(library(rgeos))
suppressMessages(library(patchwork))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
suppressMessages(library(MuMIn))
suppressMessages(library(RStoolbox))

newproj<-"+proj=utm +zone=16Q +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
E_scale <- 1
PX_EXPAND <- if_else(Sys.info()['sysname'] == 'Windows', 128, 355) #1028
make_plots <- TRUE
write_output <- if_else(Sys.info()['sysname'] == 'Windows', FALSE, TRUE)
min_weight <- 0.05

# https://www.muscardinus.be/2018/07/inlabru-bru/ - marginal plots. use of bru 

#### Set up Computer ####
save_suffix <- 'All'
COMPUTER <- if_else(Sys.info()['nodename'] == 'JASONDELL', 'laptop', 
                    if_else(Sys.info()['nodename'] == 'jdselwyn', 'Gawain', 'HPC'))
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

#### Functions ####
sfill <- function(data, where = NULL) {
  #Function to extend factor variable across mesh based on NN
  #https://github.com/fbachl/inlabru/commit/2b6bff077eeda35e471e1b3960eeb1043a93e6b0#
  
  if ( is.null(where) ) { where = data }
  # vallist = list()
  registerDoParallel(cores = 1)
  vallist <- foreach(k=1:ncol(data@data), .packages = 'sp') %dopar% {
    
    dpoints = sp::SpatialPoints(data)
    vals = data@data[,k]
    dpoints = dpoints[!is.na(vals),]
    vals = vals[!is.na(vals)]
    
    data.ow = spatstat::owin(range(coordinates(dpoints)[,1]), range(coordinates(dpoints)[,2]))
    data.ppp = spatstat::as.ppp(coordinates(dpoints), data.ow)
    where.ow = spatstat::owin(range(coordinates(where)[,1]), range(coordinates(where)[,2]))
    where.ppp = spatstat::as.ppp(coordinates(where), where.ow)
    
    nn = spatstat::nncross(where.ppp, data.ppp)[,"which"]
    vals[nn]
  }
  
  
  ret = data.frame(do.call(cbind, vallist))
  colnames(ret) = colnames(data@data)
  ret = sp::SpatialPixelsDataFrame(where, data = ret)
}

#### Read in Data ####
message(str_c('Start Reading in data: ', Sys.time()))

shoals <- st_read(str_c(INTERMEDIATE_FILES, "/Topography/shoal_mosaic_", model_choice,".shp"))

topography <- rast(str_c(INTERMEDIATE_FILES, "/Topography/reef_mosaic_", model_choice, ".tif"))


overall_stats <- read_csv(str_c(SAVE_LOCATION, '/Topography/overall_metrics_', model_choice, '.csv'),
                          col_types = cols(
                            metric = col_character(),
                            mean = col_double(),
                            sd = col_double()
                          ))

site_stats <- read_csv(str_c(SAVE_LOCATION, '/Topography/site_metrics_', model_choice, '.csv'),
                       col_types = cols(
                         .default = col_double(),
                         Site = col_character()
                       ))

message(str_c('Finish Reading in data: ', Sys.time()))

#### Preprocess data ####
## Extract site boundaries from rasters
if(!file.exists(str_c(INTERMEDIATE_FILES, '/INLA_models/', 'site.rds', sep = ''))){
  site_borders <- topography$depth %>%
    is.na %>% 
    equals(0) %>% 
    classify(matrix(c(0, NA), ncol = 2)) %>% 
    as.polygons()
  
  tmpf2 <- tempfile()
  writeVector(site_borders, tmpf2, overwrite=TRUE)
  site <- st_read(tmpf2, crs = newproj, quiet = TRUE) %>%
    rename(weight = depth) %>%
    mutate(weight = 1) %>%
    st_simplify(preserveTopology = TRUE, dTolerance = 0.05) %>%
    
    st_cast("POLYGON") %>% 
    mutate(area = st_area(geometry)) %>%
    arrange(-area) %>%
    slice(1:12) %>%
    select(-area) %>%
    
    st_set_crs(NA_character_) %>%
    as('Spatial')
  
  
  site_names <- rgeos::gCentroid(site, byid = TRUE) %>%
    as.data.frame %>%
    terra::extract(topography$site, .) %>%
    as_tibble %>%
    mutate(site = case_when(site == 0 ~ 'BZ17-0A',
                            site == 1 ~ 'BZ17-0B',
                            site == 100 ~ 'BZ17-100N',
                            site == -100 ~ 'BZ17-100S',
                            site == 10000 ~ 'BZ17-10KN',
                            site == 1000 ~ 'BZ17-1KN',
                            site == -1000 ~ 'BZ17-1KS',
                            site == 500 ~ 'BZ17-500N',
                            site == -500 ~ 'BZ17-500S',
                            site == 5000 ~ 'BZ17-5KN',
                            site == -5000 ~ 'BZ17-5KS',
                            site == -60 ~ 'BZ17-60S',
                            TRUE ~ NA_character_))
  
  site$names <- site_names$site
  
  write_rds(site, str_c(INTERMEDIATE_FILES, '/INLA_models/', 'site.rds', sep = ''), compress = 'gz')
} else {
  site <- read_rds(str_c(INTERMEDIATE_FILES, '/INLA_models/', 'site.rds', sep = ''))
}

#### Read in and process habitat data ####
if(!file.exists(str_c(INTERMEDIATE_FILES, '/INLA_models/', 'topography.rds', sep = ''))){
  
  ## Process Habitat Type data ##
  message(str_c('Start processing habitat data: ', Sys.time()))
  habitat_type <- topography %>%
    subset(c('classified_habitat_unsmoothed')) %>%
    round %>%
    raster %>%
    set_names('habitat') %>%
    as('SpatialGridDataFrame')
  habitat_type$habitat <- case_when(habitat_type$habitat == 1 ~ 'Sand',
                                    habitat_type$habitat == 2 ~ 'Reef',
                                    TRUE ~ NA_character_)
  habitat_type$habitat <- factor(habitat_type$habitat)
  crs(habitat_type) <- NA_character_  
  
  # habitat_type <- sfill(habitat_type, where = pixels(mesh, mask = FALSE, nx = PX_EXPAND, ny = PX_EXPAND))
  # habitat_type <- sfill(habitat_type)
  message(str_c('Finish processing habitat data: ', Sys.time()))
  
  ## Process Site Data ##
  message(str_c('Start processing site data: ', Sys.time()))
  site_selection <- topography %>%
    subset(c('site')) %>%
    round %>%
    raster %>%
    set_names('site') %>%
    as('SpatialGridDataFrame')
  
  site_selection$site <- case_when(site_selection$site == 0 ~ 'BZ17-0A',
                                   site_selection$site == 1 ~ 'BZ17-0B',
                                   site_selection$site == 100 ~ 'BZ17-100N',
                                   site_selection$site == -100 ~ 'BZ17-100S',
                                   site_selection$site == 10000 ~ 'BZ17-10KN',
                                   site_selection$site == 1000 ~ 'BZ17-1KN',
                                   site_selection$site == -1000 ~ 'BZ17-1KS',
                                   site_selection$site == 500 ~ 'BZ17-500N',
                                   site_selection$site == -500 ~ 'BZ17-500S',
                                   site_selection$site == 5000 ~ 'BZ17-5KN',
                                   site_selection$site == -5000 ~ 'BZ17-5KS',
                                   site_selection$site == -60 ~ 'BZ17-60S',
                                   TRUE ~ NA_character_)
  site_selection$site <- factor(site_selection$site)
  crs(site_selection) <- NA_character_  
  
  # site_order <- levels(site_selection$site)
  
  # site_selection <- sfill(site_selection, where = pixels(mesh, mask = FALSE, nx = PX_EXPAND, ny = PX_EXPAND))
  # site_selection <- sfill(site_selection)
  # site_selection$site <- site_order[site_selection$site]
  # site_selection$site <- factor(site_selection$site)
  message(str_c('Finish processing site data: ', Sys.time()))
  
  ## PCA To create fine-scale complexity metric ##
  for_pca <- topography %>%
    brick %>% 
    subset(c('rugosity', 'vectorDispersion',
             'relief', 'slope')) 
  
  complexity_pca <- rasterPCA(for_pca, nSamples = 1000000, spca = TRUE, maskCheck = TRUE)
  # complexity_pca <- rasterPCA(for_pca, nSamples = 100, spca = TRUE, maskCheck = TRUE)
  
  library(exactextractr)
  
  site_stats <- exact_extract(complexity_pca$map, st_as_sf(site), c('mean', 'stdev')) %>%
    as_tibble %>%
    rename_with(~str_c(str_extract(., 'PC[0-9]+'), '.', str_extract(., 'mean|std'))) %>%
    rename_with(~str_remove(., 't')) %>%
    mutate(Site = site$names) %>%
    left_join(site_stats, ., by = 'Site') %T>%
    write_csv(str_c(INTERMEDIATE_FILES, '/Topography/site_metrics_', model_choice, '.csv'))
  
  variance_explained <- ((complexity_pca$model$sdev^2)/sum(complexity_pca$model$sdev^2))
  
  complexity_pca_pctVar <- variance_explained %>%
    enframe %>%
    ggplot(aes(x = name, y = value)) +
    geom_col() +
    geom_text(aes(y = 1.1 * value, label = round(value, 4))) +
    scale_y_continuous(labels = scales::percent_format(1)) +
    labs(y = 'Percent Variance Explained',
         x = NULL) +
    theme_classic()
  ggsave(str_c(SAVE_LOCATION, '/Topography/complexity_pca_var.png'), plot = complexity_pca_pctVar, height = 7, width = 7)
  
  pca_load_tibble <- loadings(complexity_pca$model) %>%
    as.data.frame.matrix() %>%
    as_tibble(rownames = 'variable')
  write_csv(pca_load_tibble, str_c(SAVE_LOCATION, '/Topography/complexity_pca_loadings.csv'))
  
  pca_loadings <- pca_load_tibble %>%
    ggplot(aes(x = Comp.1, y = Comp.2)) +
    geom_segment(aes(xend = 0, yend = 0)) +
    geom_text(aes(label = variable)) +
    labs(x = str_c('PC1 (', round(variance_explained[1]*100, 2), '%)'),
         y = str_c('PC2 (', round(variance_explained[2]*100, 2), '%)')) +
    theme_classic()
  ggsave(str_c(SAVE_LOCATION, '/Topography/complexity_pca_loadings.png'), plot = pca_loadings, height = 7, width = 7)
  
  # topography$complexity <- rast(complexity_pca$map)$layer.1
  
  ## Process Topography Data ##
  message(str_c('Start processing topography data: ', Sys.time()))
  ## Isolate to topography of interest and scale/center
  topography <- topography %>%
    brick %>% 
    subset(c('depth', 'moran', 'viewshed', 'boundaryDistance')) %>% 
    c(complexity_pca$map$PC1) %>%
    stack %>%
    set_names(c(names(.)[1:4], 'complexity')) %>%
    as('SpatialGridDataFrame')
  crs(topography) <- NA_character_
  
  message(str_c('Finish preprocessing data: ', Sys.time()))
  
  #### Exploratory Analysis ####
  
  ## Pairs plot
  # pairs_plot <-  topography@data %>%
  #   as_tibble %>%
  #   na.omit %>%
  #   ggplot(aes(x = .panel_x, y = .panel_y)) + 
  #   geom_smooth(method = 'lm', formula = y ~ x) +
  #   geom_point(alpha = 0.2, shape = 16, size = 0.5) +
  #   geom_autodensity() +
  #   facet_matrix(vars(everything()), layer.diag = 3, grid.y.diag = FALSE)
  # ggsave(str_c(INTERMEDIATE_FILES, '/INLA_models/exploratory_pairs_plot.png'), plot = pairs_plot, height = 15, width = 15)
  
  ## Histograms
  histograms <- topography@data %>%
    as_tibble %>%
    na.omit %>%
    pivot_longer(cols = everything()) %>%
    ggplot(aes(x = value)) +
    geom_histogram(bins = 50) +
    facet_wrap(~name, scales = 'free')
  ggsave(str_c(SAVE_LOCATION, '/INLA_models/exploratory_histogram_plot.png'), plot = histograms, height = 14, width = 14)
  
  ## Transform and z-score scale ind vars 
  topography@data$boundaryDistance <- log(topography@data$boundaryDistance + 1, base = 2)
  topography@data$moran <- log(topography@data$moran + 1, base = 2)
  # topography@data$relief <- log(topography@data$relief + 1, base = 2)
  # topography@data$rugosity <- log(topography@data$rugosity, base = 2)
  # 
  # histograms <- topography@data %>%
  #   as_tibble %>%
  #   na.omit %>%
  #   pivot_longer(cols = everything()) %>%
  #   ggplot(aes(x = value)) +
  #   geom_histogram(bins = 50) +
  #   facet_wrap(~name, scales = 'free')
  # ggsave(str_c(SAVE_LOCATION, '/INLA_models/exploratory_histogram_plot_transformed.png'), plot = histograms, height = 14, width = 14)
  
  
  topography@data <- as.data.frame(scale(topography@data, center = TRUE, scale = TRUE))
  
  # histograms <- topography@data %>%
  #   as_tibble %>%
  #   na.omit %>%
  #   pivot_longer(cols = everything()) %>%
  #   ggplot(aes(x = value)) +
  #   geom_histogram(bins = 50) +
  #   facet_wrap(~name, scales = 'free')
  # ggsave(str_c(SAVE_LOCATION, '/INLA_models/exploratory_histogram_plot_zScored.png'), plot = histograms, height = 14, width = 14)
  
  ## Correlations among vars
  corrs <- topography@data %>%
    as_tibble %>%
    na.omit %>%
    correlate() %>%
    rearrange %>%
    shave
  
  fashion(corrs)
  
  cor_plot <- rplot(corrs, print_cor = TRUE)
  ggsave(str_c(SAVE_LOCATION, '/INLA_models/exploratory_correlation_plot.png'), plot = cor_plot, height = 7, width = 7)
  
  
  #### Extrapolate covariates to mesh buffer area ####
  message(str_c('Start filling in topography outside site: ', Sys.time()))
  
  # topography <- sfill(topography, where = pixels(mesh, mask = FALSE, nx = PX_EXPAND, ny = PX_EXPAND))
  # topography <- sfill(topography)
  message(str_c('Finish filling in topography outside site: ', Sys.time()))
  
  
  # message(str_c('Start plotting expanded universe: ', Sys.time()))
  # pdf(str_c(SAVE_LOCATION, "/Topography/expanded_reef_mosaic_", model_choice,".pdf"), height = 10, width = 10, onefile = TRUE)
  # tmp <- rast(brick(topography))
  # for(i in 1:length(names(tmp))){
  #   plot(tmp[[i]], main = names(tmp)[i])
  #   plot(site, add = TRUE)
  # }
  # plot(rast(raster(habitat_type)), main = 'Habitat Type')
  # plot(site, add = TRUE)
  # 
  # plot(rast(raster(site_selection)), main = 'Site')
  # plot(site, add = TRUE)
  # dev.off()
  # message(str_c('Finish plotting expanded universe: ', Sys.time()))
  
  write_rds(topography, str_c(INTERMEDIATE_FILES, '/INLA_models/', 'topography.rds', sep = ''), compress = 'gz')
  write_rds(habitat_type, str_c(INTERMEDIATE_FILES, '/INLA_models/', 'habitat_type.rds', sep = ''), compress = 'gz')
  write_rds(site_selection, str_c(INTERMEDIATE_FILES, '/INLA_models/', 'site_selection.rds', sep = ''), compress = 'gz')
  
} else {
  topography <- read_rds(str_c(INTERMEDIATE_FILES, '/INLA_models/', 'topography.rds', sep = ''))
  habitat_type <- read_rds(str_c(INTERMEDIATE_FILES, '/INLA_models/', 'habitat_type.rds', sep = ''))
  site_selection <- read_rds(str_c(INTERMEDIATE_FILES, '/INLA_models/', 'site_selection.rds', sep = ''))
}


#### Model Set Up ####
message(str_c('Start setting up exhaustive model search: ', Sys.time()))

#### Create Exhaustive Search ####
model_building_blocks <- names(topography) %>%
  tibble(variable = .) %>%
  # sample_n(5) %>%
  mutate(variable_component = str_c("beta.", variable, "(main = continuous_covars(x, y, '", 
                                    variable, "'), model = 'linear', mean.linear = 0, prec.linear = 0.01)", 
                                    sep = '')) %>%
  add_row(variable = '1', variable_component = "spatialSmooth(main = coordinates, model = matern)") %>%
  # add_row(variable = '1', variable_component = "Intercept") %>% #change for spatial
  add_row(variable = 'habitat', 
          variable_component = "habitat(main = habitat_type, model = 'factor_contrast')") %>%
  # add_row(variable = 'site', 
  #         variable_component = "site(main = Site, model = 'iid', hyper = list(prec = list(prior = 'pc.prec', param = c(1.5, 0.001))))") %>%
  identity()




exhaustive_model_search <- model_building_blocks %>%
  select(variable) %>% 
  filter(variable != '1') %>%
  mutate(tmp = list(rnorm(nrow(model_building_blocks) - 1))) %>%
  unnest(tmp) %>%
  mutate(row_number = rep(1:(nrow(model_building_blocks) - 1), nrow(model_building_blocks) - 1)) %>%
  pivot_wider(names_from = 'variable',
              values_from = 'tmp') %>%
  lm(row_number ~ ., data = ., na.action = na.fail) %>%
  dredge(evaluate = FALSE) %>%
  enframe %>%
  mutate(name = as.integer(name),
         model_var = map_chr(value, ~as.character(.x)[2])) %>%
  select(-value) %>%
  rename(model_code = name) %>%
  mutate(model_var = str_remove(model_var, '^.*~ ')) %>%
  
  rowwise %>%
  mutate(model_components = list(str_split(model_var, ' \\+ ') %>% unlist)) %>%
  unnest(model_components) %>%
  # filter(!(model_components == '1' & str_detect(model_var, 'habitat'))) %>%
  left_join(model_building_blocks, by = c('model_components' = 'variable')) %>%
  group_by(model_code) %>%
  summarise(model_name = unique(model_var), 
            model_formula = str_c(variable_component, collapse = ' + '), .groups = 'rowwise') %>%
  mutate(model_formula = if_else(str_detect(model_name, 'site'), 
                                 str_c(model_formula, 'Intercept', sep = ' - '),
                                 str_c(model_formula, 'Intercept', sep = ' + '))) %>%
  
  mutate(model_formula = list(formula(str_c("coordinates ~ ", model_formula)))) %>%
  ungroup %>%
  mutate(model_save_name = str_c(INTERMEDIATE_FILES, '/INLA_models/models/', 
                                 str_replace_all(model_name, ' [\\+\\-] ', '_'), 
                                 '.rds')) %>%
  
  # sample_n(4) %>% #for testing
  
  sample_frac(size = 1, replace = FALSE) %>%
  # filter(model_code == max(model_code)) %>% #only fit full model
  mutate(job_number = 1:n())
message(str_c('Finished setting up exhaustive model search: ', Sys.time()))

#### Save outputs for use in other scripts ####
write_csv(select(exhaustive_model_search, -model_formula), 
          str_c(INTERMEDIATE_FILES, '/INLA_models/models_to_fit.csv', sep = ''))
write_rds(exhaustive_model_search, 
          str_c(INTERMEDIATE_FILES, '/INLA_models/', 'models_to_fit.rds', sep = ''), compress = 'gz')

