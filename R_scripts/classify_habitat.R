#### Code to Classify sites as sand/reef ####
rm(list=ls())

#on hpc #R/gcc/64/3.5.1

#### Parameter Options ####
COMPUTER<-'HPC' #Set up some logic here to check
save_suffix <- 'All' #justJason or All
train_test_split <- 0.7 #What should the split between training and testing be?
bayes_iter <- 200L #maximum number of iterations for bayesian optimization
bayes_noImprov <- 50L #if not improvement after 30 iterations then stop optimization
initial_grid_multiplier <- 5 #initial grid will be drawn from X times the number of parameters in the model

#### To Do ####
# 2) Reasonable parameter ranges (?)

#### Libraries ####
library(sf)
library(raster)
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
library(patchwork)

options(tidymodels.dark = TRUE)

#### Functions ####
zero_to_NA<-function(x){
  NAvalue(x$Red)<-0
  NAvalue(x$Green)<-0
  NAvalue(x$Blue)<-0
  x
}

ortho_to_dem_res <- function(x, y){
  out_name <- str_replace(x, '\\.tif', '_demRez.tif')
  
  if(file.exists(out_name)){
    out_name
  } else {
    x <- stack(x)
    
    y <- raster(y)
    
    out <- projectRaster(x, y, filename = out_name)
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
  plan(multiprocess)
  
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
                             ortho_file = list.files(path = str_c(DATA_folder,'COPE_Sites/Orthomosaics',sep='/'), pattern = 'orthomosaic.tif$', full.names = TRUE),
                             training_file = list.files(path = str_c(DATA_folder, classified_habitat, sep='/'), pattern = 'shp$', full.names = TRUE),
                             shoals_file = list.files(str_c(DATA_folder,'COPE_Sites/Shoals',sep='/'), recursive = TRUE, pattern = 'shp$', full.names = TRUE)) %>%
  mutate(Site = str_extract(shoals_file, 'BZ17-[0-9]*.[ABKNS]')) %>%
  dplyr::select(Site, everything()) %>%
  
  mutate(ortho_file = future_map2_chr(ortho_file, dem_file, ortho_to_dem_res)) #Change resolution of orthomosaic to match the DEM


#### Extract and Process Data ####
#Subsample points such that each site/habitat class has an ~equal number of points. This was done to ensure points come from all shapes in proportion to their proportion of the site/habitatClass 
classification_data <- all_site_locations %>%
  
  # slice(c(1, 5, 7)) %>% #for testing
  
  #Read in and process shape files
  dplyr::select(Site, ortho_file, training_file) %>%
  mutate(train_data_full = future_map(training_file, function(x) read_sf(x) %>%
                                        filter(classStr!='flag',classStr!='GCP') %>%
                                        mutate(class=if_else(classStr=='sand','C1','C2')) %>%
                                        mutate(class=factor(class)) %>%
                                        mutate(Shape=1:nrow(.),Shape=as.character(Shape)),
                                      .progress = TRUE)) %>%
  dplyr::select(-training_file) %>%
  unnest(train_data_full) %>%
  st_as_sf %>%
  
  mutate(area = future_map_dbl(geometry, function(x) st_sfc(x, crs = oldproj) %>%
                                 st_as_sf %>%
                                 st_transform(crs = newproj) %>%
                                 st_area,
                               .progress = FALSE)) %>%
  arrange(-area) %T>% #Sort to make the parallelization more efficient hopefully
  print() %>% 
  #Could chunk larger areas into multiple smaller ones to accelerate

  # Temp for coding
  # filter(area < 5) %>%
  
  #Read in the points
  mutate(points = future_map2(ortho_file, geometry, ~stack(.x) %>%
                                dropLayer(4) %>%
                                set_names(c('Red', 'Green', 'Blue')) %>%
                                zero_to_NA %>%
                                exact_extract(x = ., 
                                              y = st_as_sf(st_sfc(.y, 
                                                                  crs = oldproj))) %>%
                                pluck(1) %>%
                                as_tibble() %>%
                                filter(coverage_fraction == 1) %>%
                                dplyr::select(-coverage_fraction),
                              .progress = TRUE)) %T>%
  print() %>%
  
  mutate(points = future_map(points, na.omit)) %>%
  mutate(n_point = future_map_int(points, nrow)) %>%
  filter(n_point > 0) %>%
  
  #Figure our number of points per shape to keep
  group_by(Site, class) %>%
  mutate(total_points = sum(n_point)) %>%
  ungroup %>%
  mutate(prop_shape = n_point/total_points,
         min_points = min(total_points),
         points_saved = floor(prop_shape * min_points))  %>%
  
  #Subsample points
  mutate(points = future_map2(points, points_saved, function(x, y) sample_n(x, size = y), .progress = TRUE)) %>%
  unnest(points) %>%
  as_tibble %>%
  dplyr::select(class, Site, Red, Green, Blue) %>%
  mutate(across(where(is.character), as.factor)) %>%

  #Split into training and testing
  initial_split(prop = train_test_split, strata = class) %T>%
  write_rds(str_c(INTERMEDIATE_FILES, '/Habitat_Classification/', save_suffix, '_training.rds', sep = ''), compress = 'gz')

# classification_data <- read_rds(str_c(INTERMEDIATE_FILES, '/Habitat_Classification/', save_suffix, '_training.rds', sep = ''))

bind_rows(Train = training(classification_data),
          Test = testing(classification_data),
          .id = 'data_class') %>%
  group_by(data_class, Site, class) %>%
  summarise(n = n(), .groups = 'drop') %>%
  mutate(class = case_when(class == 'C1' ~ 'Sand',
                           class == 'C2' ~ 'Reef',
                           TRUE ~ 'Error')) %>%
  pivot_wider(names_from = c(data_class, class), values_from = n) %T>%
  write_csv(str_c(SAVE_LOCATION, '/Habitat_Classification/', save_suffix, '_training_testing_split.csv', sep = ''))

# classification_data <- training(classification_data) %>%
#   bind_rows(testing(classification_data)) %>%
#   sample_frac(0.01) %>%
#   initial_split(prop = 0.7, strata = class)

#### Preprocess ####
preprocess_data <- training(classification_data) %>%
  recipe(class ~ .) %>%
  step_YeoJohnson(all_numeric()) %>%
  step_normalize(all_numeric()) %>%
  step_pca(all_numeric(), num_comp = tune()) %>%
  step_dummy(Site)

# prep(preprocess_data) %>% juice

# preprocess_data_noPCA <- training(classification_data) %>%
#   recipe(class ~ .) %>%
#   # step_string2factor(Site) %>%
#   # step_naomit(all_predictors(), all_outcomes(), skip = TRUE) %>%
#   step_YeoJohnson(all_numeric()) %>%
#   step_normalize(all_numeric()) %>%
#   step_dummy(Site) 
# prepared_training <- prep(preprocess_data_noPCA)

mtry_max <- nlevels(training(classification_data)$Site) + 3

#### Model Set Up ####
#### Random Forest ####
RF_model <- rand_forest(mode = 'classification', mtry = tune(), trees = tune(), min_n = tune()) %>%
  set_engine("ranger")

RF_wflow <- workflow() %>%
  add_recipe(preprocess_data %>% finalize_recipe(tibble(num_comp = 3L))) %>%
  add_model(RF_model)

RF_param <- parameters(RF_wflow) %>%
  # finalize(x = juice(prepared_training) %>% select(-class)) #%>%
  # update(num_comp = num_comp(c(0L, 3L))) %>%
  update(mtry = mtry(c(1, mtry_max)))

#### k-Nearest Neighbor ####
KNN_model <- nearest_neighbor(mode = 'classification', neighbors = tune(), weight_func = tune()) %>%
  set_engine("kknn")

KNN_wflow <- workflow() %>%
  add_model(KNN_model) %>%
  add_recipe(preprocess_data)

KNN_param <- parameters(KNN_wflow) %>%
  update(num_comp = num_comp(c(0L, 3L)))

#### Boosted Random Forest ####
boostRF_model <- boost_tree(mode = 'classification', mtry = tune(), trees = tune(), min_n = tune(), 
           tree_depth = tune(), learn_rate = tune(), loss_reduction = tune(), sample_size = tune()) %>%
  set_engine("xgboost")

boostRF_wflow <- workflow() %>%
  add_model(boostRF_model) %>%
  add_recipe(preprocess_data %>% finalize_recipe(tibble(num_comp = 3L)))

boostRF_param <- parameters(boostRF_wflow) %>%
  # finalize(x = juice(prepared_training) %>% select(-class))
  # update(num_comp = num_comp(c(0L, 3L))) %>%
  update(mtry = mtry(c(1, mtry_max)))
  

#### Logistic Regression ####
GLM_model <- logistic_reg(mode = 'classification') %>%
  set_engine("glm")

GLM_wflow <- workflow() %>%
  add_model(GLM_model) %>%
  add_recipe(preprocess_data)

#GLM - no param for default
GLM_param <- parameters(GLM_wflow) %>%
  update(num_comp = num_comp(c(0L, 3L)))

#### Decision Tree ####
CART_model <- decision_tree(mode = 'classification', cost_complexity = tune(), tree_depth = tune(), min_n = tune()) %>%
  set_engine("rpart")

CART_wflow <- workflow() %>%
  add_model(CART_model) %>%
  add_recipe(preprocess_data)

CART_param <- parameters(CART_wflow) %>%
  update(num_comp = num_comp(c(0L, 3L)))

#### Support Vector Machine - radial basis function ####
RBF_model <- svm_rbf(mode = 'classification', cost = tune(), rbf_sigma = tune()) %>%
  set_engine("kernlab")

RBF_wflow <- workflow() %>%
  add_model(RBF_model) %>%
  add_recipe(preprocess_data)

RBF_param <- parameters(RBF_wflow) %>%
  update(num_comp = num_comp(c(0L, 3L)))

#### Support Vector Machine - polynomial basis function ####
POLY_model <- svm_poly(mode = 'classification', cost = tune(), degree = tune(), scale_factor = tune()) %>%
  set_engine("kernlab")

POLY_wflow <- workflow() %>%
  add_model(POLY_model) %>%
  add_recipe(preprocess_data)

POLY_param <- parameters(POLY_wflow) %>%
  update(num_comp = num_comp(c(0L, 3L)))

#### Mars Model ####
MARS_model <- mars(mode = 'classification', num_terms = tune(), prod_degree = tune(), prune_method = tune()) %>%
  set_engine("earth")

MARS_wflow <- workflow() %>%
  add_model(MARS_model) %>%
  add_recipe(preprocess_data %>% finalize_recipe(tibble(num_comp = 3L)))

MARS_param <- parameters(MARS_wflow) %>%
  # finalize(x = juice(prepared_training) %>% select(-class)) %>%
  update(prune_method = prune_method(values = str_subset(values_prune_method, 'cv', negate = TRUE))) %>%
  # update(num_comp = num_comp(c(0L, 3L))) %>%
  update(num_terms = num_terms(c(2L, mtry_max)))
  

#### Flexible Disciminant Analysis ####
fda_model <- discrim_flexible(mode = 'classification', num_terms = tune(), 
                                      prod_degree = tune(), prune_method = tune()) %>%
  set_engine('earth')

fda_wflow <- workflow() %>%
  add_model(fda_model) %>%
  add_recipe(preprocess_data %>% finalize_recipe(tibble(num_comp = 3L)))

fda_param <- parameters(fda_wflow) %>%
  update(num_terms = num_terms(c(2L, mtry_max))) %>%
  # finalize(x = juice(prepared_training) %>% select(-class)) %>%
  # update(num_comp = num_comp(c(0L, 3L))) %>%
  update(prune_method = prune_method(values = str_subset(values_prune_method, 'cv', negate = TRUE)))

#### Linear Discriminant Analysis ####
lda_model <- discrim_linear(mode = 'classification') %>% 
  set_engine('MASS')

lda_wflow <- workflow() %>%
  add_model(lda_model) %>%
  add_recipe(preprocess_data)

#lda - no param for default
lda_param <- parameters(lda_wflow) %>%
  update(num_comp = num_comp(c(0L, 3L)))

#### Regularized Discriminant Analysis ####
rda_model <- discrim_regularized(mode = 'classification', frac_common_cov = tune(), 
                                        frac_identity = tune()) %>%
  set_engine('klaR')

rda_wflow <- workflow() %>%
  add_model(rda_model) %>%
  add_recipe(preprocess_data)

rda_param <- parameters(rda_wflow) %>%
  update(num_comp = num_comp(c(0L, 3L)))

#### Naive Bayes ####
nb_model <- naive_Bayes(mode = 'classification', smoothness = tune(), Laplace = tune()) %>%
  set_engine('klaR')

nb_wflow <- workflow() %>%
  add_model(nb_model) %>%
  add_recipe(preprocess_data)

nb_param <- parameters(nb_wflow) %>%
  update(num_comp = num_comp(c(0L, 3L)))

#### C5 Model ####
c5_model <- C5_rules(mode = 'classification', trees = tune(), min_n = tune()) %>%
  set_engine('C5.0')

c5_wflow <- workflow() %>%
  add_model(c5_model) %>%
  add_recipe(preprocess_data)

c5_param <- parameters(c5_wflow) %>%
  update(trees = trees(range = c(1L, 100L))) %>%
  update(num_comp = num_comp(c(0L, 3L)))

#### Rules based model ####
rule_model <- rule_fit(mode = 'classification', mtry = tune(), trees = tune(), min_n = tune(),
                       tree_depth = tune(), learn_rate = tune(), loss_reduction = tune(),
                       sample_size = tune()) %>%
  set_engine('xrf')

rule_wflow <- workflow() %>%
  add_model(rule_model) %>%
  add_recipe(preprocess_data %>% finalize_recipe(tibble(num_comp = 3L)))

rule_param <- parameters(rule_wflow) %>%
  # finalize(x = juice(prepared_training) %>% select(-class)) %>%
  # update(num_comp = num_comp(c(0L, 3L))) %>%
  update(mtry = mtry(c(1, mtry_max)))
  

#### Bagged MARS model ####
bagMars_model <- bag_mars(mode = 'classification', num_terms = tune(), prod_degree = tune(),
                          prune_method = tune()) %>%
  set_engine('earth')

bagMars_wflow <- workflow() %>%
  add_model(bagMars_model) %>%
  add_recipe(preprocess_data)

bagMars_param <- parameters(bagMars_wflow) %>%
  update(prune_method = prune_method(values = str_subset(values_prune_method, 'cv', negate = TRUE))) %>%
  update(num_comp = num_comp(c(0L, 3L)))

#### Bagged tree based model ####
bagTree_model <- bag_tree(mode = 'classification', cost_complexity = tune(), tree_depth = tune(),
                          min_n = tune(), class_cost = tune()) %>%
  set_engine('rpart')

bagTree_wflow <- workflow() %>%
  add_model(bagTree_model) %>%
  add_recipe(preprocess_data)

bagTree_param <- parameters(bagTree_wflow) %>%
  update(num_comp = num_comp(c(0L, 3L)))

#### Fit models to parameter grid ####
# cross_val <- bootstraps(training(classification_data), times = 3)
# cross_val <- mc_cv(training(classification_data), prop = 3/4, times = 3)
cross_val <- vfold_cv(training(classification_data), v = 10, repeats = 5, strata = 'class')
eval_metrics <- metric_set(j_index, roc_auc, accuracy, f_meas, kap, mcc) 
# eval_metrics <- metric_set(accuracy)

plan('sequential')
registerDoFuture()
cl <- makeClusterPSOCK(detectCores())
plan('cluster', workers = cl)

hold <- parallel::clusterEvalQ(cl, {library(kernlab)})
hold <- parallel::clusterEvalQ(cl, {library(discrim)})
hold <- parallel::clusterEvalQ(cl, {library(rules)})
hold <- parallel::clusterEvalQ(cl, {library(baguette)})

parameter_tuning <- tibble(model_name = c('rf', 
                                          'knn', 
                                          'boostRF', 
                                          'glm',
                                          'cart', 
                                          'rbf_svm', 
                                          'poly_svm', 
                                          'mars',
                                          'fda',
                                          'lda',
                                          'rda', 
                                          'nb', 
                                          'c5', 
                                          'rule',
                                          'bag_mars',
                                          'bag_tree'),
                           wflow = list(RF_wflow, 
                                        KNN_wflow, 
                                        boostRF_wflow, 
                                        GLM_wflow,
                                        CART_wflow, 
                                        RBF_wflow, 
                                        POLY_wflow, 
                                        MARS_wflow,
                                        fda_wflow,
                                        lda_wflow,
                                        rda_wflow, 
                                        nb_wflow, 
                                        c5_wflow, 
                                        rule_wflow,
                                        bagMars_wflow,
                                        bagTree_wflow), 
                           params = list(RF_param, 
                                         KNN_param, 
                                         boostRF_param, 
                                         GLM_param,
                                         CART_param, 
                                         RBF_param, 
                                         POLY_param, 
                                         MARS_param,
                                         fda_param,
                                         lda_param,
                                         rda_param, 
                                         nb_param, 
                                         c5_param, 
                                         rule_param,
                                         bagMars_param,
                                         bagTree_param
                           )) %>%
  filter(!model_name %in% c('poly_svm', # - always has errors about votematrix given particular combinations of hyper-params
                            'rule' # -can't install package on HPC
                            )) %>%
  
  #temp for coding 
  # sample_n(3) %>%
  # slice(1, 2, 4, 10) %>%
  
  mutate(initial_grid_size = initial_grid_multiplier * map_int(params, nrow)) %>%
  mutate(initial_grid = pmap(list(params, initial_grid_size), 
                            ~possibly(~grid_max_entropy(..1, 
                                                  size = ..2,
                                                  iter = 1000) %>%
                                        distinct,
                                      otherwise = NULL)(..1, ..2, ..3)))  %>%
  
  mutate(initial_grid = map(initial_grid, sample_n, size = 3)) %>% #Temp for coding
  
  mutate(initial_tune = pmap(list(wflow, params, initial_grid),
                             ~possibly(~tune_grid(..1, 
                                                   resamples = cross_val,
                                                   param_info = ..2,
                                                   grid = ..3,
                                                   metrics = eval_metrics, 
                                                   control = control_grid(verbose = FALSE,
                                                                          save_pred = FALSE)),
                                       otherwise = NULL)(..1, ..2, ..3)
                             )) %>% 

  mutate(tuned_bayes = pmap(list(wflow, params, initial_tune), 
                            ~possibly(~tune_bayes(..1, 
                                                  resamples = cross_val,
                                                  param_info = ..2,
                                                  initial = ..3,
                                                  iter = bayes_iter,
                                                  metrics = eval_metrics,
                                                  control = control_bayes(no_improve = bayes_noImprov,
                                                                          uncertain = 5L,
                                                                          verbose = TRUE)),
                                      otherwise = NULL)(..1, ..2, ..3)
                            )) %>%
  
  filter(!map_lgl(tuned_bayes, is.null)) %>%
  mutate(metrics = map(tuned_bayes, collect_metrics),
         best_params = map(tuned_bayes, select_best, metric = 'j_index')) %T>%
  write_rds(str_c(INTERMEDIATE_FILES, '/Habitat_Classification/', 
                  save_suffix, '_parameter_tuning.rds', sep = ''),
            compress = 'gz')

# parameter_tuning <- read_rds(str_c(INTERMEDIATE_FILES, '/Habitat_Classification/', save_suffix, '_parameter_tuning.rds', sep = ''))

#### Plot iteration accuracy and accuracy by parameter for each model ####
eval_iterations <- parameter_tuning %>%
  select(model_name, metrics) %>%
  unnest(metrics) %>%
  
  ggplot(aes(x = .iter, y = mean, ymin = mean - std_err, ymax = mean + std_err)) +
  geom_linerange() +
  geom_point() +
  facet_grid(.metric ~ model_name, scales = 'free') +
  scale_y_continuous(labels=percent) +
  ylab('Accuracy') +
  xlab('Iteration') +
  theme_bw()
ggsave(str_c(SAVE_LOCATION, '/Habitat_Classification/', 
             save_suffix, '_iteration_eval.png', sep = ''), 
       plot = eval_iterations, height = 10, width = 10)


eval_param <- parameter_tuning %>%
  select(model_name, metrics) %>%
  mutate(metric_plots = map(metrics, ~.x %>%
                              mutate(across(-c(.iter, .metric, .estimator, mean, n, std_err, .config), as.character)) %>%
                              pivot_longer(cols = -c(.iter, .metric, .estimator, mean, n, std_err, .config),
                                           names_to = 'parameter',
                                           values_to = 'value') %>%
                              nest(data = -c(parameter)) %>%
                              mutate(plots = map2(data, parameter, ~.x %>%
                              {if(all(!is.na(parse_number(.$value)))) mutate(., value = parse_number(value)) else .} %>%
                                
                                ggplot(aes(x = value, y = mean, ymin = mean - std_err, 
                                           ymax = mean + std_err)) +
                                geom_errorbar(width = 0) +
                                geom_point() +
                                geom_smooth(method = 'gam') + #formula = 'y ~ s(x, bs = "cs")'
                                # geom_smooth(method = 'loess', formula = 'y ~ x') +
                                facet_grid(rows = vars(.metric), scales = 'free') +
                                labs(x = .y,
                                     y = 'Evaluation Metric') +
                                theme_bw()))),
         metric_plots = map2(metric_plots, model_name, ~wrap_plots(.x$plots) + plot_annotation(title = .y)))

pdf(str_c(SAVE_LOCATION, '/Habitat_Classification/', 
          save_suffix, '_param_eval.pdf', sep = ''),
    onefile = TRUE, height = 15, width = 15)
for(i in 1:nrow(eval_param)){
  print(eval_param$metric_plots[[i]])
}
dev.off()

#### Choose best parameters of each model - and fit to full training data ####
best_parameter_fits <- parameter_tuning %>%
  
  #Fit to full training data
  mutate(wflow = map2(wflow, best_params, finalize_workflow)) %>%
  select(model_name, best_params, wflow) %T>%
  write_rds(str_c(SAVE_LOCATION, '/Habitat_Classification/', 
                  save_suffix, '_top_param.rds', sep = ''))

# best_parameter_fits <- read_rds(str_c(INTERMEDIATE_FILES, '/Habitat_Classification/', save_suffix, '_fitted_training_best_parameters.rds', sep = ''))


#### Model assessment ####
hold <- parallel::clusterEvalQ(cl, {library(recipes)})

eval_met <- metric_set(roc_auc, sens, spec, accuracy, bal_accuracy, f_meas, j_index, kap, detection_prevalence, mcc, npv, ppv, precision, recall)

tuned_model_results <- best_parameter_fits %>%
  mutate(test_fits = map(wflow, ~last_fit(.x, classification_data,
                                          metrics = eval_met))) %>%
  select(-wflow) %>%
  unnest(test_fits) %>%
  filter(!map_lgl(.metrics, is.null),
         !map_lgl(.predictions, is.null)) %>%
  mutate(confusion = map(.predictions, ~conf_mat(.x, truth = class, estimate = .pred_class)),
         mcnemar = map(confusion, ~tidy(.x) %$% matrix(value, nrow = 2) %>% mcnemar.test %>% tidy %>% 
                                select(statistic, p.value) %>% rename_with(~str_c('mcnemar_', .)))) %>%
  select(model_name, best_params, .workflow, .predictions, confusion, .metrics, mcnemar) %>%
  mutate(.metrics = map(.metrics, ~select(.x, .metric, .estimate) %>% 
                          pivot_wider(names_from = .metric, values_from = .estimate))) %>%
  unnest(c(.metrics, mcnemar)) %T>%
  write_rds(path = str_c(SAVE_LOCATION, '/Habitat_Classification/', save_suffix, '_topParam_model_assessment.rds', sep = ''),
            compress = 'xz')

tuned_model_results %>%
  select(-.workflow, -.predictions) %>%
  write_rds(str_c(SAVE_LOCATION, '/Habitat_Classification/', save_suffix, '_topParam_model_assessment_small.rds', sep = ''))

colour_choices <- map(c('Reds', 'Greens', 'Blues', 'Oranges', 'Purples'), ~scales::brewer_pal(palette = .x)(4) %>% 
                        extract(-1)) %>%
  unlist %>% 
  sample(14) %>%
  set_names(c('rf', 'knn', 'boostRF', 'glm', 'cart', 'rbf_svm', 'mars', 'fda', 'lda', 
              'rda', 'nb', 'c5', 'bag_mars', 'bag_tree'))

model_comparison <- tuned_model_results %>%
  select(-best_params, -.workflow, -.predictions, -confusion, -mcnemar_p.value) %>%
  pivot_longer(cols = -model_name,
               names_to = 'metric',
               values_to = 'value') %>%
  mutate(model_name_sorted = reorder_within(model_name, value, metric)) %>%
  
  ggplot(aes(y = model_name_sorted, x = value, colour = model_name)) +
  geom_text(aes(label = model_name), show.legend = FALSE) +
  facet_wrap(~metric, scales = 'free') +
  scale_y_reordered() +
  labs(x = NULL, 
       y = NULL) +
  theme_classic() +
  theme(axis.text = element_text(colour = 'black'),
        strip.background = element_blank(),
        panel.background = element_rect(colour = 'black', size = 1))
ggsave(str_c(SAVE_LOCATION, '/Habitat_Classification/', save_suffix, '_model_selection.png', sep = ''), 
       plot = model_comparison, height = 15, width = 15)

#### Get accuracy stats per site ####
site_by_model_assessment <- tuned_model_results %>%
  select(model_name, .predictions) %>%
  mutate(.predictions = map(.predictions, ~bind_cols(.x, testing(classification_data) %>%
                                                      select(Site)) %>%
                              nest(.predictions = -c(Site)))) %>%
  unnest(.predictions) %>%
  mutate(confusion = map(.predictions, ~conf_mat(.x, truth = class, estimate = .pred_class)),
         roc_auc = map_dbl(.predictions, ~roc_auc(.x, class, .pred_C1) %>% pull(.estimate)),
         sens = map_dbl(.predictions, ~sens(.x, class, .pred_class) %>% pull(.estimate)),
         spec = map_dbl(.predictions, ~spec(.x, class, .pred_class) %>% pull(.estimate)),
         accuracy = map_dbl(.predictions, ~accuracy(.x, class, .pred_class) %>% pull(.estimate)),
         bal_accuracy = map_dbl(.predictions, ~bal_accuracy(.x, class, .pred_class) %>% pull(.estimate)),
         f_meas = map_dbl(.predictions, ~f_meas(.x, class, .pred_class) %>% pull(.estimate)),
         j_index = map_dbl(.predictions, ~j_index(.x, class, .pred_class) %>% pull(.estimate)),
         kap = map_dbl(.predictions, ~kap(.x, class, .pred_class) %>% pull(.estimate)),
         detection_prevalence = map_dbl(.predictions, ~detection_prevalence(.x, class, .pred_class) %>% pull(.estimate)), 
         mcc = map_dbl(.predictions, ~mcc(.x, class, .pred_class) %>% pull(.estimate)), 
         npv = map_dbl(.predictions, ~npv(.x, class, .pred_class) %>% pull(.estimate)), 
         ppv = map_dbl(.predictions, ~ppv(.x, class, .pred_class) %>% pull(.estimate)), 
         precision = map_dbl(.predictions, ~precision(.x, class, .pred_class) %>% pull(.estimate)), 
         recall = map_dbl(.predictions, ~recall(.x, class, .pred_class) %>% pull(.estimate)),
         mcnemar = map(confusion, ~tidy(.x) %$% matrix(value, nrow = 2) %>% mcnemar.test %>% tidy %>% select(statistic, p.value) %>% rename_with(~str_c('mcnemar_', .)))) %>%
  unnest(mcnemar) %>%
  select(-.predictions) %T>%
  write_rds(str_c(SAVE_LOCATION, '/Habitat_Classification/', save_suffix, '_topParam_model_bySite.rds', sep = ''),
            compress = 'gz')


site_plots <- site_by_model_assessment %>%
  select(-confusion, -mcnemar_p.value) %>%
  pivot_longer(cols = -c(model_name, Site),
               names_to = 'metric',
               values_to = 'value') %>%
  nest(data = -metric) %>%
  # mutate(data = map(data, ~mutate(.x, Site = reorder_within(Site, value, model_name)))) %>%
  mutate(plot = map(data, ~ggplot(.x, aes(x = Site, y = value, colour = model_name)) +
                      # geom_point(position = position_dodge(width = 0.1)) +
                      geom_text(aes(label = model_name), position = position_dodge(width = 0.1), show.legend = FALSE) +
                      # facet_wrap(~model_name, scales = 'free') +
                      scale_colour_manual(values = colour_choices) +
                      # scale_shape_manual() +
                      labs(x = NULL, 
                           y = NULL) +
                      theme_classic()))


site_plot2 <- site_by_model_assessment %>%
  select(-confusion, -mcnemar_p.value) %>%
  pivot_longer(cols = -c(model_name, Site),
               names_to = 'metric',
               values_to = 'value') %>%
  mutate(Site = str_remove(Site, 'BZ17-')) %>%
  
  ggplot(aes(x = Site, y = value, colour = model_name)) +
  # geom_point(position = position_dodge(width = 0.1)) +
  geom_text(aes(label = model_name), position = position_dodge(width = 0.2), show.legend = FALSE) +
  # facet_wrap(~model_name, scales = 'free') +
  scale_colour_manual(values = colour_choices) +
  # scale_shape_manual() +
  facet_wrap(~metric, scales = 'free') +
  labs(x = NULL, 
       y = NULL) +
  theme_classic() +
  theme(axis.text = element_text(colour = 'black'),
        strip.background = element_blank(),
        panel.background = element_rect(colour = 'black', size = 1))
ggsave(str_c(SAVE_LOCATION, '/Habitat_Classification/', save_suffix, '_site_eval.png', sep = ''),
       plot = site_plot2, height = 15, width = 15)


#### Write out all fit models ####
#Need to choose which model to use based on graph/table
print(str_c('Starting to fit all models to the final dataset', Sys.time(), sep = ': '))
top_models <- tuned_model_results %>%
  mutate(save_file_name = str_c(SAVE_LOCATION, '/Habitat_Classification/', save_suffix, "_", model_name, '_completeModel.rds', sep = '')) %>%
  mutate(model_fit = map2(.workflow, save_file_name, ~fit(.x, training(classification_data)) %>%
                            butcher(verbose = TRUE) %>%
                            write_rds(path = .y, compress = 'xz')))

parallel::stopCluster(cl)
plan(sequential)