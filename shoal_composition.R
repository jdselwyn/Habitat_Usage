rm(list = ls())

library(sf)
library(terra)
library(raster)
library(tidyverse)
library(magrittr)
library(readxl)
library(janitor)
library(brms)


#### Read in Data ####
shoal_data <- list.files("~/Coryphopterus/Maps/COPE_Sites/Shoals", pattern = '.shp$', recursive = TRUE, full.names = TRUE) %>%
  tibble(file = .) %>%
  mutate(Site = str_extract(file, 'BZ17-[0-9ABNSK]+'),
         Site = factor(Site, levels = str_c('BZ17', c('10KN', '5KN', '1KN', '500N', '100N', 
                                                      '0A', '0B', '60S', '100S', '500S', '1KS', '5KS'),
                                            sep = '-'))) %>%
  mutate(shoals = map(file, st_read, quiet = TRUE)) %>%
  select(-file) %>%
  unnest(shoals) %>%
  filter(Shoal.Size > 0) %>%
  mutate(Nut = map_chr(Nut, ~str_split(.x, '-', simplify = TRUE) %>% sort %>% str_to_lower %>% str_c(collapse = '-')),
         Nut = str_replace(Nut, 'orage', 'orange'))

shoal_data %>%
  filter(Y.N == 'Y') %>%
  group_by(Site) %>%
  summarise(n = n())


individual_data <- unlist(str_split('cdcccccccnnnciiiccclcccccicccc','')) %>%
  tibble(col_type = .) %>%
  mutate(col_type = case_when(col_type == 'c' ~ 'text',
                              col_type == 'd' ~ 'numeric',
                              col_type == 'n' ~ 'numeric',
                              col_type == 'i' ~ 'numeric',
                              col_type == 'l' ~ 'text')) %$%
  read_excel('~/Coryphopterus/Master COPE collection database.xlsx',sheet = 1, 
             col_types = col_type, na = c('N/A','NA', "?")) %>% 
  clean_names %>%
  dplyr::select(tube_label, year, site, cloud, sl_mm, tl_mm) %>%
  dplyr::rename(ID = tube_label) %>%
  left_join(read_csv("~/Coryphopterus/Bioinformatics/Mitochondrial Mapping/mitochondrial_blast_results.csv"), by = 'ID') %>%
  left_join(read_csv('~/Coryphopterus/Bioinformatics/Species_Split/Results/DAPC_assignments.csv'), by = 'ID') %>%
  left_join(read_csv("~/Coryphopterus/Bioinformatics/Species_Split/Results/Structure_assignments.csv"), by = 'ID') %>%
  mutate(cloud = map_chr(cloud, ~str_split(.x, '-', simplify = TRUE) %>% sort %>% str_to_lower %>% str_c(collapse = '-')),
         cloud = str_remove(cloud, '^[ab]-'))

pca_model_loadings <- read_csv('../Results/Topography/complexity_pca_loadings.csv')

topography <- list.files('../Intermediate_Files/Topography', pattern = 'tif$', full.names = TRUE) %>%
  tibble(habitat = .) %>%
  mutate(Site = str_extract(habitat, 'BZ17-[0-9ABKNS]+'))

#### Which cluster goes to which species ####
individual_data %>%
  filter(!is.na(blast_species_hit), !is.na(dapc_assignment)) %>%
  count(blast_species_hit, dapc_assignment)

individual_data <- individual_data %>%
  mutate(dapc_assignment = if_else(dapc_assignment == 'Cluster 1', 'Coryphopterus hyalinus', 'Coryphopterus personatus')) 


#### Get habitat data for shoal centroids ####
shoal_composition <- inner_join(shoal_data, individual_data, by = c('Site' = 'site', 'Nut' = 'cloud')) %>%
  filter(!is.na(dapc_assignment)) %>%
  # mutate(dapc_assignment = if_else(dapc_assignment == 'Cluster 1', 'chya', 'cpers')) %>%
  mutate(dapc_assignment = str_c('c', str_extract(dapc_assignment, 'hya|pers'), sep = '')) %>%
  group_by(Site, Nut, Shoal.Size, dapc_assignment) %>%
  summarise(n = n(), .groups = 'drop') %>%
  pivot_wider(names_from = 'dapc_assignment',
              values_from = 'n', 
              values_fill = 0L) %>%
  mutate(total = chya + cpers) %>%
  left_join(dplyr::select(shoal_data, Site, Nut, geometry), by = c("Site", "Nut")) %>%
  st_as_sf %>%
  st_transform(crs = '+proj=utm +zone=16 +datum=WGS84 +units=m +no_defs') %>%
  
  nest(shoals = -c(Site)) %>%
  # inner_join(global_model_results, by = 'Site') %>%
  inner_join(topography, by = 'Site') %>%
  
  mutate(across(c(habitat), ~map2(., shoals, ~stack(.x) %>% extract(as(.y, 'Spatial'))))) %>%
  mutate(across(c(habitat), ~map(., as_tibble))) %>%
  unnest(c(shoals, habitat)) %>%
  dplyr::select(Site, Nut, Shoal.Size, chya, cpers, total, 
                depth, boundaryDistance, moran, slope, relief, rugosity, 
                vectorDispersion, viewshed, classified_habitat_unsmoothed) %>%
  mutate(to_transform = cbind(slope, relief, rugosity, vectorDispersion),
         coarse_complexity = -1 * moran,
         depth = -1 * depth,
         habitat = if_else(classified_habitat_unsmoothed == 2, 'Reef', 'Sand')) %>%
  dplyr::select(-slope, -relief, -rugosity, -vectorDispersion, -moran, -classified_habitat_unsmoothed) %>%
  mutate(fine_complexity = to_transform %*% matrix(pca_model_loadings$Comp.1, nrow = 4, ncol = 1) %>% as.numeric) %>%
  dplyr::select(-to_transform) %>%
  mutate(shoal = str_c(Site, Nut, sep = '_') %>% factor,
         Site = factor(Site),
         habitat = factor(habitat),
         across(c(Shoal.Size, depth, boundaryDistance, viewshed, coarse_complexity), ~scale(.) %>% as.numeric)) %>%
  dplyr::select(Site, shoal, habitat, chya, cpers, total, Shoal.Size, everything(), -Nut) 


#### Shoal Composition Model ####
get_prior(cpers | trials(total) ~ Shoal.Size + habitat + depth + boundaryDistance + viewshed + coarse_complexity + fine_complexity + 
            (1 | Site) + (1 | shoal),
          family = 'binomial',
          data = shoal_composition)


x <- seq(0, 10, length.out = 100)

plot(dgamma(x, shape = 1, rate = 1) ~ x, type = 'b', col = 'black')
lines(dgamma(x, shape = 2, rate = 1) ~ x, type = 'b', col = 'blue')
lines(dgamma(x, shape = 2, rate = 2) ~ x, type = 'b', col = 'red')
lines(dt(x, df = 3) ~ x, type = 'b', col = 'purple')

composition_model <- brm(cpers | trials(total) ~ habitat + Shoal.Size + depth + boundaryDistance + viewshed + coarse_complexity + fine_complexity + 
                           (1 | Site) + (1 | shoal),
                         family = 'binomial',
                         data = shoal_composition,
                         prior = c(prior(normal(0, 10), class = Intercept),
                                   prior(normal(0, 10), class = b),
                                   prior(gamma(2, 2), class = sd)),
                         
                         backend = 'cmdstanr',
                         cores = 4,
                         chains = 4,
                         iter = 2000,
                         warmup = 1000,
                         
                         sample_prior = 'yes',
                         algorithm = 'sampling',
                         file = '../Intermediate_Files/BRMS_models/shoal_composition_model',
                         file_refit = 'on_change')

#### Model Diagnostics ####
plot(composition_model) #Trace plots for MCMC
summary(composition_model, priors = TRUE) #Check Rhat and Bulk/Tail ESS
pp_check(composition_model) #Check observed data is similar to posterior samples

#### Model Interpretation ####
bayes_R2(composition_model)

parnames(composition_model)

bayesplot::mcmc_areas(as.matrix(composition_model),
                      regex_pars =  c('Shoal.Size', 'depth', 'boundaryDistance', 
                                      'coarse_complexity', 'fine_complexity', 'viewshed',
                                      'habitat')) +
  geom_vline(xintercept = 0, linetype = 'dashed')
bayesplot::mcmc_areas(as.matrix(composition_model), regex_pars = c('viewshed', 'habitat')) +
  geom_vline(xintercept = 0, linetype = 'dashed')


hypothesis(composition_model, c('Shoal.Size < 0', 'depth < 0', 
                                'boundaryDistance > 0', 'viewshed > 0', 
                                'coarse_complexity > 0', 'fine_complexity < 0',
                                'habitatSand < 0'), alpha = 0.05) 


