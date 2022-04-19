library(terra)
library(tidyverse)


all_rasters <- list.files(path = '~/Coryphopterus/Habitat Association (Paper Y - PhD Ch. Y)/Results/INLA_models/model_Overall.9.7.21',
           pattern = 'tif$', full.names = TRUE) %>%
  tibble(file = .) %>%
  mutate(site = str_extract(file, 'BZ17-[0-9ABSNK]+')) %>%
  filter(str_detect(file, 'unsmoothed')) %>%
  rowwise %>%
  mutate(preds = list(rast(file))) %>%
  select(-file) %>%
  mutate(mean_pred = list(preds$mean),
         cv_pred = list(preds$cv),
         mean_log = list(log(mean_pred)))

all_rasters$preds[[1]]$mean

plot((log(all_rasters$preds[[1]]$q0.75) - log(all_rasters$preds[[1]]$q0.25)) / log(all_rasters$preds[[1]]$mean))

tmp <- (log(all_rasters$preds[[1]]$q0.75) - log(all_rasters$preds[[1]]$q0.25)) / log(all_rasters$preds[[1]]$mean) 

plot(tmp > 2000)

r <- tmp 
minimum <- -10
maximum <- 10

raster_filter <- function(r, minimum, maximum){
  intermediate <- ifel(r <= maximum & r >= minimum, r, NA)
  
  mask(r, intermediate)
}


plot(raster_filter(tmp, -1, 0.1))

tmo2 <- tmp %>%
  values() %>%
  as.numeric()


tmp3 <- tmo2[!is.nan(tmo2)]

length(tmp3[tmp3 > 0])


plot(log(tmp))

plot(all_rasters$preds[[1]])

plot(log(all_rasters$preds[[1]]$sd) / log(all_rasters$preds[[1]]$mean))

all_rasters$preds[[1]]$cv

all_rasters$cv_pred[[1]] %>% plot

vector_qualities_mean <- all_rasters %>%
  mutate(mean_vec = list(values(mean_log))) %>%
  select(mean_vec) %>%
  unnest(mean_vec) %>%
  filter(!is.nan(mean_vec)) %>%
  pull(mean_vec) %>%
  as.numeric()

vector_cv <- all_rasters %>%
  mutate(mean_vec = list(values(cv_pred))) %>%
  select(mean_vec) %>%
  unnest(mean_vec) %>%
  filter(!is.nan(mean_vec)) %>%
  pull(mean_vec) %>%
  as.numeric()


tmp <- all_rasters %>%
  mutate(quality = list((mean_log - min(vector_qualities_mean)) / (max(vector_qualities_mean) - min(vector_qualities_mean))),
         quality_cv = list((cv_pred - min(vector_cv)) / (max(vector_cv) - min(vector_cv)))) %>%
  select(site, quality, quality_cv)

for(i in 1:nrow(tmp)){
  writeRaster(tmp$quality[[i]], str_c('~/Coryphopterus/Dispersal/Data/', tmp$site[i], '_quality.tif'))
}



tmp$quality[[2]] %>% plot
tmp$quality_cv[[2]] %>% plot

r <- rast('~/Coryphopterus/Dispersal/Data/OLD/BZ17-0A_quality.tif')
plot(r)
plot(tmp$quality[[1]])
