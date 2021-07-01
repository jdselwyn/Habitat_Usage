rm(list = ls())

library(sf)
library(tidyverse)
library(magrittr)
library(readxl)
library(janitor)
library(broom)
select <- dplyr::select

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

#### Which cluster goes to which species ####
individual_data %>%
  filter(!is.na(blast_species_hit), !is.na(dapc_assignment)) %>%
  count(blast_species_hit, dapc_assignment)

individual_data <- individual_data %>%
  mutate(dapc_assignment = if_else(dapc_assignment == 'Cluster 1', 'Coryphopterus hyalinus', 'Coryphopterus personatus')) 

#### Size difference btw species ####
individual_data %>%
  filter(!is.na(dapc_assignment)) %>%
  select(ID, dapc_assignment, sl_mm, tl_mm) %>%
  pivot_longer(cols = ends_with('mm'),
               names_to = 'measure',
               values_to = 'length') %>%
  group_by(measure) %>%
  summarise(            n_chya = sum(dapc_assignment == 'Coryphopterus hyalinus'),
                        n_cpers = sum(dapc_assignment == 'Coryphopterus personatus'),
                        t_test = list(t.test(length ~ dapc_assignment)),
                        .groups = 'drop') %>%
  rowwise %>%
  mutate(t_test = list(tidy(t_test))) %>%
  unnest(t_test) %>%
  rename_with(~str_replace_all(., c('1$' = '_chya', '2$' = '_cpers')))


#### NEED TO SORT OUT WHY THESE DISAGREE!!!! ####
## Sorted out - some typos and some that were collected but not sequenced
shoal_data$Nut %>% unique
shoal_data$Site %>% unique

shoal_data %>%
  filter(!is.na(Nut), Y.N == 'Y') %>%
  group_by(Site, Nut) %>%
  summarise(n = n()) %>%
  summarise(n = n()) %>%
  mutate(Site = as.character(Site)) %>%
  arrange(Site)

individual_data %>%
  filter(year == 2017) %>%
  group_by(site, cloud) %>%
  summarise(n = n()) %>%
  summarise(n = n()) %>%
  mutate(site = as.character(site)) %>%
  arrange(site)

individual_data %>%
  filter(year == 2017) %>%
  distinct(site, cloud)


full_join(
  
  shoal_data %>%
    filter(Y.N == 'Y') %>%
    distinct(Site, Nut) %>%
    mutate(source = 'shoal'),
  
  individual_data %>%
    filter(year == 2017) %>%
    distinct(site, cloud) %>%
    mutate(source = 'individual'),
  
  by = c('Site' = 'site', 'Nut' = 'cloud')
  
) %>% 
  arrange(Site) %>%

  write_csv('tmp.csv')

#### Join shoals and genetics ####
a<-1;b<-1
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
  mutate(prop_cpers = (cpers+a)/(total+a+b),
         cpers_lwr.95 = qbeta(0.025, cpers+a, total-cpers+b),
         cpers_lwr.50 = qbeta(0.25, cpers+a, total-cpers+b),
         cpers_upr.50 = qbeta(0.75, cpers+a, total-cpers+b),
         cpers_upr.95 = qbeta(0.975, cpers+a, total-cpers+b)) %>%
  left_join(dplyr::select(shoal_data, Site, Nut, geometry), by = c("Site", "Nut")) %>%
  st_as_sf

shoal_composition %>%
  mutate(prop_cpers = cpers / total) %>%
  select(Site, Nut, geometry, Shoal.Size, prop_cpers) %>%
  mutate(prop_chya = 1 - prop_cpers) %>%
  st_write('../Results/Shoal Composition/shoal_composition.shp', 
           delete_dsn = file.exists('../Results/Shoal Composition/shoal_composition.shp'))


percent_cpers_overall <- shoal_composition %>%
  as_tibble %>%
  group_by(Site) %>%
  summarise(n_shoal = n(), 
            across(c(cpers, chya, total), sum),
            .groups = 'drop') %>%
  mutate(site_cpers = (cpers+a)/(total+a+b),
         site_cpers_lwr.95 = qbeta(0.025, cpers+a, total-cpers+b),
         site_cpers_lwr.50 = qbeta(0.25, cpers+a, total-cpers+b),
         site_cpers_upr.50 = qbeta(0.75, cpers+a, total-cpers+b),
         site_cpers_upr.95 = qbeta(0.975, cpers+a, total-cpers+b)) %>%
  mutate(Site = factor(Site, levels = str_c('BZ17', c('10KN', '5KN', '1KN', '500N', '100N', 
                                                      '0A', '0B', '60S', '100S', '500S', '1KS', '5KS'),
                                            sep = '-'))) 

percent_cpers_overall %>%
  select(Site, site_cpers, site_cpers_lwr.95, site_cpers_upr.95) %>%

  arrange(Site) %>%
  mutate(chya_perc = 1 - site_cpers,
         chya_lwr = 1 - site_cpers_upr.95,
         chya_upr = 1 - site_cpers_lwr.95) %>%
  mutate(across(where(is.numeric), ~.*100),
         across(where(is.numeric), ~round(., 1) %>% as.character)) %>%
  mutate(cpers = str_c(site_cpers, ' (', site_cpers_lwr.95, ' - ', site_cpers_upr.95, ')'),
         chya = str_c(chya_perc, ' (', chya_lwr, ' - ', chya_upr, ')')) %>%
  select(Site, chya, cpers) %>%
  write_csv('tmp_table.csv')


percent_cpers_overall %>%
  summarise(across(chya:total, sum)) %>%
  mutate(site_chya = (chya+a)/(total+a+b),
         site_chya_lwr.95 = qbeta(0.025, chya+a, total-chya+b),
         site_chya_upr.95 = qbeta(0.975, chya+a, total-chya+b))


a_shoal <- 1; b_shoal<-1; a_site <- 1; b_site<-1 # 
difference_from_site <- shoal_composition %>%
  dplyr::select(Site:total) %>%
  left_join(percent_cpers_overall %>%
              dplyr::select(Site, cpers, total) %>%
              dplyr::rename(site_cpers = cpers, site_total = total)) %>%
  mutate(mean_diff_site = (a_shoal+cpers)/(a_shoal+b_shoal+total)-(a_site+site_cpers)/(a_site+b_site+site_total),
         sd_diff_site = sqrt(((a_shoal+cpers)*(b_shoal+total-cpers))/((a_shoal+b_shoal+total)^2*(a_shoal+b_shoal+total+1)) + ((a_site+site_cpers)*(b_site+site_total-site_cpers))/((a_site+b_site+site_total)^2*(a_site+b_site+site_total+1))),
         lwr_95_site_diff = qnorm(0.025,mean_diff_site,sd_diff_site),
         upr_95_site_diff = qnorm(0.975,mean_diff_site,sd_diff_site),
         different = case_when(lwr_95_site_diff < 0 & upr_95_site_diff < 0 ~ TRUE,
                               lwr_95_site_diff > 0 & upr_95_site_diff > 0 ~ TRUE,
                               TRUE ~ FALSE))

sum(!difference_from_site$different)/nrow(difference_from_site)

library(tidytext)
shoal_composition %>%
  as_tibble %>%
  left_join(difference_from_site %>%
              select(Site, Nut, different),
            by = c('Site', 'Nut')) %>%
  mutate(Nut = reorder_within(Nut, prop_cpers, Site)) %>%
  ggplot(aes(x = prop_cpers, y = Nut, xmin = cpers_lwr.95, xmax = cpers_upr.95, colour = different)) +
  geom_vline(data = percent_cpers_overall, aes(xintercept = site_cpers_lwr.95), colour = 'black') +
  geom_vline(data = percent_cpers_overall, aes(xintercept = site_cpers_upr.95), colour = 'black') +
  geom_linerange(linetype = 'dotted') +
  geom_linerange(aes(xmin = cpers_lwr.50, xmax = cpers_upr.50)) +
  geom_point(aes(size = total)) +
  scale_y_reordered() +
  facet_wrap(~Site, scales = 'free_y')

shoal_differences <- shoal_composition %>%
  as_tibble %>%
  left_join(difference_from_site %>%
              select(Site, Nut, different),
            by = c('Site', 'Nut')) %>%
  mutate(Nut = reorder_within(Nut, prop_cpers, Site)) %>%
  mutate(Site = factor(Site, levels = str_c('BZ17', c('10KN', '5KN', '1KN', '500N', '100N', 
                                                      '0A', '0B', '60S', '100S', '500S', '1KS', '5KS'),
                                            sep = '-'))) %>%
  mutate(different = if_else(different, 'Different from Background Reef', NA_character_))


percent_cpers_overall <- mutate(percent_cpers_overall,
                                Site = LETTERS[as.integer(Site)])

shoal_differences <- mutate(shoal_differences,
                            Site = LETTERS[as.integer(Site)])


shoal_differences %>%
  ggplot(aes(y = Site)) +
  geom_linerange(data = percent_cpers_overall, aes(xmin = site_cpers_lwr.95, xmax = site_cpers_upr.95), size = 10, alpha = 0.5) +
  # geom_rect(data = percent_cpers_overall, aes(xmin = Site-0.5, xmax = Site+0.5, ymin = site_cpers_lwr.95, ymax = site_cpers_upr.95)) +
  geom_linerange(aes(x = prop_cpers, xmin = cpers_lwr.95, 
                     xmax = cpers_upr.95, group = Nut, colour = different), 
                  position = position_dodge(0.5), show.legend = FALSE) +
  geom_point(aes(x = prop_cpers, group = Nut, colour = different), 
                 position = position_dodge(0.5), show.legend = TRUE) +
  theme_classic() +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_y_discrete(limits = rev) +
  scale_color_manual(values = c('Different from Background Reef' = 'red'), 
                     na.value = 'black',
                     breaks = 'Different from Background Reef') +
  labs(x = 'Percent C. personatus',
       y = NULL,
       colour = NULL) +
  theme(axis.text = element_text(colour = 'black'),
        legend.position = 'top')
 # ggsave('../Manuscript/Figures/shoal_composition.png', width = 7, height = 6)
ggsave('../Manuscript/Figures/Figure S1.png', width = 7, height = 6)

shoal_differences %>%
  filter(!is.na(different)) %>%
  filter(total == 20) %>% pull(Nut) %>% as.character()

