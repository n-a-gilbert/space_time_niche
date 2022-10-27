# data files to run this code cannot be shared, since
# geographic coordinates are necessary to check for spatial
# autocorrelation; these coordinates are not allowed to be
# shared to protect privacy of program volunteers

library(here)
library(tidyverse)
library(tidybayes)
library(brms)
library(moranfast)

# read in and prepare data used to fit the model
setwd(here::here("data"))

dat <- read_csv("tte_data_v04_with_abundance_camera_active.csv")

mass_antag <- read_csv("pair_antagonism_new_old_v02.csv") %>% 
  mutate(pair = paste(sp1, sp2, sep = "-")) %>% 
  dplyr::select(pair, antag = new_antag)

keeper_pairs <- dat %>% 
  group_by(xcoord, ycoord) %>% 
  mutate(id = cur_group_id()) %>% 
  filter(tte < 43800) %>% 
  mutate(pair = paste(sp1, sp2, sep = "-")) %>%
  ungroup() %>% 
  group_by(pair) %>%
  count(id) %>% 
  mutate(n = 1) %>% 
  summarise(nsite = sum(n)) %>% 
  filter(nsite > 49) %>% 
  pull(pair)

d <- dat %>%
  mutate(pair = paste(sp1, sp2, sep = "-"),
         y = tte + 1, 
         x = as.numeric(scale(ghm)), 
         x2 = x*x,
         rai = as.numeric( scale ( log1p( ntot / active )))) %>% 
  filter(tte < 43800) %>% 
  filter(pair %in% keeper_pairs) %>% 
  arrange(pair) %>% 
  left_join(mass_antag) %>% 
  dplyr::select(pair,
                xcoord, ycoord, season,
                y, x, x2,
                antag, rai)

# model output
setwd(here::here("results/revision"))
mod <- readRDS("ghm5k_antag_v01.rds")

# add posterior residuals to the data used to fit the model
residuals <- d %>%
  add_residual_draws(mod)

# going to evaluate spatial autocorrelation in residuals at
### 4 levels, starting from fine and going to broad

### 1) pair x season x site level;
#-#-#-#-# that is, for each location, calculate mean residual
#-#-#-#-# for each pair-season combination

### 2) pair x site level;
#-#-#-#-# that is, for each location, calculate mean residual
#-#-#-#-# for each pair

### 3) aggregated to antagonism class;
#-#-#-#-# that is, for each location, calculate mean residual
#-#-#-#-# for each antagonism class

### 4) fully aggregated;
#-#-#-#-# that is, for each location, calculate mean residual
#-#-#-#-# for across all pairs / seasons 

### 1) pair x season x site level;
rmean_pair_season <- residuals %>%
  group_by(pair, season, xcoord, ycoord) %>%
  summarise(residual = mean(.residual)) %>% 
  mutate(pair_season = paste(pair, season, sep = "_")) %>% 
  group_by(xcoord, ycoord) %>%
  mutate(site_id = cur_group_id()) %>%
  ungroup() %>%
  group_by(pair_season) %>%
  mutate(n_sites = n_distinct(site_id)) %>%
  filter(n_sites > 2)

pair <- unique(rmean_pair_season$pair_season)
pair_pval <- list()
for(i in 1:length(pair)) {
  pair_pval[i] <- moranfast(rmean_pair_season[rmean_pair_season$pair_season == pair[i], ]$residual,
                            rmean_pair_season[rmean_pair_season$pair_season == pair[i], ]$xcoord,
                            rmean_pair_season[rmean_pair_season$pair_season == pair[i], ]$ycoord)[[4]]
  
  print(paste("finished", i, "of", length(pair)))
}

# 5/140 (3.57%) of pair-season combinations with
# evidence of spatial autocorrelation in residuals
cbind(
  pair = pair,
  moranpval = unlist(pair_pval)) %>% 
  as_tibble() %>% 
  filter(moranpval < 0.05)

rmean_pair <- residuals %>% 
  group_by(pair, xcoord, ycoord) %>% 
  summarise(residual = mean(.residual)) %>% 
  group_by(xcoord, ycoord) %>%
  mutate(site_id = cur_group_id()) %>%
  ungroup() %>%
  group_by(pair) %>%
  mutate(n_sites = n_distinct(site_id)) %>%
  filter(n_sites > 2)

pair <- unique(rmean_pair$pair)
pair_pval <- list()
for(i in 1:length(pair)) {
  pair_pval[i] <- moranfast(rmean_pair[rmean_pair$pair == pair[i], ]$residual,
                            rmean_pair[rmean_pair$pair == pair[i], ]$xcoord,
                            rmean_pair[rmean_pair$pair == pair[i], ]$ycoord)[[4]]
  
  print(paste("finished", i, "of", length(pair)))
}

# 5 / 74 (6.76%) of pairs with evidence of 
# spatial autocorrelation in residuals
cbind(
  pair = pair,
  moranpval = unlist(pair_pval)) %>% 
  as_tibble() %>% 
  filter(moranpval < 0.05)

rmean_antag <- residuals %>% 
  group_by(antag, xcoord, ycoord) %>% 
  summarise(residual = mean(.residual))%>% 
  # mutate(pair_season = paste(pair, season, sep = "_")) %>% 
  group_by(xcoord, ycoord) %>%
  mutate(site_id = cur_group_id()) %>%
  ungroup() %>%
  group_by(antag) %>%
  mutate(n_sites = n_distinct(site_id)) %>%
  filter(n_sites > 2)

antag <- unique(rmean_antag$antag)
antag_pval <- list()
for(i in 1:length(antag)) {
  antag_pval[i] <- moranfast(rmean_antag[rmean_antag$antag == antag[i], ]$residual,
                            rmean_antag[rmean_antag$antag == antag[i], ]$xcoord,
                            rmean_antag[rmean_antag$antag == antag[i], ]$ycoord)[[4]]
  
  print(paste("finished", i, "of", length(antag)))
}

# no evidence of spatial autocorrelation for anatagonism levels
cbind(
  antag = antag,
  moranpval = unlist(antag_pval)) %>% 
  as_tibble() %>% 
  filter(moranpval < 0.05)

rmean <- residuals %>% 
  group_by(xcoord, ycoord) %>% 
  summarise(residual = mean(.residual)) 

# no evidence of spatial autocorrelation at global level
moranfast(rmean$residual, 
          rmean$xcoord, 
          rmean$ycoord)