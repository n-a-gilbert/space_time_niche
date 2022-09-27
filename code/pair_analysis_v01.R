library(brms)
library(here)
library(tidyverse)

setwd(here::here("data"))

load("supplement_pairwise_data_v01.RData")

# VARIABLES
### pair: the species pair (74 unique pairs)
### antag: antagonism classification (1 = low, 2 = medium, 3 = high)
### season: winter / summer
### period: which 3-month data collection period (8 total)
### tte: the time between successive detections (minutes) of the pair
### ghm: composite human disturbance variable, mean value within 5 km of camera
### ghm5k_sd: standard deviation of composite disturbance variable within 5 km of camera
### imperv: impervious cover (%) within 5 km of camera
### crop: cropland (%) within 5 km of camera
### ntot: total number of photos of the two species of the pair for a given period
### active: number of days the camera was active for a given period

glimpse(d)

# standardize predictors
d <- d %>% 
  mutate(
    ghm = as.numeric(scale(ghm)),
    ghm5k_sd = as.numeric(scale(ghm5k_sd)),
    crop = as.numeric(scale(crop)), 
    imperv = as.numeric(scale(imperv)),
    rai = as.numeric( scale ( log1p( ntot / active )))) # relative abundance index

# PRIMARY MODEL
setwd(here::here("results"))
m1 <- brm(
  data    = d,
  family  = Gamma(link = "log"),
  formula = tte ~ 1 + ghm + mo(antag) + ghm:mo(antag) + ( 1 + ghm |  pair / season),
  prior   = c(prior(normal(7.5, 2), class = Intercept),
              prior(normal(0, 0.5), class = b),
              prior(exponential(1), class = sd),
              prior(lkj(2), class = cor),
              prior(exponential(1), class = shape)),
  iter    = 4000,
  warmup  = 2000,
  chains  = 3,
  cores   = 3,
  thin    = 2,
  file = "primary_model_v01")

# HETEROGENEITY MODEL
m2 <- brm(
  data    = d,
  family  = Gamma(link = "log"),
  formula = tte ~ 1 + ghm + ghm5k_sd + mo(antag) + ghm:mo(antag) + ghm5k_sd:mo(antag) + ghm:ghm5k_sd:mo(antag) + ghm:ghm5k_sd + ( 1 + ghm + ghm5k_sd |  pair / season),
  prior   = c(prior(normal(7.5, 2), class = Intercept),
              prior(normal(0, 0.5), class = b),
              prior(exponential(1), class = sd),
              prior(lkj(2), class = cor),
              prior(exponential(1), class = shape)),
  iter    = 4000,
  warmup  = 2000,
  chains  = 3,
  cores   = 3,
  thin    = 2,
  file = "heterogeneity_model_v01")

# CROPLAND MODEL
m3 <- brm(
  data    = d,
  family  = Gamma(link = "log"),
  formula = tte ~ 1 + crop + mo(antag) + crop:mo(antag) + ( 1 + crop |  pair / season),
  prior   = c(prior(normal(7.5, 2), class = Intercept),
              prior(normal(0, 0.5), class = b),
              prior(exponential(1), class = sd),
              prior(lkj(2), class = cor),
              prior(exponential(1), class = shape)),
  iter    = 4000,
  warmup  = 2000,
  chains  = 3,
  cores   = 3,
  thin    = 2,
  file = "cropland_model_v01")

# IMPERVIOUS MODEL
m4 <- brm(
  data    = d,
  family  = Gamma(link = "log"),
  formula = tte ~ 1 + imperv + mo(antag) + imperv:mo(antag) + ( 1 + imperv |  pair / season),
  prior   = c(prior(normal(7.5, 2), class = Intercept),
              prior(normal(0, 0.5), class = b),
              prior(exponential(1), class = sd),
              prior(lkj(2), class = cor),
              prior(exponential(1), class = shape)),
  iter    = 4000,
  warmup  = 2000,
  chains  = 3,
  cores   = 3,
  thin    = 2,
  file = "imperv_model_v01")

# ABUNDANCE MODEL
m5 <- brm(
  data    = d,
  family  = Gamma(link = "log"),
  formula = tte ~ 1 + ghm + rai + mo(antag) + ghm:mo(antag) + rai:mo(antag) + ghm:rai:mo(antag) + ghm:rai + ( 1 + ghm |  pair / season),
  prior   = c(prior(normal(7.5, 2), class = Intercept),
              prior(normal(0, 0.5), class = b),
              prior(exponential(1), class = sd),
              prior(lkj(2), class = cor),
              prior(exponential(1), class = shape)),
  iter    = 4000,
  warmup  = 2000,
  chains  = 3,
  cores   = 3,
  thin    = 2,
  file = "abundance_model_v01")