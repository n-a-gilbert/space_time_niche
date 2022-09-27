library(tidyverse)
library(here)
library(brms)

setwd(here::here("data"))

# id = site x season id key
# ghm = mean of composite human disturbance metric within 5 km of camera
# connectance - proportion of species that co-occur our of total possible
# size - the number of species (nodes) in the network
# res - temporal resolution used to define network (i.e. time threshold between detections to be considered co-occuring)
d <- read_csv("network_connectance_v01.csv") %>% 
  mutate(x = as.numeric(scale(ghm)))

# zero-one inflated models are appropriate because connectance can take on values of 0 and 1
# we model these components by network size, since connectances of 1 or 0 are
# much more likely with smaller networks
d %>% 
  mutate(zero = ifelse(connectance == 0, "Zero", 
                       ifelse(connectance == 1, "One", "Mid"))) %>% 
  mutate(zero = factor(zero, 
                       levels = c(
                         "One",
                         "Mid", 
                         "Zero"))) %>% 
  mutate(res = factor(res, 
                      levels = c(
                        "10min",
                        "1hr", 
                        "1day", 
                        "1week", 
                        "1month"))) %>% 
  ggplot(aes(x = factor(size), fill = zero)) +
  geom_bar() + 
  facet_wrap(~res)

# run the model separately for the different temporal resolutions
setwd(here::here("results"))

m_10 <-
  brm(
    formula = bf(
      connectance ~ x + season + x:season, 
      zoi ~ size, 
      coi ~ size),
    data = filter(d, res == "10min"),
    family = zero_one_inflated_beta(),
    prior = c(
      prior(normal(0, 1), class = Intercept),
      prior(normal(0, 0.5), class = b),
      prior(exponential(1), class = phi)),
    chains = 3,
    cores = 3,
    file = "connectance_10minutes_v01")

m_60 <-
  brm(
    formula = bf(
      connectance ~ x + season + x:season, 
      zoi ~ size, 
      coi ~ size),
    data = filter(d, res == "1hr"),
    family = zero_one_inflated_beta(),
    prior = c(
      prior(normal(0, 1), class = Intercept),
      prior(normal(0, 0.5), class = b),
      prior(exponential(1), class = phi)),
    chains = 3,
    cores = 3,
    file = "connectance_60minutes_v01")

m_1440 <-
  brm(
    formula = bf(
      connectance ~ x + season + x:season, 
      zoi ~ size, 
      coi ~ size),
    data = filter(d, res == "1day"),
    family = zero_one_inflated_beta(),
    prior = c(
      prior(normal(0, 1), class = Intercept),
      prior(normal(0, 0.5), class = b),
      prior(exponential(1), class = phi)),
    chains = 3,
    cores = 3,
    file = "connectance_1440minutes_v01")

m_10080 <-
  brm(
    formula = bf(
      connectance ~ x + season + x:season, 
      zoi ~ size, 
      coi ~ size),
    data = filter(d, res == "1week"),
    family = zero_one_inflated_beta(),
    prior = c(
      prior(normal(0, 1), class = Intercept),
      prior(normal(0, 0.5), class = b),
      prior(exponential(1), class = phi)),
    chains = 3,
    cores = 3,
    file = "connectance_10080minutes_v01")

m_43200 <- 
  brm(
    formula = bf(
      connectance ~ x + season + x:season, 
      zoi ~ size, 
      coi ~ size),
    data = filter(d, res == "1month"),
    family = zero_one_inflated_beta(),
    prior = c(
      prior(normal(0, 1), class = Intercept),
      prior(normal(0, 0.5), class = b),
      prior(exponential(1), class = phi)),
    chains = 3,
    cores = 3,
    file = "connectance_43200minutes_v01")