library(here)
library(tidyverse)
library(brms)

setwd(here::here("data"))

p <- read_csv("network_antag_proportions_v02.csv") %>% 
  group_by(antag, res) %>% 
  mutate(index_id = cur_group_id()) %>% 
  arrange(index_id) %>% 
  mutate(x = as.numeric(scale(ghm))) 

index_id_key <- p %>% 
  dplyr::select(antag, res, index_id) %>% 
  distinct(.)

setwd(here::here("results/proportion_model_results"))

for(i in 1:max(p$index_id)){
  
  d <- p %>%
    filter(index_id == i)
  
  antag_lab <- unique(d$antag)
  res_lab <- unique(d$res)
  
  mod <- brm(
    formula = bf(
      prop ~ x + season + x:season,
      zoi ~ size, 
      coi ~ size),
    data = d,
    family = zero_one_inflated_beta(),
    prior = c(
      prior(normal(0, 1), class = Intercept),
      prior(normal(0, 0.5), class = b),
      prior(exponential(1), class = phi)),
    chains = 3,
    cores = 3,
    file = paste0("antag", antag_lab, "_proportions_", res_lab, "minutes_v01"))
  
  print(paste("Finished", i, "of", max(p$index_id)))
}
