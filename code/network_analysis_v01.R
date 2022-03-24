library(readr)
library(dplyr)
library(tibble)
library(here)
library(network)
library(tidyr)
library(stringr)
library(brms)

setwd(here::here("data"))

d <- readr::read_csv("space_time_pair_ttd_data_v01.csv")

# create an ID column to indicate unique site - sampling period combinations
# separate the pair column to create network
# clean up some species names
d <- d %>% 
  dplyr::group_by(site, period) %>% 
  dplyr::mutate(id = dplyr::cur_group_id()) %>% 
  tidyr::separate(pair, into = c("sp1", "sp2")) %>% 
  dplyr::mutate(sp1 = stringr::str_replace(sp1, "SQUIRRELSANDCHIPMUNKS", "SQUIRRELS"),
                sp1 = stringr::str_replace(sp1, "FOXRED", "RED FOX"), 
                sp1 = stringr::str_replace(sp1, "SKUNKSTRIPED", "SKUNK"),
                sp1 = stringr::str_replace(sp1, "SNOWSHOEHARE", "HARE"),
                sp2 = stringr::str_replace(sp2, "SQUIRRELSANDCHIPMUNKS", "SQUIRRELS"),
                sp2 = stringr::str_replace(sp2, "FOXRED", "RED FOX"), 
                sp2 = stringr::str_replace(sp2, "SKUNKSTRIPED", "SKUNK"),
                sp2 = stringr::str_replace(sp2, "SNOWSHOEHARE", "HARE"))  

den <- list()
prop <- list()

for(i in 1:length(unique(d$id))){
  
  # grab data for each site - period combo
  filt <- d %>%
    filter(id == i) %>% 
    dplyr::select(id, sp1, sp2, antag) %>% 
    distinct(.)
  
  sp1 <- tolower(filt$sp1)
  sp2 <- tolower(filt$sp2)
  netmat <- cbind(sp1, sp2)
  
  # create network
  net <- network::network(netmat, matrix.type = "edgelist", directed = FALSE)
  set.edge.attribute(net, "antag", filt$antag)
  
  # store results in lists
  # first, density of network (propotion of total observed links)
  den[[i]] <- network.density(net)
  # now, proportion of links in each network for each antag ranking
  prop[[i]] <- as.data.frame(table(get.edge.attribute(net, attrname = 'antag')) /
                               length(get.edge.attribute(net, attrname = 'antag'))) %>% 
    setNames(., c("antag", "prop")) %>% 
    add_column(id = unique(filt$id))
  
}

density_df <-
  cbind(id = 1:length(unique(d$id)),
                    den = do.call(rbind, den)) %>%
  tibble::as_tibble() %>% 
  dplyr::rename(density = V2) %>% 
  dplyr::full_join(d) %>% 
  dplyr::select(id, x, x2, season, density) %>% 
  dplyr::distinct(.)

prop_df <-
  expand.grid(id = 1:length(unique(d$id)),
            antag = 1:3) %>% 
  dplyr::mutate(antag = as.factor(antag)) %>%
  dplyr::full_join(do.call(rbind, prop)) %>% 
  dplyr::mutate(prop = ifelse(is.na(prop), 0, prop),
                antag = as.numeric(antag)) %>% 
  dplyr::full_join(dplyr::select(d, -antag)) %>% 
  dplyr::select(id, x, x2, season, antag, prop) %>% 
  dplyr::distinct(.) %>% 
  tibble::as_tibble(.)

# Save files if you want (used for plotting)
# setwd(here::here("data"))
# write_csv(density_df, "space_time_network_density_v01.csv")
# write_csv(prop_df, "space_time_network_proportions_v01.csv")

density_df <- read_csv("space_time_network_density_v01.csv")
prop_df <- read_csv("space_time_network_proportions_v01.csv")

# network density analysis
# ART a few minutes
density_model <-  
  brms::brm(density ~ x + season + season:x,
      data = density_df,
      family = zero_one_inflated_beta(),
      prior = c(
        prior(normal(0, 1), class = Intercept),
        prior(normal(0, 0.5), class = b),
        prior(exponential(1), class = phi)),
      chains = 3,
      cores = 3)

# network proportion models
# each takes a few minutes to run
prop1_model <- 
  brms::brm(prop ~ x + season + x:season,
      data = filter(prop_df, antag == 1),
      family = zero_one_inflated_beta(),
      prior = c(
        prior(normal(0, 1), class = Intercept),
        prior(normal(0, 0.5), class = b),
        prior(exponential(1), class = phi)),
      chains = 3,
      cores = 3)

prop2_model <- 
  brms::brm(prop ~ x + season + x:season,
            data = filter(prop_df, antag == 2),
            family = zero_one_inflated_beta(),
            prior = c(
              prior(normal(0, 1), class = Intercept),
              prior(normal(0, 0.5), class = b),
              prior(exponential(1), class = phi)),
            chains = 3,
            cores = 3)

prop3_model <- 
  brms::brm(prop ~ x + season + x:season,
            data = filter(prop_df, antag == 3),
            family = zero_one_inflated_beta(),
            prior = c(
              prior(normal(0, 1), class = Intercept),
              prior(normal(0, 0.5), class = b),
              prior(exponential(1), class = phi)),
            chains = 3,
            cores = 3)

setwd(here::here("results"))
save(density_model, 
     prop1_model, 
     prop2_model, 
     prop3_model, 
     file = "network_analysis_results_v01.RData")