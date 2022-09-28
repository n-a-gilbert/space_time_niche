# code to create co-occurrence networks for each camera site
# and calculate network connectance and proportion of each antagonism classfication

library(tidyverse)
library(here)
library(network)

setwd(here::here("data"))

# site_id - id for site 
# id - site x period id
# sp1 - species 1
# sp2 - species 2
# antag - antagonism (1 = low, 2 = medium, 3 = high)
# ghm - human disturbance composite metric (mean within 5 km of cam)
# season - summer/winter
# tte - time between detections of the 2 species
d <- read_csv("data_to_create_networks_v01.csv")

res <-  c(10, 60, 1440, 10080, 43200)

d2 <- tibble(
  site_id = rep(d$site_id, length(res)),
  period = rep(d$period, length(res)), 
  season = rep(d$season, length(res)),
  sp1 = rep(d$sp1, length(res)),
  sp2 = rep(d$sp2, length(res)), 
  tte = rep(d$tte, length(res)), 
  antag = rep(d$antag, length(res)),
  res = rep(res, each = nrow(d))) %>% 
  group_by(site_id, period, res) %>% 
  mutate(id = cur_group_id())

prop <- list(list())
size <- list(list())
connectance <- list(list())

for(i in 1:max(d2$id)){
  
  tot <- d2 %>% 
    ungroup(.) %>% 
    filter(id == i)
  
  filt <- tot %>% 
    filter(tte < res) %>% 
    group_by(sp1, sp2, antag) %>% 
    summarise(weight = mean(tte)) %>% 
    filter(!is.na(antag))
  
  if(nrow(filt) == 0){
    
    prop[[i]] <- tibble(
      antag = factor(c(1, 2, 3)), 
      prop = c(0, 0, 0), 
      id = i)
    
    size[[i]] <- tibble(
      size = 0, 
      id = i)
    
    connectance[i] <- 0
    
  } else {
    
    sp1 <- tolower(filt$sp1)
    sp2 <- tolower(filt$sp2)
    
    totsp <- length( unique(c(tot$sp1, tot$sp2)))
    
    netmat <- cbind(sp1, sp2)
    
    net <- network::network(netmat, matrix.type = "edgelist", directed = FALSE)
    
    network::set.edge.attribute(net, "antag", filt$antag)
    
    sz <- network::network.size(net)
    
    size[[i]] <- tibble(
      size = network::network.size(net),
      id = i)
    
    key <- tibble(
      antag = factor(c(1, 2, 3)),
      id = i)
    
    prop[[i]] <-
      as.data.frame(table(network::get.edge.attribute(net, attrname = 'antag')) /
                      length(network::get.edge.attribute(net, attrname = 'antag'))) %>%
      setNames(., c("antag", "prop")) %>%
      add_column(id = i) %>% 
      full_join(key) %>% 
      mutate(prop = ifelse(is.na(prop), 0, prop))
    
    if( totsp == sz ) {
      connectance[[i]] <- tibble(
        connectance = network.density(net),
        id = i)
      
    } else {
      
      nsp_to_add <- totsp - sz
      
      add.vertices(net, nsp_to_add)
      
      connectance[[i]] <- tibble(
        connectance = network.density(net),
        id = i)
    }
  }
}

site_info <- d %>% 
  ungroup() %>% 
  dplyr::select(site_id, ghm) %>% 
  distinct(.)

info <- d2 %>% 
  dplyr::select(site_id, period, season, res, id) %>% 
  distinct(.) %>% 
  left_join(site_info)

n_sp <- do.call(rbind, size)

all <- do.call(rbind, prop) %>% 
  left_join(info) %>% 
  left_join(n_sp) %>% 
  dplyr::select(site_id, season, period, res, antag, prop, ghm, size)

info_ord <- info %>% ungroup() %>% arrange(id)

con <- do.call(rbind, connectance) %>% 
  dplyr::select(connectance) %>% 
  cbind(info_ord) %>% 
  as_tibble() %>% 
  dplyr::select(site_id, season, period, res, connectance, ghm)

setwd(here::here("data"))
# write_csv(con, "network_connectance_v01.csv")
# write_csv(all, "network_antagonism_proportions_v02.csv")

