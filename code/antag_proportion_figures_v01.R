# visualize antagonism proportion predictions from models
# part of figure 5 (panel c), fig. s9


library(here)
library(tidyverse)
library(brms)
library(tidybayes)
library(MetBrewer)

setwd(here::here("results/proportion_model_results"))

df <- list.files(pattern = "v01.rds") %>% 
  map(readRDS)

setwd(here::here("data"))

res_key <- tibble(
  res = c(10, 60, 1440, 10080, 43200), 
  name = c("10 minutes", "1 hour", "1 day", "1 week", "1 month")
) %>% 
  mutate(name = factor(name, 
                       levels = c(
                         "10 minutes", 
                         "1 hour", 
                         "1 day", 
                         "1 week", 
                         "1 month")))

p <- read_csv("network_antag_proportions_v02.csv") %>% 
  group_by(antag, res) %>% 
  mutate(index_id = cur_group_id()) %>% 
  arrange(index_id) %>% 
  mutate(x = as.numeric(scale(ghm))) %>% 
  full_join(res_key) %>% 
  mutate(antag_name = ifelse(antag == 1, "Low antagonism", 
                             ifelse(antag == 2, "Medium antagonism", 
                                    "High antagonism"))) %>% 
  mutate(antag_name = factor(antag_name, 
                             levels = c(
                               "Low antagonism", 
                               "Medium antagonism", 
                               "High antagonism")))


index_id_key <- p %>% 
  dplyr::select(antag, res, index_id) %>% 
  distinct(.)

pred <- expand.grid(x = seq(from = min(p$x),
                            to = max(p$x),
                            by = 0.1),
                    season = c("summer", "winter"),
                    size = 6)

center <- attr(scale(p$ghm), "scaled:center")
sc <- attr(scale(p$ghm), "scaled:scale")

pred_list <- list(list())
for(i in 1:max(index_id_key$index_id)){
  
  d <- p %>% 
    filter(index_id == i)
  
  antag_lab <- unique(d$antag)
  res_lab <- unique(d$res)
  
  pred_list[[i]] <- add_linpred_draws(
    pred, 
    df[[i]], 
    scale = "response", 
    n = 1000 ) %>% 
    ungroup() %>% 
    group_by(season, x) %>% 
    summarise(mean = mean(.value), 
              l95 = quantile(.value, c(0.025)), 
              u95 = quantile(.value, c(0.975))) %>% 
    ungroup() %>% 
    mutate(ghm = x * sc + center) %>% 
    add_column(antag = antag_lab,
               res = res_lab)
}

all <- do.call(rbind, pred_list) %>% 
  full_join(res_key) %>% 
  mutate(antag_name = ifelse(antag == 1, "Low antagonism", 
                             ifelse(antag == 2, "Medium antagonism", 
                                    "High antagonism"))) %>% 
  mutate(antag_name = factor(antag_name, 
                             levels = c(
                               "Low antagonism", 
                               "Medium antagonism", 
                               "High antagonism")))


pal <- MetPalettes$Hiroshige[[1]][c(1, 10)]

ggplot() +
  geom_point(
    data = p,
    # data = filter(p, name == "1 day"),
             aes(x = ghm,
                 y = prop),
             alpha = 0.1,
             color = "gray60") +
  geom_ribbon(
    data = all,
    # data =  filter(all, name == "1 day"), 
              aes(x = ghm, 
                  ymin = l95, 
                  ymax = u95, 
                  fill = season), 
              color = NA, 
              alpha = 0.4) +
  geom_line(
    data = all, 
    # data =  filter(all, name == "1 day"),
            aes(x = ghm, 
                y = mean,
                color = season), 
            size = 1.5, 
            alpha = 1) +
  facet_grid(name ~ antag_name) +
  theme_classic() +
  labs(x = "Human disturbance (5km)", 
       y = "Proportion of links in network",
       color = "Season", 
       fill = "Season") +
  scale_color_manual(values = pal,
                     labels = c("Summer", "Winter")) +
  scale_fill_manual(values = pal,
                    labels = c("Summer", "Winter")) +
  theme(
    axis.line = element_line(size = 0.1,
                             color = "gray20"),
    axis.ticks = element_line(size = 0.1, 
                              color = "gray20"), 
    axis.text = element_text(size = 10, 
                             color = "black"), 
    # strip.text.y = element_blank(),
    strip.text.x = element_text(size = 11),
    axis.title = element_text(size = 11, 
                              color = "black"), 
    legend.title = element_blank(),
    legend.text = element_text(size = 9, 
                               color = "black"),
    legend.background = element_blank(),
    strip.background = element_rect(color = NA), 
    panel.background = element_rect(color = NA, 
                                    fill = "white"))

