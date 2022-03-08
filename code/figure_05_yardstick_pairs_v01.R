library(here)
library(tidyverse)
library(brms)
library(tidybayes)
library(wesanderson)

setwd(here::here("results"))
mod <- readRDS("pair_analysis_results_v01.rds")

setwd(here::here("data"))
d <- readr::read_csv("space_time_pair_ttd_data_v01.csv")

yardsticks <- c(
  "COYOTE-DEER",
  "COTTONTAIL-COYOTE", 
  "OPOSSUM-RACCOON",
  "COYOTE-FOXRED",
  "DEER-WOLF",
  "RACCOON-SKUNKSTRIPED",
  "COTTONTAIL-FOXRED",
  "COYOTE-WOLF",
  "OPOSSUM-SKUNKSTRIPED")

yardstick_pred <- expand.grid(pair = yardsticks,
                              season = c("winter", "summer"),
                              x = seq(from = min(d$x), to = max(d$x), by = 0.2)) %>% 
  tibble::as_tibble() %>% 
  dplyr::mutate(x2 = x*x) %>% 
  dplyr::left_join(dplyr::distinct(dplyr::select(d, pair, antag))) %>% 
  dplyr::filter(! (pair == "BEAR-BOBCAT" & season == "winter")) %>% 
  dplyr::filter(! (pair == "BEAR-PORCUPINE" & season == "winter")) %>% 
  dplyr::filter(! (pair == "BEAR-RACCOON" & season == "winter")) %>% 
  dplyr::filter(! (pair == "BEAR-TURKEY" & season == "winter")) %>%   
  tidybayes::add_linpred_draws(mod, 
                               scale = "response", 
                               n = 1000) %>% 
  dplyr::mutate(pair = str_replace(pair, "SQUIRRELSANDCHIPMUNKS", "SQUIRRELS"),
                pair = str_replace(pair, "FOXRED", "RED FOX"), 
                pair = str_replace(pair, "SKUNKSTRIPED", "SKUNK"),
                pair = str_replace(pair, "SNOWSHOEHARE", "HARE"),
                pair = str_replace(pair, "CRANESANDHILL", "CRANE")) %>% 
  dplyr::mutate(pair = tolower(pair)) %>% 
  dplyr::mutate(.value = .value / 1440 ,
                season = stringr::str_to_sentence(season)) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(pair, season, x) %>% 
  dplyr::summarise(mean = mean(.value), 
                   lower = quantile(.value, c(0.025)), 
                   upper = quantile(.value, c(0.975)))

pal <- wesanderson::wes_palette("Zissou1")[c(5, 1)]

ggplot(data = yardstick_pred, 
       aes(x = x, y = mean, color = season)) +
  facet_wrap(~pair) + 
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = season), 
              color = NA, alpha = 0.45) + geom_line(size = 1.5) +
  scale_color_manual("Season", values = pal) +
  scale_fill_manual("Season", values = pal) +
  theme_classic()  +
  scale_x_continuous(breaks = c(-1.98, 2.36), 
                     labels = c("Low", "High")) + 
  labs(x = "Human disturbance (5km)",
       y = "Time between detections (days)") +
  theme(axis.line = element_line(size = 0.1,
                                 color = "gray20"),
        axis.ticks = element_line(size = 0.1, 
                                  color = "gray20"), 
        axis.text = element_text(size = 11, 
                                 color = "black"), 
        axis.title = element_text(size = 12, 
                                  color = "black"), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 11, 
                                   color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(size = 12, 
                                  color = "black"),
        panel.background = element_rect(color = NA, 
                                        fill = "gray95"),
        legend.position = "bottom",
        legend.margin = margin(0, 0, 0, 0), 
        legend.box.margin = margin(-5, -10, -5, -10))