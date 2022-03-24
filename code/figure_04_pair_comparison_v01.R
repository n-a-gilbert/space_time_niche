library(here)
library(tidyverse)
library(brms)
library(tidybayes)
library(MetBrewer)

setwd(here::here("results"))
mod <- readRDS("pair_analysis_results_v01.rds")

setwd(here::here("data"))
d <- read_csv("space_time_pair_ttd_data_v01.csv")

# get predictions "globally" for each antagonism ranking
global <- expand.grid(
  x = c(min(d$x), max(d$x)), 
  antag = unique(d$antag)) %>% 
  # dplyr::mutate(x2 = x*x) %>% 
  tidybayes::add_linpred_draws(mod,
                               re_formula = NA, 
                               scale = "response", 
                               n = 1000) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::mutate(x = ifelse(x == min(x), "Low", "High")) %>% 
  dplyr::mutate(value = .value / 1440) %>% # convert minutes to days
  dplyr::select(x, antag, draw = .draw, value) %>% 
  tidyr::pivot_wider(names_from = x, values_from = value) %>% 
  dplyr::mutate(difference = Low - High) %>% 
  dplyr::mutate(antag = ifelse(antag == 1, "P(Lethal Encounter): Low", 
                               ifelse(antag == 2, "P(Lethal Encounter): Medium", 
                                      "P(Lethal Encounter): High"))) %>% 
  mutate(antag = factor(antag,
                        levels = c(
                          "P(Lethal Encounter): High",
                          "P(Lethal Encounter): Medium",
                          "P(Lethal Encounter): Low"))) %>% 
  dplyr::group_by(antag) %>% 
  dplyr::summarise(mean = mean(difference), 
                   lower95 = quantile(difference, c(0.025)), 
                   upper95 = quantile(difference, c(0.975)),
                   lower68 = quantile(difference, c(0.160)), 
                   upper68 = quantile(difference, c(0.840)))

# for each pair, calculate the difference in simulated detection time between low- and high-disturbance sites
pair <- expand.grid(
  pair = unique(d$pair), 
  x = c(min(d$x), max(d$x))) %>% 
  # dplyr::mutate(x2 = x*x) %>% 
  dplyr::left_join( dplyr::distinct(dplyr::select(d, pair, antag))) %>%
  tibble::as_tibble() %>% 
  tidybayes::add_linpred_draws(mod, 
                               re_formula = ~(1 + x | pair), 
                               scale = "response", 
                               n = 1000) %>% 
  dplyr::ungroup(.) %>%
  dplyr::select(pair, antag, x, draw = .draw, value = .value) %>% 
  dplyr::mutate(x = ifelse(x == min(x), "Low", "High")) %>% 
  dplyr::mutate(value = value / 1440) %>% # convert minutes to days
  tidyr::pivot_wider(names_from = x, values_from = value) %>% 
  dplyr::mutate(difference = Low - High) %>% 
  dplyr::mutate(antag = ifelse(antag == 1, "P(Lethal Encounter): Low", 
                               ifelse(antag == 2, "P(Lethal Encounter): Medium", 
                                      "P(Lethal Encounter): High"))) %>% 
  dplyr::group_by(pair, antag) %>% 
  dplyr::summarise(mean = mean(difference),
                   lower95 = quantile(difference, c(0.025)), 
                   upper95 = quantile(difference, c(0.975)),
                   lower68 = quantile(difference, c(0.160)), 
                   upper68 = quantile(difference, c(0.840))) %>%
  dplyr::arrange(mean) %>% 
  dplyr::ungroup(.) %>% 
  # clean up some names
  dplyr::mutate(pair = str_replace(pair, "SQUIRRELSANDCHIPMUNKS", "SQUIRRELS"),
                pair = str_replace(pair, "FOXRED", "RED FOX"), 
                pair = str_replace(pair, "SKUNKSTRIPED", "SKUNK"),
                pair = str_replace(pair, "SNOWSHOEHARE", "HARE"),
                pair = str_replace(pair, "CRANESANDHILL", "CRANE")) %>% 
  dplyr::mutate(pair = tolower(pair)) %>% 
  # factor pair so they are arranged in order of effect size on y-axis
  dplyr::mutate(pair = factor(pair, levels = unique(pair))) %>% 
  mutate(antag = factor(antag,
                        levels = c(
                          "P(Lethal Encounter): High",
                          "P(Lethal Encounter): Medium",
                          "P(Lethal Encounter): Low")))

pal <-  MetBrewer::MetPalettes$Demuth[[1]][c(8, 4, 2)]

ggplot() +
  geom_vline(xintercept = 0, color = "red") + 
  geom_rect(data = global,
            aes(xmin = lower95, 
                xmax = upper95,
                fill = antag), 
            ymin = -Inf, 
            ymax = Inf,
            color = NA, 
            alpha = 0.2) +
  scale_fill_manual("P(Lethal Encounter)",
                    values = pal, 
                    labels = c("High", "Medium", "Low")) +
  geom_errorbar(data = pair, 
                aes(y = pair, 
                    xmin = lower68, 
                    xmax = upper68,
                    color = antag),
                width = 0, size = 1) + 
  geom_errorbar(data = pair, 
                aes(y = pair, 
                    xmin = lower95,
                    xmax = upper95,
                    color = antag), 
                width = 0, size = 0.5) +
  geom_point(data = pair, 
             aes(y = pair, 
                 x = mean, 
                 color = antag)) +
  scale_color_manual("P(Lethal Encounter)",
                     values = pal, 
                     labels = c("High", "Medium", "Low")) + 
  facet_grid(antag ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Difference in detection time between\n low- and high-disturbance sites (days)") + 
  theme(axis.title.y = element_blank(),
        axis.line = element_line(size = 0.1, 
                                 color = "gray20"), 
        axis.ticks = element_line(size = 0.1, 
                                  color = "gray20"),
        panel.background = element_rect(color = NA, 
                                        fill = "gray90"),
        panel.grid.minor.x = element_blank(),
        axis.text = element_text(size = 8.5, 
                                 color = "black"), 
        axis.title = element_text(size = 10, 
                                  color = "black"),
        legend.title = element_text(size = 10, 
                                    color = "black"), 
        legend.text = element_text(size = 8.5, 
                                   color = "black"),
        strip.background = element_blank(), 
        strip.text = element_blank(),
        legend.position = "bottom",
        legend.justification = "left",
        legend.margin = margin(0, 0, 0, -50), 
        legend.box.margin = margin(-5, -10, -5, -10))
