library(here)
library(tidyverse)
library(tidybayes)
library(MetBrewer)

setwd(here::here("results"))
mod <- readRDS("pair_analysis_results_v01.rds")

setwd(here::here("data"))
d <- read_csv("space_time_pair_ttd_data_v01.csv")

global <- expand.grid(x = seq(from = min(d$x), 
                              to = max(d$x),
                              by = 0.1), 
                      antag = unique(d$antag)) %>% 
  tibble::as_tibble() %>% 
  tidybayes::add_linpred_draws(mod, 
                               re_formula = NA, 
                               scale = "response", 
                               n = 1000) %>% 
  mutate(.value = .value/ 1440) %>% 
  mutate(antag_lab = ifelse(antag < 2 , "Low",
                            ifelse(antag > 2, "High", "Medium"))) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(x, y = .value, antag_lab) %>% 
  dplyr::mutate(antag_lab = factor(antag_lab, levels = c("High", "Medium", "Low"))) %>% 
  dplyr::group_by(antag_lab, x) %>% 
  dplyr::summarise(mean = mean(y), 
                   lower = quantile(y, c(0.025)), 
                   upper = quantile(y, c(0.975)))

pal <-  MetBrewer::MetPalettes$Demuth[[1]][c(8, 4, 2)]
xvals <- c(min(global$x) + 0.25, max(global$x) - 0.25)

ggplot(data = global, aes(x = x, y = mean, color = antag_lab, fill = antag_lab)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), color = NA, alpha = 0.4) +
  geom_line(size = 1.5) +
  scale_color_manual("P(Lethal Encounter)",
                     values = pal) +
  scale_fill_manual("P(Lethal Encounter)",
                    values = pal) +
  scale_x_continuous(breaks = xvals,
                     labels = c("Low", "High")) +
  labs(y = "Time between detections (days)", 
       x = "Human disturbance (5km)") + 
  theme_classic()  +
  theme(legend.position = c(0.825, 0.8),
        axis.line = element_line(size = 0.1,
                                 color = "gray20"),
        axis.ticks = element_line(size = 0.1, 
                                  color = "gray20"), 
        axis.text = element_text(size = 10, 
                                 color = "black"), 
        axis.title = element_text(size = 11, 
                                  color = "black"), 
        legend.title = element_text(size = 9, 
                                    color = "black"),
        legend.text = element_text(size = 9, 
                                   color = "black"),
        legend.background = element_blank(),
        panel.background = element_rect(color = NA, 
                                        fill = "gray95")) 
