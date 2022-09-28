# visualize predicted difference in time-between-detections for low- and high
# disturbance landscapes, at pair level
# figure 3

library(tidyverse)
library(brms)
library(here)
library(MetBrewer)
library(tidybayes)

setwd(here::here("results"))

res <- readRDS("ghm5k_antag_v01.rds")

setwd(here::here("data"))

load("supplement_pairwise_data_v01.RData")

d <- d %>% 
  mutate(ghm  = as.numeric(scale(ghm)),
         rai = as.numeric( scale (log1p (ntot / active ))),
         ghm5k_sd = as.numeric(scale(ghm5k_sd))) %>% 
  mutate(x = ghm)

# calculate global comparison for each antag level to put in on the background
global_comparison <- 
  expand.grid(x = c(min(d$x), max(d$x)), 
              antag = unique(d$antag)) %>% 
  # mutate(x2 = x*x) %>% 
  add_fitted_draws(res, 
                   re_formula = NA,
                   scale = "response", 
                   n = 1000) %>% 
  ungroup(.) %>% 
  mutate(x = ifelse(x == min(x), "Low", "High")) %>% 
  mutate(value = .value / 1440) %>% 
  dplyr::select(x, antag, draw = .draw, value) %>% 
  pivot_wider(names_from = x, values_from = value) %>% 
  mutate(difference = Low - High) %>% 
  mutate(antag = ifelse(antag == 1, "P(Lethal Encounter): Low", 
                        ifelse(antag == 2, "P(Lethal Encounter): Medium", 
                               "P(Lethal Encounter): High"))) %>% 
  mutate(antag = factor(antag,
                        levels = c(
                          "P(Lethal Encounter): High",
                          "P(Lethal Encounter): Medium",
                          "P(Lethal Encounter): Low"))) %>% 
  group_by(antag) %>% 
  summarise(mean_diff = mean(difference), 
            lower = quantile(difference, c(0.025)), 
            upper = quantile(difference, c(0.975)),
            sd_diff = sd(difference))

# for each pair, calculate difference in simulated TTE between low- and high-disturbance sites
plot_comp <-
  expand.grid( x = c(min(d$x), max(d$x)) ,
               pair = unique(d$pair) ) %>% 
  as_tibble() %>% 
  full_join(distinct(select(d, c(pair, antag)))) %>% 
  add_linpred_draws(res,
                    re_formula = ~(1 + x | pair),
                    scale = "response", 
                    n = 1000) %>% 
  ungroup(.) %>% 
  dplyr::select(pair, antag, x, draw = .draw, value = .value) %>% 
  mutate(x = ifelse(x == min(x), "Low", "High")) %>% 
  mutate(value = value / 1440) %>%
  pivot_wider(names_from = x, values_from = value) %>% 
  mutate(difference = Low - High) %>% 
  group_by(pair) %>% 
  summarise(mean_diff = mean(difference),
            sd_diff = sd(difference),
            lower95 = quantile(difference, c(0.025)), 
            upper95 = quantile(difference, c(0.975)),
            lower68 = quantile(difference, c(0.160)), 
            upper68 = quantile(difference, c(0.840))) %>%
  arrange(mean_diff) %>% 
  full_join(distinct(select(d, c(pair, antag)))) %>% 
  mutate(antag = ifelse(antag == 1, "P(Lethal Encounter): Low", 
                        ifelse(antag == 2, "P(Lethal Encounter): Medium", 
                               "P(Lethal Encounter): High"))) %>% 
  ungroup(.) %>% 
  mutate(pair = str_replace(pair, "SQUIRRELSANDCHIPMUNKS", "SQUIRRELS"),
         pair = str_replace(pair, "FOXRED", "RED FOX"), 
         pair = str_replace(pair, "SKUNKSTRIPED", "SKUNK"),
         pair = str_replace(pair, "SNOWSHOEHARE", "HARE"),
         pair = str_replace(pair, "CRANESANDHILL", "CRANE")) %>% 
  mutate(pair = tolower(pair)) %>% 
  # factor pair so they are arranged in order of effect size on y-axis
  mutate(pair = factor(pair, levels = unique(pair))) %>% 
  mutate(antag = factor(antag,
                        levels = c(
                          "P(Lethal Encounter): High",
                          "P(Lethal Encounter): Medium",
                          "P(Lethal Encounter): Low")))

# pal <- viridisLite::inferno(n = 3, begin = 0.7, end = 0.15)
pal <-  MetBrewer::MetPalettes$Demuth[[1]][c(8, 4, 2)]

ggplot() +
  geom_vline(xintercept = 0, color = "red") + 
  geom_rect(data = global_comparison,
            aes(xmin = lower, 
                xmax = upper,
                fill = antag), 
            ymin = -Inf, 
            ymax = Inf,
            color = NA, 
            alpha = 0.2) +
  scale_fill_manual("Antagonism",
                    values = pal, 
                    labels = c("High", "Medium", "Low")) +
  geom_errorbar(data = plot_comp, 
                aes(y = pair, 
                    xmin = lower68, 
                    xmax = upper68,
                    color = antag),
                width = 0, size = 1) + 
  geom_errorbar(data = plot_comp, 
                aes(y = pair, 
                    xmin = lower95,
                    xmax = upper95,
                    color = antag), 
                width = 0, size = 0.5) +
  geom_point(data = plot_comp, 
             aes(y = pair, 
                 x = mean_diff, 
                 color = antag)) +
  scale_color_manual("Antagonism",
                     values = pal, 
                     labels = c("High", "Medium", "Low")) + 
  # facet_wrap(~antag, scales = "free_y", ncol = 1) +
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
        legend.margin = margin(0, 0, 0, -40), 
        legend.box.margin = margin(-5, -10, -5, -10))