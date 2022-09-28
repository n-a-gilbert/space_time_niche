# code to visualize predictions from model
# and recreate fig 2

library(here)
library(tidyverse)
library(brms)
library(tidybayes)

setwd(here::here("results"))

# primary model
m1 <- readRDS("ghm5k_antag_v01.rds")
# model wtih ghm heterogeneity
m2 <- readRDS("ghm5k_mean_sd_antag_v01.rds")
# model with rai
m3 <- readRDS("ghm5k_antag_rai_v01.rds")

# prep data
setwd(here::here("data"))

load("supplement_pairwise_data_v01.RData")

center <- attr(scale(d$ghm), "scaled:center")
sc <- attr(scale(d$ghm), "scaled:scale")

d <- d %>% 
  mutate(ghm  = as.numeric(scale(ghm)),
         rai = as.numeric( scale (log1p (ntot / active ))),
         ghm5k_sd = as.numeric(scale(ghm5k_sd))) %>% 
  mutate(x = ghm)

predm1 <- expand.grid(
  antag = unique(d$antag), 
  x = seq(from = min(d$x, na.rm = TRUE), to = max(d$x, na.rm = TRUE), by = 0.1)) %>% 
  add_linpred_draws(model = m1,
                    re_formula = NA, 
                    scale = "response", 
                    n = 1000) %>% 
  mutate(.value = .value/ 1440) %>% 
  mutate(antag_lab = ifelse(antag < 2 , "Low",
                            ifelse(antag > 2, "High", "Medium"))) %>% 
  mutate(antag_lab = factor(antag_lab, levels = c("Low", "Medium", "High"))) %>%
  ungroup() %>% 
  dplyr::select(x, antag = antag_lab, value = .value) %>% 
  group_by(antag, x) %>% 
  summarise(mean = mean(value), 
            lower = quantile(value, c(0.025)), 
            upper = quantile(value, c(0.975))) %>% 
  mutate(ghm = x * sc + center)

pal <-  MetBrewer::MetPalettes$Demuth[[1]][c(2, 4, 8)]

predm2 <- expand.grid(
  antag = unique(d$antag), 
  x = seq(from = min(d$x),
          to = max(d$x), by = 0.1),
  ghm5k_sd = c(min(d$ghm5k_sd),
               max(d$ghm5k_sd))) %>% 
  add_linpred_draws(model = m2,
                    re_formula = NA, 
                    scale = "response", 
                    n = 1000) %>% 
  mutate(.value = .value/ 1440) %>% 
  mutate(antag_lab = ifelse(antag < 2 , "Low",
                            ifelse(antag > 2, "High", "Medium"))) %>% 
  mutate(antag_lab = factor(antag_lab, levels = c("Low", "Medium", "High"))) %>%
  ungroup() %>% 
  dplyr::select(x, x_sd = ghm5k_sd, antag = antag_lab, value = .value) %>% 
  group_by(antag, x, x_sd) %>% 
  summarise(mean = mean(value), 
            lower = quantile(value, c(0.025)), 
            upper = quantile(value, c(0.975))) %>% 
  mutate(ghm = x * sc + center,
         het = ifelse(x_sd < 0, "Low", "High"),
         het = factor(het, levels = c("Low", "High")))

predm3 <-
  expand.grid(
    antag = unique(d$antag), 
    x = seq(from = min(d$x),
            to = max(d$x), by = 0.1),
    rai = c(quantile(d$rai, c(0.10)),
            mean(d$rai),
            quantile(d$rai, c(0.90)))) %>% 
  add_linpred_draws(model = m3,
                    re_formula = NA, 
                    scale = "response", 
                    n = 1000) %>% 
  mutate(.value = .value/ 1440) %>% 
  mutate(antag_lab = ifelse(antag < 2 , "Low",
                            ifelse(antag > 2, "High", "Medium"))) %>% 
  mutate(antag_lab = factor(antag_lab, levels = c("Low", "Medium", "High"))) %>%
  ungroup() %>% 
  dplyr::select(x, rai, antag = antag_lab, value = .value) %>% 
  group_by(antag, x, rai) %>% 
  summarise(mean = mean(value), 
            lower = quantile(value, c(0.025)), 
            upper = quantile(value, c(0.975))) %>% 
  mutate(ghm = x * sc + center,
         scenario = ifelse(rai == min(rai), "Abundance: Low",
                           ifelse(rai == max(rai), "Abundance: High", 
                                  "Abundance: Mean"))) %>% 
  dplyr::select(-rai) %>% 
  ungroup()


test <- predm2 %>% 
  ungroup() %>% 
  dplyr::select(ghm, scenario = het, antag, mean:upper) %>% 
  mutate(scenario = ifelse(scenario == "Low", "Heterogeneity: Low", "Heterogeneity: High")) %>% 
  full_join( add_column(select(ungroup(predm1), c(ghm, antag, mean:upper)),
                        scenario = "Primary")) %>% 
  full_join(predm3)

scenario_key <- test %>% 
  group_by(scenario) %>% 
  slice(1) %>% 
  dplyr::select(scenario) %>% 
  add_column(label = paste0("(", c("f", "d", "e", "c", "b", "a"), ")"))

( panela <-
    test %>% 
    full_join(scenario_key) %>% 
    mutate(scenario_name = paste(label, scenario, sep = "  "),
           scenario_name = factor(scenario_name,
                                  levels = c(
                                    "(a)  Primary", 
                                    "(b)  Heterogeneity: Low", 
                                    "(c)  Heterogeneity: High",
                                    "(d)  Abundance: Low", 
                                    "(e)  Abundance: Mean", 
                                    "(f)  Abundance: High"))) %>% 
    filter(scenario == "Primary") %>% 
    ggplot(aes(x = ghm, y = mean, color = antag)) + 
    # facet_wrap(~scenario_name) +
    geom_ribbon(aes(ymin = lower, 
                    ymax = upper, 
                    fill = antag),
                color = NA, 
                alpha = 0.3) +
    geom_line(size = 1) +
    scale_color_manual("Antagonism",
                       values = pal#,
                       # labels = c("Low", "Medium", "High")
    ) +
    scale_fill_manual("Antagonism",
                      values = pal#,
                      # labels = c("Low", "Medium", "High")
    ) +
    labs(y = "Time between detections (days)", 
         x = "Mean human disturbance (5km)",
         title = "(a)") + 
    theme_classic()  +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 9, color = "black"),
      axis.line = element_line(size = 0.2,
                               color = "gray20"),
      axis.ticks = element_line(size = 0.2, 
                                color = "gray20"), 
      axis.text = element_text(size = 8, 
                               color = "black"), 
      axis.title = element_text(size = 9, 
                                color = "black"), 
      legend.title = element_text(size = 9, 
                                  color = "black"),
      legend.text = element_text(size = 8, 
                                 color = "black"),
      legend.background = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = 9, color = "black"),
      legend.box.just = "center",
      panel.background = element_rect(color = NA, 
                                      fill = "white"),
      legend.key.size = unit(0.4, "cm"),
      legend.margin = margin(0, 0, 0, 0), 
      legend.box.margin = margin(-5, 0, 0,0)) +
    guides(color  = guide_legend(reverse = T), 
           fill = guide_legend(reverse = T)) )


scenario_key2 <- test %>% 
  group_by(scenario) %>% 
  slice(1) %>% 
  dplyr::select(scenario) %>% 
  filter(! (scenario == "Primary" | scenario == "Abundance: Mean")) %>% 
  add_column(label = paste0("(", c("e", "d", "c", "b"), ")"))


( panelb <-
    
    test %>% 
    full_join(scenario_key2) %>% 
    mutate(scenario_name = paste(label, scenario, sep = "  "),
           scenario_name = factor(scenario_name,
                                  levels = c(
                                    # "(a)  Primary", 
                                    "(b)  Heterogeneity: Low", 
                                    "(c)  Heterogeneity: High",
                                    "(d)  Abundance: Low", 
                                    # "(e)  Abundance: Mean", 
                                    "(e)  Abundance: High"))) %>% 
    # filter(scenario == "Primary") %>%
    # test %>% 
    filter(! (scenario == "Primary" | scenario == "Abundance: Mean")) %>%
    # mutate(scenario = factor(scenario, 
    #                          levels = c(
    #                            "Heterogeneity: Low", 
    #                            "Heterogeneity: High", 
    #                            "Abundance: Low", 
    #                            "Abundance: High"))) %>% 
    
    ggplot(aes(x = ghm, y = mean, color = antag)) + 
    facet_wrap(~scenario_name) +
    geom_ribbon(aes(ymin = lower, 
                    ymax = upper, 
                    fill = antag),
                color = NA, 
                alpha = 0.3) +
    geom_line(size = 1) +
    scale_color_manual("Antagonism",
                       values = pal#,
                       # labels = c("Low", "Medium", "High")
    ) +
    scale_fill_manual("Antagonism",
                      values = pal#,
                      # labels = c("Low", "Medium", "High")
    ) +
    labs(
      # title = "(b)",
      y = "Time between detections (days)", 
      x = "Mean human disturbance (5km)") + 
    theme_classic()  +
    theme(
      legend.position = "bottom",
      axis.line = element_line(size = 0.2,
                               color = "gray20"),
      axis.ticks = element_line(size = 0.2, 
                                color = "gray20"), 
      axis.text = element_text(size = 7.5, 
                               color = "black"), 
      axis.title = element_text(size = 8.5,
                                color = "black"),
      legend.title = element_text(size = 9, 
                                  color = "black"),
      legend.text = element_text(size = 8, 
                                 color = "black"),
      legend.background = element_blank(),
      strip.background = element_blank(),
      # plot.title = element_text(size = 9, color = "black"),
      strip.text = element_text(size = 8, color = "black"),
      legend.box.just = "center",
      panel.background = element_rect(color = NA, 
                                      fill = "white"),
      legend.key.size = unit(0.4, "cm"),
      legend.margin = margin(0, 0, 0, 0), 
      legend.box.margin = margin(-5, 0, 0,0)) +
    guides(color  = guide_legend(reverse = T),
           fill = guide_legend(reverse = T)) )


library(patchwork)
panela + panelb + plot_layout(ncol = 1, heights = 1)

setwd(here::here("figures/revision"))
ggsave("figure_02_v03.png", 
       width = 4,
       height = 5, 
       units = "in", 
       dpi = 300)