# code to make network connectance figures
#panel b of fig 5 and fig. s8

library(here)
library(tidyverse)
library(tidybayes)
library(MetBrewer)
library(lemon)

setwd(here::here("results/connectance_model_results"))

m10 <- readRDS("connectance_10minutes_v02.rds")
m60 <- readRDS("connectance_60minutes_v02.rds")
m1440 <- readRDS("connectance_1440minutes_v02.rds")
m10080 <- readRDS("connectance_10080minutes_v02.rds")
m43200 <- readRDS("connectance_43200minutes_v02.rds")

setwd(here::here("data"))

d <- read_csv("network_connectance_v01.csv") %>% 
  mutate(x = as.numeric(scale(ghm))) %>% 
  mutate(res = ifelse(res == "1hr", "1 hour",
                      ifelse(res == "10min", "10 minutes", 
                             ifelse(res == "1day", "1 day", 
                                    ifelse(res == "1week", "1 week", "1 month"))))) %>% 
  mutate(res = factor(res, 
                      levels = c(
                        "10 minutes", 
                        "1 hour", 
                        "1 day", 
                        "1 week", 
                        "1 month")))

center <- attr(scale(d$ghm), "scaled:center")
sc <- attr(scale(d$ghm), "scaled:scale")

pred <- expand.grid(x = seq(from = min(d$x),
                            to = max(d$x),
                            by = 0.1),
                    season = c("summer", "winter"),
                    size = 6)

m10_pred <- 
  add_linpred_draws(
    pred, 
    m10, 
    scale = "response", 
    n = 1000
  ) %>% 
  ungroup() %>% 
  group_by(season, x) %>% 
  summarise(mean = mean(.value), 
            l95 = quantile(.value, c(0.025)), 
            u95 = quantile(.value, c(0.975))) %>% 
  ungroup() %>% 
  mutate(ghm = x * sc + center) %>% 
  add_column(res = "10 minutes")

m60_pred <- 
  add_linpred_draws(
    pred, 
    m60, 
    scale = "response", 
    n = 1000
  ) %>% 
  ungroup() %>% 
  group_by(season, x) %>% 
  summarise(mean = mean(.value), 
            l95 = quantile(.value, c(0.025)), 
            u95 = quantile(.value, c(0.975))) %>% 
  ungroup() %>% 
  mutate(ghm = x * sc + center) %>% 
  add_column(res = "1 hour")

m1440_pred <- 
  add_linpred_draws(
    pred, 
    m1440, 
    scale = "response", 
    n = 1000
  ) %>% 
  ungroup() %>% 
  group_by(season, x) %>% 
  summarise(mean = mean(.value), 
            l95 = quantile(.value, c(0.025)), 
            u95 = quantile(.value, c(0.975))) %>% 
  ungroup() %>% 
  mutate(ghm = x * sc + center) %>% 
  add_column(res = "1 day")

m10080_pred <- 
  add_linpred_draws(
    pred, 
    m10080, 
    scale = "response", 
    n = 1000
  ) %>% 
  ungroup() %>% 
  group_by(season, x) %>% 
  summarise(mean = mean(.value), 
            l95 = quantile(.value, c(0.025)), 
            u95 = quantile(.value, c(0.975))) %>% 
  ungroup() %>% 
  mutate(ghm = x * sc + center) %>% 
  add_column(res = "1 week")

m43200_pred <- 
  add_linpred_draws(
    pred, 
    m43200, 
    scale = "response", 
    n = 1000
  ) %>% 
  ungroup() %>% 
  group_by(season, x) %>% 
  summarise(mean = mean(.value), 
            l95 = quantile(.value, c(0.025)), 
            u95 = quantile(.value, c(0.975))) %>% 
  ungroup() %>% 
  mutate(ghm = x * sc + center) %>% 
  add_column(res = "1 month")

pred_all <- full_join(m10_pred, m60_pred) %>% 
  full_join(m1440_pred) %>% 
  full_join(m10080_pred) %>% 
  full_join(m43200_pred) %>% 
  mutate(res = factor(res, 
                      levels = c(
                        "10 minutes", 
                        "1 hour", 
                        "1 day", 
                        "1 week", 
                        "1 month")))

pal <- MetPalettes$Hiroshige[[1]][c(1, 10)]

# p <-
  ggplot() + 
  geom_point(data = d,
    # data = filter(d, res == "1 day"), 
             aes(x = ghm, 
                 y = connectance#, 
                 # color = season
                 ), 
             color = "gray60",
             alpha = 0.1) +
  geom_ribbon(data = pred_all,
    # data = filter(pred_all, res == "1 day"), 
              aes(x = ghm, 
                  ymin = l95, 
                  ymax = u95, 
                  fill = season), 
              color = NA, 
              alpha = 0.4) +
  geom_line(data = pred_all,
    # data = filter(pred_all, res == "1 day"),
            aes(x = ghm, 
                y = mean, 
                color = season), 
            size = 1.5) +
  facet_wrap(~res) +
  theme_classic()  +
  labs(x = "Human disturbance (5km)", 
       y = "Network connectance",
       color = "Season",
       fill = "Season") +
  scale_color_manual(values = pal,
                     labels = c("Summer", "Winter")) +
  scale_fill_manual(values = pal,
                    labels = c("Summer", "Winter")) +
  theme(
    legend.position = "bottom",
    axis.line = element_line(size = 0.1,
                             color = "gray20"),
    axis.ticks = element_line(size = 0.1, 
                              color = "gray20"), 
    axis.text = element_text(size = 10, 
                             color = "black"), 
    axis.title = element_text(size = 11, 
                              color = "black"), 
    # strip.text = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 9, 
                               color = "black"),
    legend.background = element_blank(),
    strip.background = element_rect(color = NA), 
    panel.background = element_rect(color = NA, 
                                    fill = "white"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(-5, 0, 0, 0)) #+ 
  