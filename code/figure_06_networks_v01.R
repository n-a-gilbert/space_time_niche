library(here)
library(tidyverse)
library(tidybayes)
library(wesanderson)
library(network)
library(ggnetwork)
library(viridisLite)
library(patchwork)
library(cowplot)

# color palette for seasons: blue for winter, red for summer
pal <- wesanderson::wes_palette("Zissou1")[c(5, 1)]

# color palette for antagonism
antagpal <- viridisLite::inferno(n = 3, begin = 0.15, end = 0.7)

# gather materials (data, model results) to create plots

setwd(here::here("results"))

load("network_analysis_results_v01.RData")

setwd(here::here("data"))

netdata <- readr::read_csv("df_to_make_example_networks_v01.csv") %>% 
  dplyr::mutate(sp1 = tolower(sp1), 
         sp2 = tolower(sp2)) %>% 
  dplyr::mutate(tte = 1 / log(tte)) # transform average time to detection for visualization

d <- readr::read_csv("space_time_network_density_v01.csv")

p <- readr::read_csv("space_time_network_proportions_v01.csv") %>% 
  dplyr::mutate(antag_lab = ifelse(antag == 3, "P(Lethal Encounter): High", 
                            ifelse(antag == 2, "P(Lethal Ecounter): Medium", 
                                   "P(Lethal Ecounter): Low")))

# panel a - example of a network
# separate data for the two examples (low and high disturbance)
d1 <- netdata %>% 
  filter(ghm_class == "Low")

d2 <- netdata %>% 
  filter(ghm_class == "High")

net1 <- network::network(cbind(d1$sp1, d1$sp2), matrix.type = "edgelist", directed = FALSE)
set.edge.attribute(net1, "antag", d1$antag)
set.edge.attribute(net1, "tte", d1$tte)
netplot1 <- ggnetwork::ggnetwork(net1)

( low <- ggplot(netplot1, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(aes( color = factor(antag,
                                   levels = c(1, 2, 3)),
                    size = tte)) +
    scale_size_continuous(range = c(1, 5)) +
    guides(size = "none") +
    scale_color_manual("P(Lethal Encounter)",
                       labels = c("Low", "Medium", "High"),
                       values = antagpal) +
    geom_nodes(color = "gray30", size = 4) +
    geom_nodelabel(aes(label = vertex.names),
                   fill = "gray90",
                   size = 3.5) +
    scale_x_continuous(limits = c(-0.1, 1.1)) +
    theme_void() + 
    theme(legend.position = "none") )

net2 <- network::network(cbind(d2$sp1, d2$sp2), matrix.type = "edgelist", directed = FALSE)
network::set.edge.attribute(net2, "antag", d2$antag)
network::set.edge.attribute(net2, "tte", d2$tte)
netplot2 <- ggnetwork::ggnetwork(net2)

( high <- ggplot(netplot2, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(aes( color = factor(antag,
                                   levels = c(1, 2, 3)),
                    size = tte)) +
    scale_size_continuous(range = c(1, 5)) +
    guides(size = "none") +
    scale_color_manual("P(Lethal Encounter)",
                       labels = c("Low", "Medium", "High"),
                       values = antagpal) +
    geom_nodes(color = "gray30", size = 4) +
    geom_nodelabel(aes(label = vertex.names),
                   fill = "gray90",
                   size = 3.5) +
    scale_x_continuous(limits = c(-0.1, 1.1)) +
    theme_void() + 
    theme(legend.position = "none") )

( for_legend <- ggplot(netplot2, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(aes( color = factor(antag,
                                   levels = c(1, 2, 3)),
                    size = tte)) +
    scale_size_continuous(range = c(1, 5)) +
    guides(size = "none") +
    scale_color_manual("P(Lethal Encounter)",
                       labels = c("Low", "Medium", "High"),
                       values = antagpal) +
    geom_nodes(color = "gray30", size = 4) +
    geom_nodelabel(aes(label = vertex.names),
                   fill = "gray90",
                   size = 3.5) +
    scale_x_continuous(limits = c(-0.1, 1.1)) +
    theme_void() + 
    theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 2))) )

legend <- cowplot::get_legend(for_legend)

# Panel A - example networks
( low + high ) / legend + patchwork::plot_layout(heights = c(10, 1))

# create some synthetic data to visualize model predictions
pred_data <- expand.grid(
  x = seq(from = min(d$x),
          to = max(d$x), 
          by = 0.1), 
  season = c("summer", "winter")
) %>% 
  dplyr::mutate(x2 = x*x)

# use model to predict (average values) for the synethic data
# use 1000 draws, then summarise mean and 95% CI for plotting
density_pred <- tidybayes::add_linpred_draws(
  pred_data, 
  density_model, 
  scale = "response", 
  n = 1000
) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::select(season, x, density = .value) %>% 
  dplyr::group_by(season, x) %>% 
  dplyr::summarise(mean = mean(density), 
                   lower = quantile(density, c(0.025)), 
                   upper = quantile(density, c(0.975)))

# network density figure - panel b
( density_plot <- ggplot() +
    geom_point(data = d,
               aes(x = x, y = density, color = season), 
               size = 3, 
               alpha = 0.1) + 
    geom_ribbon(data = density_pred, 
                aes(x = x, ymin = lower, ymax = upper, fill = season), 
                color = NA, 
                alpha = 0.4) +
    geom_line(data = density_pred,
              aes(x = x, y = mean, color = season), 
              size = 2) +
    scale_color_manual("Season", 
                       values = pal,
                       labels = c("Summer", "Winter")) + 
    scale_fill_manual("Season", 
                      values = pal,
                      labels = c("Summer", "Winter")) + 
    scale_x_continuous("Human disturbance (5km)",
                       breaks = c(-1.73, 2.97), 
                       labels = c("Low", "High")) +
    labs(x = "Human disturbance (5km)", 
         y = "Network density") +
    theme_classic() +
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
          panel.background = element_rect(color = NA, 
                                          fill = "gray95"),
          legend.position = "bottom") )


# add predicted average values for proportion of each antagonism ranking in network

# low
prop_antag1 <- tidybayes::add_linpred_draws(
  pred_data, 
  prop1_model, 
  scale = "response", 
  n = 1000) %>% 
  dplyr::ungroup(.) %>% 
  tibble::add_column(antag = 1) %>% 
  dplyr::select(antag, x, season, draw = .draw, prop = .value) %>% 
  dplyr::group_by(antag, season, x) %>%
  dplyr::summarise(mean = mean(prop), 
            lower = quantile(prop, c(0.025)), 
            upper = quantile(prop, c(0.975)))

# medium
prop_antag2 <- tidybayes::add_linpred_draws(
  pred_data, 
  prop2_model, 
  scale = "response", 
  n = 1000) %>% 
  dplyr::ungroup(.) %>% 
  tibble::add_column(antag = 2) %>% 
  dplyr::select(antag, x, season, draw = .draw, prop = .value) %>% 
  dplyr::group_by(antag, season, x) %>%
  dplyr::summarise(mean = mean(prop), 
            lower = quantile(prop, c(0.025)), 
            upper = quantile(prop, c(0.975)))

#high
prop_antag3 <- tidybayes::add_linpred_draws(
  pred_data, 
  prop3_model, 
  scale = "response", 
  n = 1000) %>% 
  dplyr::ungroup(.) %>% 
  tibble::add_column(antag = 3) %>% 
  dplyr::select(antag, x, season, draw = .draw, prop = .value) %>% 
  dplyr::group_by(antag, season, x) %>%
  dplyr::summarise(mean = mean(prop), 
            lower = quantile(prop, c(0.025)), 
            upper = quantile(prop, c(0.975)))

# join them all together
prop_all <- dplyr::full_join(prop_antag1, prop_antag2) %>% 
  dplyr::full_join(prop_antag3) %>% 
  dplyr::mutate(antag_lab = ifelse(antag == 3, "P(Lethal Encounter): High", 
                            ifelse(antag == 2, "P(Lethal Ecounter): Medium", 
                                   "P(Lethal Ecounter): Low")))

# plot showing proprtions and model predictions  - panel c
( prop_plot <- ggplot() + 
    geom_point(data = p, 
               aes(x = x, y = prop, color = season), 
               size = 3, 
               alpha = 0.025) +
    geom_ribbon(data = prop_all, 
                aes(x = x, ymin = lower, ymax = upper, fill = season), 
                color = NA, 
                alpha = 0.3) +
    geom_line(data = prop_all, 
              aes(x = x, y = mean, color = season),
              size = 2) + 
    facet_wrap(~antag_lab) +
    scale_color_manual("Season",
                       values = pal,
                       labels = c("Summer", "Winter")) +
    scale_fill_manual("Season",
                      values = pal,
                      labels = c("Summer", "Winter")) +
    scale_x_continuous("Human disturbance (5km)",
                       breaks = c(-1.73, 2.97), 
                       labels = c("Low", "High")) +
    ylab("Proportion of links in network") +
    theme_classic() +
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
          panel.background = element_rect(color = NA, 
                                          fill = "gray95"),
          strip.background = element_blank()) )
