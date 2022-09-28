# code to create network visualization
# panel a of figure 5
library(here)
library(tidyverse)
library(ggnetwork)
library(network)

setwd(here::here("figures"))

d <- read_csv("df_to_make_example_networks_v01.csv") %>% 
  mutate(sp1 = tolower(sp1), 
         sp2 = tolower(sp2)) %>%
  mutate(tte = tte / 1440)
  
d1 <- d %>% 
  filter(ghm_class == "Low")

d2 <- d %>% 
  filter(ghm_class == "High")

net <- network(cbind(d1$sp1, d1$sp2), matrix.type = "edgelist", directed = FALSE)

set.edge.attribute(net, "antag", d1$antag)
set.edge.attribute(net, "tte", d1$tte)

p <- ggnetwork(net)

net2 <- network(cbind(d2$sp1, d2$sp2), matrix.type = "edgelist", directed = FALSE)

set.edge.attribute(net2, "antag", d2$antag)
set.edge.attribute(net2, "tte", d2$tte)

p2 <- ggnetwork(net2)

both <- p %>% 
  add_column(disturbance = "Low") %>% 
  full_join(add_column(p2, disturbance = "High")) %>% 
  mutate(disturbance = factor(disturbance, 
                              levels = c(
                                "Low", "High"
                              ))) %>% 
  mutate(tte_new = log1p(1 / tte)) %>% 
  group_by(disturbance) %>% 
  mutate(density = ifelse(disturbance == "Low", 0.31, 0.56),
         tte_mean = round(mean(tte, na.rm = TRUE), 1)) %>% 
  mutate(dist_lab = paste(disturbance, "Human Disturbance\n Connectance:", density,
                          "\nMean time:", tte_mean, "days")) %>% 
  mutate(dist_lab = factor(dist_lab, 
                           levels = unique(dist_lab)))

pal <-  MetBrewer::MetPalettes$Demuth[[1]][c(2, 4, 8)]

ggplot(both, aes(x = x, y = y, xend = xend, yend = yend)) +
  facet_wrap(~dist_lab) + 
  geom_edges(aes(color = factor(antag, 
                                levels = c(1, 2, 3)), 
                 size = tte_new)) +
  scale_size_continuous(range = c(1, 5)
                        ) +
  guides(size = "none") +
  scale_color_manual("Antagonism",
                     labels = c("Low", "Medium", "High"),
                     values = pal) +
  geom_nodes(color = "gray30", size = 3) +
  geom_nodelabel(aes(label = vertex.names),
                 fill = "gray90", 
                 size = 3) +
  scale_x_continuous(limits = c(-0.1, 1.1)) +
  theme_void() + 
  theme(legend.position = "bottom",
        strip.text = element_text(size = 10, color = "black"),
        legend.title = element_text(size = 10, color = "black"), 
        legend.text = element_text(size = 10, color = "black")) +
  guides(color = guide_legend(override.aes = list(size = 3)))
