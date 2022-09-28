library(overlap)
library(tidyverse)
library(here)

setwd(here::here("data"))

# this has two data frames

# pair_key is a key for the pairs to be considered
# detection_times contains the time of day (converted to suntime) of each species, 
# for the 20% most-disturbed (high) and least-disturbed (low) landscapes
# raw data cannot be shared, since geographic coordinates are private and are
# needed to convert the deteciton times to suntime
# detections of each species have been filtered so each is at least 30 mins from each other
load("detection_time_DAPs_v01.RData")

all_thin <- detection_times %>% 
  rename(lo_hi = disturbance)

# list to store results in
res <- list(list())

for(i in 1:nrow(pair_key)){

  # skip if only 1 member of pair is detected
  checkpoint <- all_thin %>%
    filter(sp == pair_key[[i, "sp1"]] | sp == pair_key[[i, "sp2"]]) %>%
    dplyr::select(sp) %>%
    distinct()

  if (nrow(checkpoint) < 2) next

  pair <- all_thin %>%
    filter(sp == pair_key[[i, "sp1"]] | sp == pair_key[[i, "sp2"]]) %>%
    dplyr::select(season, lo_hi, sp, suntime) %>%
    group_by(sp, lo_hi) %>%
    mutate(id = row_number()) %>%
    pivot_wider(names_from = sp, values_from = suntime) %>%
    setNames(.,
             c("season", "lo_hi", "id", "sp1", "sp2")) %>%
    ungroup() %>%
    group_by(lo_hi, season) %>%
    mutate(nsp1 = sum(!is.na(sp1)),
           nsp2 = sum(!is.na(sp2))) %>%
    filter(nsp1 > 0 & nsp2 > 0)

  # 1: summer, high ghm
  summer_hi <- pair %>%
    filter(season == "summer") %>%
    filter(lo_hi == "High")

  if(nrow(summer_hi) == 0){

    dhat_sum_hi = tibble(
      measure = NA,
      value = NA,
      season = "summer",
      lo_hi = "High",
      sp1 = pair_key[[i, "sp1"]],
      sp2 = pair_key[[i, "sp2"]],
      nsp1 = NA,
      nsp2 = NA)

  } else {

    dhat_sum_hi <-
      as_tibble( overlapEst(pull(filter(summer_hi, !is.na(sp1)), sp1),
                            pull(filter(summer_hi, !is.na(sp2)), sp2)),
                 rownames = "measure") %>%
      add_column(season = "summer",
                 lo_hi = "High",
                 sp1 = pair_key[[i, "sp1"]],
                 sp2 = pair_key[[i, "sp2"]],
                 nsp1 = unique(summer_hi$nsp1),
                 nsp2 = unique(summer_hi$nsp2))
  }

  # 2 summer, low ghm
  summer_lo <- pair %>%
    filter(season == "summer") %>%
    filter(lo_hi == "Low")

  if(nrow(summer_lo) == 0){
    dhat_sum_lo = tibble(
      measure = NA,
      value = NA,
      season = "summer",
      lo_hi = "Low",
      sp1 = pair_key[[i, "sp1"]],
      sp2 = pair_key[[i, "sp2"]],
      nsp1 = NA,
      nsp2 = NA)

  } else {

    dhat_sum_lo <-
      as_tibble( overlapEst(pull(filter(summer_lo, !is.na(sp1)), sp1),
                            pull(filter(summer_lo, !is.na(sp2)), sp2)),
                 rownames = "measure") %>%
      add_column(season = "summer",
                 lo_hi = "Low",
                 sp1 = pair_key[[i, "sp1"]],
                 sp2 = pair_key[[i, "sp2"]],
                 nsp1 = unique(summer_lo$nsp1),
                 nsp2 = unique(summer_lo$nsp2))
  }

  # 3 winter, hi ghm
  winter_hi <- pair %>%
    filter(season == "winter") %>%
    filter(lo_hi == "High")

  if(nrow(winter_hi) == 0){
    dhat_win_hi <- tibble(
      measure = NA,
      value = NA,
      season = "winter",
      lo_hi = "High",
      sp1 = pair_key[[i, "sp1"]],
      sp2 = pair_key[[i, "sp2"]],
      nsp1 = NA,
      nsp2 = NA)
  } else {

    dhat_win_hi <-
      as_tibble( overlapEst(pull(filter(winter_hi, !is.na(sp1)), sp1),
                            pull(filter(winter_hi, !is.na(sp2)), sp2)),
                 rownames = "measure") %>%
      add_column(season = "winter",
                 lo_hi = "High",
                 sp1 = pair_key[[i, "sp1"]],
                 sp2 = pair_key[[i, "sp2"]],
                 nsp1 = unique(winter_hi$nsp1),
                 nsp2 = unique(winter_hi$nsp2))
  }

  # 4 winter, lo ghm
  winter_lo <- pair %>%
    filter(season == "winter") %>%
    filter(lo_hi == "Low")

  if(nrow(winter_lo) == 0){
    dhat_win_lo <- tibble(
      measure = NA,
      value = NA,
      season = "winter",
      lo_hi = "Low",
      sp1 = pair_key[[i, "sp1"]],
      sp2 = pair_key[[i, "sp2"]],
      nsp1 = NA,
      nsp2 = NA)
  } else {

    dhat_win_lo <-
      as_tibble( overlapEst(pull(filter(winter_lo, !is.na(sp1)), sp1),
                            pull(filter(winter_lo, !is.na(sp2)), sp2)),
                 rownames = "measure") %>%
      add_column(season = "winter",
                 lo_hi = "Low",
                 sp1 = pair_key[[i, "sp1"]],
                 sp2 = pair_key[[i, "sp2"]],
                 nsp1 = unique(winter_lo$nsp1),
                 nsp2 = unique(winter_lo$nsp2))

  }

  res[[i]] <- full_join(dhat_sum_hi, dhat_sum_lo) %>%
    full_join(dhat_win_hi) %>%
    full_join(dhat_win_lo)

  print(paste("Finished", i, "of", nrow(pair_key)))

}

res <- do.call(rbind, res)

# pretty colors
pal <- MetBrewer::MetPalettes$Hiroshige[[1]][c(1, 10)]

# do.call(rbind, res) %>%
res %>% 
  filter(!is.na(measure)) %>% 
  filter(measure == "Dhat4") %>% 
  mutate(pair = paste(sp1, sp2, sep = "-")) %>% 
  group_by(pair, season) %>% 
  filter(nsp1 > 50 & nsp2 > 50) %>%
  # filter(season == "winter") %>% 
  mutate(n = n()) %>% 
  filter(n > 1) %>% 
  dplyr::select(pair, value, season, lo_hi) %>% 
  ungroup() %>% 
  group_by(pair) %>% 
  mutate(avg_value = mean(value)) %>% 
  arrange(avg_value) %>% 
  ungroup() %>% 
  mutate(season = str_to_sentence(season)) %>% 
  mutate(pair = str_replace(pair, "SQUIRRELSANDCHIPMUNKS", "SQUIRRELS"),
         pair = str_replace(pair, "FOXRED", "RED FOX"), 
         pair = str_replace(pair, "SKUNKSTRIPED", "SKUNK"),
         pair = str_replace(pair, "SNOWSHOEHARE", "HARE"),
         pair = str_replace(pair, "CRANESANDHILL", "CRANE")) %>% 
  mutate(pair = tolower(pair)) %>% 
  mutate(pair = factor(pair, levels = unique(pair))) %>% 
  pivot_wider(names_from = lo_hi, values_from = value) %>% 
  mutate(bigger = ifelse(Low > High, "Low", "High")) %>% 
  ggplot()  +
  aes(y = pair) +
  facet_wrap(~season) +
  geom_segment(aes(x = Low, xend = High, 
                   yend = pair, color = bigger)) +
  geom_point(aes(x = High), color = pal[1], size = 3, alpha = 0.7) +
  geom_point(aes(x = Low), color = pal[2], size = 3, alpha = 0.7) +
  scale_color_manual("Disturbance",     
                     values = pal,
                     labels = c("High", 
                                "Low")) +
  theme_minimal() +
  labs(x = "Daily activity pattern overlap") +
  theme(legend.position = "bottom",
        axis.title.y = element_blank()) +
  guides(color = guide_legend(override.aes = list(alpha = 1,
                                                  size = 1)))

# create plot of activity pattern overlap for raccoon-deer
# (figure 4)

# first make a dataframe for annotation
anno <- res %>% 
  filter(!is.na(measure)) %>% 
  filter(measure == "Dhat4") %>% 
  mutate(pair = paste(sp1, sp2, sep = "-")) %>% 
  group_by(pair, season) %>% 
  filter(nsp1 > 50 & nsp2 > 50) %>%
  # filter(season == "winter") %>% 
  mutate(n = n()) %>% 
  filter(n > 1) %>% 
  dplyr::select(pair, value, season, lo_hi) %>% 
  ungroup() %>% 
  group_by(pair) %>% 
  mutate(avg_value = mean(value)) %>% 
  arrange(avg_value) %>% 
  ungroup() %>% 
  mutate(season = str_to_sentence(season)) %>% 
  mutate(pair = str_replace(pair, "SQUIRRELSANDCHIPMUNKS", "SQUIRRELS"),
         pair = str_replace(pair, "FOXRED", "RED FOX"), 
         pair = str_replace(pair, "SKUNKSTRIPED", "SKUNK"),
         pair = str_replace(pair, "SNOWSHOEHARE", "HARE"),
         pair = str_replace(pair, "CRANESANDHILL", "CRANE")) %>% 
  mutate(pair = tolower(pair)) %>% 
  filter(pair == "deer-raccoon") %>% 
  add_column(x = 3, 
             y = 0.35) %>% 
  mutate(label = paste("Overlap:", round(value, 2)),
         lo_hi = paste(lo_hi, "Disturbance"),
         lo_hi = factor(lo_hi, 
                        levels = c(
                          "Low Disturbance", 
                          "High Disturbance"))) 

# get data for racoon and deer
rac_deer <- all_thin %>% 
  mutate(season = str_to_sentence(season)) %>% 
  mutate(sp = str_replace(sp, "SQUIRRELSANDCHIPMUNKS", "SQUIRRELS"),
         sp = str_replace(sp, "FOXRED", "RED FOX"), 
         sp = str_replace(sp, "SKUNKSTRIPED", "SKUNK"),
         sp = str_replace(sp, "SNOWSHOEHARE", "HARE"),
         sp = str_replace(sp, "CRANESANDHILL", "CRANE")) %>% 
  mutate(sp = tolower(sp)) %>%
  filter(sp == "deer" | sp == "raccoon") %>% 
  mutate(lo_hi = paste(lo_hi, "Disturbance"),
         lo_hi = factor(lo_hi, 
                        levels = c(
                          "Low Disturbance", 
                          "High Disturbance"
                        ))) 
# plot (figure 4)
ggplot() + 
  geom_density(data = rac_deer,
               aes(x = suntime, 
                   color = sp, 
                   fill = sp),
               alpha = 0.5) + 
  facet_grid(cols = vars(lo_hi),
             rows = vars(season)) +
  theme_minimal() +
  geom_text(data = anno, aes(x = x, y = y, label = label),
            size = 3.25) +
  scale_x_continuous(breaks = c(1.9, 4.2), 
                     labels = c("Sunrise", "Sunset")) +
  scale_color_manual(values = pal) + 
  scale_fill_manual(values = pal) +
  labs(x = "Time of day",
       y = "Activity") +
  theme(axis.text.y = element_blank(),
        strip.text = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 11, color = "black"),
        axis.text = element_text(size = 10, color = "black"), 
        legend.text = element_text(size = 11, color = "black"),
        axis.ticks.x = element_line(color = "black"),
        legend.title = element_blank(), 
        legend.position = "bottom",
        panel.grid = element_blank(),
        legend.margin = margin(0, 0, 0, 0), 
        legend.box.margin = margin(-5, 0, 0,0))
