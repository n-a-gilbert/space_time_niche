library(readr)
library(brms)
library(here)

setwd(here::here("data"))

# pair-level data
# column definitions:
# site = ID key for the camera location
# pair = the two species forming the pair
# antag = ranked p(lethal encounter) between members of pair; 1 = low, 2 = medium, 3 = high
# period = which 3-month sampling period
# season = whether data is from Dec-Feb (winter) or Jun-Aug (summer)
# x = average human disturbance within 5-km buffer around camera locations. z-standardized
# x2 = just x squared 
# ttd = time-to-detection - the time (minutes) between successive detections of the two pair species @ a camera

d <- readr::read_csv("space_time_pair_ttd_data_v01.csv")

# this model will take awhile to run
# will depend on machine, but on a lab server, w/ parallel chains, it took ~ 3 days
start <- Sys.time()
mod <- brms::brm(data    = d,
                 family  = Gamma(link = "log"), 
                 formula = ttd ~ 1 + x + mo(antag) + mo(antag):x + ( 1 + x |  pair / season ),
                 prior   = c(prior(normal(7.5, 2), class = Intercept), 
                             prior(normal(0, 0.5), class = b),
                             prior(exponential(1), class = sd), 
                             prior(lkj(2), class = cor),
                             prior(exponential(1), class = shape)), 
                 iter = 4000, warmup = 2000, chains = 3, cores = 3, thin = 2)
end <- Sys.time()
end - start