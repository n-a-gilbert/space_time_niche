# Data and code from Gilbert et al, "Human disturbance compresses the spatiotemporal niche"

Processed data (raw data cannot be shared since location information from Snapshot Wisconsin volunteers is private) is located in the "data" folder. The code folder contains scripts to run models and to reproduce figures.

## code
pair_analysis_v01.R - code to run the pair models
network_density_analysis_v01.R - code to run the network connectance models
network_proportion_analysis_v01.R - code to run the network proportion (of each antagonism level) models
create_networks_v01.R - code demonstrating how networks were created and statistics derived from them
calculate_activity_pattern_overlap_v01.R - code demonstrating how activity pattern overlap was calculated; also contains code for Fig. 4 and Fig. S7
figure_02.R - code to create figure 2
figure_03.R - code to create figure 3
antagonism_proportion_figures_v01.R - code to create plots of model predictions for antagonsim proportions in network - used to create Fig. 5c and Fig. S9
network_connectance_figure_v01.R - code to create plots of model predictions for network connectance - used to create Fig. 5b and Fig. S8
figure_05c_example_network.R - code to create Fig. 5c visualizing an example network

## data
detection_time_DAPS_v01.RData - data used to calculate daily activity pattern overlap. Pre-processed to be thinned so that detections of a species are at least 30 minutes from another (standard practice for such analyses) and detection times converted to suntime to account for seasonal variation in daylength
supplement_pairwise_data_v01.RData - data for pairwise analysis
network_antag_proportions_v02.csv - network stats: proportion of each antagonism classification in each network
network_connectance_v01.csv - network stat: connectance of each network
df_to_make_example_networks_v01.csv - dataframe to create an example network visualization

## results
connectance_model_results - folder with model output of connectance models (5 models for the different temporal scales)
proportion_model_results - folder with model output of antagonism proportion models (15 models; 5 temporal scales x 3 antagonism rankings)
ghm5k_antag_rai_v01.RDS - model object for relative abundance model
ghm5k_antag_v01.RDS - model object for primary model
ghm5k_mean_sd_antag_v01.RDS - model object for heterogeneity model

