# Human disturbance compresses the spatiotemporal niche

### [Neil A. Gilbert](https://gilbertecology.com), Jennifer L. Stenglein, [Jonathan N. Pauli](https://pauli.russell.wisc.edu/), & [Benjamin Zuckerberg](https://zuckerberg.russell.wisc.edu/)

#### Please contact the first author for questions about the code or data: Neil A. Gilbert (neil.allen.gilbert@gmail.com)
__________________________________________________________________________________________________________________________________________

## Abstract:  
Human disturbance may fundamentally alter the way that species interact, a prospect that remains poorly understood. We investigated whether anthropogenic landscape modification increases or decreases co-occurrence—a prerequisite for species interactions—within wildlife communities. Using 4 years of data from >2,000 camera traps across a human disturbance gradient in Wisconsin, USA, we considered 74 species pairs (classifying pairs as low-, medium-, or high-antagonism to account for different interaction types) and used the time between successive detections of pairs as a measure of their co-occurrence probability and to define co-occurrence networks. Pairs averaged 6.1 [95% CI: 5.3, 6.8] days between detections in low-disturbance landscapes (e.g., national forests) but 4.1 [3.5, 4.7] days between detections in high-disturbance landscapes, such as those dominated by urbanization or intensive agriculture. Co-occurrence networks showed higher connectance (i.e., a larger proportion of the possible co-occurrences) and greater proportions of low-antagonism pairs in disturbed landscapes. Human-mediated increases in species abundance (possibly via resource subsidies) appeared more important than behavioral mechanisms (e.g., changes in daily activity timing) in driving these patterns of compressed co-occurrence in disturbed landscapes. The spatiotemporal compression of species co-occurrences in disturbed landscapes likely strengthens interactions like competition, predation, and infection unless species can avoid each other at fine spatiotemporal scales. Regardless, human-mediated increases in co-occurrence with—and hence increased exposure to—predators or competitors might elevate stress levels in individual animals, with possible cascading effects across populations, communities, and ecosystems.

## Repository Directory

### [code](./code): Contains code for modeling and generating figures
* pair_analysis_v01.R - code to run the pair models
* network_density_analysis_v01.R - code to run the network connectance models
* network_proportion_analysis_v01.R - code to run the network proportion (of each antagonism level) models
* create_networks_v01.R - code demonstrating how networks were created and statistics derived from them
* calculate_activity_pattern_overlap_v01.R - code demonstrating how activity pattern overlap was calculated; also contains code for Fig. 4 and Fig. S7
* figure_02.R - code to create figure 2
* figure_03.R - code to create figure 3
* antagonism_proportion_figures_v01.R - code to create plots of model predictions for antagonsim proportions in network - used to create Fig. 5c and Fig. S9
* network_connectance_figure_v01.R - code to create plots of model predictions for network connectance - used to create Fig. 5b and Fig. S8
* figure_05c_example_network.R - code to create Fig. 5c visualizing an example network

### [data](./data): Contains data, processed to remove geographic coordinates to protect privacy of volunteers
* detection_time_DAPS_v01.RData - data used to calculate daily activity pattern overlap. Pre-processed to be thinned so that detections of a species are at least 30 minutes from another (standard practice for such analyses) and detection times converted to suntime to account for seasonal variation in daylength
* supplement_pairwise_data_v01.RData - data for pairwise analysis
* network_antag_proportions_v02.csv - network stats: proportion of each antagonism classification in each network
* network_connectance_v01.csv - network stat: connectance of each network
* df_to_make_example_networks_v01.csv - dataframe to create an example network visualization

### [results](./results): Contains model output
* connectance_model_results - folder with model output of connectance models (5 models for the different temporal scales)
* proportion_model_results - folder with model output of antagonism proportion models (15 models; 5 temporal scales x 3 antagonism rankings)
* ghm5k_antag_rai_v01.RDS - model object for relative abundance model
* ghm5k_antag_v01.RDS - model object for primary model
* ghm5k_mean_sd_antag_v01.RDS - model object for heterogeneity model

