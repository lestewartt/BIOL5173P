##### Modelling the Potential Effect of Drug Treatment SARS-CoV-2 Transmission #####

GUID: 2915700s 
Supervisor: Dr. Chris Illingworth 

## This page covers the expected output from each script and how to process them for further analysis ##

#### File #1 -- MATLAB SCRIPT 'Fitting_Models.m' 
Output consists of 5 SVG files: 
     -Plot describing UDCA-induced changes in ACE2 expression over time 
     -Plot describing the fit of the gamma distribution across each defined time period 
     -Lung data fit to a linear model 
     -Bronchi data fit to a linear model 
     -Bottleneck data fit to a negative-binomial model 

    Relevant parameters to downstream analyses 
    Expression model:
     -'initial_params': table describing the parameters for a gamma distribution to each individual day. Day0 is the first column. 
     -'phat1_2': variable describing the parameters for a gamma distribution fit to day 1 and 2 data 
     -'phat3_5': variable describing the parameters for a gamma distribution fit to day 3, 4 and 5 data 
    Virus Exposure Model:
     -'phatgam': variable describing the parameters for a gamma distribution fit to augmented data reflecting the secondary attack rate 


#### File #2 -- PYTHON SCRIPT 'Model_Integration.py' 
Primary output provides the odds of infection relative to the duration of UDCA treatment: 
    odds_ratio_1_2
    odds_ratio_3_5


#### File #3 -- PYTHON SCRIPT 'Hospital_UDCA_Simulation.py' 
This script's output is as follows: 
    -A new directory, titled 'hospital_infection_sim'
    -Each simulation will produce a folder within this directory that is appropriately tagged according to the simulation number (e.g. 'simulation1')
    -Within each simulation folder there will be:
        -a gephi network plot with the subnetworks of infection clusters 
        -A csv file containing information for each infection cluster in the simulation ('ward_infection_clusters')
        -A master dataframe containing the aggregated patient and healthcare worker data ('master_df')
        -A csv file containing the ward infection cluster data ('infection_network_df')

To visualize results in Gephi:
    -Gephi (free, open-sourced software) must be downloaded and installed
    -Load the gephi plot within the Gephi application (e.g. File>Open>gephi_66)
    Appearance
        -Nodes can be coloured according to infection odds, red indicating highest odds and blue indicating lowest 
        -The size of nodes can be adjusted relative to the degree of interactions they have with other nodes
            -Select ranking --> degree --> range from 10 to 40 
    Layout
        -Force Atlas can be used to group nodes according to interactions 
        -Default settings are preserved except: 
            -Repulsion Strength: 1200 
            -Attraction Strength: 6.0 
            -Select Attraction Distrib.
            -Select Adjust By Sizes 
    Statistics
        -Community Detection --> Modularity 
            Set resolution to 0.5 
    Filters 
        -Attributes --> Partition --> ModularityClass
        -Select desired networks from the Query box directly below
To export Gephi results:
    -Preview 
        -Presets: Default curved 
        -Node Labels:
            -Deselect proportional size 
            -Arial 18 Bold 
        Edges:
            -Thickness: 3.5 
        Edge Arrows:
            -Size: 3.0 
    -Export as SVG

#### File #4 -- PYTHON SCRIPT 'Processing_Sim_Results.py' 
Output of this script is as follows: 
    -New folder titled 'hospital_infection_sim_summary'
    -this folder contains:
        -summary information across all 180 simulations ('summary_results')
        -summary stats for all 180 simulations ('summary_results_stats')
        -summary stats for each set of community infection rates (e.g. 'subset_high_comm_inf_rate_stats')
        -a file indicating the proportion of community infection in each simulation ('comm_inf_stats')
        -an svg file of a histogram plot comparing the distribution of percent reduction in infection odds across each set of community infection rates ('percent_reduction_dist')


