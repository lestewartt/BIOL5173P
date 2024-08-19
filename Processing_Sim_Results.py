'''
GUID: 2915700s
Supervised by: Dr. Chris Illingworth 

Processing Hospital Infection Data 
-This script generates summary information for the data produced for each simulation in the hospital infection simulation script 

To run this script, please update the 'main_directory' variable with the path to the 'hospital_infection_sim' folder

Output of this script is as follows: 
    -New folder titled 'hospital_infection_sim_summary'
    -this folder contains:
        -summary information across all 180 simulations 
        -summary stats for all 180 simulations 
        -summary stats for each set of community infection rates 
        -a file indicating the proportion of community infection in each simulation 
        -an svg file of a histogram plot comparing the distribution of percent reduction in infection odds across each set of community infection rates

NOTE: This script will not run normally without input of all 180 simulations
    -To run a subset of the simulations, comment out code from line 278 onwards 
        '''

###Processing results from all 180 simulations###
import os
from natsort import natsorted
import re
import pandas as pd 
import logging 
import scipy.stats
import matplotlib.pyplot as plt

#Setting up a logger to be used in error handling/providing helpful info for user guidance 
logger = logging.getLogger()
#Logging level
logger.setLevel(logging.ERROR)
#Create stream handler and defining the format
sh = logging.StreamHandler()
#Type of error, the date and time, and the associated message 
sh.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - %(message)s'))
logger.addHandler(sh)

#Feeding the hospital sim folder into this script 
#This is the only user input necessary to run 
main_directory = 'C:/Users/2915700s/OneDrive/UoG Masters/Project_Start/hospital_infection_sim'

try:
    #Creating a new output folder for summary information 
    output_folder = 'hospital_infection_sim_summary'
    os.makedirs(output_folder)
    logger.info(f'Folder {output_folder} sucessfully created. Located in: {os.getcwd()}')
except FileExistsError:
    logger.error('File already exists. Rename folder or continue with existing file.')

#Creating empty lists to store the desired files from each simulation folder
subdir_path_networks = []
subdir_path_infection_clusters = []
def pulling_files(primary_directory):
    '''
    Function that pulls the infection network and ward infection cluster files for downstream analysis

    Parameter
        primary_directory = path to the directory containing the 180 simulation folders 
    Return
        lists containing the infection network and ward infection cluster files for each simulation  
    '''
    #os walk is a common funciton used for parsing through a main directory and several subdirectories 
    #subdir_path contains the actual file path to the subdirectories
    #dirs has the name of the subdirectories
    #files is the name of every file within those subdirectories 
    for subdir, dirs, files in os.walk(primary_directory):
        for file in files:
            #Using a regular expression to match the file name 
            #This first file is for the infection_network_df file 
            #Leveraging the naming conventions used in the simulation script 
            match = re.match(r'[a-z]{9}_[a-z]{7}_[a-z]{2}_\d+\.csv$', file)
            if match: 
                #Joins the path to the subdirectory with the file name to create a path to this exact file
                full_path_networks = os.path.join(subdir, file)
                #Appends to the list outside the loop 
                subdir_path_networks.append(full_path_networks)
            #This is a match to the ward infection clusters file 
            match2 = re.search(r'([a-z]{4}_[a-z]{9}_[a-z]{8}_)((\d+))(.csv$)', file)
            if match2:
                #Creating a path to the desired file from this simulation folder
                full_path_infection_clusters = os.path.join(subdir, file)
                #Appending to the list outside the loop 
                subdir_path_infection_clusters.append(full_path_infection_clusters)
            #Using this for calculating the number of patient zero cases are coming in from the community 

#Calling the above function
pulling_files(main_directory)

#Sorting the files in these lists so that they are no longer in lexicographic order 
sorted_subdir_network_data = natsorted(subdir_path_networks)
sorted_infection_cluster_data = natsorted(subdir_path_infection_clusters)
#Creating empty lists to store the paths to the simulation files that correspond to one of the three sets of conditions 
sim_set_1 = []
sim_set_2 = [] 
sim_set_3 = []
def sorting_simulations(subdir_network_data):
    '''
    Function that sorts the infection network files according to the grouping of community infection rate applied 

    Parameter
        subdir_network_data = sorted file paths to the files containing information on the infection networks in each simulation
    Return
        sim_set_1, sim_set_2, sim_set_3 = lists of the file paths grouped by rate of community infection 
    '''
    #Iterating over the list of infecton network files 
    for path in subdir_network_data: 
        #This investigates the results across one simulation == 1 hospital 
        #The regular expression conditions were grouped by parentheses so the 2nd match group, which contains the 
        #simulation number, can be isolated and used for downstream output file naming conventions 
        sim_number = re.split(r'([a-z]{9}_[a-z]{7}_[a-z]{2})_(\d+)\.csv$', path)[2]
        #Must be an integer to be passed through the conditionals below
        sim_number = int(sim_number)
        #Grouping the simulation data based on the 3 sets of unique community infection rates 
        #1 Low community infection rate: 1-20, 61-80, 121-140
        if (1 <= sim_number <= 20) or (61 <= sim_number <= 80) or (121 <= sim_number <= 140):
            sim_set_1.append(path)
        #2 Medium community infection rate: 21-40, 81-100, 141-160
        elif (21 <= sim_number <= 40) or (81 <= sim_number <= 100) or (141 <= sim_number <= 160):
            sim_set_2.append(path)
        #3 High community infection rate: 41-60, 101-120, 161-180 
        elif (41 <= sim_number <= 60) or (101 <= sim_number <= 120) or (161 <= sim_number <= 180):
            sim_set_3.append(path)

#Calling the above function
sorting_simulations(sorted_subdir_network_data)

def processing_infection_clusters(infection_cluster_data):
    '''
    Function that processes the ward infection cluster data for analysis

    Parameter
        infection_cluster_data = path to the file containing the infection cluster data for each simulation 
    Return
        infection_cluster_df = dataframe with the proportion of infection clusters that were initiated by a community infection event for each simulation 
    '''
    #Initiating an empty dataframe to store the community infection stats for each simulation 
    infection_cluster_df = pd.DataFrame(columns=['SimulationNumber', 'PropCommInf'])
    #Iterating over the ward infection cluster files  
    for path in infection_cluster_data:
        #Again, grouping the file name to isolate the simulation number 
        sim_number = re.split(r'([a-z]{4}_[a-z]{9}_[a-z]{8}_)((\d+))(.csv$)', path)[2]
        #Opening the ward infection cluster file as a pandas dataframe 
        ward_cluster_df = pd.read_csv(path, sep='\t')
        #This dataframe contains information indicating if the origin of the infection cluster was from the community or within the hospital
        #Starting an empty tracker to count how many of the infection clusters were initiated by a community infection in each simulation
        community_infections = 0
        #Iterating over the dataframe by row 
        for _, row in ward_cluster_df.iterrows():
            #Pulling the column that indicates where the infection originated
            outbreak_origin = row['OutbreakOrigin']
            #Community infections are marked as coming from the community or incoming from the emergency department, so this conditional checks for these instances 
            if outbreak_origin == 'Comm_' or outbreak_origin == 'ED':
                #If this infection cluster passes the conditional, it is added to the count of community infections 
                community_infections += 1 
        #Calculating the number of infection clusters recorded in this simulation 
        num_networks = len(ward_cluster_df)
        #Calculating the proportion of community infections in this simulation 
        proportion_comm_infections = community_infections / num_networks
        #print(f'{sim_number}: {proportion_comm_infections}')
        #Adding to the dataframe 
        infection_cluster_df.loc[len(infection_cluster_df)] = {'SimulationNumber': sim_number, 'PropCommInf': proportion_comm_infections}

    return infection_cluster_df

#Calling the above function
infection_cluster_df = processing_infection_clusters(sorted_infection_cluster_data)

#Exporting the community infection data to the newly created folder 
output_path = os.path.join(output_folder, 'comm_inf_stats.csv')   
infection_cluster_df.to_csv(output_path, sep='\t')
#print(infection_cluster_df.describe())
    

def summary_subsets(network_file_paths):
    '''
    Function that processes the infection netowrk dataframe for analysis 

    Parameter
        network_file_paths = file paths to the infection network data 
    Return
        infection_network_analysis = pandas dataframe containing summary data for each simulation 
            SimulationNumber = The simulation number the row of data is associated with 
            PercentReduction = % reduction in infection odds across ALL cases (pt and hcw)
            PropCasesReduced = proportion of individuals who were both directly and indirectly impacted by UDCA (pt and hcw)
            InitialInf = cumulative initial infection odds before UDCA intervention 
            AdjustedInf = new cumulative infection odds after UDCA intervention 
            NumHCWs = number of healthcare workers involved in the infection networks for that simulation 
            NumPTs = number of patients involved in the infection networks for that simulation 
            TreatmentEffect = the count of individuals impacted (directly or indirectly) by UDCA 
            PropImpactedPTs = proportion of pts across all infection networks in that simulation who were impacted by UDCA 
            PropImpactedHCWs = proportion of hcws across all infection networks in that simulation who were impacted by UDCA 
            PTPercentReduction = % reduction in infection odds across PATIENT cases 
    '''
    #Initializing an empty dataframe to store this information 
    infection_network_analysis = pd.DataFrame(columns=['SimulationNumber', 'PercentReduction', 'PropCasesReduced', 'InitialInf', 
                                                       'AdjustedInf', 'NumHCWs', 'NumPTs', 'TreatmentEffect', 'PropImpactedPTs', 
                                                       'PropImpactedHCWs', 'PTPercentReduction'])
    for path in network_file_paths: 
        #This investigates the results across one simulation == 1 hospital 
        #The regular expression conditions were grouped by parentheses so the 2nd match group, which contains the 
        #simulation number, can be isolated and used for downstream output file naming conventions 
        sim_number = re.split(r'([a-z]{9}_[a-z]{7}_[a-z]{2})_(\d+)\.csv$', path)[2]
        #Opening the file as a pandas dataframe 
        infection_df = pd.read_csv(path, sep='\t')
        #Creating empty columns to be populated in this dataframe -- this information will be called outside the loop for analysis across the entire dataframe (i.e. the entire simulation) 
        infection_df['NumPTs'] = float()
        infection_df['TreatedPts'] = float()
        infection_df['PropReduced'] = float()

        #Iterating over each infection network in the simulation
        for index, row in infection_df.iterrows():
            #Storing the columns as variables 
            network_size = row['NetworkSize']
            treated_count = row['TreatedCount']
            reduction = row['Reduction']
            hcw_count = row['NumHCWs']
            treatment_impact = row['TreatmentEffect']
            hcw_impacted = row['TreatedHCWs']

            #Cant get the indirect counter to work right so this is my solution (pressed for time)
            indirect_effects = int(treatment_impact) - int(treated_count)
            #Proportion of reduced cases
            prop_reduced = int(treatment_impact) / int(network_size)
            #Getting the number of patients in the network 
            pt_count = int(network_size) - int(hcw_count)
            #Getting the number of patients who were impacted (directly and indirectly) by distribution of UDCA
            pt_impacted = int(treatment_impact) - int(hcw_impacted)

            #Adding the derived information to the previously created columns 
            infection_df.at[index, 'NumPTs'] = pt_count
            infection_df.at[index, 'TreatedPts'] = pt_impacted
            infection_df.at[index, 'PropReduced'] = prop_reduced

        #The total infection odds pre treatment is also synonymous to the total number of individuals from the simulation included 
        total_odds_pre_treatment = infection_df['NetworkSize'].sum()
        total_odds_post_treatment = infection_df['Reduction'].sum()
        difference_in_reduction = total_odds_pre_treatment - total_odds_post_treatment
        #Percent reduction of infection odds across the entire hospital 
        percent_reduction = (difference_in_reduction/total_odds_pre_treatment) * 100
        #Calculating the total number of healthcare workers in the simulation 
        num_hcws = infection_df['NumHCWs'].sum()
        #This is not only the number of patients in the network but also the intial infection odds of all patients 
        num_pts = (infection_df['NetworkSize'].sum()) - (infection_df['NumHCWs'].sum()) 
        #Percent reduction in patients 
        total_pt_odds_post_treatment = infection_df['PtReduction'].sum()
        difference_in_pt_reduction = num_pts - total_pt_odds_post_treatment
        pt_percent_reduction = (difference_in_pt_reduction / num_pts) * 100 
        #Sum of individuals who were impacted either directly or indirectly 
        treatment_effect = infection_df['TreatmentEffect'].sum()
        impacted_patients = infection_df['TreatedPts'].sum()
        impacted_hcws = infection_df['TreatedHCWs'].sum()
        #Prop of impacted patients 
        prop_impacted_patients = impacted_patients / num_pts
        #Prop of impacted healthcare workers 
        prop_impacted_hcws = impacted_hcws / num_hcws
        #Proportion of cases being reduced across the entire hospital 
        prop_reduced = (infection_df['TreatmentEffect'].sum()) / (infection_df['NetworkSize'].sum())
        #Adding the corresponding information for the current simulation to the dataframe 
        infection_network_analysis.loc[len(infection_network_analysis)] = {'SimulationNumber': sim_number, 'PercentReduction': percent_reduction, 
                                                                           'PropCasesReduced': prop_reduced, 'InitialInf': total_odds_pre_treatment, 
                                                                           'AdjustedInf': total_odds_post_treatment, 'NumHCWs': num_hcws, 
                                                                           'NumPTs': num_pts, 'TreatmentEffect': treatment_effect, 
                                                                           'PropImpactedPTs': prop_impacted_patients, 'PropImpactedHCWs': prop_impacted_hcws, 
                                                                           'PTPercentReduction': pt_percent_reduction}
    return infection_network_analysis 

#Calling the above function and exporting the results for all 180 simulations 
all_data = summary_subsets(sorted_subdir_network_data)
output_path = os.path.join(output_folder,'summary_results.csv')   
all_data.to_csv(output_path, sep='\t')

#Calling the above function and storing the results for the descriptive statistics 
all_data_summ = all_data.describe()
output_path = os.path.join(output_folder,'summary_results_stats.csv')   
all_data_summ.to_csv(output_path, sep='\t')

#Calling the above function 
subset_summary1 = summary_subsets(sim_set_1)
summ1 = subset_summary1.describe()
#print(summ1)
#Adding a column to track the 'condition' applied to this set of simulations for downstream processing 
subset_summary1['Condition'] = 1
#Storing the descriptive statistics for the first rate of community infection (low)
output_path = os.path.join(output_folder, 'subset_low_comm_inf_rate_stats.csv')   
summ1.to_csv(output_path, sep='\t')

#Calling the above function 
subset_summary2 = summary_subsets(sim_set_2)
summ2 = subset_summary2.describe()
#print(summ2)
#Adding a column to track the 'condition' applied to this set of simulations for downstream processing 
subset_summary2['Condition'] = 2
#Storing the descriptive statistics for the first rate of community infection (low)
output_path = os.path.join(output_folder, 'subset_medium_comm_inf_rate_stats.csv')   
summ2.to_csv(output_path, sep='\t')

#Calling the above function 
subset_summary3 = summary_subsets(sim_set_3)
summ3 = subset_summary3.describe()
#print(summ3)
#Adding a column to track the 'condition' applied to this set of simulations for downstream processing 
subset_summary3['Condition'] = 3
#Storing the descriptive statistics for the first rate of community infection (low)
output_path = os.path.join(output_folder, 'subset_high_comm_inf_rate_stats.csv')   
summ3.to_csv(output_path, sep='\t')

simulations_of_interest_set = ['6', '16', '66', '76', '126', '136', '26', '36', '86', '96', '146', '156', '46', '56', '106', '116', '166', '176']
parameter_subset = all_data[all_data['SimulationNumber'].isin(simulations_of_interest_set)]
output_path = os.path.join(output_folder, 'subset_parameter_set6.csv')   
parameter_subset.to_csv(output_path, sep='\t')

#Now to compare parameter sets, focused on the set including simulation 66 
simulations_of_interest_set1 = ['6', '16', '66', '76', '126', '136']
#simulations_of_interest_set1 = ['27', '37', '87', '97', '147', '157']
#simulations_of_interest_set1 = ['4', '14', '64', '74', '124', '134']
#simulations_of_interest_set1 = ['1', '11', '61', '71', '121', '131']
parameter_subset_1 = subset_summary1[subset_summary1['SimulationNumber'].isin(simulations_of_interest_set1)]
#print(parameter_subset_1)
param6_subset1_stats = parameter_subset_1.describe()
output_path = os.path.join(output_folder, 'param_set6_low_comm_stats.csv')   
param6_subset1_stats.to_csv(output_path, sep='\t')

#Pulling the data for parameter set 6
simulations_of_interest_set2 = ['26', '36', '86', '96', '146', '156']
#simulations_of_interest_set2 = ['27', '37', '87', '97', '147', '157']
#simulations_of_interest_set2 = ['24', '34', '84', '94', '144', '154']
#simulations_of_interest_set2 = ['21', '31', '81', '91', '141', '151']
parameter_subset_2 = subset_summary2[subset_summary2['SimulationNumber'].isin(simulations_of_interest_set2)]
#print(parameter_subset_2)
#print(parameter_subset_2.describe())
param6_subset2_stats = parameter_subset_2.describe()
output_path = os.path.join(output_folder, 'param_set6_med_comm_stats.csv')   
param6_subset2_stats.to_csv(output_path, sep='\t')

#Pulling the data for parameter set 6
simulations_of_interest_set3 = ['46', '56', '106', '116', '166', '176']
#simulations_of_interest_set3 = ['47', '57', '107', '117', '167', '177']
#simulations_of_interest_set3 = ['44', '54', '104', '114', '164', '174']
#simulations_of_interest_set3 = ['41', '51', '101', '111', '161', '171']
parameter_subset_3 = subset_summary3[subset_summary3['SimulationNumber'].isin(simulations_of_interest_set3)]
#print(parameter_subset_3)
#print(parameter_subset_3.describe())

param6_subset3_stats = parameter_subset_3.describe()
output_path = os.path.join(output_folder, 'param_set6_high_comm_stats.csv')   
param6_subset3_stats.to_csv(output_path, sep='\t')


all_dfs = [subset_summary1, subset_summary2, subset_summary3]
#Making sure all the columns are the same
for df in all_dfs:
    df.columns = ['SimulationNumber', 'PercentReduction', 'PropCasesReduced', 'InitialInf', 'AdjustedInf', 'NumHCWs', 'NumPTs', 
                  'TreatmentEffect', 'PropImpactedPTs', 'PropImpactedHCWs', 'PTPercentReduction', 'Condition']
#Merging these dataframes with the newly added 'condition' column
merged = pd.concat(all_dfs, ignore_index=True)
#print(merged.groupby("Condition")['PTPercentReduction'].describe())

#print(merged.groupby("Condition")['PropImpactedPTs'].describe())

## NOTE: ONE WAY ANOVA TO COMPARE THE MEAN % REDUCTION ACROSS THE 3 CONDITIONS
#Shapiro-Wilks test for normality in each subset 
normality_p_red_1 = scipy.stats.shapiro(subset_summary1['PTPercentReduction'])
#print(normality_p_red_1)
#ShapiroResult(statistic=0.9456226461013568, pvalue=0.009744750633332614)
#non-normal dist
normality_p_red_2 = scipy.stats.shapiro(subset_summary2['PTPercentReduction'])
#print(normality_p_red_2)
#ShapiroResult(statistic=0.9731736070270203, pvalue=0.2080906291925918)
#normal dist
normality_p_red_3 = scipy.stats.shapiro(subset_summary3['PTPercentReduction'])
#print(normality_p_red_3)
#ShapiroResult(statistic=0.9892472825563939, pvalue=0.8764704763208155)
#normal dist 

#Plotting the above data in a histogram to visually assess the distribution 
#alpha makes them partially transparent 
#blue 
plt.hist(subset_summary1['PTPercentReduction'], alpha=0.5, label='low rate', color='#495af2', edgecolor='#323ea8')
#red
plt.hist(subset_summary3['PTPercentReduction'], alpha=0.5, label='high rate', color='#d64545', edgecolor ='#ab1616')
#green
plt.hist(subset_summary2['PTPercentReduction'], alpha=0.5, label='medium rate', color='#338a36', edgecolor='#044d07')
#Adding a legend 
plt.legend(loc='upper right')
#Adding a title
plt.title('Percent Reduction Distributions Across Community Infection Rates')
plt.xlabel('Percent Reduction in Patient Infection Odds')
plt.ylabel('Frequency')
output_path = os.path.join(output_folder, 'percent_reduction_dist.svg')   
plt.savefig(output_path)
#plt.show()

#Testing for the variance
levene_percent_reduction = scipy.stats.levene(subset_summary1['PTPercentReduction'], subset_summary2['PTPercentReduction'], subset_summary3['PTPercentReduction'])
#print(levene_percent_reduction)
#LeveneResult(statistic=0.1327043343464391, pvalue=0.8758110219542314)
#p value > 0.05 == variances are not significantly different

#Going to test as if the data is normally distributed and homoskedastic --> one way ANOVA 
#This tests the means of two or more independent groups 
ANOVA_percent_reduction = scipy.stats.f_oneway(subset_summary1['PTPercentReduction'], subset_summary2['PTPercentReduction'], subset_summary3['PTPercentReduction'])
#print(ANOVA_percent_reduction)
#F_onewayResult(statistic=0.18125359512689757, pvalue=0.8343784198095198)
#p value is very high, indicating failure to reject the null hypothesis 
#--> no statistically significant difference between the mean percent reduction of the 3 groups 

#To account for the one shapiro result for condition 2 stating its non normal, we will also use a non parametric test to compare with the above results 
kruskal_percent_reduction = scipy.stats.kruskal(subset_summary1['PTPercentReduction'], subset_summary2['PTPercentReduction'], subset_summary3['PTPercentReduction'])
#print(kruskal_percent_reduction)
#KruskalResult(statistic=0.7954450583179096, pvalue=0.6718484201477672)
#kruskal results substantiate that there is no significant difference 




