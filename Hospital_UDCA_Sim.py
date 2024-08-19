'''
GUID: 2915700s
Supervised by: Dr. Chris Illingworth 

Hospital Simulation Infection Tracking Code
This script accepts the 180 simulations from the provided hospital simulation data (Evans et al., 2024)

NOTE: To run this script, please update the 'pt_folder' and 'hcw_folder' file paths (lines 41 & 43) to the appropriate directory where this data is stored.
This script's output is as follows: 
    -A new directory, titled 'hospital_infection_sim
    -Each simulation will produce a folder within this directory that is appropriately tagged according to the simulation number
    -Within each simulation folder there will be a gephi network plot with the subnetworks of infection clusters
    -A csv file containing information for each infection cluster in the simulation
    -A master dataframe containing the aggregated patient and healthcare worker data 
    -A csv file containing the ward infection cluster data 
Processing of all 180 simulations takes approximately 20 minutes 
    
Comment on University of Glasgow's AI Guidance and Policy:

Artificial intelligence should be used appropriately and responsibly, and it should not be used as a replacement for independent thinking.
All use of ChatGPT was as an information tool and not for directly generating code. ChatGPT was used to assist in debugging, highlighted by **

Statement of Academic Integrity: This submission is a product of this student's work and approach to solving the designated task.  
'''

#Importing necessary libraries 
import pandas as pd 
import networkx as nx 
import logging
import os
from natsort import natsort

#Setting up a logger to be used in error handling/providing helpful info for user guidance 
logger = logging.getLogger()
#Logging level
logger.setLevel(logging.WARNING)
#Create stream handler and defining the format
sh = logging.StreamHandler()
#Type of error, the date and time, and the associated message 
sh.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - %(message)s'))
logger.addHandler(sh)

#File path to patient simulation data folder
pt_folder = 'C:/Users/2915700s/OneDrive/UoG Masters/Project_Start/patient_sim_data'
#File path to healthcare worker simulation data folder
hcw_folder = 'C:/Users/2915700s/OneDrive/UoG Masters/Project_Start/hcw_sim_data'

#List of file names in the directory 
pt_files = os.listdir(pt_folder)
#Ensures the files are not sorted according to python's nature, which is lexicographic order
pt_files = natsort.natsorted(pt_files)
#List of file names in the directory 
hcw_files = os.listdir(hcw_folder)
#Ensures the files are not sorted according to python's nature, which is lexicographic order
hcw_files = natsort.natsorted(hcw_files)

#Ensuring the pt and hcw folders contain the same number of simulations to be correctly paired together 
try:
    len(pt_files) == len(hcw_files)
except:
    logger.warning(f'Discrepancy in paired data. Please check contents of simulated data folders.')

#Ensures the files are processed in pairs 
paired_sim_data = zip(pt_files, hcw_files)

'''NOTE: TRACKING INFECTION CLUSTERS'''

def add_to_ward_detected_infections(ward_outbreaks, ward, outbreak_start, last_infection_date, outbreak_indvidual, outbreak_origin):
    '''
    This is a function called within the 'detected_infection_events' function, 
    used to add infection cluster information to the dataframe

    Parameters
        ward_outbreaks = dataframe containing information on each recorded infection cluster 
        ward = the ward in which this infection cluster originated
        outbreak_start = the date of the detected infection event 
        last_infection_date = the last recorded infection event in this cluster 
        outbreak_individual = the individual who was detected
        outbreak_origin = indicates how the outbreak_individual became infected 
    Return
        An updated row to the dataframe tracking information on infection clusters in the hospital
    '''
    #Calculates the number of days the infection network lasted
    days_since_outbreak = (last_infection_date - outbreak_start).days
    #Returns the above as '__ days'
    outbreak_duration = pd.Timedelta(days=days_since_outbreak)
    #Adds the new row to the dataframe
    ward_outbreaks.loc[len(ward_outbreaks)] = {'WardNumber': ward, 
                                               'Start_Date': outbreak_start, 
                                               'End_Date': last_infection_date, 
                                               'OutbreakDuration': outbreak_duration, 
                                               'OutbreakIndividual':outbreak_indvidual, 
                                               'OutbreakOrigin': outbreak_origin}

 
def detected_infection_events(infection_df, wards_of_interest):
    '''
    Function to track and record outbreaks across wards with detected infections 

    Parameters
        infection_df = dataframe containing patient data from the hospital simulation 
        wards_of_interest = list of the wards in which a detected infection case was recorded 
    Return
        A dataframe storing all the above parameters to the corresponding infection cluster 
    '''
    #Initializing a dataframe to store the outbreak information 
    ward_outbreaks = pd.DataFrame(columns=['WardNumber', 'Start_Date', 'End_Date', 'OutbreakDuration', 
                                           'OutbreakIndividual', 'OutbreakOrigin'])
    #Iterating through the wards that have detected cases 
    for ward in wards_of_interest: 
        #The ward of interest is then used to filter the sorted patient dataframe made on line 105 
        ward_data = infection_df[infection_df['WardNumber'] == ward]
        #Initializing variables to be used to track relevant infection cluster information 
        #Stored in this part of the loop to account for the first ward in the list, as well as transitioning to a new ward
        #This variable stores the previous infection date passed through the loop
        last_infection_date = None
        #Stores the date of the detected infection 
        outbreak_start = None
        #Stores the ID of the individual who was detected 
        outbreak_individual = None 
        #Stores information indicating how the detected individual became infected 
        outbreak_origin = None 
        #Iterating through the list of the patients who were infected on this ward where an infection was detected 
        for _, line in ward_data.iterrows():
            #Establishing the current infection data 
            infection_date = line['Date']
            #Storing the detection status of this case 
            detection_status = line['Detected']
            #Individual being considered
            individual = line['ID']
            #How the individual was infected (i.e. patient, healthcare worker, community)
            infector = line['InfectedBy']
            #If there has been no outbreak initiated yet and this infection was detected, then this is a new outbreak 
            if detection_status == True and outbreak_start is None: 
                #Setting the start of the outbreak to the current date of infection 
                #Updating the tracking variables 
                outbreak_start = infection_date 
                last_infection_date = infection_date
                outbreak_individual = individual 
                outbreak_origin = infector
            #If this current iteration is in the middle of an ongoing outbreak, then we now determine how long it has been since the last infection on that ward
            #Conditions for a new outbreak: must exceed 14 days since the last infection and must also be detected 
            elif outbreak_start is not None and (infection_date - last_infection_date > pd.Timedelta(days=14)): 
                #If it passes the above conditions and was detected, its a new outbreak 
                if detection_status == True:
                    #Storing the previous outbreak information 
                    add_to_ward_detected_infections(ward_outbreaks, ward, outbreak_start, last_infection_date, outbreak_individual, outbreak_origin)
                    #Updating the outbreak tracking variables to reflect the new outbreak period
                    outbreak_start = infection_date 
                    last_infection_date = infection_date 
                    outbreak_individual = individual 
                #If it has been more than 14 days since the last infection, but this new infection was not detected
                #Then it is not considered a new outbreak, but the previous outbreak is over 
                if detection_status == False:
                    #Closing the previous outbreak 
                    add_to_ward_detected_infections(ward_outbreaks, ward, outbreak_start, last_infection_date, outbreak_individual, outbreak_origin)
                    #Resetting the outbreak status to indicate no current outbreak 
                    outbreak_start = None 
                    outbreak_individual = None 
                    outbreak_origin = None
                    #Updating the last infection date to the current date 
                    last_infection_date = infection_date 
            #Ensures the last infection date is updated**
            #Otherwise all recorded days are 0 
            else:  
                last_infection_date = infection_date 
        #This handles the final line in the current ward if there was an ongoing outbreak
        if outbreak_start is not None:
           add_to_ward_detected_infections(ward_outbreaks, ward, outbreak_start, last_infection_date, outbreak_individual, outbreak_origin)
    #Returns the ward outbreak tracking file
    return ward_outbreaks #to_csv('ward_infection_clusters', sep='\t')

''' NOTE: BUILDING THE NETWORK PLOT'''

def build_network(dataframe): 
    '''
    This function handles building the entire network 

    Parameter
        dataframe = master dataframe containing the aggregated patient and healthcare worker dataframes 
    Return
        Directed Networkx graph that contained node and edge information connecting the infector to the infected individual
    '''
    #Initializing an empty graph object 
    #Using Di graph because edges should be directional 
    G = nx.DiGraph()
    #Iterating over each row in the master dataframe 
    for _, row in dataframe.iterrows():
        #Naming variables 
        infected = row['ID']
        infector = row['InfectedBy']
        #Solution for repeated 'Comm_' entries denoting a case of community infection: by making the HCW name a suffix so it will be a unique node 
        #If the infection came from the emergency department, then the patient's name will be appended to the end so its a unique node
        #Modify infector if necessary
        if infector == 'Comm_' or infector == 'ED':
            infector = infector + str(infected)
        #Ensuring the nodes are unique 
        if infected not in G.nodes():
            G.add_node(infected)
        if infector not in G.nodes():
            G.add_node(infector)
        #Adding edges to connect the infector with the infected
        G.add_edge(infector, infected, date=str(row['Date']))
    return G 

'''NOTE: UPDATING INFECTION CLUSTER NETWORKS WITH ADJUSTED INFECTION ODDS ACCORDING TO UDCA TREATMENT'''

def build_dfs_tree(outbreak_origin_list, G):
    '''
    Helper function to make the depth first search trees

    Parameters
        outbreak_origin_list = list of all the unique ID's associated with the start of an infection cluster (i.e. the origin/source)
        G = networks of infection 
    Return 
        A dictionary containing all the descending nodes associated to the source node (stored as the key to the dict)
    '''
    #Initializing empty dict
    outbreak_tree = {}
    #Starts the search at the point of origin (root node), which in this case is the detected patient
    #It traverses each branch as far as possible before backtracking 
    for patient_zero in outbreak_origin_list:
        #Calling the dfs function to trace the infection networks relative to the root node 
        dfs_tree = nx.dfs_tree(G, patient_zero)
        #Storing as a dict with the root node as the key 
        outbreak_tree[patient_zero] = dfs_tree
    return outbreak_tree


def infection_odds(master_df, outbreak_origin_list, G):
    '''
    Function that applies the UDCA treatment regime to the infection clusters identified from the simulated hospital data

    Parameters
        master_df = the master dataframe containing aggregated patient and healthcare worker data 
        outbreak_origin_list = list of all the unique ID's associated with the start of an infection cluster (i.e. the origin/source)
        G = networks of infection
    Return 
        A comprehensive dataframe containing information relative to each infection cluster, such as: 
            -Source: the source node's unique ID 
            -NetworkSize: the number of nodes in the infection network (also doubles as the initial cumulative infection odds)
            -TreatedCount: the number of patients treated 
            -Reduction: the adjusted infection odds after targeted UDCA treatment 
            -NumHCWs: the number of healthcare workers in the network 
            -TreatmentEffect: the total count of individuals impacted both directly and indirectly by UDCA 
            -TreatedHCWs: the number of healthcare workers who were impacted by UDCA (always indirect)
            -PtReduction: reduction in patient infection odds 
        The DFS tree dictionary (created by the helper function above)
    '''
    #Initializing a dict to store the infection odds for each node that is iterated over 
    adjusted_infection_odds = {}
    #Creating an empty dataframe to store all the associated information for each infection network 
    infection_odds_table = pd.DataFrame(columns=['Source', 'NetworkSize', 'TreatedCount', 'Reduction', 'NumHCWs', 
                                                 'TreatmentEffect', 'TreatedHCWs', 'PtReduction'])
    #Calling the dfs_tree function to pull the nodes that are descended from each source node 
    outbreak_tree = build_dfs_tree(outbreak_origin_list, G)
    #Iterating through the source nodes and their corresponding infection networks 
    for source_node, infection_network in outbreak_tree.items():
        #Initializing variables to indicate origin of outbreak 
        outbreak_origin = False
        outbreak_start = None
        ward_origin = None 
        #This will store the reduction calculation 
        infection_sum = 0 
        #Calculating the adjusted infection odds for patients 
        pt_infection_sum = 0 
        #Tracking the number of nodes in the network 
        network_length = 0
        #Number of patients treated 
        treated = 0 
        #Count of healthcare workers in the network 
        hcw_count = 0
        #Number of individuals who are indirectly impacted --> doesnt work it tracks both direct and indirect instances..
        indirect_effects = 0 
        #Healthcare workers who were indirectly impacted 
        hcw_treated = 0 
        #The source node will always have an odds of 1 so hard coding this is not a problem 
        infection_source_odds = float(1.0)
        #Iterating over each node corresponding to the source node 
        for node in infection_network.nodes():
            #Filtering the master dataframe for the node's matching information 
            filtered_df = master_df[(master_df['ID'] == node)]
            #Incrementally adding to the network length 
            network_length += 1 
            #Can now access the information related to the node to determine their infection odds
            for index, row in filtered_df.iterrows():
                #Storing variables from the dataframe 
                infector = row['InfectedBy']
                infection_date = row['Date']
                ward_number = row['WardNumber']
                role = row['Role']
                current_inf_odds = row['InfectionOdds']
                admission = row['Admission']
                discharge = row['Discharge']
                treatment_status = row['TreatmentStatus']
                #Checking if the node is the source (because of how dfs_tree works, the first iteration will always be the source node)
                if node == source_node: 
                    #Setting the appropriate variables to correspond to the source node's data for accurate downstream calculations/tracking
                    ward_origin = row['WardNumber']
                    outbreak_start = infection_date
                    #Updating the boolean 
                    outbreak_origin = True
                    #Adding the source node information to the dict so it can be referred back to 
                    adjusted_infection_odds[source_node] = infection_source_odds
                    #Adjusting the cumulative infection odds 
                    infection_sum += infection_source_odds
                    #The source node should always be a patient but including this conditional just to be certain 
                    if role == 'Patient':
                        #Adding to the patient infection odds counter 
                        pt_infection_sum += infection_source_odds
                    continue 

                #Passing a conditional to see who is treated with UDCA 
                #Not administering the drug to the outbreak patient, they must be on the same ward as the origin, and they must also be a patient  
                if (outbreak_origin == True) and (ward_number == ward_origin) and (role == 'Patient'):
                    #Determining how many days have passed between the detected infection and this individual's case of infection
                    days_since_outbreak = (infection_date - outbreak_start).days
                    #Determining how many days have passed between the date of detected infection and the date of this individual's admission to the hospital
                    #This is relevant for tracking those who were admitted after the date of detection 
                    admission_post_outbreak = (admission - outbreak_start).days
                    #If the detected infection falls between this patient's admission and discharge date: 
                    if (admission <= outbreak_start) and (discharge >= outbreak_start):
                        #Updating their infection odds to reflect how long they have been treated with UDCA relative to when they would have become infected originally
                        if days_since_outbreak == 1 or days_since_outbreak == 2:
                            #Using .at instead of .loc to update one row at one specific column 
                            updated_odds = 0.639
                            #Updating the count of treated patients 
                            treated += 1 
                            #Storing the adjusted infection odds for downstream calculations 
                            adjusted_infection_odds[node] = updated_odds
                            #Updating the variable indicating this patient's infection odds 
                            current_inf_odds = updated_odds
                            #Updating the master dataframe with the adjusted treatment information 
                            master_df.at[index, 'TreatmentStatus'] = True
                            master_df.at[index, 'InfectionOdds'] = updated_odds
                        #if patients were infected more than 2 days after the detected infection (i.e. treated with UDCA for more than 2 days)
                        #their infection odds are reduced to reflect a longer duration of treatment than the above condition
                        elif days_since_outbreak > 2: 
                            updated_odds = 0.496
                            #Updating the count of treated patients 
                            treated += 1 
                            #Storing the adjusted infection odds for downstream calculations 
                            adjusted_infection_odds[node] = updated_odds
                            #Updating the variable indicating this patient's infection odds 
                            current_inf_odds = updated_odds
                            #Updating the master dataframe with the adjusted treatment information 
                            master_df.at[index, 'TreatmentStatus'] = True
                            master_df.at[index, 'InfectionOdds'] = updated_odds

                    #Applying treatment to patients who were admitted to the ward within 5 days of the detected infection
                    elif admission_post_outbreak <= 5:
                        #Determining the number of days they would have been treated with UDCA relative to their infection date 
                        days_since_admission = (infection_date - admission).days
                        if days_since_admission == 1 or days_since_admission == 2:
                            #If they would have become infected within 2 days 
                            updated_odds = 0.639
                            #Updating the count of treated patients 
                            treated += 1 
                            #Storing the adjusted infection odds for downstream calculations 
                            adjusted_infection_odds[node] = updated_odds
                            #Updating the variable indicating this patient's infection odds
                            current_inf_odds = updated_odds
                            #Updating the master dataframe with the adjusted treatment information 
                            master_df.at[index, 'TreatmentStatus'] = True
                            master_df.at[index, 'InfectionOdds'] = updated_odds
                        elif days_since_admission > 2: 
                            #If they would have become infected between 3 and 5 days 
                            updated_odds = 0.496
                            #Updating the count of treated patients 
                            treated += 1 
                            #Storing the adjusted infection odds for downstream calculations 
                            adjusted_infection_odds[node] = updated_odds
                            #Updating the variable indicating this patient's infection odds
                            current_inf_odds = updated_odds
                            #Updating the master dataframe with the adjusted treatment information 
                            master_df.at[index, 'TreatmentStatus'] = True
                            master_df.at[index, 'InfectionOdds'] = updated_odds
                       
                #Some healthcare workers slipped through the cracks --> didnt pass the above check but still listed as being treated... 
                #At this point, all healthcare workers should have an infection odds of 1 because none of them are being treated with UDCA.
                #Their infection odds will not change until the downstream effects of UDCA are implemented (which is conducted below)
                if role == 'HCW':
                    current_inf_odds = 1
                    hcw_count += 1 
                #If the individual was infected by the source node (this should always be the case for at least the second node in the list)
                if infector == source_node: 
                    #Updating their infection odds relative to the parent node 
                    new_infection_odds = float(infection_source_odds) * float(current_inf_odds)
                    #Adding this information to the dict so it can be referred back to 
                    adjusted_infection_odds[node] = new_infection_odds
                    #Updating the original dataframe (this is necesary for exporting the trees to gephi)
                    master_df.at[index, 'InfectionOdds'] = new_infection_odds
                    #Updating the infection sum variable 
                    infection_sum += new_infection_odds
                    #Updating the adjusted patient odds 
                    if role == 'Patient': 
                        pt_infection_sum += new_infection_odds
                    #If the individual was not treated by UDCA but they still have reduced infection odds 
                    if treatment_status == False and new_infection_odds != 1:
                        #Track as an indirect effect 
                        indirect_effects += 1 
                        #Tracking the number of healthcare workers that were impacted  
                        if role == 'HCW':
                            hcw_treated += 1
                    else:
                        pass
                else:
                    #Pulling the infection odds of the parent node from the dict using .get()
                    parent_inf_odds = adjusted_infection_odds.get(infector)
                    try:
                        #Calculating the current node's infection odds relative to the parent 
                        new_infection_odds = float(parent_inf_odds) * float(current_inf_odds)
                    except TypeError as e:
                        logger.warning(f'failed to compute infection odds for {node}: {e}')
                    #Storing in the dict for later referencing 
                    adjusted_infection_odds[node] = new_infection_odds
                    #Updating the dataframe and infection sum variable 
                    master_df.at[index, 'InfectionOdds'] = new_infection_odds
                    infection_sum += new_infection_odds
                    #Tracking the adjusted patient infection odds 
                    if role == 'Patient': 
                        pt_infection_sum += new_infection_odds
                    #If the individual was not directly administered UDCA but their likelihood of infection was reduced
                    if treatment_status == False and new_infection_odds != 1:
                        indirect_effects += 1 
                        #Tracking the number of healthcare workers that were impacted 
                        if role == 'HCW':
                            hcw_treated += 1 
                    else: 
                        pass
        #Updating the total reduced infection odds once all the nodes branching from this source have been visited
        infection_odds_table.loc[len(infection_odds_table)] = {'Source': source_node, 'NetworkSize': network_length, 
                                                               'TreatedCount': treated, 'Reduction': infection_sum, 
                                                               'NumHCWs': hcw_count, 'TreatmentEffect': indirect_effects, 
                                                               'TreatedHCWs': hcw_treated, 'PtReduction': pt_infection_sum}

    return infection_odds_table, outbreak_tree


'''NOTE: Body of the code continued -- just calling functions and exporting things to gephi'''

def parallel_processing_files(paired_files, pt_folder, hcw_folder):
    '''
    Primary function of the script that nests the above functions to track infection clusters, create infection networks, and apply UDCA treatment

    Parameters
        paired_files = the zipped patient and healthcare file names paired by simulation number 
        pt_folder = path to the folder containing the patient simulated data
        hcw_folder = path to the folder containing the healthcare worker simulated data 
    Return 
        master data = contains the aggregated patient and healthcare worker information 
        ward infection cluster data = contains information stored from the infection cluster function 
        infection network data = contains information stored from the infection odds function 
    '''
    #Serves as a naming convention for output files 
    simulation_count = 0 
    #Iterating over the grouped patient and healthcare worker file names 
    for file_name1, file_name2 in paired_files:
        #Joining the filepath with the name of the files in the loop 
        file_1_path = os.path.join(pt_folder, file_name1)
        file_2_path = os.path.join(hcw_folder, file_name2)
        try:
            #Reading in as a pandas dataframe and uncompressing the files 
            pt_df = pd.read_csv(file_1_path, compression='gzip')
            hcw_df = pd.read_csv(file_2_path, compression='gzip')
            #Incrementally adding to the sim count as files are processed 
            simulation_count += 1 
        except FileNotFoundError: 
            logger.error(f'Error loading file. Please check {file_1_path} and {file_2_path}')
            raise SystemExit(1)

        ''' NOTE: Pre-Processing Data '''

        #Adding explicit columns that designate the individual's 'role' in the simulation 
        pt_df['Role'] = 'Patient'
        hcw_df['Role'] = 'HCW'
        #Adapting the ID's to solve the problem of duplicate IDs that was found between patients and healthcare workers...
        pt_df['ID'] = 'Patient_' + pt_df['ID'].astype(str)
        hcw_df['ID'] = 'HCW_' + hcw_df['ID'].astype(str)
        #Converting to datetime format (important for downstream calculations)
        pt_df['Date'] = pd.to_datetime(pt_df['Date'])
        hcw_df['Date'] = pd.to_datetime(hcw_df['Date'])

        #Selecting the earliest and latest recorded date of the patients' hospital visit
        pt_df_hostpital_stay = pt_df.groupby('ID').agg(Admission=('Date','min'), Discharge=('Date','max')).reset_index()
        #Filtering for instances in which there was an infection 
        infection_statuses = ['INFECTED_A', 'INFECTED_S', 'INFECTED_P']
        pt_df = pt_df[pt_df['Infection.Status.CT'].isin(infection_statuses)]
        hcw_df = hcw_df[hcw_df['Infection.Status.CT'].isin(infection_statuses)]
        #Setting baselines -- none have been treated yet and the table now only contains infection events 
        pt_df['InfectionOdds'] = float(1)
        pt_df['TreatmentStatus'] = False
        hcw_df['InfectionOdds'] = float(1)
        hcw_df['TreatmentStatus'] = False

        #Patient data frame is a bit different from the healthcare worker dataframe, where if it was a community infection it was just left blank
        #Using .loc to select rows where community infection is set to true and the infectedby/person columns are either blank or set to NA.
        #Any rows selected will be replaced with 'Comm_' or 'Comm' to match the HCW dataframe 
        pt_df.loc[(pt_df['Community.Infection'] == True) & (pt_df['InfectedBy'].isna() | (pt_df['InfectedBy'] == '')),'InfectedBy'] = 'Comm_'
        pt_df.loc[(pt_df['Community.Infection'] == True) & (pt_df['Person'].isna() | (pt_df['Person'] == '')),'Person'] = 'Comm'


        ''' NOTE: Below Code Pertains to Infection Clusters in Wards '''
        #Filtering to find detected cases of infection --> note: the pt_df now only contains recorded infections 
        #Using the pt_df to not lose detection events from aggregating the pt information 
        detected_master_df = pt_df[pt_df['Detected'] == True]
        #Finding the earliest detected case in each ward 
        earliest_ward_detections = detected_master_df.loc[detected_master_df.groupby('WardNumber')['Date'].idxmin()]
        #Earliest_ward_detections.to_csv('earliest_ward_detections', sep='\t')
        #Filtering for the wards that have detected infections 
        wards_of_interest = earliest_ward_detections['WardNumber'].unique()
        #Only patients are given ward assignments (and they will be the individuals given UDCA, so using the patient dataframe for tracking)

        #Original method of filtering to make the IDs unique created a problem with duplicate HCW IDs. 
        #Aggregating the tables solved this problem and also creates a more flexible way of viewing the data across different entries for one individual,
        #which can be done easily below by adjusting the data that is stored upon aggregation. 

        #Sorting dataframe by ID and Date 
        pt_aggregated = pt_df.sort_values(by=['ID', 'Date'])
        #Filtering for the columns I want to keep for downstream analysis 
        patient_aggregated = pt_aggregated[['ID', 'Date', 'Infection.Status.CT', 'Community.Infection', 'WardNumber', 'InfectedBy', 
                                            'Detected', 'InfectionOdds', 'TreatmentStatus', 'Role']]
        #Grouping by ID to make it a unique value, and then adding the information for the other columns according to my needs --> have to convert numerical values to strings (i.e. infection odds)
        #Note: originally filtered for unique wards and saw that some patients were moved to different wards: 'WardNumber':lambda x: ', '.join(map(str, x.unique()))
        pt_aggregated = pt_aggregated.groupby('ID').agg({'Date':lambda x: x.min(), 'Infection.Status.CT':lambda x: ', '.join(x.unique()), 
                                                         'Community.Infection':'first', 'WardNumber':'first', 'InfectedBy': 'first', 
                                                         'Detected': 'first', 'InfectionOdds':lambda x: ', '.join(map(str, x.unique())), 
                                                         'TreatmentStatus': 'first', 'Role': 'first'})
        #Its important to reset the index after using the groupby and merge functions to maintain continuity 
        pt_aggregated = pt_aggregated.reset_index()
        #Adding the columns that mark the beginning and end of the patient's hospital stay 
        pt_aggregated = pt_aggregated.merge(pt_df_hostpital_stay[['ID', 'Admission', 'Discharge']], on='ID', how='left').reset_index()

        #Same order of operations as above but adapted to the hcw df as they are slightly different in content
        hcw_aggregated = hcw_df.sort_values(by=['ID', 'Date'])
        hcw_aggregated = hcw_aggregated[['ID', 'Date', 'Infection.Status.CT', 'Community.Infection', 'InfectedBy', 'Detected', 'InfectionOdds', 
                                         'TreatmentStatus', 'Role']]
        hcw_aggregated = hcw_aggregated.groupby('ID').agg({'Date':lambda x: x.min(), 'Infection.Status.CT':lambda x: ', '.join(x.unique()), 
                                                           'Community.Infection':'first', 'InfectedBy': 'first', 'Detected': 'first', 
                                                           'InfectionOdds':lambda x: ', '.join(map(str,x.unique())), 'TreatmentStatus': 'first', 
                                                           'Role': 'first'})
        hcw_aggregated = hcw_aggregated.reset_index()

        #Calling the above functions 
        outbreak_wards = detected_infection_events(pt_df, wards_of_interest)
        #Pulling the list of IDs of the individuals linked to initiating an infection cluster 
        outbreak_starter = outbreak_wards['OutbreakIndividual'].tolist()

        #Ensuring the index's are correct because its creating problems when updating the infection odds based on index      
        merged_aggregated = pd.concat([pt_aggregated, hcw_aggregated])
        #The indices become a problem after concatenating the two dataframes, and using the unique identifier (the ID) ensures the index is unique for each row
        #Drop set to false ensures the original ID column is not removed (critical for downstream processes)
        merged_aggregated = merged_aggregated.set_index('ID', drop=False)

        #calling the function to build the networks of infection 
        G = build_network(merged_aggregated)
        #Calling the function that creates the infection cluster dataframe and the depth first search networks of the infection clusters 
        infection_networks, dfs_networks = infection_odds(merged_aggregated, outbreak_starter, G)
            
        #Making a new graph for the dfs_plots 
        dfs_plot = nx.DiGraph()
        #To join the dfs_trees i have to convert it to a list of tuples..
        for root, network in dfs_networks.items():
            #Iterating through the dfs dict and pulling the nodes and the their edges to build a network plot 
            for node in network.nodes():
                dfs_plot.add_node(node)
            for edge in network.edges():
                infector, infected = edge
                #The edge moves from top down so the first individual is the infector 
                dfs_plot.add_edge(infector, infected)
                #Pulling information for the infector and infected individuals to add attributes for processing in Gephi 
                filtered_df_infector = merged_aggregated[(merged_aggregated['ID'] == infector)]
                filtered_df_infected = merged_aggregated[(merged_aggregated['ID'] == infected)]
                for _, row in filtered_df_infector.iterrows():
                    #Pulls the node associated with the infector from the network plot and adds their infection odds 
                    dfs_plot.nodes[infector]['InfectionOdds'] = row['InfectionOdds']
                    #If this person is a patient then their ward number is added as well 
                    if row['Role'] == 'Patient':
                        dfs_plot.nodes[infector]['WardNumber'] = row['WardNumber']
                    #Adds their date of infection -- converting to a string because pandas timestamps are not supported objects in the NetworkX writer 
                    dfs_plot.nodes[infector]['Date'] = str(row['Date'])
                #Same process as above but for the infected individual 
                for _, row in filtered_df_infected.iterrows():
                    #Infection odds is added
                    dfs_plot.nodes[infected]['InfectionOdds'] = row['InfectionOdds']
                    #Checks if they are listed as a patient
                    if row['Role'] == 'Patient':
                        dfs_plot.nodes[infector]['WardNumber'] = row['WardNumber']
                    #Adds their date of infection 
                    dfs_plot.nodes[infected]['Date'] = str(row['Date'])

        try:
            #Setting the name of the new folder that will store all the simulation data 
            output_folder = 'hospital_infection_sim'
            #The subfolder name will correspond to its simulation number 
            sim_subfolder = f'simulation{simulation_count}'
            #Creates a path to the directory in which the data can be added 
            full_path = os.path.join(output_folder,sim_subfolder)
            os.makedirs(full_path)
            logger.info(f'Folder {sim_subfolder} sucessfully created. Located in: {os.getcwd()}')
        except FileExistsError:
            logger.error('File already exists. Rename folder or continue with existing file.')

        #To be exported
        #Joins the directory path to the name of the file to be exported 
        output_path = os.path.join(full_path, f'master_df_{simulation_count}.csv')  
        #Converting the pandas df to a csv file  
        merged_aggregated.to_csv(output_path, sep='\t')
        #Ward infection cluster data 
        output_path = os.path.join(full_path, f'ward_infection_clusters_{simulation_count}.csv')
        outbreak_wards.to_csv(output_path, sep='\t')
        #Information about UDCA treatment in each infection cluster 
        output_path = os.path.join(full_path, f'infection_network_df_{simulation_count}.csv')   
        infection_networks.to_csv(output_path, sep='\t')
        #The dfs network tree exported as a gephi plot
        #Gephi was used rather than networkx in its native form because there are too many nodes/edges present to have a visually appealing output
        #without any further processing. Gephi allows for further filtering and refinment of these subnetwork plots 
        output_path = os.path.join(full_path, f'gephi_{simulation_count}.gexf')           
        nx.write_gexf(dfs_plot, output_path)

#Calling the above function 
parallel_processing_files(paired_sim_data, pt_folder, hcw_folder)

