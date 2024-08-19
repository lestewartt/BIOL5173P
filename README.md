# BIOL5173P
Master's Thesis Project; GUID: 2915700s

##### Modelling the Potential Effect of Drug Treatment on SARS-CoV-2 Transmission #####

##### GUID: 2915700s 
##### Supervisor: Dr. Chris Illingworth 
##### A project submitted in partial fulfilment of the requirements for the MSc Bioinformatics Degree at The University of Glasgow

#### Project Aim: Each of the scripts in this repository serve to quantify the effect of UDCA on SARS-CoV-2 infection odds and simulate these effects in a hospital environment based on a targeted treatment regime. The results of which are then processed and analysed to assess epidemiological conditions that may dictate the success rate of UDCA. 


## Installation 

First, clone this repository:
<!-- start:code block -->
git clone https://github.com/lestewartt/BIOL5173P.git
<!-- end:code block -->
Or download as a zip file

## Once installed, the script should be run in the following order: 

# File #1 -- MATLAB SCRIPT 'Fitting_Models.m' 

This portion of the project derives the parameters that will be used to build a model of how UDCA affects viral transmission. 
This model has three components: expression model, viral entry model, and the viral exposure model 

Files can be accessed from the 'Test_Data/MATLAB_Data' folder

To run this script, please assign data files to the variables corresponding to their model use:

Expression model 
Data describing ACE2 expression levels on days 0, 1, 2, 3, 4, and 5 of UDCA treatment.  
From Brevini et al., 2022
<!-- start:code block -->
responses = readmatrix('Responses.txt');
<!-- end:code block -->

Virus entry model 
Data describing the levels of SARS-CoV-2 infection in treated and untreated 
lung and bronchial tissue
From Brevini et al., 2022
<!-- start:code block -->
file1 = 'lung_master.dat';
lung_master = readtable(file1, "Delimiter", "\t");
file2 = 'bronch_master.dat';
bronch_master = readtable(file2, "Delimiter", "\t");
<!-- end:code block -->

Virus exposure model
Data describing transmission bottlenecks for SARS-CoV-2 
From Lythgoe et al., 2021
<!-- start:code block -->
file1 = 'Lythgoe_bottlenecks.dat';
bottleneck_data = readtable(file1, "Delimiter", "\t");
file2 = 'augmented_data.dat';
augmented_data = readtable(file2, "Delimiter", "\t");
<!-- end:code block -->


# File #2 -- PYTHON SCRIPT 'Model_Integration.py' 

This part of the project deals with integrating the models that were fit in MATLAB
    1. ACE2 Expression
    2. Virus Entry
    3. Virus Exposure

This script can be run without direct user input 
Exceptions include cases in which ACE2_expression_params and virus_exposure_params are to be updated with new parameter values 


<img width="580" alt="Screenshot 2024-08-18 at 7 38 46 PM" src="https://github.com/user-attachments/assets/35e00730-e4ac-4e3e-a1ed-2100be19b8c7">

Flowchart created using Whimsical: https://whimsical.com/flowcharts


# File #3 -- PYTHON SCRIPT 'Hospital_UDCA_Simulation.py' 

This script accepts the 180 simulations from the provided hospital simulation data (Evans et al., 2024)
To run this script, please update the 'pt_folder' and 'hcw_folder' file paths (lines 41 & 43) to the appropriate directory where this data is stored.
This script can be run successfully as long as there is one pair of pt and hcw data input

Files can be accessed from the 'Simulation_Test_Data' folder

File path to patient simulation data folder
<!-- start:code block -->
pt_folder = 'C:/Users/2915700s/OneDrive/UoG Masters/Project_Start/patient_sim_data'
<!-- end:code block -->
File path to healthcare worker simulation data folder
<!-- start:code block -->
hcw_folder = 'C:/Users/2915700s/OneDrive/UoG Masters/Project_Start/hcw_sim_data'
<!-- end:code block -->

<img width="414" alt="image" src="https://github.com/user-attachments/assets/12424463-30a3-4858-a01e-50c7268c1334">


Flowchart created using Whimsical: https://whimsical.com/flowcharts


# File #4 -- PYTHON SCRIPT 'Processing_Sim_Results.py' 

This script generates summary information for the data produced for each simulation in the hospital infection simulation script 
To run this script, please update the 'main_directory' variable with the path to the 'hospital_infection_sim' folder
Note, this script was designed to process all 180 simulation results. If running using only a single subset of hospital data, code from line 278 onwards should be commented out.
Run using all 3 simulations provided to simulate the effect with all 180 simulations 

Files can be accesed from the 'Processing_UDCA_Model' folder

<!-- start:code block -->
main_directory = 'C:/Users/2915700s/OneDrive/UoG Masters/Project_Start/hospital_infection_sim'
<!-- end:code block -->

# Program Executables Manual 
Describes the output of each script and how to further process results for analysis 
