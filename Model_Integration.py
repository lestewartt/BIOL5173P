import numpy as np 
from scipy.stats import gamma, poisson 

'''This part of the project deals with integrating the models that were fit in MATLAB
Integrating information from three models:
    1. ACE2 Expression ---> MLE parameters fit to a Gamma distribution (a,b)
        ACE2(t)
        *the distibution of ACE2 expression has been grouped by days so there are 3 sets of parameters: 0, 1-2, 3-5
    2. Virus Entry ---> We have fit both sigmoidal and linear models, and used BIC scores to substantiate the linear model
        Ventry(t) = x(ACE2(t))
    3. Virus Exposure --> Gamma~Poisson mixed model to predict the probability of infection based on various levels of exposure 
        Vexposure - Gamma(alpha,theta)
'''
#MODEL 1 - phatgam parameters 
#These are the Maximum Likelihood Estimation (MLE) parameters for the gamma distribution fit to ACE2 expression data
ACE2_expression_params = {'day0': {'shape': 3.3082, 'scale': 0.0344}, 'day1-2': {'shape': 2.2333, 'scale': 0.0222}, 'day3-5': {'shape': 1.0817, 'scale': 0.0310}}

#MODEL 2 - parameters for the linear model estimating viral entry based on the ACE2 expression levels 
virus_entry_lung_params = {'slope': 30.84, 'intercept': 0.01}
virus_entry_bronch_params = {'slope': 2506.23, 'intercept': -0.12}

#MODEL 3 - parameters for the mixed model fit to the viral exposure data 
virus_exposure_params = {'alpha': 0.1191, 'theta': 3.6439} 

#Calculting the mean reduction in ACE2 expression, by time
'''This was used to validate if the below code logic was sound, as there was some variation in online resources as to the
calculations for finding the mean of a gamma distribution. The results here and below were the same, which suggests that 
multiplying the shape and scale is indeed the correct method. '''
mean_exp_0 = gamma.mean(ACE2_expression_params['day0']['shape'], loc=0, scale=ACE2_expression_params['day0']['scale'])
mean_exp_1_2 = gamma.mean(ACE2_expression_params['day1-2']['shape'], loc=0, scale=ACE2_expression_params['day1-2']['scale'])
mean_exp_3_5 = gamma.mean(ACE2_expression_params['day3-5']['shape'], loc=0, scale=ACE2_expression_params['day3-5']['scale'])

'''These changes in ACE2 pull through directly to effective changes in exposure to the virus, 
as we are using a linear model ( y=kx ) for the relationship between ACE2 expression and virus entry.'''
def viral_entry(slope, ACE2_levels, intercept):
    '''
    Calculates the probability of viral entry based on a linear model relationship with ACE2 expression

    Parameters: 
        slope = slope of the linear model 
        ACE2_levels = ACE2 expression levels 
        intercept = intercept of the linear model 
    '''
    p_viral_entry = slope * ACE2_levels + intercept  
    return p_viral_entry

def exposure_levels(params):
    '''
    Generates exposure levels fit to a gamma distribution 

    MLE parameters obtained from MATLAB code
        alpha: shape parameter 
        theta: scale parameter 

    x: number of exposure levels to be generated 
    '''
    #creating an array from 0.001 to 0.999 with a step of 0.001
    x = np.linspace(0.001,0.999,1000)
    #gamma ppf is the inverse cdf (this is the equivalent to the gaminv function used in Matlab to obtain the gi variable)
    exposures = gamma.ppf(x, a=params['alpha'], scale=params['theta'])
    return exposures 

def mean_reduction_calc(ACE2_exp_params, exposure_params):
    '''
    Calculates the mean reduction in infection probability and viral entry probability for each time period based on 
    ACE2 expression levels for a population treated with UDCA and a population left untreated (i.e. control group). 

    Parameters 
        ACE2_exp_params = dict containing the ACE2 expression parameters for each time period, fit to a gamma distribution 
        exposure_params = dict containing the parameters for the mixed model of virus exposure data 
    '''
    #initializing empty dictionaries to store outputs
    mean_reductions = {}
    probability_untreated = {}
    probability_treated = {}

    #setting to none to be replaced by calculation below
    untreated_exp = None
    #setting to true to iterate through the first time period, day0
    day_0 = True 

    #generating exposure levels for downstream calculations using the previously defined function 
    #compared these exposure values to the ones generated in MATLAB and they are both essentially the same
    exposure = exposure_levels(exposure_params)

    #iterating over each time period and the corresponding parameters 
    for time_period, param in ACE2_exp_params.items():
        #calculating the mean ACE2 expression for each time period (same method as the scipy.gamma function)
        if day_0 is True:
            #day0 is the metric for the untreated population 
            untreated_exp = param['shape'] * param['scale']
            #once this value is extracted for day0 the boolean is changed to false so subsequent iterations are passed through the condition below
            day_0 = False

            #creating an empty list to store the probabilities of infection in an untreated population 
            prob_untreated = []
            #iterating over each exposure level created from the above function, stored in the 'exposure' variable
            for e in exposure: 
                #using the poisson distribution to calculate the probability of infection at each exposure level
                #set to zero to calculate the probability of no infection events  
                prob = 1 - poisson.cdf(0, e)
                #adding the results to the initialized list outside the loop 
                prob_untreated.append(prob)
            #calculating the mean probability of infection at each time point and storing the results in the corresponding dict outside the loop
            probability_untreated[time_period] = np.mean(prob_untreated)
        
        #calculating the mean reduction in ACE2 for each time period relative to the baseline (untreated population)
        #subsequent days will use the untreated value as obtained by the day0 calculations 
        else:
            #calculating the mean ACE2 expression 
            mean_exp = param['shape'] * param['scale']
            #dividing the mean expression for the current time period by the mean expression for day 0 (untreated) to determine the mean reduction in ACE2 expression
            mean_reduction = mean_exp / untreated_exp
            #storing the results in a dict
            mean_reductions[f'mean_reduction_{time_period}'] = mean_reduction
            
            #creating another empty list but for the UDCA treated population 
            prob_treated = []
            #iterating over each exposure level generated 
            for e in exposure: 
                '''using the poisson distribution calculation but adjusting the exposure level, e, to prob
                of viral entry at that time point (adjusted for the mean reduction in ACE2 expression).
                set to zero to calculate the probability of no infection events. 
                models the efficacy of UDCA treatment on reducing the risk of infection '''
                adjusted_exposure = e * mean_reduction
                prob = 1 - poisson.cdf(0, adjusted_exposure)
                #adding the results to the list outside the loop 
                prob_treated.append(prob)
            #calculating the mean probability of infection at each time point and storing the results in dict outside the loop
            probability_treated[time_period] = np.mean(prob_treated)
            
    return mean_reductions, probability_untreated, probability_treated 

#calling the above function and defining variable names to the output 
mean_reductions, probability_untreated, probability_treated = mean_reduction_calc(ACE2_expression_params, virus_exposure_params)

#calculating the infection odds after being treated with UDCA for 1-2 and 3-5 days 
odds_ratio_1_2 = probability_treated.get('day1-2') / probability_untreated.get('day0')
odds_ratio_3_5 = probability_treated.get('day3-5') / probability_untreated.get('day0')

#printing results to terminal 
print(f'Mean Reduction in ACE2 Expression: {mean_reductions}')
print(f'Probability of infection for untreated population: {probability_untreated}')
print(f'Probability of infection for UDCA treated population: {probability_treated}')
print(f'Odds of COVID infection after 1-2 days of UDCA treatment: {odds_ratio_1_2}')
print(f'Odds of COVID infection after 3-5 days of UDCA treatment: {odds_ratio_3_5}')

'''
Results of the above: 
ACE2 expression untreated: 0.11380208
Mean Reduction in ACE2 Expression: {'mean_reduction_day1-2': 0.4356621601292349, 'mean_reduction_day3-5': 0.2946580589739661}
Probability of infection for untreated population: {'day0': 0.16680215971783008}
Probability of infection for UDCA treated population: {'day1-2': 0.1066592432729963, 'day3-5': 0.08278282036496551}
Odds of COVID infection after 1-2 days of UDCA treatment: 0.6394356251347453
Odds of COVID infection after 3-5 days of UDCA treatment: 0.49629345630179245
'''






