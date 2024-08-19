
%%%%%%%%%%%%%%%%%%% PART 1: MODELLING %%%%%%%%%%%%%%%%%% 
%GUID: 2915700s
%Supervised by: Dr. Chris Illingworth 

% This portion of the project derives the parameters that will be used 
% to build a model of how UDCA affects viral transmission. 
% This model has three components: expression model, viral entry model, and
% the viral exposure model 

%To run this script, please assign data files to the variables 
%corresponding to their model use:

%Expression model 
%Data describing ACE2 expression levels on days 0, 1, 2, 3, 4, and 5 of UDCA treatment.  
%From Brevini et al., 2022
responses = readmatrix('Responses.txt');

% Virus entry model 
%Data describing the levels of SARS-CoV-2 infection in treated and untreated 
%lung and bronchial tissue
%From Brevini et al., 2022
file1 = 'lung_master.dat';
lung_master = readtable(file1, "Delimiter", "\t");
file2 = 'bronch_master.dat';
bronch_master = readtable(file2, "Delimiter", "\t");

% Virus exposure model
%Data describing transmission bottlenecks for SARS-CoV-2 
%From Lythgoe et al., 2021
file1 = 'Lythgoe_bottlenecks.dat';
bottleneck_data = readtable(file1, "Delimiter", "\t");
file2 = 'augmented_data.dat';
augmented_data = readtable(file2, "Delimiter", "\t");

% Comment on University of Glasgow's AI Guidance and Policy:
% 
% Artificial intelligence should be used appropriately and responsibly, 
% and it should not be used as a replacement for independent thinking.
% All use of ChatGPT was as an information tool and not for directly generating code.
% ChatGPT was used to assist in debugging, highlighted by **
% 
% Statement of Academic Integrity: This submission is a product of this student's work and 
% approach to solving the designated task. 

%%%%%%%%%%%%%%%%%%% ACE2 EXPRESSION MODEL %%%%%%%%%%%%%%%%%% 
%File of temporal data describing changes in ACE2 expression levels during
%a course of treatment with UDCA 

%Isolating the column that records the number of days the data was recorded
time = responses(:,1);
%Deleting the first column so the data alone can be processed downstream 
responses(:,1) = [];
%Naming the columns to correspond to the sample individual -- necessary for
%ensuring the legend in future plots is formatted well
columns = {'Sample1', 'Sample2', 'Sample3', 'Sample4', 'Sample5', 'Sample6'}; 
%making a table so i can extract the variable names -- this is important to
%iterate over each column rather than by row 
response_data = array2table(responses, 'VariableNames',columns);

%%%%%%%%%%%%%%%%%%% DATA VISUALIZATION

%Plotting data -- using a multi-line timeseries plot to visualize changes
%in ACE2 expression over time for each sample 
VarNames=response_data.Properties.VariableNames;

%Assigning a unique colour to each sample to be used for plotting
column_num = length(VarNames);
colormap = turbo(column_num); 

figure;
%The below command allows for more than one plot to be added without
%deleting the existing plot, thus creating a multi-line plot. 
hold on;
%Establishing that i will be iterating over each data entry by column 
for i=1:length(VarNames)
    %Pulling the data corresponding to the current iterated column 
    col_data = response_data{:,i};
    %Pulling a unique color to be used for plotting this data from the
    %colormap 
    current_color = colormap(i,:);
    %Plotting the data accordingly
    plot(time, col_data, '-pentagram', 'Color', current_color, 'LineWidth', 1.5);
end 
%Once all the data has been added to the plot, further details can be
%implemented -- cant do this in the loop has to be added at the end 
title('TimeSeries ACE2 Expression');
xlabel('Time (days)');
ylabel('ACE2 Expression');
%Ensuring the x-axis does not exceed the range of the temporal data
xlim([0, 5]);
%Populating the legend with the sample names previously defined 
legend(VarNames); 
hold off; 
%Plots will be stored as an svg file to the current working directory 
%gcf stands for 'get current figure'
saveas(gcf,'TimeSeriesPlot.svg');

%%%%%%%%%%%%%%%%%%% FIT A STATISTICAL MODEL

%The data must now be iterated over by day rather than by sample to fit to
%a gamma distribution 
%gamfit requires the days to be in columns not rows, accomplished with an
%apostrophe 
transposed_responses = responses';
%Initalizing empty cells for estimated params for each day alone to be stored 
initial_params = cell(1,6); 
for i = 1:6
    %Selecting one column (same system as above but slightly different for
    %an array than a table --> i.e. () instead of {}) 
    day_data = transposed_responses(:, i);
    %Saving the estimated parameters for each day to the respective table 
    [initial_params{i}] = gamfit(day_data);
end

%This next section groups the days into different time periods and
%produces the gamma distribution parameters for each to be used to calculate the 
%negative log likelihood as a means of comparison 

%Since day 0 is the measure of ACE2 levels before UDCA treatment, it is
%being isolated as a control measure 
initial_params_day0 = initial_params{:,1};
day0_data = transposed_responses(:,1);

%Pulling the corresponding data for the days to be grouped and fitting them
%to a gamma distribution 
day0to1_data = [transposed_responses(:,1);transposed_responses(:,2)];
[phat0_1, pci0_1] = gamfit(day0to1_data);

day0to2_data = [transposed_responses(:,1);transposed_responses(:,2);transposed_responses(:,3)];
[phat0_2, pci0_2] = gamfit(day0to2_data);

day1to2_data = [transposed_responses(:,2);transposed_responses(:,3)];
[phat1_2, pci1_2] = gamfit(day1to2_data);

day1to3_data = [transposed_responses(:,2);transposed_responses(:,3);transposed_responses(:,4)];
[phat1_3, pci1_3] = gamfit(day1to3_data);

day3to5_data = [transposed_responses(:,4);transposed_responses(:,5);transposed_responses(:,6)];
[phat3_5, pci3_5] = gamfit(day3to5_data);

day3to4_data = [transposed_responses(:,4);transposed_responses(:,5)];
[phat3_4, pci3_4] = gamfit(day3to4_data);

day4to5_data = [transposed_responses(:,5);transposed_responses(:,6)];
[phat4_5, pci4_5] = gamfit(day4to5_data);

initial_params_day5 = initial_params{:,6};
day5_data = transposed_responses(:,6);

day1to5_data = [transposed_responses(:,2);transposed_responses(:,3);transposed_responses(:,4);transposed_responses(:,5);transposed_responses(:,6)];
[phat1_5, pci1_5] = gamfit(day1to5_data);

%Using the function gamlike to compute negative log likelihood using the
%parameters for the corresponding time periods 
%function for the pdf of a gamma distribution ==> y = gampdf(x,a,b)
%gamlike is a utility function for the MLE of a gamma distribution 
    %it calculates the negative gamma log-likelihood, so minimizing gamlike
    %using fminsearch is the same as maximizing the likelihood 
nlogL0 = gamlike(initial_params_day0, day0_data);
logL0 = -nlogL0;

nlogL0_2 = gamlike(phat0_2, day0to2_data);
logL_0_2 = -nlogL0_2;

nlogL0_1 = gamlike(phat0_1, day0to1_data);
logL_0_1 = -nlogL0_1;

nlogL1_2 = gamlike(phat1_2, day1to2_data);
logL_1_2 = -nlogL1_2;

nlogL1_3 = gamlike(phat1_3, day1to3_data);
logL_1_3 = -nlogL1_3;

nlog3_5 = gamlike(phat3_5, day3to5_data);
logL_3_5 = -nlog3_5;

nlog3_4 = gamlike(phat3_4, day3to4_data);
logL_3_4 = -nlog3_4;

nlogL4_5 = gamlike(phat4_5, day4to5_data);
logL4_5 = -nlogL4_5;

nlogL5 = gamlike(initial_params_day5, day5_data);
logL5 = -nlogL5;

nlogL1_5 = gamlike(phat1_5, day1to5_data);
logL1_5 = -nlogL1_5;

%Another metric for comparison is bayesion information criterion (BIC),
%There are variations of this equation but the version here takes the log
%likelihood calculated above 
function bic = bayesioninfo(params)
    %L is the log likelihood 
    L = params(1);
    %n is the sample size
    n = params(2);
    %k is the number of parameters 
    k = params(3);
    bic = -2*L + k*log(n);
end

%Setting variables for the parameters according to each time period: the
%log likelihood, sample size, and number of parameters 
params0 = [logL0, 6, 2];
bic0 = bayesioninfo(params0);

params1_3 = [logL_1_3, 18, 2];
bic1_3 = bayesioninfo(params1_3);

params1_2 = [logL_1_2, 12, 2];
bic1_2 = bayesioninfo(params1_2);

params3_5 = [logL_3_5, 18, 2];
bic3_5 = bayesioninfo(params3_5);

params3_4 = [logL_3_4, 12, 2];
bic3_4 = bayesioninfo(params3_4);

params4_5 = [logL4_5, 12, 2];
bic4_5 = bayesioninfo(params4_5);

params5 = [logL5, 6, 2];
bic5 = bayesioninfo(params5);

params1_5 = [logL1_5, 30, 2];
bic1_5 = bayesioninfo(params1_5);

%the optimal groupings are selected according to the returned log
%likelihood and BIC values 

%%%%%%%%%%%%%%%%%%% PLOTTING PDFs OF THE GAMMA DISTRIBUTION 

%Pulling the shape and scale parameters for each chosen time period 
shape0 = initial_params_day0(1);
scale0 = initial_params_day0(2);

shape3_5 = phat3_5(1);
scale3_5 = phat3_5(2);

shape1_2 = phat1_2(1);
scale1_2 = phat1_2(2);

% Define the range of x-values for plotting the PDF
x_range3_5 = linspace(min(day3to5_data), max(day3to5_data), 100); 

% Calculate the PDF values for the entire range
pdf_day3_5 = gampdf(x_range3_5, shape3_5, scale3_5);

% % % Plot the PDF
% figure;
% hold on;
% %plotting the PDF (model prediction) 
% plot(x_range3_5, pdf_day3_5, 'LineWidth', 2, 'Color', '#A2142F');
% %plotting the empirical data using kernel denstiy estimate 
% kde_estimate = ksdensity(day3to5_data, x_range3_5);
% plot(x_range3_5, kde_estimate, 'LineWidth', 2, 'Color', '#000000');
% %another visualization of the empirical data is via histogram 
% h = histogram(day3to5_data);
% %modifies the display of the histogram to a PDF 
% h.Normalization = 'pdf';
% %changes the colour of the histogram 
% h.FaceColor = '#0072BD';
% %formatting the title and axis names of the plot 
% title('Gamma Distribution Fit for Days 3 to 5');
% xlabel('ACE2 Expression Level');
% ylabel('Density');
% %applying a legend 
% legend('Gamma PDF' ,'KDE', 'Location','best')
% hold off;
% %saving the plot as an SVG file 
% %gcf stands for 'get current figure'
% saveas(gcf,'pdf_day3_5.svg');


% Define the range of x-values for plotting the PDF
x_range1_2 = linspace(min(day1to2_data), max(day1to2_data),100); 

% Calculate the PDF values for the entire range
pdf_day1_2 = gampdf(x_range1_2, shape1_2, scale1_2);

% %Plot the PDF
% figure;
% hold on;
% %plotting the PDF (model prediction) 
% plot(x_range1_2, pdf_day1_2, 'LineWidth', 2, 'Color', '#A2142F');
% %plotting the empirical data using kernel denstiy estimate 
% kde_estimate = ksdensity(day1to2_data, x_range1_2);
% plot(x_range1_2, kde_estimate, 'LineWidth', 2, 'Color', '#000000');
% %another visualization of the empirical data is via histogram 
% h = histogram(day1to2_data);
% %modifies the display of the histogram to a PDF 
% h.Normalization = 'pdf';
% %changes the colour of the histogram 
% h.FaceColor = '#0072BD';
% %formatting the title and axis names of the plot
% title('Gamma Distribution Fit for Days 1 to 2');
% xlabel('ACE2 Expression Level');
% ylabel('Density');
% %applying a legend 
% legend('Gamma PDF' ,'KDE', 'Location','best')
% hold off;
% %saving as an SVG plot 
% %gcf stands for 'get current figure'
% saveas(gcf,'pdf_day1_2.svg');



% Define the range of x-values for plotting the PDF
x_range0 = linspace(min(day0_data), max(day0_data), 100); 

% Calculate the PDF values for the entire range
pdf_day0 = gampdf(x_range0, shape0, scale0);


% % Plot the PDF
% figure;
% hold on;
% %plotting the PDF (model prediction) 
% plot(x_range0, pdf_day0, 'LineWidth', 2, 'Color', '#A2142F');
% %plotting the empirical data using kernel denstiy estimate 
% kde_estimate = ksdensity(day0_data, x_range0);
% plot(x_range0, kde_estimate, 'LineWidth', 2, 'Color','#000000');
% %another visualization of the empirical data is via histogram 
% h = histogram(day0_data);
% %modifies the display of the histogram to a PDF 
% h.Normalization = 'pdf';
% %changes the colour of the histogram 
% h.FaceColor = '#0072BD';
% %formatting the title and axis names of the plot
% title('Gamma Distribution Fit for Day 0');
% xlabel('ACE2 Expression Level');
% ylabel('Density');
% %applying a legend 
% legend('Gamma PDF' ,'KDE', 'Location','best')
% hold off;
% %saving as an SVG file 
% %gcf stands for 'get current figure'
% saveas(gcf,'pdf_day0.svg');


%%plotting them all together%% 
%storing together to be accessed in the loop below 
all_data = {day0_data, day1to2_data, day3to5_data}; 
x_range_all = {x_range0, x_range1_2, x_range3_5};
all_pdfs = {pdf_day0, pdf_day1_2, pdf_day3_5};

time_periods = {'Day_0', 'Day1-2', 'Day3-5'};

pdf_colour = {'#DD7208', '#0A8728', '#0303C0'};
empirical_data_colour = {'#DD2f08', '#14C003', '#1B9CCD'};

%initializing a figure 
figure;
hold  on;
%storing the key for each plot to be input to the legend 
legend_entries = {};
%iterating over the time period data 
for i=1:length(time_periods)
    %pulling the current info for the time period being iterated over 
    current_x_range = x_range_all{i};
    current_data = all_data{i};
    current_pdf = all_pdfs{i}; 
    %pulling associated colors for this time period (defined above) 
    current_pdf_colour = pdf_colour{i};
    current_empirical_data_colour = empirical_data_colour{i};
    %plotting the pdf for the current time period
    plot(current_x_range, current_pdf, 'LineWidth', 2, 'Color', current_pdf_colour);
    %as each time period is iterated over, the legend is updated
    %accordingly 
    legend_entryPDF = sprintf('%s-PDF', time_periods{i});
    legend_entries = [legend_entries, legend_entryPDF];
    %plotting the empirical data 
    kde_estimate = ksdensity(current_data, current_x_range);
    plot(current_x_range, kde_estimate, 'LineWidth', 2, 'Color', current_empirical_data_colour);
    %the legend will have the corresponding time period and associated
    %color to the predicted and empirical data 
    legend_entryKDE = sprintf('%s-KDE', time_periods{i});
    legend_entries = [legend_entries, legend_entryKDE]; 

end
title('Gamma Distribution Fit Across All Time Periods');
xlabel('ACE2 Expression Level');
ylabel('Probability Density'); 
legend(legend_entries, 'Location', 'best');
hold off;
%gcf stands for 'get current figure'
saveas(gcf,'PDF_KDE_all_time_periods.svg');


%%%%%%%%%%%%%%%%%%% VIRUS ENTRY MODEL %%%%%%%%%%%%%%%%%%%%%% 
%This portion of the script takes the lung and bronchial data and fits them
%to linear and logistic models to identify which best characterizes the
%relationship between changes in ACE2 expression and changes in the level
%of SARS-CoV-2 infection

%Converting to numerical format with just numerical data (i.e. no sample
%type info) 
lung_data = [double(lung_master.ACE2), double(lung_master.COVID)];
bronch_data = [double(bronch_master.Var2), double(bronch_master.Var3)];

%Isolating the columns by ACE2 and SARS-CoV-2 data 
lung_ACE2 = lung_master.ACE2;
lung_COVID = lung_master.COVID;

bronch_ACE2 = bronch_master.Var2;
bronch_COVID = bronch_master.Var3;

%Defining the data as numerical input 
x_lung = double(lung_ACE2);
y_lung = double(lung_COVID);

x_bronch = double(bronch_ACE2);
y_bronch = double(bronch_COVID);

%Using the curve fitter app shows that the lung data best suits logarithmic and
%linear models, whereas bronch data also fits sigmoidal 

%%%%%%%%%%%%%%%%%%% LINEAR MODEL 

%Fitting the data to a linear model using fitlm 
lung_linear = fitlm(x_lung, y_lung);
bronch_linear = fitlm(x_bronch, y_bronch); 

%plotting the results of the linear model for lung and bronchial data 

%disp(lung_linear);
figure;
plot(lung_linear);
xlabel('ACE2 Levels')
ylabel('COVID Levels')
title('Lung Data Linear Fit')
%gcf stands for 'get current figure'
saveas(gcf,'lung_linear_plot.svg');


%disp(bronch_linear);
figure;
plot(bronch_linear);
xlabel('ACE2 Levels')
ylabel('COVID Levels')
title('Bronchi Data Linear Fit')
%gcf stands for 'get current figure'
saveas(gcf,'bronch_linear_plot.svg');

%%%%%%%%%%%%%%%%%%% LOGISTIC MODEL 

%sigmoid function 
function yi = sigmoid(params, x)
    %upper limit 
    L = params(1);
    %growth rate 
    k = params(2);
    %inflection point (i.e. midpoint) 
    ip = params(3); 
    %sigmoidal equation
    yi = L./(1 + exp(-k.*(x-ip)));
end

%fitting to a logistic model using nlinfit
%L is defined as the maximum of the y axis, the growth rate is commonly
%denoted as 1**, and the inflection point is the mean of the x axis
l_params = [max(y_lung), 1, mean(x_lung)];
[l_beta, l_R, l_J, l_CovB, l_MSE] = nlinfit(x_lung, y_lung, @sigmoid, l_params);

b_params = [max(y_bronch), 1, mean(x_bronch)];
[b_beta, b_R, b_J, b_CovB, b_MSE] = nlinfit(x_bronch, y_bronch, @sigmoid, b_params);


%%%%%%%%%%%%%%%%%%% CALCULATING BIC 

%While the linear model function output does provide the loglikelihood,
%nlinfit does not, so we are calculating the BIC differently from the previous model, 
% using the residual sum of squares 
RSS_lung = sum(l_R.^2);
RSS_bronch = sum(b_R.^2);

%finding BIC with RSS 
function bic = bicrss(params)
    %sample size
    n = params(1);
    %residual sum of squares
    rss = params(2);
    %number of parameters 
    k = params(3); 
    bic = n.*log(rss./n) + k.*log(n);
end

%Calculating the BIC for the logistic models 
lung_sig_bic_rss_params = [8,RSS_lung,2];
lung_sig_bic_rss = bicrss(lung_sig_bic_rss_params);

bronch_sig_bic_rss_params = [8,RSS_bronch,2];
bronch_sig_bic_rss = bicrss(bronch_sig_bic_rss_params); 

%Calculating the BIC for the linear models  
linear_residuals_lung = lung_linear.Residuals.Raw;
RSS_lung_linear = sum(linear_residuals_lung.^2);

linear_res_bronch = bronch_linear.Residuals.Raw;
RSS_bronch_linear = sum(linear_res_bronch.^2);

lung_lin_bic_rss_params = [8,RSS_lung_linear,2];
lung_lin_bic_rss = bicrss(lung_lin_bic_rss_params);

bronch_lin_bic_rss_params = [8,RSS_bronch_linear, 2];
bronch_lin_bic_rss = bicrss(bronch_lin_bic_rss_params); 

%%%%%%%%%%%%%%%%%%% VALIDATING MODEL ACCURACY    
%Created a test to validate the accuracy of the linear and sigmoidal models
%at estimating parameters. Defined known parameters and 

%initializing a table to store the below results 
linsim_results = table;
for i = 1:50
    %true parameters
    n = 8;
    m = 0.3;
    b = 1;
    st_dev = 0.2; 
    %creating 8 random x data points from 0 to 10
    x = linspace(0,10,8); 
    %creating noise to simulate real data (st. dev is kept small)**
    noise = randn(size(x)) * st_dev; 
    %creating y data points to fit a linear model 
    y = m.*x + b + noise; 
    %fitting the generated data to a linear model using fitlm 
    fitsim_lin = fitlm(x,y);
    %isolating the resulting intercept and slope for evaluation 
    intercept = fitsim_lin.Coefficients.Estimate(1); 
    slope = fitsim_lin.Coefficients.Estimate(2);  
    %mse will help identify the results' closeness to the true parameters 
    mse = fitsim_lin.MSE;
    %adding the results of each iteration to the table established outside
    %the loop 
    addrow = table(i, intercept, slope, mse, 'VariableNames', {'#', 'b', 'm', 'mse'});
    linsim_results = [linsim_results; addrow];
end


%establishing a table for results to be stored 
sigsim_results = table;
for i = 1:50 
    %true parameters 
    L = 2;
    k = 5;
    ip = 1;
    %generating 8 random points from 0 to 1.8 (generated a sigmoidal plot
    %to ensure these metrics created a sensical sigmoidal curve) 
    x_sig = linspace(0,1.8,8); 
    st_dev_sig = 0.01;
    %creating values for noise based around the values of x **
    noise_sig = randn(size(x_sig)) * st_dev_sig;
    %generating y values with the added noise to simulate environmental data 
    sig = L./(1 + exp(-k.*(x_sig-ip)));
    y_noise_sig = sig + noise_sig; 
    %mimicking my model's parameter estimation model to test accuracy 
    sim_params = [max(y_noise_sig), 1, mean(x_sig)];
    [sim_beta, sim_R, sim_J, sim_CovB, sim_MSE] = nlinfit(x_sig, y_noise_sig, @sigmoid, sim_params);
    %generating the predicted y values based on the estimated parameters 
    y_sim_pred = sigmoid(sim_beta, x_sig);
    %calculating mse for accuracy report 
    RSS = sum(sim_R.^2);
    mse = RSS/8;

    %adding results to the table 
    addrow = table(sim_beta(1), sim_beta(2), sim_beta(3), sim_MSE, mse, 'VariableNames', {'L', 'k', 'ip','sim_mse', 'calc_mse'});
    sigsim_results = [sigsim_results; addrow];
end


%%%%%%%%%%%%%%%%%%% VIRUS EXPOSURE MODEL %%%%%%%%%%%%%%%%%%%%%% 
%The final model takes the bottleneck data to fit a model describing levels of exposure 
% to the SARS-CoV-2 virus in household situations

%A review of household transmission studies has found a secondary attack rate of 18.9% 
% i.e. for every 189 transmission events there are 811 cases of non-transmission.  
%We will need to include these cases of non-transmission in our model.
%Getting the count of transmission events 
transmission_events = height(bottleneck_data);
%Using the above information to determine the degree to which the data
%needs to be augmented 
nontransmission_events_needed = round((transmission_events * 811) / 189);
%converting the augmented data table to an array 
augmented_data = table2array(augmented_data);

%testing for overdispersion -- relevant for using a negative binomial
%distribution
mean_count = mean(augmented_data);
var_count = var(augmented_data);
dispersion = var_count / mean_count;
%if dispersion ratio is > 1 then the data is overdispersed --> very
%dispersion: 3.64

%%%%%%%%%%%%%%%%%%% NEGATIVE BINOMIAL MODEL 

%nbinfit returns the MLEs and confidence intervals for the parameters of
%the negative binomial distribution 
%phat contains the MLE for the parameters 
[phat,pci] = nbinfit(augmented_data);

%visualising the model's goodness of fit by plotting the CDF against the
%empirical data 
figure;
hold on;
%CDF of a negative binomial distribution 
plot(0:10,nbincdf(0:10,phat(1),phat(2)),".-", 'LineWidth',2, 'Color','#A2142F');
%plotting the actual data as a histogram 
h = histogram(augmented_data);
%modifying the display to fit a cdf 
h.Normalization = 'cdf';
%changing the colour 
h.FaceColor = '#0072BD';
%formatting the name of the title and x and y axis 
title('Transmission Bottleneck Negative Binomial Fit');
xlabel("Transmission Trials");
ylabel("Cumulative Probability");
hold off; 
%saving as an SVG file 
%gcf stands for 'get current figure'
saveas(gcf,'nbincdf_plot_transmissionbottleneck.svg');

%The negative binomial model describes the transmission bottlenecks we would expect in a household 
%situation.  We want to pull out the details of this, learning the distribution of 
%exposure levels.  For this we need an explicit Gamma-Poisson model.

%gamma dist model: phatgam = MLE params 
[phatgam, pcigam] = gamfit(augmented_data);
 
%creating a table of probability values to be input into the inverse CDF
%function below 
p = table; 
count = 0; 
for i = 1:999
    %starting from 0.001 to 999
    count = count + 0.001;
    addrow = table(count, 'VariableNames', {'probabilities'});
    p = [p; addrow];
end 
p = table2array(p);

%Using the poisson density function (PDF) for each gi to
%calculate the probabilities at each exposure level 
%Find the mean probability of observation for any x calculated across all of the 999 exposure levels gi.

%calculates the log likelihood of observing the data under a mixed
%gamma-poisson model to evaluate the goodness of fit. 
function loglike_mixedmodel = gamma_poisson_loglike(data, phatgam, p)
    %shape parameter
    a = phatgam(1); 
    %scale parameter
    b = phatgam(2);

    % x = gaminv(p,a,b) returns the icdf of the gamma distribution with shape
    % parameter a and the scale parameter b, evaluated at the values in p.
    % Using the parameters found in gamfit from the augmented data 
    %The output, contained in gi, are the levels of exposure 
    gi = gaminv(p, a, b);

    % For each value in the inverse CDF (gi), find the
    % probability of observing a bottleneck of any given size x, based on a
    % Poisson model with rate gi.  

    %making a table to match the dimensions of the exposure and bottleneck data
    probabilities = zeros(length(data), length(gi));
    %iterating over the exposure levels 
    for i = 1:length(gi)
        %this is the current exposure level 
        exposure = gi(i); 
        %iterating over the observed bottlenecks from the augmented data  
        for x = 1:length(data)
            %this is the current bottleneck size 
            bottleneck = data(x); 
            %pdf and adding results to the probabilities table 
            %must specify the row and column in which to store the results 
            probabilities(x, i) = poisspdf(bottleneck, exposure);
        end
    end

    %Find the mean probability of observation for any x calculated across all of the 999 
    %exposure levels gi.
    %Have to set it to 2 so the dimensions sum across the columns 
    sum_prob = sum(probabilities, 2);
    mean_prob = sum_prob / length(gi);

    %log likelihood 
    loglike_mixedmodel = sum(log(mean_prob)); 
end
%calling the above function to produce the loglikelihood 
mixed_loglike = gamma_poisson_loglike(augmented_data, phatgam, p);
