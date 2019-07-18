import operational_sep_quantities as sep
import compare_data_model as compare

__version__ = "0.4"
__author__ = "Katie Whitman"
__maintainer__ = "Katie Whitman"
__email__ = "kathryn.whitman@nasa.gov"

"""This and supporting codes are found in the respository:
        https://github.com/ktindiana/operational-sep

   EDIT THE INFORMATION STARTING BELOW LINE 23 IN THIS CODE TO MATCH YOUR MODEL
   AND THE AVAILABLE EVENT DATES.
   Indicated by !!!! EDIT HERE !!!! and !!!! END EDITS !!!!

   IN operational_sep_quantities.py:
   You also need to edit variables around lines 38 - 52 to be specific for your
   model. EDIT:
        user_col
        user_delim
        user_energy_bins
"""

#!!!!!!!!!!!!!!!!!!! EDIT HERE !!!!!!!!!!!!!!!!!!!!!!
model_start_dates = ['2012-03-07','2012-05-17','2017-09-10']
model_end_dates = ['2012-03-14','2012-05-20','2017-09-16']
#Specify model file names associated with the SEP events as a list
#---expect files to be in data folder---
model_file_names = ['MODEL/sample_model_mar2012.csv',
                    'MODEL/sample_model_may2012.csv',
                    'MODEL/sample_model_sep2017.csv']

your_model_name = 'MODEL'
model_flux_type = 'integral'

#List all SEP dates available for your model as a single string (comma
#separated, no spaces), e.g.:   '2012-03-07,2012-05-17'
#May or may not be the same as model_start_dates. Should be date on the day
#the SEP flux crossed threshold.
#Plots will be created for these dates
sep_dates = '2012-03-07,2012-05-17,2017-09-10'

#Choose whether to run operational_sep_quantities for model or data.
run_model = True
run_data = False
#!!!!!!!!!!!!!!!!!!! END EDITS !!!!!!!!!!!!!!!!!!!!!!


#---RUN YOUR MODEL FOR SELECTED SHINE SEP EVENTS---
if run_model:
    Nsep = len(model_start_dates)
    for i in range(Nsep):
        start_date = model_start_dates[i]
        end_date = model_end_dates[i]
        experiment = 'user'
        model_name = your_model_name
        user_file = model_file_names[i]
        showplot = False #set False if don't want to see plots
        detect_prev_event = True #probably doesn't hurt, turn off if need to
        threshold = '100,1' #default; modify to add threshold to 10,10 and 100,1
        #CALCULATE SEP INFO AND OUTPUT RESULTS
        try:
            sep.run_all(start_date, end_date, experiment, model_flux_type,
                model_name, user_file, showplot, detect_prev_event, threshold)
        except SystemExit:
            continue



#---RUN ALL GOES AND SEPEM EXPERIMENTS FOR ALL SHINE SEP EVENTS---
#--------------------DO NOT MODIFY--------------------------------
#(although you can if you want to try different events)
sep_start_dates = ['2012-03-07','2012-05-17','2017-09-10']
sep_end_dates = ['2012-03-14','2012-05-20','2017-09-16']
experiments = ['GOES-13','GOES-13','GOES-15','GOES-15','SEPEM']
flux_types = ['integral','differential','integral','differential',
                'differential']

if run_data:
    Nsep = len(sep_start_dates)
    Nexp = len(experiments)
    for i in range(Nsep):
        for j in range(Nexp):
            start_date = sep_start_dates[i]
            end_date = sep_end_dates[i]
            experiment = experiments[j]
            flux_type = flux_types[j]
            model_name = '' #leave as default
            user_file = '' #leave as default
            showplot = False #set True if want to see plots
            detect_prev_event = True #doesn't hurt for these events
            threshold = '100,1' #default; modify to add threshold to 10,10 and 100,1
            #CALCULATE SEP INFO AND OUTPUT RESULTS
            try:
                sep.run_all(start_date, end_date, experiment, flux_type, model_name,
                    user_file, showplot, detect_prev_event, threshold)
            except SystemExit:
                continue



#MAKE PLOTS COMPARING MEASURMENTS TO MODEL
thresholds = ['10,10','100,1'] #<--- Modify if model only predicts one threshold
showplot = False
for str_threshold in thresholds:
    compare.run_all(sep_dates, your_model_name, model_flux_type, str_threshold,
                    showplot)
