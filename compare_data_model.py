import matplotlib.pyplot as plt
import math
import numpy as np
import sys
#import urllib2
import re
import calendar
import datetime
import argparse
from datetime import timedelta
from calendar import monthrange
import os
import matplotlib.mlab as mlab
import matplotlib.cm as cm
import matplotlib.dates as mdates
import seaborn as sns

__version__ = "0.1"
__author__ = "Katie Whitman"
__maintainer__ = "Katie Whitman"
__email__ = "kathryn.whitman@nasa.gov"

#Compare results output by operational_sep_quantities.py
#for the SHINE challenge SEP events.


def check_files(filename):
    """Check if a file exists. Graceful exit if doesn't."""
    exists = os.path.isfile(filename)
    if not exists:
        print("NOTE: File " + filename + " does not exist. Continuing.")
        return False

    return True


def calculate_time_difference(key1, dict1, key2, dict2):
    """Calculate the amount of time between two times in two dictionaries.
       key1 - key to get time 1, e.g. 'proton_start_time'
       dict1 - dictionary with values
       key2 - key to get time 2, e.g. 'proton_end_time'
       dict2 - dictionary with values

       Get the duration of an SEP event by entering the same dictionary for
       dict1 and dict2, but specify different keys.
       Calculate how different data or model answers are for the same value
       by inputting two different dictionaries, but specifying the same value
       for key1 and key2.

       Time difference is given in hours and found by (time2-time1).
    """
    time1 = dict1[key1]
    time2 = dict2[key2]
    #We will save in time units of hours
    time_difference = (time2 - time1).total_seconds()/(60.*60.)

    return time_difference


def read_proton_info(threshold, sep_date, expmt_keys, experiments, flux_types):
    """Read in the proton-related information for all SEP events. Use the
       sep_keys (i.e. dates) to read in the SEP files which were output from
       operational_sep_quantities.py.
       For two threshold crossings (>10 MeV, 10 pfu and >100 MeV, 1 pfu):
            SEP start time, peak time, rise time, end time, duration
            peak flux, >10 MeV fluence, >100 MeV fluence
       SEP filenames like sep_values_GOES-13_differential_2012_3_7.csv
    """
    Nexp = len(expmt_keys)
    proton_dict = {}
    proton_keys = ['proton_start_time','proton_peak_flux','proton_peak_time',
            'proton_rise_time','proton_end_time','proton_duration',
            'fluence_gt10','fluence_gt100']
    date_cols = [2,4,6]
    float_cols = [3,8,9]
    threshold_found = False
    for j in range(Nexp):
        year = sep_date.year
        month = sep_date.month
        day = sep_date.day
        fname = 'shine/sep_values_' + experiments[j] + '_' + flux_types[j] \
                + '_'  + str(year) + '_' + str(month) + '_' + str(day) + '.csv'
        success = check_files(fname)
        #If experiment file doesn't exist
        if not success:
            #empty dictionaries
            proton_dict.update({expmt_keys[j]: {}})
            continue

        #Read in the file
        infile = open(fname,"r")
        #SKIP LINES starting with hash #
        for line in infile:
            line = line.lstrip()
            if line[0] == "#":
                continue
            else:
                #Break up line into different column by splitting at commas
                row = line.split(",")
                #Check if row matches the requested threshold
                if row[0] == threshold[0] and row[1] == threshold[1]:
                    threshold_found = True
                    row_dict = {}
                    for i in range(2,len(row)):
                        val = row[i]  #string
                        for col in date_cols: #if column should be a date
                            if i == col:
                                val = datetime.datetime.strptime(row[i],
                                                       "%Y-%m-%d %H:%M:%S")
                        for col in float_cols: #if column should be a number
                            if i == col:
                                val = float(row[i])

                        if proton_keys[i-2] == 'proton_rise_time': #recalculate
                            val = calculate_time_difference(\
                                    'proton_start_time',row_dict,
                                    'proton_peak_time',row_dict)
                        if proton_keys[i-2] == 'proton_duration': #recalculate
                            val = calculate_time_difference(\
                                    'proton_start_time',row_dict,
                                    'proton_end_time',row_dict)
                        #string all the columns together to make one dict per row
                        row_dict.update({proton_keys[i-2] : val})

                    proton_dict.update({expmt_keys[j]: row_dict})

        #if request threshold was not crossed or calculated for this experiment
        if not threshold_found:
            proton_dict.update({expmt_keys[j]: {}})
            continue

    return proton_keys, proton_dict



def read_fluence_info(threshold, sep_date, expmt_keys, experiments, flux_types):
    """Read in the proton-related information for all SEP events. Use the
       sep_keys (i.e. dates) to read in the SEP files which were output from
       operational_sep_quantities.py.
       For two threshold crossings (>10 MeV, 10 pfu and >100 MeV, 1 pfu):
            SEP start time, peak time, rise time, end time, duration
            peak flux, >10 MeV fluence, >100 MeV fluence
       SEP filenames like sep_values_GOES-13_differential_2012_3_7.csv
    """
    Nexp = len(expmt_keys)
    dict = {}
    fluence_keys = ['Emid','fluence']
    threshold_found = False
    thresh_label = 'gt' + threshold[0][1:]
    for j in range(Nexp):
        year = sep_date.year
        month = sep_date.month
        day = sep_date.day
        fname = 'shine/fluence_' + experiments[j] + '_' + flux_types[j] + '_' \
                + thresh_label + '_'  + str(year) + '_' + str(month) + '_' \
                + str(day) + '.csv'
        success = check_files(fname)
        #If experiment file doesn't exist
        if not success:
            #empty dictionaries
            dict.update({expmt_keys[j]: {}})
            continue

        #Read in the file
        infile = open(fname,"r")
        Emid = []
        fluence = []
        #SKIP LINES starting with hash #
        for line in infile:
            line = line.lstrip()
            if line[0] == "#":
                continue
            else:
                #Break up line into different column by splitting at commas
                row = line.split(",")
                if flux_types[j] == "differential":
                    Emid.append(float(row[1]))
                    fluence.append(float(row[3]))
                if flux_types[j] == "integral":
                    Emid.append(float(row[0]))
                    fluence.append(float(row[1]))

        dict.update({expmt_keys[j]: {'Emid': Emid, 'fluence': fluence}})

    return fluence_keys, dict



def reference_time_difference(keys,ref_key,select_key,dict):
    """Since we are dealing with 2D dictionaries, keys refers to the top
       level set of keys. e.g. list of experiments and flux_types
       ref_key indicates which of the top level set of keys should be used to
       define the reference time. e.g. "GOES-13_integral"
       select_key refers to the 2nd level information that you are interested
       in, e.g. "proton_start_time"
       Returns time difference in hours. Times prior to reference are
       negative. Times later than reference are positive.
    """
    N = len(keys)
    time_difference = []
    for i in range(N):
        if not dict[keys[i]]: #empty dictionary
            time_difference.append([])
            continue
        time_diff = calculate_time_difference(select_key, dict[ref_key],
                    select_key, dict[keys[i]])
        time_difference.append(time_diff)


    return time_difference


def reference_ratio(keys,ref_key,select_key,dict):
    """Since we are dealing with 2D dictionaries, keys refers to the top
       level set of keys. e.g. list of experiments and flux_types
       ref_key indicates which of the top level set of keys should be used to
       define the reference time. e.g. "GOES-13_integral"
       select_key refers to the 2nd level information that you are interested
       in, e.g. "proton_start_time"
       Returns ratio of values w.r.t. reference = val/ref.
       > 1 indicates higher than reference. < 1 indicates lower than reference.
    """
    N = len(keys)
    ratios = []
    for i in range(N):
        if not dict[keys[i]]: #empty dictionary
            ratios.append([])
            continue
        ratio= dict[keys[i]][select_key]/dict[ref_key][select_key]
        ratios.append(ratio)

    return ratios


def make_array_from_dict(keys, select_key, dict):
    """Loop through higher level keys, e.g. experiments or sep events.
       Pull out the selected info, designated in select_key.
                e.g. 'proton_peak_flux'
       Put that info into an array that can be plotted.
    """
    N = len(keys)
    list = []
    for i in range(N):
        if not dict[keys[i]]:
            list.append([])
        else:
            val = dict[keys[i]][select_key]
            list.append(val)

    return list






if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--SEPEventDate", type=str, default='', help=("Start "
                    "date of SEP event in YYYY-MM-DD or a sequence of "
                    "SEP dates separated by commas with no spaces, "
                    "e.g. 2012-03-07,2012-05-17"))
    parser.add_argument("--ModelName", type=str, default='', help=("Name of "
                "your model. Please use operational_sep_quantities.py "
                "to create results for your model and then rename them with "
                "names e.g."
                "sep_values_ModelName_differential_2012_5_17.csv,  "
                "fluence_ModelName_integral_gt10_2012_5_17.csv, and"
                "fluence_ModelName_integral_gt100_2012_5_17.csv"))
    parser.add_argument("--FluxType", type=str, choices=['integral',
            'differential'], default='differential',help=("Are your model "
                    "results derived from differential or "
                    "integral fluxes? (Default is differential)"))
    parser.add_argument("--Threshold", type=str, default="100,1",
            help=("The energy and flux thresholds (written as 100,1 with no "
                    "spaces) which will be used to plot the event. "
                    "e.g. 100,1 indicates >100 MeV fluxes crossing 1 pfu "
                    "(1/[cm^2 s sr])."))
    parser.add_argument("--showplot",help="Set flag to display plots", \
                    action="store_true")


    args = parser.parse_args()
    sep_dates = args.SEPEventDate.strip().split(",")
    model_name = args.ModelName
    model_flux_type = args.FluxType
    threshold = args.Threshold.strip().split(",")
    threshold[0] = '>'+threshold[0]
    showplot = args.showplot

    sep_keys = []
    for i in range(len(sep_dates)):
        date = datetime.datetime.strptime(sep_dates[i], "%Y-%m-%d")
        sep_keys.append(date)
    Nsep = len(sep_keys)
    experiments = ['GOES-13','GOES-13','GOES-15','GOES-15','SEPEM','SEPEM_BG',
                model_name]
    flux_types = ['integral','differential','integral','differential',
            'differential','differential',model_flux_type]
    Nexp = len(experiments) #number of comparisons for plots
    expmt_keys = [] #keys defining the experiment + flux type combos
    for i in range(Nexp):
        expmt_keys.append(experiments[i] + ' ' + flux_types[i])

    #READ IN SEP PROTON FLUX FOR SPECIFIED SEP EVENTS AND ALL EXPERIMENT
    #AND MODEL COMBINATIONS
    proton_dict = {}
    for i in range(Nsep):
        proton_keys, proton_one = read_proton_info(threshold,\
                sep_keys[i], expmt_keys, experiments, flux_types)
        proton_dict.update({sep_keys[i]:proton_one})

    ref_key = 'GOES-13 integral'  #choose reference data set
    #RATIO of peak flux and time difference w.r.t. GOES-13
    ref_time_diff_all = []
    ref_peak_ratio_all = []
    for i in range(Nsep):
        ref_time_diff = reference_time_difference(expmt_keys,ref_key,
                    'proton_peak_time',proton_dict[sep_keys[i]])
        ref_peak_ratio = reference_ratio(expmt_keys,ref_key,'proton_peak_flux',
                    proton_dict[sep_keys[i]])
        ref_time_diff_all.append(ref_time_diff)
        ref_peak_ratio_all.append(ref_peak_ratio)


    #FLUENCE COMPARISON
    fluence_dict = {}
    for i in range(Nsep):
        fluence_keys, fluence_one = read_fluence_info(threshold, sep_keys[i],
                        expmt_keys, experiments, flux_types)
        fluence_dict.update({sep_keys[i]:fluence_one})


    #-----------PLOTTING-----------
    if showplot:
        thresh_label = threshold[0] + ' MeV, ' + threshold[1] + ' pfu'
        sns.set()
        colors = cm.rainbow(np.linspace(0, 1, Nexp))
        #PLOT OF PROTON PEAK TIME DIFFERENCE VS PROTON PEAK FLUX RATIO
        for i in range(Nsep):
            plt.figure(figsize=(8,5))
            ax = plt.subplot(111)
            #plt.grid(which="both", axis="both")
            plt.title(str(sep_keys[i])[0:10] + ' SEP Peak ('+thresh_label+ ')')
            plt.xlabel('Hours from GOES-13 Integral Start Time')
            plt.ylabel('Ratio: Proton ' + threshold[0] + ' MeV\n Peak '
                        'Flux/GOES-13 Integral')
            for j in range(Nexp):
                if not ref_peak_ratio_all[i][j]:
                    continue #skip experiment combo if doesn't exist
                ax.plot(ref_time_diff_all[i][j], ref_peak_ratio_all[i][j],
                        'bo', color=colors[j],
                        label=(experiments[j] + ", " + flux_types[j]))
            chartBox = ax.get_position()
            ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.75,
                             chartBox.height])
            ax.legend(loc='upper center', bbox_to_anchor=(1.25, 0.95))


        #PLOT OF EVENT-INTEGRATED FLUENCE SPECTRA
        for i in range(Nsep):
            plt.figure(figsize=(8,6))
            ax = plt.subplot(111)
            #plt.grid(which="both", axis="both")
            plt.title(str(sep_keys[i])[0:10] + ' Event-Integrated Fluence '
                        'Spectrum ('+thresh_label+ ')')
            plt.xlabel('Energy [MeV]')
            if model_flux_type == "integral":
                plt.ylabel('Integral Flux 1/[cm^2 s sr]')
            if model_flux_type == "differential":
                plt.ylabel('Differential Flux 1/[MeV cm^2 s sr]')
            for j in range(Nexp):
                if flux_types[j] != model_flux_type:
                    continue  #only plot same type of flux
                if not fluence_dict[sep_keys[i]][expmt_keys[j]]:
                    continue #no values for experiment
                energy_bins = fluence_dict[sep_keys[i]][expmt_keys[j]]['Emid']
                spectrum = fluence_dict[sep_keys[i]][expmt_keys[j]]['fluence']
               # if experiments[j] == model_name:
               #     colors = 'red'
               # else:
                #    colors = 'blue'

                ax.plot(energy_bins, spectrum,'bo', color = colors[j],
                        label=(experiments[j] + ", " + flux_types[j]))
            chartBox = ax.get_position()
            ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.75,
                             chartBox.height])
            ax.legend(loc='upper center', bbox_to_anchor=(1.25, 0.95))
            plt.xscale("log")
            plt.yscale("log")


        #BAR CHART OF SEP EVENT PEAK FLUXES - ALL EXPERIMENTS
        for i in range(Nsep):
            peak_fluxes = make_array_from_dict(expmt_keys, 'proton_peak_flux',
                proton_dict[sep_keys[i]])
            plt.figure(figsize=(8,5))
            colors = []
            labels = []
            for j in range(Nexp):
                labels.append(experiments[j] + '\n' + flux_types[j])
                #set colors for bars
                if experiments[j] == model_name:
                    colors.append('red')
                else:
                    colors.append('blue')
            for j in range(Nexp-1,-1,-1):
                if peak_fluxes[j] == []:
                    del peak_fluxes[j]
                    del labels[j]
                    del colors[j]
                #else:
                #    #Convert to percent
                 #   ref_peak_ratio_all[i][j] = ref_peak_ratio_all[i][j]*100.

            y_pos = np.arange(len(labels))
            plt.bar(y_pos, peak_fluxes, color=colors, align='center', alpha=0.5)
            plt.xticks(y_pos, labels)
            plt.ylabel(threshold[0] + ' MeV Peak Flux 1/[cm^2 s sr]')
            plt.title(str(sep_keys[i])[0:10] +' Event ' + threshold[0]
                    +' MeV Peak Flux (' + thresh_label + ')')


        #BAR CHART OF SEP EVENT PEAK FLUXES - ERROR BAR PLOT W/GOES-13 REFERENCE
        for i in range(Nsep):
            peak_fluxes = make_array_from_dict(expmt_keys, 'proton_peak_flux',
                proton_dict[sep_keys[i]])
            plt.figure(figsize=(8,5))
            colors = []
            labels = []
            yerrl = 1e9
            yerrh = -1e9
            selected_peak_flux = []
            for j in range(Nexp):
                #Pull out only the reference value
                if expmt_keys[j] == ref_key:
                    ref_peak_flux = peak_fluxes[j]
                    labels.append(experiments[j] + ' ' + flux_types[j] +
                            '\n with data range')
                    colors.append('blue')
                #set colors for bars
                if experiments[j] == model_name:
                    model_peak_flux = peak_fluxes[j]
                    labels.append(experiments[j] + '\n' + flux_types[j])
                    colors.append('red')
                #Get the range of data values
                if peak_fluxes[j] != [] and experiments[j] != model_name:
                    val = peak_fluxes[j] - ref_peak_flux
                    if val < yerrl:
                        yerrl = val
                    if val > yerrh:
                        yerrh = val

            print('yerrl = ' + str(yerrl) + ' yerrh = ' + str(yerrh))
            y_pos = np.arange(len(labels))
            plt.bar(y_pos, [ref_peak_flux,model_peak_flux],
                    yerr = [[abs(yerrl),0],[yerrh,0]], color=colors, align='center',
                    alpha=0.5,capsize=10)
            plt.xticks(y_pos, labels)
            plt.ylabel(threshold[0] + ' MeV Peak Flux 1/[cm^2 s sr]')
            plt.title(str(sep_keys[i])[0:10] +' Event ' + threshold[0] +' MeV Peak '
                            'Flux (' + thresh_label + ')')



        #BAR CHART OF SEP EVENT RISE TIMES
        for i in range(Nsep):
            rise_time = make_array_from_dict(expmt_keys, 'proton_rise_time',
                proton_dict[sep_keys[i]])
            plt.figure(figsize=(8,5))
            colors = []
            labels = []
            for j in range(Nexp):
                labels.append(experiments[j] + '\n' + flux_types[j])
                #set colors for bars
                if experiments[j] == model_name:
                    colors.append('red')
                else:
                    colors.append('blue')
            for j in range(len(rise_time)-1,-1,-1):
                if rise_time[j] == []:
                    del rise_time[j]
                    del labels[j]
                    del colors[j]

            y_pos = np.arange(len(labels))
            plt.bar(y_pos, rise_time, align='center', color=colors, alpha=0.5)
            plt.xticks(y_pos, labels)
            plt.ylabel('Rise Time [Hours]')
            plt.title(str(sep_keys[i])[0:10] +' Event Rise Time ('
                        + thresh_label + ')')


        #BAR CHART OF SEP EVENT DURATIONS
        for i in range(Nsep):
            durations = make_array_from_dict(expmt_keys, 'proton_duration',
                proton_dict[sep_keys[i]])
            plt.figure(figsize=(8,5))
            colors = []
            labels = []
            for j in range(Nexp):
                labels.append(experiments[j] + '\n' + flux_types[j])
                #set colors for bars
                if experiments[j] == model_name:
                    colors.append('red')
                else:
                    colors.append('blue')
            for j in range(len(durations)-1,-1,-1):
                if durations[j] == []:
                    del durations[j]
                    del labels[j]
                    del colors[j]

            y_pos = np.arange(len(labels))
            plt.bar(y_pos, durations, align='center', color=colors,alpha=0.5)
            plt.xticks(y_pos, labels)
            plt.ylabel('Duration [Hours]')
            plt.title(str(sep_keys[i])[0:10] +' Event Duration ('
                        + thresh_label + ')')


        #BAR CHART OF SEP FLUENCE for >10 MeV
        for i in range(Nsep):
            fluence = make_array_from_dict(expmt_keys,
                    'fluence_gt10',proton_dict[sep_keys[i]])
            plt.figure(figsize=(8,5))
            colors = []
            labels = []
            for j in range(Nexp):
                labels.append(experiments[j] + '\n' + flux_types[j])
                #set colors for bars
                if experiments[j] == model_name:
                    colors.append('red')
                else:
                    colors.append('blue')
            for j in range(len(fluence)-1,-1,-1):
                if fluence[j] == []:
                    del fluence[j]
                    del labels[j]
                    del colors[j]

            y_pos = np.arange(len(labels))
            plt.bar(y_pos, fluence, align='center', color=colors,alpha=0.5)
            plt.xticks(y_pos, labels)
            plt.ylabel('>10 MeV Fluence [cm^-2]')
            plt.title(str(sep_keys[i])[0:10] +' Event-Integrated >10 MeV '
                        'Fluence ('+ thresh_label + ')')


        #BAR CHART OF SEP FLUENCE for >100 MeV
        for i in range(Nsep):
            fluence = make_array_from_dict(expmt_keys,
                    'fluence_gt100',proton_dict[sep_keys[i]])
            plt.figure(figsize=(8,5))
            colors = []
            labels = []
            for j in range(Nexp):
                labels.append(experiments[j] + '\n' + flux_types[j])
                #set colors for bars
                if experiments[j] == model_name:
                    colors.append('red')
                else:
                    colors.append('blue')
            for j in range(len(fluence)-1,-1,-1):
                if fluence[j] == []:
                    del fluence[j]
                    del labels[j]
                    del colors[j]

            y_pos = np.arange(len(labels))
            plt.bar(y_pos, fluence, align='center', color=colors,alpha=0.5)
            plt.xticks(y_pos, labels)
            plt.ylabel('>100 MeV Fluence [cm^-2]')
            plt.title(str(sep_keys[i])[0:10] +' Event-Integrated >100 MeV '
                        'Fluence (' + thresh_label + ')')


        #BAR CHART OF START TIMES WITH RESPECT TO GOES-13 INTEGRAL
        for i in range(Nsep):
            ref_time_diff = reference_time_difference(expmt_keys,ref_key,
                        'proton_start_time',proton_dict[sep_keys[i]])
            plt.figure(figsize=(8,5))
            colors = []
            labels = []
            for j in range(Nexp):
                labels.append(experiments[j] + '\n' + flux_types[j])
                #set colors for bars
                if experiments[j] == model_name:
                    colors.append('red')
                else:
                    colors.append('blue')
            for j in range(len(ref_time_diff)-1,-1,-1):
                if ref_time_diff[j] == []:
                    del ref_time_diff[j]
                    del labels[j]
                    del colors[j]

            y_pos = np.arange(len(labels))
            plt.bar(y_pos, ref_time_diff, align='center', color=colors,alpha=0.5)
            plt.xticks(y_pos, labels)
            plt.ylabel('Time from GOES-13 Integral Start Time [Hours]')
            plt.title(str(sep_keys[i])[0:10] +' Event Start Times' +
                        ' (' + thresh_label + ')')

        #BAR CHART OF PEAK TIMES WITH RESPECT TO GOES-13 INTEGRAL
        for i in range(Nsep):
            ref_time_diff = reference_time_difference(expmt_keys,ref_key,
                        'proton_peak_time',proton_dict[sep_keys[i]])
            plt.figure(figsize=(8,5))
            colors = []
            labels = []
            for j in range(Nexp):
                labels.append(experiments[j] + '\n' + flux_types[j])
                #set colors for bars
                if experiments[j] == model_name:
                    colors.append('red')
                else:
                    colors.append('blue')
            for j in range(len(ref_time_diff)-1,-1,-1):
                if ref_time_diff[j] == []:
                    del ref_time_diff[j]
                    del labels[j]
                    del colors[j]

            y_pos = np.arange(len(labels))
            plt.bar(y_pos, ref_time_diff, align='center', color=colors,alpha=0.5)
            plt.xticks(y_pos, labels)
            plt.ylabel('Time from GOES-13 Integral Peak Time [Hours]')
            plt.title(str(sep_keys[i])[0:10] +' Event Peak Times' +
                        ' (' + thresh_label + ')')




        #TIMELINE PLOT OF SEP START TIMES
        levels = np.tile([-5, 5, -3, 3, -1, 1],int(np.ceil(Nexp/6)))[:Nexp]
        for i in range(Nsep):
            start_times = make_array_from_dict(expmt_keys, 'proton_start_time',
               proton_dict[sep_keys[i]])
            fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
            ax.set(title=str(sep_keys[i])[0:10] + ' Event Start Times ('
                        + thresh_label + ')')
            labels = []
            for j in range(Nexp):
                labels.append(experiments[j] + '\n' + flux_types[j])
            for j in range(len(start_times)-1,-1,-1):
                if start_times[j] == []: #missing experiment
                    del start_times[j]
                    del labels[j]
            if start_times == []:
                continue
            markerline, stemline, baseline = ax.stem(start_times,
                levels[:len(start_times)], linefmt='C3-', basefmt='k-')
            plt.setp(markerline, mec="k", mfc="w", zorder=3)
            markerline.set_ydata(np.zeros(len(start_times)))
            vert = np.array(['top', 'bottom'])[(levels > 0).astype(int)]
            for d, l, r, va in zip(start_times, levels, labels, vert):
                ax.annotate(r, xy=(d, l), xytext=(-3, np.sign(l)*3),
                    textcoords="offset points", va=va, ha="right")

            #SET TIME LABEL POSITIONS to 2 hour divisions
            ax.get_xaxis().set_major_locator(mdates.HourLocator(interval=1))
            ax.get_xaxis().set_major_formatter(mdates.DateFormatter("%b-%d\n "
                            "%H:%M"))
            plt.setp(ax.get_xticklabels(), rotation=30, ha="right")

            ax.get_yaxis().set_visible(False)
            for spine in ["left", "top", "right"]:
                ax.spines[spine].set_visible(False)
            ax.margins(y=0.1)



        #TIMELINE PLOT OF SEP PEAK TIMES
        levels = np.tile([-5, 5, -3, 3, -1, 1],int(np.ceil(Nexp/6)))[:Nexp]
        for i in range(Nsep):
            peak_times = make_array_from_dict(expmt_keys, 'proton_peak_time',
               proton_dict[sep_keys[i]])
            fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
            ax.set(title=str(sep_keys[i])[0:10] + ' Event Peak Times ('
                        + thresh_label + ')')
            labels = []
            for j in range(Nexp):
                labels.append(experiments[j] + '\n' + flux_types[j])
            for j in range(len(peak_times)-1,-1,-1):
                if peak_times[j] == []: #missing experiment
                    del peak_times[j]
                    del labels[j]
            if peak_times == []:
                continue
            markerline, stemline, baseline = ax.stem(peak_times,
                levels[:len(peak_times)], linefmt='C3-', basefmt='k-')
            plt.setp(markerline, mec="k", mfc="w", zorder=3)
            markerline.set_ydata(np.zeros(len(peak_times)))
            vert = np.array(['top', 'bottom'])[(levels > 0).astype(int)]
            for d, l, r, va in zip(peak_times, levels, labels, vert):
                ax.annotate(r, xy=(d, l), xytext=(-3, np.sign(l)*3),
                    textcoords="offset points", va=va, ha="right")

            #SET TIME LABEL POSITIONS to 2 hour divisions
            ax.get_xaxis().set_major_locator(mdates.HourLocator(interval=1))
            ax.get_xaxis().set_major_formatter(mdates.DateFormatter("%b-%d\n "
                        "%H:%M"))
            plt.setp(ax.get_xticklabels(), rotation=30, ha="right")

            ax.get_yaxis().set_visible(False)
            for spine in ["left", "top", "right"]:
                ax.spines[spine].set_visible(False)
            ax.margins(y=0.1)


        # #TIMELINE PLOT OF SEP END TIMES
        # levels = np.tile([-5, 5, -3, 3, -1, 1],int(np.ceil(Nexp/6)))[:Nexp]
        # for i in range(Nsep):
        #     end_times = make_array_from_dict(expmt_keys, 'proton_end_time',
        #             proton_gt10_dict[sep_keys[i]])
        #     fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
        #     ax.set(title=str(sep_keys[i])[0:10] + " SEP End Times")
        #     labels = []
        #     for k in range(Nexp):
        #         labels.append(experiments[k] + '\n' + flux_types[k])
        #     for j in range(len(end_times)-1,-1,-1):
        #         if end_times[j] == []:
        #             del end_times[j]
        #             del labels[j]
        #     markerline, stemline, baseline = ax.stem(end_times,
        #         levels[:len(end_times)], linefmt='C3-', basefmt='k-')
        #     plt.setp(markerline, mec="k", mfc="w", zorder=3)
        #     markerline.set_ydata(np.zeros(len(end_times)))
        #     vert = np.array(['top', 'bottom'])[(levels > 0).astype(int)]
        #     for d, l, r, va in zip(end_times, levels, labels, vert):
        #         ax.annotate(r, xy=(d, l), xytext=(-3, np.sign(l)*3),
        #             textcoords="offset points", va=va, ha="right")
        #
        #     #SET TIME LABEL POSITIONS to 2 hour divisions
        #     ax.get_xaxis().set_major_locator(mdates.HourLocator(interval=2))
        #     ax.get_xaxis().set_major_formatter(mdates.DateFormatter("%b-%d\n "
        #                     "%H:%M"))
        #     plt.setp(ax.get_xticklabels(), rotation=30, ha="right")
        #
        #     ax.get_yaxis().set_visible(False)
        #     for spine in ["left", "top", "right"]:
        #         ax.spines[spine].set_visible(False)
        #     ax.margins(y=0.1)

        print('If plots are empty, then no data was available for requested '
            'thresholds, or no data files were present.')
        plt.show()
