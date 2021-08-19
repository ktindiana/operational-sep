from library import read_datasets as datasets
from library import global_vars as vars
from library import ccmc_json_handler as ccmc_json
from library import derive_background as bgsub
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
import os
import wget
from calendar import monthrange
import urllib.request
import csv
from dateutil.parser import parse
import scipy.integrate
from numpy import exp
import array as arr
import pandas as pd
import scipy
from scipy import signal
from statistics import mode

__version__ = "3.0"
__author__ = "Katie Whitman"
__maintainer__ = "Katie Whitman"
__email__ = "kathryn.whitman@nasa.gov"

################ CHANGE LOG ####################
#CHANGES in 0.4: allows users to specify user model or experiment name for
#output files
#CHANGES in 0.4: user passes filename through an argument if selects "user" for
#experiment rather than setting filename at the top of the code
#Changes in 0.5: Added ability to download and plot SOHO/ERNE data for the
#>10 MeV threshold calculations. -- INCOMPLETE
#Changes in 0.5: Indicate peak flux selection on plots.
#Changes in 0.6: SEP event considered to start only after 3 consecutive points.
#   Start time set to the first of the three points.
#Changes in 0.6: Added the --TwoPeaks flag to allow users to indicate that
#   the event has an initial threshold crossing for just a few points,
#   drops below threshold, then continues the increase to the main part of
#   the event. This is phenomenalogical addition after finding a few events
#   with this problem.
#Changes in 0.7: The main function now returns flags to indicate various
#   things that may help a user identify if an event has been poorly timed.
#   These return flags are generated with the idea that this could would
#   be run in batch mode for multiple events without looking at the output.
#   The flags would help a user decide after the fact if event timing should be
#   checked more closely. The flags include:
#       - Event starts at first time point on start date (i.e. previous
#           event already ongoing)
#       - Event is less than 12 hours long (i.e. start time might indicate an
#           initial rise above threshold, but miss main event. TwoPeaks flag
#           was created to handle events like this.)
#       - >100 MeV, 1 pfu threshold start is more than 12 hours after >10 MeV,
#           10 pfu threshold start (e.g. there might be multiple events in a
#           row and the first event is only >10 MeV and second events has >100)
#       - Event ends on very last time point indicating that it ended with the
#           time range rather than dropping below the threshold
#   A flag was added with the "UMASEP" option. When this flag is selected, the
#   code finds all information for four energy channels and thresholds used by
#   UMASEP: >10 MeV, 10 pfu; >100 MeV, 1 pfu; >30 MeV, 1 pfu; >50 MeV, 1 pfu
#   The proton flux in each of these channels is reported for specific times
#   after threshold crossing (Ts) and added to the sep_values_ output file for
#   times Ts + 3, 4, 5, 6, and 7 hours
#   Returns start date of SEP event found so that know the names of the output
#   files produced by this code.
#   Added a saveplot option to automatically write plots to file with a
#   unique filename.
#Changes in 0.7: Adding SEP event onset peak, i.e. the peak associated with the
#   initial acceleration at the sun. Versions 0.6 and previous find the absolute
#   peak in the entire time period between onset and end. For lower energy
#   channels, that peak often corresponds to the ESP portion of the event.
#   Now including a calculation of onset peak.
#Changes in 1.0: Adding support for SOHO/COSTEP/EPHIN differential proton
#   flux data. The 4 energy channels span from 4.3 - 53 MeV.
#   The user must input a threshold in differential flux to calculate all
#   SEP quantities. The operational thresholds cannot be applied.
#   Add ability for user to specify a differential flux threshold to define SEP
#   event.
#Changes in 1.1: Added ability for users to input multiple thresholds. If
#   user uses differential fluxes, then may request a combination of thresholds
#   based off of differential flux bins and estimated integral fluxes.
#Changes in 1.2: Adding options to:
#       Choose corrected or uncorrected GOES fluxes.
#       Choose to apply Bruno (2017) or Sandberg et al. (2014) effective
#       energies to GOES uncorrected data.
#2020-11-13, Changes in 2.0: Complete restructuring of code. A library directory has
#   been created and some of the subroutines originally in
#   operational_sep_quantities is now in read_datasets.py. Global variables,
#   such as directory names and the information users must input to run their
#   data sets, is not in global_vars.py.
#   A code has been written to perform a background subtraction of the SEP
#   flux using a time period specified by the user - presumably spanning a day
#   or more prior to the SEP event. The separation of SEP flux and background,
#   followed by a background subtraction of the SEP flux is in
#   derive_background.py.
#2020-11-24, Changes in 2.1: Added ability to write out json file in CCMC SEP Scoreboard
#   format. More information about those files can be found at
#   https://ccmc.gsfc.nasa.gov/challenges/sep.php#format
#   If a threshold isn't crossed, the code now calculates the maximum flux
#   during the full time period (start date to end date) and the time and saves
#   it in the onset peak and onset date. This info won't be reported in the
#   csv files, but will be reported in the json.
#2021-01-06, Changes in 2.2: DHConsultancy (Daniel Heynderickx) sent out a beta
#   version of the SEPEM RSDv3 data sets. Added here as a native data set.
#2021-01-12, Changes in 2.3: Added functionality to read in REleASE data.
#2021-01-19, Changes in 2.4: If time resolution of data set is >15 minutes,
#   relax three point requirement to exceed or fall below a threshold.
#   For finding onset time, restricted onset to fall between event start and
#   start + 18 hours or start and end time, which ever is shorter.
# 2021-01-25, Changes in 2.4.1: Changed the logic when converting differential
#   flux to integral flux (from_differential_to_integral_flux). Previously,
#   if one bin (bin[i]) had non-zero flux and the next bin (bin[i+1]) had zero
#   flux, the bin[1+1] was set to a value of 1e-15 and then an integral was
#   calculate from bin[i] to 1e-15 to get the flux contribution. This gave
#   strange results. Now, if bin[i] or bin[i+1] has a flux of zero or None, no
#   flux is added to the integral flux.
#2021-02-25, Changes in 2.5: Changing the all_integral_fluences array to
#   all_threshold_fluences. In previous versions, all_integral_fluences
#   contained fluences for integral energy channels only, i.e. >10, >100 MeV.
#   If a user input a differential channel thresold, e.g. for an energy bin
#   of 5 - 7 MeV, then this all_integral_fluxes would hold a NaN value for
#   the fluence corresponding to the threshold array element in
#   all_integral_fluences.
#   Instead, we will now have all_threshold_fluences. The fluence will
#   will correspond to the fluence in each energy channel, whether it
#   is integral or differential. If the user specifies a differential energy
#   channel (e.g. 5 - 7 MeV) with a threshold, then the corresponding element in
#   all_threshold_fluences will have the fluence in the 5 - 7 MeV bin.
#2021-02-28, Changes in 2.5.1: Adjusting the onset peak estimation to fix bugs.
#   Reorganized some of the error checking into subroutines to clean up run_all.
#2021-03-02, Changes in 2.5.2: Added resorting of bins if they are in "reverse"
#   order. For example, SEPMOD produces flux files for energies in order
#   of 1000, 750, 500, ... MeV bins. This doesn't matter for integral fluxes,
#   but from_differential_to_integral flux expects differential bins
#   in increasing energy order. Added sort_bin_order to ensure that
#   energy bins are always in increasing order of effective energies.
#!***!2021-03-03, Changes in 2.5.3: Changed threshold crossing logic.
#   Previously,required flux to exceed threshold (>) for 3 consecutive points
#   to get event start. Changed to flux >= threshold for 3 consecutive points.
#2021-05-17, changes in 2.5.4: Nothing has changed in the
#   operational_sep_quantities.py code, but the ccmc_json_handler.py
#   code was modified to v0.4 and changed the format of the outpuer
#   json file a little bit. I want to mark this change within this code
#   as well and explicitly state the date when it happened.
#2021-05-27, changes in 2.5.5: Added the NoInterp flag to allow users
#   to specify that negative flux or None values should be set to None.
#   If NoInterp is not set, the default behavior is to fill in bad data
#   points using linear interpolation in time. This may not be desired
#   for model output as it inherently changes the model predictions.
#   Zeroes are always treated as valid values and are not replaced.
#   If there are gaps in the time steps, the code does NOT try to
#   fill in the gaps. It will only perform interpolation for time
#   steps present in the input data set.
#   Additionally, changed the way in which the time resolution was
#   calculated in calculate_fluence. Previously, the time resolution
#   was determined by the time difference between the first and
#   second time points. Now the difference in time is calculated for
#   every consecutive set of time points and the most common value
#   for the difference is used as the time resolution. This is done
#   in case there are time gaps in the data set.
#2021-05-28, changes in 2.6: Discovered an unintended bug in the
#   calculate_onset_peak code. I had intended to normalize the
#   deriv variable by the changing flux value to minimize
#   the derivative for jumps and dips when the flux is already
#   elevated. Instead, I had simply normalized by the very first
#   flux value in the array.
#   deriv = smooth_flux[i][j] - smooth_flux[i][j-nwin]
#   Previously: run_deriv[i].append(deriv/smooth_flux[i][0])
#   Changed to: run_deriv[i].append(deriv/smooth_flux[i][j-nwin])
#   Now using two versions of the derivative to get the onset peak.
#   One is the derivative divided by smooth_flux[i][j-nwin], which
#   really highlights the intial rise. The second is deriv divided
#   by a constant normalization factor (the first non-zero flux in
#   the array). This is used to escape local minima.
#   I intended to always use deriv/smooth_flux[i][j-nwin], but
#   I made a typo in v0.8 that ended up dividing by a constant
#   factor. I feel that combining the two concepts works best.
#!!!!2021-07-19, changes in 3.0: Went up to the next integer version!!!
#   number because this version has reconciled all differences
#   with the CCMC SEP Scoreboard JSON format. i.e. this version
#   produces JSON files and supporting output files in exactly
#   the format required by the SEP Scoreboard. The format
#   is specified at https://ccmc.gsfc.nasa.gov/challenges/sep.php.
#       - JSON in CCMC format
#       - cleaning up code and writing a lot more small subroutines.
#       - changed calculate_fluence to multiply by 4pi to remove
#            units of sr^-1
#       - added global units variables at top of code
#       - added ability for user to set units for a user-read file
#           in global_vars.py
#       - when threshold is not crossed, maximum flux in time period
#           is saved in peak_flux rather than onset_peak, as before
#   2021-08-06: Changing onset peak flux definition to be more
#   independent of applied threshold. For events with maxmimum flux
#   values less than 500/energy_channel, search for the onset peak
#   starting 12 hours before the flux threshold is crossed.
#   Onset peak will still only be derived if thresholds
#   are crossed.
########################################################################

#See full program description in all_program_info() below
datapath = vars.datapath
outpath = vars.outpath
plotpath = vars.plotpath
badval = vars.badval #bad data points will be set to this value; must be negative

#####UNITS#####
#If a user data set is read in, the units in global_vars.py will
#supercede these units.
energy_units = "MeV"
flux_units_integral = "pfu"
fluence_units_integral = "cm^-2"
flux_units_differential = "MeV^-1*cm^-2*s^-1*sr^-1"
fluence_units_differential = "MeV^-1*cm^-2"

######FOR USER DATA SETS######
#(expect the first (0th) column contains date in YYYY-MM-DD HH:MM:SS format)
#Identify columns containing fluxes you want to analyze
user_col = vars.user_col
#DELIMETER between columns; for whitespace separating columns, use " " or ""
user_delim = vars.user_delim
#DEFINE ENERGY BINS associated with user file and columns specified above as:
user_energy_bins = vars.user_energy_bins
############################

#FILENAME(s) containing user input fluxes - WILL BE SET THROUGH ARGUMENT
#Can be list of files that are continuous in time
 #      e.g. user_fname = ["file1.txt","file2.txt"]
user_fname = ['tmp.txt']


def all_program_info(): #only for documentation purposes
    """ Program description for operational_sep_quantities.py v3.0.
    
    Formatting for Sphinx web-based documentation, which can be viewed
    within the docs/index.html directory or online at:
    https://ktindiana.github.io/operational-sep/index.html
    
    This program will calculate various useful pieces of operational
    information about SEP events from GOES-08, -10, -11, -12, -13, -14, -15
    data and the SEPEM (RSDv2 and RSDv3) dataset.

    SEP event values are always calculated for threshold definitions:
        
        * >10 MeV exceeds 10 pfu
        * >100 MeV exceed 1 pfu

    The user may add multiple additional thresholds through the command line.
    This program will check if data is already present in a 'data' directory. If
    not, GOES or EPHIN data will be automatically downloaded from the web. SEPEM
    (RSDv2 and RSDv3) data must be downloaded by the user and unzipped inside
    the 'data' directory. Because the SEPEM data set is so large (every 5
    minutes from 1974 to 2015 for RSDv2 and to 2017 for RSDv3), the program will
    break up the data into yearly files for faster reading.
    
    --NoInterp: Data sets are checked for bad data point (negative or None value
    fluxes) and the default behavior is to fill in those bad data points by
    performing a linear interpolation with time. This choice was made to
    calculate more accurate event-intergrated fluence values from data.
    Interpolation with time is not appropriate for model predictions, as
    it will inherently change the prediction or may not be desired by
    the user for the data set. Turn off interpolation with time by
    setting the --NoInterp flag (or nointerp=True). If the interpolation
    is turned off, negative flux values will be set to None.
    Zeroes are always treated as valid values and are not replaced.
    If there are gaps in the time steps, the code does NOT try to
    fill in the gaps. It will ONLY perform interpolation for time
    steps present in the input data set. i.e. gaps in time are not
    interpolated, only time steps with negative or None flux values.

    The values calculated here are important for space radiation operations:
       
       * Onset time, i.e. time to cross thresholds
       * Onset peak intensity
       * Onset peak time
       * Maximum intensity
       * Time of maximum intensity
       * Rise time (onset to peak)
       * End time, i.e. fall below 0.85*threshold for 3 points (15 mins for GOES)
       * Duration
       * Event-integrated fluences
       * Proton fluxes at various times after threshold crossing (UMASEP option)

    UNITS: User may choose differential proton fluxes (e.g. [MeV s sr cm^2]^-1)
    or integral fluxes (e.g. [s sr cm^2]^-1 or pfu). Default units are:
    MeV, cm, s, sr
    Units are specified in library/global_vars.py.
    The program has no internal checks or requirements on units - EXCEPT FOR
    THE THRESHOLD DEFINITIONS OF >10, 10 and >100, 1.
    If you convert those thresholds in the main program to your units,
    you should be able to generate consistent results.
    Currently there are no features to change units automatically. SEE MORE
    about units in the USER INPUT DATA section below.

    OPTIONS: User may specify various options, that currently only apply to
    GOES data:
        
        * Choose corrected or uncorrected GOES fluxes.
        * Choose to apply Bruno (2017) or Sandberg et al. (2014) effective
          energies to GOES uncorrected data.
    
    .. code-block::
    
        --options uncorrected
        --options uncorrected,S14,Bruno2017 (recommend using background subtraction)

    BACKGROUND SUBTRACTION: Users may choose to perform a background
    subtraction by specifying:
    
    .. code-block::
    
        --SubtractBG --BGStartDate YYYY-MM-DD --BGEndDate YYYY-MM-DD
    
    The user should look at the data and select an appropriate time frame
    prior to the event when the background is calm and well-defined. If
    performing background subtraction, the mean background will be subtracted
    from the fluxes in the SEP time frame (StartDate to EndDate). Plots
    showing the mean background level, the background flux only, and the
    background-subtracted SEP fluxes will be created by derive_background to
    verify the quality of the background estimation and subtraction.

    If a previous event is ongoing and the specified time period starts with a
    threshold already crossed, you may try to set the --DetectPreviousEvent
    flag. If the flux drops below threshold before the next event starts, the
    program will identify the second event. This will only work if the
    threshold is already crossed for the very first time in your specified
    time period, and if the flux drops below threshold before the next event
    starts.

    If the event has an initial increase above threshold for a few points, falls
    below threshold, then continues to increase above threshold again, you
    may try to use the --TwoPeak feature to capture the full duration of the
    event. The initial increase above threshold must be less than a day. An
    example of this scenario can be seen in >100 MeV for 2011-08-04.

    A flag was added with the "UMASEP" option. When this flag is used
    (--UMASEP), the code finds all information for four energy channels and
    thresholds used by UMASEP: >10 MeV, 10 pfu; >100 MeV, 1 pfu; >30 MeV, 1 pfu;
    >50 MeV, 1 pfu. The proton flux in each of these channels is reported for
    multiple times after threshold crossing (Ts). The applied time delays are
    as follows:
        
        * >10 MeV - Ts + 3, 4, 5, 6, 7 hours
        * >30 MeV - Ts + 3, 4, 5, 6, 7 hours
        * >50 MeV - Ts + 3, 4, 5, 6, 7 hours
        * >100 MeV - Ts + 3, 4, 5, 6, 7 hours
    
    --spase_id: If you know the appropriate spase_id for the your model or
    experiment, you may specify it here to be filled in to the json file
    for CCMC's SEP Scoreboard.
    
    Note about EVENT END TIME: Currently, the code calculates the end of the
    event as the first time that the flux drops below threshold*endfac for
    three consecutive points. The endfac variable is specified in
    library/global_vars.py.
    Mimicking a SRAG code, endfac is set to 0.85 (85% of threshold). This
    may be changed in the future to update this to more NOAA SWPC-like logic
    to define the end of an SEP event.
    
        Note about ONSET PEAK: In operational_sep_quantities.py version 3.0,
    the onset peak may be found up to 12 hours prior to the threshold crossing
    that determines the SEP event start time. In previous versions, the onset
    peak could only be found at or after the threshold crossing time. For events
    that crossed threshold by a small amount, this often meant that the
    actual onset peak as identified by eye occurred prior to threshold crossing
    at a lower flux value. In an effort to derive the onset peak from the shape
    of the flux time profile independent of the applied threshold value, the code
    will search for the onset peak 12 hours earlier for events with lower
    flux levels. The documentation for
    operational_sep_quantities.calculate_onset_peak() has further details.
    
    RUN CODE FROM COMMAND LINE (put on one line), e.g.:
    
    .. code-block::
    
        python3 operational_sep_quantities.py --StartDate 2012-05-17
        --EndDate "2012-05-19 12:00:00" --Experiment GOES-13
        --FluxType integral --showplot --saveplot

    RUN CODE FROM COMMAND FOR USER DATA SET (put on one line), e.g.:
    
    .. code-block::
    
        python3 operational_sep_quantities.py --StartDate 2012-05-17
        --EndDate '2012-05-19 12:00:00' --Experiment user --ModelName MyModel
        --UserFile MyFluxes.txt --FluxType integral --showplot

    RUN CODE FROM COMMAND LINE AND PERFORM BACKGROUND SUBTRACTION AND APPLY
    Sandberg et al. (2014) and Bruno (2017) effective energies to the GOES bins.
    (note: cannot bg-subtract GOES integral fluxes), e.g.:
    
    .. code-block::
        
        python3 operational_sep_quantities.py --StartDate 2012-05-17
        --EndDate '2012-05-19 12:00:00' --Experiment GOES-13
        --FluxType differential  --showplot --options uncorrected,S14,Bruno2017
        --SubtractBG --BGStartDate 2012-05-10 --BGEndDate --2012-05-17

    RUN CODE IMPORTED INTO ANOTHER PYTHON PROGRAM, e.g.:
    
    .. code-block::
    
        import operational_sep_quantities as sep
        start_date = '2012-05-17'
        end_date = '2012-05-19 12:00:00'
        experiment = 'GOES-13'
        flux_type = 'integral'
        spase_id = ''
        model_name = '' #if experiment is user, set model_name to describe data set
        user_file = '' #if experiment is user, specify filename containing fluxes
        showplot = True  #Turn to False if don't want to see plots
        saveplot = False #turn to true if you want to save plots to file
        options = '' #various options: S14, Bruno2017, uncorrected
        doBGSub = False #Set true if want to perform background subtraction
        bgstart_date = "2012-05-10" #Dates used to estimate mean background if
        bgend_date = "2012-05-17"   #doBGSub is set to True
        detect_prev_event = True  #Helps if previous event causes high intensities
        two_peaks = False  #Helps if two increases above threshold in one event
        umasep = False #Set to true if want UMASEP values (see explanation above)
        threshold = '' #Add a threshold to 10,10 and 100,1: '30,1' or '4.9-7.3,0.01'
        nointerp = False #Default False; set to True to stop linear interpolatin in time

        sep_year, sep_month,sep_day, jsonfname = sep.run_all(start_date, \
            end_date, experiment, flux_type, model_name, user_file,\
            spase_id, showplot, saveplot, detect_prev_event,  \
            two_peaks, umasep, threshold, options, doBGSub, bgstart_date, \
            bgend_date,nointerp)

    Set the desired directory locations for the data and output at the beginning
    of the program in datapath and outpath. Defaults are 'data' and 'output'.

    In order to calculate the fluence, the program determines time_resolution
    (seconds) by finding the difference between every consecutive set of
    time points in the data set. The most common difference is identified as
    the time resolution. This method should find an accurate time resolution
    even if there are gaps in the time steps.
    If the time steps in the data set are truly irregular, the user will
    have to manually set the time resolution inside the subroutine
    calculate_fluence.

    OUTPUT: This program outputs 3 to 4 files, 1 per defined threshold plus
    a summary file containing all of the values calculated for each threshold.
    A file named as e.g. fluence_GOES-13_differential_gt10_2012_3_7.csv contains
    the event-integrated fluence for each energy channel using the specified
    threshold (gt10) to determine start and stop times.
    A file named as e.g. sep_values_GOES-13_differential_2012_3_7.csv contains
    start time, peak flux, etc, for each of the defined thresholds.

    The program write to file the >10 MeV and >100 MeV time series for the
    date range input by the user. If the original data were integral fluxes,
    then the output files simply contain the >10 and >100 MeV time series from
    the input files. If the original data were differential fluxes, then the
    estimated >10 and >100 MeV fluxes are output as time series.
    
    The program also writes any time profile to file if a threshold was applied
    to it. Each energy channel is written to an independent file, named according
    to CCMC's SEP Scoreboard naming conventions. These accompany a json file
    that contains all quantities calculated for the SEP event. There is also a csv
    file that contains most of the same values recorded in the json file.
    
    Example output for values derived from SEPMOD forecasts for the 2021-05-29
    SEP event with additional threshold >10, 0.001 pfu, >100, 0.0001 pfu,
    >30 MeV, 1 pfu, >50 MeV, 1 pfu, and >60 MeV, 0.079 pfu. Some of the files
    below are only created if a threshold was crossed. A default run would produce
    only 10 and 100 MeV files, sep_values_*.csv and json file:
    
    * fluence_SEPMOD_RT_60min_integral_gt10_2021_5_29.csv
    * fluence_SEPMOD_RT_60min_integral_gt10.0_2021_5_29.csv
    * fluence_SEPMOD_RT_60min_integral_gt30.0_2021_5_29.csv
    * fluence_SEPMOD_RT_60min_integral_gt50.0_2021_5_29.csv
    * fluence_SEPMOD_RT_60min_integral_gt60.0_2021_5_29.csv
    * fluence_SEPMOD_RT_60min_integral_gt100.0_2021_5_29.csv
    * integral_fluxes_SEPMOD_RT_60min_integral_2021_5_29.csv
    * sep_values_SEPMOD_RT_60min_integral_2021_5_29.csv
    * SEPMOD_RT_60min_integral.2021-05-29T000000Z.2021-08-06T172140Z.10.0MeV.txt
    * SEPMOD_RT_60min_integral.2021-05-29T000000Z.2021-08-06T172140Z.10MeV.txt
    * SEPMOD_RT_60min_integral.2021-05-29T000000Z.2021-08-06T172140Z.30.0MeV.txt
    * SEPMOD_RT_60min_integral.2021-05-29T000000Z.2021-08-06T172140Z.50.0MeV.txt
    * SEPMOD_RT_60min_integral.2021-05-29T000000Z.2021-08-06T172140Z.60.0MeV.txt
    * SEPMOD_RT_60min_integral.2021-05-29T000000Z.2021-08-06T172140Z.100.0MeV.txt
    * SEPMOD_RT_60min_integral.2021-05-29T000000Z.2021-08-06T172140Z.100MeV.txt
    * SEPMOD_RT_60min_integral.2021-05-29T000000Z.2021-08-06T172140Z.json
    
    The json and txt files listed above would be the ones that would be appropriate
    to read into the CCMC SEP Scoreboard or to pass to the SEP validation code
    being developed in conjunction with the SEP Scoreboard. The csv files are legacy
    files, but may also be easier for some users to read.
    

    USER INPUT DATA SETS: Users may input their own data set. For example, if an
    SEP modeler would like to feed their own intensity time series into this
    code and calculate all values in exactly the same way they were calculated
    for data, it is possible to do that. Default flux units are
    1/[MeV cm^2 s sr] or 1/[cm^2 s sr] and energy channels in MeV for the default
    thresholds to be correct. You can specify your units in the
    library/global_vars.py file.
    
    You can use any units, as long as you are consistent with energy units in
    energy channel/bin definition and in fluxes and you MODIFY THE THRESHOLD
    VALUES TO REFLECT YOUR UNITS. If you want to use different units, but
    still have the correct operational definitions, you need to modify these
    lines in define_thresholds() below:
    
        * energy_thresholds = [10,100] #MeV; flux for particles of > this MeV
        * flux_thresholds = [10,1] #pfu; exceed this level of intensity
    
    NOTE: The first column in your flux file is assumed to be time in format
    YYYY-MM-DD HH:MM:SS. IMPORTANT FORMATTING!!
    NOTE: The flux file may contain header lines that start with a hash #,
    including blank lines.
    NOTE: Any bad or missing fluxes must be indicated by a negative value.
    NOTE: Put your flux file into the "datapath" directory. Filenames will be
    relative to this path.
    NOTE: Please use only differential or integral channels. Please do not mix
    them. You may have one integral channel in the last bin, as this is the way
    HEPAD works and the code has been written to include that HEPAD >700 MeV
    bin along with lower differential channels.

    USER VARIABLES: The user must modify the following variables in
    library/global_vars.py:
    
        :user_col: identify columns in your file containing fluxes to analyze;
                even if your delimeter is white space, consider the date-time
                column as one single column. SET IN library/global_vars.py.
        :user_delim: delimeter between columns, e.g. " " or ","   Use " " for
                any amount of whitespace. SET IN library/global_vars.py.
        :user_energy_bins: define your energy bins at the top of the code in the
                variable user_energy_bins. Follow the format in the subroutine
                define_energy_bins. SET IN library/global_vars.py.
        :user_fname: specify the name of the file containing the fluxes
                through an argument in the command line. --UserFile  The
                user_fname variable will be updated with that filename. ARGUMENT
        :time_resolution: the program determines time_resolution
                (seconds) by finding the difference between every consecutive
                set of time points in the data set. The most common difference
                is identified as the time resolution. This method should find
                an accurate time resolution even if there are gaps in the
                time steps.. AUTOMATICALLY DETERMINED.
                
    Running the code for a user-input file may look like the example below.
    Note that the --UserFile location is with respect to the "data" directory
    inside the operational-sep directory:
    
    .. code-block::
    
        python3 operational_sep_quantities.py --StartDate 2021-05-29 --EndDate 2021-06-05 --Experiment user --UserFile SEPMOD/Scoreboard/SEPMOD.20210529_000000.20210529_165133.20210529_133005_geo_integral_tseries_timestamped_60min.txt --ModelName SEPMOD_RT_60min --showplot --Threshold "10,0.001;100,0.0001;30,1;50,1;60,0.079" --spase_id "spase://CCMC/SimulationModel/SEPMOD" --FluxType integral
    
    VALUES SPECIFIED IN library/global_vars.py:
    
        :datapath: directory containing data, 'data'
        :outpath: directory for program output, 'output'
        :plotpath: directory for saving plots, 'plots'
        :listpath: directory for lists (for run_multi_sep.py)
        :badval: will set any bad data points to this value
        :endfac: multiplicative factor to define threshold for
                end of event; threshold*endfac (default 0.85)
        :nsigma: number of sigma to define SEP versus background
                flux in background subtraction routine
        :version: if you are running a model or data set, allows you
                to enter a version number
        :user_col: array defining flux columns (0 is always datetime)
        :user_delim: delimeter used to separate the columns in the time
                profile file that you will read in
        :user_energy_bins: energy bins associated with the columns
                specified in user_col
        :energy_units: e.g. "MeV'
        :flux_units_integral: e.g. "pfu"
        :fluence_units_integral: e.g. "cm^-2"
        :flux_units_differential: e.g. "MeV^-1*cm^-2*s^-1*sr^-1" (CCMC format)
        :fluence_units_differential: e.g. "MeV^-1*cm^-2" (CCMC format)
        
        (setting the units here will make correct units on plots and in json
        file, but doesn't change operational threshold values; must be done
        accordingly by hand)
    """


def from_differential_to_integral_flux(experiment, min_energy, energy_bins,
                fluxes, options):
    """ If user selected differential fluxes, convert to integral fluxes to
        caluculate operational threshold crossings (>10 MeV protons exceed 10
        pfu, >100 MeV protons exceed 1 pfu).
        Assume that the measured fluxes correspond to the center of the energy
        bin and use power law interpolation to extrapolate integral fluxes
        above user input min_energy.
        The intent is to calculate >10 MeV and >100 MeV fluxes, but leaving
        flexibility for user to define the minimum energy for the integral flux.
        An integral flux will be provided for each timestamp (e.g. every 5 mins).
       
        INPUTS:
        
        :experiment: (string)
        :min_energy: (float) - bottom energy for integral flux calculation
        :energy_bins: (float 1xn array) - bins for each energy channel
        :fluxes: (float nxm array) - fluxes with time for each energy channel
        :options: (string array) - effective energies for GOES channels
        
        OUTPUTS:
        
        :integral_fluxes: (float 1xm array) - estimate integral flux for >min_energy
            (Returns all zero values if no energy bins above min_energy)
            
    """
    print('Converting differential flux to integral flux for >'
            + str(min_energy) + 'MeV.')
    nbins = len(energy_bins)
    nflux = len(fluxes[0])
    #Check requested min_energy inside data energy range
    if min_energy < energy_bins[0][0] or min_energy >= energy_bins[nbins-1][0]:
        print('The selected minimum energy ' + str(min_energy) + ' to create'
                ' integral fluxes is outside of the range of the data: '
                + str(energy_bins[0][0]) + ' - ' + str(energy_bins[nbins-1][0]))
        print('Setting all >'+ str(min_energy) + ' fluxes to zero.')
        integral_fluxes = [0]*nflux
        return integral_fluxes

    #Calculate bin center in log space for each energy bin
    bin_center = []
    for i in range(nbins):
        if energy_bins[i][1] != -1:
            centerE = math.sqrt(energy_bins[i][0]*energy_bins[i][1])
        else:
            centerE = -1
        bin_center.append(centerE)

    #The highest energy EPEAD bin overlaps with all of the HEPAD bins
    #For this reason, excluding the highest energy EPEAD bin in
    #integral flux estimation, 110 - 900 MeV
    if (experiment == "GOES-13" or experiment == "GOES-14" or
        experiment == "GOES-15") and "Bruno2017" not in options:
        remove_bin = -1
        for i in range(nbins):
            if energy_bins[i][0] == 110 and energy_bins[i][1] == 900:
                remove_bin = i
        if remove_bin == -1:
            sys.exit("Attempting to remove 110 - 900 MeV bin for "
                    + experiment + '. Cannot locate bin. Please check '
                    'define_energy_bins to see if GOES-13 to GOES-15 bins '
                    'include 110 - 900 MeV. if not please comment out this '
                    'section in from_differential_to_integral_flux.')
        fluxes = np.delete(fluxes,remove_bin,0)
        energy_bins = np.delete(energy_bins,remove_bin,0)
        bin_center = np.delete(bin_center,remove_bin,0)
        nbins = nbins-1

    #integrate by power law interpolation; assume spectrum the shape of a
    #power law from one bin center to the next. This accounts for the case
    #where the minimum energy falls inside of a bin or there are overlapping
    #energy bins or gaps between bins.
    #An energy bin value of -1 (e.g. [700,-1]) indicates infinity - already an
    #integral channel. This happens for HEPAD. If the integral channel is the
    #last channel, then the flux will be added. If it is a middle bin, it will
    #be skipped.
    integral_fluxes = []
    integral_fluxes_check = []
    for j in range(nflux):  #flux at each time
        sum_flux = 0
        ninc = 0 #number of energy bins included in integral flux estimate
        for i in range(nbins-1):
            if bin_center[i+1] < min_energy:
                continue
            else:
                if energy_bins[i][1] == -1 or energy_bins[i+1][1] == -1:
                    #bin is already an integral flux, e.g. last HEPAD bin
                    continue
                if fluxes[i,j] < 0 or fluxes[i+1,j] < 0: #bad data
                    sys.exit('Bad flux data value of ' + str(fluxes[i,j])
                            + ' and ' + str(fluxes[i+1,j]) +
                            ' found for bin ' + str(i) + ','
                            + str(j) + '. This should not happen. '
                            + 'Did you call check_for_bad_data() first?')

                #if fluxes[i,j] == 0 and fluxes[i+1,j] == 0: #add 0 flux
                #    ninc = ninc + 1
                #    continue

                if fluxes[i,j] == 0 or fluxes[i+1,j] == 0: #add 0 flux
                    ninc = ninc + 1
                    continue


                F1 = fluxes[i,j]
                F2 = fluxes[i+1,j]
                if fluxes[i,j] == 0 or fluxes[i,j] == None: #sub very small value for interpolation
                    sys.exit('from_differential_to_integral_flux: found bin '
                            'flux of zero. Should not happen here, bin [i,j] ['
                            + str(i) + ',' + str(j) + '].' )
                    #F1 = 1e-15
                if fluxes[i+1,j] == 0 or fluxes[i+1,j] == None:
                    sys.exit('from_differential_to_integral_flux: found bin '
                            'flux of zero. Should not happen here, bin [i+1,j] '
                            '['+ str(i+1) + ',' + str(j) + '].' )
                    #F2 = 1e-15 #Set to very small number for interpolation
                logF1 = np.log(F1)
                logF2 = np.log(F2)
                logE1 = np.log(bin_center[i])
                logE2 = np.log(bin_center[i+1])
                endE = bin_center[i+1]
                if i+1 == nbins-1:
                    endE = energy_bins[nbins-1][1] #extend to edge of last bin

                f = lambda x:exp(logF1
                            + (np.log(x)-logE1)*(logF2-logF1)/(logE2-logE1))
                startE = max(bin_center[i],min_energy)
                fint = scipy.integrate.quad(f,startE,endE)
                if math.isnan(fint[0]):
                    print("from_differential_to_integral_flux: flux integral"
                        "across bins is NaN. Setting to zero. Bin values are "
                        + str(F1) + ' and ' + str(F2))
                    fint = [0]
                if fint[0] < 1e-10:
                    fint = [0]
                sum_flux = sum_flux + fint[0]
                ninc = ninc + 1

        #if last bin is integral, add (HEPAD)
        if energy_bins[nbins-1][1] == -1 and fluxes[nbins-1,j] >= 0:
            sum_flux = sum_flux + fluxes[nbins-1,j]
            ninc = ninc + 1

        if ninc == 0:
            sum_flux = -1
        integral_fluxes.append(sum_flux)

    return integral_fluxes



def extract_integral_fluxes(fluxes, experiment, flux_type, flux_thresholds,
            energy_thresholds, energy_bins, options):
    """ Select or create the integral fluxes that correspond to the desired
        energy thresholds.
        If the user selected differential fluxes, then the
        differential fluxes will be converted to integral fluxes with the
        minimum energy defined by the set energy thresholds.
        If the user selected integral fluxes, then the channels corresponding to
        the desired energy thresholds will be identified.
        If an energy channel with a specified is not present, will return a
        None value for the integral fluxes, and remove that threshold in
        run_all. (e.g. if >100 MeV fluxes were not in the data or model
        fluxes, return a None value for integral fluxes and remove the
        >100 MeV, 1 pfu threshold.)
        
        INPUTS:
        
        :fluxes: (float nxm array) - flux time profile for each energy channel
        :experiment: (string)
        :flux_type: (string) - integral or differential
        :flux_thresholds: (float 1xp array) - flux threshold values
        :energy_thresholds: (float 1xp array) - energy channels to which the
            flux threshold values are applied
        :energy_bins: (float 1xn array) - bins for each energy channel in fluxes
        :options: (string array) - bg subtraction, effective energies for
            GOES channels
            
        OUTPUTS:
        
        :integral_fluxes: (float array) - integral fluxes associated with each
            energy channel in energy_thresholds; >energy_thresholds[i]
            
    """
    nthresh = len(flux_thresholds)
    nenergy = len(energy_bins)

    #Calculate integral fluxes corresponding the the energy_thresholds
    integral_fluxes = []
    if flux_type == "differential":
        for i in range(nthresh):
            integral_flux = from_differential_to_integral_flux(experiment, \
                            energy_thresholds[i], energy_bins, fluxes, options)
            if i == 0:
                integral_fluxes = [np.array(integral_flux)]
            else:
                integral_fluxes = np.concatenate((integral_fluxes,
                                        [np.array(integral_flux)]))
    if flux_type == "integral":
        for i in range(nthresh):
            found_thresh = False
            for j in range(nenergy):
                if energy_bins[j][0] == energy_thresholds[i]:
                    found_thresh = True
                    if len(integral_fluxes) == 0:
                        integral_fluxes = [np.array(fluxes[j,:])]
                    else:
                        integral_fluxes = np.concatenate((integral_fluxes,
                                                [np.array(fluxes[j,:])]))
            #if the energy channel is not present in the input fluxes
            if not found_thresh:
                print("Didn't find energy threshold for " \
                + str(energy_thresholds[i]) + ", " + str(flux_thresholds[i]))
                if len(integral_fluxes) == 0:
                    integral_fluxes = [None]*len(fluxes[0,:])
                else:
                    integral_fluxes = np.concatenate((integral_fluxes,
                                        [np.array([0]*len(fluxes[0,:]))]))

        if len(integral_fluxes) == 0:
            sys.exit("There were no integral fluxes with the requested energy "
                    "thresholds. Exiting.")
        if len(integral_fluxes) != len(energy_thresholds):
            sys.exit("Integral fluxes do not exist for some of the energy "
                    "thresholds. Exiting.")

    return integral_fluxes


def calculate_threshold_crossing(energy_threshold,flux_threshold,dates,fluxes):
    """ Calculate the time that a threshold is crossed.
        Operational thresholds used by the NASA JSC Space Radiation Analysis
        Group to determine actions that should be taken during an SEP event are:
        
            * >10 MeV proton flux exceeds 10 pfu (1/[cm^2 s sr])
            * >100 MeV proton flux exceeds 1 pfu (1/[cm^2 s sr])
            
        An SEP event is considered to start if 3 consecutive points are
        above threshold. The start time is set to the first point that crossed
        threshold.
        
        Values are calcualated for any additional thresholds provided by user.
        
        This program also calculates peak flux and time of peak flux for time
        period specified by user. It then calculates rise time by comparing:
        rise time = time of peak flux - threshold crossing time
        If no thresholds are crossed during specified time period, peak time
        and rise time will be 0.
        The event end time is much more subjective. For this program, I follow
        the code that officially calls the end of an event when SRAG supports
        ISS operations. When the flux has three consecutive points (15 minutes
        of GOES data) below 0.85*threshold, then the event is ended. The end time
        is corrected back to the first of the three data points. This end
        criteria has no physics basis. It simply ensures that a new event will
        not erroneously begin within the SRAG SEP alarm code due to flux
        fluctuations around threshold as the SEP flux decays away slowly.
        If the time period specified by the user ends before the event falls
        below 0.85*threshold, then the event_end_time is set to the last time
        in the specified time range. The program will give a message indicating
        that the user may want to extend the requested time frame in order to
        better-capture the end of the event.
        
        INPUTS:
        
        :energy_threshold: (float) - energy channel for which the
            flux threshold value should be applied
        :flux_threshold: (float) - flux threshold value
        :dates: (datetime 1xn array) - dates associated with the flux time profile
        :fluxes: (float 1xn array) - flux with time for the single energy channel
                associated with energy_threshold
        
        OUTPUTS:
        
        :crossing_time: (datetime)
        :peak_flux: (float) - maximum flux value between start and end time
        :peak_time: (datetime)
        :rise_time: (timedelta) - (peak_time - crossing_time)
        :event_end_time: (datetime)
        :duration: (timedelta) - (event_end_time - crossing_time)
        
    """
    #flux_threshold, e.g. 10 pfu (1/[cm^2 s sr])
    #energy_threshold, e.g. 10 MeV --> integral flux for >10 MeV protons
    #dates contain all of the datetimes for each data point
    #fluxes are a 1D array of integral fluxes (not multiple energy channels)
    print('Calculating threshold crossings and SEP event characteristics.')

    ndates = len(dates)
    end_threshold = vars.endfac*flux_threshold
                    #endfac = 0.85 used by SRAG operators;
                    #Can specify value in library/global_vars.py

    threshold_crossed = False
    event_ended = False
    peak_flux = 0
    peak_time = 0
    crossing_time = 0 #define in case threshold not crossed
    rise_time = 0 #define in case threshold not crossed
    event_end_time = 0
    duration = 0
    npoints = 3 #require 3 points above threshold
    tdiff = determine_time_resolution(dates)
    tdiff = tdiff.total_seconds()/(60) #time resolution of data set
    if tdiff > 15:
        npoints = 1 #time resolution >15 mins, require one point above threshold

    for i in range(ndates):
        if not threshold_crossed:
            if(fluxes[i] >= flux_threshold):
                start_counter = 0
                if i+(npoints-1) < ndates:
                    for ii in range(npoints):
                        if fluxes[i+ii] >= flux_threshold:
                            start_counter = start_counter + 1
                if start_counter == npoints:
                    crossing_time = dates[i]
                    threshold_crossed = True
        if threshold_crossed and not event_ended:
            if (fluxes[i] >= end_threshold):
                end_counter = 0  #reset if go back above threshold
            if (fluxes[i] <= end_threshold): #flux drops below 85% threshold
                end_counter = end_counter + 1
                if end_counter == npoints: #N consecutive points triggers end
                    event_ended = True
                    event_end_time = dates[i-(npoints-1)] #correct back time steps

    #In case that date range ended before fell before threshold,
    #use the last time in the file
    if crossing_time != 0 and event_end_time == 0:
        event_end_time = dates[ndates-1]
        print("!!!!File ended before SEP event ended for >"
            + str(energy_threshold) + ", " + str(flux_threshold) + " pfu! "
            "Using the last time in the date range as the event end time. "
            "Extend your date range to get an improved estimate of the event "
            "end time and duration.")

    if crossing_time != 0:
        #Find peak flux within event start and stop time
        for i in range(ndates):
            if dates[i] >= crossing_time and dates[i] <= event_end_time:
                if fluxes[i] > peak_flux:
                    peak_flux = fluxes[i]
                    peak_time = dates[i]
        if peak_time != 0:
            rise_time = peak_time - crossing_time
            duration = event_end_time - crossing_time

    return crossing_time,peak_flux,peak_time,rise_time,event_end_time,duration


def consistent_peak_fluxes(crossing_time, onset_peak, onset_date,
        peak_flux, peak_time):
    """ If a threshold is not crossed, then the max flux is saved
        as zero and the onset peak is saved as the maximum value
        in the time range considered. This is technically
        not correct.
        Swapping the two values to more correct.
        
        INPUTS (for n SEP events):
        
        :crossing_time: (datetime 1xn array) - value will be zero
            if no threshold was crossed
        :onset_peak: (float 1xn array) - value of onset peak for each
            SEP event, max flux if no thresholds crossed
        :onset_date: (datetime 1xn array) - time of onset peaks
        :peak_flux: (float 1xn array) - values of maximum flux for
            each SEP event, zero if no threshold crossed
        :peak_time: (datetime 1xn array) - value of peak flux, zero if
            threshold not crossed
            
        OUTPUTS (same as above, but values switched around):
        
        :onset_peak: value set to zero if no threshold crossed
        :onset_date: value set to zero if no threshold crossed
        :peak_flux: value of max flux during time period
        :peak_time: time of max flux
        
    """
    if len(crossing_time) != len(onset_peak):
        sys.exit("consistent_peak_flux: Different lengths for "
                "threshold crossings and onset peaks. Exiting.")

    for i in range(len(crossing_time)):
        if crossing_time[i] == 0:
            peak_flux[i] = onset_peak[i]
            peak_time[i] = onset_date[i]
            onset_peak[i] = 0
            onset_date[i] = 0
    
    return onset_peak, onset_date, peak_flux, peak_time


def check_bin_exists(threshold, energy_bins):
    """ If a user specifies a threshold for a differential energy bin, check
        to see that the energy bin is in the requested data set.
        threshold is a list of strings.
        Expect threshold[0] = lowedge-highedge, threshold[1] = flux value
        
        INPUTS:
        
        :threshold: (string) - threshold input by user in string format
        :energy_bins: (float 2xn array) - energy bins for n energy channels
        
        OUTPUTS:
        
        :Boolean: indicating specified bins exist, or exit
        
    """
    edges = threshold[0].split("-")
    lowedge = float(edges[0])
    highedge = float(edges[1])
    for bin in energy_bins:
        if lowedge == bin[0] and highedge == bin[1]:
            return True

    sys.exit("check_bin_exists: The user-defined threshold (" + threshold[0] \
            + ',' + threshold[1] + "cannot be found in the experiment's energy "
            "bins: " + str(energy_bins))



def get_energy_bin_index(energy_threshold, energy_bins):
    """ Find the index associated with the integral or differential energy bin, check
        to see that the energy bin is in the requested data set.
        threshold is a list of strings.
        For integral, expect:
            
            * energy_threshold = lowedge
        
        For differential, expect:
            
            * energy_threshold[0] = lowedge
            * energy_threshold[1] = highedge
            
        INPUTS:
        
        :energy_threshold: (float or 1x2 array)
        :energy_bins: (float 2xn array) - energy bins for n energy channels
        
        OUTPUTS:
        
        :index: (integer) - index of energy bin associated with energy_threshold
        
    """
    
    is_array = isinstance(threshold, list) #number or array?
    
    if not is_array: #integral
        for index in range(len(energy_bins)):
            if energy_threshold == energy_bins[index][0]:
                return index
    
    if is_array: #differential
        lowedge = threshold[0]
        highedge = threshold[1]
        for index in range(len(energy_bins)):
            if lowedge == energy_bins[index][0] \
                and highedge == energy_bins[index][1]:
                return index

    sys.exit("get_energy_bin_index: Could not identify bin that was "
            "requested: " + str(energy_threshold) +". Exiting.")
    
    
def determine_time_resolution(dates):
    """ The time resolution is found by taking the difference between
        every consecutive data point. The most common difference is
        taken as the time resolution. Even if the data set has gaps,
        if there are enough consecutive time points in the observational
        or model output, the correct time resolution should be identified.
        
        INPUTS:
        
        :dates: (datetime 1xm array) - dates associated with flux time profile
        
        OUTPUTS:
        
        :time_resolution: (time delta object)
        
    """
    ndates = len(dates)
    time_diff = [a - b for a,b in zip(dates[1:ndates],dates[0:ndates-1])]
    time_resolution = mode(time_diff)
    return time_resolution


def calculate_fluence(dates, flux):
    """ This subroutine sums up all of the flux in the 1xn array "flux". The
        "dates" and "flux" arrays input here should reflect only the intensities
        between the SEP start and stop times, determined by the subroutine
        calculate_threshold_crossing. The subroutine does not differentiate between
        differential or integral flux.
        
        The extract_date_range subroutine is used prior to calling this one to
        make the dates and fluxes arrays covering only the SEP time period.
        The flux will be multiplied by time_resolution and summed for all of the
        data between the start and end times. Negative or bad flux values
        should have been set to None or interpolated with check_bad_data().
        None values will be skipped.
        
        The time resolution is found by taking the difference between
        every consecutive data point. The most common difference is
        taken as the time resolution. Even if the data set has gaps,
        if there are enough consecutive time points in the observational
        or model output, the correct time resolution should be identified.
        
        Fluence units will be 1/[MeV cm^2] for GOES or SEPEM differential
        fluxes or 1/[cm^2] for integral fluxes.
        
        INPUTS:
        
        :flux: (float 1xn array) - intensity time series for a single energy bin
            or single integral channel (1D array).
        :dates: (datetime 1xn array) - datetimes that correspond to the fluxes
        
        OUTPUTS:
        
        :fluence: (float) - sum of all the flux values in flux
        
    """
    ndates = len(dates)
    time_resolution = determine_time_resolution(dates)
    
    #print("calculate_fluence: Identified a time resolution of "
    #        + str(time_resolution.total_seconds()) + " seconds.")
    
    fluence = 0
    for i in range(ndates):
        if flux[i] == None: continue
        if flux[i] >= 0:  #0 flux ok for models
                fluence = fluence + flux[i]*time_resolution.total_seconds()
        else:
            sys.exit('Bad flux data value of ' + str(flux[i]) +
                    ' found for bin ' + str(i) + ', '
                    + str(dates[i]) + '. This should not happen. '
                    + 'Did you call check_for_bad_data() first?')
                    
    fluence = fluence*4.0*math.pi #multiply 4pi steradians
    return fluence


def get_fluence_spectrum(experiment, flux_type, options, doBGSub,
                model_name, energy_threshold,
                flux_threshold, sep_dates, sep_fluxes, energy_bins,
                diff_thresh, save_file):
    """ Calculate the fluence spectrum for each of the energy channels in the
        user selected data set. If the user selected differential fluxes, then
        the fluence values correspond to each energy bin. If the user selected
        integral fluxes, then the fluence values correspond to each integral bin.
        Writes fluence values to file according to boolean save_file.
        If the user input a threshold for a differential energy bin,
        is_diff_thresh will be true and the filename will not contain "gt".
        
        INPUTS:
        
        :experiment: (string)
        :flux_type: (string) - integral or differential
        :options: (string array) - bg subtraction, effective energies for
            GOES channels
        :model_name: (string) - name of model or user experiment, if relevant
        :energy_threshold: (float) - energy channel to which the
            flux threshold value is applied
        :flux_threshold: (float) - flux threshold value
        :sep_dates: (datetime 1xm array) - dates trimmed between SEP start and
            end times
        :sep_fluxes: (float nxm array) - flux time profiles for each energy channel
            trimmed between SEP start and end times
        :energy_bins: (float 1xn array) - bins for each energy channel
        :diff_thresh: (boolean) - indicates if the energy_threshold and
            flux_threshold refer to a differential channel (True)
        :save_file: (boolean) - set True to save fluence values to file
        
        OUTPUTS:
        
        :fluence: (float 1xn array) - fluence value in each energy bin
        :energies: (float 1xn array) - bin centers for each energy bin
        
    """
    nenergy = len(energy_bins)
    fluence = np.zeros(shape=(nenergy))
    energies = np.zeros(shape=(nenergy))
    for i in range(nenergy):
        #Multiplied by 4pi sr (units of e.g. 1/[cm^2] or 1/[MeV cm^2])
        fluence[i] = calculate_fluence(sep_dates, sep_fluxes[i,:])
        if energy_bins[i][1] != -1:
            energies[i] = math.sqrt(energy_bins[i][0]*energy_bins[i][1])
        else:
            energies[i] = energy_bins[i][0]

    if save_file:
        #Write fluence to file
        year = sep_dates[0].year
        month = sep_dates[0].month
        day = sep_dates[0].day
        mod1 = 'gt'
        mod2 = '>'
        unit = flux_units_integral #'pfu'
        if diff_thresh:
            mod1 = ''
            mod2 = 'differential energy bin with low edge '
            unit = flux_units_differential #'1/[MeV cm^2 s sr]'
        modifier = ''
        if options[0] != '':
            for opt in options:
                modifier = modifier + '_' + opt
        if doBGSub:
            modifier = modifier + '_bgsub'
        foutname = outpath + '/fluence_' + str(experiment) + modifier + '_' \
                    + str(flux_type) + '_' + mod1 + str(energy_threshold) \
                    + '_' +str(year) + '_' + str(month) + '_' + str(day) +'.csv'
        if experiment == 'user' and model_name != '':
            foutname = outpath + '/fluence_' + model_name + modifier + '_' \
                        + str(flux_type) + '_' + mod1 + str(energy_threshold) \
                        + '_' +str(year) + '_' + str(month) + '_' + str(day) \
                        +'.csv'
        fout = open(foutname,"w+")

        fout.write('\"#Event defined by ' + mod2 + str(energy_threshold) \
                    + ' '+ energy_units + ', ' + str(flux_threshold) \
                    +' '+ unit + '; start time '
                    + str(sep_dates[0]) + ', end time '
                    + str(sep_dates[len(sep_dates)-1]) + '\"\n')
        if flux_type == "differential":
            fout.write("#Elow,Emid,Ehigh,Fluence " +
                        fluence_units_differential + "\n")
        if flux_type == "integral":
            fout.write("#>Elow,Fluence " + fluence_units_integral + "\n")

        for i in range(nenergy):
            if flux_type == "differential":
                fout.write("{0},{1},{2},{3}\n".format(energy_bins[i][0],
                        energies[i], energy_bins[i][1], fluence[i]))
            if flux_type == "integral":
                fout.write("{0},{1}\n".format(energy_bins[i][0], fluence[i]))
        fout.close()

    return fluence, energies


def calculate_event_info(energy_thresholds,flux_thresholds,dates,
                integral_fluxes, detect_prev_event, two_peaks, diff_thresh):
    """ Applies energy and flux thresholds to calculate SEP event quantities:
            
            * Threshold crossing time (onset)
            * Peak Flux in date range specified by user
            * Time of Peak Flux
            * Rise Time (onset to peak flux)
            * Event End Time (below 0.85*threshold for 3 data points, e.g 15 min)
            * Duration (onset to end)
            
        If the detect_prev_event flag is set to true and the threshold is crossed
        on the first time in the specified date range, this indicates that fluxes
        were already high due to a previous event. The code will look for the
        flux to drop below threshold and then increase above threshold again
        during the specified time period. If this flag is not set, then the code
        simply take the first threshold crossing as the start of the SEP event.
        
        If the two_peaks flag is set to true, the code will use the first
        identified threshold crossing as the start time. A second threshold
        crossing will be allowed. The event end will be determined as the drop
        below threshold after the second threshold crossing. The duration will
        be the end time - first threshold crossing time. The peak will be
        determined as the largest flux value found in both threshold crossings.
        The peak and rise times will be calculated from the first threshold
        crossing time.
        
        NOTE: integral_fluxes is a "historic" name and may actual contain
        integral or differential fluxes. They correspond to the flux channels
        for which a threshold has been applied.
        
        INPUTS:
        
        :energy_thresholds: (float 1xn array) - energy channels for which a
            threshold is applied
        :flux_thresholds: (float 1xn array) - flux thresholds to apply
        :dates: (datetime 1xm array)
        :integral_fluxes: (float nxm array) - fluxes only associated with the
            energy channels in energy_thresholds
        :detect_prev_event: (boolean)
        :two_peaks: (boolean)
        :diff_thresh: (boolean) - indicates if the energy_threshold and
            flux_threshold refer to a differential channel (True)
            
        OUTPUTS:
        
        :crossing_time: (datetime 1xn array) - crossing times for n applied thresholds
        :peak_flux: (float 1xn array) - max flux value between start and end time
        :peak_time: (datetime 1xn array)
        :rise_time: (timedelta 1xn array) - (peak_time - crossing_time)
        :event_end_time: (datetime 1xn array)
        :duration: (timedelta 1xn arrat) - (event_end_time - crossing_time)
        
    """
    nthresh = len(flux_thresholds)
    crossing_time = []
    peak_flux = []
    peak_time = []
    rise_time = []
    event_end_time = []
    duration = []
    for i in range(nthresh):
        ct,pf,pt,rt,eet,dur = calculate_threshold_crossing(energy_thresholds[i],
                        flux_thresholds[i],dates,integral_fluxes[i])
        if detect_prev_event and ct == dates[0]:
            print("Threshold may have been high due to previous event."
                "Recalculating event info for remaining time period in data "
                "set.")
            last_date = dates[len(dates)-1]
            tmp_dates, tmp_fluxes = datasets.extract_date_range(eet,last_date,
                                dates,integral_fluxes)
            ct,pf,pt,rt,eet,dur = calculate_threshold_crossing(\
                            energy_thresholds[i],flux_thresholds[i],
                            tmp_dates,tmp_fluxes[i])

        if dur != 0 and two_peaks:
            if dur < timedelta(days=1):
                print("User specified that event has two peaks. Extending "
                    "event to second decrease below threshold.")
                last_date = dates[len(dates)-1]
                tmp_dates, tmp_fluxes = datasets.extract_date_range(eet,
                                    last_date,dates,integral_fluxes)
                ct2,pf2,pt2,rt2,eet2,dur2 = calculate_threshold_crossing(\
                                energy_thresholds[i],flux_thresholds[i],
                                tmp_dates,tmp_fluxes[i])

                #Only apply changes if another threshold crossing was found and
                #new event end time is not the last point in the file
                if dur2 != 0:
                    if eet2 < dates[len(dates)-1]:
                        dur = eet2 - ct #adjust duration to include both peaks
                        eet = eet2 #adjust event end time
                        if pf2 > pf: #determine which peak has the highest flux
                            pf = pf2
                            pt = pt2
                            rt = pt2 - ct

        crossing_time.append(ct)
        peak_flux.append(pf)
        peak_time.append(pt)
        rise_time.append(rt)
        event_end_time.append(eet)
        duration.append(dur)
        mod = '>'
        units = flux_units_integral #'pfu'
        if diff_thresh:
            mod = ''
            units = flux_units_differential #'1/[MeV cm^2 s sr]'
        print(
               'Flux      Threshold    Time Crossed         Peak Flux'
                + '            Peak Time' + '            Rise Time'
                + '  End Time' + '            Duration'
             )
        print(
                mod + str(energy_thresholds[i]) + ' '
                + energy_units + '   '
                + str(flux_thresholds[i]) + ' ' + units + '       '
                + str(crossing_time[i]) + '  '+ str(peak_flux[i]) + '   '
                + str(peak_time[i]) + '  ' + str(rise_time[i]) + '    '
                + str(event_end_time[i]) + '  ' + str(duration[i])
             )
    return crossing_time,peak_flux,peak_time,rise_time,event_end_time,duration


def calculate_onset_peak(experiment, energy_thresholds, dates, integral_fluxes,
                crossing_time, event_end_time, showplot):
    """ Calculate the peak associated with the initial SEP onset. This subroutine
        searches for the rollover that typically occurs after the SEP onset.
        The peak value will be specified as the flux value at the rollover
        location.
        If the event peak is larger than 500/energy channel (e.g. for >10
        MeV, if peak is above 50 pfu), then the rollover will be found
        starting from the threshold crossing time.
        If the event is a "weaker" event with a maximum below the value above,
        the rollover location will be identified starting from 12 hours (buffer)
        prior to threshold crossing time. The onset peak time may be BEFORE
        THE EVENT START TIME with this logic, but is more independent of the
        level of the threshold and more linked to the behavior of the flux
        time profile.
        The onset peak may provide a more physically appropriate comparison
        with models.
        If code cannot identify onset peak, it will return a value of 0 on
        the date None.
        
        INPUTS:
        
        :experiment: (string) e.g. GOES-13
        :energy_thresholds: (float 1xn array) - energy channels for which thresholds
            are applied
        :dates: (datetime 1xm array) - dates associated with flux time profile
        :integral_fluxes: (float nxm array) - fluxes for each energy channel for
            which a threshold is applied; each is the same length as dates
        :crossing_time: (datetime 1xn array) - threshold crossing times for each energy
            channel for which a threshold is applied
        :event_end_time: (datetime 1xn array) - end times for each energy channel for which
            a threshold is applied
        :showplot: (bool)
        
        OUTPUTS:
        
        :onset_date: (datetime 1xn array) - time of onset peak
        :onset_peak: (float 1xn array) - flux value of onset peak
        
    """
    nthresh = len(energy_thresholds)
    smooth_flux = [[]]*nthresh
    for i in range(nthresh):
        #tried out a smoothing algorithm, but abandoned it
        smooth_flux[i] = integral_fluxes[i]


    run_deriv = [[]]*nthresh
    run_deriv_norm = [[]]*nthresh
    nwin = 8 #Number of points away for calculating derivative, 5 min data
    time_resolution = determine_time_resolution(dates)
    if time_resolution > datetime.timedelta(minutes=5):
        nwin = 1
    for i in range(nthresh):
        if crossing_time[i] == 0: #no event, no threshold crossed
            continue
        zeroes = False #Models may output zero flux
        if 0 in smooth_flux[i]: zeroes = True
        if None in smooth_flux[i]: zeroes = True
        
        #Normalize the derivative by the flux value
        #If there are zeroes or None values in the array,
        #then normalize by a constant flux value
        #Find the first positive flux value in the array
        #Most likely smooth_flux[i][0]
        normindx = np.nonzero(smooth_flux[i])
        normfac = np.amax(smooth_flux[i][normindx[0]])

        run_deriv[i] = [0] #normalized by constant factor
        run_deriv_norm[i] = [0] #normalized by changing flux
        for j in range(1,nwin):
            deriv = smooth_flux[i][j] - smooth_flux[i][0]
            run_deriv[i].append(deriv/normfac)
            #run_deriv_norm[i].append(deriv/normfac)
            if not zeroes:
                #normalize difference by the flux
                run_deriv_norm[i].append(deriv/smooth_flux[i][0])
            else:
                #use difference directly
                run_deriv_norm[i].append(deriv/normfac)
                
        for j in range(nwin,len(smooth_flux[i])):
            deriv = smooth_flux[i][j] - smooth_flux[i][j-nwin]
            run_deriv[i].append(deriv/normfac)
            #run_deriv_norm[i].append(deriv/normfac)
            if not zeroes:
                run_deriv_norm[i].append(deriv/smooth_flux[i][j-nwin])
            else:
                run_deriv_norm[i].append(deriv/normfac)


    #Search for onset peak
    #Peak likely occurs near first time derivative goes from generally positive
    #to zero
    onset_date = [[]]*nthresh
    onset_peak = [[]]*nthresh
    for i in range(nthresh):
        if crossing_time[i] == 0:
            #threshold not crossed
            #Find the max flux in the time period and time
            onset_peak[i] = np.max(smooth_flux[i])
            peak_index = np.where(smooth_flux[i] == np.amax(smooth_flux[i]))
            onset_date[i] = dates[peak_index[0][0]]
            continue
        #If all entries are zero flux (bg-sub or SEPEMv3)
        is_all_zero = np.all((smooth_flux[i] == 0))
        if is_all_zero:
            #No non-zero fluxes
            onset_peak[i] = 0.
            onset_date[i] = None
            continue

        #Get value of maximum positive derivative in first 18 hours
        #record where deriv first goes negative
        #For an event that just barely crosses threshold, look for
        #onset peak prior to threshold crossing
        #For very large event, begin looking after threshold crossing
        buffer = 0
        if np.amax(smooth_flux[i]) < 500/energy_thresholds[i]:
            buffer = 12
        
        index_cross = 0
        index_stp = 0
        last_date = crossing_time[i] + datetime.timedelta(hours=18)
        if last_date > event_end_time[i]:
            last_date = event_end_time[i]
        for j in range(len(dates)):
            if dates[j] <= (crossing_time[i] - datetime.timedelta(hours=buffer)):
                index_cross = j
            if dates[j] <= last_date:
                index_stp = j

        #Max value of the normalized derivative in first 18 hours
        max_index =  np.argmax(run_deriv_norm[i][index_cross:index_stp]) + index_cross
        
        
#        #Find first index where derivative drops below a small negative number
#        index_neg = next((i for i, x in enumerate(run_deriv_norm[i][max_index:index_stp+1]) if x < -0.01), None)
#        if index_neg == None: #No negative value found
#                print("calculate_onset_peak: Could not locate onset peak for "
#                    +  str(energy_thresholds[i])
#                    + " MeV. Setting onset peak and date to None.")
#                onset_peak[i] = None
#                onset_date[i] = None
#                continue #go to next threshold
#        else:
#            index_neg = index_neg + max_index
#            #Find the maximum flux in the range where the onset peak should be
#            onset_peak[i] = max(integral_fluxes[i][max_index:index_neg+1])
#            onset_index = np.argmax(integral_fluxes[i][max_index:index_neg+1])
#            onset_date[i] = dates[max_index + onset_index]
#            print("Found onset peak for " + str(energy_thresholds[i]) \
#                + " MeV: " + str(onset_peak[i]) + ", Onset peak time: " \
#                + str(onset_date[i]))
#            #Need index where derivative first drops below zero
#            #for checks below
#            index_zero = next((i for i, x in enumerate(run_deriv_norm[i][max_index:index_stp+1]) if x < 0), None)
#            index_zero = index_zero + max_index
#
#        #Find the next positive value and then next negative value
#        #Idea is to get the next bump in derivative and compare
#        #this next rise with the previous one to see if it is
#        #similar.
#        #Keep searching every dip as needed
#        index_neg_next = index_neg
#        while index_neg_next < index_stp:
#            index_pos = next((i for i, x in enumerate(run_deriv_norm[i][index_neg_next:index_stp+1]) if x > 0), None)
#            if index_pos == None:
#                break
#            else:
#                index_pos = index_pos + index_neg_next
#
#            index_neg_next = next((i for i, x in enumerate(run_deriv_norm[i][index_pos:index_stp+1]) if x < -0.01), None)
#            if index_neg_next == None:
#                break
#            else:
#                index_neg_next = index_neg_next + index_pos
#
#            #Now have indices bounding a "bump" of positive derivative
#            #Find the max value and location of that value
#            max_index_next =  np.argmax(run_deriv_norm[i][index_pos:index_neg_next+1]) + index_pos
#            npts = max_index_next - index_pos
#            if npts == 0: continue
#            #Go the same number of points backwards in time to compare the
#            #delta in the intial derivative with the new "bump" in the
#            #derivative to see if the flux is continuing to increase
#            #at a similar rate
#            ave_deriv_init = sum(run_deriv_norm[i][index_zero - npts:index_zero+1])/npts
#            ave_deriv_next = sum(run_deriv_norm[i][index_pos:max_index_next+1])/npts
#            deriv_diff = abs(ave_deriv_init - ave_deriv_next)/ave_deriv_init
#            print("ave_deriv_init: "+str(ave_deriv_init)+" ave_deriv_next: "\
#                +str(ave_deriv_next)+" deriv_diff: "+str(deriv_diff)\
#                +" n pts: " + str(npts) \
#                +" starting at date " + str(dates[index_pos]))
#
#
#            if deriv_diff <= 0.1 \
#                or (ave_deriv_next>0.1 and ave_deriv_next>ave_deriv_init):
#                onset_peak[i] = max(integral_fluxes[i][index_neg:index_neg_next])
#                onset_index = np.argmax(integral_fluxes[i][index_neg:index_neg_next])
#                onset_date[i] = dates[index_neg + onset_index]
#                print("Recalculated onset peak for " + str(energy_thresholds[i]) \
#                        + " MeV: " + str(onset_peak[i]) + ", Onset peak time: " \
#                        + str(onset_date[i]))
#            else:
#                break #don't keep looking for onset peak
        
        
        
        
        
        #Find where derivative falls below zero after the max
        index_neg = 0
        first_neg = False
        deriv_thresh = -0.04
        for j in range(max_index,index_stp+1):
            if run_deriv_norm[i][j] < deriv_thresh and not first_neg:
                index_neg = j
                first_neg = True
        while not first_neg:
            deriv_thresh = round(deriv_thresh + 0.01, 2) #negative value
            for j in range(max_index,index_stp+1):
                if run_deriv_norm[i][j] < deriv_thresh and not first_neg:
                    index_neg = j
                    first_neg = True
            if deriv_thresh > 0:
                print("calculate_onset_peak: Could not locate onset peak for "
                    +  str(energy_thresholds[i]) + " " + energy_units
                    + ". Setting onset peak and date to None.")
                onset_peak[i] = None
                onset_date[i] = None
                first_neg = True #exit loop

        if onset_peak[i] == None:
            continue
        #Find the maximum flux in the range where the onset peak should be
        onset_peak[i] = max(integral_fluxes[i][max_index:index_neg])
        onset_index = np.argmax(integral_fluxes[i][max_index:index_neg])
        onset_date[i] = dates[max_index + onset_index]
        print("Found onset peak for " + str(energy_thresholds[i]) + " "
            + energy_units + ": " + str(onset_peak[i])
            + ", Onset peak time: " + str(onset_date[i]))

        #Check and see if the event continues rising at a similar rate to the
        #onset peak. If so, it's likely a little dip on the way up.
        #Check the average derivative 1 hour prior to the onset peak
        #and the average 1 hour after the peak. If similar, likely event is
        #continuing to rise.
        time_res = determine_time_resolution(dates)
        npts = math.ceil(50.*60./time_res.total_seconds())
        stpt = index_neg - npts
        if stpt < 0: stpt = 0
        endpt = index_neg + npts
        if endpt >= index_stp: #don't check past 18 hours into event
            continue
        if endpt >= len(dates): endpt = len(dates)-1
        #Use both derivatives to compare pre and post rise
        deriv_ave_pre = sum(run_deriv[i][stpt:index_neg])/len(run_deriv[i][stpt:index_neg])
        deriv_ave_post =  sum(run_deriv[i][index_neg:endpt])/len(run_deriv[i][index_neg:endpt])
        deriv_diff = (deriv_ave_pre - deriv_ave_post)/deriv_ave_pre
        deriv_diff = abs(deriv_diff)

        deriv_ave_pre_norm = sum(run_deriv_norm[i][stpt:index_neg])/len(run_deriv_norm[i][stpt:index_neg])
        deriv_ave_post_norm =  sum(run_deriv_norm[i][index_neg:endpt])/len(run_deriv_norm[i][index_neg:endpt])
        deriv_diff_norm = (deriv_ave_pre_norm - deriv_ave_post_norm)/deriv_ave_pre_norm
        deriv_diff_norm = abs(deriv_diff_norm)

        #print("deriv_ave_pre: "+str(deriv_ave_pre)+" deriv_ave_post: "\
        #    +str(deriv_ave_post)+" deriv_diff: "+str(deriv_diff)\
        #    +" starting at date " + str(dates[index_neg]))
        #print("deriv_ave_pre_norm: "+str(deriv_ave_pre_norm)+" deriv_ave_post_norm: "\
         #   +str(deriv_ave_post_norm)+" deriv_diff_norm: "+str(deriv_diff_norm)\
         #   +" starting at date " + str(dates[index_neg]))
        #if the two differ by 10% or less or deriv is bigger
        #going forward
        if deriv_diff <= 0.1 or deriv_diff_norm <= 0.1\
            or (deriv_ave_post>0.1 and deriv_ave_post>deriv_ave_pre)\
            or (deriv_ave_post_norm>0.1 and deriv_ave_post_norm>deriv_ave_pre_norm):
#        if deriv_diff_norm <= 0.1\
#            or (deriv_ave_post_norm>0.1 and deriv_ave_post_norm>deriv_ave_pre_norm):
            #Max value of the normalized derivative in first 18 hours
            max_index =  np.argmax(run_deriv_norm[i][index_neg:index_stp]) + index_neg
            #Find where derivative falls below zero after the max
            index_neg2 = index_neg
            first_neg = False
            for j in range(max_index,index_stp+1):
                if run_deriv_norm[i][j] < 0 and not first_neg:
                    index_neg2 = j
                    first_neg = True
            #Find the maximum flux in the range where the onset peak should be
            if index_neg2 > index_neg:
                check_peak = max(integral_fluxes[i][index_neg:index_neg2])
                if check_peak > onset_peak[i]:
                    onset_peak[i] = max(integral_fluxes[i][index_neg:index_neg2])
                    onset_index = np.argmax(integral_fluxes[i][index_neg:index_neg2])
                    onset_date[i] = dates[index_neg + onset_index]
                    print("Recalculated onset peak for " + str(energy_thresholds[i]) + " " + energy_units
                        + ": " + str(onset_peak[i]) + ", Onset peak time: "
                        + str(onset_date[i]))


    if showplot:
        for i in range(nthresh):
            if crossing_time[i] == 0:
                continue

            fig, ax = plt.subplots(figsize=(9, 4))
            plt.title("Onset Peak Derivation for " + str(energy_thresholds[i]) \
                    + " MeV")
            ax.plot_date(dates,integral_fluxes[i],'-',label=experiment)
            ax.plot_date(dates,smooth_flux[i],'-',color="red",label="Smoothed")
            if onset_peak[i] != None:
                ax.plot_date(onset_date[i],onset_peak[i],'o',color="black")
            plt.yscale("log")
            ax.set_xlabel("Date")
            ax.set_ylabel("Flux")

            ax2 = ax.twinx()
            ax2.plot_date(dates,run_deriv[i],'-',color="blue",\
                        label="Derivative")
            ax2.plot_date(dates,run_deriv_norm[i],'-',color="green",\
                        label="Normalized Derivative")
            ax2.axhline(0,color='red',linestyle=':')
            ax2.set_ylabel("Derivative")

    return onset_date, onset_peak


def calculate_umasep_info(energy_thresholds,flux_thresholds,dates,
                integral_fluxes, crossing_time):
    """ Uses the integral fluxes (either input or estimated from differential
        channels) and all the energy and flux thresholds set in the main program
        to calculate SEP event quantities specific to the UMASEP model.
            Flux at threshold crossing time + 3, 4, 5, 6, 7 hours
            
        INPUTS:
        
        :energy_thresholds: (float 1xn array) - energy channels for which thresholds
            are applied
        :flux_thresholds: (float 1xn array) - flux thresholds that are applied
        :dates: (datetime 1xm array) - dates associated with flux time profile
        :integral_fluxes: (float nxm array) - fluxes for each energy channel for
            which a threshold is applied; each is the same length as dates
        :crossing_time: (datetime 1xn array) - threshold crossing times for each energy
            channel for which a threshold is applied
            
        OUTPUTS:
        
        :proton_delay_times: (datetime nx5 array) - times 3, 4, 5, 6, 7 hours
            after crossing time for n thresholds
        :proton_flux: (float nx5 array) - value of flux at each delay time and for
            each threshold
        
    """
    nthresh = len(flux_thresholds)
    proton_flux = []
    delays = [datetime.timedelta(hours=3), datetime.timedelta(hours=4),
                    datetime.timedelta(hours=5), datetime.timedelta(hours=6),
                    datetime.timedelta(hours=7)]
    ndelay = len(delays)
    delay_times = []
    proton_delay_times = []  #actual time point corresponding to flux

    #Match the correct time delay with the correct threshold
    for i in range(nthresh):
        if crossing_time[i] == 0:
            delay_times.append(0)
            continue
        all_delays = []
        for delay in delays:
            all_delays.append(crossing_time[i] + delay)
            #Make sure that delayed time doesn't exceed input time range
            if crossing_time[i] + delay > dates[len(dates)-1]:
                sys.exit("An UMASEP delayed time (Ts+3, 4, 5, 6, 7 hrs) "
                        "exceeded the user's input time range. Please rerun "
                        "and extend end time.")

        delay_times.append(all_delays) #all delays for a given threshold

    for i in range(nthresh):
        save_flux = [0]*ndelay
        save_dates = [0]*ndelay
        if delay_times[i] == 0:
            proton_delay_times.append(0)
            proton_flux.append(0)
            continue
        for k in range(ndelay):
            save_index = -1
            for j in range(len(dates)):
                if dates[j] <= delay_times[i][k]:
                    save_index = j

            #GET FLUX AT DELAYED TIME WITH 10 MINUTE AVERAGE
            #May choose to modify if input data set has something other than
            #5 minute time cadence.
            if save_index == -1: #should not happen, unless dates has no length
                sys.exit("Did not find an appropriate UMASEP flux point. "
                        "Exiting.")
            if save_index == 0: #also should be no way for this to happen
                save_flux[k] = (integral_fluxes[i][save_index] + \
                            integral_fluxes[i][save_index +1])/2.
            else:
                save_flux[k] = (integral_fluxes[i][save_index] + \
                            integral_fluxes[i][save_index - 1])/2.
            save_dates[k] = dates[save_index]

        proton_flux.append(save_flux)
        proton_delay_times.append(save_dates)

    return proton_delay_times, proton_flux



def report_threshold_fluences(experiment, flux_type, model_name,
                energy_thresholds, energy_bins, sep_dates, sep_fluxes):
    """ Report fluences for specified thresholds, typically >10, >100 MeV.
        These values are interesting to use for comparison with literature and
        for quantifying event severity.
        
        Assumes that energy_thresholds are referring to integral energy
        channels.
        Assumes that sep_fluxes are integral or estimated integral
        fluxes for each energy channel in energy_thresholds.
        If the input fluxes were differential, then this subroutine only
        works correctly for integral energy_thresholds and matching sep_fluxes
        that were estimated for each integral channel in energy_thresholds.
        
        INPUTS:
        
        :experiment: (string)
        :flux_type: (string) - integral or differential
        :model_name: (string) - name of model or experiment if experiment = "user"
        :energy_thresholds: (float 1xn array) - integral energy channels for
            which a threshold has been applied
        :energy_bins: (float 2xm array) - energy bins for input fluxes
        :sep_dates: (datetime 1xp array) - time profile dates trimmed to
            start and stop of SEP event
        :sep_fluxes: (float mxp array) - flux time profiles for each energy
            channel trimmed to start and end dates
            
        OUTPUTS:
        
        :integral_fluence: (float 1xn array) - event fluence for all
            integral energy channels for which a threshold was applied
    """
    tmp_energy_bins = []
    nthresh = len(energy_thresholds)
    ndates = len(sep_dates)
    sep_integral_fluxes = np.zeros(shape=(nthresh,ndates))
    for i in range(nthresh):
        tmp_energy_bins.append([energy_thresholds[i],-1])
        if flux_type == "differential":
            #integral fluxes were estimated for only the threshold energies
            #so the indices in the fluxes match the order of the thresholds
            sep_integral_fluxes[i,:] = sep_fluxes[i,:]
        if flux_type == "integral":
            #Pull out the integral fluxes from the correct energy channels
            #corresponding to the thresholds that were applied
            #This is selecting from all integral channels in the input flux file
            for j in range(len(energy_bins)):
                if energy_bins[j][0] == energy_thresholds[i]:
                    sep_integral_fluxes[i,:] = sep_fluxes[j,:]

    #Generate integral fluence spectrum for only the integral channels
    #If input fluxes in sep_fluxes are differential, convert them to
    #estimated integral fluxes in get_fluence_spectrum
    #that were specified with thresholds and in the same index order
    #Returns fluences multiplied by 4pi
    #get_fluence_spectrum(experiment, flux_type, options, doBGSub,
    #            model_name, energy_threshold,
    #            flux_threshold, sep_dates, sep_fluxes, energy_bins,
    #            diff_thresh, save_file)
    #No output file is written, so only the sep_dates, sep_integral_fluxes,
    #and tmp_energy_bins arguments are used for anything
    integral_fluence, integral_energies = get_fluence_spectrum(experiment,
                    "integral", '', False, #Filler values for filename b/c file not saved
                    model_name, 0, 0,sep_dates, sep_integral_fluxes,
                    tmp_energy_bins, False, False) #diff_thresh; savefile

    return integral_fluence


def save_integral_fluxes_to_file(experiment, flux_type, options, doBGSub,
        model_name, energy_thresholds, crossing_time, dates, integral_fluxes):
    """Output the time series of integral fluxes to a file. If the input
        data set was in integral channels, then this file will contain exactly
        the same values in the time series.
        If the input data set was in differential energy bins, then this file
        contains the estimated integral fluxes calculated in this program.
        
        Writes out csv file with dates in first column and integral fluxes
        for which a threshold was applied in the remaining columns. Writes
        fluxes for full time period specified by user.
        Any differential channels for which a threshold was applied are not
        included in this file.
        
        INPUTS:
        
        :experiment: (string)
        :flux_type: (string) - integral or differential
        :options: (string array) - S14, Bruno2017, uncorrected options for GOES data
        :doBGSub: (boolean) - indicates if background subtration is performed
        :model_name: (string) - name of model or user experiment, if relevant
        :energy_thresholds: (float 1xn array) - energy channels for which
            flux thresholds are applied
        :crossing_time: (datetime 1xn array) - start times of sep event for each
            energy channel (energy_thresolds) for which a threshold was applied
        :dates: (datetime 1xm array) - dates for flux time profile
        :integral_fluxes: (float nxm array) - flux time profiles for each energy channel
            for which a threshold was applied; assumed to be (estimated)integral fluxes
        
        OUTPUTS:
        
        No outputs except output file named e.g.
            integral_fluxes_GOES-13_differential_2012_3_7.csv
    """
    nthresh = len(energy_thresholds)
    ndates = len(dates)
    year = 0
    month = 0
    day = 0
    for i in range(nthresh):
        if crossing_time[i] != 0 and year == 0:
            year = crossing_time[i].year
            month = crossing_time[i].month
            day = crossing_time[i].day
    if year == 0:
        print("No thresholds were crossed during this time period. "
                "Integral flux time profiles not written to file. Exiting.")
        return

    #e.g. integral_fluxes_GOES-13_differential_2012_3_7.csv
    modifier = ''
    if options[0] != '':
        for opt in options:
            modifier = modifier + '_' + opt
    if doBGSub:
        modifier = modifier + '_bgsub'

    foutname = outpath + '/integral_fluxes_' + experiment + modifier + '_' \
                 + flux_type + '_' + str(year) + '_' + str(month) \
                 + '_' + str(day) + '.csv'
    if experiment == 'user' and model_name != '':
        foutname = outpath + '/integral_fluxes_' + model_name + modifier + '_' \
                     + flux_type + '_' + str(year) + '_' + str(month) \
                     + '_' + str(day) + '.csv'
    print('Writing integral flux time series to file --> ' + foutname)
    fout = open(foutname,"w+")
    if flux_type == "integral":
        fout.write('#Integral fluxes in units of '
                    + flux_units_integral + '\n')
    if flux_type == "differential":
        fout.write('#Estimated integral fluxes in units of '
                    +flux_units_differential + '\n')

    fout.write('#Columns headers indicate low end of integral channels in '
                    + energy_units + ';'
                    ' e.g. >10 ' + energy_units +'\n')
    fout.write('#Date')
    for thresh in energy_thresholds: #build header
        fout.write(',' + str(thresh))
    fout.write('\n')
    for i in range(ndates):
        fout.write(str(dates[i]))
        for j in range(nthresh):
            fout.write(',' + str(integral_fluxes[j][i]))
        fout.write('\n')

    fout.close()




def print_values_to_file(experiment, flux_type, options, doBGSub,
                model_name, startdate, energy_thresholds,
                flux_thresholds, crossing_time, onset_peak, onset_date,
                peak_flux, peak_time, rise_time, event_end_time, duration,
                threshold_fluences, is_diff_thresh, umasep, umasep_times,
                umasep_fluxes):
    """ Write all calculated values to file for all thresholds. Event-integrated
        fluences for >10, >100 MeV (and user-defined threshold) will also be
        included. Writes out file with name e.g.
        output/sep_values_experiment_fluxtype_YYYY_M_D.csv
        If the UMASEP option was selected, add on the proton values calculated
        at the UMASEP Ts + Xhr time points.
        is_diff_thresh indicates whether the user input a differential
        threshold. If so, the user threshold bin(s) will refer to differential
        fluxes.
        
        Writes all values out to a csv file. These values and more are also
        written out to json file in a different subroutine.
        
        INPUTS:
        
        :experiment: (string)
        :flux_type: (string) - integral or differential
        :options: (string array) - can be S14, Bruno2017, uncorrected and apply
            to GOES data
        :doBGSub: (boolean) - indicates if background subtraction performed
        :model_name: (string) - model or experiment name if experiment = "user"
        :startdate: (datetime) - start date of time period entered by user
        :energy_thresholds: (float 1xn array) - all energy channels for which
            a flux threshold was applied (both integral and differential)
        :flux_thresholds: (float 1xn array) - flux thresholds for each of the
            energy channels in energy_thresholds
        :crossing_time: (datetime 1xn array) - SEP event start time for each
            applied threshold
        :onset_peak: (float 1xn array) - onset peak for each applied threshold
        :onset_date: (datetime 1xn array) - onset peak time for each applied
            threshold
        :peak_flux: (float 1xn array) - maximum flux for each applied threshold
        :peak_time: (datetime 1xn array) - time of maximum flux for each applied
            threshold
        :rise_time: (timedelta 1xn array) - start time to max flux
        :event_end_time: (datetime 1xn array) - end time for each applied threshold
        :duration: (timedelta 1xn array) - end time minus start time
        :threshold_fluences: (float 1xn array) - event fluence for each energy
            channel in energy_thresholds
        :is_diff_thresh: (bool 1xn array) - indicates if threshold applied to
            integral or differential channel, e.g. if channel in energy_thresholds
            is an integral or differential channel
        :umasep: (boolean) - indicate if user called UMASEP flag
        :umasep_times: (datetime 5xn array) - times 3, 4, 5, 6, 7, 9 hours after
            crossing time for each applied threshold
        :umasep_fluxes: (float 5xn) - fluxes at each of those times for each
            applied threshold
            
        OUTPUTS:
        
        :year: (integer) - year of either SEP event or startdate (if no
            thresholds crossed)
        :month: (integer) - month of either SEP event or startdate (if no
            thresholds crossed)
        :day: (integer) - day of either SEP event or startdate (if no
            thresholds crossed)
        :Output file: named e.g.
            output/sep_values_experiment_fluxtype_YYYY_M_D.csv
        
        
    """
    nthresh = len(energy_thresholds)
    numa = 0
    if umasep:
        for i in range(nthresh):
            if umasep_times[i] == 0:
                continue
            if len(umasep_times[i]) > numa:
                numa = len(umasep_times[i])

    year = 0
    month = 0
    day = 0
    for i in range(nthresh):
        if crossing_time[i] != 0 and year == 0:
            year = crossing_time[i].year
            month = crossing_time[i].month
            day = crossing_time[i].day
    if year == 0:
        year = startdate.year
        month = startdate.month
        day = startdate.day
        print("No thresholds were crossed during this time period. "
            "Using starting date of time period for json file and not "
            "writing out a csv file.")
        return year, month, day, False

    modifier = ''
    if options[0] != '':
        for opt in options:
            modifier = modifier + '_' + opt
    if doBGSub:
        modifier = modifier + '_bgsub'

    foutname = outpath + '/sep_values_' + experiment + modifier + '_' \
                + flux_type + '_' + str(year) + '_' + str(month) + '_' \
                + str(day) +'.csv'
    if experiment == 'user' and model_name != '':
        foutname = outpath + '/sep_values_' + model_name + modifier + '_' \
                    + flux_type + '_' + str(year) + '_' + str(month) + '_' \
                    + str(day) +'.csv'
    print('Writing SEP values to file --> ' + foutname)
    fout = open(foutname,"w+")
    #Write header
    fout.write('#For thresholds that depend on integral fluxes (annotated '
            'with >) Flux Threshold, Onset Peak Flux, and Max Peak Flux '
            'have units of [' + flux_units_integral
            +'] and Bin Fluence has units of ['
            + fluence_units_integral + '].\n')
    fout.write('#For thresholds that depend on differential fluxes (no '
            '>) Flux Threshold, Onset Peak Flux and Max Peak Flux '
            'have units of [' + flux_units_differential
            + '] and Bin Fluence has units of ['
            + fluence_units_differential + '].\n')
    if flux_type == "differential":
        fout.write('#For thresholds that depend on integral fluxes (annotated '
                    'with >) - differential fluxes were converted '
                    'to integral fluxes. Onset Peak Flux and Max Flux are '
                    'estimated integral flux values.\n')
        if True in is_diff_thresh:
            fout.write('#The bottom threshold(s) are differential thresholds '
                        'using the energy bin with low edge specified in the '
                        'Energy Threshold column and differential flux in the '
                        'Flux Threshold column to define the SEP quantities.\n')
            
    fout.write('#Energy Threshold [' + energy_units
            + '],Flux Threshold,Start Time,Onset Peak Flux,Onset Time,'
            'Max Flux,Max Time,Rise Time,End Time,Duration, Bin Fluence')

    if umasep:
        for jj in range(numa):
            fout.write(',UMASEP Delay [hr],Flux ['
                        + flux_units_integral + ']')
    fout.write('\n')
    nthresh = len(energy_thresholds)

    for i in range(nthresh):
        if crossing_time[i] == 0: #no threshold crossed
            continue
        if is_diff_thresh[i]:
            fout.write(str(energy_thresholds[i]) + ',')
        else:
            fout.write('>' + str(energy_thresholds[i]) + ',')
        fout.write(str(flux_thresholds[i]) + ',')
        fout.write(str(crossing_time[i]) + ',')
        fout.write(str(onset_peak[i]) + ',')
        fout.write(str(onset_date[i]) + ',')
        fout.write(str(peak_flux[i]) + ',')
        fout.write(str(peak_time[i]) + ',')
        str_rise_time = str(rise_time[i]).split(',')
        for k in range(len(str_rise_time)):
            fout.write(str_rise_time[k] + ' ')
        fout.write(',')
        fout.write(str(event_end_time[i]) + ',')
        str_duration = str(duration[i]).split(',')
        for k in range(len(str_duration)):
            fout.write(str_duration[k] + ' ')
        fout.write(',' + str(threshold_fluences[i])) #units of flux*time
        if umasep:
            for jj in range(numa):
                fout.write(',' + str(umasep_times[i][jj] - crossing_time[i]) \
                    + ',' + str(umasep_fluxes[i][jj]))
        fout.write('\n')

    fout.close()
    return year, month, day, True


def error_check_options(experiment, flux_type, options, doBGSub):
    """ Make sure the selected options make sense for the experiment.
        
        INPUTS:
        
        :experiment: (string)
        :flux_type: (string) - integral or differential
        :options: (string array) - various options applied to GOES data
        :doBGSub: (boolean) - indicates if background subtraction to be
            performed
        
        OUTPUTS:
        
        no outputs by system exit if error found
    """
    if "S14" in options and experiment[0:4] != "GOES":
        sys.exit("Sandberg et al. (2014) effective energies (S14) may only "
            "be applied to GOES data.")
    if "S14" in options and "uncorrected" not in options:
        sys.exit("Sandberg et al. (2014) effective energies (S14) may only be "
            "applied to GOES uncorrected fluxes. Please add "
            "\"uncorrected\" to options.")
    if "uncorrected" in options and flux_type == "integral":
        sys.exit("The uncorrected option cannot be used with integral fluxes. "
                "Please remove this option and run again. Exiting.")
    if "uncorrected" in options and experiment[0:4] != "GOES":
        sys.exit("The uncorrected option may only be specified for GOES "
                "differential fluxes. Exiting.")
    if "Bruno2017" in options and (experiment != "GOES-13" and \
        experiment != "GOES-15"):
        sys.exit("Bruno2017 effective energies may only be appied to GOES-13 "
                "or GOES-15 fluxes. Exiting.")
    if "S14" in options and (experiment == "GOES-13" or \
        experiment == "GOES-15"):
        print("Sandberg et al. (2014) effective energies found for GOES-11 "
            "will be applied to channels P2-P7. Continuing.")
    if "S14" in options and "Bruno2017" in options:
        print("Sandberg et al. (2014) effective energies from GOES-11 will be "
            "applied to P2-P5. Bruno (2017) effective energies will be applied "
            "to P6-P11.")
    if doBGSub and experiment[0:4] == "GOES" and flux_type == "integral":
        sys.exit("Do not perform background subtraction on GOES integral "
                "fluxes. Integral fluxes have already been derived by "
                "applying corrections for cross-contamination and removing "
                "the instrument background levels.")
    if ("uncorrected" in options or "S14" in options or "Bruno2017" in options)\
                and experiment[0:4] != "GOES":
        sys.exit("The options you have selected are only applicable to GOES "
                "data. Please remove these options and run again: "
                "uncorrected, S14, or Bruno2017.")



def error_check_inputs(startdate, enddate, experiment, flux_type, is_diff_thresh):
    """ Check that all of the user inputs make sense and fall within bounds.
        
        INPUTS:
        
        :startdate: (datetime) - start of time period entered by user
        :enddate: (datetime) - end of time period entered by user
        :experiment: (string) - name of experiment specifed by user
        :flux_type: (string) - integral or differential
        :is_diff_thresh: (bool 1xn array) - where n indicates the number of
            thresholds input by the user, e.g. "30,1;50,1" n=2
            Indicates if the user-input thresholds apply to integral or
            differential channels
            
        OUTPUTS:
        
        None, but system exit if error found
    """
    #CHECKS ON INPUTS
    if (enddate < startdate):
        sys.exit('End time before start time! Enter a valid date range. '
                'Exiting.')
                
    if flux_type == "":
        sys.exit('User must indicate whether input flux is integral or '
                'differential. Exiting.')

    if (experiment == "SEPEM" and flux_type == "integral"):
        sys.exit('The SEPEM (RSDv2) data set only provides differential fluxes.'
            ' Please change your FluxType to differential. Exiting.')

    if (experiment == "SEPEMv3" and flux_type == "integral"):
        sys.exit('The SEPEM (RSDv3) data set only provides differential fluxes.'
            ' Please change your FluxType to differential. Exiting.')

    if ((experiment == "EPHIN" or experiment == "EPHIN_REleASE") \
        and flux_type == "integral"):
        sys.exit('The SOHO/EPHIN data set only provides differential fluxes.'
            ' Please change your FluxType to differential. Exiting.')

    for diff_thresh in is_diff_thresh:
        if diff_thresh and flux_type == "integral":
            sys.exit('The input flux type is specified as integral, but you '
                    'have requested a threshold in a differential energy bin. '
                    'Flux must be differential to impelement a threshold on a '
                    'differential energy bin. Exiting.')

    sepem_end_date = datetime.datetime(2015,12,31,23,55,00)
    if(experiment == "SEPEM" and (startdate > sepem_end_date or
                   enddate > sepem_end_date)):
        sys.exit('The SEPEM (RSDv2) data set only extends to '
                  + str(sepem_end_date) +
            '. Please change your requested dates. Exiting.')

    sepemv3_end_date = datetime.datetime(2017,12,31,23,55,00)
    if(experiment == "SEPEMv3" and (startdate > sepemv3_end_date or
                   enddate > sepemv3_end_date)):
        sys.exit('The SEPEM (RSDv3) data set only extends to '
                  + str(sepemv3_end_date) +
            '. Please change your requested dates. Exiting.')



def sort_bin_order(all_fluxes, energy_bins):
    """Check the order of the energy bins. Usually, bins go from
        low to high energies, but some modelers or users may
        go in reverse order. Usually expect:
        [[10,20],[20,30],[30,40]]
        But user may instead input files with fluxes in order of:
        [[30,40],[20,30],[10,20]]
        
        This subroutine will reorder the fluxes and energy_bins
        to go in increasing order. If differential fluxes were input,
        this reordering will ensure that integral fluxes are
        estimated properly.
        
        INPUTS:
        
        :all_fluxes: (float nxm array) - fluxes for n energy channels
            and m time points
        :energy_bins: (float 2xn array) - energy bins for each of the
            energy channels
            
        OUTPUTS:
        
        :sort_fluxes: (float nxm array) - same as above, but sorted so
            that the lowest energy channel is first and highest is last
        :sort_bins: (float 2xn array) - same as above, but sorted
        
    """

    nbins = len(energy_bins)
    #Rank energy bins in order of lowest to highest effective
    #energies
    eff_en = []
    for i in range(nbins):
        if energy_bins[i][1] == -1:
            eff_en.append(energy_bins[i][0])
        else:
            midpt = math.sqrt(energy_bins[i][0]*energy_bins[i][1])
            eff_en.append(midpt)
            
    eff_en_np = np.array(eff_en)
    sort_index = np.argsort(eff_en_np) #indices in sorted order
    
    sort_fluxes = np.array(all_fluxes)
    sort_bins = []
    for i in range(nbins):
        sort_fluxes[i] = all_fluxes[sort_index[i]]
        sort_bins.append(energy_bins[sort_index[i]])
    
    return sort_fluxes, sort_bins


####TOOLS#####
def get_input_thresholds(str_thresh):
    """ Get user specified input thresholds,
        
        INPUTS:
        :str_thresh: (string) input by user with thresholds
            separated by semi-colons
        
        OUTPUTS:
        
        :input_threshold: (2xn floats) array as
            [[binLow, flux threshold],
             [bin Low, fluxthreshold],
             ...]
        
        :is_diff_thresh: (bool array)
    """
    input_threshold = []
    is_diff_thresh = []
    if str_thresh == "":
        return input_threshold, is_diff_thresh
    
    nin_thresh = len(str_thresh)
    is_diff_thresh = [False]*nin_thresh #True if differential flux threshold
    for i in range(nin_thresh):
        if "-" in str_thresh[i][0]:
            is_diff_thresh[i] = True
            thresh0 = str_thresh[i][0].split("-")
            input_threshold.append([float(thresh0[0]), float(str_thresh[i][1])])
            print("Found differential threshold " + str_thresh[i][0])
        else:
            input_threshold.append([float(str_thresh[i][0]), \
                                    float(str_thresh[i][1])])
                                    
    return input_threshold, is_diff_thresh


def str_to_datetime(date):
    """ String date to datetime
        
        INPUTS:
        
        :date: (string) - date as "YYYY-MM-DD" or "YYYY-MM-DD HH:MM:SS"
        
        OUTPUTS:
        
        :dt: (datetime) - datetime conversion of date
    """
    if len(date) == 10: #only YYYY-MM-DD
        date = date  + ' 00:00:00'
    dt = datetime.datetime.strptime(date, "%Y-%m-%d %H:%M:%S")
    return dt


def read_in_flux_files(experiment, flux_type, user_file, model_name, startdate,
        enddate, str_startdate, str_enddate, str_bgstartdate, str_bgenddate,
        options, doBGSub, nointerp, showplot, saveplot):
    """ Read in the appropriate data or user files. Performs
        background subtraction, if requested. Trims to dates
        between start time and end time. Interpolates bad
        points with linear interpolation in time.
        
        INPUTS:
        
        :experiment: (string)
        :flux_type: (string) - integral, differential
        :user_file: (string) - file containing user's flux time profiles
        :model_name: (string) - model name or experiment if experiment = "user"
        :startdate: (datetime) - start date of time period entered by user
        :enddate: (datetime) - end date of time period entered by user
        :str_startdate: (string) - start date in string form
        :str_enddate: (string) - end date in string form
        :str_bgstartdate: (string) - start date of time period to be used for
            background subtraction
        :str_bgenddate: (string) - end date of time period to be used or
            background subtraction
        :options: (string array) - options that could be applied
        :doBGSub: (bool) - indicate if background subtraction to be performed
        :nointerp: (bool) - indicates if user doesn't want to do linear
            interpolation in time for negative of bad flux values
        :showplot: (bool)
        :saveplot: (bool) - save plots automatically to "plots" directory
        
        OUTPUTS:
        
        :dates: (datetime 1xm array) - times in flux time profile trimmed
            between startdate and enddate
        :fluxes: (numpy float nxm array) - fluxes for n energy channels and m
            time steps; these are background subtracted fluxes if background
            subtraction was selected.
        :energy_bins: (array nx2 for n thresholds)
        
    """
    
    filenames1, filenames2, filenames_orien = datasets.check_data(startdate,
                                    enddate, experiment, flux_type, user_file)
                                    
    #read in flux files
    if experiment != "user":
        all_dates, all_fluxes, west_detector =datasets.read_in_files(experiment,
                    flux_type, filenames1, filenames2, filenames_orien, options)
    if experiment == "user":
        all_dates, all_fluxes = datasets.read_in_user_files(filenames1)
        west_detector = []
    #Define energy bins
    energy_bins = datasets.define_energy_bins(experiment, flux_type, \
                                west_detector, options)
    
    all_fluxes, energy_bins = sort_bin_order(all_fluxes, energy_bins)
    
    #IF REQUESTED BACKGROUND SUBTRACTION
    if doBGSub:
        #sepfluxes are background subtracted fluxes
        bgfluxes, sepfluxes, bgdates = bgsub.derive_background(str_startdate, \
                    str_enddate, str_bgstartdate, str_bgenddate, experiment, \
                    flux_type, model_name,user_file, showplot, saveplot,options)
        #Extract the date range specified by the user
        dates, fluxes = datasets.extract_date_range(startdate, enddate,
                                    bgdates, sepfluxes)

    #NO BACKGROUND SUBTRACTION
    if not doBGSub:
        #Extract the date range specified by the user
        dates, fluxes = datasets.extract_date_range(startdate, enddate,
                                all_dates, all_fluxes)
    
    if nointerp:
        #Set bad data points (negative flux or None) to None
        fluxes = datasets.check_for_bad_data(dates,fluxes,energy_bins,nointerp)
    else:
        #Remove bad data points (negative flux or None) w/ linear interp in time
        fluxes = datasets.check_for_bad_data(dates,fluxes,energy_bins)
        
    if len(dates) <= 1:
        sys.exit("The specified start and end dates were not present in the "
                "specified input file. Exiting.")
    
    return dates, fluxes, energy_bins
      
      
      
def define_thresholds(input_threshold, is_diff_thresh, str_thresh, energy_bins,
        umasep):
    """ Two operational thresholds plus any thresholds
        specified by the user.
        
        INPUTS:
        
        :input_threshold: (float 2xn array) - thresholds specified by user e.g.
            [[30,1],[50,1]]
        :is_diff_thresh: (bool 1xn array) - indicates whether the input threshold
            is a differential channel
        :str_thresh: (string) - threshold as input by user in command line e.g.
            "30,1;50,1;4.9-8.3,10"
            
        OUTPUTS:
        
        :energy_thresholds: (float 1x(n+2) array) - e.g. [10,100,30,50,4.9]
            energy channels for which a threshold will be applied. Includes
            the two operational thresholds plus any thresholds entered by
            the user
        :flux_thresholds: (float 1x(n+2) array) - e.g. [10,1,1,1,10]
            flux thresholds to be applied to the energy channels in
            energy_thresholds. Includes the two operational thresholds
            plus any thresholds entered by the user
        
    """
    energy_thresholds = [10,100] #MeV; flux for particles of > this MeV
    flux_thresholds = [10,1] #pfu; exceed this level of intensity
    
    if input_threshold != []:
        nin_thresh = len(input_threshold)
        for i in range(nin_thresh):
            if not is_diff_thresh[i]:
                energy_thresholds.append(input_threshold[i][0])
                flux_thresholds.append(input_threshold[i][1])
            #Check if user entered differential threshold
            #Don't append to thresholds; deal with these later
            if is_diff_thresh[i]:
                check_bin_exists(str_thresh[i],energy_bins)

    if umasep: #add two additional thresholds
        energy_thresholds.append(30)
        flux_thresholds.append(1) #Used for UMASEP-30 development
        energy_thresholds.append(50)
        flux_thresholds.append(1) #Used for UMASEP-50 development
                 
    return energy_thresholds, flux_thresholds
                                    

def calculate_integral_fluences(experiment, flux_type, options,
        model_name, doBGSub, startdate, enddate, energy_bins,
        energy_thresholds, flux_thresholds, dates, fluxes,
        integral_fluxes,
        crossing_time, event_end_time, all_threshold_fluences,
        all_fluence, all_energies):
    """ Calculate fluence values for all integral energy thresholds
        between start and end times determined in each channel.
        This subroutine is only called for the thresholds applied
        to integral channels.
        
        INPUTS:
        
        :experiment: (string)
        :flux_type: (string) - integral, differential
        :options: (string array) - options that could be applied
        :model_name: (string) - model name or experiment if experiment = "user"
        :doBGSub: (bool) - indicate if background subtraction to be performed
        :startdate: (datetime) - start date of time period entered by user
        :enddate: (datetime) - end date of time period entered by user
        :energy_bins: (float 2xp array) - energy bins for all of the flux channels
        :energy_thresholds: (float 1xn array) - energy channels for which a
            threshold is applied (at this point, only the integral ones)
        :flux_thresholds: (float 1xn array) - flux thresholds applied to
            the energy channels in energy_thresholds
        :dates: (datetime 1xm array) - time points for the flux time profiles
        :fluxes: (float pxm array) - flux time profiles for p energy bins and m
            time steps
        :integral_fluxes: (float nxm array) - flux time profiles of integral
            or estimated integral fluxes for n energy channels and m time steps
        :crossing_time: (datetime 1xn array) - SEP event start times for each
            energy channel with an applied threshold
        :event_end_time: (datetime 1xn array) - SEP event end times for each
            energy channel with an applied threshold
        :all_threshold_fluences: (float 1xn array) - fluence values for integral
            or estimated integral fluxes in the n energy channels for which
            thresholds were applied
        :all_fluence: (float nxp array) - fluence spectrum in the original
            integral or differential energy channels associated with
            energy_bins. Fluence calculated between crossing_time and
            event_end_time, so this is the fluence spectrum associated with
            the events defined by the thresholds applied to each energy channel
        :all_energies: (float 1xp array) - effective energies for energy_bins,
            defined by sqrt(low bin edge*high bin edge)
        
        
        OUTPUTS:
        
        Arrays filled in for the integral thresholds:
        
        :all_threshold_fluences: (float 1xn array) - fluence values for integral
            or estimated integral fluxes in the n energy channels for which
            thresholds were applied
        :all_fluence: (float nxp array) - fluence spectrum in the original
            integral or differential energy channels associated with
            energy_bins. Fluence calculated between crossing_time and
            event_end_time, so this is the fluence spectrum associated with
            the events defined by the thresholds applied to each energy channel
        :all_energies: (float 1xp array) - effective energies for energy_bins,
            defined by sqrt(low bin edge*high bin edge)
        
    """
    nthresh = len(energy_thresholds)
    for i in range(nthresh):
        #If no threshold was crossed during specified date range
        if crossing_time[i] == 0:
            print("The >" + str(energy_thresholds[i]) + " "
                    + energy_units + " " + energy_units + " threshold was "
                     "not crossed during the specified date range. No SEP "
                     "event. Continuing.")
            continue

        #Extract the original fluxes trimmed for the SEP start and stop times
        sep_dates, sep_fluxes = datasets.extract_date_range(crossing_time[i],
                             event_end_time[i],dates,fluxes)

        #Calculate fluence spectrum for the SEP event.
        #Fluence spectrum will be of integral fluxes if the original
        #data set read in was integral; fluence spectrum will be of
        #differential fluxes if original data set was differential
        print('=====Calculating event fluence for event defined by >'
                + str(energy_thresholds[i]) + ' ' + energy_units + ', for '
                + str(crossing_time[i]) + ' to ' + str(event_end_time[i]))
        if crossing_time[i] == event_end_time[i]:
            sys.exit("Event start and end time the same (did you set "
            "--DetectPreviousEvent? May not work in this case). Exiting.")

        fluence, energies = get_fluence_spectrum(experiment, flux_type,
                         options, doBGSub,
                         model_name, energy_thresholds[i], flux_thresholds[i],
                         sep_dates, sep_fluxes, energy_bins, False, True)
                         #diff_thresh; savefile
                         #Only thresholds applied to integral flux
                         #channels are specified so far
                         #The False indicates the threshold is associated
                         #with integral channels and has nothing to do
                         #with whether the fluxes themselves are integral
                         #or differential

        all_fluence[i] = fluence #in native units of experiment
        all_energies[i] = energies

        #Always calculate fluences for integral fluxes >10, >100 MeV and
        #any user input thresholds applied to integral channels
        #The fluences produced here are only for integral flux channels
        #that have had a threshold applied and the fluence is only
        #for that channel (not a spectrum as calculated above)
        if flux_type == "integral":
            integral_fluence = report_threshold_fluences(experiment, flux_type,
                        model_name, energy_thresholds, energy_bins,
                        sep_dates, sep_fluxes)
            
        if flux_type == "differential":
            #Extract the estimated integral fluxes in the SEP event date range
            sep_integral_dates, sep_integral_fluxes = \
                            datasets.extract_date_range(\
                            crossing_time[i],event_end_time[i],
                            dates, integral_fluxes)

            integral_fluence = report_threshold_fluences(experiment, flux_type,
                         model_name, energy_thresholds, energy_bins,
                         sep_integral_dates, sep_integral_fluxes)

        #integral_fluence produced by report_threshold_fluences has fluences in the
        #same index order as the energies specified in energy_thresholds
        all_threshold_fluences[i] = integral_fluence[i]

    return all_threshold_fluences, all_fluence, all_energies



def append_differential_thresholds(energy_thresholds, flux_thresholds,
        options, str_thresh, is_diff_thresh, plt_energy, plt_flux,
        input_threshold, energy_bins,
        dates, fluxes, detect_prev_event, two_peaks,
        umasep, umasep_times, umasep_fluxes,
        all_threshold_fluences, all_fluence, all_energies,
        crossing_time, peak_flux, peak_time, rise_time, event_end_time,
        duration, onset_date, onset_peak, integral_fluxes):
    """ Add the threshold crossing information for differential channels
        as specified by the user.
        
        Definitions:
        
        * n integral energy thresholds (>10, >100 + user input integral thresholds)
        * m energy bins associated with the input flux files
        * p user input thresholds (integral + differential, define d differential channels for which a threshold is applied n-2 = integral channels, (n-2)+d = p)
        * q time points in the flux time profiles (between user input start and end dates)
        
        INPUTS (filled in up to thresholds applied to integral channels):
        
        :energy_thresholds: (float 1xn array) - energy channels for which a
            threshold is applied (at this point, only the integral ones)
        :flux_thresholds: (float 1xn array) - flux thresholds applied to
            the energy channels in energy_thresholds
        :str_thresh: (string 1xp array) - p user input thresholds broken into
            an array (need for generating plot labels)
        :is_diff_thresh: (bool 1xp array) - indicates which of the user
            thresholds are applied to differential channels
        :plt_energy: (float 1xn array) - holds energy thresholds for the
            integral channels at this point
        :plt_flux: (float 1xn array) - hold the flux thresholds for the
            integral channels at this point
        :input_threshold: (float 2xp array) - p input thresholds in form
            [binlowedge, flux threshold]
        :energy_bins: (float 2xm array) - m energy bins for each input
            flux channel
        :dates: (datetime 1xq array) - q time points for the flux time profile
        :fluxes: (float mxq array) - flux time profiles for m energy bins and
            q time points
        :detect_prev_event: (bool) - option for finding start of event
        :two_peaks: (bool) - option for extending event length
        :umasep: (boolean) - indicate if user called UMASEP flag
        :umasep_times: (datetime 5xn array) - times 3, 4, 5, 7, 9 hours after
            crossing time for each applied threshold
        :umasep_fluxes: (float 5xn) - fluxes at each of those times for each
            applied threshold
        :all_threshold_fluences: (float 1xn array) - fluence values for integral
            or estimated integral fluxes in the n energy channels for which
            thresholds were applied
        :all_fluence: (float nxm array) - fluence spectrum in the original
            integral or differential energy channels associated with
            energy_bins. Fluence calculated between crossing_time and
            event_end_time, so this is the fluence spectrum associated with
            the events defined by the thresholds applied to each energy channel
        :all_energies: (float 1xm array) - effective energies for energy_bins,
            defined by sqrt(low bin edge*high bin edge)
        :crossing_time: (datetime 1xn array) - SEP event start time for each
            applied threshold
        :peak_flux: (float 1xn array) - maximum flux for each applied threshold
        :peak_time: (datetime 1xn array) - time of maximum flux for each applied
            threshold
        :rise_time: (timedelta 1xn array) - start time to max flux
        :event_end_time: (datetime 1xn array) - end time for each applied threshold
        :duration: (timedelta 1xn array) - end time minus start time
        :onset_peak: (float 1xn array) - onset peak for each applied threshold
        :onset_date: (datetime 1xn array) - onset peak time for each applied
            threshold
        :integral_fluxes: (float nxq array) - flux time profiles of integral
            or estimated integral fluxes for n energy channels and m time steps
       
        
        OUTPUTS:
        
        (thresholds applied to differential channels have been
        appended; define d as the number of differential channels for
        which the user applied a threshold.
        n+d = total number of thresholds applied to energy channels):
        
        :energy_thresholds: (float 1x(n+d) array)
        :flux_thresholds: (float 1x(n+d) array)
        :plot_diff_thresh: (bool 1x(n+d) array)
        :plt_energy: (float 1x(n+d) array)
        :plt_flux: (float 1x(n+d) array)
        :all_threshold_fluences: (float 1x(n+d) array)
        :all_fluence: (float (n+d)xm array)
        :all_energies: (float 1xm array)
        :crossing_time: (datetime 1x(n+d) array)
        :peak_flux: (float 1x(n+d) array)
        :peak_time: (datetime 1x(n+d) array)
        :rise_time: (timedelta 1x(n+d) array)
        :event_end_time: (datetime 1x(n+d) array)
        :duration: (timedelta 1x(n+d) array)
        :onset_date: (datetime 1x(n+d) array)
        :onset_peak: (float 1x(n+d) array)
        :integral_fluxes: (float (n+d)xq array)
        
    """
    #plot_diff_threshold will indicate whether all thresholds applied
    #are integral or differential for all channels
    #So far, is_diff_thresh only held entries related to the
    #thresholds input by the user. plot_diff_thresh includes the native
    #>10 MeV, 10 pfu and >100 MeV, 1 pfu thresholds
    #Because all the integral channel thresholds were already grouped
    #together in energy_thresholds, fill plot_diff_thresh for all
    #integral channel thresholds, then tack on the differential channels
    #at the end
    nthresh = len(energy_thresholds) #only integral channel thresholds so far
    plot_diff_thresh = [False]*nthresh  #integral thresholds

    #If a differential threshold was specified, calculate everything with
    #correct differential channel
    nin_thresh = len(input_threshold)
    for i in range(nin_thresh):
        if is_diff_thresh[i]:
            plot_diff_thresh.append(True) #diff tacked onto end
            plt_energy.append(str_thresh[i][0])
            plt_flux.append(str_thresh[i][1])
            print("=====INFORMATION FOR DIFFERENTIAL BIN=======")
            #identify which bin is the correct energy bin
            svbin = 0
            for k in range(len(energy_bins)):
                if energy_bins[k][0] == input_threshold[i][0]:
                    svbin = k
            #Expect an nxm array, but looking at one flux profile at
            #a time. Hack it by putting the single threshold info
            #into arrays
            energy_thresh = [input_threshold[i][0]] #calculate_onset_peak needs arrays
            flux_thresh = [input_threshold[i][1]]
            in_flx = np.array([fluxes[svbin]])
            ct,pf,pt,rt,eet,dur=calculate_event_info(energy_thresh,\
                        flux_thresh, dates, in_flx, detect_prev_event,two_peaks,
                        is_diff_thresh[i])
            if ct[0] == 0:
                print("The energy bin " + str_thresh[i][0] + " "
                        + energy_units +
                        " threshold was not crossed during the specified date "
                        "range. No SEP event. Continuing.")
                fl = np.array([0]*len(energy_bins))
                en = np.array([0]*len(energy_bins))
                bin_fl = 0
                pf = [np.amax(fluxes[svbin])]
                pf_idx = np.where(fluxes[svbin] == np.amax(fluxes[svbin]))
                pt = [dates[pf_idx[0][0]]]
                op = pf #onset peak = max flux when no threshold crossed
                od = pt #time of max flux
                
            else:
                #Extract the original fluxes for the SEP start and stop times
                sep_d, sep_f = datasets.extract_date_range(ct[0],eet[0],dates,
                                fluxes)
                fl, en = get_fluence_spectrum(experiment, flux_type,
                                 options, doBGSub,
                                 model_name, input_threshold[i][0],
                                 input_threshold[i][1], sep_d, sep_f,
                                 energy_bins, is_diff_thresh[i], True) #savefile
                bin_fl = fl[svbin] #fluence for bin associated with threshold
                
                od,op=calculate_onset_peak(experiment, energy_thresh,
                            dates, in_flx, ct, eet, showplot)

                if umasep:
                    umasep_t, umasep_f = calculate_umasep_info(\
                                [input_threshold[i][0]],[input_threshold[i][1]],
                                dates, [fluxes[svbin]], ct)
                    umasep_times.append(umasep_t[0])
                    umasep_fluxes.append(umasep_f[0])

            #user-specified differential thresholds will be tacked onto the end
            all_fluence = np.append(all_fluence, [fl], axis=0) #in native units of experiment
            all_energies = np.append(all_energies, [en], axis=0)
            crossing_time.append(ct[0])
            peak_flux.append(pf[0])
            peak_time.append(pt[0])
            rise_time.append(rt[0])
            event_end_time.append(eet[0])
            duration.append(dur[0])
            onset_date.append(od[0])
            onset_peak.append(op[0])
            energy_thresholds.append(input_threshold[i][0])
            flux_thresholds.append(input_threshold[i][1])
            all_threshold_fluences.append(bin_fl) #differential bin fluence for threshold
            integral_fluxes = np.append(integral_fluxes, [fluxes[svbin]], \
                                    axis=0)

    return energy_thresholds, flux_thresholds, plot_diff_thresh, plt_energy,\
        plt_flux, all_threshold_fluences, all_fluence, all_energies,\
        crossing_time, peak_flux, peak_time, rise_time, event_end_time,\
        duration, onset_date, onset_peak, integral_fluxes

 
def write_zulu_time_profile(filename, dates, fluxes):
    """ Write out the time profile with the date in the
        first column as the ISO standard and flux in the
        second column as:
        
        YYYY-MM-DDTHH:MM:SSZ    Float
        
        INPUTS:
        
        :Filename: (string) - name of file to write
        :date: (datetime 1xn array) - list of dates
        :fluxes: (float 1xn array) - corresponding fluxes
        
        OUTPUTS:
        
        None but writes output file with filename
    """
    fname = outpath + "/" + filename
    outfile = open(fname, "w")
    for i in range(len(dates)):
        zdate = ccmc_json.make_ccmc_zulu_time(dates[i])
        outfile.write(zdate + "    " + str(fluxes[i]) + "\n")
        
    outfile.close()
    print("write_zulu_time_profile: Wrote file --> " + fname)



def write_info_to_file(experiment, flux_type, options, doBGSub,
        energy_bins, model_name, spase_id, startdate, enddate,
        energy_thresholds, flux_thresholds, dates, integral_fluxes,
        crossing_time, onset_peak, onset_date, peak_flux, peak_time,
        rise_time, event_end_time, duration, all_threshold_fluences,
        all_fluence, plot_diff_thresh,
        umasep, umasep_times, umasep_fluxes):
    """ Write information to csv and json file.
        Writes all of the derived information for all the
        energy-threshold combinations to file.
        
        INPUTS:
        
        The inputs here are the same as the outputs produced by
        append_differential_thresholds(). Please see the list
        of those outputs for detailed explanations of these inputs.
        
        OUTPUTS:
        
        No outputs except files are written:
            
            * fluence spectra files for each energy channel for which a threshold was applied
            * integral fluxes in csv file for which a threshold was applied
            * flux time profiles in txt files for each energy channel for which a threshold was applied
            * json file in CCMC format containing blocks of information for each energy channel for which a threshold was applied
            * csv file containing all timing and peak flux values
            
    """
    
    #Save all calculated values for all threshold definitions to csv file
    sep_year, sep_month, sep_day, IsCrossed = print_values_to_file(experiment,
                    flux_type, options, doBGSub,
                    model_name, startdate, energy_thresholds, flux_thresholds,
                    crossing_time, onset_peak, onset_date, peak_flux, peak_time,
                    rise_time, event_end_time, duration, all_threshold_fluences,
                    plot_diff_thresh, umasep, umasep_times, umasep_fluxes)
    
    
    #SAVE TO JSON FILE
    type = "observations"
    if experiment == "user":
        type = "model"
    template = ccmc_json.read_in_json_template(type)

    modifier = ''
    if options[0] != '':
        options = sorted(options) #to always make consistent filenames
        for opt in options:
            modifier = modifier + '_' + opt
    if doBGSub:
        modifier = modifier + '_bgsub'
        options.append("BGSubtracted")
        
        
    #Get issue time of forecast (now)
    now = datetime.datetime.now()
    issue_time = ccmc_json.make_ccmc_zulu_time(now)
    
    #Generate filenames in CCMC preferred format
    #issue time in CCMC format YYYY-MM-DDTHH:MM:SSZ
    zstdate = ccmc_json.make_ccmc_zulu_time(startdate)
    #For native data set which is an observation, don't include issue time
    #in the filename since nothing should change in the measurements
    if experiment != "user":
        fnameprefix = experiment + "_" + flux_type + modifier + "." + zstdate.replace(":","")
    
    #For a user data set or model, include issue time
    if experiment == "user":
        fnameprefix = experiment + "_" + flux_type + modifier + "." + zstdate.replace(":","") + "." + issue_time.replace(":","")
    if experiment == 'user' and model_name != '':
            fnameprefix = model_name + "_" + flux_type + modifier + "." + zstdate.replace(":","") + "." + issue_time.replace(":","")
            
    jsonfname = outpath +'/' + fnameprefix + ".json"
    
    #filenames for time profiles
    proffnames = []
    for j in range(len(energy_thresholds)):
        energy = energy_thresholds[j]
        if (not plot_diff_thresh[j] and flux_type == "differential")\
            or flux_type == "integral":
            profname = fnameprefix + "." + str(energy) \
                        + energy_units + ".txt"
        else:
            bin = ccmc_json.find_energy_bin(energy, energy_bins)
            profname = fnameprefix + "." + str(bin[0]) + "-" + str(bin[1])\
                        + energy_units + ".txt"
        proffnames.append(profname)
    
    ##### WRITE JSON FILE #######
    filled_json = ccmc_json.fill_json(template, issue_time,
                    experiment, flux_type,
                    energy_bins, model_name, spase_id, startdate, enddate,
                    options, energy_thresholds, flux_thresholds, crossing_time,
                    onset_peak, onset_date, peak_flux, peak_time, rise_time,
                    event_end_time, duration, all_threshold_fluences,
                    plot_diff_thresh, all_fluence,
                    umasep, umasep_times, umasep_fluxes, proffnames,
                    energy_units, flux_units_integral, fluence_units_integral,
                    flux_units_differential, fluence_units_differential)
    
    filled_json = ccmc_json.clean_json(filled_json,experiment)
    isgood = ccmc_json.write_json(filled_json, jsonfname)
    if not isgood:
        print("WARNING: ccmc_json_handler: write_json could not write your " \
                "file "+ str(jsonfname))
    
    ##### WRITE INDIVIDUAL FLUX FILES #####
    #Note that integral_fluxes actually contains time profiles for all
    #the energy channels for which a threshold was applied, including
    #for any differential channels
    for j in range(len(energy_thresholds)):
        profnm = proffnames[j]
        write_zulu_time_profile(profnm, dates, integral_fluxes[j])
        
    #IF NO THRESHOLDS CROSSED, EXIT PROGRAM. ONLY PLOTTING REMAINS.
    if not IsCrossed:
        sys.exit("No thresholds were crossed during this time period. "
                "Max flux has been written to json file in onset peak in file "
                + jsonfname + ". Exiting. ")


    return sep_year, sep_month, sep_day, jsonfname


######## MAIN PROGRAM #########
def run_all(str_startdate, str_enddate, experiment, flux_type, model_name,
        user_file, spase_id, showplot, saveplot, detect_prev_event,
        two_peaks, umasep, str_thresh, options, doBGSub, str_bgstartdate,
        str_bgenddate, nointerp=False):
    """"Runs all subroutines and gets all needed values. Takes the command line
        arguments as input. Code may be imported into other python scripts and
        run using this routine.
        
        INPUTS:
        
        :str_startdate: (string) - user input start date "YYYY-MM-DD" or
            "YYYY-MM-DD HH:MM:SS"
        :str_enddate: (string) - user input end date "YYYY-MM-DD" or
            "YYYY-MM-DD HH:MM:SS"
        :experiment: (string) - "GOES-08" up to "GOES-15", "SEPEM", "SEPEMv3",
            "EPHIN", "EPHIN_REleASE", or "user"
        :flux_type: (string) - "integral" or "differential" indicates the type
            of flux to read in
        :model_name: (string) - If model is "user", set model_name to describe
            your model or data set (e.g. MyModel), otherwise set to ''.
        :user_file: (string) - Default is ''. If "user" is selected for experiment,
            specify name of flux file.
        :spase_id: (string) - Default is ''. If you know the spase_id of you
            model or experiment, enter it here for the json file
        :showplot: (bool) - Set to True to show plots when run
        :saveplot: (bool) - Set to True to automatically save plots to the
            plots directory when run
        :detect_prev_event: (bool) - option for finding start of event
        :two_peaks: (bool) - option for extending event length
        :umasep: (boolean) - call flag to run code for specific time related
            to the UMASEP model
        :str_thresh: (string) - user-input thresholds in the format "30,1"
            for >30 MeV exceeds 1 pfu, "4-7,0.01" for 4-7 MeV differential
            channel exceeds 0.01.  "30,1;4-7,0.01" multiple thresholds
            separated by semi-colon.
        :nointerp: (boolean) - set to true to fill in negative fluxes with None
            value rather than filling in via linear interpolation in time
        
        OUTPUTS:
        
        :sep_year: (integer) - year of start of SEP event or year of beginning
            of user-input time period (if no threshold crossed)
        :sep_month: (integer) - month of start of SEP event or year of beginning
            of user-input time period (if no threshold crossed)
        :sep_day: (integer) - day of start of SEP event or year of beginning
            of user-input time period (if no threshold crossed)
        :jsonfname: (string) - name of saved json file (could include issue time,
            which is the time that the program is run, so absolutely
            need to pass file name back as argument)
        :and generates multiple plots:
        
    """
    #Check for empty dates
    if (str_startdate == "" or str_enddate == ""):
        sys.exit('You must enter a valid date range. Exiting.')

    #PROCESS INPUTS
    options = options.split(";")
    user_fname[0] = user_file #input as argument, default is 'tmp.txt'
    
    #Break up thresholds into a string array
    input_threshold = []
    is_diff_thresh = []
    if str_thresh != "":
        str_thresh = str_thresh.strip().split(";")
        for kk in range(len(str_thresh)):
            str_thresh[kk] = str_thresh[kk].strip().split(",")
        input_threshold, is_diff_thresh = get_input_thresholds(str_thresh)
    
    startdate = str_to_datetime(str_startdate)
    enddate = str_to_datetime(str_enddate)

    #PERFORM CHECKS AND VALIDATE INPUTS
    error_check_options(experiment, flux_type, options, doBGSub)
    error_check_inputs(startdate, enddate, experiment, flux_type, is_diff_thresh)
    datasets.check_paths()
    
    #SET UNITS
    if experiment == "user":
        global energy_units
        global flux_units_integral
        global fluence_units_integral
        global flux_units_differential
        global fluence_units_differential
        
        energy_units = vars.energy_units
        flux_units_integral = vars.flux_units_integral
        fluence_units_integral = vars.fluence_units_integral
        flux_units_differential = vars.flux_units_differential
        fluence_units_differential = vars.fluence_units_differential
    
    #READ IN FLUXES AND ENERGY BINS, BG SUBTRACT, INTERPOLATE
    dates, fluxes, energy_bins = read_in_flux_files(experiment,
        flux_type, user_file, model_name, startdate,
        enddate, str_startdate, str_enddate, str_bgstartdate,
        str_bgenddate, options, doBGSub, nointerp, showplot, saveplot)
   
    #Define ONLY INTEGRAL thresholds to use for start and end of event
    #Will tack on differential thresholds later on
    #At this point, energy_thresholds and flux_thresholds correspond
    #ONLY to threshold applied to integral flux channels
    energy_thresholds, flux_thresholds = define_thresholds(input_threshold,
                is_diff_thresh, str_thresh, energy_bins, umasep)
    
    #Estimate or select integral fluxes corresponding the energy_thresholds
    #Pull out or estimate only the integral flux channels for which a
    #threshold will be applied
    integral_fluxes = extract_integral_fluxes(fluxes, experiment, flux_type,
                    flux_thresholds, energy_thresholds, energy_bins, options)

    #Calculate SEP event quantities for energy and flux threshold combinations
    #integral fluxes are used to define event start and stop
    crossing_time,peak_flux,peak_time,rise_time,event_end_time,duration=\
        calculate_event_info(energy_thresholds,flux_thresholds, dates,
                    integral_fluxes, detect_prev_event, two_peaks, False)
                    #Always integral fluxes here, so False

    #Calculate onset peak for all thresholds
    onset_date, onset_peak = calculate_onset_peak(experiment, energy_thresholds,
                dates, integral_fluxes, crossing_time, event_end_time, showplot)

    #Calculate times used in UMASEP
    umasep_times =[]
    umasep_fluxes=[]
    if umasep:
        umasep_times, umasep_fluxes = calculate_umasep_info(energy_thresholds,
                        flux_thresholds, dates, integral_fluxes, crossing_time)

    #Arrays for event-integrated fluences for all thresholds, both
    #integral and differential
    nthresh = len(energy_thresholds)
    all_threshold_fluences = [0]*nthresh #fluences corresponding to >10, >100 MeV and other thresholded bins
    all_fluence = np.zeros(shape=(nthresh,len(energy_bins))) #fluence spectrum
    all_energies = np.zeros(shape=(nthresh,len(energy_bins))) #corresponding energy bin centers

    #------INTEGRAL THRESHOLD FLUENCES------
    #Fill in fluence arrays for the integral channels
    all_threshold_fluences, all_fluence, all_energies = \
        calculate_integral_fluences(experiment, flux_type, options,
            model_name, doBGSub, startdate, enddate, energy_bins,
            energy_thresholds, flux_thresholds, dates, fluxes,
            integral_fluxes, crossing_time, event_end_time,
            all_threshold_fluences,
            all_fluence, all_energies)
    
    #Save to file the integral fluxes for all integral channels
    #where thresholds were applied - multiple fluxes in a single file
    save_integral_fluxes_to_file(experiment, flux_type, options, doBGSub,
                model_name, energy_thresholds, crossing_time, dates,
                integral_fluxes)

    #------APPEND INFO FOR DIFFERENTIAL THRESHOLDS------
    #I'm sorry. This is a mess. I will rewrite eventually
    plot_diff_thresh = [False]*nthresh
    plt_energy = []  #string thresholds for plot labels
    plt_flux = []
    for i in range(nthresh):
        plt_energy.append(str(energy_thresholds[i]))
        plt_flux.append(str(flux_thresholds[i]))
    
    if True in is_diff_thresh:
        energy_thresholds, flux_thresholds, plot_diff_thresh,\
        plt_energy, plt_flux,\
        all_threshold_fluences, all_fluence, all_energies,\
        crossing_time, peak_flux, peak_time, rise_time, event_end_time,\
        duration, onset_date, onset_peak, integral_fluxes\
            = append_differential_thresholds(energy_thresholds,
                flux_thresholds, options,
                str_thresh, is_diff_thresh, plt_energy, plt_flux,
                input_threshold, energy_bins,
                dates, fluxes, detect_prev_event, two_peaks,
                umasep, umasep_times, umasep_fluxes,
                all_threshold_fluences, all_fluence, all_energies,
                crossing_time, peak_flux, peak_time, rise_time, event_end_time,
                duration, onset_date, onset_peak, integral_fluxes)

    #If threshold isn't crossed, the maximum flux was stored in onset
    #peak and peak_flux was saved as 0.
    #Transfer onset_peak to peak_flux and set onset_peak to zero.
    onset_peak, onset_date, peak_flux, peak_time = consistent_peak_fluxes(\
                crossing_time, onset_peak, onset_date, peak_flux, peak_time)
    #####################################################################
    #Write information to csv and json files
    #Note that, at this point in the code, integral_fluxes contains
    #both the integral flux time profiles for which thresholds were applied
    #and any differential flux channel time profiles for which
    #thresholds were applied
    sep_year, sep_month, sep_day, jsonfname\
        = write_info_to_file(experiment, flux_type, options, doBGSub,
        energy_bins, model_name, spase_id,  startdate, enddate,
        energy_thresholds, flux_thresholds, dates, integral_fluxes,
        crossing_time, onset_peak, onset_date, peak_flux, peak_time,
        rise_time, event_end_time, duration, all_threshold_fluences,
        all_fluence, plot_diff_thresh,
        umasep, umasep_times, umasep_fluxes)


    
    #####################################################################
    #===============PLOTS==================
    if saveplot or showplot:
        #Plot selected results
        #Event definition from integral fluxes
        if flux_type == "differential":
            print("Generating figure of estimated integral fluxes with threshold "
                   "crossings.")
        if flux_type == "integral":
            print("Generating figure of integral fluxes with threshold crossings.")

        #Additions to titles and filenames according to user-selected options
        modifier = ''
        title_mod = ''
        if "uncorrected" in options:
            modifier = modifier + '_uncorrected'
            title_mod = title_mod + 'uncorrected '
        if doBGSub:
            modifier = modifier + '_bgsub'
            title_mod = title_mod + 'BG-subtracted '
        if "S14" in options:
            modifier = modifier + '_S14'
            title_mod = title_mod + 'S14 '
        if "Bruno2017" in options:
            modifier = modifier + '_Bruno2017'
            title_mod = title_mod + 'Bruno2017 '

        #plot integral fluxes (either input or estimated)
        nthresh = len(flux_thresholds)
        figname = str(sep_year) + '_' + str(sep_month)+ '_' + str(sep_day) \
                + '_' + experiment + '_' + flux_type + modifier \
                + '_' + 'Event_Def'
        if experiment == 'user' and model_name != '':
            figname = str(sep_year) + '_' + str(sep_month)+ '_' + str(sep_day) \
                    + '_' + model_name + '_' + flux_type + modifier \
                    + '_' + 'Event_Def'
        if umasep or nthresh > 4:
            fig = plt.figure(figname,figsize=(9,9))
        else:
            fig = plt.figure(figname,figsize=(9,7))
        for i in range(nthresh):
            if crossing_time[i] == 0:
                continue
            data_label = (experiment + ' >'+ plt_energy[i] + ' '
                            + energy_units)
            plot_title = 'Threshold crossings for ' + experiment + '\n ' \
                            + title_mod + ' ' + flux_type + ' Fluxes '
            if experiment == 'user' and model_name != '':
                data_label = (model_name + ' >' + plt_energy[i] + ' '
                            + energy_units)
                plot_title = 'Threshold crossings for ' + model_name + '\n ' \
                                + title_mod + ' ' + flux_type + ' Fluxes '

            if flux_type == 'differential':
                data_label = (experiment + ' Estimated >' + plt_energy[i] \
                                + ' ' + energy_units)
                if experiment == 'user' and model_name != '':
                    data_label = (model_name + ' Estimated >' + plt_energy[i] \
                                + ' ' + energy_units)

            if plot_diff_thresh[i]: #differential threshold tacked on to end
                data_label = (experiment + ' ' + plt_energy[i] + ' '
                                + energy_units)
                if experiment == 'user' and model_name != '':
                    data_label = (model_name + ' ' + plt_energy[i] + ' '
                                + energy_units)
            ax = plt.subplot(nthresh, 1, i+1)
            #Don't want to plot negative values, particularly in background-subtracted plots
            if doBGSub:
                maskfluxes = np.ma.masked_where(integral_fluxes[i] <0, \
                                integral_fluxes[i])
                plt.plot_date(dates,maskfluxes,'-',label=data_label)
            else:
                plt.plot_date(dates,integral_fluxes[i],'-',label=data_label)

            plt.axvline(crossing_time[i],color='black',linestyle=':')
            plt.axvline(event_end_time[i],color='black',linestyle=':',
                        label="Start, End")
            plt.axhline(flux_thresholds[i],color='red',linestyle=':',
                        label="Threshold")
            if onset_peak[i] != None and onset_date[i] != None:
                plt.plot_date(onset_date[i],onset_peak[i],'o',color="black",
                        label="Onset Peak")
            if peak_time[i] != None and peak_flux[i] != None:
                plt.plot_date(peak_time[i],peak_flux[i],'ro',mfc='none',
                        label="Max Flux")
            if umasep:
                for k in range(len(umasep_times[i])):
                    plt.plot_date(umasep_times[i][k],umasep_fluxes[i][k],'bo')

            plt.xlabel('Date')
            plt.ylabel('Integral Flux\n' + '[' + flux_units_integral + ']')
            plt.suptitle(plot_title)
            if plot_diff_thresh[i]:
                plt.ylabel('Differential Flux\n' + '['
                        + flux_units_differential + ']')
            if sum(integral_fluxes[i]) > 0:
                plt.yscale("log")
            #ymin = max(1e-6, min(integral_fluxes[i]))
            # plt.ylim(ymin, peak_flux[i]+peak_flux[i]*.2)
            ax.legend(loc='upper right')
        if saveplot:
            fig.savefig(plotpath + '/' + figname + '.png')
        if not showplot:
            plt.close(fig)


        #All energy channels in specified date range with event start and stop
        print("Generating figure of fluxes in original energy bins. Any bad data "
              "points were interpolated. Lines indicate event start and stop for "
              "thresholds.")
        #Plot all channels of user specified data
        figname = str(sep_year) + '_' + str(sep_month)+ '_' + str(sep_day) \
                + '_' + experiment + '_' + flux_type + modifier \
                + '_' + 'All_Bins'
        if experiment == 'user' and model_name != '':
            figname = str(sep_year) + '_' + str(sep_month)+ '_' + str(sep_day) \
                    + '_' + model_name + '_' + flux_type + modifier \
                    + '_' + 'All_Bins'
        fig = plt.figure(figname,figsize=(9,4))
        ax = plt.subplot(111)
        nbins = len(energy_bins)
        for i in range(nbins):
            legend_label = ""
            if energy_bins[i][1] != -1:
                legend_label = str(energy_bins[i][0]) + '-' \
                               + str(energy_bins[i][1]) + ' ' + energy_units
            else:
                legend_label = '>'+ str(energy_bins[i][0]) + ' ' + energy_units

            if doBGSub:
                maskfluxes = np.ma.masked_where(fluxes[i] <0, fluxes[i])
                ax.plot_date(dates,maskfluxes,'-',label=legend_label)
            else:
                ax.plot_date(dates,fluxes[i],'-',label=legend_label)

        colors = ['black','red','blue','green','cyan','magenta','violet',\
                'orange','brown']
        for j in range(len(energy_thresholds)):
            if crossing_time[j] == 0:
                continue
            line_label = '>' + plt_energy[j] + ' MeV, ' \
                        + plt_flux[j] + ' pfu'
            if plot_diff_thresh[j]: #tacked on to end
                line_label = (plt_energy[j] + ' MeV, ' + plt_flux[j] + \
                            '\n' + flux_units_differential)
            ax.axvline(crossing_time[j],color=colors[j],linestyle=':',
                        label=line_label)
            ax.axvline(event_end_time[j],color=colors[j],linestyle=':')
        if flux_type == "integral":
            plt.ylabel('Integral Flux [' + flux_units_integral + ']')
            plt.title(experiment + ' '+ title_mod + '\n'\
                        + "Integral Energy Bins with Threshold Crossings")
            if experiment == 'user' and model_name != '':
                plt.title(model_name + ' '+ title_mod + '\n'\
                        + "Integral Energy Bins with Threshold Crossings")
        if flux_type == "differential":
            plt.ylabel('Flux [' + flux_units_differential + ']')
            plt.title(experiment + ' ' + title_mod + '\n'\
                        + "Differential Energy Bins with Threshold Crossings")
            if experiment == 'user' and model_name != '':
                plt.title(model_name + ' ' + title_mod + '\n' \
                        + "Differential Energy Bins with Threshold Crossings")
        plt.xlabel('Date')
        plt.yscale("log")
        chartBox = ax.get_position()
        ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.85,
                         chartBox.height])
        ax.legend(loc='upper center', bbox_to_anchor=(1.17, 1.05))
        if saveplot:
            fig.savefig(plotpath + '/' +figname + '.png')
        if not showplot:
            plt.close(fig)


        #Event-integrated fluence for energy channels
        print("Generating figure of event-integrated fluence spectrum.")
        #Plot fluence spectrum summed between SEP start and end dates
        figname = str(sep_year) + '_' + str(sep_month)+ '_' + str(sep_day) \
                + '_' + experiment + '_' + flux_type + modifier \
                + '_' + 'Fluence'
        if experiment == 'user' and model_name != '':
            figname = str(sep_year) + '_' + str(sep_month)+ '_' + str(sep_day) \
                    + '_' + model_name + '_' + flux_type + modifier \
                    + '_' + 'Fluence'
        fig = plt.figure(figname,figsize=(6,5))
        ax = plt.subplot(111)
        markers = ['o','P','D','v','^','<','>','*','d','+']
        for j in range(len(energy_thresholds)):
            if crossing_time[j] == 0:
                continue
            legend_label = '>' + plt_energy[j] + ' ' + energy_units + ', '\
                        + plt_flux[j] + ' ' + flux_units_integral
            if plot_diff_thresh[j]: #tacked on to end
                legend_label = (plt_energy[j] + ' ' + energy_units + ',\n'
                            + plt_flux[j] + '\n' + flux_units_differential)
            ax.plot(all_energies[j,:],all_fluence[j,:],markers[j],
                    color=colors[j],mfc='none',label=legend_label)
        plt.grid(which="both", axis="both")
        plt.title(experiment + ' ' + title_mod + '\n Event-Integrated Fluences '
                    'for All Event Definitions')
        if experiment == 'user' and model_name != '':
            plt.title(model_name + ' ' + title_mod + '\n Event-Integrated '
                    'Fluences for All Event Definitions')
        plt.xlabel('Energy [' + energy_units +']')
        if flux_type == "integral":
            plt.ylabel('Integral Fluxes [' + flux_units_integral + ']')
        if flux_type == "differential":
            plt.ylabel('Flux [' + flux_units_differential + ']')
        plt.xscale("log")
        plt.yscale("log")
        ax.legend(loc='upper right')
        if saveplot:
            fig.savefig(plotpath + '/' + figname + '.png')
        if not showplot:
            plt.close(fig)

    return sep_year, sep_month, sep_day, jsonfname


if __name__ == "__main__":
    #INPUTS:
    #   Start and end dates of SEP event
    #   Experiment is the source of the dataset - GOES or SEPEM:
    #       The energy bins and input file format will be determined
    #       by the 1) selected experiment, 2) SEP dates, and 3) FluxType
    #   FluxType can be differential or integral. If differential is specified
    #       and ThresholdType is chosen as 0 (>10 MeV exceeds 10 pfu, >100 MeV
    #       exceeds 1 pfu), then differential bins will be converted to
    #       integral flux
    #   showplot: set this flag to plot the SEP fluxes or background-
    #       subtracted SEP fluxes; also plot fluence
    parser = argparse.ArgumentParser()
    parser.add_argument("--StartDate", type=str, default='',
            help=("Start date in YYYY-MM-DD or \"YYYY-MM-DD HH:MM:SS\""
                   " with quotes"))
    parser.add_argument("--EndDate", type=str, default='',
            help=("End date in YYYY-MM-DD or \"YYYY-MM-DD HH:MM:SS\""
                    " with quotes"))
    parser.add_argument("--Experiment", type=str, choices=['GOES-08',
            'GOES-10', 'GOES-11', 'GOES-12', 'GOES-13', 'GOES-14', 'GOES-15',
            'SEPEM', 'SEPEMv3', 'EPHIN', 'EPHIN_REleASE', 'user'],
            default='', help="Enter name of spacecraft or dataset")
    parser.add_argument("--FluxType", type=str, choices=['integral',
            'differential'], default='',
            help=("Do you want to use integral or differential fluxes?"))
    parser.add_argument("--ModelName", type=str, default='', help=("If you "
            "chose user for experiment, specify the name of the model or "
            "experiment that you are analyzing (no spaces)."))
    parser.add_argument("--UserFile", type=str, default='tmp.txt', help=("If "
            "you chose user for experiment, specify the filename containing "
            "the fluxes. Specify energy bins and delimeter in code at top. "
            "Default is tmp.txt."))
    parser.add_argument("--spase_id", type=str, default='', help=("If your "
            "model or data source has an associated Spase ID, specify here."))
    parser.add_argument("--Threshold", type=str, default="",
            help=("Additional energy and flux threshold which will be used to "
                    "define the event. To define an integral flux threshold: "
                    "write 100,1 with no spaces; e.g. 100,1 indicates >100 MeV "
                    "fluxes crossing 1 pfu (1/[cm^2 s sr])."
                    "To define a differential flux threshold: write 25-40.9,0.01"
                    "with no spaces; e.g. energy bin "
                    "low edge-high edge,threshold (1/[MeV/n cm^2 s sr])). "
                    "Multiple thresholds may be entered separated by a "
                    "semi-colon with no spaces and surrounded by quotes, "
                    "e.g. \"30,1;50,1;25-40.9,0.001\""))
    parser.add_argument("--options", type=str, default='', help=("You "
            "may specify a series of options as a semi-colon separated list "
            "surrounded by quotations.\n"
            "\"uncorrected\" for GOES uncorrected differential fluxes with "
            "nominal GOES energy bins,\n "
            "\"S14\" to apply Sandberg et al. (2014) "
            "effective energies to GOES uncorrected fluxes for P2-P6,\n "
            "\"Bruno2017\" to apply Bruno (2017) effective energies to GOES-13 "
            "or GOES-15 P6-P11 channels for either corrected or uncorrected "
            "GOES fluxes. \n"
            "If both S14 and Bruno2017 are "
            "specified for GOES-13 or GOES-15, S14 bins will be applied to "
            "P2-P5 and Bruno2017 bins will be applied to P6-P11 for uncorrected "
            "fluxes.\n"
            "e.g. \"uncorrected;S14;Bruno2017\""))
    parser.add_argument("--NoInterp",
            help=("Do not fill in negative or missing fluxes via "
                    "linear interpolation in time. Set as None values "
                    "instead."), action="store_true")
    parser.add_argument("--SubtractBG",
            help="Set to calculate the background and subtract from the "\
                "SEP flux. Must define start and end dates for the background.",
                action="store_true")
    parser.add_argument("--BGStartDate", type=str, default='',
            help=("Start date in YYYY-MM-DD or \"YYYY-MM-DD HH:MM:SS\""
                   " with quotes to define the background time period."))
    parser.add_argument("--BGEndDate", type=str, default='',
            help=("End date in YYYY-MM-DD or \"YYYY-MM-DD HH:MM:SS\""
                    " with quotes to define the background time period."))
    parser.add_argument("--showplot",
            help="Flag to display plots", action="store_true")
    parser.add_argument("--saveplot",
            help="Flag to save plots to file", action="store_true")
    parser.add_argument("--DetectPreviousEvent",
            help=("Flag to indicate that the threshold is crossed at the "
                   "first point due to a previous event."), action="store_true")
    parser.add_argument("--TwoPeaks",
            help=("Flag to indicate that the event exceeds threshold (usually "
                    "for a couple of points), then drops below threshold "
                    "before increasing again to the true peak of the event."),
                    action="store_true")
    parser.add_argument("--UMASEP",
            help=("Flag to calculate flux values and thresholds specific to "
                "the UMASEP model. Thresholds for >10, >30, >50, >100 MeV and "
                "flux values at various time periods after crossing "
                "thresholds."), action="store_true")


    args = parser.parse_args()

    str_startdate = args.StartDate
    str_enddate = args.EndDate
    experiment = args.Experiment
    flux_type = args.FluxType
    model_name = args.ModelName
    user_file = args.UserFile
    spase_id = args.spase_id
    str_thresh = args.Threshold
    doBGSub = args.SubtractBG
    str_bgstartdate = args.BGStartDate
    str_bgenddate = args.BGEndDate
    showplot = args.showplot
    saveplot = args.saveplot
    detect_prev_event = args.DetectPreviousEvent
    two_peaks = args.TwoPeaks
    umasep = args.UMASEP
    options = args.options
    nointerp = args.NoInterp



    sep_year, sep_month, sep_day, jsonfname = run_all(str_startdate,
        str_enddate, experiment,
        flux_type, model_name, user_file, spase_id, showplot, saveplot,
        detect_prev_event, two_peaks, umasep, str_thresh, options, doBGSub,
        str_bgstartdate, str_bgenddate, nointerp)

    if showplot: plt.show()
