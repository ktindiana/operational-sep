import operational_sep_quantities as sep
from importlib import reload
import argparse
import csv
import datetime
import logging
import sys
import os
import asciitable

__version__ = "0.3"
__author__ = "Katie Whitman"
__maintainer__ = "Katie Whitman"
__email__ = "kathryn.whitman@nasa.gov"

#Changes in 0.2: Modified so that output list files will indicate when an
#   observation or flux did not exceed a certain threshold for a given SEP
#   event. Added a column specifying SEP date to sep_list
#2021-01-14, Changes in 0.3: Made consistent with operational_sep_quantities.py
#   v2.3 which includes background subtraction and various energy bin options.
#   Added more fields to list file to allow better specification of each data
#   set.

"""This and supporting codes are found in the respository:
        https://github.com/ktindiana/operational-sep

    This code will run operational_sep_quantities.py for multiple SEP events
    and record various quality flags for each event. The quality flags can
    be used to identify events that may need to be rerun with different options
    or times:
    FirstStart - start on first time point, likely previous ongoing event
    LastEnd - ended on last time point, need to increase time range
    ShortEvent - event less than 12 hours, need to check if actually continues
    LateHundred - >100 MeV crossed 24 hrs after >10 MeV, may be different events

    The SEP dates will be read in from a csv file with the columns:
    SEP Start Date - YYYY-MM-DD or YYYY-MM-DD HH:MM:SS
    Experiment - GOES-08 up to GOES-15 or SEPEM
    Ndays - number of days to end date (start time + Ndays)
    Flags - Options are blank, TwoPeak, or DetectPreviousEvent
    Example:
    2012-01-23 00:00:00,GOES-13,5,
    2012-01-27 00:00:00,GOES-13,5,DetectPreviousEvent
    2012-03-07 00:00:00,GOES-13,6,
    2012-03-13 12:00:00,GOES-13,5,

    Information for each SEP date is aggregated and sorted into separate lists
    for each threshold range. The output file columns are:
    Start Time,Peak Flux [pfu],Peak Time,End Time,Ts + 3hr,Ts + 4hr,Ts + 5hr,
    Ts + 6hr,Ts + 7hr
"""
datapath = 'data'
outpath = 'output'
listpath = 'lists'

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger('sep')
logging.getLogger("matplotlib").setLevel(logging.WARNING)

############## SET INPUTS ##################
showplot = False
saveplot = True
detect_prev_event_default = False #Set to true if get FirstStart flag
two_peaks_default = False #Set to true if get ShortEvent flag
############## END INPUTS #################

#READ IN SEP DATES FROM A CSV FILE
#ALSO READ IN ASSOCIATED EXPERIMENT FOR EACH DATE
def read_sep_dates(sep_filename):
    ''' Reads in a csv list file of SEP events. List must have the format:
        StartDate, Enddate, Experiment, FluxType, Flags,,,options,bgstartdate,
            bgenddate
        If the experiment is 'user', indicating a user-input flux file, then
        the file must have the format:
        StartDate, Enddate, Experiment, FluxType, Flags, Model Name,
            User Filename, options, bgstartdate, bgenddate

        Flags may be: TwoPeak, DetectPreviousEvent, SubtractBG
        options may be: "S14,Bruno2017,uncorrected"
    '''
    print('Reading in file ' + sep_filename)
    start_dates = [] #row 0
    end_dates = []
    experiments = [] #row 1, e.g. GOES-11, GOES-13, GOES-15, SEPEM, user
    flux_types = [] #row 3
    flags = [] #row 4
    model_names = [] #row 5
    user_files = [] #row 6
    options = [] #row 7
    bgstartdate = [] #row 8
    bgenddate = [] #row 9

    with open(sep_filename) as csvfile:
        readCSV = csv.reader(csvfile, delimiter=',')
        #Define arrays that hold dates
        for row in readCSV:
            if len(row[0]) > 10:
                stdate = datetime.datetime.strptime(row[0][0:19],
                                            "%Y-%m-%d %H:%M:%S")
            if len(row[0]) == 10:
                stdate = datetime.datetime.strptime(row[0][0:10],
                                            "%Y-%m-%d")
            if len(row[1]) > 10:
                enddate = datetime.datetime.strptime(row[1][0:19],
                                            "%Y-%m-%d %H:%M:%S")
            if len(row[1]) == 10:
                enddate = datetime.datetime.strptime(row[1][0:10],
                                            "%Y-%m-%d")
            start_dates.append(str(stdate))
            end_dates.append(str(enddate))
            experiments.append(row[2])
            flux_types.append(row[3])

            if len(row) > 4:
                flags.append(row[4])
            else:
                flags.append('')

            if len(row) > 5:
                model_names.append(row[5])
            else:
                model_names.append('')

            if len(row) > 6:
                user_files.append(row[6])
            else:
                user_files.append('')

            if len(row) > 7:
                options.append(row[7])
            else:
                options.append('')

            if len(row) > 8:
                bgstartdate.append(row[8])
            else:
                bgstartdate.append('')

            if len(row) > 9:
                bgenddate.append(row[9])
            else:
                bgenddate.append('')


            if row[1] == 'user':
                if len(row) < 7:
                    sys.exit("For a user file, you must specify model name and "
                            "input filename in the list.")


    return start_dates, end_dates, experiments, flux_types, flags, model_names,\
        user_files, options, bgstartdate, bgenddate


def write_sep_lists(sep_year, sep_month, sep_day, experiment, flux_type,
                    model_name, umasep, threshold):
    '''Reads in sep_values_* files output by operational_sep_quantities.
        Selected information is taken and sorted into lists for each threshold
        definition. Output is then an SEP list with associated quantities for
        each threshold.

        If the operational thresholds or the user input threshold is not
        present in the sep_values_ file, then values of None will be saved
        for all columns for that date. This indicates in the output list file
        that the model or observations did not cross those thresholds.
    '''
    #Put together thresholds that might be found inside sep_values files
    op_thresh = [['10','10'],['100','1']]
    if threshold != '100,1':
        thresh = threshold.split(';')
        for th in thresh:
            th = th.split(",")
            op_thresh.append([str(float(th[0])),str(float(th[1]))])
    is_thresh = [False]*len(op_thresh)

    infile = outpath + '/' + 'sep_values_' + experiment + '_' + flux_type \
                + '_' + str(sep_year) + '_' + str(sep_month) + '_' \
                + str(sep_day) + '.csv'
    if experiment == 'user' and model_name != '':
        infile = outpath + '/' + 'sep_values_' + model_name + '_' + flux_type \
                    + '_' + str(sep_year) + '_' + str(sep_month) + '_' \
                    + str(sep_day) + '.csv'
    isgood = os.path.isfile(infile)
    if not isgood:
        print('Cannot open file ' + infile)
        return False
    data = open(infile).read()
    lines = data.split('\n')
    #Pick out columns to extract and save to SEP list
    #start time, onset peak, onset time, peak flux, peak time, end time
    #If UMASEP, then all delayed proton values
    cols = [2,3,4,5,6,8]
    umacols = [15,17,19,21,23]

    for line in lines:
        if line == '': continue
        if line[0] == '#': continue
        row = line.split(',')
        thresh = [row[0][1:len(row[0])], row[1]]
        threshfile = listpath + '/' +'sep_list_' + str(thresh[0]) + 'MeV_' \
                    + str(thresh[1]) + 'pfu.csv'
        isgood = os.path.isfile(threshfile)
        if not isgood:
            print('Cannot open file ' + threshfile)
            continue
        #Check which threshold are and are not included in the sep_values_ file
        for i in range(len(op_thresh)):
            if thresh[0] == op_thresh[i][0] and thresh[1] == op_thresh[i][1]:
                is_thresh[i] = True

        fin = open(threshfile,'a')
        if experiment == 'user' and model_name != '':
            fin.write(model_name + ',')
        if experiment != 'user':
            fin.write(experiment + ',')
        date = '{0:d}-{1:02d}-{2:02d}'.format(sep_year, sep_month,sep_day)
        fin.write(date + ',')
        for col in cols:
            fin.write(str(row[col]) + ',')

        if umasep:
            for ucol in umacols:
                fin.write(str(row[ucol]) + ',')

        fin.write('\n')
        fin.close()

    #If a threshold is not present in a file for the requested date, then
    #thresholds were not crossed. Record None for those values to indicate
    #the event did not cross threshold.
    for i in range(len(op_thresh)):
        if not is_thresh[i]:   #False, missing threshold
            threshfile = listpath + '/' +'sep_list_' + str(op_thresh[i][0])  \
                        + 'MeV_' + str(op_thresh[i][1]) + 'pfu.csv'
            isgood = os.path.isfile(threshfile)
            if not isgood:
                print('Cannot open file ' + threshfile)
                continue

            fin = open(threshfile,'a')
            if experiment == 'user' and model_name != '':
                fin.write(model_name + ',')
            if experiment != 'user':
                fin.write(experiment + ',')
            date = '{0:d}-{1:02d}-{2:02d}'.format(sep_year, sep_month,sep_day)
            fin.write(date + ',')
            for col in cols:
                fin.write('None' + ',')

            if umasep:
                for ucol in umacols:
                    fin.write('None' + ',')

            fin.write('\n')
            fin.close()

    return True


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--Filename", type=str, default='tmp.csv', \
            help=("Name of csv file containing list of SEP start dates."
            "Default is tmp.csv."))
    parser.add_argument("--OutFilename", type=str, default='lists/out.csv', \
            help=("Name of csv file containing list of SEP dates with "
                "flags indicating status after run with "
                "operational_sep_quantities.py. Default is lists/out.csv."))
    parser.add_argument("--Threshold", type=str, default="100,1",
            help=("An additional energy and flux threshold (written as 100,1 "
                    "with no spaces) which will be used to define the event. "
                    "e.g. 100,1 indicates >100 MeV fluxes crossing 1 pfu "
                    "(1/[cm^2 s sr]). Default = '100,1'"))
    parser.add_argument("--UMASEP",
            help=("Flag to calculate flux values and thresholds specific to "
                "the UMASEP model. Thresholds for >10, >30, >50, >100 MeV and "
                "flux values at 3, 4, 5, 6, 7 hours after "
                "crossing thresholds."), action="store_true")

    args = parser.parse_args()
    sep_filename = args.Filename
    outfname = args.OutFilename
    threshold = args.Threshold
    umasep = args.UMASEP

    #Define filenames for aggregated SEP lists
    #MAKES NEW FILE EVERY TIME THIS CODE IS RUN
    #Always >10 MeV, 10 pfu and >100 MeV, 1 pfu thresholds
    open(listpath + '/' 'sep_list_10MeV_10pfu.csv','w+').close()
    open(listpath + '/' 'sep_list_100MeV_1pfu.csv','w+').close()
    if umasep:
        open(listpath + '/' 'sep_list_30MeV_1pfu.csv','w+').close()
        open(listpath + '/' 'sep_list_50MeV_1pfu.csv','w+').close()
    #Additional threshold
    if threshold != '100,1':
        thresh = threshold.split(';')
        for th in thresh:
            th = th.split(",")
            open(listpath + '/' 'sep_list_' + str(float(th[0])) + 'MeV_' \
                + str(float(th[1])) + 'pfu.csv','w+').close()

    #READ IN SEP DATES AND experiments
    start_dates, end_dates, experiments, flux_types, flags, model_names, \
        user_files, options, bgstart, bgend = read_sep_dates(sep_filename)

    #Prepare output file listing events and flags
    fout = open(outfname,"w+")
    fout.write('#Experiment,SEP Date,FirstStart,LastEnd,ShortEvent,LateHundred,'
                + 'Exception\n')

    #---RUN ALL SEP EVENTS---
    Nsep = len(start_dates)
    print('Read in ' + str(Nsep) + ' SEP events.')
    for i in range(Nsep):
        start_date = start_dates[i]
        end_date = end_dates[i]
        experiment = experiments[i]
        flux_type = flux_types[i]
        flag = flags[i]
        model_name = model_names[i]
        user_file = user_files[i]
        option = options[i]
        bgstartdate = bgstart[i]
        bgenddate = bgend[i]

        flag = flag.split(';')
        detect_prev_event = detect_prev_event_default
        two_peaks = two_peaks_default
        doBGSub = False
        if "DetectPreviousEvent" in flag:
            detect_prev_event = True
        if "TwoPeak" in flag:
            two_peaks = True
        if "SubtractBG" in flag:
            doBGSub = True

        print('\n-------RUNNING SEP ' + start_date + '---------')
        #CALCULATE SEP INFO AND OUTPUT RESULTS TO FILE
        try:
            FirstStart, LastEnd, ShortEvent, LateHundred, sep_year, sep_month, \
            sep_day = sep.run_all(start_date, end_date, experiment, flux_type,
                model_name, user_file, showplot, saveplot, detect_prev_event,
                two_peaks, umasep, threshold, option, doBGSub, bgstartdate,
                bgenddate)

            sep_date = datetime.datetime(year=sep_year, month=sep_month,
                            day=sep_day)
            if experiment == 'user' and model_name != '':
                fout.write(model_name + ',')
            if experiment != 'user':
                fout.write(experiment + ',')
            fout.write(str(sep_date)+','+ str(FirstStart)+',' +str(LastEnd)\
                        +','+str(ShortEvent) + ',' + str(LateHundred) + ', ')
            fout.write('\n')

            success=write_sep_lists(sep_year, sep_month, sep_day, experiment,
                        flux_type, model_name, umasep, threshold)
            if not success:
                print('Could not write values to file for ' + str(sep_year) \
                    + '/' + str(sep_month) + '/' + str(sep_day))

            sep = reload(sep)

        except SystemExit as e:
            # this log will include traceback
            logger.exception('operational_sep_quantities failed with exception')
            # this log will just include content in sys.exit
            logger.error(str(e))
            fout.write(str(start_date) +', , , , ,' + '\"' + str(e) + '\"' )
            fout.write('\n')
            sep = reload(sep)
            continue

    fout.close()
