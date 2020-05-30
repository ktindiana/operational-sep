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

__version__ = "0.6"
__author__ = "Katie Whitman"
__maintainer__ = "Katie Whitman"
__email__ = "kathryn.whitman@nasa.gov"

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

#See full program description in all_program_info() below
datapath = 'data'
outpath = 'output'
badval = -1 #bad data points will be set to this value; must be negative

######FOR USER DATA SETS######
#(expect the first column contains date in YYYY-MM-DD HH:MM:SS format)
#Identify columns containing fluxes you want to analyze
user_col = arr.array('i',[1,2,3,4,5])

#DELIMETER between columns; for whitespace separating columns, use " " or ""
user_delim = ","

#DEFINE ENERGY BINS associated with user file and columns specified above as:
#   [[Elow1,Ehigh1],[Elow2,Ehigh2],[Elow3,Ehigh3],etc]
#Use -1 in the second edge of the bin to specify integral channel (infinity):
#   [[Elow1,-1],[Elow2,-1],[Elow3,-1],etc]
user_energy_bins = [[10,-1],[30,-1],[40,-1],[50,-1],[100,-1]]
#EPREM [[10,-1],[30,-1],[40,-1],[50,-1],[100,-1]]
#SEPMOD [[300,-1],[100,-1],[60,-1],[50,-1],[30,-1],[10,-1]]
############################

#FILENAME(s) containing user input fluxes - WILL BE SET THROUGH ARGUMENT
#Can be list of files that are continuous in time
 #      e.g. user_fname = ["file1.txt","file2.txt"]
user_fname = ['tmp.txt']


def all_program_info(): #only for documentation purposes
    """This program will calculate various useful pieces of operational
    information about SEP events from GOES-08, -10, -11, -12, -13, -14, -15
    data and the SEPEM (RSDv2) dataset.

    SEP event values are always calculated for threshold definitions:
        >10 MeV exceeds 10 pfu
        >100 MeV exceed 1 pfu

    The user may add an additional threshold through the command line.
    This program will check if data is already present in a 'data' directory. If
    not, GOES data will be automatically downloaded from NOAA ftp site. SEPEM
    (RSDv2) data must be downloaded by the user and unzipped inside the 'data'
    directory. Because the SEPEM data set is so large (every 5 minutes from 1974
    to 2015), the program will break up the data into yearly files for faster
    reading.

    The values calculated here are important for space radiation operations:
       Onset time, i.e. time to cross thresholds
       Peak intensity
       Time of peak intensity
       Rise time (onset to peak)
       End time, i.e. fall below 0.85*threshold for 3 points (15 mins for GOES)
       Duration
       Event-integrated fluences

    User may choose differential proton fluxes (e.g. [MeV s sr cm^2]^-1) or
    integral fluxes (e.g. [s sr cm^2]^-1). The program has no internal checks or
    requirements on units - EXCEPT FOR THE THRESHOLD DEFINITIONS OF >10, 10
    and >100, 1. If you change those thresholds in the main program accordingly,
    you should be able to use other units. Also, all of the plots and messages
    refer to MeV, pfu, and cm.

    If a previous event is ongoing and the specified time period starts with a
    threshold already crossed, you may try to set the --DetectPreviousEvent
    flag. If the flux drops below threshold before the next event starts, the
    program will identify the second event. This will only work if the
    threshold is already crossed for the very first time in your specified
    time period, and if the flux drops below threshold before the next event
    starts.

    If the event has an initial increase above threshold for a few points, falls
    below threshold, then continues to increase above threshold again, you
    may try to use the --TwoPeaks feature to capture the full duration of the
    event. The initial increase above threshold must be less than a day. An
    example of this scenario can be seen in >100 MeV for 2011-08-04.

    RUN CODE FROM COMMAND LINE (put on one line), e.g.:
    python3 operational_sep_quantities.py --StartDate 2012-05-17
        --EndDate '2012-05-19 12:00:00' --Experiment GOES-13
        --FluxType integral --showplot

    RUN CODE FROM COMMAND FOR USER DATA SET (put on one line), e.g.:
    python3 operational_sep_quantities.py --StartDate 2012-05-17
        --EndDate '2012-05-19 12:00:00' --Experiment user --ModelName MyModel
        --UserFile MyFluxes.txt --FluxType integral --showplot

    RUN CODE IMPORTED INTO ANOTHER PYTHON PROGRAM, e.g.:
    import operational_sep_quantities as sep
    start_date = '2012-05-17'
    end_date = '2012-05-19 12:00:00'
    experiment = 'GOES-13'
    flux_type = 'integral'
    model_name = '' #if experiment is user, set model_name to describe data set
    user_file = '' #if experiment is user, specify filename containing fluxes
    showplot = True  #Turn to False if don't want to see plots
    detect_prev_event = True  #Helps if previous event causes high intensities
    two_peaks = False  #Helps if two increases above threshold in one event
    threshold = '100,1' #default; modify to add a threshold to 10,10 and 100,1
    sep.run_all(start_date, end_date, experiment, flux_type, model_name,
        user_file, showplot, detect_prev_event, two_peaks, threshold)

    Set the desired directory locations for the data and output at the beginning
    of the program in datapath and outpath. Defaults are 'data' and 'output'.

    In order to calculate the fluence, the program determines time_resolution
    (seconds) from two (fairly random) data points at the start of the SEP
    event. GOES and SEPEM data sets have a time resolution of 5 minutes. If the
    user wishes to use a data set with measurements at irregular times, then the
    subroutine calculate_fluence should be modified.

    OUTPUT: This program outputs 3 to 4 files, 1 per defined threshold plus
    a summary file containing all of the values calculated for each threshold.
    A file named as e.g. fluence_GOES-13_differential_gt10_2012_3_7.csv contains
    the event-integrated fluence for each energy channel using the specified
    threshold (gt10) to determine start and stop times.
    A file named as e.g. sep_values_GOES-13_differential_2012_3_7.csv contains
    start time, peak flux, etc, for each of the defined thresholds.

    Added functionality to ouput the >10 MeV and >100 MeV time series for the
    date range input by the user. If the original data were integral fluxes,
    then the output files simply contain the >10 and >100 MeV time series from
    the input files. If the original data were differential fluxes, then the
    estimated >10 and >100 MeV fluxes are output as time series.

    USER INPUT DATA SETS: Users may input their own data set. For example, if an
    SEP modeler would like to feed their own intensity time series into this
    code and calculate all values in exactly the same way they were calculated
    for data, it is possible to do that. Fluxes should be in units of
    1/[MeV cm^2 s sr] or 1/[cm^2 s sr] and energy channels in MeV for the plot
    labels to be correct. You can use any units, as long as you are consistent
    with energy units in energy channel/bin definition and in fluxes and you
    MODIFY THE THRESHOLD VALUES TO REFLECT YOUR UNITS. You may then want to
    modify plot labels accordingly if not using MeV and cm.
    NOTE: The first column in your flux file is assumed to be time in format
    YYYY-MM-DD HH:MM:SS. IMPORTANT FORMATTING!!
    NOTE: The flux file may contain header lines that start with a hash #,
    including blank lines.
    NOTE: Any bad or missing fluxes must be indicated by a negative value.
    NOTE: Put your flux file into the "datapath" directory.
    NOTE: Please use only differential or integral channels. Please do not mix
    them. You may have one integral channel in the last bin, as this is the way
    HEPAD works and the code has been written to include that HEPAD >700 MeV
    bin along with lower differential channels.

    USER VARIABLES: The user must modify the following variables at the very
    top of the code (around line 30):
        user_col - identify columns in your file containing fluxes to analyze;
                even if your delimeter is white space, consider the date-time
                column as one single column. SET AT TOP OF CODE.
        user_delim - delimeter between columns, e.g. " " or ","   Use " " for
                any amount of whitespace. SET AT TOP OF CODE.
        user_energy_bins - define your energy bins at the top of the code in the
                variable user_energy_bins. Follow the format in the subroutine
                define_energy_bins. SET AT TOP OF CODE.
        user_fname - specify the name of the file containing the fluxes
                through an argument in the command line. --UserFile  The
                user_fname variable will be updated with that filename. ARGUMENT
        time_resolution - will be calculated using two time points in your file;
                if you have irregular time measurements, calculate_fluence()
                must be modified/rewritten. AUTOMATICALLY DETERMINED.
    """


def check_paths():
    """Check that the paths that hold the data and output exist. If not, create.
    """
    print('Checking that paths exist: ' + datapath + ' and ' + outpath)
    if not os.path.isdir(datapath):
        print('check_paths: Directory containing fluxes, ' + datapath +
        ', does not exist. Creating.')
        os.mkdir(datapath);
    if not os.path.isdir(outpath):
        print('check_paths: Directory to store output information, ' + outpath
            + ', does not exist. Creating.')
        os.mkdir(outpath);


def make_yearly_files(filename):
    """Convert a large data set into yearly files."""
    print('Breaking up the SEPEM data into yearly data files. (This could '
            + 'take a while, but you will not have to do it again.)')
    fnamebase = filename.replace('.csv','')  #if csv file
    fnamebase = fnamebase.replace('.txt','')  #if txt file

    with open(datapath + '/' + filename) as csvfile:
        readCSV = csv.reader(csvfile, delimiter=',')
        has_header = csv.Sniffer().has_header(csvfile.readline())
        if has_header:
            next(readCSV)  # Skip single header row.
        ncol = len(next(readCSV))
        csvfile.seek(0) #back to beginning of file
        if has_header:
            header = csvfile.readline()  # save header row.

        check_year = 0
        for row in readCSV:
            date = datetime.datetime.strptime(row[0][0:18],
                                            "%Y-%m-%d %H:%M:%S")
            year = date.year
            if check_year != year:
                if check_year != 0:
                    outfile.close()
                    outfname = fnamebase + '_' + str(year) + '.csv'
                    outfile = open(datapath + '/' + outfname,'w+')
                    check_year = year

                if check_year == 0:
                    outfname = fnamebase + '_' + str(year) + '.csv'
                    outfile = open(datapath + '/' + outfname,'w+')
                    if has_header:
                        outfile.write(header)
                    check_year = year

            outfile.write(','.join(row))
            outfile.write('\n')

    outfile.close()
    csvfile.close()
    return


def check_data(startdate, enddate, experiment, flux_type):
    """Check that the files containing the data are in the data directory. If
        the files for the requested dates aren't present, they will be
        downloaded from the NOAA website. For SEPEM (RSDv2) data, if missing,
        the program prints the URL from which the data can be downloaded and
        unzipped manually.
        The RSDv2 data set is very large and takes a long time to read as a
        single file. This program will generate files containing fluxes for
        each year for faster reading.
    """
    print('Checking that the requested data is present on your computer.')
    styear = startdate.year
    stmonth = startdate.month
    stday = startdate.day
    endyear = enddate.year
    endmonth = enddate.month
    endday = enddate.day

    #Array of filenames that contain the data requested by the User
    filenames1 = []  #SEPEM, eps, or epead
    filenames2 = []  #hepad
    filenames_orien = []  #orientation flag for GOES-13+

    #If user wants to use own input file (filename defined at top of code)
    if experiment == "user":
        nuser = len(user_fname)
        for i in range(nuser):
            exists = os.path.isfile(datapath + '/' + user_fname[i])
            if exists:
                filenames1.append(user_fname[i])
            if not exists:
                sys.exit("You have selected to read a user-input input file "
                "with filename " + datapath + '/' + user_fname[i]
                + ". This file is not found! Exiting.")

        return filenames1, filenames2, filenames_orien

    #SEPEM data set is continuous, but too long; prefer yearly files
    #Check if user has yearly files; if not:
        #check if user has original SEPEM, then create yearly files
        #Otherwise alert user to download data set and try again
    if (experiment == "SEPEM"):
        year = styear
        while (year <= endyear):
            fname = 'SEPEM_H_reference_' + str(year) + '.csv'
            exists = os.path.isfile(datapath + '/' + fname)
            if exists:
                filenames1.append(fname)
                year = year + 1
            if not exists:
                full_exists = os.path.isfile(datapath
                                    + '/SEPEM_H_reference.txt')
                if not full_exists:
                    sys.exit("Please download and unzip the RSDv2 data set."
                        " You may download the file at"
                        " http://sepem.eu/help/SEPEM_RDS_v2-00.zip for full "
                        "fluxes or http://sepem.eu/help/SEPEM_RDS_v2-00.zip "
                        "for ESA background-subtracted fluxes.")
                if full_exists:
                    #Break up SEPEM data set into yearly files
                    print('The SEPEM (RSDv2) is more tractable when breaking'
                        + ' into yearly data files. Producing yearly files.')
                    make_yearly_files('SEPEM_H_reference.txt')
                    year = styear

        return filenames1, filenames2, filenames_orien

    #GOES data is stored in monthly data files
    get_years = []
    get_months = []
    test_year = styear
    test_month = stmonth
    test_date = datetime.datetime(year=test_year, month=test_month, day=1)
    while (test_date <= enddate):
        get_years.append(test_year)
        get_months.append(test_month)
        test_month = test_month + 1
        if (test_month > 12):
            test_month = 1
            test_year = test_year + 1
        test_date = datetime.datetime(year=test_year, month=test_month, day=1)

    NFILES = len(get_months)  #number of data files to download

    #Set correct file prefix for data files
    if experiment == "GOES-08":
        prefix1 = 'g08_eps_5m_'
        prefix2 = 'g08_hepad_5m_'
        satellite = 'goes08'

    if experiment == "GOES-10":
        prefix1 = 'g10_eps_5m_'
        prefix2 = 'g10_hepad_5m_'
        satellite = 'goes10'

    if experiment == "GOES-11":
        prefix1 = 'g11_eps_5m_'
        prefix2 = 'g11_hepad_5m_'
        satellite = 'goes11'

    if experiment == "GOES-12":
        prefix1 = 'g12_eps_5m_'
        prefix2 = 'g12_hepad_5m_'
        satellite = 'goes12'

    if experiment == "GOES-13":
        prefix2 = 'g13_hepad_ap_5m_'
        prefix_orien = 'g13_epead_orientation_flag_1m_'
        satellite = 'goes13'
        if flux_type == "differential":
            prefix1 = 'g13_epead_p17ew_5m_'
        if flux_type == "integral":
            prefix1 = 'g13_epead_cpflux_5m_'

    if experiment == "GOES-14":
        prefix2 = 'g14_hepad_ap_5m_'
        prefix_orien = 'g14_epead_orientation_flag_1m_'
        satellite = 'goes14'
        if flux_type == "differential":
            prefix1 = 'g14_epead_p17ew_5m_'
        if flux_type == "integral":
            prefix1 = 'g14_epead_cpflux_5m_'

    if experiment == "GOES-15":
        prefix2 = 'g15_hepad_ap_5m_'
        prefix_orien = 'g15_epead_orientation_flag_1m_'
        satellite = 'goes15'
        if flux_type == "differential":
            prefix1 = 'g15_epead_p17ew_5m_'
        if flux_type == "integral":
            prefix1 = 'g15_epead_cpflux_5m_'

    #for every month that data is required, check if file is present or
    #needs to be downloaded.
    for i in range(NFILES):
        year = get_years[i]
        month = get_months[i]
        last_day = calendar.monthrange(year,month)[1]
        date_suffix = '%i%02i01_%i%02i%02i' % (year,month,year,month,
                        last_day)
        fname1 = prefix1 + date_suffix + '.csv'
        exists1 = os.path.isfile(datapath + '/' + fname1)
        fname2 = prefix2 + date_suffix + '.csv'
        exists2 = os.path.isfile(datapath + '/' + fname2)
        if (experiment == "GOES-13" or experiment == "GOES-14"
            or experiment == "GOES-15"):
            fname_orien = prefix_orien + date_suffix + '_v1.0.0.csv'
            exists_orien = os.path.isfile(datapath + '/' + fname_orien)
            filenames_orien.append(fname_orien)

        filenames1.append(fname1)
        filenames2.append(fname2)

        if not exists1: #download file if not found on your computer
            url = ('https://satdat.ngdc.noaa.gov/sem/goes/data/avg/' +
                '%i/%02i/%s/csv/%s' % (year,month,satellite,fname1))
            print('Downloading GOES data: ' + url)
            try:
                urllib.request.urlopen(url)
                wget.download(url, datapath + '/' + fname1)
            except urllib.request.HTTPError:
                sys.exit("Cannot access file at " + url +
                ". Please check that selected spacecraft covers date range.")


        if not exists2: #download file if not found on your computer
            url = ('https://satdat.ngdc.noaa.gov/sem/goes/data/avg/' +
               '%i/%02i/%s/csv/%s' % (year,month,satellite,fname2))
            print('Downloading GOES data: ' + url)
            try:
                urllib.request.urlopen(url)
                wget.download(url, datapath + '/' + fname2)
            except urllib.request.HTTPError:
                sys.exit("Cannot access file at " + url +
               ". Please check that selected spacecraft covers date range.")

        if (experiment == "GOES-13" or experiment == "GOES-14"
            or experiment == "GOES-15"):
            if not exists_orien: #download file if not found on your computer
                url = ('https://satdat.ngdc.noaa.gov/sem/goes/data/avg/' +
                   '%i/%02i/%s/csv/%s' % (year,month,satellite,fname_orien))
                print('Downloading GOES data: ' + url)
                try:
                    urllib.request.urlopen(url)
                    wget.download(url, datapath + '/' + fname_orien)
                except urllib.request.HTTPError:
                    sys.exit("Cannot access orientation file at " + url +
                   ". Please check that selected spacecraft covers date range.")

    return filenames1, filenames2, filenames_orien


def find_goes_data_dimensions(filename):
    """Input open csv file of GOES data. Identifies the start of the data by
       searching for the string 'data:', then returns the number of header
       rows and data rows present in the file.
    """
    with open(datapath + '/' + filename) as csvfile:
        #GOES data has very large headers; figure out where the data
        #starts inside the file and skip the required number of lines
        nhead = 0
        for line in csvfile:
            nhead = nhead + 1
            if 'data:' in line: #location of line before column headers
                break
        nhead = nhead + 1 #account for column header line
        csvfile.readline() #proceed to column header line
        readCSV = csv.reader(csvfile, delimiter=',')
        nrow = len(csvfile.readlines()) #count lines of data
        #print('\nThere are ' + str(nhead) + ' header lines and '
        #    + str(nrow) + ' rows of data in ' + filename)

    csvfile.close()
    return nhead, nrow


def get_west_detector(filename, dates):
    """For GOES-13+, identify which detector is facing west from the
       orientation flag files. Get an orientation for each data point.
       EPEAD orientation flag. 0: A/W faces East and B/E faces West.
       1: A/W faces West and B/E faces East. 2: yaw-flip in progress.
    """
    nhead, nrow = find_goes_data_dimensions(filename)
    orien_dates = []
    orientation = []

    with open(datapath + '/' + filename) as orienfile:
        #GOES data has very large headers; figure out where the data
        #starts inside the file and skip the required number of lines
        readCSV = csv.reader(orienfile, delimiter=',')
        for k in range(nhead):
            next(readCSV)  #to line with column headers

        for row in readCSV:
            date = datetime.datetime.strptime(row[0][0:18],
                                            "%Y-%m-%d %H:%M:%S")
            orien_dates.append(date)
            orientation.append(float(row[1]))

    orienfile.close()

    #orientation data is in 1 minute intervals while flux data is in 5
    #minute intervals. Identify W facing detector for flux times.
    #Assume that time stamps match every 5 minutes.
    date_index = 0
    west_detector = []
    for i in range(nrow):
        if orien_dates[i] == dates[date_index]:
            if orientation[i] == 0:
                west_detector.append("B")
            if orientation[i] == 1:
                west_detector.append("A")
            if orientation[i] == 2:
                west_detector.append("Flip")
            date_index = date_index + 1
            if date_index == len(dates):
                break

    #print('There were ' + str(len(dates)) + ' input dates and there are ' +
    #        str(len(west_detector)) + ' detector orientations.')
    return west_detector


def read_in_files(experiment, flux_type, filenames1, filenames2,
                filenames_orien):
    """Read in the appropriate data files with the correct format. Return an
       array with dates and fluxes. Bad flux values (any negative flux) are set
       to -1. Format is defined to work with the files downloaded directly from
       NOAA or the RSDv2 (SEPEM) website as is.
       The fluxes output for the GOES-13+ satellites are always from the
       westward-facing detector (A or B) by referring to the orientation flags
       provided in the associated orientation file. Data taken during a yaw
       flip (orientation flag = 2) are excluded and fluxes are set to -1.
       Note that the EPS detectors on GOES-08 and -12 face westward. The
       EPS detector on GOES-10 faces eastward. GOES-11 is a spinning satellite.
    """
    print('Reading in data files for ' + experiment + '.')
    NFILES = len(filenames1)

    if experiment == "SEPEM":
        for i in range(NFILES):
            print('Reading in file ' + datapath + '/' + filenames1[i])
            with open(datapath + '/' + filenames1[i]) as csvfile:
                readCSV = csv.reader(csvfile, delimiter=',')
                has_header = csv.Sniffer().has_header(csvfile.readline())
                if has_header:
                    next(readCSV)  # Skip single header row.
                ncol = len(next(readCSV))
                csvfile.seek(0) #back to beginning of file
                if has_header:
                    next(readCSV)  # Skip header row.
                nrow = len(csvfile.readlines())
                #print('There are ' + str(ncol) + ' columns and ' + str(nrow) +
                #    ' rows of data in ' + filenames1[i])

                #Define arrays that hold dates and fluxes
                dates = []
                fluxes = np.zeros(shape=(ncol-1,nrow))

                csvfile.seek(0) #back to beginning of file
                if has_header:
                    next(readCSV)  # Skip header row.

                count = 0
                for row in readCSV:
                    date = datetime.datetime.strptime(row[0][0:18],
                                                    "%Y-%m-%d %H:%M:%S")
                    dates.append(date)
                    for j in range(1,ncol):
                        flux = float(row[j])
                        if flux < 0:
                            flux = badval
                        fluxes[j-1][count] = flux
                    count = count + 1
            #If reading in multiple files, then combine all data into one array
            #SEPEM currently only has one file, but making generalized
            if i==0:
                all_fluxes = fluxes
                all_dates = dates
            else:
                all_fluxes = np.concatenate((all_fluxes,fluxes),axis=1)
                all_dates = all_dates + dates

        return all_dates, all_fluxes

    #All GOES data
    if experiment != "SEPEM":
        if (experiment == "GOES-08" or experiment == "GOES-10" or
            experiment == "GOES-11" or experiment == "GOES-12"):
            if flux_type == "differential":
                #CORRECTED CHANNELS
                columns = [18,19,20,21,22,23] #eps
                hepad_columns = [1,2,3,4]
                #UNCORRECTED channels
                #columns = [5,6,7,8,9,10] #eps
                #hepad_columns = [1,2,3,4]
            if flux_type == "integral":
                columns = [25,26,27,28,29,30] #eps
                hepad_columns = [4]
        if (experiment == "GOES-13" or experiment == "GOES-14" or
            experiment == "GOES-15"):
            if flux_type == "differential":
                #CORRECTED CHANNELS
                columns = [16,24,32,40,48,56] #epead, default A detector
                columnsB = [12,20,28,36,44,52] #epead, B detector
                hepad_columns = [9,12,15,18]
                #UNCORRECTED CHANNELS
                #columns = [15,23,31,39,47,55] #epead, default A detector
                #columnsB = [11,19,27,35,43,51] #epead, B detector
            if flux_type == "integral":
                #ONLY CORRECTED CHANNELS AVAILABLE
                columns = [18,20,22,24,26,28] #epead, default A detector
                columnsB = [4,6,8,10,12,14] #epead, B detector
                hepad_columns = [18]

        ncol = len(columns)
        nhcol = len(hepad_columns)
        totcol = ncol + nhcol

        #Read in fluxes from files
        for i in range(NFILES):
            #FIRST set of files for lower energy eps or epead
            nhead, nrow = find_goes_data_dimensions(filenames1[i])
            dates = []
            fluxes = np.zeros(shape=(totcol,nrow))
            print('Reading in file ' + datapath + '/' + filenames1[i])
            with open(datapath + '/' + filenames1[i]) as csvfile:
                #GOES data has very large headers; figure out where the data
                #starts inside the file and skip the required number of lines
                readCSV = csv.reader(csvfile, delimiter=',')
                for k in range(nhead):
                    next(readCSV)  #to start of data

                #Get dates first; need dates to identify spacecraft orientation
                #for GOES-13+
                for row in readCSV:
                    date = datetime.datetime.strptime(row[0][0:18],
                                                    "%Y-%m-%d %H:%M:%S")
                    dates.append(date)

                if (experiment == "GOES-13" or experiment == "GOES-14"
                    or experiment == "GOES-15"):
                    west_detector = get_west_detector(filenames_orien[i], dates)


                #Go back and get fluxes
                count = 0
                csvfile.seek(0)
                for k in range(nhead):
                    next(readCSV)  #to start of data
                for row in readCSV:
                    for j in range(ncol):
                        flux = float(row[columns[j]])
                        #Account for orientation
                        if (experiment == "GOES-13" or experiment == "GOES-14"
                            or experiment == "GOES-15"):
                            if west_detector[count] == "B":
                                flux = float(row[columnsB[j]])
                            if west_detector[count] == "Flip":
                                flux = badval
                        if flux < 0:
                            flux = badval
                        fluxes[j][count] = flux
                    count = count + 1
            csvfile.close()


            #SECOND set of files for higher energy hepad
            nhead, nrow = find_goes_data_dimensions(filenames2[i])
            with open(datapath + '/' + filenames2[i]) as csvfile:
                readCSV = csv.reader(csvfile, delimiter=',')
                for k in range(nhead):
                    next(readCSV)  #to start of data

                count = 0
                for row in readCSV:
                    for j in range(nhcol):
                        flux = float(row[hepad_columns[j]])
                        if flux < 0:
                            flux = badval
                        fluxes[ncol+j][count] = flux
                    count = count + 1
            csvfile.close()

            #If reading in multiple files, then combine all data into one array
            if i==0:
                all_fluxes = fluxes
                all_dates = dates
            else:
                all_fluxes = np.concatenate((all_fluxes,fluxes),axis=1)
                all_dates = all_dates + dates

        return all_dates, all_fluxes



def read_in_user_files(filenames1):
    """Read in file containing flux time profile information that was
       specified by the user.
       The first column MUST contain the date in YYYY-MM-DD HH:MM:SS
       format. The remaining flux columns to be read in are specified by the
       user in the variable user_col at the very beginning of this program.
       The date column should always be considered column 0, even if you used
       whitespace as your delimeter. The code will consider the date format
       YYYY-MM-DD HH:MM:SS as one column even though it contains whitespace.
       Any number of header lines are allowed, but they must be indicated by #
       at the very beginning, including empty lines.
       Be sure to add the energy bins associated with your flux columns in the
       subroutine define_energy_bins under the "user" is statement.
    """
    print('Reading in user-specified files.')
    NFILES = len(filenames1)
    ncol = len(user_col) #include column for date
    for i in range(NFILES):
        print('Reading in ' + datapath + '/' + filenames1[i])
        with open(datapath + '/' + filenames1[i]) as csvfile:
            #Count header lines indicated by hash #
            nhead = 0
            for line in csvfile:
                line = line.lstrip()
                if line[0] == "#":
                    nhead = nhead + 1
                else:
                    break
            #number of lines containing data
            nrow = len(csvfile.readlines())+1

            #Define arrays that hold dates and fluxes
            dates = []
            fluxes = np.zeros(shape=(ncol,nrow))

            csvfile.seek(0) #back to beginning of file
            for k in range(nhead):
                csvfile.readline()  # Skip header rows.

            if user_delim == " " or user_delim == "":
                for j in range(len(user_col)):
                    #date takes two columns if separated by whitespace
                    #adjust the user input columns to account for this
                    user_col[j] = user_col[j] + 1

            count = 0
            for line in csvfile:
                if line == " " or line == "":
                    continue
                if user_delim == " " or user_delim == "":
                    row = line.split()
                    str_date = row[0][0:10] + ' ' + row[1][0:8]
                if user_delim != " " and user_delim != "":
                    row = line.split(user_delim)
                    str_date = row[0][0:18]

                date = datetime.datetime.strptime(str_date,
                                                "%Y-%m-%d %H:%M:%S")
                dates.append(date)
                for j in range(len(user_col)):
                   # print("Read in flux for column " + str(user_col[j]) + ': '
                   #     + str(date) + ' ' + row[user_col[j]])
                    flux = float(row[user_col[j]])
                    if flux < 0:
                        flux = badval
                    fluxes[j][count] = flux
                count = count + 1
        #If reading in multiple files, then combine all data into one array
        #SEPEM currently only has one file, but making generalized
        if i==0:
            all_fluxes = fluxes
            all_dates = dates
        else:
            all_fluxes = np.concatenate((all_fluxes,fluxes),axis=1)
            all_dates = all_dates + dates

    return all_dates, all_fluxes


def define_energy_bins(experiment,flux_type):
    """Define the energy bins for the selected spacecraft or data set.
       If the user inputs their own file, they must set the user_energy_bins
       variable at the top of the code.
    """
    #use corrected proton flux for GOES eps or epead; include hepad
    #-1 indicates infinity ([700, -1] means all particles above 700 MeV)
    if experiment == "SEPEM":
        energy_bins = [[5.00,7.23],[7.23,10.46],[10.46,15.12],[15.12,21.87],
                       [21.87,31.62],[31.62,45.73],[45.73,66.13],
                       [66.13,95.64],[95.64,138.3],[138.3,200.0],
                       [200.0,289.2]]

    if experiment == "ERNEf10":
        #The top 3 channels tend to saturate and show incorrect values during
        #high intensity SEP events. For this reason, only the >10 MeV
        #integral fluxes should be derived from ERNE data.
        #f10 format, from 2 December 1996
        energy_bins = [[10.0,13.0],[14.0,17.0],[17.0,22.0],[21.0,28.0],
                       [26.0,32.0],[32.0,40.0],[41.0,51.0],
                       [51.0,67.0],[54.0,79.0],[79.0,114.0],[111.0,140.]]

    if experiment == "ERNEf40":
        #f40 format, from 19 May 2000
        energy_bins = [[10.0,13.0],[14.0,17.0],[17.0, 22.0],[21.0,28.0],
                       [26.0,32.0],[32.0,40.0],[40.0,51.0],[51.0,67.0],
                       [64.0,80.0],[80.0,101.0],[101.0,131.0]]

    if (flux_type == "integral" and experiment != "SEPEM"):
         energy_bins = [[5.0,-1],[10.0,-1],[30.0,-1],
                        [50.0,-1],[60.0,-1],[100.0,-1],[700.0,-1]]

    if (experiment == "GOES-08" or experiment == "GOES-10" or
        experiment == "GOES-11" or experiment == "GOES-12"):
        if (flux_type == "differential"):
            #files named e.g. g08_eps_5m_yyyymmdd_yyyymmdd.csv
            energy_bins = [[4.0,9.0],[9.0,15.0],[15.0,44.0],
                           [40.0,80.0],[80.0,165.0],[165.0,500.0],
                           [350.0,420.0],[420.0,510.0],[510.0,700.0],
                           [700.0,-1]]

    if (experiment == "GOES-13" or experiment == "GOES-14" or
        experiment == "GOES-15"):
        if (flux_type == "differential"):
            energy_bins = [[4.2,8.7],[8.7,14.5],[15.0,40.0],
                            [38.0,82.0],[84.0,200.0],[110.0,900.0],
                            [330.0,420.0],[420.0,510.0],[510.0,700.0],
                            [700.0,-1]]

    if experiment == "user":
        #modify to match your energy bins or integral channels
        #use -1 in the second edge of the bin for integral channel (infinity)
        energy_bins = user_energy_bins

    return energy_bins


def extract_date_range(startdate,enddate,all_dates,all_fluxes):
    """Extract fluxes only for the dates in the range specified by the user."""
    print('Extracting fluxes for dates: ' + str(startdate) + ' to '
        + str(enddate))
    ndates = len(all_dates)
    nst = 0
    nend = 0
    for i in range(ndates):
        if all_dates[i] <= startdate:
            nst = i
        if all_dates[i] <= enddate:
            nend = i
    nend = min(nend+1, ndates)  #adjust to include nend in date range
    dates = all_dates[nst:nend]
    fluxes = all_fluxes[:,nst:nend]

    return dates, fluxes


def do_interpolation(i,dates,flux):
    """If bad fluxes (flux < 0) are found in the data, find the first prior
       data point and the first following data point that have good flux values.
       Perform linear interpolation in time:
            F(t) = F1 + (t - t1)*(F2 - F1)/(t2 - t1)
       This subroutine does the calculation for a single instance of bad data
       that corresponds to array index i.
    """
    ndates = len(dates)

    #If first point is bad point, use the next good point to fill gap
    if i == 0:
        for j in range(i,ndates-1):
            if flux[j] != badval:
                postflux = flux[j]
                postdate = dates[j]
                print('First point in array is bad. The first good value after '
                    'the gap is on '+ str(dates[j]) + ' with value '
                    + str(flux[j]))
                break
        preflux = postflux

    #If last point is bad point, use the first prior good point to fill gap
    if i == ndates - 1:
        for j in range(i,-1,-1):
            if flux[j] != badval:
                preflux = flux[j]
                predate = dates[j]
                print('Last point in the array is bad. The first good value '
                    'previous to gap is on '+ str(dates[j]) + ' with value '
                    + str(flux[j]))
                break
        postflux = preflux

    #Within the flux array
    if i != 0 and i != ndates-1:
        #search for first previous good value prior to the gap
        for j in range(i,-1,-1):
            if flux[j] != badval:
                preflux = flux[j]
                predate = dates[j]
                print('The first good value previous to gap is on '
                    + str(dates[j]) + ' with value ' + str(flux[j]))
                break
            if j == 0:
                sys.exit('There is a data gap at the beginning of the '
                        'selected time period. Program cannot estimate '
                        'flux in data gap.')

        #search for first previous good value after to the gap
        for j in range(i,ndates-1):
            if flux[j] != badval:
                postflux = flux[j]
                postdate = dates[j]
                print('The first good value after to gap is on '
                    + str(dates[j]) + ' with value ' + str(flux[j]))
                break
            if j == ndates-2 and flux[j] == badval:
                if flux[ndates-1] != badval:
                    postflux = flux[ndates-1]
                    postdate = dates[ndates-1]
                else:
                    postflux = preflux
                    postdate = predate
                    print(' Bad values continue to the end of the data set. '
                        'Using the first good value previous to gap on '
                        + str(postdate) + ' with value ' + str(postflux))

    if preflux == postflux:
        interp_flux = preflux
    if preflux != postflux:
        interp_flux = preflux + (dates[i] - predate).total_seconds()\
             *(postflux - preflux)/(postdate - predate).total_seconds()
    print('Filling gap at time ' + str(dates[i])
            + ' with interpolated flux ' + str(interp_flux))
    return interp_flux


def check_for_bad_data(dates,fluxes,energy_bins):
    """Search the data for bad values (flux < 0) and fill the missing data with
       an estimate flux found by performing a linear interpolation with time,
       using the good flux values immediately surrounding the data gap.
    """
    print('Checking for bad data values and filling with linear interpolation '
            + 'with time.')
    ndates = len(dates)
    nbins = len(energy_bins)

    for j in range(ndates):  #flux at each time
        for i in range(nbins):
            if fluxes[i,j] < 0: #bad data
                #estimate flux with interpolation in time
                print()
                print('There is a data gap for time ' + str(dates[j])
                        + ' and energy ' + str(energy_bins[i][0]) + ' - '
                        + str(energy_bins[i][1]) + ' MeV.'
                        + ' Filling in missing value with linear '
                        + 'interpolation in time.')
                interp_flux = do_interpolation(j,dates,fluxes[i,:])
                fluxes[i,j] = interp_flux
    print('Finished checking for bad data.')
    print()
    return fluxes


def from_differential_to_integral_flux(experiment, min_energy, energy_bins,
                fluxes):
    """If user selected differential fluxes, convert to integral fluxes to
       caluculate operational threshold crossings (>10 MeV protons exceed 10
       pfu, >100 MeV protons exceed 1 pfu).
       Assume that the measured fluxes correspond to the center of the energy
       bin and use power law interpolation to extrapolate integral fluxes
       above user input min_energy.
       The intent is to calculate >10 MeV and >100 MeV fluxes, but leaving
       flexibility for user to define the minimum energy for the integral flux.
       An integral flux will be provided for each timestamp (e.g. every 5 mins).
    """
    print('Converting differential flux to integral flux for >'
            + str(min_energy) + 'MeV.')
    nbins = len(energy_bins)
    nflux = len(fluxes[0])
    #Check requested min_energy inside data energy range
    if min_energy < energy_bins[0][0] or min_energy >= energy_bins[nbins-1][0]:
        sys.exit('The selected minimum energy ' + str(min_energy) + ' to create'
                ' integral fluxes is outside of the range of the data: '
                + str(energy_bins[0][0]) + ' - ' + str(energy_bins[nbins-1][0]))
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
        experiment == "GOES-15"):
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
                    sys.exit('Bad flux data value found for bin ' + str(i) + ','
                            + str(j) + '. This should not happen. '
                            + 'Did you call check_for_bad_data() first?')

                if fluxes[i,j] == 0 and fluxes[i+1,j] == 0: #add 0 flux
                    ninc = ninc + 1
                    continue

                F1 = fluxes[i,j]
                F2 = fluxes[i+1,j]
                if fluxes[i,j] == 0: #sub very small value for interpolation
                    F1 = 1e-15
                if fluxes[i+1,j] == 0:
                    F2 = 1e-15
                logF1 = np.log(F1)
                logF2 = np.log(F2)
                logE1 = np.log(bin_center[i])
                logE2 = np.log(bin_center[i+1])
                endE = bin_center[i+1]

                f = lambda x:exp(logF1
                            + (np.log(x)-logE1)*(logF2-logF1)/(logE2-logE1))
                startE = max(bin_center[i],min_energy)
                fint = scipy.integrate.quad(f,startE,endE)
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
            energy_thresholds, energy_bins):
    """Select or create the integral fluxes that correspond to the desired
       energy thresholds.
       If the user selected differential fluxes, then the
       differential fluxes will be converted to integral fluxes with the
       minimum energy defined by the set energy thresholds.
       If the user selected integral fluxes, then the channels corresponding to
       the desired energy thresholds will be identified.
    """
    nthresh = len(flux_thresholds)
    nenergy = len(energy_bins)

    #Calculate integral fluxes corresponding the the energy_thresholds
    integral_fluxes = []
    if flux_type == "differential":
        for i in range(nthresh):
            integral_flux = from_differential_to_integral_flux(experiment, \
                            energy_thresholds[i], energy_bins, fluxes)
            if i == 0:
                integral_fluxes = [np.array(integral_flux)]
            else:
                integral_fluxes = np.concatenate((integral_fluxes,
                                        [np.array(integral_flux)]))
    if flux_type == "integral":
        for i in range(nthresh):
            for j in range(nenergy):
                if energy_bins[j][0] == energy_thresholds[i]:
                    if len(integral_fluxes) == 0:
                        integral_fluxes = [np.array(fluxes[j,:])]
                    else:
                        integral_fluxes = np.concatenate((integral_fluxes,
                                                [np.array(fluxes[j,:])]))
        if len(integral_fluxes) == 0:
            sys.exit("There were no integral fluxes with the energy "
                    "thresholds specified by the user. Exiting.")
        if len(integral_fluxes) != len(energy_thresholds):
            sys.exit("Integral fluxes did not exist for some of the energy "
                    "thresholds specified by the user. Exiting.")

    return integral_fluxes


def integral_threshold_crossing(energy_threshold,flux_threshold,dates,fluxes):
    """Calculate the time that a threshold is crossed.
       Operational thresholds used by the NASA JSC Space Radiation Analysis
       Group to determine actions that should be taken during an SEP event are:
            >10 MeV proton flux exceeds 10 pfu (1/[cm^2 s sr])
            >100 MeV proton flux exceeds 1 pfu (1/[cm^2 s sr])
        An SEP event is considered to start if 3 consecutive points are
        above threshold. The start time is set to the first point that crossed
        threshold.
       If the user input differential flux, then the program has converted it
       to integral fluxes to calculate the integral threshold crossings.
       If the user selected fluxes that were already integral fluxes, then
       the correct channel was identified to determine these crossings.
       GOES and SEPEM data have a 5 minute time resolution, so initial
       crossings are accurate to 5 minutes.
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
    """
    #flux_threshold, e.g. 10 pfu (1/[cm^2 s sr])
    #energy_threshold, e.g. 10 MeV --> integral flux for >10 MeV protons
    #dates contain all of the datetimes for each data point
    #fluxes are a 1D array of integral fluxes (not multiple energy channels)
    print('Calculating threshold crossings and SEP event characteristics.')

    ndates = len(dates)
    end_threshold = 0.85*flux_threshold #used by SRAG operators

    threshold_crossed = False
    event_ended = False
    peak_flux = 0
    peak_time = 0
    crossing_time = 0 #define in case threshold not crossed
    rise_time = 0 #define in case threshold not crossed
    event_end_time = 0
    duration = 0
    end_counter = 0
    for i in range(ndates):
        if not threshold_crossed:
            if(fluxes[i] > flux_threshold):
                if i+2 < ndates:
                    if fluxes[i+1] > flux_threshold and fluxes[i+2] > \
                    flux_threshold:
                        crossing_time = dates[i]
                        threshold_crossed = True
        if threshold_crossed and not event_ended:
            if (fluxes[i] > end_threshold):
                end_counter = 0  #reset if go back above threshold
            if (fluxes[i] <= end_threshold): #flux drops below 85% threshold
                end_counter = end_counter + 1
                if end_counter == 3: #three consecutive points triggers end
                    event_ended = True
                    event_end_time = dates[i-2] #correct back two time steps

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

        rise_time = peak_time - crossing_time
        duration = event_end_time - crossing_time

    return crossing_time,peak_flux,peak_time,rise_time,event_end_time,duration


def calculate_fluence(dates, flux):
    """This subroutine sums up all of the flux in the 1D array "flux". The
       "dates" and "flux" arrays input here should reflect only the intensities
       between the SEP start and stop times, determined by the subroutine
       integral_threshold_crossing. "flux" should contain the intensity
       time series for a single energy bin or single integral channel (1D
       array). The subroutine does not differentiate between differential or
       integral flux. dates contains the 1D array of datetimes that
       correspond to the flux measurements.
       The extract_date_range subroutine is used prior to calling this one to
       make the dates and fluxes arrays covering only the SEP time period.
       The flux will be multiplied by time_resolution and summed for all of the
       data between the start and end times. Data gaps are taken into account.
       Fluence units will be 1/[MeV cm^2 sr] for GOES or SEPEM differential
       fluxes or 1/[cm^2 sr] for integral fluxes.
       For SRAG operational purposes, the event start and end is determined by
       the >100 MeV fluxes, as they are most pertinent to astronaut health
       inside of the space station.
    """
    ndates = len(dates)
    time_resolution = (dates[1] - dates[0]).total_seconds()
    sum_flux = 0
    for i in range(ndates):
        if flux[i] >= 0:  #0 flux ok for models
            sum_flux = sum_flux + flux[i]*time_resolution
        else:
            sys.exit('Bad flux data value found for bin ' + str(i) + ', '
                    + str(dates[i]) + '. This should not happen. '
                    + 'Did you call check_for_bad_data() first?')
    return sum_flux


def get_fluence_spectrum(experiment, flux_type, model_name, energy_threshold,
                flux_threshold, sep_dates, sep_fluxes, energy_bins, save_file):
    """Calculate the fluence spectrum for each of the energy channels in the
       user selected data set. If the user selected differential fluxes, then
       the fluence values correspond to each energy bin. If the user selected
       integral fluxes, then the fluence values correspond to each integral bin.
       Writes fluence values to file according to boolean save_file.
    """
    nenergy = len(energy_bins)
    fluence = np.zeros(shape=(nenergy))
    energies = np.zeros(shape=(nenergy))
    for i in range(nenergy):
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
        foutname = outpath + '/fluence_' + str(experiment) + '_' \
                    + str(flux_type) + '_' + 'gt' + str(energy_threshold) \
                    + '_' +str(year) + '_' + str(month) + '_' + str(day) +'.csv'
        if experiment == 'user' and model_name != '':
            foutname = outpath + '/fluence_' + model_name + '_' \
                        + str(flux_type) + '_' + 'gt' + str(energy_threshold) \
                        + '_' +str(year) + '_' + str(month) + '_' + str(day) \
                        +'.csv'
        fout = open(foutname,"w+")
        fout.write('#Event defined by >' + str(energy_threshold) + ', '
                    + str(flux_threshold) + ' pfu; start time '
                    + str(sep_dates[0]) + ', end time '
                    + str(sep_dates[len(sep_dates)-1]) + '\n')
        if flux_type == "differential":
            fout.write("#Elow,Emid,Ehigh,Fluence 1/[MeV cm^2 sr]\n")
        if flux_type == "integral":
            fout.write("#>Elow,Fluence 1/[cm^2 sr]\n")

        for i in range(nenergy):
            if flux_type == "differential":
                fout.write("{0},{1},{2},{3}\n".format(energy_bins[i][0],
                        energies[i], energy_bins[i][1], fluence[i]))
            if flux_type == "integral":
                fout.write("{0},{1}\n".format(energy_bins[i][0], fluence[i]))
        fout.close()

    return fluence, energies


def calculate_event_info(energy_thresholds,flux_thresholds,dates,
                integral_fluxes, detect_prev_event, two_peaks):
    """Uses the integral fluxes (either input or estimated from differential
       channels) and all the energy and flux thresholds set in the main program
       to calculate SEP event quantities.
            Threshold crossing time (onset)
            Peak Flux in date range specified by user
            Time of Peak Flux
            Rise Time (onset to peak flux)
            Event End Time (below 0.85*threshold for 3 data points, e.g 15 min)
            Duration (onset to end)
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
    """
    nthresh = len(flux_thresholds)
    crossing_time = []
    peak_flux = []
    peak_time = []
    rise_time = []
    event_end_time = []
    duration = []
    for i in range(nthresh):
        ct,pf,pt,rt,eet,dur = integral_threshold_crossing(energy_thresholds[i],
                        flux_thresholds[i],dates,integral_fluxes[i])
        if detect_prev_event and ct == dates[0]:
            print("Threshold may have been high due to previous event."
                "Recalculating event info for remaining time period in data "
                "set.")
            last_date = dates[len(dates)-1]
            tmp_dates, tmp_fluxes = extract_date_range(eet,last_date,dates,
                                integral_fluxes)
            ct,pf,pt,rt,eet,dur = integral_threshold_crossing(\
                            energy_thresholds[i],flux_thresholds[i],
                            tmp_dates,tmp_fluxes[i])

        if dur != 0 and two_peaks:
            if two_peaks and dur < timedelta(days=1):
                print("User specified that event has two peaks. Extending "
                    "event to second decrease below threshold.")
                last_date = dates[len(dates)-1]
                tmp_dates, tmp_fluxes = extract_date_range(eet,last_date,dates,
                                    integral_fluxes)
                ct2,pf2,pt2,rt2,eet2,dur2 = integral_threshold_crossing(\
                                energy_thresholds[i],flux_thresholds[i],
                                tmp_dates,tmp_fluxes[i])

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
        print(
               'Flux      Threshold    Time Crossed         Peak Flux'
                + '            Peak Time' + '            Rise Time'
                + '  End Time' + '            Duration'
             )
        print(
                '>' + str(energy_thresholds[i]) + ' MeV   '
                + str(flux_thresholds[i]) + ' pfu       '
                + str(crossing_time[i]) + '  '+ str(peak_flux[i]) + '   '
                + str(peak_time[i]) + '  ' + str(rise_time[i]) + '    '
                + str(event_end_time[i]) + '  ' + str(duration[i])
             )
    return crossing_time,peak_flux,peak_time,rise_time,event_end_time,duration


def report_threshold_fluences(experiment, flux_type, model_name,
                energy_thresholds, energy_bins, sep_dates, sep_fluxes):
    """Report fluences for specified thresholds, typically >10, >100 MeV.
       These values are interesting to use for comparison with literature and
       for quantifying event severity.
       If the user has input integral channels, then this information has also
       been saved in the output fluence file. If the user selected differential
       channels, then these values come from the estimated integral fluxes.
    """
    tmp_energy_bins = []
    nthresh = len(energy_thresholds)
    ndates = len(sep_dates)
    sep_integral_fluxes = np.zeros(shape=(nthresh,ndates))
    for i in range(nthresh):
        tmp_energy_bins.append([energy_thresholds[i],-1])
        if flux_type == "differential":
            #integral fluxes were estimated for only the threshold energies
            sep_integral_fluxes[i,:] = sep_fluxes[i,:]
        if flux_type == "integral":
            #Pull out the integral fluxes from the correct energy channel
            for j in range(len(energy_bins)):
                if energy_bins[j][0] == energy_thresholds[i]:
                    sep_integral_fluxes[i,:] = sep_fluxes[j,:]

    integral_fluence, integral_energies = get_fluence_spectrum(experiment,
                    "integral",model_name, 0, 0,sep_dates, sep_integral_fluxes,
                    tmp_energy_bins, False)

    for i in range(nthresh):
        integral_fluence[i] = integral_fluence[i]*4.*math.pi
        print("Event-integrated fluence for >" + str(integral_energies[i])
            + " MeV: " + str(integral_fluence[i]) + " 1/[cm^2]")

    return integral_fluence


def save_integral_fluxes_to_file(experiment, flux_type, model_name,
        energy_thresholds, crossing_time, dates, integral_fluxes):
    """
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
        sys.exit("No thresholds were crossed during this time period. "
                "Integral flux time profiles not written to file. Exiting.")

    #e.g. integral_fluxes_GOES-13_differential_2012_3_7.csv
    foutname = outpath + '/integral_fluxes_' + experiment + '_' \
                 + flux_type + '_' + str(year) + '_' + str(month) \
                 + '_' + str(day) + '.csv'
    if experiment == 'user' and model_name != '':
        foutname = outpath + '/integral_fluxes_' + model_name + '_' \
                     + flux_type + '_' + str(year) + '_' + str(month) \
                     + '_' + str(day) + '.csv'
    print('Writing integral flux time series to file --> ' + foutname)
    fout = open(foutname,"w+")
    fout.write('#Integral fluxes in units of 1/[cm^2 s sr] \n')
    fout.write('#Columns headers indicate low end of integral channels in MeV;'
                    ' e.g. >10 MeV \n')
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



def print_values_to_file(experiment, flux_type, model_name, energy_thresholds,
                flux_thresholds, crossing_time, peak_flux, peak_time, rise_time,
                event_end_time, duration, integral_fluences):
    """Write all calculated values to file for all thresholds. Event-integrated
       fluences for >10, >100 MeV (and user-defined threshold) will also be
       included. Writes out file with name e.g.
       output/sep_values_experiment_fluxtype_YYYY_M_D.csv
    """
    nthresh = len(energy_thresholds)
    year = 0
    month = 0
    day = 0
    for i in range(nthresh):
        if crossing_time[i] != 0 and year == 0:
            year = crossing_time[i].year
            month = crossing_time[i].month
            day = crossing_time[i].day
    if year == 0:
        sys.exit("No thresholds were crossed during this time period. "
                "Sep values not written to file. Exiting.")

    foutname = outpath + '/sep_values_' + experiment + '_' \
                + flux_type + '_' + str(year) + '_' + str(month) + '_' \
                + str(day) +'.csv'
    if experiment == 'user' and model_name != '':
        foutname = outpath + '/sep_values_' + model_name + '_' \
                    + flux_type + '_' + str(year) + '_' + str(month) + '_' \
                    + str(day) +'.csv'
    print('Writing SEP values to file --> ' + foutname)
    fout = open(foutname,"w+")
    #Write header
    if flux_type == "integral":
        fout.write('#Energy Threshold [MeV],Flux Threshold [pfu],Start Time,'
            + 'Peak Flux [pfu],Peak Time,Rise Time,End Time,Duration')
    if flux_type == "differential":
        fout.write('#Energy Threshold [MeV],Flux Threshold [pfu],'
            + 'Start Time,Peak Flux 1/[MeV cm^2 s sr],Peak Time,Rise Time,'
            'End Time,Duration')
    for i in range(nthresh):
        fout.write(',Fluence >' + str(energy_thresholds[i]) + ' MeV [cm^-2]')
    fout.write('\n')
    nthresh = len(energy_thresholds)

    for i in range(nthresh):
        if crossing_time[i] == 0: #no threshold crossed
            continue
        fout.write('>' + str(energy_thresholds[i]) + ',')
        fout.write(str(flux_thresholds[i]) + ',')
        fout.write(str(crossing_time[i]) + ',')
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
        for j in range(nthresh):
            fout.write(',' + str(integral_fluences[i,j]))
        fout.write('\n')

    fout.close()


######## MAIN PROGRAM #########
def run_all(str_startdate, str_enddate, experiment, flux_type, model_name,
        user_file, showplot, detect_prev_event, two_peaks, str_thresh):
    """"Runs all subroutines and gets all needed values. Takes the command line
        areguments as input. Written here to allow code to be imported into
        other python scripts.
        str_startdate, str_enddate, experiment, flux_type are strings.
        model_name is a string. If model is "user", set model_name to describe
        your model (e.g. MyModel), otherwise set to ''.
        user_file is a string. Defaul is ''. If user is selected for experiment,
        then name of flux file is specified in user_file.
        showplot, detect_prev_event, and two_peaks are booleans.
        Set str_thresh to be '100,1' for default value or modify to add your own
        threshold.
    """

    str_thresh = str_thresh.strip().split(",")
    input_threshold = [float(str_thresh[0]), float(str_thresh[1])]

    if len(str_startdate) == 10: #only YYYY-MM-DD
        str_startdate = str_startdate  + ' 00:00:00'
    if len(str_enddate) == 10: #only YYYY-MM-DD
        str_enddate = str_enddate  + ' 00:00:00'
    startdate = datetime.datetime.strptime(str_startdate, "%Y-%m-%d %H:%M:%S")
    enddate = datetime.datetime.strptime(str_enddate, "%Y-%m-%d %H:%M:%S")

    #CHECKS ON INPUTS
    if (str_startdate == "" or str_enddate == ""):
        sys.exit('You must enter a valid date range. Exiting.')
    if (experiment == "SEPEM" and flux_type == "integral"):
        sys.exit('The SEPEM (RSDv2) data set only provides differential fluxes.'
            ' Please change your FluxType to differential. Exiting.')
    sepem_end_date = datetime.datetime(2015,12,31,23,55,00)
    if(experiment == "SEPEM" and (startdate > sepem_end_date or
                   enddate > sepem_end_date)):
        sys.exit('The SEPEM (RSDv2) data set only extends to '
                  + str(sepem_end_date) +
            '. Please change your requested dates. Exiting.')

    user_fname[0] = user_file #input as argument, default is 'tmp.txt'

    #create data and output paths if don't exist
    check_paths()

    #Check and prepare the data
    filenames1, filenames2, filenames_orien = check_data(startdate,enddate,
                                                experiment,flux_type)
    #read in flux files
    if experiment != "user":
        all_dates, all_fluxes = read_in_files(experiment, flux_type, filenames1,
                                filenames2, filenames_orien)
    if experiment == "user":
        all_dates, all_fluxes = read_in_user_files(filenames1)
    #Define energy bins
    energy_bins = define_energy_bins(experiment,flux_type)
    #Extract the date range specified by the user
    dates, fluxes = extract_date_range(startdate,enddate,all_dates,all_fluxes)
    #Remove bad data points (negative fluxes) with linear interpolation in time
    fluxes = check_for_bad_data(dates,fluxes,energy_bins)

    #Define thresholds to use for start and end of event
    energy_thresholds = [10,100] #MeV; flux for particles of > this MeV
    flux_thresholds = [10,1] #pfu; exceed this level of intensity
    if input_threshold[0] != 100 or input_threshold[1] != 1:
        energy_thresholds.append(input_threshold[0])
        flux_thresholds.append(input_threshold[1])

    #Estimate or select integral fluxes corresponding the energy_thresholds
    integral_fluxes = extract_integral_fluxes(fluxes, experiment, flux_type,
                    flux_thresholds, energy_thresholds, energy_bins)

    #Calculate SEP event quantities for energy and flux threshold combinations
    #integral fluxes are used to define event start and stop
    #Calculate SEP values for all thresholds, however:
    #Operationally, >100 MeV fluxes above 1 pfu are used to define an SEP event
    #that has radiation dose signifigance for astronauts behind shielding.
    crossing_time,peak_flux,peak_time,rise_time,event_end_time,duration=\
        calculate_event_info(energy_thresholds,flux_thresholds, dates,
                    integral_fluxes, detect_prev_event, two_peaks)

    #Calculate event-integrated fluences for all thresholds
    nthresh = len(energy_thresholds)
    all_integral_fluences = np.zeros(shape=(nthresh,nthresh)) #fluences corresponding to >10, >100 MeV
    all_fluence = np.zeros(shape=(nthresh,len(energy_bins)))
    all_energies = np.zeros(shape=(nthresh,len(energy_bins)))
    for i in range(nthresh):
        #If no threshold was crossed during specified date range
        if crossing_time[i] == 0:
            print("The >" + str(energy_thresholds[i]) + " MeV threshold was "
                     "not crossed during the specified date range. No SEP "
                     "event. Continuing.")
            continue
        #Extract the original fluxes for the SEP start and stop times
        sep_dates, sep_fluxes = extract_date_range(crossing_time[i],
                             event_end_time[i],dates,fluxes)

        #Calculate fluence spectrum for the SEP event; either differential or
        #integral as specified by the user
        print('=====Calculating event fluence for event defined by >'
                + str(energy_thresholds[i]) + ' MeV, for '
                + str(crossing_time[i]) + ' to ' + str(event_end_time[i]))
        if crossing_time[i] == event_end_time[i]:
            print("Event start and end time the same (did you set "
            "--DetectPreviousEvent? May not work in this case). Exiting.")
            return
        fluence, energies = get_fluence_spectrum(experiment, flux_type,
                         model_name, energy_thresholds[i], flux_thresholds[i],
                         sep_dates, sep_fluxes, energy_bins, True)

        all_fluence[i] = fluence
        all_energies[i] = energies

        #Always calculate fluences for integral fluxes >10, >100 MeV and
        #any user input threshold fluences
        if flux_type == "integral":
            integral_fluence = report_threshold_fluences(experiment, flux_type,
                        model_name, energy_thresholds, energy_bins, sep_dates,
                        sep_fluxes)
        if flux_type == "differential":
            #Extract the estimated integral fluxes in the SEP event date range
            sep_integral_dates, sep_integral_fluxes = extract_date_range(\
                             crossing_time[i],event_end_time[i],
                             dates, integral_fluxes)
            integral_fluence = report_threshold_fluences(experiment, flux_type,
                         model_name, energy_thresholds, energy_bins,
                         sep_integral_dates, sep_integral_fluxes)

        all_integral_fluences[i] = integral_fluence


    #Save all calculated values for all threshold definitions to file
    print_values_to_file(experiment, flux_type, model_name, energy_thresholds,
                    flux_thresholds, crossing_time, peak_flux, peak_time,
                    rise_time, event_end_time, duration, all_integral_fluences)
    #Write >10, >100 and user input threshold integral fluxes to file
    save_integral_fluxes_to_file(experiment, flux_type, model_name,
            energy_thresholds, crossing_time, dates, integral_fluxes)

    #===============PLOTS==================
    #Plot selected results
    if showplot:
        #Event definition from integral fluxes
        if flux_type == "differential":
            print("Plotting estimated integral fluxes with threshold "
                   "crossings.")
        if flux_type == "integral":
            print("Plotting integral fluxes with threshold crossings.")
        #plot integral fluxes (either input or estimated)
        nthresh = len(flux_thresholds)
        plt.figure("Event_Def",figsize=(9,6))
        for i in range(nthresh):
            if crossing_time[i] == 0:
                continue
            plt.subplot(nthresh, 1, i+1)
            plt.plot_date(dates,integral_fluxes[i],'-')
            plt.axvline(crossing_time[i],color='black',linestyle=':')
            plt.axvline(event_end_time[i],color='black',linestyle=':')
            plt.axhline(flux_thresholds[i],color='red',linestyle=':')
            plt.plot_date(peak_time[i],peak_flux[i],'ro')
            if flux_type == "integral":
                plt.title(experiment + ' Flux for >'+str(energy_thresholds[i])
                           + ' MeV', x=0.75, y=0.8)
                if experiment == 'user' and model_name != '':
                    plt.title(model_name + ' Flux for >'
                        + str(energy_thresholds[i]) + ' MeV', x=0.75, y=0.8)

            if flux_type == "differential":
                plt.title(experiment + ' Estimated Flux for >'
                           + str(energy_thresholds[i]) + ' MeV', x=0.7, y=0.8)
                if experiment == 'user' and model_name != '':
                    plt.title(model_name + ' Estimated Flux for >'
                               + str(energy_thresholds[i]) + ' MeV',
                               x=0.7, y=0.8)
            plt.xlabel('Date')
            plt.ylabel('Integral Flux 1/[cm^2 s sr]')
            plt.yscale("log")
            ymin = max(1e-4, min(integral_fluxes[i]))
           # plt.ylim(ymin, peak_flux[i]+peak_flux[i]*.2)



        #All energy channels in specified date range with event start and stop
        print("Plotting fluxes in original energy bins. Any bad data points "
              "were interpolated. Lines indicate event start and stop for "
              "thresholds.")
        #Plot all channels of user specified data
        plt.figure("All_Bins",figsize=(9,4))
        ax = plt.subplot(111)
        nbins = len(energy_bins)
        for i in range(nbins):
            legend_label = ""
            if energy_bins[i][1] != -1:
                legend_label = str(energy_bins[i][0]) + '-' \
                               + str(energy_bins[i][1]) + ' MeV'
            else:
                legend_label = '>'+ str(energy_bins[i][0]) + ' MeV'
            ax.plot_date(dates,fluxes[i],'-',label=legend_label)
        colors = ['black','red','blue']
        for j in range(len(energy_thresholds)):
            if crossing_time[j] == 0:
                continue
            line_label = '>' + str(energy_thresholds[j]) + ' MeV, ' \
                        + str(flux_thresholds[j]) + ' pfu'
            ax.axvline(crossing_time[j],color=colors[j],linestyle=':',
                        label=line_label)
            ax.axvline(event_end_time[j],color=colors[j],linestyle=':')
        if flux_type == "integral":
            plt.ylabel('Integral Flux 1/[cm^2 s sr]')
            plt.title(experiment + " Integral Energy Bins with Threshold "
                            "Crossings")
            if experiment == 'user' and model_name != '':
                plt.title(model_name + " Integral Energy Bins with Threshold "
                                "Crossings")
        if flux_type == "differential":
            plt.ylabel('Flux 1/[MeV cm^2 s sr]')
            plt.title(experiment + " Differential Energy Bins with Threshold "
                            "Crossings")
            if experiment == 'user' and model_name != '':
                plt.title(model_name + " Differential Energy Bins with Threshold "
                                "Crossings")
        plt.xlabel('Date')
        plt.yscale("log")
        chartBox = ax.get_position()
        ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.85,
                         chartBox.height])
        ax.legend(loc='upper center', bbox_to_anchor=(1.17, 0.95))


        #Event-integrated fluence for energy channels
        print("Plotting event-integrated fluence spectrum.")
        #Plot fluence spectrum summed between SEP start and end dates
        plt.figure("Fluence",figsize=(6,5))
        ax = plt.subplot(111)
        markers = ['bo','P','D']
        for j in range(len(energy_thresholds)):
            if crossing_time[j] == 0:
                continue
            legend_label = '>' + str(energy_thresholds[j]) + ' MeV, ' \
                        + str(flux_thresholds[j]) + ' pfu'
            ax.plot(all_energies[j,:],all_fluence[j,:],markers[j],
                    color=colors[j],mfc='none',label=legend_label)
        plt.grid(which="both", axis="both")
        plt.title(experiment + ' Event-Integrated Fluences for All Event '
                    'Definitions')
        if experiment == 'user' and model_name != '':
            plt.title(model_name + ' Event-Integrated Fluences for All Event '
                        'Definitions')
        plt.xlabel('Energy [MeV]')
        if flux_type == "integral":
            plt.ylabel('Integral Fluxes 1/[cm^2 sr]')
        if flux_type == "differential":
            plt.ylabel('Flux 1/[MeV cm^2 sr]')
        plt.xscale("log")
        plt.yscale("log")
        ax.legend(loc='upper right')




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
            'SEPEM', 'user'],
            default='', help="Enter name of spacecraft or dataset")
    parser.add_argument("--FluxType", type=str, choices=['integral',
            'differential'], default='differential',
            help=("Do you want to use integral or differential fluxes?"
                 " (Default is differential)"))
    parser.add_argument("--ModelName", type=str, default='', help=("If you "
            "chose user for experiment, specify the name of the model or "
            "experiment that you are analyzing (no spaces)."))
    parser.add_argument("--UserFile", type=str, default='tmp.txt', help=("If you "
            "chose user for experiment, specify the filename containing the "
            "fluxes. Specify energy bins and delimeter in code at top. "
            "Default is tmp.txt."))
    parser.add_argument("--Threshold", type=str, default="100,1",
            help=("The energy and flux thresholds (written as 100,1 with no "
                    "spaces) which will be used to define the event. "
                    "e.g. 100,1 indicates >100 MeV fluxes crossing 1 pfu "
                    "(1/[cm^2 s sr]). Default = '100,1'"))
    parser.add_argument("--showplot",
            help="Flag to display plots", action="store_true")
    parser.add_argument("--DetectPreviousEvent",
            help=("Flag to indicate that the threshold is crossed at the "
                   "first point due to a previous event."), action="store_true")
    parser.add_argument("--TwoPeaks",
            help=("Flag to indicate that the event exceeds threshold (usually "
                    "for a couple of points), then drops below threshold "
                    "before increasing again to the true peak of the event."),
                    action="store_true")


    args = parser.parse_args()

    str_startdate = args.StartDate
    str_enddate = args.EndDate
    experiment = args.Experiment
    flux_type = args.FluxType
    model_name = args.ModelName
    user_file = args.UserFile
    str_thresh = args.Threshold
    showplot = args.showplot
    detect_prev_event = args.DetectPreviousEvent
    two_peaks = args.TwoPeaks

    run_all(str_startdate, str_enddate, experiment, flux_type, model_name,
        user_file, showplot, detect_prev_event, two_peaks, str_thresh)
    if showplot: plt.show()
