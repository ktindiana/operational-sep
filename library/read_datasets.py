from library import global_vars as vars
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
import numpy as np
import matplotlib.pyplot as plt
import sys

__version__ = "0.3"
__author__ = "Katie Whitman"
__maintainer__ = "Katie Whitman"
__email__ = "kathryn.whitman@nasa.gov"

#2021-01-06, Changes in 0.2: Added SEPEMv3 data set and made changes to
#   necessary subroutines to accomodate the new data set.
#2021-02-25, Changes 0.3: Changed GOES-13, 14, 15 S14 option to include
#   S14 corrections to channels P6 and P7. Had previously only
#   applied S14 to P2 - P5 for those experiments.


datapath = vars.datapath
outpath = vars.outpath
plotpath = vars.plotpath
badval = vars.badval #bad data points will be set to this value; must be negative
user_col = vars.user_col
user_delim = vars.user_delim
user_energy_bins = vars.user_energy_bins


def check_paths():
    """Check that the paths that hold the data and output exist. If not, create.
    """
    print('Checking that paths exist: ' + datapath + ' and ' + outpath)
    if not os.path.isdir(datapath):
        print('check_paths: Directory containing fluxes, ' + datapath +
        ', does not exist. Creating.')
        os.mkdir(datapath);
    if not os.path.isdir(datapath + '/GOES'):
        print('check_paths: Directory containing fluxes, ' + datapath +
        '/GOES, does not exist. Creating.')
        os.mkdir(datapath + '/GOES');
    if not os.path.isdir(datapath + '/SEPEM'):
        print('check_paths: Directory containing fluxes, ' + datapath +
        '/SEPEM, does not exist. Creating.')
        os.mkdir(datapath + '/SEPEM');
    if not os.path.isdir(datapath + '/SEPEMv3'):
        print('check_paths: Directory containing fluxes, ' + datapath +
        '/SEPEMv3, does not exist. Creating.')
        os.mkdir(datapath + '/SEPEMv3');
    if not os.path.isdir(datapath + '/EPHIN'):
        print('check_paths: Directory containing fluxes, ' + datapath +
        '/EPHIN, does not exist. Creating.')
        os.mkdir(datapath + '/EPHIN');
    if not os.path.isdir(outpath):
        print('check_paths: Directory to store output information, ' + outpath
            + ', does not exist. Creating.')
        os.mkdir(outpath);
    if not os.path.isdir(plotpath):
        print('check_paths: Directory to store plots does not exist. Creating.')
        os.mkdir(plotpath);


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
            date = datetime.datetime.strptime(row[0][0:19],
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


def check_sepem_data(startdate, enddate, experiment, flux_type):
    """Check if SEPEM data is present on the computer. Break into yearly
        files if needed. Return SEPEM filenames for analysis.
    """
    styear = startdate.year
    stmonth = startdate.month
    stday = startdate.day
    endyear = enddate.year
    endmonth = enddate.month
    endday = enddate.day

    filenames1 = []  #SEPEM, eps, or epead

    year = styear

    if experiment == 'SEPEM':
        basenm = 'SEPEM_H_reference'
        dir = experiment

    if experiment == 'SEPEMv3':
        basenm = 'SEPEM_RDS_v3_H'
        dir = experiment

    while (year <= endyear):
        fname = basenm + '_' + str(year) + '.csv'
        exists = os.path.isfile(datapath + '/' + dir + '/' + fname)
        if exists:
            filenames1.append(dir + '/' + fname)
            year = year + 1
        if not exists:
            full_exists = os.path.isfile(datapath + '/' + dir + '/' + \
                                 '/' + basenm + '.txt')
            if not full_exists:
                if experiment == 'SEPEM':
                    sys.exit("Please download and unzip the RSDv2 data set."
                        " You may download the file at"
                        " http://sepem.eu/help/SEPEM_RDS_v2-00.zip for full "
                        "fluxes or http://sepem.eu/help/SEPEM_RDS_v2-00.zip "
                        "for ESA background-subtracted fluxes.")
                if experiment == 'SEPEMv3':
                    sys.exit('Please contact DH Consultancy for the SEPEM '
                            'RDSv3 data set. Unzip and put SEPEM_RDS_V3_H.txt '
                            'in the data/SEPEMv3 folder.')
            if full_exists:
                #Break up SEPEM data set into yearly files
                print('The SEPEM (RSDv2 and RDSv3) is more tractable when '
                        'breaking into yearly data files. '
                        'Producing yearly files.')
                make_yearly_files(dir + '/' + basenm + '.txt')
                year = styear

    return filenames1


def check_goes_data(startdate, enddate, experiment, flux_type):
    """Check that GOES data is on your computer or download it from the NOAA
        website. Return the filenames associated with the correct GOES data.
    """
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
        exists1 = os.path.isfile(datapath + '/GOES/' + fname1)
        fname2 = prefix2 + date_suffix + '.csv'
        exists2 = os.path.isfile(datapath + '/GOES/' + fname2)
        if (experiment == "GOES-13" or experiment == "GOES-14"
            or experiment == "GOES-15"):
            fname_orien = prefix_orien + date_suffix + '_v1.0.0.csv'
            exists_orien = os.path.isfile(datapath + '/GOES/' + fname_orien)
            filenames_orien.append('GOES/' + fname_orien)

        filenames1.append('GOES/' + fname1)
        filenames2.append('GOES/' + fname2)

        if not exists1: #download file if not found on your computer
            url = ('https://satdat.ngdc.noaa.gov/sem/goes/data/avg/' +
                '%i/%02i/%s/csv/%s' % (year,month,satellite,fname1))
            print('Downloading GOES data: ' + url)
            try:
                urllib.request.urlopen(url)
                wget.download(url, datapath + '/GOES/' + fname1)
            except urllib.request.HTTPError:
                sys.exit("Cannot access file at " + url +
                ". Please check that selected spacecraft covers date range.")


        if not exists2: #download file if not found on your computer
            url = ('https://satdat.ngdc.noaa.gov/sem/goes/data/avg/' +
               '%i/%02i/%s/csv/%s' % (year,month,satellite,fname2))
            print('Downloading GOES data: ' + url)
            try:
                urllib.request.urlopen(url)
                wget.download(url, datapath + '/GOES/' + fname2)
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
                    wget.download(url, datapath + '/GOES/' + fname_orien)
                except urllib.request.HTTPError:
                    sys.exit("Cannot access orientation file at " + url +
                   ". Please check that selected spacecraft covers date range.")

    return filenames1, filenames2, filenames_orien


def check_ephin_data(startdate, enddate, experiment, flux_type):
    """Check for SOHO/COSTEP/EPHIN data on your computer. If not there,
        download from http://ulysses.physik.uni-kiel.de/costep/level3/l3i/
        30 minute data will be downloaded. Intensities are in units of
        (cm^2 s sr mev/nuc)^-1
        First available date is 1995 12 8 (DOY = 342).
        The files are available in daily or yearly format.
    """
    styear = startdate.year
    stmonth = startdate.month
    stday = startdate.day
    endyear = enddate.year
    endmonth = enddate.month
    endday = enddate.day

    #Array of filenames that contain the data requested by the User
    filenames1 = []  #SEPEM, EPHIN, eps, or epead

    Nyr = endyear - styear + 1
    for year in range(styear, endyear+1):
        fname = str(year) + '.l3i'
        filenames1.append('EPHIN/' + fname)

        exists = os.path.isfile(datapath + '/EPHIN/' + fname)
        if not exists: #download file if not found on your computer
            url = ('http://ulysses.physik.uni-kiel.de/costep/level3/l3i/30min/%s'
                    % (fname))
            print('Downloading EPHIN data: ' + url)
            try:
                urllib.request.urlopen(url)
                wget.download(url, datapath + '/EPHIN/' + fname)
            except urllib.request.HTTPError:
                sys.exit("Cannot access EPHIN file at " + url +
               ". Please check that selected spacecraft covers date range.")

    return filenames1


def check_ephin_release_data(startdate, enddate, experiment, flux_type):
    """Check for SOHO/COSTEP/EPHIN data on your computer provided by the
        HESPERIA collaboration on the website
        https://www.hesperia.astro.noa.gr/index.php/results/
        real-time-prediction-tools/data-retrieval-tool
        (one long URL). Monthly files of this data set are provided in the
        public git containing these codes. Otherwise, users may go to the URL
        above, download the data sets, and use the file naming convention
        adopted here:
            data/EPHIN_REleASE/HESPERIA_SOHO_PROTON_YYYY.txt
        Intensities are in units of (cm^2 s sr mev/nuc)^-1
        First available date is 1995 12 8 (DOY = 342).
        The files are saved in yearly format.
    """
    styear = startdate.year
    stmonth = startdate.month
    stday = startdate.day
    endyear = enddate.year
    endmonth = enddate.month
    endday = enddate.day

    #Array of filenames that contain the data requested by the User
    filenames1 = []  #SEPEM, EPHIN, eps, or epead

    Nyr = endyear - styear + 1
    for year in range(styear, endyear+1):
        fname = 'HESPERIA_SOHO_PROTON_' + str(year) + '.txt'
        filenames1.append('EPHIN_REleASE/' + fname)

        exists = os.path.isfile(datapath + '/EPHIN_REleASE/' + fname)
        if not exists: #download file if not found on your computer
                sys.exit("Cannot access EPHIN file " + fname +
               ". Please check that selected spacecraft covers date range.")

    return filenames1



def check_data(startdate, enddate, experiment, flux_type, user_file):
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

    user_fname = [user_file]

    #Array of filenames that contain the data requested by the User
    filenames1 = []  #SEPEM, eps, or epead
    filenames2 = []  #hepad
    filenames_orien = []  #orientation flag for GOES-13+

    #If user wants to use own input file (filename defined as input)
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
    if (experiment == "SEPEM" or experiment == "SEPEMv3"):
        filenames1 = check_sepem_data(startdate, enddate, experiment, flux_type)
        return filenames1, filenames2, filenames_orien

    if experiment[0:4] == "GOES":
        filenames1, filenames2, filenames_orien = check_goes_data(startdate, \
                                        enddate, experiment, flux_type)
        return filenames1, filenames2, filenames_orien

    if experiment == "EPHIN":
        filenames1 = check_ephin_data(startdate, enddate, experiment, flux_type)
        return filenames1, filenames2, filenames_orien

    if experiment == "EPHIN_REleASE":
        filenames1 = check_ephin_release_data(startdate, enddate, experiment,
                                        flux_type)
        return filenames1, filenames2, filenames_orien

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
            date = datetime.datetime.strptime(row[0][0:19],
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


def read_in_sepem(experiment, flux_type, filenames1):
    """Read in SEPEM data files from the computer."""
    NFILES = len(filenames1)
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
                date = datetime.datetime.strptime(row[0][0:19],
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


def read_in_goes(experiment, flux_type, filenames1, filenames2,
                filenames_orien, options):
    """Read in GOES data from your computer. User may specify option to choose
        corrected or uncorrected GOES fluxes.
    """
    NFILES = len(filenames1)
    west_detector = [] #place holder, will be filled if needed

    if (experiment == "GOES-08" or experiment == "GOES-10" or
        experiment == "GOES-11" or experiment == "GOES-12"):
        if flux_type == "differential":
            if "corrected" in options or "uncorrected" not in options:
                #CORRECTED CHANNELS; DEFAULT
                columns = [18,19,20,21,22,23] #eps
                hepad_columns = [1,2,3,4]
            if "uncorrected" in options:
                #UNCORRECTED channels
                columns = [5,6,7,8,9,10] #eps
                hepad_columns = [1,2,3,4]
        if flux_type == "integral":
            #ONLY CORRECTED CHANNELS AVAILABLE
            columns = [25,26,27,28,29,30] #eps
            hepad_columns = [4]
    if (experiment == "GOES-13" or experiment == "GOES-14" or
        experiment == "GOES-15"):
        if flux_type == "differential":
            if "corrected" in options or "uncorrected" not in options:
                #CORRECTED CHANNELS
                columns = [16,24,32,40,48,56] #epead, default A detector
                columnsB = [12,20,28,36,44,52] #epead, B detector
                hepad_columns = [9,12,15,18]
            if "uncorrected" in options:
                #UNCORRECTED CHANNELS
                columns = [15,23,31,39,47,55] #epead, default A detector
                columnsB = [11,19,27,35,43,51] #epead, B detector
                hepad_columns = [9,12,15,18]
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
                date = datetime.datetime.strptime(row[0][0:19],
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

    return all_dates, all_fluxes, west_detector


def read_in_ephin(experiment, flux_type, filenames1):
    """Read in EPHIN files from your computer."""
    NFILES = len(filenames1)

    datecols = [0,1,2,4,5] #yr, mth, dy, hr, min
    fluxcols = [8,9,10,11]
    ncol= len(fluxcols)

    for i in range(NFILES):
        print('Reading in file ' + datapath + '/' + filenames1[i])
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

            count = 0
            for line in csvfile:
                if line == '': continue
                if line[0] == "#": continue

                row = line.split()
                yr = int(row[datecols[0]])
                mnth = int(row[datecols[1]])
                dy = int(row[datecols[2]])
                hr = int(row[datecols[3]])
                min = int(row[datecols[4]])
                date = datetime.datetime(yr,mnth,dy,hr,min,0,0)
                dates.append(date)
                for j in range(ncol):
                    flux = float(row[fluxcols[j]])
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


def read_in_ephin_release(experiment, flux_type, filenames1):
    """Read in EPHIN files from your computer."""
    NFILES = len(filenames1)

    datecols = [0] #yr, mth, dy, hr, min
    fluxcols = [1,2,3,4]
    ncol= len(fluxcols)

    for i in range(NFILES):
        print('Reading in file ' + datapath + '/' + filenames1[i])
        with open(datapath + '/' + filenames1[i]) as csvfile:
            #Count header lines indicated by hash #
            nhead = 0
            for line in csvfile:
                line = line.lstrip()
                if line == '':
                    nhead = nhead + 1
                elif line[0] == "#":
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

            count = 0
            for line in csvfile:
                if line == '': continue
                if line[0] == "#": continue

                row = line.split(';')
                date = datetime.datetime.strptime(row[0][0:19],
                                                "%Y-%m-%d %H:%M:%S")
                dates.append(date)
                for j in range(ncol):
                    flux = float(row[fluxcols[j]])
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

    print(all_fluxes)
    return all_dates, all_fluxes


def read_in_files(experiment, flux_type, filenames1, filenames2,
                filenames_orien, options):
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
    all_dates = []
    all_fluxes = []
    west_detector = []

    if experiment == "SEPEM" or experiment == "SEPEMv3":
        all_dates, all_fluxes = read_in_sepem(experiment, flux_type, filenames1)
        return all_dates, all_fluxes, west_detector

    #All GOES data
    if experiment[0:4] == "GOES":
        all_dates, all_fluxes, west_detector = read_in_goes(experiment, \
                    flux_type, filenames1, filenames2, filenames_orien, options)
        return all_dates, all_fluxes, west_detector

    if experiment == "EPHIN":
        all_dates, all_fluxes = read_in_ephin(experiment, flux_type, filenames1)
        return all_dates, all_fluxes, west_detector

    if experiment == "EPHIN_REleASE":
        all_dates, all_fluxes = read_in_ephin_release(experiment, flux_type,
                    filenames1)
        return all_dates, all_fluxes, west_detector

    return all_dates, all_fluxes, west_detector


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
                if line == '' or line == '\n':
                    nhead = nhead + 1
                elif line[0] == "#" or line[0] == '\"':
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
                    str_date = row[0][0:19]

                date = datetime.datetime.strptime(str_date,
                                                "%Y-%m-%d %H:%M:%S")
                dates.append(date)
                for j in range(len(user_col)):
                   # print("Read in flux for column " + str(user_col[j]) + ': '
                   #     + str(date) + ' ' + row[user_col[j]])
                    if row[user_col[j]] == 'n/a': #REleASE
                        flux = None
                    else:
                        flux = float(row[user_col[j]])
                        if flux < 0:
                            flux = badval
                    fluxes[j][count] = flux
                count = count + 1

        #Remove dates that have None values (because they were n/a in REleASE)
        for k in range(len(dates)-1,-1,-1):
            if None in fluxes[:,k] or np.isnan(np.sum(fluxes[:,k])):
                del dates[k]
                fluxes = np.delete(fluxes, k, 1)

        #If reading in multiple files, then combine all data into one array
        #SEPEM currently only has one file, but making generalized
        if i==0:
            all_fluxes = fluxes
            all_dates = dates
        else:
            all_fluxes = np.concatenate((all_fluxes,fluxes),axis=1)
            all_dates = all_dates + dates

    return all_dates, all_fluxes


def extract_date_range(startdate,enddate,all_dates,all_fluxes):
    """Extract fluxes only for the dates in the range specified by the user."""
    #print('Extracting fluxes for dates: ' + str(startdate) + ' to '
    #    + str(enddate))
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
            if flux[j] != badval and flux[j] != None:
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
            if flux[j] != badval and flux[j] != None:
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
            if flux[j] != badval and flux[j] != None:
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
            if flux[j] != badval and flux[j] != None:
                postflux = flux[j]
                postdate = dates[j]
                print('The first good value after to gap is on '
                    + str(dates[j]) + ' with value ' + str(flux[j]))
                break
            if j == ndates-2 and flux[j] == badval:
                if flux[ndates-1] != badval and flux[ndates-1] != None:
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


def check_for_bad_data(dates,fluxes,energy_bins,dointerp=True):
    """Search the data for bad values (flux < 0) and fill the missing data with
       an estimate flux found by performing a linear interpolation with time,
       using the good flux values immediately surrounding the data gap.
    """
    if dointerp:
        print('Checking for bad data values and filling with linear '
              'interpolation with time.')
    else:
        print('Checking for bad data values and filling with None values. ')

    ndates = len(dates)
    nbins = len(energy_bins)

    for j in range(ndates):  #flux at each time
        for i in range(nbins):
            if fluxes[i,j] == None: #bad data
                #estimate flux with interpolation in time
                if dointerp:
                    print()
                    print('There is a data gap for time ' + str(dates[j])
                            + ' and energy ' + str(energy_bins[i][0]) + ' - '
                            + str(energy_bins[i][1]) + ' MeV.'
                            + ' Filling in missing value with linear '
                            + 'interpolation in time.')
                    interp_flux = do_interpolation(j,dates,fluxes[i,:])
                    fluxes[i,j] = interp_flux
                else:
                    print()
                    print('There is a data gap for time ' + str(dates[j])
                            + ' and energy ' + str(energy_bins[i][0]) + ' - '
                            + str(energy_bins[i][1]) + ' MeV.'
                            + ' Filling in missing value with None ')
                    fluxes[i,j] = None

            elif fluxes[i,j] < 0:
                #estimate flux with interpolation in time
                if dointerp:
                    print()
                    print('There is a data gap for time ' + str(dates[j])
                            + ' and energy ' + str(energy_bins[i][0]) + ' - '
                            + str(energy_bins[i][1]) + ' MeV.'
                            + ' Filling in missing value with linear '
                            + 'interpolation in time.')
                    interp_flux = do_interpolation(j,dates,fluxes[i,:])
                    fluxes[i,j] = interp_flux
                else:
                    print()
                    print('There is a data gap for time ' + str(dates[j])
                            + ' and energy ' + str(energy_bins[i][0]) + ' - '
                            + str(energy_bins[i][1]) + ' MeV.'
                            + ' Filling in missing value with None ')
                    fluxes[i,j] = None



    print('Finished checking for bad data.')
    print()
    return fluxes



def define_energy_bins(experiment,flux_type,west_detector,options):
    """Define the energy bins for the selected spacecraft or data set.
       If the user inputs their own file, they must set the user_energy_bins
       variable in library/global_vars.py.
       User may select options to apply Sandberg et al. (2014) effective
       energies for GOES EPS by specifying "S14" and/or apply Bruno (2017)
       effective energies for GOES-13 or -15 P6, P7 and HEPAD by specifying
       "Bruno2017"
    """
    #use corrected proton flux for GOES eps or epead; include hepad
    #-1 indicates infinity ([700, -1] means all particles above 700 MeV)
    if experiment == "SEPEM":
        energy_bins = [[5.00,7.23],[7.23,10.46],[10.46,15.12],[15.12,21.87],
                       [21.87,31.62],[31.62,45.73],[45.73,66.13],
                       [66.13,95.64],[95.64,138.3],[138.3,200.0],
                       [200.0,289.2]]

    if experiment == "SEPEMv3":
        energy_bins = [[5.00,7.23],[7.23,10.46],[10.46,15.12],[15.12,21.87],
                       [21.87,31.62],[31.62,45.73],[45.73,66.13],
                       [66.13,95.64],[95.64,138.3],[138.3,200.0],
                       [200.0,289.2],[289.2,418.3],[418.3,604.9],
                       [604.9,874.7]]

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

    if experiment == "EPHIN":
        #http://ulysses.physik.uni-kiel.de/costep/level3/l3i/
        #DOCUMENTATION-COSTEP-EPHIN-L3-20181002.pdf
        energy_bins = [[4.3,7.8],[7.8,25],[25,40.9],[40.9,53]]

    if experiment == "EPHIN_REleASE":
        #This data may be downloaded by hand through a web interface at:
        #https://www.hesperia.astro.noa.gr/index.php/results/real-time-prediction-tools/data-retrieval-tool
        energy_bins = [[4.0,9.0],[9.0,15.8],[15.8,39.6],[20.0,35.5]]

    if (flux_type == "integral" and experiment[0:4] == "GOES"):
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
            if "S14" in options:
                if experiment == "GOES-08":
                    energy_bins = [[4.0,7.9],[7.4,15.0],[13.3,21.3],
                               [37.0,53.6],[91.5,113.0],[119.0,179.0],
                               [350.0,420.0],[420.0,510.0],[510.0,700.0],
                               [700.0,-1]]
                if experiment == "GOES-11":
                    energy_bins = [[5.0,7.9],[9.4,15.9],[16.7,23.2],
                               [32.5,56.4],[89.8,114.0],[120.0,186.0],
                               [350.0,420.0],[420.0,510.0],[510.0,700.0],
                               [700.0,-1]]

    if (experiment == "GOES-13" or experiment == "GOES-14" or
        experiment == "GOES-15"):
        if (flux_type == "differential"):
            energy_bins = [[4.2,8.7],[8.7,14.5],[15.0,40.0],
                            [38.0,82.0],[84.0,200.0],[110.0,900.0],
                            [330.0,420.0],[420.0,510.0],[510.0,700.0],
                            [700.0,-1]]
            if "S14" in options:
                #S14 is not specifically calibrated to these GOES, but
                #apply the GOES-11 EPS energy bins for P2 - P7
                energy_bins = [[5.0,7.9],[9.4,15.9],[16.7,23.2],
                           [32.5,56.4],[89.8,114.0],[120.0,186.0],
                           [330.0,420.0],[420.0,510.0],[510.0,700.0],
                           [700.0,-1]]
            if "Bruno2017" in options:
                #EPEAD CHANNELS P6 and P7
                if "uncorrected" in options:
                    if west_detector.count("A") >= west_detector.count("B"):
                        #A detector bins
                        if experiment == "GOES-13":
                            energy_bins[4] = [93.3,129.0]
                            energy_bins[5] = [145.2,203.9]
                        if experiment == "GOES-15":
                            energy_bins[4] = [97.5,134.8]
                            energy_bins[5] = [142.2,199.0]
                    if west_detector.count("B") > west_detector.count("A"):
                        #B detector bins
                        if experiment == "GOES-13":
                            energy_bins[4] = [92.3,127.5]
                            energy_bins[5] = [143.5,200.5]
                        if experiment == "GOES-15":
                            energy_bins[4] = [95.9,132.3]
                            energy_bins[5] = [144.6,202.3]
                if "corrected" in options or "uncorrected" not in options: #Z89 applied
                    if west_detector.count("A") >= west_detector.count("B"):
                        #A detector bins
                        if experiment == "GOES-13":
                            energy_bins[4] = [93.3,129.1]
                            energy_bins[5] = [146.7,205.6]
                        if experiment == "GOES-15":
                            energy_bins[4] = [97.9,135.3]
                            energy_bins[5] = [145.0,202.3]
                    if west_detector.count("B") > west_detector.count("A"):
                        #B detector bins
                        if experiment == "GOES-13":
                            energy_bins[4] = [92.4,127.8]
                            energy_bins[5] = [144.6,202.5]
                        if experiment == "GOES-15":
                            energy_bins[4] = [96.1,132.8]
                            energy_bins[5] = [147.3,205.7]
                if experiment == "GOES-13": #Should use with bg-subtracted flux
                    energy_bins[6] = [273.9,387.5] #P8
                    energy_bins[7] = [330.0,458.0] #P9
                    energy_bins[8] = [418.7,566.0] #P10
                    energy_bins[9] = [852.6,1081.2] #P11
                if experiment == "GOES-15":
                    energy_bins[6] = [240.4,335.6] #P8
                    energy_bins[7] = [325.3,464.6] #P9
                    energy_bins[8] = [420.4,573.1] #P10
                    energy_bins[9] = [878.6,1230.0] #P11


    if experiment == "user":
        #modify to match your energy bins or integral channels
        #use -1 in the second edge of the bin for integral channel (infinity)
        energy_bins = user_energy_bins

    return energy_bins
