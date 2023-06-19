import re
import calendar
import datetime
from datetime import timedelta
import os
import wget
from calendar import monthrange
import urllib.request
import csv
import numpy as np
import sys
import math
import netCDF4

__version__ = "1.5"
__author__ = "Katie Whitman"
__maintainer__ = "Katie Whitman"
__email__ = "kathryn.whitman@nasa.gov"


#2021-01-06, Changes in 0.2: Added SEPEMv3 data set and made changes to
#   necessary subroutines to accomodate the new data set.
#2021-02-25, Changes 0.3: Changed GOES-13, 14, 15 S14 option to include
#   S14 corrections to channels P6 and P7. Had previously only
#   applied S14 to P2 - P5 for those experiments.
#2021-09-24, Changes in 0.4: added a global_var called time_shift
#   which allows users to shift the times in user-input files by
#   time_shift number of hours. Changed in read_in_user_files and
#   added convert_decimal_hour.
#2021-11-16, Changes in 0.5: added support for GOES-16 and GOES-17
#   differential fluxes.
#2022-02-11, Changes in 0.6: Added support for GOES-R (16&17) primary
#   integral fluxes served by CCMC. These are the real time fluxes
#   from NOAA archived on the CCMC website. These are not the
#   official NOAA L2 integral fluxes. Those are not yet
#   available.
#2022-02-18, Changes in 0.7: Added checking for data/GOES-R
#   directory and will make if not present.
#2022-03-23, Changes in 0.8: Added ability to download and read GOES-16
#   SEP event file on NOAA's website in check_goesR_data. Modified
#   read_in_goesR data since the special file contains 30 days
#   of data with slightly different variable names.
#2022-05-20. Changes in 0.9: Changed SOHO/EPHIN L3 data from 30 minute
#   to 10 min data.
#2022-06-16, changes in 1.0: Added Shaowen Hu's recalibrated GOES
#   data set as a native data set in the code (SRAG1.2)
#2022-08-04, changes in 1.1: in extract_date_range, abjusted
#   the trimming so that the selected time range starts either
#   on or one point AFTER the specified start. Previously,
#   the point right before the specified start was included.
#2022-09-19, changes in 1.2: GOES-14 and GOES-15 hepad files from
#   2019-09-01 forward are missing a column. Added code to
#   read_in_goes() to change the expected columns for later dates.
#2022-11-20, changes in 1.3: Added STEREO-A and B to native data sets.
#2023-02-09, changes in 1.4: NOAA SWPC moved the location of the historical
#   GOES-15 and previous data. Updated the url in check_goes_data().
#   Updated check_goesR_data() to account for two different version
#   numbers possible in the differential files. Updated read_in_goesR()
#   to account for the different keys used to extract the flux
#   values in the different versions.
#2023-04-20: SIMPLIFIED FOR ONLY GOES-R
#2023-06-19, changes in 1.5: NOAA added a v3-0-1 format for files
#   starting in April 2023. Rewrote check_goesR to be more versatile.
#   Added checking that include v3-0-1 in read_in_goesR.

datapath = "data"
badval = -1 #bad data points will be set to this value; must be negative

########## READ IN OWN FILE ######
#IF READING IN YOUR OWN FILE. SPECIFY COLUMNS CONTAINING FLUX,
#DELIMITER SEPARATING COLUMNS, AND ENERGY BINS ASSOCIATED WITH FLUX
#COLUMNS. FIRST COLUMN IN THE FILE MUST CONTAIN DATE IN FORMAT:
#YYYY-MM-DD HH:MM:SS
#Flux columns must start at column 1 or more.
#Even if the delimeter is white space, the code will treat the
#ISO date as the 0th column.
#user_col = arr.array('i',[1,2,3,4,5,6,7,8])
#user_delim = ""
#user_energy_bins = [[750,-1],[500,-1],[300,-1],[100,-1],\
#                    [60,-1],[50,-1],[30,-1],[10,-1]]

def about_read_datasets():
    """ About read_datasets.py
        
        Subroutines that are required to read in the data sets
        native to this code or the user-specified data set.
        
        Reads GOES-08 to GOES-15, GOES-R, SOHO/EPHIN Level 3 data,
        SOHO/EPHIN data from the REleASE website, SEPEM RDSv2
        and SEPEM RDSv3 (if the user downloads and unzips the
        files into the data directory).
        
        When possible, data is pulled from online databases and
        saved on your computer in the directory specified by
        datapath in global_vars.py. The default path is "data".
        
        Data files will be stored in subdirectories named as:
        
        * data/GOES-R/
 
        
    """

def check_paths():
    """Check that the paths that hold the data and output exist. If not, create.
    """
    print('Checking that paths exist: ' + datapath)
    if not os.path.isdir(datapath):
        print('check_paths: Directory containing fluxes, ' + datapath +
        ', does not exist. Creating.')
        os.mkdir(datapath);
    if not os.path.isdir(datapath + '/GOES-R'):
        print('check_paths: Directory containing fluxes, ' + datapath +
        '/GOES-R, does not exist. Creating.')
        os.mkdir(datapath + '/GOES-R');
 


def check_goesR_data(startdate, enddate, experiment, flux_type):
    """Check that GOES-R data is on your computer or download it from the NOAA
        website. Return the filenames associated with the correct GOES data.
        GOES-R files are saved daily in cdf format.
        
        GOES-16 & 17 differential fluxes are the official science-grade product
        provided by NOAA SWPC.
        
        INPUTS:
        
        :startdate: (datetime) start of time period specified by user
        :enddate: (datetime) end of time period entered by user
        :experiment: (string) name of native experiment or "user"
        :flux_type: (string) "integral" or "differential"
        
        OUTPUTS:
        
        :filenames1: (string array) the files containing the GOES
            EPS or EPEAD data that span the desired time range
            (monthly files)
        :filenames2: (string array) the files containing the GOES
            HEPAD data that span the desired time range
        :filenames_orien: (string array) the files
            that indicate the orientation of the GOES EPS or
            EPEAD detector (so can choose westward facing detector)
        
    """
    if flux_type == "integral":
        sys.exit("check_goesR_data: This subroutine is only valid for GOES-R "
                "differential fluxes. Please set the flux_type to differential "
                "and try again.")
                
    styear = startdate.year
    stmonth = startdate.month
    stday = startdate.day
    endyear = enddate.year
    endmonth = enddate.month
    endday = enddate.day

    #Array of filenames that contain the data requested by the User
    filenames1 = []  #GOES-R

    
    #SPECIAL FILE FOR 2017-09-10 SEP EVENTS
    if experiment == "GOES-16" and styear == 2017:
        fname1 = 'se_sgps-l2-avg5m_g16_s20172440000000_e20172732355000_v2_0_0.nc'
        filenames1.append('GOES-R/' + fname1)
        exists1 = os.path.isfile(datapath + '/GOES-R/' + fname1)
        if not exists1:
            url=('https://www.ngdc.noaa.gov/stp/space-weather/satellite-data/satellite-systems/goesr/solar_proton_events/sgps_sep2017_event_data/%s' % (fname1))
            try:
                urllib.request.urlopen(url)
                wget.download(url, datapath + '/GOES-R/' + fname1)
            except urllib.request.HTTPError:
                sys.exit("Cannot access SEP event file at " + url +
               ". Please check that the url is still active.")
        
        return filenames1
    
    
    #GOES-R data is stored in daily data files
    td = enddate - startdate
    NFILES = td.days #number of data files to download
    if td.seconds > 0: NFILES = NFILES + 1

    if experiment == "GOES-16":
        prefix = 'sci_sgps-l2-avg5m_g16_'
        satellite = 'goes16'

    if experiment == "GOES-17":
        prefix = 'sci_sgps-l2-avg5m_g17_'
        satellite = 'goes17'


    #for every day that data is required, check if file is present or
    #needs to be downloaded.
    for i in range(NFILES):
        date = startdate + datetime.timedelta(days=i)
        year = date.year
        month = date.month
        day = date.day
        date_suffix = 'd%i%02i%02i' % (year,month,day)
 
        #GOES-R differential data has three possible version numbers
        file_ext = ['_v1-0-1.nc', '_v2-0-0.nc', '_v3-0-0.nc', '_v3-0-1.nc']
        
        foundfile = None
        for ext in file_ext:
            fname_data = prefix + date_suffix + ext
            exists = os.path.isfile(datapath + '/GOES-R/' + fname_data)
            if exists:
                foundfile = fname_data
            
        #Try versions
        if foundfile == None:
            for ext in file_ext:
                fname_data = prefix + date_suffix + ext
                url=('https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/%s/l2/data/sgps-l2-avg5m/%i/%02i/%s' % (satellite,year,month,fname_data))
                try:
                    urllib.request.urlopen(url)
                    wget.download(url, datapath + '/GOES-R/' + fname_data)
                    foundfile = fname_data
                    break
                except urllib.request.HTTPError:
                    foundfile = None

        if foundfile == None:
            sys.exit("Cannot access GOES-R file at " + url +
               ". Tried file versions " + str(file_ext) + ". Please check that selected spacecraft covers date range.")
  
        filenames1.append('GOES-R/' + foundfile)
        
    return filenames1


def check_goesR_RTdata(startdate, enddate, experiment, flux_type):
    """Check that GOES-R data is on your computer or download it from the NOAA
        website. Return the filenames associated with the correct GOES data.
        GOES-R real time integral files are saved daily in txt format.
        
        The GOES-16 & 17 integral fluxes are the real time product plotted
        by SWPC on a daily basis and archived at CCMC. These fluxes are
        not the official science-grade integral flux product from NOAA, as
        these are not yet available. Also, only the PRIMARY spacecraft
        integral fluxes are available and it isn't possible to choose
        between GOES-16 or GOES-17 for the integral fluxes. When the official
        integral fluxes are available, it will be possible to select the
        spacecraft.
        
        INPUTS:
        
        :startdate: (datetime) start of time period specified by user
        :enddate: (datetime) end of time period entered by user
        :experiment: (string) name of native experiment or "user"
        :flux_type: (string) "integral" or "differential"
        
        OUTPUTS:
        
        :filenames1: (string array) the files containing the GOES
            EPS or EPEAD data that span the desired time range
            (monthly files)
        :filenames2: (string array) the files containing the GOES
            HEPAD data that span the desired time range
        :filenames_orien: (string array) the files
            that indicate the orientation of the GOES EPS or
            EPEAD detector (so can choose westward facing detector)
        
    """
    if flux_type == "differential":
        sys.exit("check_goesR_RTdata: This subroutine is only valid for GOES-R "
                "integral fluxes. Please set the flux_type to integral and try "
                "again.")
    
    styear = startdate.year
    stmonth = startdate.month
    stday = startdate.day
    endyear = enddate.year
    endmonth = enddate.month
    endday = enddate.day

    #Array of filenames that contain the data requested by the User
    filenames1 = []  #GOES-R


    #GOES-R data is stored in daily data files
    td = enddate - startdate
    NFILES = td.days #number of data files to download
    if td.seconds > 0: NFILES = NFILES + 1

    #pulls primary spacecraft fluxes
    prefix = '_Gp_part_5m'


    #for every day that data is required, check if file is present or
    #needs to be downloaded.
    for i in range(NFILES):
        date = startdate + datetime.timedelta(days=i)
        year = date.year
        month = date.month
        day = date.day
        date_suffix = '%i%02i%02i' % (year,month,day)
 
        fname1 = date_suffix + prefix + '.txt'
        exists1 = os.path.isfile(datapath + '/GOES-R/' + fname1)
        filenames1.append('GOES-R/' + fname1)

        if not exists1:
            url=('https://iswa.gsfc.nasa.gov/iswa_data_tree/observation/magnetosphere/goes_p/particle/%i/%02i/%s' % (year,month,fname1))
            try:
                urllib.request.urlopen(url)
                wget.download(url, datapath + '/GOES-R/' + fname1)
            except urllib.request.HTTPError:
                sys.exit("Cannot access GOES-R file at " + url +
               ". Please check that selected spacecraft covers date range.")
  

    return filenames1




def check_data(startdate, enddate, experiment, flux_type):
    """Check that the files containing the data are in the data directory. If
        the files for the requested dates aren't present, they will be
        downloaded from the NOAA website. For SEPEM (RSDv2) data, if missing,
        the program prints the URL from which the data can be downloaded and
        unzipped manually.
        The RSDv2 data set is very large and takes a long time to read as a
        single file. This program will generate files containing fluxes for
        each year for faster reading.
        
        INPUTS:
        
        :startdate: (datetime) start of time period specified by user
        :enddate: (datetime) end of time period entered by user
        :experiment: (string) name of native experiment or "user"
        :flux_type: (string) "integral" or "differential"
        
        OUTPUTS:
        
        :filenames1: (string array) the files containing the data that
            span the desired time range (monthly files)
        
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

 
    if flux_type == "differential":
        filenames1 = check_goesR_data(startdate, enddate, experiment, flux_type)
        return filenames1
        
    if flux_type == "integral":
        filenames1 = check_goesR_RTdata(startdate, enddate, experiment, flux_type)
        return filenames1
        

    return filenames1




def read_in_goesR(experiment, flux_type, filenames1):
    """Read in GOES-R data from your computer.
        Appears that only differential channels + one >500 MeV
        integral channel are available in the files.
        
        Reading in Level 2 data.
        Flux fill value = -1e+31
        Differential flux units = protons/(cm^2 sr keV s)
        
        Time stamp is seconds since 2000-01-01 12:00:00
        
        +X and -X are stored in the same array as the 0th and 1st
        array entry... figuring out which are which
        Typicaly, -X faces West and +X faces East
        
        YawFlipFlag indicates if the spacecraft is flipped.
        Assume a value of 1 indicates the spacecraft is flipped.
        
        https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes16/l1b/docs/GOES-16_SEISS_SGPS_L1b_Provisional_Maturity_ReadMe.pdf
        "There are two SGPS sensor units mounted on each GOES-R series spacecraft, facing in the spacecraft -X and +X directions. When the spacecraft is not in the yaw-flipped configuration SGPS-X faces west and SGPS+X faces east. Each SGPS unit has three solid-state (silicon detector) telescopes T1, T2, and T3 for measuring 1-25, 25-80, and 80-500 MeV protons, respectively. All three telescopes have the same look direction (i.e., +X or -X). T1 and T2 have 60o (full cone angle) fields of view, and T3 has a 90o field of view. Each unit measures 1-500 MeV proton fluxes in 13 logarithmically spaced differential channels (P1-P10) and >500 proton flux in a single integral channel (P11). The L1b data product is one-second cadence fluxes. The channels generally register counts above backgrounds only during solar energetic particle events, except for P11 which measures galactic cosmic rays in the absence of a solar particle event."
        
        
        INPUTS:
        
        :experiment: (string) experiment name
        :flux_type: (string) integral or differential
        :filenames1: (string array) the files containing the GOES-R
            netcdf data that span the desired time range
        
            
        OUTPUTS:
        
        :all_dates: (datetime 1xm array) time points for every time in
            all the data points in the files contained in filenames1
        :all_fluxes: (float nxm array) fluxes for n energy channels and m
            time points
    
        Note that all_dates and all_fluxes will be trimmed down to the
        user time period of interest.
        
    """
    ndiff_chan = 13 #5
    conversion = 1000. #keV/MeV
    
    NFILES = len(filenames1)
    all_dates = []
    all_fluxes = []

    #GOES-R times are wrt reference time of 2000-01-01 12:00:00
    ref_date = datetime.datetime(year=2000, month=1, day=1, hour=12)
    
    #Read in fluxes from files
    for i in range(NFILES):
        infile = os.path.expanduser(datapath + "/" + filenames1[i])
        data = netCDF4.Dataset(infile)
        
        if "v3-0" in filenames1[i]:
            ntstep = len(data.variables["time"])
        else:
            ntstep = len(data.variables["L2_SciData_TimeStamp"])
        
        #13 differential channels, one integral channel
        #5 minute time steps
        fluxes = np.zeros(shape=(ndiff_chan+1,ntstep))
        
        for j in range(ntstep):
            if "v3-0" in filenames1[i]:
                time_sec = float(data.variables["time"][j].data)
            else:
                time_sec = float(data.variables["L2_SciData_TimeStamp"][j].data)
            td = datetime.timedelta(seconds=time_sec)
            date = ref_date + td
            all_dates.append(date)
            
            #Orientation flag
            if "v3-0" in filenames1[i]:
                flip_flag = data.variables["yaw_flip_flag"][j]
            else:
                flip_flag = data.variables["YawFlipFlag"][j]
            idx = 0
            if flip_flag == 1:
                idx = 1
            #if flip_flag > 1:
            #    idx = None #exclude these points because in process of flip
            
            #Extract the 13 differential channels
            for k in range(ndiff_chan):
                ###TEMP###
                #kk = k + 8
                #[288 time step, 2 +/-X, 13 energy chan]
                
                flux = data.variables["AvgDiffProtonFlux"][j][idx][k]
                if flux < 0:
                    flux = badval
                fluxes[k][j] = flux*conversion
            

            flux = data.variables["AvgIntProtonFlux"][j][idx]
            if flux < 0:
                flux = badval
            fluxes[-1][j] = flux
                
        if all_fluxes == []:
            all_fluxes = fluxes
        else:
            all_fluxes = np.concatenate((all_fluxes,fluxes),axis=1)


    return all_dates, all_fluxes



def read_in_goesR_RT(experiment, flux_type, filenames1):
    """ Read in GOES-R data from your computer.
        Read in the NOAA SWPC real time integral flux files from
        the 1 day json files for the primary GOES spacecraft.
        These files are archived on the CCMC website.
        
        The >500 MeV fluxes are in the differential files. This
        file contains fluxes up to >100 MeV.
        
        Reading in Level 2 data.
        Flux fill value = -100000.0
        Units: Particles = Protons/cm2-s-sr
        Units: Electrons = Electrons/cm2-s-sr
        
        Time stamp is YYYYY M D HHMM: 2022  1  20  0000
        
        No selection for +/-X direction. Assume always westward facing.
        
        
        INPUTS:
        
        :experiment: (string) experiment name
        :flux_type: (string) integral or differential
        :filenames1: (string array) the files containing the GOES-R
            netcdf data that span the desired time range
        
            
        OUTPUTS:
        
        :all_dates: (datetime 1xm array) time points for every time in
            all the data points in the files contained in filenames1
        :all_fluxes: (float nxm array) fluxes for n energy channels and m
            time points
    
        Note that all_dates and all_fluxes will be trimmed down to the
        user time period of interest.
        
    """
    n_chan = 6
    
    NFILES = len(filenames1)
    all_dates = []
    all_fluxes = []

    
    #Read in fluxes from files
    for i in range(NFILES):
        with open(datapath + "/" + filenames1[i]) as infile:
        
            #6 integral channels
            #5 minute time steps up to 00:00 of the next day
            #Exclude last time step at midnight of next day or will
            #end up with repeat entries for the same time.
            fluxes = np.zeros(shape=(n_chan,288))
            j = 0 #counter for number of time steps in file
        
            for row in infile:
                if row[0] == "#": continue
                if row[0] == ":": continue
                if j == 288: continue
                
                row = row.split()
                
                #Get Date
                year = int(row[0])
                month = int(row[1])
                day = int(row[2])
                hour = int(row[3][0:2])
                minute = int(row[3][2:4])
                date = datetime.datetime(year=year,month=month,day=day,hour=hour,minute=minute)
                all_dates.append(date)
                
                for k in range(n_chan):
                    flux = float(row[6+k])
                    if flux < 0:
                        flux = badval
                    fluxes[k][j] = flux
                
                j = j+1 #count dates
                    
        if all_fluxes == []:
            all_fluxes = fluxes
        else:
            all_fluxes = np.concatenate((all_fluxes,fluxes[:,0:len(all_dates)]),axis=1)

    return all_dates, all_fluxes



def read_in_files(experiment, flux_type, filenames1):
    """ Read in the appropriate data files with the correct format. Return an
        array with dates and fluxes. Bad flux values (any negative flux) are set
        to -1. Format is defined to work with the files downloaded directly from
        NOAA or the RSDv2 (SEPEM) website as is.
 
 
        INPUTS:
        
        :experiment: (string) name of native experiment or "user"
        :flux_type: (string) "integral" or "differential"
        :filenames1: (string array) the files containing the data that
            span the desired time range
        
        OUTPUTS:
        
        :all_dates: (datetime 1xm array) time points for every time in
            all the data points in the files contained in filenames1
        :all_fluxes: (float nxm array) fluxes for n energy channels and m
            time points
   
        Note that all_dates and all_fluxes will be trimmed down to the
        user time period of interest.
       
    """
    print('Reading in data files for ' + experiment + '.')
    all_dates = []
    all_fluxes = []

    if (experiment == "GOES-16" or experiment == "GOES-17") and flux_type == "differential":
        all_dates, all_fluxes = read_in_goesR(experiment, \
                    flux_type, filenames1)
        return all_dates, all_fluxes
        
    if (experiment == "GOES-16" or experiment == "GOES-17") and flux_type == "integral":
        all_dates, all_fluxes = read_in_goesR_RT(experiment, \
                    flux_type, filenames1)
        return all_dates, all_fluxes

    return all_dates, all_fluxes




def do_interpolation(i,dates,flux):
    """ If bad fluxes (flux < 0) are found in the data, find the first prior
        data point and the first following data point that have good flux values.
        Perform linear interpolation in time:
        
        F(t) = F1 + (t - t1)*(F2 - F1)/(t2 - t1)
        
        This subroutine does the calculation for a single instance of bad data
        that corresponds to array index i.
        
        INPUTS:
        
        :i: (integer) index of time point to interpolate
        :dates: (datetime 1xp array) dates for p time points within
            the date range specified by startdate and enddate
        :fluxes: (float 1xp array) flux time profiles for n energy channels
            and p time points
            
        OUTPUTS:
        
        :interp_flux: (float 1xp array) flux time profile with any negative
            flux values replaced with linear interpolation with time
       
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
    """ Search the data for bad values (flux < 0) and fill the missing data with
        an estimate flux found by performing a linear interpolation with time,
        using the good flux values immediately surrounding the data gap.
        
        INPUTS:
        
        :dates: (datetime 1xp array) dates for p time points within
            the date range specified by startdate and enddate
        :fluxes: (float nxp array) flux time profiles for n energy channels
            and p time points
        :energy_bins: (float nx2 array) energy bins associated with fluxes
        :dointerp: (bool) Set True to perform linear interpolation in time,
            otherwise will fill bad data points with None values
            
        OUTPUT:
        
        :fluxes: (float 1xp array) flux time profile with any negative
            flux values replaced with linear interpolated of None values
       
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
                    print('There is a data gap for time ' + str(dates[j])
                            + ' and energy bin ' + str(energy_bins[i][0]) + ' - '
                            + str(energy_bins[i][1]) + '.'
                            + ' Filling in missing value with linear '
                            + 'interpolation in time.')
                    interp_flux = do_interpolation(j,dates,fluxes[i,:])
                    fluxes[i,j] = interp_flux
                else:
                    print('There is a data gap for time ' + str(dates[j])
                            + ' and energy bin ' + str(energy_bins[i][0]) + ' - '
                            + str(energy_bins[i][1]) + '.'
                            + ' Filling in missing value with None ')
                    fluxes[i,j] = None

            elif fluxes[i,j] < 0:
                #estimate flux with interpolation in time
                if dointerp:
                    print('There is a data gap for time ' + str(dates[j])
                            + ' and energy bin ' + str(energy_bins[i][0]) + ' - '
                            + str(energy_bins[i][1]) + '.'
                            + ' Filling in missing value with linear '
                            + 'interpolation in time.')
                    interp_flux = do_interpolation(j,dates,fluxes[i,:])
                    fluxes[i,j] = interp_flux
                else:
                    print('There is a data gap for time ' + str(dates[j])
                            + ' and energy bin ' + str(energy_bins[i][0]) + ' - '
                            + str(energy_bins[i][1]) + '.'
                            + ' Filling in missing value with None ')
                    fluxes[i,j] = None #results in NaN value in np array

    
    print('Finished checking for bad data.')
    print()
    return fluxes



def define_energy_bins(experiment,flux_type):
    """ Define the energy bins for the selected spacecraft or data set.
        GOES-R only
        
        INPUTS:
        
        :experiment: (string) name of experiment or "user"
        :flux_type: (string) integral or differential

        
        OUTPUTS:
        
        :energy_bins: (float nx2 array) appropriate energy bins for the
            experiment specified by the user
       
    """
    #use corrected proton flux for GOES eps or epead; include hepad
    #-1 indicates infinity ([700, -1] means all particles above 700 MeV)

    if experiment == "GOES-16" or experiment == "GOES-17":
        if flux_type == "differential":
            #energy_bins = [[83.7,98.5],
            #               [99.9,118.0],[115.0,143.0],[160.0,242.0],
            #               [276.0,404.0],[500.0,-1]]
            energy_bins = [[1.02,1.86],[1.9,2.3],[2.31,3.34],
                           [3.4,6.48],[5.84,11.0],[11.64,23.27],
                           [24.9,38.1],[40.3,73.4],[83.7,98.5],
                           [99.9,118.0],[115.0,143.0],[160.0,242.0],
                           [276.0,404.0],[500.0,-1]]
        if flux_type == "integral":
            energy_bins = [[1,-1],[5,-1],[10,-1],[30,-1],[50,-1],[100,-1]]

    return energy_bins
