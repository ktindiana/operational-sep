#This is a program that shows how to read in data using read_datasets.py
#The goal of this program is to get the data into arrays that are easily 
#used for other purposes.
#Also some plotting code
#CODE EXPECTS THERE TO BE A data/ DIRECTORY
#K. Whitman 2022-02-02


from library import read_datasets as datasets
from library import global_vars as vars
import numpy as np
import matplotlib.pyplot as plt
import datetime
import argparse
import math
import sys
from lmfit import minimize, Parameters
from statistics import mode
from scipy.odr import Model, Data, RealData, ODR

#From global_vars.py
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
#From global_vars.py
#(expect the first (0th) column contains date in YYYY-MM-DD HH:MM:SS format)
#Identify columns containing fluxes you want to analyze
user_col = vars.user_col
#DELIMETER between columns; for whitespace separating columns, use " " or ""
user_delim = vars.user_delim
#DEFINE ENERGY BINS associated with user file and columns specified above as:
user_energy_bins = vars.user_energy_bins
############################



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
    


def read_in_flux_files(experiment, flux_type, user_file, model_name, startdate,
        enddate, options, nointerp):
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
        :options: (string array) - options that could be applied
        :nointerp: (bool) - indicates if user doesn't want to do linear
            interpolation in time for negative of bad flux values
        
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
    
    #Optional - if bins in backwards order or something, this can be useful
    #all_fluxes, energy_bins = sort_bin_order(all_fluxes, energy_bins)
    
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
      

def make_diff_from_int(integral_fluxes, energy_bins):
    """ Make differential fluxes from integral channels.
        Subtract each integral flux channel from the consecutive
        one to create differential bins.
    """

    new_bins = []
    nchan = len(integral_fluxes)
    nflx = len(integral_fluxes[0])
    diff_fluxes = np.zeros((nchan-1,nflx))
    diff_bins = np.zeros((nchan-1,2))
    
    for i in range(nchan-1):
        diff_bins[i][0] = energy_bins[i][0]
        diff_bins[i][1] = energy_bins[i+1][0]
        bin_width = diff_bins[i][1] - diff_bins[i][0]
        for j in range(nflx):
            diff_fluxes[i][j] = \
                (integral_fluxes[i][j] - integral_fluxes[i+1][j])/bin_width
            
    
    print("make_diff_from_int converted to energy bins " + str(diff_bins))
    
    return diff_bins, diff_fluxes



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
    if not time_diff:
        sys.exit("determine_time_resolution: Require more than 1 data point "
                "to determine time resolution. Please extend your "
                "requested time range and try again. Exiting.")
    time_resolution = mode(time_diff)
    return time_resolution



def integrate_with_time(int_time, dates, fluxes):
    """ Enter a time resolution in units of seconds and
        create spectra integrated over each
        consecutive time point.
        
        The spectra retain the initial units and are not
        fluence spectra.
        
        INPUTS:
        
        :int_time: (float or integer) desired integration time
            in units of seconds
        :dates: (1xn datetime array) dates corresponding to the fluxes
        :fluxes: (mxn float array) fluxes for m energy bins
        
        OUTPUTS:
        :int_dates: (1xp datetime array) dates for p time periods;
            the date is at the middle of the integrated time period
        :int_fluxes: (mxp float array) fluxes integrated over desired
            time frame
        
    """
    
    int_dates = []
    int_fluxes = []
    time_res = determine_time_resolution(dates).total_seconds()
    
    fluxes_np = np.array(fluxes)

    n = math.floor((dates[-1] - dates[0]).total_seconds()/int_time) + 1
    
    for i in range(n):
        start = dates[0] + datetime.timedelta(seconds=i*int_time)
        end = start + datetime.timedelta(seconds=int_time)
        mid = start + datetime.timedelta(seconds=int_time/2.)
        
        if end > dates[-1]: continue
        
        int_dates.append(mid)
        
        idx = []
        for j in range(len(dates)):
            if dates[j] >= start and dates[j] < end:
                idx.append(j)
       
        print("integrate_with_time: Summing between " + str(dates[idx[0]]) \
                + " and " + str(dates[idx[-1]]))
        sum_flux = np.sum(fluxes[:,idx[0]:idx[-1]+1],axis=1,keepdims=True)
        #retain flux units with s^-1 by dividing by number of time points
        #that were summed together
        sum_flux = sum_flux/(len(idx))
        if int_fluxes == []:
            int_fluxes = sum_flux
        else:
            int_fluxes = np.append(int_fluxes,sum_flux,axis=1)


    return int_dates, int_fluxes


def ER(beta,E):
    """ E is the particle energies or energy bins
        beta[0] is a normalization
        beta[1] is gamma
        beta[2] is rollover energy E0
    """
    y = beta[0]*E**(-beta[1])*np.exp(-E/beta[2])
    return y


def fit_ER(energy_bins, flux):
    """ Fit a flux spectrum with the Ellison-Ramatay function.
    
        INPUTS:
        
        :energy_bins: (2xn float array) energy bins
        :flux: (1xn float array) spectrum at a single time point
        
        OUTPUS:
        :beta: (ODR object) fit parameters
    """

    energies = []
    for bin in energy_bins:
        energies.append(pow(bin[0]*bin[1],0.5))
    
    energies_np = np.array(energies)
    flux_np = np.array(flux)

    data = RealData(energies_np,flux_np)
    model = Model(ER)
    odr = ODR(model, data, beta0=[100.,2.,50.])
    odr.set_job(fit_type=2)
    output = odr.run()
    output.pprint()

    yn = ER(output.beta, energies_np)

    if showplot or saveplot:
        #=========FLUX TO PHI CORRELATION FOR IONS Z > 2 ========
        title = "Fit to " + experiment + " Spectrum \n" \
                    + r'fit = ' + "({:.4f}".format(output.beta[0]) + "*E$^{" \
                    + "{:.4f}".format(output.beta[1]) + "}*(E/" + "{:.4f}".format(output.beta[2]) +")"
            
            
        fig = plt.figure(figsize=(10,7))
        ax = plt.subplot(111)
        plt.title(title)
        plt.xlabel("Energy (MeV)")
        plt.ylabel(r'Summed Flux (cm$^{-2}$ s$^{-1}$ sr$^{-1}$ MeV$^{-1}$)')
        
        ax.plot(energies,flux,"o")
        ax.plot(energies,yn,label='Fit',color="red",zorder=10)
        
        chartBox = ax.get_position()
        ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.75,
                        chartBox.height])
        ax.legend(loc='upper center', bbox_to_anchor=(1.25, 0.95),ncol=2,fontsize=10)
        plt.xscale("log")
        plt.yscale("log")
        
 #       if saveplot:
 #           figname = vars.figpath+'/'+"IntegratedFlux_Fit.png"
  #          fig.savefig(figname, dpi=300, bbox_inches='tight')

    return beta


def run_all(str_startdate, str_enddate, experiment,
        flux_type, model_name, user_file, showplot, saveplot,
        options, nointerp, int_time, int_to_diff):
        
    """ :str_startdate: (string) - user input start date "YYYY-MM-DD" or
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
        :showplot: (bool) - Set to True to show plots when run
        :saveplot: (bool) - Set to True to automatically save plots to the
            plots directory when run
        :nointerp: (boolean) - set to true to fill in negative fluxes with None
            value rather than filling in via linear interpolation in time
        :int_to_diff: (boolean) - convert integral channels to differential
            channels by subtracting the higher channels from the lower channels
    """

    #Check for empty dates
    if (str_startdate == "" or str_enddate == ""):
        sys.exit('You must enter a valid date range. Exiting.')

    #PROCESS INPUTS
    options = options.split(";")
    user_fname = [user_file] #input as argument, default is 'tmp.txt'
    
    startdate = str_to_datetime(str_startdate)
    enddate = str_to_datetime(str_enddate)
    
    #READ IN FLUXES AND ENERGY BINS, BG SUBTRACT, INTERPOLATE
    dates, fluxes, energy_bins = read_in_flux_files(experiment,
        flux_type, user_file, model_name, startdate,
        enddate, options, nointerp)
    
    int_dates, int_fluxes = integrate_with_time(int_time, dates, fluxes)
    
    #beta = fit_ER(energy_bins, int_fluxes[:,10])

                        
    nbins = len(energy_bins)
    
    #===============PLOTS==================
    if saveplot or showplot:
        #Additions to titles and filenames according to user-selected options
        modifier = ''
        title_mod = ''
        if "uncorrected" in options:
            modifier = modifier + '_uncorrected'
            title_mod = title_mod + 'uncorrected '
        if "S14" in options:
            modifier = modifier + '_S14'
            title_mod = title_mod + 'S14 '
        if "Bruno2017" in options:
            modifier = modifier + '_Bruno2017'
            title_mod = title_mod + 'Bruno2017 '
    
    
        ###PLOT ORIGINAL INPUT FLUXES WITH TIME####
        #Plot all channels of user specified data
        figname = 'Fluxes' \
                + '_' + experiment + '_' + flux_type + modifier \
                + '_' + 'All_Bins'
        if experiment == 'user' and model_name != '':
            figname = 'Fluxes' \
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

            ax.plot_date(dates,fluxes[i],'-',label=legend_label)

        if flux_type == "integral":
            plt.ylabel('Integral Flux [' + flux_units_integral + ']')
            plt.title(experiment + ' '+ title_mod)
            if experiment == 'user' and model_name != '':
                plt.title(model_name + ' '+ title_mod)
        if flux_type == "differential":
            plt.ylabel('Flux [' + flux_units_differential + ']')
            plt.title(experiment + ' ' + title_mod)
            if experiment == 'user' and model_name != '':
                plt.title(model_name + ' ' + title_mod)
        plt.xlabel('Date')
        plt.yscale("log")
        #plt.grid(which="major", axis="both", linestyle="dotted")
        chartBox = ax.get_position()
        ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.85,
                         chartBox.height])
        ax.legend(loc='upper center', bbox_to_anchor=(1.17, 1.05))
        if saveplot:
            fig.savefig(plotpath + '/' +figname + '.png')
        if not showplot:
            plt.close(fig)
    


        ####PLOT FLUXES INTEGRATED OVER SPECIFIC TIME WITH TIME#####
        #Plot all channels of user specified data
        figname = 'IntegratedFluxes' \
                + '_' + experiment + '_' + flux_type + modifier \
                + '_' + 'All_Bins'
        if experiment == 'user' and model_name != '':
            figname = 'IntegratedFluxes' \
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

            ax.plot_date(int_dates,int_fluxes[i],'.',label=legend_label)

        if flux_type == "integral":
            plt.ylabel('Time Integrated Integral Flux [' + flux_units_integral + ']')
            plt.title(experiment + ' '+ title_mod)
            if experiment == 'user' and model_name != '':
                plt.title(model_name + ' '+ title_mod)
        if flux_type == "differential":
            plt.ylabel('Time Integrated Flux [' + flux_units_differential + ']')
            plt.title(experiment + ' ' + title_mod)
            if experiment == 'user' and model_name != '':
                plt.title(model_name + ' ' + title_mod)
        plt.xlabel('Date')
        plt.yscale("log")
        #plt.grid(which="major", axis="both", linestyle="dotted")
        chartBox = ax.get_position()
        ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.85,
                         chartBox.height])
        ax.legend(loc='upper center', bbox_to_anchor=(1.17, 1.05))
        if saveplot:
            fig.savefig(plotpath + '/' +figname + '.png')
        if not showplot:
            plt.close(fig)
    
    
    return dates,fluxes,energy_bins




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
            'GOES-16', 'GOES-17','SEPEM', 'SEPEMv3', 'EPHIN', 'EPHIN_REleASE',
            'user'],
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
    parser.add_argument("--IntegrationTime", type=int, default=3600, help=("Integration time to sum flux. Must be in seconds. "
                "Default=3600 (1 hour)"))
    parser.add_argument("--IntToDiff",
            help=("Convert integral channels to differential channels"
                    "by subtracting."), action="store_true")
    parser.add_argument("--showplot",
            help="Flag to display plots", action="store_true")
    parser.add_argument("--saveplot",
            help="Flag to save plots to file", action="store_true")


    args = parser.parse_args()

    str_startdate = args.StartDate
    str_enddate = args.EndDate
    experiment = args.Experiment
    flux_type = args.FluxType
    model_name = args.ModelName
    user_file = args.UserFile
    showplot = args.showplot
    saveplot = args.saveplot
    options = args.options
    nointerp = args.NoInterp
    int_time = args.IntegrationTime
    int_to_diff = args.IntToDiff



    dates,fluxes,energy_bins = run_all(str_startdate, str_enddate, experiment,
        flux_type, model_name, user_file, showplot, saveplot,
        options, nointerp, int_time, int_to_diff)

    if showplot: plt.show()