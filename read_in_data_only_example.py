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


def modified_weibull(times, Ip, a, b):
    """ Create a Weibull for times in seconds since first
        date with params a and b.
    """
    
    weibull = []
    for t in times:
        W = Ip*(-a/b)*(t/b)**(a-1)*math.exp(-(t/b)**a)
        weibull.append(W)
        
    return weibull



def residual(fit, data):
    """ Calculate difference between fit and data."""
    
    resid = []
    for i in range(len(fit)):
        resid.append(data[i] - fit[i])
    
    return resid
    


def weibull_residual(params, *args):
    """ Caluate the residual of the Weibull fit
        compared to data.
    """
    
    pars = params.valuesdict()
    a = pars['alpha']
    b = pars['beta']
    Ip = pars['peak_intensity']
    
    #print('Generating Weibull Ip: ' + str(Ip) + ', a: ' + str(a) + ', b: ' + str(b))
    
    times = args[0]
    data = args[1]
    
    fit = modified_weibull(times, Ip, a, b)
    
    resid = residual(fit, data)
    #print(resid)
    return resid
    

def find_max_curvature(x, y):
    """ Calculate the curvature along a curve
        and find the maximum curvature location.
        
        For any three points (a,b,c) the curvature is
        2 * |(a-b) x (b-c)| / (|a-b| * |b-c| * |c-b|).
        
        https://stackoverflow.com/questions/59041924/how-can-i-go-about-finding-points-where-there-is-a-bend-cut-in-my-data
        
        INPUTS:
        
        :curve: (float 1xn array) weibull fit points
    
    """
    xarr = np.array(x)
    yarr = np.array(y)
    yderiv = yarr[1:] - yarr[:-1]
    yderiv2 = yderiv[1:] - yderiv[:-1]
    
    k_x = yderiv2/((1 + yderiv[1:]**2))**(3./2.)
    
    max_k_idx= np.argmin(k_x)
    
    #rescale the curvature to overplot
    max_y = np.max(yarr)
    max_y_idx = np.argmax(yarr)
    k_x = (np.max(yarr)/np.max(k_x))*k_x
    
    plt.figure("Curvature",figsize=(9,5))
    plt.plot(xarr,yarr,label="orig")
    plt.plot(xarr[max_k_idx+2], yarr[max_k_idx+2],"o",label="max curvature on Weibull")
    plt.plot(xarr[max_y_idx], yarr[max_y_idx],"o",label="max Weibull")
    plt.plot(xarr[2:],k_x,label="curvature")
    plt.plot(xarr[max_k_idx+2],k_x[max_k_idx],"o",label="max curvature")
    plt.legend(loc='lower left')
    plt.xlabel("Hours")
    plt.ylabel("Y")
    #plt.yscale("log")
    #plt.ylim(1e-4,1e6)

    plt.show()

    return max_k_idx+2
    


def calculate_onset_peak_from_fit(experiment, energy_bins, dates, fluxes):
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
        
        The onset peak is found by fitting a Weibull function to the SEP
        time profile. A peak time is estimated and an appropriate flux
        value at or near that peak will be identified as the peak value.
        
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
    nbins = len(energy_bins)
    
    #Convert dates into a series of times for fitting
    #Seconds since first date
    times = [((t - dates[0]).total_seconds() + 60)/(60*60) for t in dates]
    
    #Do a fit of the Weibull function for each time profile
    params_weib = Parameters()
    params_weib.add('alpha', value = -3, min = -5, max = -0.1)
    params_weib.add('beta', value = 10, min = 1, max =100)
    params_weib.add('peak_intensity', value = 100, min = 1e-3, max =1e6)
    
    for i in range(nbins):
        if energy_bins[i][0] < 10: continue
        #log_fluxes = [math.log10(f) for f in fluxes[i]]
        minimize_weib = minimize(weibull_residual, params_weib, args = [times, fluxes[i]], nan_policy= 'propagate', max_nfev=np.inf)
        
        best_pars = minimize_weib.params.valuesdict()
        best_a = best_pars['alpha']
        best_b = best_pars['beta']
        best_Ip = best_pars['peak_intensity']
        best_weibull = modified_weibull(times, best_Ip, best_a, best_b)
        print("==== >" + str(energy_bins[i][0]) + " MeV =====")
        print('Best Fit Weibull Ip: ' + str(best_Ip) + ', a: ' + str(best_a) + ', b: ' + str(best_b))
        
        max_curve_idx = find_max_curvature(times, best_weibull)
        max_curve_time = times[max_curve_idx]
        max_curve_val = best_weibull[max_curve_idx]
        
        
        max_val = np.max(best_weibull)
        max_idx = np.where(best_weibull == max_val)
        max_time = times[max_idx[0][0]]
        
        #Pull out max measured value around this identified maximum in the fit
        model_max_date = datetime.timedelta(seconds=(max_time*60*60 - 60)) + dates[0]
        max_meas = 0
        max_date = 0
        dt = datetime.timedelta(hours=1)
        for j in range(len(dates)):
            if dates[j] >= model_max_date - dt and dates[j] <= model_max_date + dt:
                if fluxes[i][j] > max_meas:
                    max_meas = fluxes[i][j]
                    max_date = times[j] #dates[j]
        
        
        plt.figure("Onset Peak",figsize=(9,5))
        label = ">" + str(energy_bins[i][0]) + " MeV"
        plt.plot(times,fluxes[i],label=label)
        label_weib = "Weibull "
        plt.plot(times,best_weibull,label="Weibull")
        #test
        rise = []
        exp_fall = []
        test_weib = []
        for t in times:
            rise.append(best_Ip*(-best_a/best_b)*(t/best_b)**(best_a-1))
            exp_fall.append(best_Ip*(-best_a/best_b)*math.exp(-(t/best_b)**best_a))
            test_weib.append(best_Ip*(-best_a/best_b)*(t/best_b)**(best_a-1)*math.exp(-(t/best_b)**best_a))
       # plt.plot(times,rise,label="rise (t/b)^(a-1)")
       # plt.plot(times,exp_fall,label="exp fall e^(-(t/b)^a)")
       # plt.plot(times,test_weib,":",label="test Weibull")
        plt.plot(max_time, max_val,"o",label="max model")
        plt.plot(max_date, max_meas,">",label="max measured")
        plt.plot(max_curve_time, max_curve_val,"D",label="max curvature")
        plt.legend(loc='lower left')
        plt.xlabel("Hours")
        plt.ylabel("Intensity")
        plt.yscale("log")
        plt.ylim(1e-4,1e6)




def run_all(str_startdate, str_enddate, experiment,
        flux_type, model_name, user_file, showplot, saveplot,
        options, nointerp, int_to_diff):
        
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
    
    exposure_time = (dates[-1] - dates[0]).total_seconds()
    
    if int_to_diff:
        energy_bins, fluxes = make_diff_from_int(fluxes, energy_bins)
    
    
#    #Code for plotting fluence with various normalization factors
#    #to try to match Kathryn Whitman's thesis
#    print("===Fluence values for time period===")
#    bin_centers = []
#    fluences = []
#    for ii in range(len(energy_bins)): #MeV
#        bin_center = math.sqrt(energy_bins[ii][0]*energy_bins[ii][1])
#        bin_centers.append(bin_center)
#        fluence = sum(fluxes[ii][:])*5.*60./(exposure_time)
#        fluences.append(fluence)
#        print(str(bin_center) + " MeV " + ", " + str(fluence) + " [MeV cm2 s sr]^-1")
#
#        for ii in range(len(energy_bins)): #GeV
#            bin_center = math.sqrt(energy_bins[ii][0]*energy_bins[ii][1])*1e-3
#            bin_centers.append(bin_center)
#            fluence = sum(fluxes[ii][:])*5.*60.*1e7/(exposure_time)
#            fluences.append(fluence)
#            print(str(bin_center) + " GeV " + ", " + str(fluence) + " [GeV m2 s sr]^-1")
#
#            flux_units_differential = "GeV^-1*m^-2*s^-1*sr^-1"
#            fluence_units_differential = "GeV^-1*m^-2"
#            energy_units = "GeV"
                        
    nbins = len(energy_bins)
    
    calculate_onset_peak_from_fit(experiment, energy_bins, dates, fluxes)
        
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
    
    
        #All energy channels in specified date range with event start and stop
        print("Generating figure of fluxes in original energy bins. Any bad data "
              "points were interpolated. Lines indicate event start and stop for "
              "thresholds.")
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
        #plt.grid(which="major", axis="both", linestyle="dotted")
        chartBox = ax.get_position()
        ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.85,
                         chartBox.height])
        ax.legend(loc='upper center', bbox_to_anchor=(1.17, 1.05))
        if saveplot:
            fig.savefig(plotpath + '/' +figname + '.png')
        if not showplot:
            plt.close(fig)
    
        if int_to_diff:
            #Event-integrated fluence for energy channels
            print("Generating figure of event-integrated fluence spectrum.")
            #Plot fluence spectrum summed between SEP start and end dates
            figname = 'Fluence' \
                    + '_' + experiment + '_' + flux_type + modifier \
                    + '_'
            if experiment == 'user' and model_name != '':
                figname = 'Fluence' \
                        + '_' + model_name + '_' + flux_type + modifier \
                        + '_'
            fig = plt.figure(figname,figsize=(6,5))
            ax = plt.subplot(111)
            markers = ['o','P','D','v','^','<','>','*','d','+','8','p','h','1','X','x']
            #for j in range(1): #len(energy_thresholds)):
               # legend_label = '>' + plt_energy[j] + ' ' + energy_units + ', '\
               #             + plt_flux[j] + ' ' + flux_units_integral
               # if plot_diff_thresh[j]: #tacked on to end
               #     legend_label = (plt_energy[j] + ' ' + energy_units + ',\n'
               #                 + plt_flux[j] + '\n' + flux_units_differential)
            ax.plot(bin_centers,fluences,"o")
            plt.grid(which="both", axis="both")
            plt.title(experiment + ' ' + title_mod + '\n Event-Integrated Average Flux')
            if experiment == 'user' and model_name != '':
                plt.title(model_name + ' ' + title_mod + '\n Event-Integrated '
                        'Average Flux')
            plt.xlabel('Energy [' + energy_units +']')
            #if flux_type == "integral":
            #    plt.ylabel('Integral Fluxes [' + flux_units_integral + ']')
            #if flux_type == "differential":
            plt.ylabel('Flux [' + flux_units_differential + ']')
            plt.xscale("log")
            plt.yscale("log")
            ax.legend(loc='upper right')
            if saveplot:
                fig.savefig(plotpath + '/' + figname + '.png')
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
    int_to_diff = args.IntToDiff



    dates,fluxes,energy_bins = run_all(str_startdate, str_enddate, experiment,
        flux_type, model_name, user_file, showplot, saveplot,
        options, nointerp, int_to_diff)

    if showplot: plt.show()
