from library import read_datasets as datasets
from library import global_vars as vars
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

__version__ = "0.1"
__author__ = "Katie Whitman"
__maintainer__ = "Katie Whitman"
__email__ = "kathryn.whitman@nasa.gov"

datapath = vars.datapath
outpath = vars.outpath
plotpath = vars.plotpath
badval = vars.badval #bad data points will be set to this value; must be negative

nsigma = vars.nsigma

######FOR USER DATA SETS######
#(expect the first (0th) column contains date in YYYY-MM-DD HH:MM:SS format)
#Identify columns containing fluxes you want to analyze
user_col = vars.user_col
#DELIMETER between columns; for whitespace separating columns, use " " or ""
user_delim = vars.user_delim
#DEFINE ENERGY BINS associated with user file and columns specified above as:
user_energy_bins = vars.user_energy_bins
############################


def about_derive_background():
    """ About derive_background.py
        
        Perform background-subtraction and return only SEP flux values.

        1. Read in a user-specified time period for the background flux
        2. Estimate the mean flux and level of variation (sigma)
        3. Define all flux above mean+nsigma*sigma (nsigma is defined in global_vars)as SEP flux and all flux below as background flux
        4. Separate the SEP fluxes and subtract the mean background value
        5. Return the total flux, the background flux, and the background-subtracted SEP fluxes
    """
    

def remove_none(flux):
    """ Takes 1D array of flux and removes None values.
        None values in a list are converted to NaN values in a numpy array.
        So check for NaN and remove.
        
        INPUTS:
        
        :flux: (float 1xn array)
        
        OUTPUTS:
        
        :clean_flux: (float 1xm array) with None or NaN values removed
        
    """
    bad_index = np.argwhere(np.isnan(flux))
    clean_flux = np.delete(flux, bad_index)
    return clean_flux


def remove_zero(flux):
    """Takes 1D array of flux and removes zero values.
    
        INPUTS:
        
        :flux: (float 1xn array)
        
        OUTPUTS:
        
        :clean_flux: (float 1xm array) with zero values removed
        
    """
    bad_index = np.argwhere(flux == 0)
    clean_flux = np.delete(flux, bad_index)
    return clean_flux


def remove_above(flux, val):
    """ Remove any flux above a specific value, val.
        
        INPUTS:
        
        :flux: (float 1xn array)
        :val: (float) lower threshold value
        
        OUTPUTS:
        
        :strip_flux: (float 1xm array) with flux > val removed
    
    """
    indices = np.nonzero(flux > val)
    strip_flux = np.delete(flux, indices)
    return strip_flux


def remove_below(flux, val):
    """ Remove any flux below a specific value, val.
    
        INPUTS:
        
        :flux: (float 1xn array)
        :val: (float) upper threshold value
        
        OUTPUTS:
        
        :strip_flux: (float 1xm array) with flux < val removed
        
    """
    indices = np.nonzero(flux < val)
    strip_flux = np.delete(flux, indices)
    return strip_flux


def separate_sep_and_background(fluxes, dates, means, sigmas):
    """ Take the input fluxes, separate them into arrays containing
        the background flux and SEP flux. Values above mean + Nsigma*sigma is
        considered SEP flux while values below are considered the background.
        Perform a background subtraction on the SEP flux by subtracting the
        mean background value.
        The input flux array is a numpy array, but the output will be a list
        for flexibility.
        Nsigma is specified in library/global_vars.py.
        
        INPUTS:
        
        :fluxes: (float nxm array) fluxes for n energy channels and m time points
        :dates: (datetime 1xm array) time points for flux time profile
        :means: (float 1xn array) mean background flux for n energy channels
        :sigmas: (float 1xn array) expected variability sigma for n energy channels
        
        OUTPUTS:
        
        :bgfluxes: (float nxm array) background fluxes for n energy channels and
            m time points
        :sepfluxes: (float nxm array) SEP fluxes for n energy channels and
            m time points
        
    """
    nflx = len(fluxes)
    sepfluxes = []
    bgfluxes = []

    for i in range(nflx):
        bgflux = [-1]
        sepflux = [-1]
        for j in range(len(fluxes[i])):
            if fluxes[i][j] <= means[i] + nsigma*sigmas[i]:
                bgflux.append(fluxes[i][j])
                sepflux.append(0)
            if fluxes[i][j] > means[i] + nsigma*sigmas[i]:
                bgflux.append(0)
                bgsubflux = fluxes[i][j] - means[i]
                if bgsubflux < 0: bgsubflux = 0
                sepflux.append(bgsubflux)
            if np.isnan(fluxes[i][j]):
                bgflux.append(None)
                sepflux.append(None)

        bgflux.pop(0)
        sepflux.pop(0)

        if not bgfluxes:
            bgfluxes = [bgflux]
        else:
            bgfluxes.append(bgflux)

        if not sepfluxes:
            sepfluxes = [sepflux]
        else:
            sepfluxes.append(sepflux)

    return bgfluxes, sepfluxes


def define_hist_bins(flux):
    """Takes a 1D numpy array of flux and defines a set of histogram bins
        between the min and max flux values equally spaced in log space.
        
        INPUTS:
        
        :flux: (float 1xm array) flux time profile for m time steps
        
        OUTPUTS:
        
        :hist_bins: (float 1x20 (nbins) array) bins for a histogram equally
            spaced in log space over 20 betweens between the minimum and
            maximum flux values in flux
        
    """
    if flux.any(0):
        flux = remove_zero(flux)
    max_val = flux.max()
    min_val = flux.min()

    nbins = 20;
    bins = []
    #Create nbins bins in log space between min_val and max_val
    logmax = math.log10(max_val)
    logmin = math.log10(min_val)
    log_bins = np.linspace(start=logmin, stop=logmax, num=nbins + 1,\
                        endpoint=True)

    for bin in log_bins:
        if not bins:
            bins = [10**bin]
        else:
            bins.append(10**bin)

    hist_bins = np.array(bins)
    return hist_bins


def create_histogram(flux, energy_bin, iteration):
    """Take a list of flux with time and generate a histogram of the values.
        NaN values are removed prior to creating histogram.
        The histogram is created with bins extending from the min flux to the
        max flux and equally spaced in log space.
        Estimate the mean by averaging the bin centers weighted by the
        frequency. Calculate the variance and take the square root to estimate
        sigma. Return the mean and sigma.
        
        INPUTS:
        
        :flux: (float 1xm array) flux time profile for m time points
        :energy_bin: (float 2x1 array) energy bin that defines the energy
            channel for the flux array (only used for plotting)
        :iteration: (integer) indicates how many times the flux has gone
            through this subroutine (only used for plotting)
            
        OUTPUTS:
        
        :hist_mean: (float) weighted mean of the histogram
        :sigma: (float) standard deviation (sqrt(variance)) of
            the histogram
        
    """
    #Bad values in the data were set to None
    #Remove None values to calculate flux distribution in a histogram
    clean_flux = remove_none(flux) #remove any None values in numpy array

    #Define a set of logarithmic bins that span the flux values
    hist_bins = define_hist_bins(clean_flux)
    hist, _ = np.histogram(clean_flux, bins=hist_bins)

    #Check if the lowest histogram bin represents a data floor (as in GOES).
    #This is indicated by the lowest energy bin containing the most values
    if hist[0] == hist.max():
        #remove lowest energy bin and recalculate histogram
        hist_bins = np.delete(hist_bins, 0)
        hist, _ = np.histogram(clean_flux, bins=hist_bins)

    centers = []
    for i in range(len(hist_bins)-1):
        center = math.sqrt(hist_bins[i]*hist_bins[i+1])
        if not centers:
            centers = [center]
        else:
            centers.append(center)
    bin_centers = np.array(centers)
    hist_mean = np.average(bin_centers, weights=hist)
    variance = np.average((bin_centers - hist_mean)**2, weights=hist)
    sigma = np.sqrt(variance)

    # Plot histogram
    # figname = 'FluxHistogram_'+ str(energy_bin[0]) + '_' + str(energy_bin[1]) \
    #             + '_it' +str(iteration)
    # fig = plt.figure(figname,figsize=(8,5))
    # ax = plt.subplot(111)
    # n, bins, patches = plt.hist(x=clean_flux, bins=hist_bins, color='#0504aa',
    #                             alpha=0.7, rwidth=0.85)
    # plt.grid(axis='y', alpha=0.75)
    # plt.xlabel('Flux')
    # plt.ylabel('Frequency')
    # plt.title('Distribution of Flux Values for ' + str(energy_bin[0]) + ' to '
    #         + str(energy_bin[1]) + ' MeV')
    # plt.text(0.8, 0.5, (r'$\mu$= %.4g, sigma= %.4g'%(hist_mean, sigma)),\
    #         horizontalalignment='center', verticalalignment='center', \
    #         transform=ax.transAxes)
    # maxfreq = n.max()
    # # Set a clean upper y-axis limit.
    # plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)


    return hist_mean, sigma



def iterate_background(fluxes, energy_bins):
    """Bin fluxes into histograms to calculate the background mean and sigma.
        Exclude fluxes above and below mean +- 3sigma and recalculate mean
        and sigma. Use thes values as the final estimates of the background flux
        and the expected level of variability in the background.
        
        INPUTS:
        
        :fluxes: (float nxm array) flux time profiles for n energy channels and
            m time points
        :energy_bins: (float nx2 array) energy bins for n energy channels
        
        OUTPUTS:
        
        :means: (float 1xn array) mean values of histogram for n energy channels
        :sigmas: (float 1xn array) sigma of histograms for n energy channels
        
    """
    means = []
    sigmas = []
    for i in range(len(fluxes)):
        strip_flux = fluxes[i]
        for it in range(2): #number of iterations
            #First iteration: Use all fluxes in the background time period to
            #estimate mean and sigma.
            mean, sigma = create_histogram(strip_flux, energy_bins[i],it)
            #exclude values above mean + 3sigma
            highval = mean + 3*sigma
            strip_flux = remove_above(strip_flux,highval)
            #exclude values below mean - 3sigma
            lowval = mean - 3*sigma
            strip_flux = remove_below(strip_flux,lowval)

        if not means:
            means = [mean]
        else:
            means.append(mean)
        if not sigmas:
            sigmas = [sigma]
        else:
            sigmas.append(sigma)

    return means, sigmas



def plot_fluxes(experiment, flux_type, options, fluxes, dates,
                energy_bins, means, sigmas, saveplot):
    """Plot fluxes with time for all of the energy bins on the same plot. The
        estimated mean background levels are plotted as dashed lines.
        Zero values are masked, which is useful when make plots of the
        background and SEP flux separately.
        
        INPUTS:
        
        :experiment: (string) name of experiment
        :flux_type: (string) "integral" or "differential"
        :options: (string array) array of options that can be applied
        :fluxes: (float nxm array) flux time profiles for n energy channels and
            m time points
        :dates: (datetime 1xm array) m time points for flux time profile
        :energy_bins: (float nx2 array) energy bins for n energy channels
        :means: (float 1xn array) mean values of histogram for n energy channels
        :sigmas: (float 1xn array) sigma of histograms for n energy channels
        :saveplot: (bool) True to save plot automatically
        
        OUTPUTS:
        
        Plot to screen and plot saved to file
        
    """
    #All energy channels in specified date range with event start and stop
    #Plot all channels of user specified data
    modifier = ''
    if options[0] != '':
        for opt in options:
            modifier = modifier + '_' + opt

    year = dates[0].year
    month = dates[0].month
    day = dates[0].day
    strdate = str(year) + '_' + str(month) + '_' + str(day)

    figname = strdate + '_Fluxes_' \
             + experiment + modifier + '_' + flux_type + '_' + 'All_Bins'

    if experiment == 'user' and model_name != '':
        figname = strdate + '_Fluxes_' \
             + model_name + modifier + '_' + flux_type + '_' + 'All_Bins'
    fig = plt.figure(figname,figsize=(9,4))
    ax = plt.subplot(111)
    nbins = len(energy_bins)
    for i in range(nbins):
        #Don't want to plot zero values, particularly in background-subtracted plots
        maskfluxes = np.ma.masked_where(fluxes[i] == 0, fluxes[i])

        legend_label = ""
        if energy_bins[i][1] != -1:
            legend_label = str(energy_bins[i][0]) + '-' \
                           + str(energy_bins[i][1]) + ' MeV'
        else:
            legend_label = '>'+ str(energy_bins[i][0]) + ' MeV'
        p = ax.plot_date(dates,maskfluxes,'-',label=legend_label)
        color = p[0].get_color()
        if i==0:
            plt.axhline(means[i],color=color,linestyle=':',\
                        label="Mean Background")
        else:
            plt.axhline(means[i],color=color,linestyle=':')
    colors = ['black','red','blue','green','cyan','magenta']
    if flux_type == "integral":
        plt.ylabel('Integral Flux 1/[cm^2 s sr]')
        plt.title(experiment + " Integral Energy Bins")
        if experiment == 'user' and model_name != '':
            plt.title(model_name + " Integral Energy Bins")
    if flux_type == "differential":
        plt.ylabel('Flux 1/[MeV cm^2 s sr]')
        plt.title(experiment + " Differential Energy Bins")
        if "uncorrected" in options:
            plt.title(experiment + " Uncorrected Differential Energy Bins")
        if experiment == 'user' and model_name != '':
            plt.title(model_name + " Differential Energy Bins")
    plt.xlabel('Date')
    plt.yscale("log")
    chartBox = ax.get_position()
    ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.85,
                     chartBox.height])
    ax.legend(loc='upper center', bbox_to_anchor=(1.17, 1.05))
    if saveplot:
        fig.savefig(plotpath + '/' +figname + '.png')



def derive_background(str_startdate, str_enddate, str_bgstartdate, \
            str_bgenddate, experiment, flux_type, model_name,user_file, \
            showplot, saveplot, options):
    """ Derive the background using fluxes in the time period between
        background start and end dates specified by the user. Derive the
        mean background value along with an expected level of variation (sigma)
        in the background. Make two separate arrays containing 1) background
        flux, 2) background-subtracted SEP fluxes.
        The fluxes will be separate by selecting fluxes above and below
        mean + Nsigma*sigma. The value of Nsigma is specified in
        library/global_vars.py.
        Return the background and background-subtracted SEP flux arrays along
        with a date array. The fluxes and dates will extend from BGStartdate to
        SEPEndDate. The fluxes will be numpy arrays and the dates are a list.
        
        INPUTS:
        
        :str_startdate: (string) starting date of SEP time period in
            "YYYY-MM-DD HH:MM:SS"
        :str_enddate: (string) ending date of SEP time period in
            "YYYY-MM-DD HH:MM:SS"
        :str_bgstartdate: (string) starting date of background time period in
            "YYYY-MM-DD HH:MM:SS"
        :str_bgenddate: (string) ending date of background time period in
            "YYYY-MM-DD HH:MM:SS"
        :experiment: (string) name of experiment
        :flux_type: (string) "integral" or "differential"
        :model_name: (string) contains name of user-input model or
            data set if experiment = "user"
        :user_file: (string) filename containing user input data
        :showplot: (bool) True to print plot to screen
        :saveplot: (bool) True to save plot automatically
        :options: (string array) array of options that can be applied
        
        OUTPUTS:
        
        :bgfluxes: (float nxm array) background fluxes for n energy channels
            and m time points
        :sepfluxes: (float nxm array) background-subtracted SEP fluxes for n
            energy channels and m time points
        :dates: (datetime 1xm array) m time points extending from
            str_bgstartdate to str_enddate containing background flux time period
            and SEP flux time period
        
    """
    if len(str_startdate) == 10: #only YYYY-MM-DD
        str_startdate = str_startdate  + ' 00:00:00'
    if len(str_enddate) == 10: #only YYYY-MM-DD
        str_enddate = str_enddate  + ' 00:00:00'
    startdate = datetime.datetime.strptime(str_startdate, "%Y-%m-%d %H:%M:%S")
    enddate = datetime.datetime.strptime(str_enddate, "%Y-%m-%d %H:%M:%S")

    if len(str_bgstartdate) == 10: #only YYYY-MM-DD
        str_bgstartdate = str_bgstartdate  + ' 00:00:00'
    if len(str_bgenddate) == 10: #only YYYY-MM-DD
        str_bgenddate = str_bgenddate  + ' 00:00:00'
    bgstartdate = datetime.datetime.strptime(str_bgstartdate,
                        "%Y-%m-%d %H:%M:%S")
    bgenddate = datetime.datetime.strptime(str_bgenddate, "%Y-%m-%d %H:%M:%S")

    #CHECKS ON INPUTS
    if (str_startdate == "" or str_enddate == ""):
        sys.exit('You must enter a valid date range. Exiting.')

    if (enddate < startdate):
        sys.exit('End time before start time! Enter a valid date range. '
                'Exiting.')

    if (experiment == "SEPEM" and flux_type == "integral"):
        sys.exit('The SEPEM (RSDv2) data set only provides differential fluxes.'
            ' Please change your FluxType to differential. Exiting.')

    if (experiment == "EPHIN" and flux_type == "integral"):
        sys.exit('The SOHO/EPHIN data set only provides differential fluxes.'
            ' Please change your FluxType to differential. Exiting.')

    sepem_end_date = datetime.datetime(2015,12,31,23,55,00)
    if(experiment == "SEPEM" and (startdate > sepem_end_date or
                   enddate > sepem_end_date)):
        sys.exit('The SEPEM (RSDv2) data set only extends to '
                  + str(sepem_end_date) +
            '. Please change your requested dates. Exiting.')

    if experiment[0:4] == "GOES" and flux_type == "integral":
        sys.exit("Do not perform background subtraction on GOES integral "
                "fluxes. Integral fluxes have already been derived by "
                "applying corrections for cross-contamination and removing "
                "the instrument background levels.")

    if experiment[0:4] == "GOES" and "uncorrected" not in options:
        print("Warning: GOES corrected fluxes have already been derived by "
                "applying corrections for cross-contamination and removing "
                "the instrument and GCR background levels up to channel P6. "
                "Please be sure it makes sense to perform a background "
                "subtraction of this data (e.g. for HEPAD energies). "
                "Otherwise, please add --options uncorrected to perform "
                "background subtracion on GOES uncorrected fluxes. Continuing.")

    if experiment[0:4] == "GOES" and "uncorrected" in options:
        print("Note: Background-subtraction of uncorrected GOES fluxes "
                "does not remove the effects of spurious increases in the low "
                "energy channels due to cross-talk from the high energy "
                "channels, particularly at the onset of well-connected, "
                "intense SEP events. It also does not remove any contamination "
                "due to particles entering the GOES detectors from the sides.")

    #create paths if don't exist
    datasets.check_paths()

    #Check and prepare the data
    filenames1, filenames2, filenames_orien = datasets.check_data(bgstartdate,
                                    enddate, experiment, flux_type, user_file)

    #read in flux files
    if experiment != "user":
        all_dates, all_fluxes, west_detector =datasets.read_in_files(experiment,
                    flux_type, filenames1, filenames2, filenames_orien, options)
    if experiment == "user":
        all_dates, all_fluxes = datasets.read_in_user_files(filenames1)
        west_detector = []

    #Extract the date range specified by the user
    dates, fluxes = datasets.extract_date_range(bgstartdate,enddate,all_dates,
                    all_fluxes)
    if len(dates) <= 1:
        sys.exit("The specified start and end dates were not present in the "
                "specified input file. Exiting.")

    #Get energy bins associated with the fluxes
    energy_bins = datasets.define_energy_bins(experiment, flux_type, \
                                west_detector, options)
    #Remove bad data points (negative fluxes) with linear interpolation in time
    #set bad values to None rather than perform a linear interpolation in time
    dointerp = False
    fluxes = datasets.check_for_bad_data(dates,fluxes,energy_bins,dointerp)

    #Pull out the fluxes in the time period to be used for calculating
    #background
    bg_dates, bg_fluxes = datasets.extract_date_range(bgstartdate, bgenddate,
                        dates, fluxes)
    print('Calculating background with data from ' + str(bgstartdate) + ' to '
        + str(bgenddate) + '.')

    means, sigmas = iterate_background(bg_fluxes, energy_bins)
    bgfluxes, sepfluxes = separate_sep_and_background(fluxes, dates,\
                     means, sigmas)
    #To get into a consistent format as used in operational_sep_quantities.py
    bgfluxes = np.array(bgfluxes)
    sepfluxes = np.array(sepfluxes)

    if showplot or saveplot:
        plot_fluxes('Total_'+experiment, flux_type, options, fluxes, dates, energy_bins,
                    means, sigmas, saveplot)
        plot_fluxes('BackgroundFlux_'+experiment, flux_type, options, bgfluxes,
                    dates, energy_bins, means, sigmas, saveplot)
        plot_fluxes('BGSubSEPFlux_'+experiment, flux_type, options, sepfluxes,
                    dates, energy_bins, means, sigmas, saveplot)

    return bgfluxes, sepfluxes, dates





if __name__ == "__main__":
    #Fluxes > mean background + 3sigma will be considered SEP fluxes
    #The the mean background value will be subtracted away and any
    #flux < mean + 3sigma will be set to zero.
    parser = argparse.ArgumentParser()
    parser.add_argument("--SEPStartDate", type=str, default='',
            help=("Start date of SEP event in YYYY-MM-DD."))
    parser.add_argument("--SEPEndDate", type=str, default='',
            help=("End date of SEP event in YYYY-MM-DD."))
    parser.add_argument("--BGStartDate", type=str, default='',
            help=("Start date to calculate background in YYYY-MM-DD."))
    parser.add_argument("--BGEndDate", type=str, default='',
            help=("End date to calculate background in YYYY-MM-DD."))
    parser.add_argument("--Experiment", type=str, choices=['GOES-08',
            'GOES-10', 'GOES-11', 'GOES-12', 'GOES-13', 'GOES-14', 'GOES-15',
            'SEPEM', 'EPHIN', 'user'],
            default='', help="Enter name of spacecraft or dataset")
    parser.add_argument("--FluxType", type=str, choices=['integral',
            'differential'], default='differential',
            help=("Do you want to use integral or differential fluxes?"
                 " (Default is differential)"))
    parser.add_argument("--ModelName", type=str, default='', help=("If you "
            "chose user for experiment, specify the name of the model or "
            "experiment that you are analyzing (no spaces)."))
    parser.add_argument("--UserFile", type=str, default='tmp.txt', help=("If "
            "you chose user for experiment, specify the filename containing "
            "the fluxes. Specify energy bins and delimeter in code at top. "
            "Default is tmp.txt."))
    parser.add_argument("--options", type=str, default='', help=("You "
            "may specify a series of options as a comma separated list "
            "surrounded by quotations.\n"
            "\"uncorrected\" for GOES uncorrected differential fluxes with "
            "nominal GOES energy bins."))
    parser.add_argument("--showplot",
            help="Flag to display plots", action="store_true")
    parser.add_argument("--saveplot",
            help="Flag to save plots to file", action="store_true")


    args = parser.parse_args()

    str_startdate = args.SEPStartDate
    str_enddate = args.SEPEndDate
    str_bgstartdate = args.BGStartDate
    str_bgenddate = args.BGEndDate
    experiment = args.Experiment
    flux_type = args.FluxType
    model_name = args.ModelName
    user_file = args.UserFile
    showplot = args.showplot
    saveplot = args.saveplot
    options = args.options

    bgfluxes, sepfluxes, dates = derive_background(str_startdate, str_enddate, \
                str_bgstartdate, str_bgenddate, experiment, flux_type, \
                model_name,user_file, showplot, saveplot, options)

    if showplot: plt.show()
