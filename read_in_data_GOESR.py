#This is a program that shows how to read in data using read_datasets.py
#The goal of this program is to get the data into arrays that are easily 
#used for other purposes.
#Also some plotting code
#CODE EXPECTS THERE TO BE A data/ DIRECTORY
#K. Whitman 2022-02-02


from library import read_datasets_GOESR as datasets
import numpy as np
import matplotlib.pyplot as plt
import datetime
import argparse
import math
import sys



#####UNITS#####
#If a user data set is read in, the units in global_vars.py will
#supercede these units.
energy_units = "MeV"
flux_units_integral = "pfu"
fluence_units_integral = "cm^-2"
flux_units_differential = "MeV^-1*cm^-2*s^-1*sr^-1"
fluence_units_differential = "MeV^-1*cm^-2"


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



def make_plots(experiment, flux_type, dates, fluxes, energy_bins,
            showplot, saveplot):
    """ Plot the fluxes for all energy bins in the specified
        time range.
    """

    #All energy channels in specified date range with event start and stop
    print("Generating figure of fluxes in original energy bins. Any bad data "
          "points were interpolated. Lines indicate event start and stop for "
          "thresholds.")
    #Plot all channels of user specified data
    figname = 'Fluxes' \
            + '_' + experiment + '_' + flux_type\
            + '_' + 'All_Bins'
    if experiment == 'user' and model_name != '':
        figname = 'Fluxes' \
                + '_' + model_name + '_' + flux_type\
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
        plt.title(experiment + '\n'\
                    + "Integral Energy Bins")

    if flux_type == "differential":
        plt.ylabel('Flux [' + flux_units_differential + ']')
        plt.title(experiment + '\n'\
                    + "Differential Energy Bins")

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


    if showplot: plt.show()



def read_in_flux_files(experiment, flux_type, str_startdate,
        str_enddate, nointerp, showplot, saveplot):
    """ Read in the appropriate data or user files. Performs
        background subtraction, if requested. Trims to dates
        between start time and end time. Interpolates bad
        points with linear interpolation in time.
        
        INPUTS:
        
        :experiment: (string)
        :flux_type: (string) - integral, differential
        :startdate: (datetime) - start date of time period entered by user
        :enddate: (datetime) - end date of time period entered by user
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
    
    if (str_startdate == "" or str_enddate == ""):
        sys.exit('You must enter a valid date range. Exiting.')

    startdate = str_to_datetime(str_startdate)
    enddate = str_to_datetime(str_enddate)
    
    datasets.check_paths()
    
    #Check if data is present on your computer, if not, download
    filenames1 = datasets.check_data(startdate, enddate, experiment, flux_type)
    print(filenames1)
    #read in daily flux files
    dates, fluxes = datasets.read_in_files(experiment, flux_type, filenames1)

    #Define energy bins
    energy_bins = datasets.define_energy_bins(experiment, flux_type)
        
    #Set bad data points (negative flux or None) to None
    dointerp = not nointerp
    fluxes = datasets.check_for_bad_data(dates,fluxes,energy_bins,dointerp)
 

    if len(dates) <= 1:
        sys.exit("The specified start and end dates were not present in the "
                "specified input file. Exiting.")
                
    if showplot or saveplot:
        make_plots(experiment, flux_type, dates, fluxes, energy_bins,
            showplot, saveplot)
    
    return dates, fluxes, energy_bins
      



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
    parser.add_argument("--Experiment", type=str,
            choices=['GOES-16', 'GOES-17'],
            default='GOES-16', help="Enter name of spacecraft or dataset")
    parser.add_argument("--FluxType", type=str, choices=['integral',
            'differential'], default='',
            help=("Do you want to use integral or differential fluxes?"))
    parser.add_argument("--NoInterp",
            help=("Do not fill in negative or missing fluxes via "
                    "linear interpolation in time. Set as None values "
                    "instead."), action="store_true")
    parser.add_argument("--showplot",
            help="Flag to display plots", action="store_true")
    parser.add_argument("--saveplot",
            help="Flag to save plots to file", action="store_true")


    args = parser.parse_args()

    str_startdate = args.StartDate
    str_enddate = args.EndDate
    experiment = args.Experiment
    flux_type = args.FluxType
    showplot = args.showplot
    saveplot = args.saveplot
    nointerp = args.NoInterp

    dates,fluxes,energy_bins = read_in_flux_files(experiment, flux_type, str_startdate,
        str_enddate, nointerp, showplot, saveplot)

