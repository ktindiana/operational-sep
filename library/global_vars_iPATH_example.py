import array as arr

#Global variables used by operational_sep_quantities, read_datasets, and
#derive_background.

datapath = 'data'
outpath = 'output'
plotpath = 'plots'
listpath = 'lists'
badval = -1 #bad data points will be set to this value; must be negative
endfac = 0.85 #factor multiplied by flux threshold to determine end of event
errval = "Value Not Found"  #alternative value to None to indicate value not present

###FOR BACKGROUND SUBTRACTION###
#derive_background.py calculates the mean background plus an expected level
#of variation (sigma).
#SEP flux is flux > mean + nsigma*sigma
#Background flux is flux <= mean + nsigma*sigma
#The value of nsigma affects how the event onset is captured because it
#determines the flux level that is considered background at the beginning of the
#event
nsigma = 2.0
#################################


########FOR USER DATA SETS#######
version = "" #Enter the version number of your model or data set

#####(expect the first (0th) column contains date in YYYY-MM-DD HH:MM:SS format)
#Identify columns containing fluxes you want to analyze
user_col = arr.array('i',[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25])


#####DELIMETER between columns; for whitespace separating columns, use ""
user_delim = ","

#####DEFINE ENERGY BINS associated with user file and columns specified above as:
#   [[Elow1,Ehigh1],[Elow2,Ehigh2],[Elow3,Ehigh3],etc]
#Use -1 in the second edge of the bin to specify integral channel (infinity):
#   [[Elow1,-1],[Elow2,-1],[Elow3,-1],etc]
user_energy_bins = [[0.1,0.1],[0.14677993,0.14677993],[0.21544347,0.21544347],\
                    [0.31622777,0.31622777],[0.46415888,0.46415888],\
                    [0.68129207,0.68129207],[1,1],[1.4677993,1.4677993],\
                    [2.1544347,2.1544347],[3.1622777,3.1622777],[4.6415888,4.6415888],\
                    [6.8129207,6.8129207],[10,10],[14.677993,14.677993],\
                    [21.544347,21.544347],[31.622777,31.622777],[46.415888,46.415888],\
                    [68.129207,68.129207],[100,100],[146.77993,146.77993],\
                    [215.44347,215.44347],[316.22777,316.22777],[464.15888,464.15888],\
                    [681.29207,681.29207],[1000,1000]]


#####UNITS
#Set the units associated with your data
#Make sure that fluence units = flux units*s*sr
#If your data set is in units other than pfu or
#MeV^-1*cm^-2*s^-1*sr^-1, the two operational thresholds
#will not be equivalent to >10 MeV, 10 pfu and >100 MeV, 1 pfu
#--energy
energy_units = "MeV"
#--integral
flux_units_integral = "pfu" #pfu = cm^-2*s^-1*sr^-1
fluence_units_integral = "cm^-2"
#--differential
flux_units_differential = "MeV^-1*cm^-2*s^-1*sr^-1"
fluence_units_differential = "MeV^-1*cm^-2"
################################


def about_global_vars():
    """ About global_vars.py
    
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
