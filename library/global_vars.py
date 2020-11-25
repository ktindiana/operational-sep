import array as arr

#Global variables used by operational_sep_quantities, read_datasets, and
#derive_background.

datapath = 'data'
outpath = 'output'
plotpath = 'plots'
badval = -1 #bad data points will be set to this value; must be negative

email = "kathryn.whitman@nasa.gov"  #Your email for output JSON files

###FOR BACKGROUND SUBTRACTION###
#derive_background.py calculates the mean background plus an expected level
#of variation (sigma).
#SEP flux is flux > mean + nsigma*sigma
#Background flux is flux <= mean + nsigma*sigma
#The value of nsigma affects how the event onset is captured
nsigma = 1.5
#################################


########FOR USER DATA SETS#######
version = "vSHINE2020" #Enter the version number of your model or data set

#(expect the first (0th) column contains date in YYYY-MM-DD HH:MM:SS format)
#Identify columns containing fluxes you want to analyze
user_col = arr.array('i',[1,2])

#DELIMETER between columns; for whitespace separating columns, use " " or ""
user_delim = ","

#DEFINE ENERGY BINS associated with user file and columns specified above as:
#   [[Elow1,Ehigh1],[Elow2,Ehigh2],[Elow3,Ehigh3],etc]
#Use -1 in the second edge of the bin to specify integral channel (infinity):
#   [[Elow1,-1],[Elow2,-1],[Elow3,-1],etc]
user_energy_bins = [[10,-1],[60,-1]]
#SEPEM_H_GOES13 P3 - P7 [[4,9],[12,23],[26,38],[40,73],[100,142],[160,242]]
#EPREM [[10,-1],[30,-1],[40,-1],[50,-1],[100,-1]]
#SEPMOD [[750,-1],[500,-1],[300,-1],[100,-1],\
#                    [60,-1],[50,-1],[30,-1],[10,-1]]
#SPARX [[10,-1],[60,-1]]
################################
