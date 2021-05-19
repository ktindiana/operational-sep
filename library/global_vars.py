import array as arr

#Global variables used by operational_sep_quantities, read_datasets, and
#derive_background.

datapath = '/Users/kwhitman/Documents/Programs/OpSEP/operational-sep/data'
outpath = '/Users/kwhitman/Documents/Programs/OpSEP/operational-sep/output'
plotpath = '/Users/kwhitman/Documents/Programs/OpSEP/operational-sep/plots'
listpath = '/Users/kwhitman/Documents/Programs/OpSEP/operational-sep/lists'
badval = -1 #bad data points will be set to this value; must be negative

email = "kathryn.whitman@nasa.gov"  #Your email for output JSON files

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
version = "2021-01-11" #Enter the version number of your model or data set

#(expect the first (0th) column contains date in YYYY-MM-DD HH:MM:SS format)
#Identify columns containing fluxes you want to analyze
user_col = arr.array('i',[6,7,8,9])
#   REleASE-30 [2,3,4,5]
#   REleASE-60 [6,7,8,9]
#   REleASE-90 [10,11,12,13]
#   SPARX [1,2]
#   SEPMOD [1,2,3,4,5,6,7,8]

#DELIMETER between columns; for whitespace separating columns, use " " or ""
user_delim = ";"

#DEFINE ENERGY BINS associated with user file and columns specified above as:
#   [[Elow1,Ehigh1],[Elow2,Ehigh2],[Elow3,Ehigh3],etc]
#Use -1 in the second edge of the bin to specify integral channel (infinity):
#   [[Elow1,-1],[Elow2,-1],[Elow3,-1],etc]
user_energy_bins = [[4,9],[9,15.8],[15.8,39.6],[20.0,35.5]]
#   SEPEM_H_GOES13 P3 - P7 [[4,9],[12,23],[26,38],[40,73],[100,142],[160,242]]
#   EPREM [[10,-1],[30,-1],[40,-1],[50,-1],[100,-1]]
#   SEPMOD integral [[750,-1],[500,-1],[300,-1],[100,-1],\
#                    [60,-1],[50,-1],[30,-1],[10,-1]]
#   SEPMOD differential [[1000,1000],[750,750],[500,500],[300,300],[100,100],\
#                    [60,60],[50,50],[30,30],[10,10]]
#   SEPMOD IMP8 differential [[72,72],[53,53],[36,36],[26,26],[17,17],\
#                    [8.6,8.6],[5.1,5.1],[2.6,2.6],[1.2,1.2]]
#   SPARX [[10,-1],[60,-1]]
#   ASPECS [[10,-1],[30,-1]]
#   RELeASE [[4,9],[9,15.8],[15.8,39.6],[20.0,35.5]]
################################
