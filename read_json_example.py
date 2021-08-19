from library import ccmc_json_handler as ccmc_json
from library import keys
#from library import global_vars as vars
import argparse
import csv
import datetime
import sys
import os


__version__ = "0.1"
__author__ = "Katie Whitman"
__maintainer__ = "Katie Whitman"
__email__ = "kathryn.whitman@nasa.gov"


def about_read_json_example():
    """This and supporting codes are found in the respository:
            
        https://github.com/ktindiana/operational-sep

        This code gives and example of how to read values from the
        json files produced by operational_sep_quantities.py.
        The json files are read using keys.py and ccmc_json_handler.py
        and the info in global_vars.py (in those subroutines).
        
        The json files contain block organized by energy channel.
        There is one block per energy channel.
        If multiple thresholds were applied to a single energy
        channel, the information for the multiple thresholds
        is stored in fields that are arrays (e.g. event_lengths,
        fluences, fluence_spectra, threshold_crossings, etc)
        
        ccmc_json_handler.py provides three different way to extract
        information from a specific block in the json file:
        
        1. return_json_value_by_energy(injson, value, energy_channel={}, index=0)
        2. return_json_value_by_index(injson, value, channel_index=0, index=0)
        3. return_json_value_by_threshold(injson, value, energy_channel={}, threshold=0)
        
        Method 1 requires the user to provide which energy channel and
        the appropriate block will be found that way. If the desired value
        is inside an array, index indicates which array element to extract.
        
        Method 2 identifies the desired block by index. If you know the
        energy channel you want is stored in the 2nd block, then
        channel_index = 1
        
        If the desired value
        is inside an array, index indicates which array element to extract.
        
        Method 3 uses an energy channel and threshold combination to uniquely
        identify the desired value.
        
        For values that are not stored in arrays, index is optional.
        
        For any values that are in the "top" level of the json file
        above the "forecasts" or "observations" field, such
        as "short_name" or "options", all of the
        methods will work equally well to extract that info.
    """



def print_json_values(filename):
    """ Print to screen all of the values found in a json file
        using 3 different methods.
    """
    
    injson = ccmc_json.read_in_json(filename)
    id_all = keys.id_all #array containing every possible value in
                         #the json file
    
    print("==============\n")
    print("Method 1: Identify the desired energy channel block by "
            "specifying the energy channel:")
    energy_channel = {'min':10, 'max':-1, 'units':"MeV"}
    print("energy_channel: " + str(energy_channel))
    for id_value in id_all:
        val = ccmc_json.return_json_value_by_energy(injson, id_value, energy_channel)
        
        print("ID: " + id_value + ", Value: " + str(val))
        
    print("==============\n")
    print("Method 2: Identify the desired energy channel block by "
            "specifying the index of the block:")
    channel_index = 0 #operational_sep_quantities.py always put >10MeV first
    
    for id_value in id_all:
        channel_index = 0 #operational_sep_quantities.py always put >10MeV first
        val = ccmc_json.return_json_value_by_index(injson, id_value, channel_index)
        
        print("ID: " + id_value + ", Value: " + str(val))
        
    
    print("================\n")
    print("Method 2: Identify the desired energy channel block by "
            "specifying the energy_channel and threshold:")
    energy_channel = {'min':10, 'max':-1, 'units':"MeV"}
    threshold = 10
    print("energy_channel: " + str(energy_channel) + ", Threshold: "
            + str(threshold))
            
    for id_value in id_all:
        val = ccmc_json.return_json_value_by_threshold(injson, id_value,\
                    energy_channel, threshold)
        
        print("ID: " + id_value + ", Value: " + str(val))
    
    ####ADDITIONAL EXAMPLES NOT PRINTED TO SCREEN#######
    ####LOOK AT keys.py TO GET THE INDIVIDUAL VALUE IDENTIFIERS###
    #To extract the onset peak flux for the 10MeV channel:
    energy_channel = {'min':10, 'max':-1, 'units':"MeV"}
    onset = ccmc_json.return_json_value_by_energy(injson, keys.id_peak_intensity,\
                    energy_channel)
                    
    #To extract the start time for the 10MeV channel with 10 pfu
    #threshold applied:
    threshold = 10
    start_time = ccmc_json.return_json_value_by_threshold(injson, \
                    keys.id_event_length_start_time,\
                    energy_channel,threshold)
                    
    #You may also directly enter the string identifier.
    #I saved them in a variable in keys.py in case any of them
    #ever changed. Then I would only have to change a single variable
    #rather than every instance of them throughout all the codes.
    start_time = ccmc_json.return_json_value_by_threshold(injson, \
                    "event-length-start-time",\
                    energy_channel,threshold)
                    
    
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--Filename", type=str, default='', \
            help=("Full path and name of json file to read in."))
    
    args = parser.parse_args()
    filename = args.Filename
    
    print_json_values(filename)
    

