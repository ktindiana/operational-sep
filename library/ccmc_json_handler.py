import json
import calendar
import datetime
from datetime import timedelta
import copy
import zulu
from library import global_vars as vars
from library import keys
import os

__version__ = "1.1"
__author__ = "Katie Whitman"
__maintainer__ = "Katie Whitman"
__email__ = "kathryn.whitman@nasa.gov"


#2021-01-12 changes in v0.2: Added fluence to json file. Added threshold
#   crossing entry for end of event.
#2021-02-11 changes in 0.3: Changed logic so that if integral flux was
#   derived from differential flux, and the flux did not cross the
#   threshold applied to the integral channel, then all_clear_boolean = True.
#   Previously, the code checked if the energy bin existed and set
#   all_clear_boolean to None if not. This produced the wrong value in the
#   case where differential channels were converted to integral.
#2021-02-24, changes in v0.4: Added subroutines to read in json file and
#   convert zulu time to datetime.
#   Modified fill_json to account for change in threshold fluences array in
#   v2.5 of operational_sep_quantities.py.
#2021-05-18, changes in 0.5: Realized that there were some differences
#   between the json files produced here and the CCMC format. CCMC defines
#   the fluences field as an array allowing multiple fluences, so the fluences
#   field was changed to an array here.
#   CCMC also has "event_lengths" as an array that allows the specification
#   of multiple start and end times. Changed the field name from "event_length"
#   to "event_lengths" and made it an array. Only one entry to
#   "event_lengths" is created by this program.
#2021-07-20, changes in 1.0: Went up an integer in version number because
#   making major changes to exactly match CCMC JSON format for the SEP
#   Scoreboard. Version 1.0 of this file works with v3.0 of
#   operational_sep_quantities.py, which has been modified for the same
#   purpose.
#   Removing contacts field from json due to NASA privacy rules.
#   Changed zulu time to keep seconds.
#   If multiple thresholds applied to a single energy channel,
#   values will be saved in one energy channel block as array
#   entries (previously writing separate blocks for unique
#   energy-threshold combinations.
#   Removing any attempt to differentiate between values not
#   filled out because they are not supported by the model
#   or values not filled out because there was no forecast.
#   Will have to apply that kind of information downstream.
#   JSON fields that have not been filled in are removed
#   from the final output json file.
#2021-08-18, changes in 1.1: Added subroutines from
#   validation_json_handler.py in the validation project
#   to access entries in the json file in a variety of ways.
#   This code, coupled with keys.py, will allow users
#   to access any values in the json files produced by
#   operational_sep_quantities.py

version = vars.version

def about_ccmc_json_handler():
    """ ABOUT ccmc_json_handler.py
    
        This module places all derived information into json and supporting
        files that follow the format specified for CCMC's SEP Scoreboard.
        https://ccmc.gsfc.nasa.gov/challenges/sep.php#format
        
        Reads in templates and fills in SEP event quantities organized
        by one block for each energy channel. If multiple thresholds
        were applied to the same energy channel, the derived information
        will be stored in a single block.
        
        Template for observations: observations_template.json
        
        Template for model output: model_template.json
        
        ADDITIONALLY:
        
        ccmc_jason_handler.py combined with library/keys.py
        provides all of the information needed to read the JSON
        files in the CCMC SEP Scoreboard format or produced
        by operational_sep_quantities.py.
        
        The user provides a unique identifier as listed in keys.py
        or the json visual scheme document on the CCMC webpage:
        https://ccmc.gsfc.nasa.gov/challenges/sep.php#format
        The CCMC identifiers are followed as closely as possible, with some
        additions specific to the operational_sep_quantities.py code and
        this SEP model validation effort.
        
        In some cases, the user must also provide an index value or
        a threshold value to uniquely identify values that are stored
        in arrays in the json file, These include 'fluences', 'fluence_spectra',
        'events_lengths', 'threshold_crossings','probabilities'
        
        In this code, a key_chain will be requested from keys.py. The
        key_chain contains the unique location of the desired value, but
        in a list of strings and integers (if an index is needed).
        validation_jason_handler.py extracts the unique value by
        extracting each sub-dictionary in an iterative manner until
        the unique value is accessed.
        
        For example, if key_chain = ['event_lengths',1,'start_time'], the
        value will be extracted from the full json dictionary (dict) as a
        loop over the key_chain, pulling out a subdictionary each time until
        the unique value is found, effectively as follows:
        
        .. code-block:: python
        
            subdict1 = dict['sep_forecast_submission']['forecasts']['event_lengths']
            subdict2 = subdict1[1]
            desired_val = subdict2['start_time']
            
        
    """


def read_in_json_template(type):
    """Read in appropriate json file template for model or observations.
        
        INPUTS:
        
        :type: (string) = "model" or "observations"
        
        OUTPUTS:
        
        :template: (dict) - appropiate template is read in and returned
        
    """
    if type != "model" and type != "observations":
        sys.exit("json_handler: read_in_template: type may be \"model\" "
                "or \"observations\". You entered " + str(type))

    if type == "model":
        with open('library/model_template.json') as f:
            template=json.load(f)

    if type == "observations":
        with open('library/observations_template.json') as f:
            template=json.load(f)

    return template


def make_ccmc_zulu_time(dt):
    """ Make a datetime string in the format YYYY-MM-DDTHH:MM:SSZ
        
        INPUTS:
        
        :dt: (datetime)
        
        OUTPUTS:
        
        :zuludate: (string) in the format YYYY-MM-DDTHH:MM:SSZ
    
    """
    if dt == None:
        return None
    if dt == 0:
        return 0

    zdt = zulu.create(dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second)
    stzdt = str(zdt)
    stzdt = stzdt.split('+00:00')
    zuludate = stzdt[0] + "Z"
    return zuludate
 
 
def zulu_to_time(zt):
    """ Convert Zulu time to datetime
    
        INPUTS:
        
        :zt: (string) - date in the format "YYYY-MM-DDTHH:MM:SSZ"
        
        OUTPUTS:
        
        :dt: (datetime)
        
    """
    #Time e.g. 2014-01-08T05:05:00Z or 2012-07-12T22:25Z
    if zt == '':
        return ''
    if zt == None:
        return None
    if zt == 0:
        return 0
  
    strzt = zt.split('T')
    strzt[1] = strzt[1].strip('Z')
    n = strzt[1].split(':')
    stdt = strzt[0] + ' ' + strzt[1]

    if len(n) == 2:
        dt = datetime.datetime.strptime(stdt, '%Y-%m-%d %H:%M')
    if len(n) == 3:
        dt = datetime.datetime.strptime(stdt, '%Y-%m-%d %H:%M:%S')
    return dt


def find_energy_bin(lowedge, energy_bins):
    """ Identify the energy bin and return low and high edges.
    
        INPUTS:
        
        :lowedge: (float) - low edge of an energy bin
        :energy_bins: (float 2xn array) - energy bins
        
        OUTPUTS:
        
        :bin: (float 2x1 array) - single energy bin
    
    """
    bin = []
    for i in range(len(energy_bins)):
        if lowedge == energy_bins[i][0]:
            bin = energy_bins[i]

    if not bin:
        print("ccmc_json_handler: find_energy_bin could not identify "
                "requested bin.")

    return bin


def id_unique_energy_channels(energy_thresholds):
    """ There may be multiple energy channel - flux threshold
        combinations. Identify the number of unique energy
        channels, e.g.
        
            * >10 MeV, 10 pfu (operations) and
            * >10 MeV, 0.001 pfu (SEPMOD testing)
        
        INPUTS:
        
        :energy_thresholds: (float 1xn array)- array containing
            all the energy channels for which a threshold
            has been applied
        
        OUTPUTS:
        
        :len(unique): (integer) - number of unique energy channels
        :unique: (float 1xm array) - the m unique energy channels to
            which a threshold was applied
        
    """
    unique = []
    for energy in energy_thresholds:
        if energy not in unique:
            unique.append(energy)
    
    return len(unique), unique

###############FILL AND WRITE JSONS##############
def fill_json(template, issue_time, experiment, flux_type, energy_bins,
                model_name, spase_id, startdate, enddate, options,
                energy_thresholds, flux_thresholds, crossing_time,
                onset_peak, onset_date, peak_flux, peak_time, rise_time,
                event_end_time, duration, all_threshold_fluences,
                diff_thresh, all_fluence,
                umasep, umasep_times, umasep_fluxes, profile_filenames,
                energy_units, flux_units_integral, fluence_units_integral,
                flux_units_differential, fluence_units_differential):
    """ Add all the appropriate values to the json template for model or
        observations.
        
        The inputs here are mainly the same as the outputs in
        operational_sep_quantities.append_differential_thresholds()
    """
    #For now, assume a user-input file is model output
    #This is not generic and should be modified in the future
    if experiment == "user":
        key = keys.model_main
        type_key = keys.model_type
        win_key = keys.model_win
        template[key]['model']['short_name'] = model_name
        template[key]['model']['flux_type'] = flux_type
        if spase_id == "":
            template[key]['model']['spase_id'] = 'spase://CCMC/SimulationModel/' \
                                            + model_name + "/" + version
        else:
            template[key]['model']['spase_id'] = spase_id
                        
    else:
        key = keys.obs_main
        type_key = keys.obs_type
        win_key = keys.obs_win
        template[key][keys.obs_exp]['short_name'] = experiment
        template[key][keys.obs_exp]['flux_type'] = flux_type
        if spase_id != "":
            template[key]['observatory']['spase_id'] = spase_id

    template[key]['options'] = options
    template[key]['issue_time'] = issue_time


    #####FILL THRESHOLDS#############
    nthresh = len(energy_thresholds)
    
    #if template doesn't contain enough energy channel entries, add more
    #Only one block per energy channel
    nunique, energy_unique = id_unique_energy_channels(energy_thresholds)
    print("fill_json: Number of unique energy channels that have had thresholds "
        "applied " + str(nunique))
    nent = len(template[key][type_key])
    if nent < nunique:
        for j in range(nent,nunique):
            add_entry = copy.deepcopy(template[key][type_key][nent-1])
            template[key][type_key].append(add_entry)
            template[key][type_key][j]['energy_channel']['min'] = energy_unique[j]
            
    #Fill in the unique energy channels
    for i in range(nthresh):
        if not diff_thresh[i] and flux_type == "differential":
            bin = [energy_thresholds[i], -1]
        else:
            bin = find_energy_bin(energy_thresholds[i], energy_bins)
        
        energy_dict = {"min":bin[0], "max": bin[1], "units":energy_units}
            
        for j in range(nunique):
            template_dict = template[key][type_key][j]['energy_channel']
            if template_dict['min'] == energy_dict['min']:
                template[key][type_key][j]['energy_channel'].update(energy_dict)
            

    #Fill in the info for each energy channel and threshold
    tidx = 0 #index of block in json template
    for i in range(nthresh):
        if not diff_thresh[i]: #integral channel
            flux_units = flux_units_integral #"pfu"
            fluence_units = fluence_units_integral #'cm^-2'
        if diff_thresh[i]: #differential channel
            flux_units = flux_units_differential #"[MeV/n cm^2 s sr]^(-1)"
            fluence_units = fluence_units_differential #'MeV^-1*cm^-2'
            
        if not diff_thresh[i] and flux_type == "differential":
            bin = [energy_thresholds[i], -1]
        else:
            bin = find_energy_bin(energy_thresholds[i], energy_bins)
        
        
        #ENERGY CHANNEL
        energy_dict = {"min":bin[0], "max": bin[1], "units":energy_units}
        #Identify the correct block in the json template
        for j in range(nunique):
            if energy_dict == template[key][type_key][j]['energy_channel']:
                tidx = j
        
        zst = make_ccmc_zulu_time(startdate)
        zend = make_ccmc_zulu_time(enddate)
        template[key][type_key][tidx]['species'] = "proton"
        template[key][type_key][tidx]['location'] = "earth"
        template[key][type_key][tidx][win_key]['start_time'] = zst
        template[key][type_key][tidx][win_key]['end_time'] = zend

        #Threshold WAS crossed
        if crossing_time[i] != 0:
            zodate = make_ccmc_zulu_time(onset_date[i])
            zpdate = make_ccmc_zulu_time(peak_time[i])
            zct = make_ccmc_zulu_time(crossing_time[i])
            zeet = make_ccmc_zulu_time(event_end_time[i])
            all_clear = False
            
            #Onset Peak Flux
            onset_dict = {"intensity":onset_peak[i],"time": zodate,
                            "units":flux_units}
            template[key][type_key][tidx]['peak_intensity'].update(onset_dict)
            
            #Maximum Flux
            max_dict = {"intensity":peak_flux[i],"time": zpdate,
                            "units":flux_units}
            template[key][type_key][tidx]['peak_intensity_max'].update(max_dict)
            
            #Start and End Times
            length_dict = { "start_time": zct,  "end_time": zeet,
                            "threshold": flux_thresholds[i],
                            "threshold_units": flux_units }
            if template[key][type_key][tidx]['event_lengths'][0]['start_time']\
                == "":
                template[key][type_key][tidx]['event_lengths'][0].update(length_dict)
            else:
                template[key][type_key][tidx]['event_lengths'].append(length_dict)
            
            #Fluence for only the specific energy channel
            #Fluence order same as event_lengths
            fluence_dict = {"fluence": all_threshold_fluences[i],
                            "units": fluence_units}
            if template[key][type_key][tidx]['fluences'][0]['fluence'] == "":
                template[key][type_key][tidx]['fluences'][0].update(fluence_dict)
            else:
                template[key][type_key][tidx]['fluences'].append(fluence_dict)
            
            
            #Fluence spectrum
            spectrum = []
            for kk in range(len(energy_bins)):
                flentry = {"energy_min": energy_bins[kk][0],
                            "energy_max": energy_bins[kk][1],
                            "fluence": all_fluence[i][kk]}
                spectrum.append(flentry)
            
            fl_spec_dict = {"start_time": zct, "end_time": zeet,
                   "threshold_start":flux_thresholds[i],
                   "threshold_end":flux_thresholds[i]*vars.endfac,
                   "threshold_units":flux_units,
                   "fluence_units": fluence_units,
                   "fluence_spectrum":spectrum }
            if template[key][type_key][tidx]['fluence_spectra'][0]['start_time']\
                == "":
                template[key][type_key][tidx]['fluence_spectra'][0].update(fl_spec_dict)
            else:
                template[key][type_key][tidx]['fluence_spectra'].append(fl_spec_dict)
            

            #Threshold crossings
            cross_dict = { "crossing_time": zct,
                            "threshold": flux_thresholds[i],
                            "threshold_units": flux_units }
            if template[key][type_key][tidx]['threshold_crossings'][0]['crossing_time']\
                == "":
                template[key][type_key][tidx]['threshold_crossings'][0].update(cross_dict)
            else:
                template[key][type_key][tidx]['threshold_crossings'].append(cross_dict)
            #template[key][type_key][i]['threshold_crossings'][0]['crossing_time']\
            #                    = zct
            #template[key][type_key][i]['threshold_crossings'][0]['threshold']\
            #                    = flux_thresholds[i]
 
            #All Clear
            #Fill in All Clear if not already filled in.
            #If there was already a forecast made for a different threshold,
            #will not replace it
            if template[key][type_key][tidx]['all_clear']['all_clear_boolean']\
                == "":
                template[key][type_key][tidx]['all_clear']['all_clear_boolean'] \
                                = all_clear
                template[key][type_key][tidx]['all_clear']['threshold'] \
                                = flux_thresholds[i]
            
            #SEP Flux Time Profile
            template[key][type_key][tidx]['sep_profile'] = profile_filenames[i]


        #Threshold was NOT crossed OR doesn't exist in data set
        if crossing_time[i] == 0:
            zodate = ""
            zpdate = make_ccmc_zulu_time(peak_time[i])
            zct = ""
            zeet = ""
            all_clear = True
            template[key][type_key][tidx]['sep_profile'] = profile_filenames[i]
            
            #The >10 and >100 MeV channels are special in that they
            #are used operationally. The All Clear is thus derived
            #ONLY from the operational thresholds.
            #>10 MeV, 10 pfu and >100 MeV, 1 pfu
            #The All Clear for any other energy channel is defined
            #by whether the user threshold was crossed or not
            if energy_thresholds[i] == 10:
                if peak_flux[i] >= 10: all_clear = False
            if energy_thresholds[i] == 100:
                if peak_flux[i] >= 1: all_clear = False
                    
            max_dict = {"intensity":peak_flux[i],"time": zpdate,
                        "units":flux_units}
            template[key][type_key][tidx]['peak_intensity_max'].update(max_dict)
            
            #Fill in All Clear if not already filled in.
            #If there was already a forecast made for a different threshold,
            #will not replace it
            if template[key][type_key][tidx]['all_clear']['all_clear_boolean']\
                == "":
                template[key][type_key][tidx]['all_clear']['all_clear_boolean'] \
                                = all_clear
                template[key][type_key][tidx]['all_clear']['threshold'] \
                                = flux_thresholds[i]

    return template


def clean_json(template, experiment):
    """ Remove any fields that didn't get filled in and
        were left as empty strings.
    """
    if experiment == "user":
        key = keys.model_main
        type_key = keys.model_type
        win_key = keys.model_win
        
        spase_id = template[key]['model']['spase_id']
        if spase_id == "":
            template[key]['model'].pop('spase_id', None)
       
                        
    else:
        key = keys.obs_main
        type_key = keys.obs_type
        win_key = keys.obs_win

        spase_id = template[key]['observatory']['spase_id']
        if spase_id != "":
            template[key]['observatory'].pop('spase_id', None)

    nent = len(template[key][type_key])
    for i in range(nent):
        #Onset Peak Flux
        if template[key][type_key][i]['peak_intensity']['intensity'] == "":
                    template[key][type_key][i].pop('peak_intensity', None)
        
        #Maximum Flux
        if template[key][type_key][i]['peak_intensity_max']['intensity'] == "":
                    template[key][type_key][i].pop('peak_intensity_max', None)
        
        #Start and End Times, Fluence
        nev = len(template[key][type_key][i]['event_lengths'])
        for j in range(nev-1,-1,-1):
            if template[key][type_key][i]['event_lengths'][j]['start_time']\
                == "":
                template[key][type_key][i]['event_lengths'].pop(j)
                            
            if template[key][type_key][i]['fluences'][j]['fluence'] == "":
                template[key][type_key][i]['fluences'].pop(j)
        
        if template[key][type_key][i]['event_lengths'] == []:
            template[key][type_key][i].pop('event_lengths')
        if template[key][type_key][i]['fluences'] == []:
            template[key][type_key][i].pop('fluences')
        
        
        #Fluence spectrum
        nev = len(template[key][type_key][i]['fluence_spectra'])
        for j in range(nev-1,-1,-1):
            if template[key][type_key][i]['fluence_spectra'][j]['start_time']\
                == "":
                template[key][type_key][i]['fluence_spectra'].pop(j)
        if template[key][type_key][i]['fluence_spectra'] == []:
            template[key][type_key][i].pop('fluence_spectra')
        
        
        #Threshold crossings
        nev = len(template[key][type_key][i]['threshold_crossings'])
        for j in range(nev-1,-1,-1):
            if template[key][type_key][i]['threshold_crossings'][j]['crossing_time'] == "":
                template[key][type_key][i]['threshold_crossings'].pop(j)
        if template[key][type_key][i]['threshold_crossings'] == []:
            template[key][type_key][i].pop('threshold_crossings')
            
        #All Clear
        if template[key][type_key][i]['all_clear']['all_clear_boolean'] \
                            == "":
            template[key][type_key][i].pop('all_clear')
        
        #SEP Flux Time Profile
        if template[key][type_key][i]['sep_profile'] == "":
            template[key][type_key][i].pop('sep_profile')
            
    return template


def write_json(template, filename):
    """Write json template to json file. """
    with open(filename, "w") as outfile:
        json.dump(template, outfile)

    if not os.path.isfile(filename):
        return False

    print("Wrote SEP values to json file --> " + filename)
    return True


###########READ JSONS####################

def read_in_json(filename):
    """Read in json file """
    if not os.path.isfile(filename):
        sys.exit("validation_json_handler: could not read in file " \
                + filename + ". Exiting.")

    with open(filename) as f:
        info=json.load(f)

    return info

def return_main_key(injson):
    """ Return the highest level key in the json file, typically
        'sep_forecast_submission' or 'sep_observation_submission'.
        Possible values assigned in keys.py.
        
        INPUTS:
        
        :injson: (dictionary) - a complete json file for a model or
            observation
            
        OUTPUTS:
        
        :main_key: (string) - the highest level key in the json file
        
    """
    #check if json file is empty
    if not injson:
        print("return_type_key: JSON file is empty.")
        return vars.errval
    
    key_list = list(injson.keys())
    main_key = key_list[0] #only one at top level
    return main_key


def return_type_key(injson):
    """ Return the key associated with type (model or observations)
        where the forecasts or observed information are stored in
        the json file, typically 'forecasts' or 'observations'.
        Possible values assigned in keys.py.
        
        INPUTS:
        
        :injson: (dictionary) - a complete json file for a model or
            observation
            
        OUTPUTS:
        
        :type_key: (string) - a second level key in the json file which
            should match either the obs_type or model_type keys in keys.py
        
    """
    main_key = return_main_key(injson)
    
    key_list = injson[main_key].keys()
    if keys.model_type in key_list:
        type_key = keys.model_type
    elif keys.obs_type in key_list:
        type_key = keys.obs_type
    else:
        print("return_type_key: Could not identify type key in "
            "set of keys. Should be " + keys.model_type + " or "
            + keys.obs_type + ". " + "keys: " + str(key_list))
        return vars.errval
    
    return type_key
    

def return_nforecasts(injson):
    """ Return the number of forecasts or observations
        in the file. e.g.
        
        len(json['sep_forecast_submission']['forecasts'])
        
        INPUTS:
        
        :injson: (dictionary) - a complete json file for a model or
            observation
            
        OUTPUTS:
        
        :n: (integer) - the number of forecast or observation blocks in
            the json file
            
    """
    main_key = return_main_key(injson)
    type_key = return_type_key(injson)
    if main_key == vars.errval or type_key == vars.errval:
        return vars.errval
    
    arr = injson[main_key][type_key]
    n = len(arr)
    
    return n


def switch_model_to_obs_keys(key_chain):
    """ keys.py provides the location of each entry in the
        json files, but uses default identifiers that correspond
        to the ones used by CCMC for model forecasts.
        
        operational_sep_quantities.py outputs similar jsons for
        measurements, but since they are observations and not
        forecasts, some of the fields have been relabeled.
        
        If looking at an observations json, fix the keys
        to reflect the correct names of the fields.
        
        INPUTS:
        
        :key_chain: (array, list) list of strings and integers
        
        OUTPUTS:
        
        :key_chain: (array, list) list of strings and integers
        
    """
    for k in range(len(key_chain)):
        if key_chain[k] == keys.model_exp:
            key_chain[k] = keys.obs_exp
        if key_chain[k] == keys.model_type:
            key_chain[k] = keys.obs_type
        if key_chain[k] == keys.model_win:
            key_chain[k] = keys.obs_win
            
    return key_chain



def return_json_value_by_energy(injson, value, energy_channel={}, index=0):
    """ Return the value of a specific entry in the json file listed in
        the fields under 'observations' or 'forecasts'.
        
        Identify the desired block under 'observations' or 'forecasts'
        by matching the 'energy_channel' entry. The array index
        for the desired block will be discovered.
        
        matches input energy_channel dictionary to either
        injson['sep_forecast_submission']['forecasts'][i]['energy_channel']
        injson['sep_observation_submission']['observations'][i]['energy_channel']
        
        INPUTS:
    
        :injson: (dictionary) is the dictionary created by reading in a CCMC
            SEP Scoreboard JSON file for a single SEP event or
            forecast period
            
        :value: (string) indicates value desired and may be any of the
            unique identifiers in keys.py
        
        :energy_channel: (dictionary, optional) defines which energy channel
            from which the value is desired. All observations and forecasts
            are organized by energy channel. e.g.
            {'min': 10, 'max': 10, 'units': 'MeV'}
            
        :index: (int, optional) - for json elements that may be arrays, index
            indicates which array element to choose; if not specified, will
            return the 0th entry in the array
          
          
        OUTPUTS:
        
        :sub: (varies) sub is the final unique value found be equating the
            unique identier "value" with a key_chain and then extracting
            that specific value from the json dictionary; returns None
            if the value cannot be found or the field is empty
        
        NOTE: values that are strings in zulu time will be converted to
            datetime prior to returning
            
    """
    
    #Check that the json file isn't empty
    if not injson:
        print("return_json_value: JSON file is empty.")
        return vars.errval
        
    #Get the sequence of keys for the desired value
    key_chain = keys.get_key_chain(value, index)
    if key_chain == vars.errval:
        return vars.errval
    
    #discover main key, e.g. 'sep_forecast_submission' or
    #'sep_observation_submission'
    main_key = return_main_key(injson)
    if main_key == vars.errval: return vars.errval
    
    #DEFAULT KEY CHAINS HAVE MODEL IDENTIFIERS IN THEM. IF LOOKING AT
    #OBSERVATION JSON, REPLACE APPROPRIATE FIELDS IN key_chain
    if main_key == keys.obs_main:
        key_chain = switch_model_to_obs_keys(key_chain)

    
    ##############TOP LEVEL INFO#################
    #Check if values saved at the top level, i.e. contacts,
    #model short name, etc.
    if key_chain[0] in injson[main_key].keys():
        sub = injson[main_key][key_chain[0]]
        for key in key_chain[1:]:
            if isinstance(key,int): #array index value
                if key >= len(sub): #check that index inside range
                    sub = vars.errval
                    break
                else:
                    sub=sub[key]
            else:
                if key not in sub.keys():
                    sub = vars.errval
                    break
                if key in sub.keys():
                    sub = sub[key]
        
        if sub == vars.errval:
            print("return_json_value: Keys for requested value " \
                + value + " not in json file: " + str(key_chain))
        else:
            #Check if value is a zulu time and convert to datetime
            for key in key_chain:
                if isinstance(key,int): continue
                if 'time' in key:
                    sub = zulu_to_time(sub)
        
        return sub #extracted down to a final value
    
    ##############UNDER FORECASTS OR OBSERVATIONS#################
    #Discover type key, i.e. 'forecasts' or 'observations'
    type_key = return_type_key(injson)
    if type_key == vars.errval:
        return vars.errval
    
    #Get the key for energy channel
    energy_chan_key = keys.get_key_chain(keys.id_energy_channel) #['energy_channel']
    
    #json[main_key][type_key] returns an array of forecasts or
    #observations, generally identified uniquely via energy channel.
    #Identify desired array by energy channel.
    Ndict = len(injson[main_key][type_key])
    sub = vars.errval #subdictionaries that will be narrowed down to single value
    for i in range(Ndict):
        if energy_chan_key[0] not in injson[main_key][type_key][i]:
            continue
        else:
            #select entry with desired energy channel
            if energy_channel == injson[main_key][type_key][i][energy_chan_key[0]]:
                #Pull out subdictionary for specific energy channel
                sub_dict = injson[main_key][type_key][i]
                 
                #check if library contatining desired value is in sub_dict
                if key_chain[0] not in sub_dict.keys():
                    continue
                #Extract the value specified by the key_chain
                else:
                    if len(key_chain) == 1:
                        return sub_dict[key_chain[0]] #We have our value!
                    else: #if dictionary
                        sub = sub_dict[key_chain[0]]
                        for key in key_chain[1:]:
                            if isinstance(key,int): #array index value
                                if key >= len(sub): #check that index inside range
                                    sub = vars.errval
                                    break
                                else:
                                    sub=sub[key]
                            else:
                                if key not in sub.keys():
                                    sub = vars.errval
                                    break
                                if key in sub.keys():
                                    sub = sub[key]
    
    if sub == vars.errval:
        print("return_json_value: Keys for requested value " \
            + value + " not in json file: " + str(key_chain))
    else:
        #Check if value is a zulu time and convert to datetime
        for key in key_chain:
            if isinstance(key,int): continue
            if 'time' in key:
                sub = zulu_to_time(sub)
    
    return sub #extracted down to a final value
                    
    


def return_json_value_by_index(injson, value, channel_index=0, index=0):
    """ Return the value of a specific entry in the json file listed in
        the fields under 'observations' or 'forecasts'. Select which block
        under under 'observations' or 'forecasts' to access by specifying
        channel_index.
        
        injson['forecasts'][channel_index]
        
        Then, if the desired value within that block is an array, e.g.
        fluences or event_lengths, specify which element in the array
        should be returned using index.
        
        injson['forecasts'][channel_index]['fluences'][index]
        
        INPUTS:
    
        :injson: (dictionary) is the dictionary created by reading in a CCMC
            SEP Scoreboard JSON file for a single SEP event or
            forecast period
            
        :value: (string) indicates value desired and may be any of the
            unique identifiers in keys.py
        
        :channel_index: (int) defines which entry in the forecast or
            observation array from which the value is desired. i.e.
            injson['sep_model_submission']['forecasts'] is an array and
            channel_index specifies which forecast to choose
            
        :index: (int, optional) - for json elements UNDER 'forecasts' or
            'observations' that may be arrays, index indicates which
            array element to choose, e.g. fluences is an array and
            index indicates which element of the array to choose;
            if not specified, will return the 0th entry in the array
            
        
        OUTPUTS:
        
        :sub: (varies) sub is the final unique value found be equating the
            unique identier "value" with a key_chain and then extracting
            that specific value from the json dictionary; returns None
            if the value cannot be found or the field is empty
        
        NOTE: values that are strings in zulu time will be converted to
            datetime prior to returning
        
    """
    
    #Check that the json file isn't empty
    if not injson:
        print("return_json_value: JSON file is empty.")
        return vars.errval
        
    #Get the sequence of keys for the desired value
    key_chain = keys.get_key_chain(value, index)
    if key_chain == vars.errval:
        return vars.errval
        
    #discover main key, e.g. 'sep_forecast_submission' or
    #'sep_observation_submission'
    main_key = return_main_key(injson)
    if main_key == vars.errval: return vars.errval
    
    #DEFAULT KEY CHAINS HAVE MODEL IDENTIFIERS IN THEM. IF LOOKING AT
    #OBSERVATION JSON, REPLACE APPROPRIATE FIELDS IN key_chain
    if main_key == keys.obs_main:
        key_chain = switch_model_to_obs_keys(key_chain)
    
    ##############TOP LEVEL INFO#################
    #Check if values saved at the top level, i.e.
    #model short name, etc.
    if key_chain[0] in injson[main_key].keys():
        sub = injson[main_key][key_chain[0]]
        for key in key_chain[1:]:
            if isinstance(key,int): #array index value
                if key >= len(sub): #check that index inside range
                    sub = vars.errval
                    break
                else:
                    sub=sub[key]
            else:
                if key not in sub.keys():
                    sub = vars.errval
                    break
                if key in sub.keys():
                    sub = sub[key]
        
        if sub == vars.errval:
            print("return_json_value: Keys for requested value " \
                + value + " not in json file: " + str(key_chain))
        else:
            #Check if value is a zulu time and convert to datetime
            for key in key_chain:
                if isinstance(key,int): continue
                if 'time' in key:
                    sub = zulu_to_time(sub)
        
        return sub #extracted down to a final value
    
    ##############UNDER FORECASTS OR OBSERVATIONS#################
    #Discover type key, i.e. 'forecasts' or 'observations'
    type_key = return_type_key(injson)
    if type_key == vars.errval:
        return vars.errval
    
    #json[main_key][type_key] returns an array of forecasts or
    #observations, generally identified uniquely via energy channel
    #and threshold. Identify desired array by specified channel_index.
    #Pull out subdictionary for specific energy channel
    sub_dict = injson[main_key][type_key][channel_index]
    #check if library contatining desired value is in sub_dict
    if key_chain[0] not in sub_dict.keys():
        return vars.errval
    #Extract the value specified by the key_chain
    else:
        if len(key_chain) == 1:
            return sub_dict[key_chain[0]] #We have our value!
        else: #if dictionary
            sub = sub_dict[key_chain[0]]
            for key in key_chain[1:]:
                if isinstance(key,int): #array index value
                    if key >= len(sub): #check that index inside range
                        sub = vars.errval
                        break
                    else:
                        sub=sub[key]
                else:
                    if key not in sub.keys():
                        sub = vars.errval
                        break
                    if key in sub.keys():
                        sub = sub[key]
    
    if sub == vars.errval:
        print("return_json_value_by_index: Keys for requested value " \
            + value + " not in json file: " + str(key_chain))
    else:
        if isinstance(sub,int):
            sub = float(sub)
        #Check if value is a zulu time and convert to datetime
        for key in key_chain:
            if isinstance(key,int): continue
            if 'time' in key:
                sub = zulu_to_time(sub)
    
    return sub #extracted down to a final value


def return_json_value_by_threshold(injson, value, energy_channel={}, threshold=0):
    """ Return the value of a specific entry in the json file listed in
        the fields under 'observations' or 'forecasts'.
        
        Select the desired energy block by matching to energy_channel.
        If the desired value is inside of an array (e.g. event_lengths,
        fluences, fluence_spectra, etc), find the correct element by
        matching the applied threshold value.
        
        Match energy_channel and threshold to pull out a unique value.
        For example, to get an event start time for a 10 pfu threshold
        applied to the >10 MeV channel:
        
        * energy_channel = {"min": 10.0, "max": -1, "units": "MeV"}
        * threshold = 10
        
        injson['sep_forecast_submission']['forecasts'][i]['event_lengths'][j]['start_time']
        
        where i was found by matching energy_channel and j was found by matching
        
        injson['sep_forecast_submission']['forecasts'][i]['event_lengths'][j]['threshold']
        
        INPUTS:
    
        :injson: (dictionary) is the dictionary created by reading in a CCMC
            SEP Scoreboard JSON file for a single SEP event or
            forecast period
            
        :value: (string) indicates value desired and may be any of the
            unique identifiers in keys.py
        
        :energy_channel: (dictionary, optional) defines which energy channel
            from which the value is desired. All observations and forecasts
            are organized by energy channel. e.g.
            {'min': 10, 'max': 10, 'units': 'MeV'}
            
        :threshold: (float) - for json elements that may be arrays, identify
            the desired entry to selecting the one associated with a
            threshold value
          
          
        OUTPUTS:
        
        :sub: (varies) sub is the final unique value found be equating the
            unique identier "value" with a key_chain and then extracting
            that specific value from the json dictionary; returns None
            if the value cannot be found or the field is empty
        
        NOTE: values that are strings in zulu time will be converted to
            datetime prior to returning
            
    """
    
    #Check that the json file isn't empty
    if not injson:
        print("return_json_value: JSON file is empty.")
        return vars.errval
        
    #Get the sequence of keys for the desired value
    key_chain = keys.get_key_chain(value)
    if key_chain == vars.errval:
        return vars.errval
        
    #If channel fluence is desired, need to take a different
    #approach and match to thresholds in the event_lengths
    #entry. CCMC requirement that those two arrays match up
    fluences_key_chain = []
    if 'fluences' in key_chain:
        fluences_key_chain = key_chain
        key_chain = keys.get_key_chain(keys.id_event_length_threshold)
    
    #discover main key, e.g. 'sep_forecast_submission' or
    #'sep_observation_submission'
    main_key = return_main_key(injson)
    if main_key == vars.errval: return vars.errval
    
    #DEFAULT KEY CHAINS HAVE MODEL IDENTIFIERS IN THEM. IF LOOKING AT
    #OBSERVATION JSON, REPLACE APPROPRIATE FIELDS IN key_chain
    if main_key == keys.obs_main:
        key_chain = switch_model_to_obs_keys(key_chain)

    
    ##############TOP LEVEL INFO#################
    #Check if values saved at the top level, i.e. contacts,
    #model short name, etc.
    if key_chain[0] in injson[main_key].keys():
        sub = injson[main_key][key_chain[0]]
        for key in key_chain[1:]:
            if isinstance(key,int): #array index value
                if key >= len(sub): #check that index inside range
                    sub = vars.errval
                    break
                else:
                    sub=sub[key]
            else:
                if key not in sub.keys():
                    sub = vars.errval
                    break
                if key in sub.keys():
                    sub = sub[key]
        
        if sub == vars.errval:
            print("return_json_value: Keys for requested value " \
                + value + " not in json file: " + str(key_chain))
        else:
            if isinstance(sub,int):
                sub = float(sub)
            #Check if value is a zulu time and convert to datetime
            for key in key_chain:
                if isinstance(key,int): continue
                if 'time' in key:
                    sub = zulu_to_time(sub)
        
        return sub #extracted down to a final value
    
    ##############UNDER FORECASTS OR OBSERVATIONS#################
    #Discover type key, i.e. 'forecasts' or 'observations'
    type_key = return_type_key(injson)
    if type_key == vars.errval:
        return vars.errval
    
    #Get the key for energy channel
    energy_chan_key = keys.get_key_chain(keys.id_energy_channel) #['energy_channel']
    
    #json[main_key][type_key] returns an array of forecasts or
    #observations, generally identified uniquely via energy channel
    #and threshold. Identify desired array by energy channel.
    Ndict = len(injson[main_key][type_key])
    sub = vars.errval #subdictionaries that will be narrowed down to single value
    thresh_index = -1
    channel_index = -1
    for i in range(Ndict):
        if energy_chan_key[0] not in injson[main_key][type_key][i]:
            continue
        else:
            #select entry with desired energy channel
            if energy_channel == injson[main_key][type_key][i][energy_chan_key[0]]:
                #Pull out subdictionary for specific energy channel
                sub_dict = injson[main_key][type_key][i]
                channel_index = i
                
                #check if library contatining desired value is in sub_dict
                if key_chain[0] not in sub_dict.keys():
                    continue
                #Extract the value specified by the key_chain
                else:
                    if len(key_chain) == 1:
                        return sub_dict[key_chain[0]] #We have our value!
                    else: #if dictionary
                        sub = sub_dict[key_chain[0]]
                        #Get correct threshold id
                        id_threshold = ""
                        if 'fluences' in key_chain \
                            or 'event_lengths' in key_chain:
                            id_threshold = 'threshold'
                        if 'fluence_spectra' in key_chain:
                            id_threshold = 'threshold_start'
                        if 'threshold_crossings' in key_chain:
                            id_threshold = 'threshold'
                        if 'probabilities' in key_chain:
                            id_threshold = 'threshold'
                            
                        for key in key_chain[1:]:
                            if isinstance(key,int): #array index value
                                #In this case, the index stored in the key
                                #is a dummy value. We want to identify
                                #the correct element of the array by matching
                                #the threshold
                                #sub is an array
                                for j in range(len(sub)):
                                    if id_threshold not in sub[j]:
                                        continue
                                    thresh = sub[j][id_threshold]
                                    if thresh == threshold:
                                        thresh_index = j
                                
                                if thresh_index == -1: #check that index inside range
                                    sub = vars.errval
                                    break
                                else:
                                    sub=sub[key]
                            else:
                                if key not in sub.keys():
                                    sub = vars.errval
                                    break
                                if key in sub.keys():
                                    sub = sub[key]
    
    #IF user requested channel fluence, then the sub currently
    #holds the value for event start time, since the event_lengths
    #and fluences arrays must match up and the threshold value is
    #only stored in the event_lengths array.
    #Get the value for the fluence
    if fluences_key_chain != [] and channel_index != -1 \
        and thresh_index != -1:
        sub = return_json_value_by_index(injson, value, channel_index,
                thresh_index)
    
    if sub == vars.errval:
        print("return_json_value: Keys for requested value " \
            + value + " not in json file: " + str(key_chain))
    else:
        if isinstance(sub,int):
            sub = float(sub)
        #Check if value is a zulu time and convert to datetime
        for key in key_chain:
            if isinstance(key,int): continue
            if 'time' in key:
                sub = zulu_to_time(sub)
    
    return sub #extracted down to a final value
