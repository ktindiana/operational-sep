import json
import calendar
import datetime
from datetime import timedelta
import copy
import zulu
from library import global_vars as vars
from library import keys
import os

__version__ = "1.0"
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

version = vars.version

def read_in_json_template(type):
    """Read in appropriate json file template for model or observations.
        type = model or observations
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
    """Make a datetime string in the format YYYY-MM-DDTHH:MM:SSZ"""
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
    """Convert Zulu time to datetime"""
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
    """Identify the energy bin and return low and high edges."""
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
        channels.
        e.g. >10 MeV, 10 pfu (operations) and
             >10 MeV, 0.001 pfu (SEPMOD testing)
        
        input:
        energy_thresholds (float array)- array containing
            all the energy channels for which a threshold
            has been applied
        
        output:
        (integer) - number of unique energy channels
    """
    unique = []
    for energy in energy_thresholds:
        if energy not in unique:
            unique.append(energy)
    
    return len(unique), unique


def fill_json(template, issue_time, experiment, flux_type, energy_bins,
                model_name, spase_id, startdate, enddate, options,
                energy_thresholds, flux_thresholds, crossing_time,
                onset_peak, onset_date, peak_flux, peak_time, rise_time,
                event_end_time, duration, all_threshold_fluences,
                diff_thresh, all_fluence,
                umasep, umasep_times, umasep_fluxes, profile_filenames,
                energy_units, flux_units_integral, fluence_units_integral,
                flux_units_differential, fluence_units_differential):
    """Add all the appropriate values to the json template for model or
        observations.
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


def read_in_json(filename):
    """Read in json file """
    if not os.path.isfile(filename):
        sys.exit("validation_json_handler: could not read in file " \
                + filename + ". Exiting.")

    with open(filename) as f:
        info=json.load(f)

    return info
