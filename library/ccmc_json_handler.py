import json
import calendar
import datetime
from datetime import timedelta
from datetime import datetime
import copy
import zulu
from library import global_vars as vars
import os

__version__ = "0.3"
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

email = vars.email
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
    """Make a datetime string in the format YYYY-MM-DDTHH:MMZ"""
    if dt == None:
        return None

    zdt = zulu.create(dt.year, dt.month, dt.day, dt.hour, dt.minute)
    stzdt = str(zdt)
    stzdt = stzdt.split(':00+00:00')
    zuludate = stzdt[0] + "Z"
    return zuludate


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


def fill_json(template, experiment, flux_type, energy_bins,
                model_name, startdate, enddate, options, energy_thresholds,
                flux_thresholds, crossing_time, onset_peak, onset_date,
                peak_flux, peak_time, rise_time, event_end_time, duration,
                all_integral_fluences, diff_thresh, all_energies, all_fluence,
                umasep, umasep_times, umasep_fluxes, profile_filename):
    """Add all the appropriate values to the json template for model or
        observations.
    """
    #For now, assume a user-input file is model output
    #This is not generic and should be modified in the future
    if experiment == "user":
        key = 'sep_forecast_submission'
        type_key = 'forecasts'
        win_key = 'prediction_window'
        template[key]['model']['short_name'] = model_name
        template[key]['model']['spase_id'] = 'spase://CCMC/SimulationModel/' \
                                            + model_name + "/" + version
    else:
        key = 'sep_observation_submission'
        type_key = 'observations'
        win_key = 'observation_window'
        template[key]['observatory']['short_name'] = experiment
        template[key]['observatory']['flux_type'] = flux_type


    template[key]['contacts'][0]['name'] = "operational_sep_quantities"
    template[key]['contacts'][0]['email'] = email
    template[key]['options'] = options
    now = datetime.now()
    znow = make_ccmc_zulu_time(now)
    template[key]['issue_time'] = znow

    #####FILL THRESHOLDS#############
    nthresh = len(energy_thresholds)
    #if template doesn't contain enough threshold entries, add more
    nent = len(template[key][type_key])
    if nent < nthresh:
        for j in range(nent,nthresh):
            add_entry = copy.deepcopy(template[key][type_key][nent-1])
            template[key][type_key].append(add_entry)

    #Fill in the info for each threshold
    for i in range(nthresh):
        #integral channel
        if not diff_thresh[i]:
            template[key][type_key][i]['energy_channel']['min'] \
                            = energy_thresholds[i]
            template[key][type_key][i]['energy_channel']['max'] = -1
            template[key][type_key][i]['peak_intensity']['units'] = "pfu"
            template[key][type_key][i]['peak_intensity_esp']['units'] \
                             = "pfu"
            template[key][type_key][i]['event_length']['units'] = "pfu"
            template[key][type_key][i]['threshold_crossings'][0]['threshold_units']\
                            = "pfu"
            template[key][type_key][i]['threshold_crossings'][1]['threshold_units']\
                            = "pfu"
            template[key][type_key][i]['all_clear']['threshold_units'] \
                            = "pfu"

        #differential channel
        if diff_thresh[i]:
            bin = find_energy_bin(energy_thresholds[i], energy_bins)
            template[key][type_key][i]['energy_channel']['min'] = bin[0]
            template[key][type_key][i]['energy_channel']['max'] = bin[1]
            template[key][type_key][i]['peak_intensity']['units'] \
                            = "[MeV/n cm^2 s sr]^(-1)"
            template[key][type_key][i]['peak_intensity_esp']['units'] \
                            = "[MeV/n cm^2 s sr]^(-1)"
            template[key][type_key][i]['event_length']['units'] \
                            = "[MeV/n cm^2 s sr]^(-1)"
            template[key][type_key][i]['threshold_crossings'][0]['threshold_units']\
                            = "[MeV/n cm^2 s sr]^(-1)"
            template[key][type_key][i]['threshold_crossings'][1]['threshold_units']\
                            = "[MeV/n cm^2 s sr]^(-1)"
            template[key][type_key][i]['all_clear']['threshold_units'] \
                            = "[MeV/n cm^2 s sr]^(-1)"


        zst = make_ccmc_zulu_time(startdate)
        zend = make_ccmc_zulu_time(enddate)
        template[key][type_key][i]['energy_channel']['units'] = "MeV"
        template[key][type_key][i]['species'] = "proton"
        template[key][type_key][i]['location'] = "earth"
        template[key][type_key][i][win_key]['start_time'] = zst
        template[key][type_key][i][win_key]['end_time'] = zend

        #Threshold WAS crossed
        if crossing_time[i] != 0:
            zodate = make_ccmc_zulu_time(onset_date[i])
            zpdate = make_ccmc_zulu_time(peak_time[i])
            zct = make_ccmc_zulu_time(crossing_time[i])
            zeet = make_ccmc_zulu_time(event_end_time[i])
            all_clear = False

            template[key][type_key][i]['peak_intensity']['intensity'] \
                                = onset_peak[i]
            template[key][type_key][i]['peak_intensity']['time'] \
                                = zodate
            template[key][type_key][i]['peak_intensity_esp']['intensity'] \
                                = peak_flux[i]
            template[key][type_key][i]['peak_intensity_esp']['time'] \
                                = zpdate
            template[key][type_key][i]['event_length']['start_time'] = zct
            template[key][type_key][i]['event_length']['end_time'] = zeet
            template[key][type_key][i]['event_length']['threshold'] \
                                = flux_thresholds[i]
            template[key][type_key][i]['fluence']['fluence_value'] \
                                = all_integral_fluences[i][i]
            template[key][type_key][i]['fluence']['units'] = 'cm^-2*sr^-1'
            template[key][type_key][i]['fluence_spectrum']['start_time'] = zct
            template[key][type_key][i]['fluence_spectrum']['end_time'] = zeet
            template[key][type_key][i]['fluence_spectrum']['energy_bins'] \
                                = all_energies[i].tolist()
            template[key][type_key][i]['fluence_spectrum']['fluences'] \
                                = all_fluence[i].tolist()
            flunits = 'cm^-2*sr-1'
            if diff_thresh[i]:
                flunits = 'MeV^-1*cm^-2*sr-1'
            template[key][type_key][i]['fluence_spectrum']['units'] = flunits

            template[key][type_key][i]['threshold_crossings'][0]['crossing_time']\
                                = zct
            template[key][type_key][i]['threshold_crossings'][0]['threshold']\
                                = flux_thresholds[i]
            template[key][type_key][i]['threshold_crossings'][1]['crossing_time']\
                                = zeet
            template[key][type_key][i]['threshold_crossings'][1]['threshold']\
                                = flux_thresholds[i]*0.85
            template[key][type_key][i]['all_clear']['all_clear_boolean'] \
                                = all_clear
            template[key][type_key][i]['all_clear']['threshold'] \
                                = flux_thresholds[i]
            template[key][type_key][i]['sep_profile'] = profile_filename


        #Threshold was NOT crossed OR doesn't exist in data set
        ######RETHINK LOGIC HERE FOR DIFFERENTIAL FLUXES THAT WERE USED
        #####TO ESTIMATE INTEGRAL FLUXES###############################
        if crossing_time[i] == 0:
            #Check if channel exists in data; if not, output None values
            bin = find_energy_bin(energy_thresholds[i], energy_bins)
            if not bin: #bin doesn't exist and no prediction can be made
                zodate = None
                zpdate = None
                zct = None
                zeet = None
                all_clear = None
                template[key][type_key][i]['sep_profile'] = None

                if flux_type == "differential" and not diff_thresh[i]:
                    zodate = ""
                    zpdate = ""
                    zct = ""
                    zeet = ""
                    all_clear = True
                    template[key][type_key][i]['sep_profile'] = profile_filename

            else:  #bin does exist and equivalent to all clear for channel
                zodate = make_ccmc_zulu_time(onset_date[i])
                zpdate = ""
                zct = ""
                zeet = ""
                all_clear = True
                template[key][type_key][i]['sep_profile'] = profile_filename


            template[key][type_key][i]['peak_intensity']['intensity'] \
                                = onset_peak[i]
            template[key][type_key][i]['peak_intensity']['time'] \
                                = zodate
            template[key][type_key][i]['peak_intensity_esp']['intensity'] \
                                = -999
            template[key][type_key][i]['peak_intensity_esp']['time'] \
                                = zpdate
            template[key][type_key][i]['event_length']['start_time'] = zct
            template[key][type_key][i]['event_length']['end_time'] = zeet
            template[key][type_key][i]['event_length']['threshold'] \
                                = flux_thresholds[i]
            template[key][type_key][i]['threshold_crossings'][0]['crossing_time']\
                                = zct
            template[key][type_key][i]['threshold_crossings'][0]['threshold']\
                                = flux_thresholds[i]
            template[key][type_key][i]['threshold_crossings'][1]['crossing_time']\
                                = zeet
            template[key][type_key][i]['threshold_crossings'][1]['threshold']\
                                = flux_thresholds[i]*0.85
            template[key][type_key][i]['all_clear']['all_clear_boolean'] \
                                = all_clear
            template[key][type_key][i]['all_clear']['threshold'] \
                                = flux_thresholds[i]

    return template


def write_json(template, filename):
    """Write json template to json file. """
    with open(filename, "w") as outfile:
        json.dump(template, outfile)

    if not os.path.isfile(filename):
        return False

    print("Wrote SEP values to json file --> " + filename)
    return True
