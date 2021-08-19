from library import global_vars as vars

__version__ = "0.3"
__author__ = "Katie Whitman"
__maintainer__ = "Katie Whitman"
__email__ = "kathryn.whitman@nasa.gov"

#2021-07-23, Changes in 0.2: Remove the contacts fields due to
#   NASA privacy rules.
#   Updating to follow updated CCMC JSON format. This version
#   is consistent with version 3.0+ of operational_sep_quantities.py
#2021-08-18, Changes in 0.3: added identifiers for all of the possible
#   fields in the json files produced by operational_sep_quantities.py
#   or that may be produced by models for the CCMC SEP Scoreboard.
#   This code matches the same version of the code in the validation
#   project library. Must keep both of these version matched and
#   up-to-date with each other.

def about_keys():
    """ About keys.py
    
        This version of keys.py is maintained in the operational-sep
        project.
    
        Define the keys used to call different variables in the CCMC
        model and observation JSON files.
        
        This library defines the keys that are different
        between the observations json and the model output json.
        
        It also defines all the keys that are shared between the different
        type of json files.
        
        This program coupled with ccmc_json_handler.py will access
        all of the possible model outputs in the CCMC SEP Scoreboard json files
        and the files produced by operational_sep_quantities.py.
        https://github.com/ktindiana/operational-sep
    """


#KEYS FOR OBSERVATIONS
obs_main = 'sep_observation_submission'
obs_exp = 'observatory'
obs_type = 'observations'
obs_win = 'observation_window'

#KEYS FOR MODELS
model_main = 'sep_forecast_submission'
model_exp = 'model'
model_type = 'forecasts'
model_win = 'prediction_window'

#The remaining keys in the JSON files are shared between model and observation

#IDS linked to KEYS for accessing
#info in CCMC SEP Scoreboard JSON files
#In main header info of json file
#id_contacts = "contacts"
#id_contact_name = "contact-name" #removed because PII
#id_contact_email = "contact-email"
id_short_name = "short-name"
id_spase_id = "spase-id"
id_issue_time = "issue-time"
id_options = "options"

###UNDER 'forecasts' or 'observations'
id_energy_channel = "energy-channel"
id_energy_min = "energy-min"
id_energy_max = "energy-max"
id_energy_units = "energy-units"
id_species = "species"
id_location = "location"
id_prediction_window = "prediction-window"
id_prediction_window_start = "prediction-window-start"
id_prediction_window_end = "prediction-window-end"
id_peak_intensity = "peak-intensity"
id_peak_intensity_units = "peak-intensity-units"
id_peak_intensity_uncertainty = "peak-intensity-uncertainty"
id_peak_intensity_uncertainty_low = "peak-intensity-uncertainty-low"
id_peak_intensity_uncertainty_high = "peak-intensity-uncertainty-high"
id_peak_intensity_time = "peak-intensity-time"
id_peak_intensity_esp = "peak-intensity-esp"
id_peak_intensity_esp_units = "peak-intensity-esp-units"
id_peak_intensity_esp_time = "peak-intensity-esp-time"
id_peak_intensity_max = "peak-intensity-max"
id_peak_intensity_max_units = "peak-intensity-max-units"
id_peak_intensity_max_uncertainty = "peak-intensity-max-uncertainty"
id_peak_intensity_max_uncertainty_low = "peak-intensity-max-uncertainty-low"
id_peak_intensity_max_uncertainty_high = "peak-intensity-max-uncertainty-high"
id_peak_intensity_max_time = "peak-intensity-max-time"
id_fluences = "fluences"
id_fluence = "fluence"
id_fluence_units = "fluence-units"
id_fluence_uncertainty_low = "fluence-uncertainty-low"
id_fluence_uncertainty_high = "fluence-uncertainty-high"
id_fluence_spectra = "fluence-spectra"
id_fluence_spectrum_start_time = "fluence-spectrum-start-time"
id_fluence_spectrum_end_time = "fluence-spectrum-end-time"
id_fluence_spectrum_threshold_start = "fluence-spectrum-threshold-start"
id_fluence_spectrum_threshold_end = "fluence-spectrum-threshold-end"
id_fluence_spectrum_threshold_units = "fluence-spectrum-threshold-units"
id_fluence_spectrum_fluence_units = "fluence-spectrum-fluence-units"
id_fluence_spectrum = "fluence-spectrum"
#id_fluence_spectrum_energy_min = "fluence-spectrum-energy-min"
#id_fluence_spectrum_energy_max = "fluence-spectrum-energy-max"
#id_fluence_spectrum_fluence = "fluence-spectrum-fluence"
id_event_lengths = "event-lengths"
id_event_length_start_time = "event-length-start-time"
id_event_length_end_time = "event-length-end-time"
id_event_length_threshold = "event-length-threshold"
id_event_length_threshold_units = "event-length-threshold-units"
id_threshold_crossings = "threshold-crossings"
id_thresh_crossing_time = "thresh-crossing-time"
id_thresh_uncertainty = "thresh-uncertainty"
id_crossing_threshold = "crossing-threshold"
id_crossing_threshold_units = "crossing-threshold-units"
id_probabilities = "probabilities"
id_probability = "probability"
id_prob_uncertainty = "prob-uncertainty"
id_prob_threshold = "prob-threshold"
id_prob_threshold_units = "prob-threshold-units"
id_all_clear = "all-clear"
id_all_clear_threshold = "all-clear-threshold"
id_all_clear_threshold_units = "all-clear-threshold-units"
id_all_clear_probability_threshold = "all-clear-probability-threshold"
id_sep_profile = "sep-profile"

id_all = [id_short_name,id_spase_id,id_issue_time,id_options,
            id_energy_channel,id_energy_min,id_energy_max,id_energy_units,
            id_species,id_location,id_prediction_window,id_prediction_window_start,
            id_prediction_window_end,id_peak_intensity,id_peak_intensity_units,
            id_peak_intensity_uncertainty,id_peak_intensity_uncertainty_low,
            id_peak_intensity_uncertainty_high,id_peak_intensity_time,
            id_peak_intensity_esp,id_peak_intensity_esp_units,
            id_peak_intensity_esp_time,id_fluences,id_fluence,id_fluence_units,
            id_fluence_uncertainty_low,id_fluence_uncertainty_high,
            id_fluence_spectra,id_fluence_spectrum_start_time,id_fluence_spectrum_end_time,
            id_fluence_spectrum_threshold_start,id_fluence_spectrum_threshold_end,
            id_fluence_spectrum_threshold_units,id_fluence_spectrum_fluence_units,
            id_fluence_spectrum,
            id_event_lengths,id_event_length_start_time,id_event_length_end_time,
            id_event_length_threshold,id_event_length_threshold_units,
            id_threshold_crossings,id_thresh_crossing_time,id_thresh_uncertainty,
            id_crossing_threshold,id_crossing_threshold_units,id_probabilities,
            id_probability,id_prob_uncertainty,id_prob_threshold,
            id_prob_threshold_units,id_all_clear,id_all_clear_threshold,
            id_all_clear_threshold_units,id_all_clear_probability_threshold,
            id_sep_profile]


def get_key_chain(value, index=0):
    """ Return the series of arrays that points to specific values
        in a CCMC JSON file compatible with the SEP Scoreboard.
        Combine this module, keys.py, with validation_json_handler.py
        to read CCMC SEP Scoreboard JSONs or JSONs produced by
        operational_sep_quantities.py.
        
        INPUTS:
        
        :value: (string) is a specifier that indicates the value desired.
        :index: (integer) is an optional input that indicates which index
            is desired if value is stored in an array.
            
        OUTPUTS:
        
        :key_chain: (list, array) contains the list of keys needed to access
            the desired value in the json file
            
        ADDITIONAL INFO:
        
        The string id to specify for value is the same as the ones listed at
        the top of this code and also specified in the json visual scheme document
        on the CCMC webpage: https://ccmc.gsfc.nasa.gov/challenges/sep.php#format
        These identifiers are followed as closely as possible, with some
        additions specific to the operational_sep_quantities.py code and
        this SEP model validation effort.
        
        e.g. "peak-intensity" requests the chain of entries that will
        give the location of the onset peak intensity in the json file.
        key_chain = ['peak_intensity']['intensity']
        
        e.g. peak-intensity-units
        key_chain = ['peak_intensity']['units']
            
            
        Some of the json entries are contained within arrays. index
        specifies which element of the array to select.
        e.g. ['fluences'] = [{'fluence': , 'units': }]
            
        index applies to 'fluences', 'fluence_spectra', 'events_lengths',
            'threshold_crossings','probabilities', 'contacts'
        
        The series of string arrays output as the key_chain are with respect to the
        ['forecasts'] or ['observations'] field. Model output will have a
        'forecasts' field while observational output will have an 'observations'
        field. An example of a model JSON for a time profile model is shown below.
        The key_chain returned will be the key_chain following the 'forecasts'
        or 'observations' field.
        
        .. code-block::
        
            {"sep_forecast_submission": {"notes": [{"note": "produced by operational_sep_quantities.py"}], "model": {"short_name": "SEPMOD_RT_60min", "spase_id": "spase://CCMC/SimulationModel/SEPMOD", "flux_type": "integral"}, "options": [""], "issue_time": "2021-08-06T17:21:40Z", "mode": "historical",
            "forecasts": [
                {"energy_channel": {"min": 10, "max": -1, "units": "MeV"},
                    "species": "proton",
                    "location": "earth",
                    "prediction_window": {"start_time": "2021-05-29T00:00:00Z", "end_time": "2021-06-05T00:00:00Z"},
                    "peak_intensity": {"intensity": 35.53, "units": "pfu", "time": "2021-05-29T09:30:00Z"},
                    "peak_intensity_max": {"intensity": 35.53, "units": "pfu", "time": "2021-05-29T09:30:00Z"},
                    "event_lengths": [{"start_time": "2021-05-29T08:30:00Z", "end_time": "2021-05-29T22:30:00Z", "threshold": 10, "threshold_units": "pfu"}, {"start_time": "2021-05-29T05:30:00Z", "end_time": "2021-05-31T12:30:00Z", "threshold": 0.001, "threshold_units": "pfu"}],
                    "fluences": [{"fluence": 12647177.403957747, "units": "cm^-2"}, {"fluence": 15077470.358017407, "units": "cm^-2"}],
                    "fluence_spectra": [
                        {"start_time": "2021-05-29T08:30:00Z", "end_time": "2021-05-29T22:30:00Z", "threshold_start": 10, "threshold_end": 8.5, "threshold_units": "pfu", "fluence_units": "cm^-2",
                        "fluence_spectrum":[
                            {"energy_min": 10, "energy_max": -1, "fluence": 12647177.403957747},
                            {"energy_min": 30, "energy_max": -1, "fluence": 1065408.567939319},
                            {"energy_min": 50, "energy_max": -1, "fluence": 458921.84421709884},
                            {"energy_min": 60, "energy_max": -1, "fluence": 365385.8238410022},
                            {"energy_min": 100, "energy_max": -1, "fluence": 191179.73597861468},
                            {"energy_min": 300, "energy_max": -1, "fluence": 14478.720894452352},
                            {"energy_min": 500, "energy_max": -1, "fluence": 4452.008754706922},
                            {"energy_min": 750, "energy_max": -1, "fluence": 1175.6965656540053}
                            ]},
                        {"start_time": "2021-05-29T05:30:00Z", "end_time": "2021-05-31T12:30:00Z", "threshold_start": 0.001, "threshold_end": 0.00085, "threshold_units": "pfu", "fluence_units": "cm^-2",
                        "fluence_spectrum": [
                            {"energy_min": 10, "energy_max": -1, "fluence": 15077470.358017407},
                            {"energy_min": 30, "energy_max": -1, "fluence": 1248155.4091932934},
                            {"energy_min": 50, "energy_max": -1, "fluence": 534833.8687969703},
                            {"energy_min": 60, "energy_max": -1, "fluence": 425533.462349928},
                            {"energy_min": 100, "energy_max": -1, "fluence": 222352.03527148295},
                            {"energy_min": 300, "energy_max": -1, "fluence": 16663.572499087855},
                            {"energy_min": 500, "energy_max": -1, "fluence": 5105.88749185624},
                            {"energy_min": 750, "energy_max": -1, "fluence": 1346.7982515392907}
                            ]}
                        ],
                    "threshold_crossings": [{"crossing_time": "2021-05-29T08:30:00Z", "threshold": 10, "threshold_units": "pfu"}, {"crossing_time": "2021-05-29T05:30:00Z", "threshold": 0.001, "threshold_units": "pfu"}],
                    "all_clear": {"all_clear_boolean": false, "threshold": 10, "threshold_units": ""},
                    "sep_profile": "SEPMOD_RT_60min_integral.2021-05-29T000000Z.2021-08-06T172140Z.10.0MeV.txt"}
                    ]
                }
            }
        
    """
    key_chain = []
    ##########################################################################
    ###ALL VALUES FALL UNDER 'sep_forecast_submission' OR 'sep_observations_submission'
    ##########################################################################
    
    #if value == id_contacts:
    #    key_chain = ['contacts']                       #points to (-->) an array of dictionaries
    #if value == id_contact_name:
    #    key_chain = ['contacts',index,'name']          #--> string
    #if value == id_contact_email:
    #    key_chain = ['contacts',index,'email']         #--> string
        
    if value == id_short_name:
        key_chain = ['model','short_name']
    if value == id_spase_id:
        key_chain = ['model','spase_id']
    if value == id_issue_time:
        key_chain = ['issue_time']                      #--> string zulu time
    if value == id_options:                             #--> array of strings
        key_chain = ['options']
    ###########################################################################
    ####ALL BELOW FALL WITHIN BOTH 'observations' OR 'forecasts' DICTIONARIES
    if value == id_energy_channel:
        key_chain = ['energy_channel']                  #--> library
    if value == id_energy_min:
        key_chain = ['energy_channel','min']            #--> float
    if value == id_energy_max:
        key_chain = ['energy_channel','max']            #--> float
    if value == id_energy_units:
        key_chain = ['energy_channel','units']          #--> string
     
    if value == id_species:
        key_chain = ['species']                         #--> string
    if value == id_location:
        key_chain =['location']                         #--> string
    
    if value == id_prediction_window:
        key_chain = ['prediction_window']               #--> library
    if value == id_prediction_window_start:
        key_chain = ['prediction_window','start_time']  #--> string zulu time
    if value == id_prediction_window_end:
        key_chain = ['prediction_window','end_time']    #--> string zulu time
        
    if value == id_peak_intensity:
        key_chain = ['peak_intensity','intensity']      #--> float
    if value == id_peak_intensity_units:
        key_chain = ['peak_intensity','units']          #--> string
    if value == id_peak_intensity_uncertainty:
        key_chain = ['peak_intensity','uncertainty']    #--> float
    if value == id_peak_intensity_uncertainty_low:
        key_chain = ['peak_intensity','uncertainty_low'] #--> float
    if value == id_peak_intensity_uncertainty_high:
        key_chain = ['peak_intensity','uncertainty_high'] #--> float
    if value == id_peak_intensity_time:
        key_chain = ['peak_intensity','time']           #--> string zulu time
        
    if value == id_peak_intensity_esp:
        key_chain = ['peak_intensity_esp','intensity']  #--> float
    if value == id_peak_intensity_esp_units:
        key_chain = ['peak_intensity_esp','units']      #--> string
    if value == id_peak_intensity_esp_time:
        key_chain = ['peak_intensity_esp','time']       #--> string zulu time
        
    if value == id_peak_intensity_max:
        key_chain = ['peak_intensity_max','intensity'] #--> float
    if value == id_peak_intensity_max_units:
        key_chain = ['peak_intensity_max','units']     #--> string
    if value == id_peak_intensity_max_uncertainty:
        key_chain = ['peak_intensity_max','uncertainty'] #--> float
    if value == id_peak_intensity_max_uncertainty_low:
        key_chain = ['peak_intensity_max','uncertainty_low'] #--> float
    if value == id_peak_intensity_max_uncertainty_high:
        key_chain = ['peak_intensity_max','uncertainty_high'] #--> float
    if value == id_peak_intensity_max_time:
        key_chain = ['peak_intensity_max','time']      #--> string zulu time

    if value == id_fluences:
        key_chain = ['fluences']                         #--> array of dictionaries
    if value == id_fluence:
        key_chain = ['fluences',index,'fluence']        #--> float
    if value == id_fluence_units:
        key_chain = ['fluences',index,'units']          #--> string
    if value == id_fluence_uncertainty_low:
        key_chain = ['fluences',index,'uncertainty_low'] #--> float
    if value == id_fluence_uncertainty_high:
        key_chain = ['fluences',index,'uncertainty_high'] #--> float
        
    if value == id_fluence_spectra:                     #--> array of dictionaries
        key_chain = ['fluence_spectra']
    if value == id_fluence_spectrum_start_time:         #--> string zulu time
        key_chain = ['fluence_spectra',index,'start_time']
    if value == id_fluence_spectrum_end_time:           #--> string zulu time
        key_chain = ['fluence_spectra',index,'end_time']
    if value == id_fluence_spectrum_threshold_start:    #--> float
        key_chain = ['fluence_spectra',index,'threshold_start']
    if value == id_fluence_spectrum_threshold_end:      #--> float
        key_chain = ['fluence_spectra',index,'threshold_end']
    if value == id_fluence_spectrum_threshold_units:    #--> string
        key_chain = ['fluence_spectra',index,'threshold_units']
    if value == id_fluence_spectrum_fluence_units:      #--> string
        key_chain = ['fluence_spectra',index,'fluence_units']
    if value == id_fluence_spectrum:                    #--> array of dicts
        key_chain = ['fluence_spectra',index,'fluence_spectrum']
    #Would need a second index value to uniquely specify values below
    #Leave it so that return the complete spectrum dictionary
    #if value == id_fluence_spectrum_energy_min = "fluence-spectrum-energy-min"
    #if value == id_fluence_spectrum_energy_max = "fluence-spectrum-energy-max"
    #if value == id_fluence_spectrum_fluence = "fluence-spectrum-fluence"
    
    if value == id_event_lengths:
        key_chain = ['event_lengths']                    #--> array of dictionaries
    if value == id_event_length_start_time:
        key_chain = ['event_lengths',index,'start_time']#--> string zulu time
    if value == id_event_length_end_time:
        key_chain = ['event_lengths',index,'end_time']  #--> string zulu time
    if value == id_event_length_threshold:
        key_chain = ['event_lengths',index,'threshold'] #--> float
    if value == id_event_length_threshold_units:
        key_chain = ['event_lengths',index,'threshold_units'] #--> string
        
    if value == id_threshold_crossings:
        key_chain = ['threshold_crossings']              #--> array of dictionaries
    if value == id_thresh_crossing_time:
        key_chain = ['threshold_crossings',index,'crossing_time'] #--> string zulu time
    if value == id_thresh_uncertainty:
        key_chain = ['threshold_crossings',index,'uncertainty'] #--> hours
    if value == id_crossing_threshold:
        key_chain = ['threshold_crossings',index,'threshold']  #--> float
    if value == id_crossing_threshold_units:
        key_chain = ['threshold_crossings',index,'threshold_units']  #--> string
        
    if value == id_probabilities:
        key_chain = ['probabilities']                    #--> array of dictionaries
    if value == id_probability:
        key_chain = ['probabilities',index,'probability_value'] #--> float 0 to 1
    if value == id_prob_uncertainty:
        key_chain = ['probabilities',index,'uncertainty']#--> float
    if value == id_prob_threshold:
        key_chain = ['probabilities',index,'threshold']  #--> float
    if value == id_prob_threshold_units:
        key_chain = ['probabilities',index,'threshold_units']  #--> string
        
    if value == id_all_clear:
        key_chain = ['all_clear','all_clear_boolean']   #--> bool
    if value == id_all_clear_threshold:
        key_chain = ['all_clear','threshold']           #--> float
    if value == id_all_clear_threshold_units:
        key_chain = ['all_clear','threshold_units']     #--> string
    if value == id_all_clear_probability_threshold:
        key_chain = ['all_clear','probability_threshold']   #--> float
    
    
    if value == id_sep_profile:
        key_chain = ['sep_profile']                      #--> string of filename
        
    if key_chain == []:
        print("get_key_chain: Could not find requested value " + value)
        key_chain = vars.errval
        
    return key_chain
        

    
    
    

