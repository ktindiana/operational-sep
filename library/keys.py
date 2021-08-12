from library import global_vars as vars

__version__ = "0.2"
__author__ = "Katie Whitman"
__maintainer__ = "Katie Whitman"
__email__ = "kathryn.whitman@nasa.gov"

#Define the keys used to call different variables in the CCMC and observation
#JSON files

#2021-07-23, Changes in 0.2: Remove the contacts fields due to
#   NASA privacy rules.
#   Updating to follow updated CCMC JSON format. This version
#   is consistent with version 3.0+ of operational_sep_quantities.py

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

#IDS linked to KEYS for accessing
#info in CCMC SEP Scoreboard JSON files
#In main header info of json file
id_short_name = "short-name"
id_spase_id = "spase-id"
id_issue_time = "issue-time"

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
#####ADD FLUENCE SPECTRA
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

id_all = [  id_short_name,
            id_spase_id,id_issue_time,
            id_energy_channel,id_energy_min,id_energy_max,id_energy_units,
            id_species,id_location,id_prediction_window,
            id_prediction_window_start,id_prediction_window_end,
            id_peak_intensity,id_peak_intensity_units,
            id_peak_intensity_uncertainty,id_peak_intensity_uncertainty_low,
            id_peak_intensity_uncertainty_high,id_peak_intensity_time,
            id_peak_intensity_esp,id_peak_intensity_esp_units,
            id_peak_intensity_esp_time,
            id_peak_intensity_max,id_peak_intensity_max_units,
            id_peak_intensity_max_uncertainty,
            id_peak_intensity_max_uncertainty_low,
            id_peak_intensity_max_uncertainty_high,id_peak_intensity_max_time,
            id_fluences,id_fluence,id_fluence_units,
            id_fluence_uncertainty_low,id_fluence_uncertainty_high,
            id_event_lengths,
            id_event_length_start_time,id_event_length_end_time,
            id_event_length_threshold,id_event_length_threshold_units,
            id_threshold_crossings,id_thresh_crossing_time,id_thresh_uncertainty,
            id_crossing_threshold,id_crossing_threshold_units,id_probabilities,
            id_probability,id_prob_uncertainty,id_prob_threshold,
            id_prob_threshold_units,id_all_clear,id_all_clear_threshold,
            id_all_clear_threshold_units,id_all_clear_probability_threshold,
            id_sep_profile]


#The remaining keys in the JSON files are shared between model and observation

def get_key_chain(value, index=0):
    """ Return the series of arrays that points to a specific values
        in a CCMC JSON file.
        
        value (string) is a specifier that indicates the value desired.
            e.g. peak-intensity requests the chain of entries that will
            results in the location of the peak intensity in the json file.
            key_chain = ['peak_intensity']['intensity']
            
            e.g. peak-intensity-units
            key_chain = ['peak_intensity']['units']
            
            value is any of the string specifiers used in the CCMC
            JSON helper script and listed in the "Corresponding
            command line argument for helper script" column in the table:
            https://ccmc.gsfc.nasa.gov/challenges/sepinfo/sepscoreboard_visual_schema_v20201112.pdf
            (check for most up-to-date version of the Visual Schema at https://ccmc.gsfc.nasa.gov/challenges/sep.php#format)
            These identifiers are followed as closely as possible, with some
            additions.
            
            
        index (integer) is an optional input; Some of the json entries
            are contained within arrays. index specifies which element
            of the array to select.
            e.g. ['fluences'] = [{'fluence': , 'units': }]
            index applies to 'fluences', 'events_lengths', 'threshold_crossings',
                'probabilities', 'contacts'
        
        The series of string arrays output as the key_chain are with respect to the
        ['forecasts'] or ['observations'] field. Model output will have a
        'forecasts' field while observational output will have an 'observations'
        field. An example of a model JSON for a time profile model is shown below.
        The key_chain returned will be the key_chain following the 'forecasts'
        or 'observations' field.
        
        {"sep_observation_submission": {
            "notes": [ { "note": "produced by operational_sep_quantities.py"} ],
            "observatory": {"short_name": "GOES-13", "spase_id": "", "flux_type": "differential"},
            "options": [""],
            "issue_time": "2021-07-23T16:25Z",
            "mode": "measurement",
            "observations": [
                {
                "energy_channel": {"min": 10, "max": -1, "units": "meV"},
                "species":"proton",
                "location": "earth",
                "observation_window": {"start_time": "2011-06-07T00:00Z", "end_time": "2011-06-12T00:00Z"},
                "peak_intensity": {"intensity": 55.11014913571528, "units": "pfu", "time": "2011-06-07T10:25Z"},
                "peak_intensity_max": {"intensity": 75.67780719646345, "units": "pfu", "time": "2011-06-07T18:20Z"},
                "event_lengths": [{"start_time": "2011-06-07T07:50Z", "end_time": "2011-06-08T20:25Z", "threshold": 10, "threshold_units": "pfu"}],
                "fluences": [{"fluence": 50523588.099411, "units": "cm^-2"}],
                "fluence_spectra": [{"start_time": "2011-06-07T07:50Z", "end_time": "2011-06-08T20:25Z", "threshold_start": 10, "threshold_end": 8.5, "threshold_units": "pfu", "fluence_units": "cm^-2",
                    "fluence_spectrum": [
                        {"energy_min": 4.2, "energy_max": 8.7, "fluence": 9471121.886108486},
                        {"energy_min": 8.7, "energy_max": 14.5, "fluence": 2390875.507650991},
                        {"energy_min": 15.0, "energy_max": 40.0, "fluence": 1309160.620161789},
                        {"energy_min": 38.0, "energy_max": 82.0, "fluence": 221348.9742486922},
                        {"energy_min": 84.0, "energy_max": 200.0, "fluence": 21820.507114220098},
                        {"energy_min": 110.0, "energy_max": 900.0, "fluence": 4491.899843708999},
                        {"energy_min": 330.0, "energy_max": 420.0, "fluence": 3346.346055318471},
                        {"energy_min": 420.0, "energy_max": 510.0, "fluence": 1892.0212727928365},
                        {"energy_min": 510.0, "energy_max": 700.0, "fluence": 887.989987960335},
                        {"energy_min": 700.0, "energy_max": -1, "fluence": 209.68393033655906}]}],
                "threshold_crossings": [{"crossing_time": "2011-06-07T07:50Z", "threshold": 10, "threshold_units": "pfu"}],
                "all_clear": {"all_clear_boolean": false, "threshold": 10, "threshold_units": "pfu"},
                "sep_profile": "GOES-13_differential.2011-06-07T0000Z.10meV.txt"}
            ]
            }
        }
        
    """
    key_chain = []
    ##################################################################################
    ###ALL VALUES FALL UNDER 'sep_forecast_submission' OR 'sep_observations_submission'
    ##################################################################################
        
    if value == id_short_name:
        key_chain = ['model','short_name']
    if value == id_spase_id:
        key_chain = ['model','spase_id']
    if value == id_issue_time:
        key_chain = ['issue_time']                      #--> string zulu time
    
    ##################################################################################
    ####ALL BELOW ARE UNDER 'observations' OR 'forecasts'
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
        key_chain = ['peak_intensity_max'],['intensity'] #--> float
    if value == "peak-intensity-max-units":
        key_chain = ['peak_intensity_max','units']     #--> string
    if value == "peak-intensity-max-uncertainty":
        key_chain = ['peak_intensity_max','uncertainty'] #--> float
    if value == "peak-intensity-max-uncertainty-low":
        key_chain = ['peak_intensity_max','uncertainty_low'] #--> float
    if value == "peak-intensity-max-uncertainty-high":
        key_chain = ['peak_intensity_max','uncertainty_high'] #--> float
    if value == "peak-intensity-max-time":
        key_chain = ['peak_intensity_max','time']      #--> string zulu time

    if value == id_fluences:
        key_chain = ['fluences']                         #--> array of libraries
    if value == id_fluence:
        key_chain = ['fluences',index,'fluence']        #--> float
    if value == id_fluence_units:
        key_chain = ['fluences',index,'units']          #--> string
    if value == id_fluence_uncertainty_low:
        key_chain = ['fluences',index,'uncertainty_low'] #--> float
    if value == id_fluence_uncertainty_high:
        key_chain = ['fluences',index,'uncertainty_high'] #--> float
    
    if value == id_event_lengths:
        key_chain = ['event_lengths']                    #--> array of libraries
    if value == id_event_length_start_time:
        key_chain = ['event_lengths',index,'start_time']#--> string zulu time
    if value == id_event_length_end_time:
        key_chain = ['event_lengths',index,'end_time']  #--> string zulu time
    if value == id_event_length_threshold:
        key_chain = ['event_lengths',index,'threshold'] #--> float
    if value == id_event_length_threshold_units:
        key_chain = ['event_lengths',index,'threshold_units'] #--> string
        
    if value == id_threshold_crossings:
        key_chain = ['threshold_crossings']              #--> array of libraries
    if value == id_thresh_crossing_time:
        key_chain = ['threshold_crossings',index,'crossing_time'] #--> string zulu time
    if value == id_thresh_uncertainty:
        key_chain = ['threshold_crossings',index,'uncertainty'] #--> hours
    if value == id_crossing_threshold:
        key_chain = ['threshold_crossings',index,'threshold']  #--> float
    if value == id_crossing_threshold_units:
        key_chain = ['threshold_crossings',index,'threshold_units']  #--> string
        
    if value == id_probabilities:
        key_chain = ['probabilities']                    #--> array of libraries
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
        

    
    
    
