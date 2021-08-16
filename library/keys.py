from library import global_vars as vars

__version__ = "0.2"
__author__ = "Katie Whitman"
__maintainer__ = "Katie Whitman"
__email__ = "kathryn.whitman@nasa.gov"

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

def about_keys():
    """ About keys.py
    
        Define the keys used to call different variables in the CCMC
        and observation JSON files. These are the keys that are different
        between the observations json and the model output json.
    """
