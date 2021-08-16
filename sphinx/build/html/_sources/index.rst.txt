.. operational_sep_quantities documentation master file, created by
   sphinx-quickstart on Fri Aug 13 16:21:38 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to operational_sep_quantities's documentation
=====================================================

operational_sep_quantities CURRENT VERSION 3.0

2021-08-13

**operational_sep_quantities** is a python code to derive quantities from solar energetic particle (SEP) event flux time profiles to use for forecasting, assessment, and validation of SEP models. SEP flux time profiles are assessed using thresholds applied to energy channels. The operational thresholds >10 MeV exceeds 10 pfu and >100 MeV exceeds 1 pfu are always applied. Users may enter additional threshold definitions. If a threshold is crossed, the program derives event start and end time, onset peak flux and time, maximum flux and time, and event-integrated fluence. 

Derived values are output to both csv files and the JSON format and accompanying text files used by CCMC's SEP Scoreboard.

This code has a few "native" data sets, including GOES-08 to GOES-15, SOHO/EPHIN, SOHO/EPHIN data in the energy bins provided on the REleASE website, SEPEM RDSv2.0, and SEPEM RDSv3.0 (when it becomes publicly available). Users may input their own time profile data sets of model output.

This code works with both integral and differential fluxes. Data gaps (time steps with negative with negative fluxes, e.g. -999) are replaced with estimated values via linear interpolation in time. 

Users may choose to perform background-subtraction by specifying an appropriate time period containing the background. For GOES data, users may choose between corrected or uncorrected fluxes, or to apply the Sandberg et al. (2014) or Bruno (2017) effective energies determined for the EPS, EPEAD, and HEPAD detectors.

This project is under active development.
 
.. toctree::
   
   :maxdepth: 2
   :caption: Contents:
   
   repository
   usage
   overall_description
   modules


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
