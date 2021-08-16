Usage
=====

operational_sep_quantities.py VERSION 3.0

2021-08-13


.. _installation:

Installation
------------

Codes may be downloaded from repository: https://github.com/ktindiana/operational-sep 

These codes are written in python3, so if necessary, install python3 on your system.
If you do not have the following libraries, download:

.. code-block:: console

   pip3 install matplotlib
   pip3 install wget
   pip3 install scipy
   pip3 install certifi

In Windows, you may need to run python3 from the Anaconda environment. This may require "pip install" or "conda install" rather than "pip3 install" example above.

Additional libraries may be needed.


Running the Code
----------------

Note that these examples were developed using a Mac. The exact call may be slightly different on a Windows machine.

Run the code for an SEP event with the "native" GOES-13 integral fluxes and show the plots on screen:

.. code-block:: console

   python3 operational_sep_quantities.py --StartDate 2012-05-17 --EndDate '2012-05-19 12:00:00' --Experiment GOES-13 --FluxType integral --showplot


Run the code for an SEP event using a user-input time profile data set. First, modify the appropriate values in library/global_vars.py. Then call the code e.g. as below. Note that the filename containing the user data is with respect to the "data" directory. So a file called MyFluxes.txt should be in data/MyFluxes.txt.

.. code-block:: console

   python3 operational_sep_quantities.py --StartDate 2012-05-17 --EndDate '2012-05-19 12:00:00' --Experiment user --ModelName MyModel --UserFile MyFluxes.txt --FluxType integral --showplot


Run the code to perform background subtraction and apply the Sandberg et al. (2014) and Bruno (2017) effective energies to the GOES bins. (note: cannot bg-subtract GOES integral fluxes), e.g.:

.. code-block:: console

   python3 operational_sep_quantities.py --StartDate 2012-05-17 --EndDate '2012-05-19 12:00:00' --Experiment GOES-13 --FluxType differential  --showplot --options uncorrected,S14,Bruno2017 --SubtractBG --BGStartDate 2012-05-10 --BGEndDate --2012-05-17


Run the code as an imported module:

.. code-block:: console

   import operational_sep_quantities as sep
   start_date = '2012-05-17'
   end_date = '2012-05-19 12:00:00'
   experiment = 'GOES-13'
   flux_type = 'integral'
   spase_id = ''
   model_name = '' #if experiment is user, set model_name to describe data set
   user_file = '' #if experiment is user, specify filename containing fluxes
   showplot = True  #Turn to False if don't want to see plots
   saveplot = False #turn to true if you want to save plots to file
   options = '' #various options: S14, Bruno2017, uncorrected
   doBGSub = False #Set true if want to perform background subtraction
   bgstart_date = "2012-05-10" #Dates used to estimate mean background if
   bgend_date = "2012-05-17"   #doBGSub is set to True
   detect_prev_event = True  #Helps if previous event causes high intensities
   two_peaks = False  #Helps if two increases above threshold in one event
   umasep = False #Set to true if want UMASEP values (see explanation above)
   threshold = '' #Add a threshold to 10,10 and 100,1: '30,1' or '4.9-7.3,0.01'
   nointerp = False #Default False; set to True to stop linear interpolation in time

   sep_year, sep_month,sep_day, jsonfname = sep.run_all(start_date, \
        end_date, experiment, flux_type, model_name, user_file,\
        spase_id, showplot, saveplot, detect_prev_event,  \
        two_peaks, umasep, threshold, options, doBGSub, bgstart_date, \
        bgend_date,nointerp)
