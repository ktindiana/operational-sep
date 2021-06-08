# operational-sep
**All files needed to run operational_sep_quantities.py:**\
operational_sep_quantities.py\
library/ccmc_json_handler.py\
library/derive_background.py\
library/global_vars.py\
library/keys.py\
library/model_template.json\
library/observations_template.json\
library/read_datasets.py\

Calculate solar energetic particle (SEP) proton flux quantities relevant to space radiation operations. The goal is for this code to be a robust, user-friendly code written in python3. Works for only one SEP event at a time.  Check back for updates, as this code will be modified as it is tested for a variety of SEP events. Please send bug reports and feedback to kathryn.whitman@nasa.gov. 

**Note that changes for each version of the code are logged and described in comments at the top of operational_sep_quantities.py and the other supporting python codes.**

This code outputs derived quantities into csv files and the JSON format used by CCMC for the SEP Scoreboard.

This code was originally developed in support of the SHINE 2019 SEP modeling challenge session to assist SEP modelers in calculating and reporting the quantities in this code. The code continues to be developed and expanded for the purpose of SEP event analysis and SEP model validation.

Users may specify any of the native data sets or supply their own flux time profiles. If differential fluxes are used, the >10 and >100 MeV integral fluxes will be estimated from the available differential channels.

Data sets that are native to this code currently include GOES-08 to GOES-15 corrected or uncorrected fluxes, SEPEM RSDv2, and SOHO/EPHIN Level3 data. Native data sets will continue to be expanded and various calibration schemes from the literature will be made available. 

An option has been added to perform a background subtracion of the SEP flux. The code for this is contained in derive_background.py, which is automatically called from operational_sep_quantities.py if background subtraction is specified.

For GOES fluxes, users have the option to apply Sandberg et al. (2014) effective energies for EPS or EPEAD. The effective energies may be applied to GOES-08 to GOES-11. More recent GOES spacecraft will apply the GOES-11 effective energies. Users may also apply Bruno (2017) effective energies to GOES-13 or GOES-15 channels P6 - P11.

The code calculates:
- Values based on the operational thresholds: >10 MeV exceeds 10 pfu, >100 MeV exceeds 1 pfu
- Values based on user input thresholds which can be applied to both integral or differential channels
- Start Time when thresholds are crossed
- End Time: when flux falls below 0.85 x threshold values
- Onset Peak Flux - the onset peak is defined to be the flux value at the location that the initial intensity rise turns over and becomes more gradual; a preliminary algorithm has been developed to locate the onset peak in an automated fashion and is still under development.
- Time of Onset Peak
- Maximum Flux for >10 MeV, >100 MeV and any other channel for which a threshold has been applied: Maximum flux value between the start and end times; this may be the ESP for lower energy channels
- Time of Maximum Flux
- Duration - defined as end time - start time
- Event fluence (total integrated event intensity) for >10 MeV, >100 MeV in cm-2
- Event fluence spectrum in [cm-2 sr]
- Plots of time profiles for >10 MeV, >100 MeV, any user defined-thresholds and fluence spectrum
- Plots of time profiles and their derivatives as part of the identification of the onset peak

## Run code from command line as, e.g.:
python3 operational_sep_quantities.py --StartDate 2012-05-17 --EndDate 2012-05-20 --Experiment GOES-13 --FluxType integral --showplot

**For a start time other than midnight:**\
python3 operational_sep_quantities.py --StartDate "2012-01-27 16:00:00" --EndDate 2012-02-02 --Experiment GOES-13 --FluxType integral --showplot

**OR if the flux is already above threshold at the very first point due to a previous event, then falls below threshold prior to the start of the desired event:**\
python3 operational_sep_quantities.py --StartDate 2012-05-17 --EndDate 2012-05-20 --Experiment GOES-13 --FluxType integral --showplot --DetectPreviousEvent

**OR if the event has an initial increase above threshold for a few points, drops below threshold, then increases again above threshold for the remainder of the event:**\
python3 operational_sep_quantities.py --StartDate 2012-05-17 --EndDate 2012-05-20 --Experiment GOES-13 --FluxType integral --showplot --TwoPeaks

**Input a user-defined threshold of >30 MeV exceeds 1 pfu:**\
python3 operational_sep_quantities.py --StartDate 2012-05-17 --EndDate 2012-05-20 --Experiment GOES-13 --FluxType integral --showplot --Threshold 30,1

**Input multiple user-defined threshold of >30 MeV exceeds 1 pfu, >50 MeV exceeds 1pfu, >5 MeV exceeds 100 pfu:**\
python3 operational_sep_quantities.py --StartDate 2012-05-17 --EndDate 2012-05-20 --Experiment GOES-13 --FluxType integral --showplot --Threshold "30,1;50,1;5,100"

**Input a user-defined threshold for a differential channel. Must specify both edges of the bin. The thresholds specifies that the flux in the 40.9 - 53 energy bin exceeds 0.001 [MeV-1 cm-2 s-1 sr-1]:**\
python3 operational_sep_quantities.py --StartDate 2012-05-17 --EndDate 2012-05-20 --Experiment EPHIN --FluxType differential --showplot --Threshold "40.9-53,0.001"

**Perform background subtraction and apply Sandberg et al. (2014) and Bruno (2017) effective energies to the GOES bins. (note: cannot bg-subtract GOES integral fluxes), e.g.:**\
    python3 operational_sep_quantities.py --StartDate 2012-05-17 --EndDate "2012-05-19 12:00:00" --Experiment GOES-13 --FluxType differential  --showplot --options uncorrected,S14,Bruno2017 --SubtractBG --BGStartDate 2012-05-10 --BGEndDate --2012-05-17


## Run code from command line for user-input file as, e.g.:
python3 operational_sep_quantities.py --StartDate 2012-05-17 --EndDate "2012-05-19 12:00:00" --Experiment user --ModelName MyModel --UserFile MyFluxes.txt --FluxType integral --showplot

## Import code and run as, e.g.:
    import operational_sep_quantities as sep
    start_date = '2012-05-17'
    end_date = '2012-05-19 12:00:00'
    experiment = 'GOES-13'
    flux_type = 'integral'
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
    threshold = '100,1' #default; modify to add a threshold to 10,10 and 100,1

    FirstStart, LastEnd, ShortEvent, LateHundred, sep_year, sep_month, \
    sep_day, jsonfname = sep.run_all(start_date, end_date, experiment, flux_type, \
        model_name, user_file, showplot, saveplot, detect_prev_event,  \
        two_peaks, umasep, threshold, options, doBGSub, bgstart_date, \
        bgend_date)


## Full documentation for operational_sep_quantities.py is located in the documentation folder:
### NAME
    operational_sep_quantities

### AUTHOR
    Katie Whitman

### FILE
operational-sep/operational_sep_quantities.py

### DATA
    __email__ = 'kathryn.whitman@nasa.gov'
    __maintainer__ = 'Katie Whitman'

### VERSION
V2.6

### FUNCTIONS
#### all_program_info()
    This program will calculate various useful pieces of operational
    information about SEP events from GOES-08, -10, -11, -12, -13, -14, -15
    data, SOHO/EPHIN Level3 and the SEPEM (RDSv2 and RDSv3) datasets.
    
    SEP event values are always calculated for threshold definitions:
        >10 MeV exceeds 10 pfu
        >100 MeV exceed 1 pfu

    The user may add multiple additional thresholds through the command line.
    This program will check if data is already present in a 'data' directory. If
    not, GOES or EPHIN data will be automatically downloaded from the web. SEPEM
    (RSDv2 and RSDv3) data must be downloaded by the user and unzipped inside
    the 'data' directory. Because the SEPEM data set is so large (every 5
    minutes from 1974 to 2015 for RSDv2 and to 2017 for RSDv3), the program will
    break up the data into yearly files for faster reading.

    Data sets are checked for bad data point (negative or None value fluxes)
    and the default behavior is to fill in those bad data points by
    performing a linear interpolation with time. This choice was made to
    calculate more accurate event-intergrated fluence values from data.
    Interpolation with time is not appropriate for model predictions, as
    it will inherently change the prediction or may not be desired by
    the user for the data set. Turn off interpolation with time by
    setting the --NoInterp flag (or nointerp=True). If the interpolation
    is turned off, negative flux values will be set to None.
    Zeroes are always treated as valid values and are not replaced.
    If there are gaps in the time steps, the code does NOT try to
    fill in the gaps. It will ONLY perform interpolation for time
    steps present in the input data set. i.e. gaps in time are not
    interpolated, only time steps with negative or None flux values.

    The values calculated here are important for space radiation operations:
        Onset time, i.e. time to cross thresholds
        Onset peak intensity
        Onset peak time
        Maximum intensity
        Time of maximum intensity
        Rise time (onset to peak)
        End time, i.e. fall below 0.85*threshold for 3 points (15 mins for GOES)
        Duration
        Event-integrated fluences
        Proton fluxes at various times after threshold crossing (UMASEP option)

    User may choose differential proton fluxes (e.g. [MeV s sr cm^2]^-1) or
    integral fluxes (e.g. [s sr cm^2]^-1 or pfu). The program has no internal
    checks orrequirements on units - EXCEPT FOR THE THRESHOLD DEFINITIONS
    OF >10, 10 and >100, 1. If you convert those thresholds in the main program
    to your units, you should be able to generate consistent results.
    Also, all of the plots and messages refer to MeV, pfu, and cm. Change those
    labels everywhere if you choose other units. Currently no features to change
    units automatically.

    User may specify various options, that currently only apply to GOES data:
        Choose corrected or uncorrected GOES fluxes.
        Choose to apply Bruno (2017) or Sandberg et al. (2014) effective
        energies to GOES uncorrected data.
        --options uncorrected
        --options uncorrected,S14,Bruno2017 (recommend using background subtraction)

    Users may choose to perform a background subtraction by specifying:
        --SubtractBG --BGStartDate YYYY-MM-DD --BGEndDate YYYY-MM-DD
    The user should look at the data and select an appropriate time frame
    prior to the event when the background is calm and well-defined. If
    performing background subtraction, the mean background will be subtracted
    from the fluxes in the SEP time frame (StartDate to EndDate). Plots
    showing the mean background level, the background flux only, and the
    background-subtracted SEP fluxes will be created by derive_background to
    verify the quality of the background estimation and subtraction.

    If a previous event is ongoing and the specified time period starts with a
    threshold already crossed, you may try to set the --DetectPreviousEvent
    flag. If the flux drops below threshold before the next event starts, the
    program will identify the second event. This will only work if the
    threshold is already crossed for the very first time in your specified
    time period, and if the flux drops below threshold before the next event
    starts.

    If the event has an initial increase above threshold for a few points, falls
    below threshold, then continues to increase above threshold again, you
    may try to use the --TwoPeak feature to capture the full duration of the
    event. The initial increase above threshold must be less than a day. An
    example of this scenario can be seen in >100 MeV for 2011-08-04.

    A flag was added with the "UMASEP" option. When this flag is used
    (--UMASEP), the code finds all information for four energy channels and
    thresholds used by UMASEP: >10 MeV, 10 pfu; >100 MeV, 1 pfu; >30 MeV, 1 pfu;
    >50 MeV, 1 pfu. The proton flux in each of these channels is reported for
    multiple times after threshold crossing (Ts). The applied time delays are
    as follows:
        >10 MeV - Ts + 3, 4, 5, 6, 7 hours
        >30 MeV - Ts + 3, 4, 5, 6, 7 hours
        >50 MeV - Ts + 3, 4, 5, 6, 7 hours
        >100 MeV - Ts + 3, 4, 5, 6, 7 hours

    RUN CODE FROM COMMAND LINE (put on one line), e.g.:
    python3 operational_sep_quantities.py --StartDate 2012-05-17 --EndDate "2012-05-19 12:00:00" --Experiment GOES-13 --FluxType integral --showplot --saveplot

    RUN CODE FROM COMMAND FOR USER DATA SET (put on one line), e.g.:
    python3 operational_sep_quantities.py --StartDate 2012-05-17 --EndDate "2012-05-19 12:00:00" --Experiment user --ModelName MyModel --UserFile MyFluxes.txt --FluxType integral --showplot

    RUN CODE FROM COMMAND LINE AND PERFORM BACKGROUND SUBTRACTION AND APPLY
    Sandberg et al. (2014) and Bruno (2017) effective energies to the GOES bins.
    (note: cannot bg-subtract GOES integral fluxes), e.g.:
    python3 operational_sep_quantities.py --StartDate 2012-05-17 --EndDate "2012-05-19 12:00:00" --Experiment GOES-13 --FluxType differential  --showplot --options uncorrected,S14,Bruno2017 --SubtractBG --BGStartDate 2012-05-10 --BGEndDate --2012-05-17

    RUN CODE IMPORTED INTO ANOTHER PYTHON PROGRAM, e.g.:
    import operational_sep_quantities as sep
    start_date = '2012-05-17'
    end_date = '2012-05-19 12:00:00'
    experiment = 'GOES-13'
    flux_type = 'integral'
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
    threshold = '100,1' #default; modify to add a threshold to 10,10 and 100,1

    FirstStart, LastEnd, ShortEvent, LateHundred, sep_year, sep_month, \
    sep_day, jsonfname = sep.run_all(start_date, end_date, experiment, flux_type, \
    model_name, user_file, showplot, saveplot, detect_prev_event,  \
    two_peaks, umasep, threshold, options, doBGSub, bgstart_date, \
    bgend_date)

    Set the desired directory locations for the data and output at the beginning
    of the program in datapath and outpath. Defaults are 'data' and 'output'.

    In order to calculate the fluence, the program determines time_resolution
    (seconds) by finding the difference between every consecutive set of
    time points in the data set. The most common difference is identified as
    the time resolution. This method should find an accurate time resolution
    even if there are gaps in the time steps.
    If the time steps in the data set are truly irregular, the user will
    have to manually set the time resolution inside the subroutine
    calculate_fluence.

    OUTPUT: This program outputs 3 to 4 files, 1 per defined threshold plus
    a summary file containing all of the values calculated for each threshold.
    A file named as e.g. fluence_GOES-13_differential_gt10_2012_3_7.csv contains
    the event-integrated fluence for each energy channel using the specified
    threshold (gt10) to determine start and stop times.
    A file named as e.g. sep_values_GOES-13_differential_2012_3_7.csv contains
    start time, peak flux, etc, for each of the defined thresholds.

    The program write to file the >10 MeV and >100 MeV time series for the
    date range input by the user. If the original data were integral fluxes,
    then the output files simply contain the >10 and >100 MeV time series from
    the input files. If the original data were differential fluxes, then the
    estimated >10 and >100 MeV fluxes are output as time series.

    USER INPUT DATA SETS: Users may input their own data set. For example, if an
    SEP modeler would like to feed their own intensity time series into this
    code and calculate all values in exactly the same way they were calculated
    for data, it is possible to do that. Fluxes should be in units of
    1/[MeV cm^2 s sr] or 1/[cm^2 s sr] and energy channels in MeV for the plot
    labels to be correct. You can use any units, as long as you are consistent
    with energy units in energy channel/bin definition and in fluxes and you
    MODIFY THE THRESHOLD VALUES TO REFLECT YOUR UNITS. You may then want to
    modify plot labels accordingly if not using MeV and cm.
    NOTE: The first column in your flux file is assumed to be time in format
    YYYY-MM-DD HH:MM:SS. IMPORTANT FORMATTING!!
    NOTE: The flux file may contain header lines that start with a hash #,
    including blank lines.
    NOTE: Any bad or missing fluxes must be indicated by a negative value.
    NOTE: Put your flux file into the "datapath" directory. Filenames will be
    relative to this path.
    NOTE: Please use only differential or integral channels. Please do not mix
    them. You may have one integral channel in the last bin, as this is the way
    HEPAD works and the code has been written to include that HEPAD >700 MeV
    bin along with lower differential channels.

    USER VARIABLES: The user must modify the following variables in
    library/global_vars.py:
        user_col - identify columns in your file containing fluxes to analyze;
            even if your delimeter is white space, consider the date-time
            column as one single column. SET IN library/global_vars.py.
        user_delim - delimeter between columns, e.g. " " or ","   Use " " for
            any amount of whitespace. SET IN library/global_vars.py.
        user_energy_bins - define your energy bins at the top of the code in the
            variable user_energy_bins. Follow the format in the subroutine
            define_energy_bins. SET IN library/global_vars.py.
        user_fname - specify the name of the file containing the fluxes
            through an argument in the command line. --UserFile  The
            user_fname variable will be updated with that filename. ARGUMENT
        time_resolution - the program determines time_resolution
            (seconds) by finding the difference between every consecutive
            set of time points in the data set. The most common difference
            is identified as the time resolution. This method should find
            an accurate time resolution even if there are gaps in the
            time steps.. AUTOMATICALLY DETERMINED.
