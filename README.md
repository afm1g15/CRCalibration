# CRCalibration
A repository to hold cosmic ray calibration work at the University of Sheffield on the DUNE experiment

The work in this directory is written for the CR sample construced by the University of Sheffield DUNE calibration group

---------------------------------------------------------------------------------------------------------

## Directory structure
- `srcs` directory contains all utilities and main macros for the studies
- `run` directory contains files which are setup to build the relevant sources for the main macro you wish to run
- `config` directory contains the configuration files for various versions of each analysis
- `docs` directory containing the doxygen for the codebase
- `python` directory containing scripts to run additional studies, such as writing the MPV dE/dx from theory

## Running
To run the macro for the sample contents study, use

    bash$ root -l -b run/run_macroName.C
    root [0] functionName("config/configName.txt")

    
where 
- `macroName` is the name of the run script associated to the relevant macro in `srcs`
- `functionName` is the name of the main function in the `srcs` macro
- `configName.txt` is the appropriate settings file for this analysis 

