# CRCalibration
A repository to hold my cosmic ray calibration work at the University of Sheffield on the DUNE experiment

The work in this directory is written for the CR sample construced by the University of Sheffield DUNE calibration group

There also exists a version of this data on the DUNE GPVM machines, but the file list will need to be tweaked accordingly

---------------------------------------------------------------------------------------------------------

# Sample Content Studies
First studies will look at the sample contents 

## Directory structure
- `srcs` directory contains all utilities and main macros for the studies
- `run_` files are setup to build the relevant sources for the main macro you wish to run 

## Running
To run the macro for the sample contents study, use

    bash$ root -l -b run_fileContentStudies.C
    root [0] functionName("config/configName.txt")

    where functionName is the name of the main function in the srcs file
    and configName is the appropriate settings file for this analysis 

where:

- `inList`: A .txt or .list file which houses the list of input files to read
- `nFiles`: Number of such files to read over

