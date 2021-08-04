# CRCalibration
A repository to hold my cosmic ray calibration work at the University of Sheffield on the DUNE experiment

---------------------------------------------------------------------------------------------------------

# Sample Content Studies
First studies will look at the sample contents 

## Directory structure
- `srcs` directory contains all utilities and main macros for the studies
- `run_` files are setup to build the relevant sources for the main macro you wish to run 

## Running
To run the macro for the sample contents study, use

    bash$ root -l -b run_fileContentStudies.C
    root [0] fileContentStudies(inList,nFiles)

where:

- `inList`: A .txt or .list file which houses the list of input files to read
- `nFiles`: Number of such files to read over

