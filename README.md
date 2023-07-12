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

## Setup

Depending on the form of the analysis chosen, there are relevant setup scripts.
You will have to modify srcs/Setup/SetupUtils.h and srcs/Setup/Utilities.h to include the correct script

- For a DUNE FD analysis, with samples generated using the AnalysisTree in duneana/AnaTree,
  use the standard setup script: srcs/Setup/Setup.h
- For a PDDUNE analysis, with samples generated using the AnalysisTree in protoduneana/singlephase/Analysis/PDSPAnalyzer,
  use the standard setup script: srcs/Setup/PDSPAnaSetup.h

## Running

To run the macro for the sample contents study, use

    bash$ root -l -b run/run_macroName.C
    root [0] functionName("config/configName.txt")

    
where 

- `macroName` is the name of the run script associated to the relevant macro in `srcs`
- `functionName` is the name of the main function in the `srcs` macro
- `configName.txt` is the appropriate settings file for this analysis 

## Main scripts

- `srcs/ActivityStudies.cpp` produces the initial plots for the CR muon absolute energy calibration
- `srcs/DeltaStudies.cpp` produces the same plots tuned using cuts on delta-rays
- `srcs/SliceAndFit.cpp` slices up the 2D plots from the above two scripts and determines the functional form of dQ/dx w.r.t the desired parameter

Each of the above have corresponding configuration files with similar names located here: `config/...` or here `config/slices/v09410002/..`

The remaining scripts were largely written for analysis interest while the calibration procedure was developed.

## Contact

Rhiannon S. Jones (she/her)

University of Sheffield

Email: r.s.jones@sheffield.ac.uk

Date first setup: July 2021
