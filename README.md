# CRCalibration
A repository to hold my cosmic ray calibration work at the University of Sheffield on the DUNE experiment

---------------------------------------------------------------------------------------------------------

# Sample Content Studies
  First studies will look at the sample contents 

## Running
  To run the sample contents study macro use

  bash$ root -l -b run_fileContentStudies.C
  root [0] fileContentStudies(inList,nFiles)

  where:

    - inList: A .txt or .list file which houses the list of input files to read

    - nFiles: Number of such files to read over

