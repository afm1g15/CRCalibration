#ifndef READFILES
#define READFILES

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "Utilities.h"

namespace calib{

  /*
   * @brief Function to read in multiple files from list
   *        and chain the TTrees
   *
   * @param fL     Input file with list of files to read
   * @param nFiles Number of files to read from the input list
   * @param chain  Output chain to return
   */
  void ReadFiles(const char *fileList, int &nFiles, TChain *chain);

  /** 
   * @brief Read a CSV input list and fill relevant vectors
   *
   * @param inputList    The input file containing the CSV list
   * @param nFiles        The number of files to read, nFiles = -1: All 
   * @param prodFiles     The list of files to read TTrees from 
   * @param prodLabels    The label for each production, for naming files
   * @param prodTeXLabels The labels for each production in TeX format, for legends and tables
   *
   */
  void ReadCSV(const std::string &inputList, 
               const int &nFiles, 
               std::vector<std::string> &prodFiles, 
               std::vector<std::string> &prodLabels, 
               std::vector<std::string> &prodTeXLabels);


} //calib (namespace)
#endif
