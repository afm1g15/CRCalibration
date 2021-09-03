#ifndef READFILES
#define READFILES

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"

#include <string>
#include <vector>

#include "Utilities.h"

namespace calib{

  /*
   * @brief Function to read in multiple files from list
   *        and chain the TTrees
   *
   * @param fL  Input file with list of files to read
   * @param chain  Output chain to return
   */
  void ReadFile(const char *fileList, int &nFiles, TChain *chain);

} //calib (namespace)
#endif
