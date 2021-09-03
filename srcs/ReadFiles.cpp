#include "ReadFiles.h"

namespace calib{
  
  //-----------------------------------------------------------------------------------------------------------------
  void ReadFile(const char *fileList, int &nFiles, TChain *chain){

    std::ifstream fL(fileList);
    std::string line;
    if(fL){
      // Determine if we want to use all files
      bool useAll = false;
      if(nFiles == -1) useAll = true;

      int nfile = 0;
      while(getline(fL,line) && ((!useAll && nfile < nFiles) || useAll)){
        if(!((nfile+1) % 10))
          std::cout << " Opening file number: " << nfile+1 << std::endl;
        if (line != ""){
          chain->Add(line.c_str());
          if(!((nfile+1) % 10))
            std::cout << " Added " << chain->GetEntries() << " entries so far" << std::endl;

        }
        nfile++;
      }
      // Update the number of files to reflect the total if we've asked for that via nFiles = -1
      if(nFiles == -1)
        nFiles = nfile;
      fL.close();
    }
    std::cout << "-----------------------------------------------------------" << std::endl;
    std::cout << " Read in the files and constructed the chain" << std::endl;
  } // ReadFile

  //-----------------------------------------------------------------------------------------------------------------
} // calib (namespace)
