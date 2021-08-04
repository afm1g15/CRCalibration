#include "readFiles.h"

namespace calib{
  
  //-----------------------------------------------------------------------------------------------------------------
  void ReadFile(const char *fileList, const int &nFiles, TChain *chain){

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
            std::cout << " Number of entries (events) added so far: " << chain->GetEntries() << std::endl;

        }
        nfile++;
      }
      fL.close();
    }
    std::cout << "-----------------------------------------------------------" << std::endl;
    std::cout << " Read in the files and constructed the chain" << std::endl;
  } // ReadFile

  //-----------------------------------------------------------------------------------------------------------------
} // calib (namespace)
