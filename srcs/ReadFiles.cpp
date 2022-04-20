#include "ReadFiles.h"

namespace calib{
  
  //-----------------------------------------------------------------------------------------------------------------
  
  void ReadFiles(const char *fileList, int &nFiles, TChain *chain){

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
  
  // --------------------------------------------------------------------------------------------------------------------------------------------------

  void ReadCSV(const std::string &inputList, 
               const int &nFiles, 
               std::vector<std::string> &prodFiles, 
               std::vector<std::string> &prodLabels, 
               std::vector<std::string> &prodTeXLabels,
               const bool verbose){
  
    // Read each line of the CSV and then setup the EventProcessor as normal for each one
    std::ifstream iL(inputList);
    std::string line;
    if(iL){
      // Determine if we want to use all files
      bool useAll = false;
      if(nFiles == -1) useAll = true;

      int nfile = 0;
      while(std::getline(iL,line) && ((!useAll && nfile < nFiles) || useAll)){
        if(!((nfile+1) % 10))
          std::cout << " Opening file number: " << nfile+1 << std::endl;
        // Check that the line isn't commented out (with #)
        if (line.at(0) == '#') {
          if(verbose)
            std::cout << __LINE__ << ". The line: " << line << " has been commented out, therefore skipping." << std::endl << std::endl;
          continue;
        }
        if (line != ""){ // Make sure the line isn't empty

          // Remove whitespace and replace tildas with whitespace
          std::string newLine(line);
          newLine.erase(remove(newLine.begin(), newLine.end(), ' '), newLine.end());
          if(verbose){
            std::cout << __LINE__ << ". Removed whitespace, string was: " << line << std::endl;
            std::cout << " and is now: " << newLine << std::endl << std::endl;
          }
          FindReplace(newLine, "~", " ", verbose);
          if(verbose){
            std::cout << __LINE__ << ". Replaced tildas with whitespace, string is now: " << newLine << std::endl << std::endl;
          }
          
          // Now separate the comma-separated string and push to vectors
          std::istringstream ssLine(newLine);
          std::vector<std::string> elements; 
          while (ssLine){
            std::string el;
            if(!std::getline(ssLine,el,',')) 
              break;
            elements.push_back(el);
          }
          // From elements, entry
          //    0: The label
          //    1: The TeX label
          //    2: The input file
          prodLabels.push_back(elements.at(0));
          prodTeXLabels.push_back(elements.at(1));
          prodFiles.push_back(elements.at(2));

        } // Breaking down each line
      } // Read line 
    } // Read list 
  } // ReadCSV

  // --------------------------------------------------------------------------------------------------------------------------------------------------
} // calib (namespace)
