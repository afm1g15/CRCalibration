#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TMath.h"

#include <iostream>
#include <iomanip>
#include <string>

void stoppingPowerPlots(const char *input_list){

  // Input list is a .txt or .list file with a root ana file per line
  // Read these in and chain the analysistree TTrees
  // 
  // First, setup the input file list for reading
  std::string line;
  std::ifstream inlist(input_list);

  // Then setup the TChain for writing
  TChain anachain("analysistree/anatree");

  // Finally setup the histograms, counters and any other variables to add to
  TH2D *h_dedx_resrange = new TH2D("h_dedx_resrange","",100,0,40,100,0,40);

  // And run the code
  if(inlist){
    unsigned int nfiles = 0;
    while(getline(inlist,line) && nfiles < 20){
      std::cout << " Opening file number: " << nfiles+1 << std::endl;
      if (line != ""){
        anachain.Add(line.c_str());
        std::cout << " Number of entries (events) added so far: " << anachain.GetEntries() << std::endl;

      }
      nfiles++;
    }
    inlist.close();
  }
}
