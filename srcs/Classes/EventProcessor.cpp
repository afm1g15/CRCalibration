#include "EventProcessor.h"

namespace calib{
  // ------------------------------------------------------------------------------------------------------------------
  EventProcessor::EventProcessor(const std::vector<TString> &allowedBranches, const std::string &inputList, int nFiles) :
  fAllowed(allowedBranches),
  fInputList(inputList.c_str()),
  fFiles(nFiles) {}
  // ------------------------------------------------------------------------------------------------------------------
  void EventProcessor::Initialize(const bool multiList) {
  
    // If we want to initialise a single file, do so
    // If not, do the usual TChain
    if(!multiList){
      TFile *f = new TFile(fInputList);
      TChain *t = static_cast<TChain*>(f->Get("analysistree/anatree"));
      anatree *evt =  new anatree(t);
      
      // Setup member variables
      fEvent = evt;
      fTree  = t;
    }
    else{
      // First, setup the TChain for writing
      TChain *anachain = new TChain("analysistree/anatree");

      // Input list is a .txt or .list file with a root ana file per line
      // Read these in and chain the analysistree TTrees
      // Then, read in the file list and fill the chain
      ReadFiles(fInputList, fFiles, anachain);

      // Finally, allocate the contents to an anatree object
      anatree* evt = new anatree(anachain);

      // Setup member variables
      fEvent = evt;
      fTree  = anachain;
    }

    // Now setup the allowed branches
    fTree->SetBranchStatus("*", 0);
    AnaTree::AllowBranches(fTree, fAllowed);
    fTree->SetMakeClass(1);

    // Get the events to loop over
    int nEvents = fTree->GetEntries();
    std::cout << std::endl;
    std::cout << " Total number of events in the chain: " << nEvents << std::endl;
    std::cout << "-----------------------------------------------------------" << std::endl;    

  } // Initialize
  // ------------------------------------------------------------------------------------------------------------------
  void EventProcessor::Finalize() {
  } // Finalize
  // ------------------------------------------------------------------------------------------------------------------
  int EventProcessor::GetFiles() const{
    return fFiles;
  }
  // ------------------------------------------------------------------------------------------------------------------
  anatree *EventProcessor::GetEvents() const{
    return fEvent;
  }
  // ------------------------------------------------------------------------------------------------------------------
  TChain *EventProcessor::GetTree() const{
    return fTree;
  }
  // ------------------------------------------------------------------------------------------------------------------
  bool EventProcessor::SelectEvent(anatree *evt) const{
    return true;
  } // Select Event
  // ------------------------------------------------------------------------------------------------------------------
  bool EventProcessor::SelectTrack(double length, int iTrk) const{
    //return (evt->trklen_pandoraTrack[iTrk] >= 300.); // at least 3-m long track (based on distributions) 
    return (length >= 200.); // at least 2-m (historical value) 
  } // Select Track
  // ------------------------------------------------------------------------------------------------------------------
  bool EventProcessor::SelectHit(anatree *evt, int iTrk, int iPlane, int iHit) const{
    return true;
  } // Select Hit
  // ------------------------------------------------------------------------------------------------------------------
  int EventProcessor::WhichTPC(double x, bool PD) {
    
    unsigned int i = 0;
    unsigned int nTPCs = 4;
    std::vector<float> APAs, CPAs;
    if(PD){
      nTPCs = 2;
      APAs.insert(APAs.begin(), PD_APA_X_POSITIONS, PD_APA_X_POSITIONS + std::size(PD_APA_X_POSITIONS));
      CPAs.insert(CPAs.begin(), PD_CPA_X_POSITIONS, PD_CPA_X_POSITIONS + std::size(PD_CPA_X_POSITIONS));
    }
    else{
      APAs.insert(APAs.begin(), APA_X_POSITIONS, APA_X_POSITIONS + std::size(APA_X_POSITIONS));
      CPAs.insert(CPAs.begin(), CPA_X_POSITIONS, CPA_X_POSITIONS + std::size(CPA_X_POSITIONS));
    }
    
    for (; i < nTPCs; ++i) {
      int iapa = (i+1)/2;
      int icpa = i/2;
      double xapa = APAs.at(iapa);
      double xcpa = CPAs.at(icpa);
      if ( (x > xapa && x <= xcpa) || (x <= xapa && x > xcpa) )
        break;
    }

    return i;
  }
  // ------------------------------------------------------------------------------------------------------------------
}// calib
