#include "EventProcessor.h"

namespace calib{
  // ------------------------------------------------------------------------------------------------------------------
  EventProcessor::EventProcessor(const std::vector<TString> &allowedBranches, const std::string &inputList, int &nFiles) :
  fAllowed(allowedBranches),
  fInputList(inputList.c_str()),
  fFiles(nFiles) {}
  // ------------------------------------------------------------------------------------------------------------------
  void EventProcessor::Initialize() {
  
    // First, setup the TChain for writing
    TChain *anachain = new TChain("analysistree/anatree");

    // Input list is a .txt or .list file with a root ana file per line
    // Read these in and chain the analysistree TTrees
    // Then, read in the file list and fill the chain
    ReadFile(fInputList, fFiles, anachain);

    // Finally, allocate the contents to an anatree object
    anatree* evt = new anatree(anachain);

    // Setup member variables
    fEvent = evt;
    fTree  = anachain;

    // Now setup the allowed branches
    anachain->SetBranchStatus("*", 0);
    AnaTree::AllowBranches(anachain, fAllowed);
    anachain->SetMakeClass(1);
    
    // Get the events to loop over
    int nEvents = anachain->GetEntries();
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
  bool EventProcessor::SelectTrack(anatree *evt, int iTrk) const{
    return (evt->trklen_pandoraTrack[iTrk] >= 300.); // at least 3-m long track (based on distributions) 
    //return (evt->trklen_pandoraTrack[iTrk] >= 200.); // at least 2-m (historical value) 
  } // Select Track
  // ------------------------------------------------------------------------------------------------------------------
  bool EventProcessor::SelectHit(anatree *evt, int iTrk, int iPlane, int iHit) const{
    return true;
  } // Select Hit
  // ------------------------------------------------------------------------------------------------------------------
  int EventProcessor::WhichTPC(double x) {
    int i = 0;
    for (; i < 4; ++i) {
      int iapa = (i+1)/2;
      int icpa = i/2;
      double xapa = APA_X_POSITIONS[iapa];
      double xcpa = CPA_X_POSITIONS[icpa];
      if ( (x > xapa && x < xcpa) || (x < xapa && x > xcpa) )
        break;
    }

    return i;
  }
  // ------------------------------------------------------------------------------------------------------------------
}// calib
