/************************************************************************
 *
 * A Macro to plot various distrubtions in order to 
 * understand the contents of the CR sample for 
 * through-going muons
 *
 * Example file list located here:
 *   /home/jones/work/cosmics/LArSoft-v08_50_00/work/files/anafiles.list
 *
 * Parameters to look into:
 *   - Muon momentum
 *   - ThetaXZ
 *   - ThetaYZ
 *   - Y vs Z position
 *   - X vs Y position
 *   - Entry plane
 *   - Exit plane
 *   - Hit distribution in YZ
 *   - Hist distribution in XY
 *
 * Ultimately, compare with Viktor and Praveen's analogous studies
 *
 *************************************************************************/

#include "setup.h"
#include "EventProcessor.h"

using namespace calib;

// Allowed branches to read from the tree
std::vector<TString> allowed = {
   "run",
   "event",
   "ntracks_pandoraTrack",
   "ntrkhits_pandoraTrack",
   "trkdqdx_pandoraTrack",
   "trkdedx_pandoraTrack",
   "trkxyz_pandoraTrack",
   "trklen_pandoraTrack",
 };

int fileContentStudies(const int n = -1, const char *input_list="/home/jones/work/cosmics/LArSoft-v08_50_00/work/files/anafiles.list"){

  // First, setup timing information so we can monitor the run
  time_t rawtime;
  std::cout << "-----------------------------------------------------------" << std::endl;
  GetTime(rawtime);
  std::cout << "-----------------------------------------------------------" << std::endl;

  // Then setup the histograms, counters and any other variables to add to
  TH2D *h_dedx_resrange = new TH2D("h_dedx_resrange","",100,0,40,100,0,40);

  // Setup TTree from input file list
  std::cout << " Reading files and filling tree" << std::endl;
  
  EventProcessor evtProc(allowed, input_list, n);
  evtProc.Initialize();

  // Now setup the tree and event objects to work with
  TChain *tree = evtProc.GetTree();
  anatree *evt = evtProc.GetEvents();
  
  // Start of analysis (loop over chain and events
  std::cout << " Running analysis..." << std::endl;

  // Now loop over the events
  unsigned int nEvts = tree->GetEntries();
  for(unsigned int i = 0; i < nEvts; ++i){
    tree->GetEntry(i);
    if(!evtProc.SelectEvent(evt)) continue;
    unsigned int nTrks = evt->ntracks_pandoraTrack;
    std::cout << " Event: " << i << "/" << nEvts << std::endl;
    std::cout << " Number of tracks : " << nTrks << std::endl;
  }// Event loop

  // End of script
  std::cout << " ...finished analysis" << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;
  time_t rawtime_end;
  GetTime(rawtime_end);
  std::cout << "-----------------------------------------------------------" << std::endl;
  GetTotalTime(rawtime, rawtime_end);
  std::cout << "-----------------------------------------------------------" << std::endl;
  
  return 0;
}
