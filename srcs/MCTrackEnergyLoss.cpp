/************************************************************************
 *
 * Script to run the absolute energy scale calibration procedure using 
 * reconstructed quantites only
 *
 * First defined: August 2022
 * Author:        R. S. Jones
 *
 * Example file list located here:
 *   /home/jones/work/cosmics/LArSoft-v08_50_00/work/files/v09_41_00_02_files.list
 *
 *************************************************************************/

#include "Classes/EventProcessor.h"
#include "Classes/ConfigReader.h"

using namespace calib;
using namespace cppsecrets;

// Allowed branches to read from the tree
std::vector<TString> allowed = {
   "run",
   "event",
   "geant_list_size",
   "inTPCActive",
   "TrackId",
   "pdg",
   "Mother",
   "StartE_tpcAV",
   "EndE_tpcAV",
   "StartPointx_tpcAV",
   "StartPointy_tpcAV",
   "StartPointz_tpcAV",
   "EndPointx_tpcAV",
   "EndPointy_tpcAV",
   "EndPointz_tpcAV",
   "StartPointx",
   "StartPointy",
   "StartPointz",
   "EndPointx",
   "EndPointy",
   "EndPointz",
   "no_mctracks",
   "mctrk_TrackId",
   "mctrk_pdg",
   "nmctrksteps",
   "mctrk_dEdx",
   "mctrk_dQdx",
   "mctrk_Process",
 };

// A translation list from plane labels to longer labels for plotting
std::map<std::string, std::string> planeLabels = {
  {"h0", "APA 1"},
  {"h1", "CPA 1"},
  {"h2", "APA 2"},
  {"h3", "CPA 2"},
  {"h4", "APA 3"},
  {"t",  "Top"},
  {"bo", "Bot."},
  {"f",  "Fro."},
  {"ba", "Back"},
};

typedef std::vector<Plane> PlaneList;

int mctrackEnergyLoss(const char *config){

  // First, setup timing information so we can monitor the run
  time_t rawtime;
  std::cout << "-----------------------------------------------------------" << std::endl;
  GetTime(rawtime);
  std::cout << "-----------------------------------------------------------" << std::endl;

  //------------------------------------------------------------------------------------------
  //                                    Configure
  //------------------------------------------------------------------------------------------
  // Create object of the classConfigReader
  // Parse the configuration file
  // Dump map on the console after parsing it
  ConfigReader* p = ConfigReader::getInstance();
  p->parseFile(config);
  std::cout << " Variables from configuration file: " << std::endl;
  p->dumpFileValues();
  std::cout << "-----------------------------------------------------------" << std::endl;

  // Get configuration variables
  int n = -1;
  int thru = 0;
 
  int nBinsFromPeak   = -1; // How many bins to traverse either side of the peak in the fit
  int nBinsFromPeakL  = -1; // How many bins to traverse left of the peak in the fit
  int nBinsFromPeakR  = -1; // How many bins to traverse right of the peak in the fit
  
  std::string input_list = "";
  std::string location="";
  std::string tag="";
  std::vector<double> minx_fid, miny_fid, minz_fid;
  std::vector<double> maxx_fid, maxy_fid, maxz_fid;
  std::vector<double> minx_av, miny_av, minz_av;
  std::vector<double> maxx_av, maxy_av, maxz_av;

  p->getValue("InputList",       input_list);
  p->getValue("Location",        location);
  p->getValue("Tag",             tag);
  p->getValue("NFiles",          n);
  p->getValue("Thru",            thru);
  p->getValue("NBinsFromPeak",   nBinsFromPeak);
  p->getValue("NBinsFromPeakL",  nBinsFromPeakL);
  p->getValue("NBinsFromPeakR",  nBinsFromPeakR);
  p->getValue("MinXFid",         minx_fid);
  p->getValue("MinYFid",         miny_fid);
  p->getValue("MinZFid",         minz_fid);
  p->getValue("MaxXFid",         maxx_fid);
  p->getValue("MaxYFid",         maxy_fid);
  p->getValue("MaxZFid",         maxz_fid);
  p->getValue("MinXAV",          minx_av);
  p->getValue("MinYAV",          miny_av);
  p->getValue("MinZAV",          minz_av);
  p->getValue("MaxXAV",          maxx_av);
  p->getValue("MaxYAV",          maxy_av);
  p->getValue("MaxZAV",          maxz_av);

  // Get the active and fiducial geometry objects
  Geometry fiducial(minx_fid,miny_fid,minz_fid,maxx_fid,maxy_fid,maxz_fid,true);
  Geometry active(minx_av,miny_av,minz_av,maxx_av,maxy_av,maxz_av,false);
  PlaneList extPlanes = active.GetExternalPlaneList();
  PlaneList allPlanes = active.GetPlaneList();
  PlaneList intPlanes = active.GetInternalPlaneList(allPlanes,extPlanes);
  PlaneList fidExtPlanes = fiducial.GetExternalPlaneList();
  PlaneList fidAllPlanes = fiducial.GetPlaneList();

  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Total number of planes in the active volume of the DUNE SP module: " << allPlanes.size() << std::endl;
  std::cout << " Consisting of " << extPlanes.size() << " external planes and " << intPlanes.size() << " internal planes" << std::endl; 
  std::cout << "-----------------------------------------------------------" << std::endl;
 
  // Sort out the file tag
  if(tag != "")
    tag = "_"+tag;

  // Open a text file to print the outputs to
  ofstream ofile;
  ofile.open((location+"/statistics"+tag+".txt").c_str());
  ofile << " --------------------------------------------------------- " << std::endl;
  ofile << std::endl;

  if((nBinsFromPeak+nBinsFromPeakL+nBinsFromPeakR) < 0){
    std::cerr << " Error: Fitting a number of bins from the peak but haven't defined by how many bins" << std::endl;
    std::exit(1);
  }
  else if(nBinsFromPeak > 0 && (nBinsFromPeakL+nBinsFromPeakR) > 0){
    std::cerr << " Error: Fitting a number of bins from the peak but have set both the number (symmetric) and the L&R values (asymmetric)" << std::endl;
    std::exit(1);
  }
  else if(nBinsFromPeak < 0 && (nBinsFromPeakL+nBinsFromPeakR) > 0){
    if(nBinsFromPeakL < 0 || nBinsFromPeakR < 0){
      std::cerr << " Error: Setting the asymmetric binning from the peak but have only set one value (L or R)" << std::endl;
      std::exit(1);
    }
  }
  else{
    nBinsFromPeakL = nBinsFromPeak;
    nBinsFromPeakR = nBinsFromPeak;
  }

  // If we've set L & R set the average to be the general nBins
  if(nBinsFromPeak < 0){
    nBinsFromPeak = (nBinsFromPeakL+nBinsFromPeakR)/2.;
  }
  std::cout << " Number of bins to fit from the peak L : " << nBinsFromPeakL << ", R : " << nBinsFromPeakR << std::endl;
  
  //--------------------------------------------------------------------------------- ---------
  //                                    Initialise
  //--------------------------------------------------------------------------------- ---------

  // Setup TTree from input file list
  std::cout << " Reading files and filling tree" << std::endl;

  EventProcessor evtProc(allowed, input_list, n);
  evtProc.Initialize();

  // Now setup the tree and event objects to work with
  TChain *tree = evtProc.GetTree();
  anatree *evt = evtProc.GetEvents();
  
  // Start of analysis (loop over chain and events
  std::cout << " Running analysis..." << std::endl;

  // Then setup the histograms, counters and any other variables to add to
  // Setup histograms
  TH1D *h_dedx_mctrk  = new TH1D("h_dedx_mctrk", "",200,0.2,6);
  TH1D *h_dedx_geant4 = new TH1D("h_dedx_geant4","",200,0.2,6);

  // Setup counters
  
  // Total number of events
  unsigned int nEvts = tree->GetEntries();
  
  // Now loop over the events
  unsigned int iIt = 1;
  std::cout << " |";
  for(unsigned int iEvt = 0; iEvt < nEvts; ++iEvt){
    tree->GetEntry(iEvt);
    if(!evtProc.SelectEvent(evt)) continue;
    int nGeant = evt->geant_list_size;
    int nMCT   = evt->no_mctracks;
    
    // Print the processing rate
    double evtFrac  = iEvt/static_cast<double>(nEvts);
    if(std::abs(0.1*iIt-evtFrac) < 1e-5){
      std::cout << " --- " << std::setprecision(2) << evtFrac*100 << " %";
      std::cout.flush();
      iIt++;
    }

    // Loop over geant tracks to plot things
    for(int iG4 = 0; iG4 < nGeant; ++iG4){

      TVector3 vtx(evt->StartPointx[iG4],evt->StartPointy[iG4],evt->StartPointz[iG4]);
      TVector3 end(evt->EndPointx[iG4],evt->EndPointy[iG4],evt->EndPointz[iG4]);
      
      // Check the particle enters the TPC volume
      if(!evt->inTPCActive[iG4]) continue;

      TVector3 vtxAV(evt->StartPointx_tpcAV[iG4],evt->StartPointy_tpcAV[iG4],evt->StartPointz_tpcAV[iG4]);
      TVector3 endAV(evt->EndPointx_tpcAV[iG4],evt->EndPointy_tpcAV[iG4],evt->EndPointz_tpcAV[iG4]);

      // If these don't match, the TPC start and end point and general start and end point are not same, 
      bool throughGoing = true;
      
      float dx = abs(endAV.X()-evt->EndPointx[iG4])+abs(vtxAV.X()-evt->StartPointx[iG4]);
      float dy = abs(endAV.Y()-evt->EndPointy[iG4])+abs(vtxAV.Y()-evt->StartPointy[iG4]);
      float dz = abs(endAV.Z()-evt->EndPointz[iG4])+abs(vtxAV.Z()-evt->StartPointz[iG4]);
      
      // If these match, the TPC end point and general end point are the same, therefore the particle stops
      if(dx+dy+dz < 1e-10) throughGoing = false;

      int pdg           = evt->pdg[iG4];
      int id            = evt->TrackId[iG4];
      float lengthAV    = (endAV-vtxAV).Mag();
      float energyDepAV = std::abs(evt->EndE[iG4]-evt->Eng[iG4])/lengthAV; // [GeV/cm]

      // Get the muon
      if(abs(pdg) != 13) continue;
    
      // Check if it's through-going
      if(thru != throughGoing) continue;

      // Make sure it's primary
      if(evt->Mother[iG4] != 0) continue;

      h_dedx_geant4->Fill(energyDepAV);

      // Now do a similar thing with the MCTracks
      for(int iMCT = 0; iMCT < nMCT; ++iMCT){

        // Make sure we're looking at the current G4 track
        int MCTid = evt->mctrk_TrackId[iMCT];
        if(MCTid != id) continue;

        // Add in checks for through-going and energy of the tracks
        // Check the PDG code
        // Check that the track enters the TPC
        // Is there a way to associate mctracks to geant4 objects?
        // Now loop over steps
        int nSteps = evt->nmctrksteps[iMCT];
        for(int iSteps = 0; iSteps < nSteps; ++iSteps){
          h_dedx_mctrk->Fill(evt->mctrk_dEdx[iMCT][iSteps]*100);
        } // iSteps
      } // iMCT
    } // iG4
  }// Event loop
  std::cout << " --- 100 % --- |" << std::endl;

  // Setup things for the analysis 
  TFile *f = new TFile((location+"/calibration_histograms"+tag+".root").c_str(), "RECREATE");
  TCanvas *c = new TCanvas("c","",900,900);
  SetCanvasStyle(c, 0.12,0.06,0.06,0.12,0,0,0);
  f->cd();

  float maxy = std::max(h_dedx_mctrk->GetMaximum(),h_dedx_geant4->GetMaximum());

  SetHistogramStyle1D(h_dedx_mctrk,"True dE/dx [MeV/cm]","Rate");
  
  h_dedx_mctrk->GetYaxis()->SetRangeUser(0,1.1*maxy);
  h_dedx_mctrk->Draw("hist");
  h_dedx_geant4->Draw("hist same");
  h_dedx_mctrk->SetLineWidth(3);
  h_dedx_geant4->SetLineWidth(3);
  h_dedx_mctrk->SetLineColor(pal.at(0));
  h_dedx_geant4->SetLineColor(pal.at(1));
  h_dedx_mctrk->SetLineStyle(2);
  h_dedx_geant4->SetLineStyle(3);

  c->SaveAs((location+"/true_dedx"+tag+".root").c_str());
  c->SaveAs((location+"/true_dedx"+tag+".png").c_str());
  c->Write();
  c->Clear();


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
