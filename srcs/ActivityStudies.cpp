/************************************************************************
 * 
 * A macro to understand the truth-level hit and track distributions for 
 * dE/dx calibration studies
 *
 *
 * Example file list located here:
 *   /home/jones/work/cosmics/LArSoft-v08_50_00/work/files/anafiles.list
 *
 *
 *************************************************************************/

#include "EventProcessor.h"
#include "ConfigReader.h"

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
   "ntracks_pandoraTrack",
   "trkId_pandoraTrack",
   "trkidtruth_pandoraTrack",
   "trkpdgtruth_pandoraTrack",
   "trkg4id_pandoraTrack",
   "ntrkhits_pandoraTrack",
   "trkdqdx_pandoraTrack",
   "trkdedx_pandoraTrack",
   "trkke_pandoraTrack",
   "trkxyz_pandoraTrack",
   "trkstartx_pandoraTrack",
   "trkstarty_pandoraTrack",
   "trkstartz_pandoraTrack",
   "trkendx_pandoraTrack",
   "trkendy_pandoraTrack",
   "trkendz_pandoraTrack",
   "trklen_pandoraTrack",
   "no_hits_stored",
   "hit_tpc",
   "hit_plane",
   "hit_charge",
   "hit_energy",
   "hit_nelec",
   "hit_trueX",
   "hit_trkid",
   "hit_peakT",
   "hit_startT",
   "hit_endT",
   "StartPointx",
   "StartPointy",
   "StartPointz",
   "EndPointx",
   "EndPointy",
   "EndPointz",
   "StartPointx_tpcAV",
   "StartPointy_tpcAV",
   "StartPointz_tpcAV",
   "EndPointx_tpcAV",
   "EndPointy_tpcAV",
   "EndPointz_tpcAV",
   "NumberDaughters",
   "P",
   "Eng"
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
     
int activityStudies(const char *config){

  // First, setup timing information so we can monitor the run
  time_t rawtime;
  std::cout << "-----------------------------------------------------------" << std::endl;
  GetTime(rawtime);
  std::cout << "-----------------------------------------------------------" << std::endl;

  //------------------------------------------------------------------------------------------
  //                                    Configure
  //------------------------------------------------------------------------------------------
  // Create object of the class ConfigReader
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
  std::string input_list = "";
  std::string location="";
  std::string tag="";
  std::vector<double> minx_fid, miny_fid, minz_fid;
  std::vector<double> maxx_fid, maxy_fid, maxz_fid;
  std::vector<double> minx_av, miny_av, minz_av;
  std::vector<double> maxx_av, maxy_av, maxz_av;

  p->getValue("InputList", input_list);
  p->getValue("Location",  location);
  p->getValue("Tag",       tag);
  p->getValue("NFiles",    n);
  p->getValue("Thru",      thru);
  p->getValue("MinXFid",   minx_fid);
  p->getValue("MinYFid",   miny_fid);
  p->getValue("MinZFid",   minz_fid);
  p->getValue("MaxXFid",   maxx_fid);
  p->getValue("MaxYFid",   maxy_fid);
  p->getValue("MaxZFid",   maxz_fid);
  p->getValue("MinXAV",    minx_av);
  p->getValue("MinYAV",    miny_av);
  p->getValue("MinZAV",    minz_av);
  p->getValue("MaxXAV",    maxx_av);
  p->getValue("MaxYAV",    maxy_av);
  p->getValue("MaxZAV",    maxz_av);

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
  n = evtProc.GetFiles();
  
  // Start of analysis (loop over chain and events
  std::cout << " Running analysis..." << std::endl;

  //TFile *fSake = new TFile((location+"/dummy"+tag+".root").c_str(), "RECREATE");

  // Then setup the histograms, counters and any other variables to add to
  // Setup histograms
  // Truth-level track quantities
  TH1D *h_length               = new TH1D("h_length","",100,0,2.2e3);   // Length of the muons
  TH1D *h_mom                  = new TH1D("h_mom","",100,0,2000);       // Momentum of the muons [GeV]
  TH1D *h_energy               = new TH1D("h_energy","",100,1e-1,1e5);       // Energy of the muons [GeV]
  TH1D *h_energy_long          = new TH1D("h_energy_long","",100,1e-1,1e5);       // Energy of the muons [GeV]
  TH1D *h_nDaughters           = new TH1D("h_nDaughters","",101,-0.5,100.5); // Number of muon daughters
  TH1D *h_reco_eng             = new TH1D("h_reco_eng","",100,1e-3,15);       // Energy of the muons [GeV]
  TH1D *h_reco_eng_long        = new TH1D("h_reco_eng_long","",100,1e-3,15);       // Energy of the muons [GeV]
  TH1D *h_reco_eng_highy       = new TH1D("h_reco_eng_highy","",100,1e-3,15);       // Energy of the muons [GeV]
  TH1D *h_reco_eng_long_highy  = new TH1D("h_reco_eng_long_highy","",100,1e-3,15);       // Energy of the muons [GeV]
  TH1D *h_reco_len             = new TH1D("h_reco_len","",100,0,1800);       // Length of the muons [cm]
  TH1D *h_eDep_0               = new TH1D("h_eDep_0","",100,0,10000);       // eDep of the muons [GeV]
  TH1D *h_eDep_1               = new TH1D("h_eDep_1","",100,0,10000);       // eDep of the muons [GeV]
  TH1D *h_eDep_2               = new TH1D("h_eDep_2","",100,0,10000);       // eDep of the muons [GeV]
  TH1D *h_eDepPerL_0           = new TH1D("h_eDepPerL_0","",100,0,8); // Energy deposition per unit length
  TH1D *h_eDepPerL_1           = new TH1D("h_eDepPerL_1","",100,0,8); // Energy deposition per unit length
  TH1D *h_eDepPerL_2           = new TH1D("h_eDepPerL_2","",100,0,8); // Energy deposition per unit length
  TH1D *h_eDepPerL_BP          = new TH1D("h_eDepPerL_BP","",100,0,8); // Energy deposition per unit length
  TH1D *h_qDepPerL_0           = new TH1D("h_qDepPerL_0","",100,0,400); // Charge deposition per unit length
  TH1D *h_qDepPerL_1           = new TH1D("h_qDepPerL_1","",100,0,400); // Charge deposition per unit length
  TH1D *h_qDepPerL_2           = new TH1D("h_qDepPerL_2","",100,0,400); // Charge deposition per unit length
  TH1D *h_qDepPerL_BP          = new TH1D("h_qDepPerL_BP","",100,0,400); // Charge deposition per unit length
  TH1D *h_dEdx_E_2_6_10        = new TH1D("h_dEdx_E_2_6_10","",80,0.2,10); // Energy deposition vs energy between Emu = 8 and 12 GeV
  TH1D *h_dEdx_E_2_50_54       = new TH1D("h_dEdx_E_2_50_54","",80,0.2,10); // Energy deposition vs energy between Emu = 50 and 54 GeV
  TH1D *h_dEdx_E_2_278_282     = new TH1D("h_dEdx_E_2_278_282","",60,0.2,10); // Energy deposition vs energy between Emu = 278 and 282 GeV
  TH1D *h_dEdx_E_2_548_552     = new TH1D("h_dEdx_E_2_548_552","",50,0.2,10); // Energy deposition vs energy between Emu = 548 and 552 GeV
  TH1D *h_eDepPerL_E_2_6_10    = new TH1D("h_eDepPerL_E_2_6_10","",80,0,8); // Energy deposition vs energy between Emu = 8 and 12 GeV
  TH1D *h_eDepPerL_E_2_50_54 = new TH1D("h_eDepPerL_E_2_50_54","",80,0,8); // Energy deposition vs energy between Emu = 50 and 54 GeV
  TH1D *h_eDepPerL_E_2_278_282 = new TH1D("h_eDepPerL_E_2_278_282","",60,0,8); // Energy deposition vs energy between Emu = 278 & 282 GeV
  TH1D *h_eDepPerL_E_2_548_552 = new TH1D("h_eDepPerL_E_2_548_552","",50,0,8); // Energy deposition vs energy between Emu = 548 & 552 GeV
  TH1D *h_true_mus             = new TH1D("h_true_mus","",101,-0.5,100.5); // True muon multiplicity
  TH1D *h_true_primary_mus     = new TH1D("h_true_primary_mus","",12,-0.5,11.5); // True muon multiplicity
  TH1D *h_reco_mus             = new TH1D("h_reco_mus","",101,-0.5,100.5); // Reco muon multiplicity
  TH1D *h_reco_long_mus        = new TH1D("h_reco_long_mus","",12,-0.5,11.5); // Reco muon multiplicity
  TH1D *h_reco_long_highy_mus  = new TH1D("h_reco_long_highy_mus","",12,-0.5,11.5); // Reco muon multiplicity
  TH1D *h_nHitsPerL_BP         = new TH1D("h_nHitsPerL_BP","",100,0,2); // Number of hits per unit length

  TH2D *h_E_nDaught      = new TH2D("h_E_nDaught","",200,0,800,51,-0.5,50.5); // Number of muon daughters per unit energy
  TH2D *h_eDep_nDaught_0 = new TH2D("h_eDep_nDaught_0","",51,-0.5,50.5,100,0,5); // Number of muon daughters vs energy depositions per unit length
  TH2D *h_eDep_nDaught_1 = new TH2D("h_eDep_nDaught_1","",51,-0.5,50.5,100,0,5); // Number of muon daughters vs energy depositions per unit length
  TH2D *h_eDep_nDaught_2 = new TH2D("h_eDep_nDaught_2","",51,-0.5,50.5,100,0,5); // Number of muon daughters vs energy depositions per unit length
  TH2D *h_dEdx_nDaught_0 = new TH2D("h_dEdx_nDaught_0","",51,-0.5,50.5,100,0.2,4); // Number of muon daughters vs energy depositions per unit length
  TH2D *h_dEdx_nDaught_1 = new TH2D("h_dEdx_nDaught_1","",51,-0.5,50.5,100,0.2,4); // Number of muon daughters vs energy depositions per unit length
  TH2D *h_dEdx_nDaught_2 = new TH2D("h_dEdx_nDaught_2","",51,-0.5,50.5,100,0.2,5); // Number of muon daughters vs energy depositions per unit length
  TH2D *h_eDep_E_0 = new TH2D("h_eDep_E_0","",100,4,1e3,100,0,8); // Number of muon daughters vs energy depositions per unit length
  TH2D *h_eDep_E_1 = new TH2D("h_eDep_E_1","",100,4,1e3,100,0,8); // Number of muon daughters vs energy depositions per unit length
  TH2D *h_eDep_E_2 = new TH2D("h_eDep_E_2","",100,4,1e3,100,0,8); // Number of muon daughters vs energy depositions per unit length
  TH2D *h_eDep_E_BP = new TH2D("h_eDep_E_BP","",100,4,1e3,100,0,8); // Number of muon daughters vs energy depositions per unit length
  TH2D *h_qDep_E_0 = new TH2D("h_qDep_E_0","",100,4,1e3,100,0,400); // Number of muon daughters vs charge depositions per unit length
  TH2D *h_qDep_E_1 = new TH2D("h_qDep_E_1","",100,4,1e3,100,0,400); // Number of muon daughters vs charge depositions per unit length
  TH2D *h_qDep_E_2 = new TH2D("h_qDep_E_2","",100,4,1e3,100,0,400); // Number of muon daughters vs charge depositions per unit length
  TH2D *h_qDep_E_BP = new TH2D("h_qDep_E_BP","",100,4,1e3,100,0,400); // Number of muon daughters vs charge depositions per unit length
  TH2D *h_reco_eDep_E_0 = new TH2D("h_reco_eDep_E_0","",100,3e-1,10,100,0.01,8); // Number of muon daughters vs energy depositions per unit length
  TH2D *h_reco_eDep_E_1 = new TH2D("h_reco_eDep_E_1","",100,3e-1,10,100,0.01,8); // Number of muon daughters vs energy depositions per unit length
  TH2D *h_reco_eDep_E_2 = new TH2D("h_reco_eDep_E_2","",100,3e-1,10,100,0.01,6); // Number of muon daughters vs energy depositions per unit length
  TH2D *h_dEdx_E_0 = new TH2D("h_dEdx_E_0","",100,4,1e3,100,0.2,10); // Energy deposition vs energy
  TH2D *h_dEdx_E_1 = new TH2D("h_dEdx_E_1","",100,4,1e3,100,0.2,10); // Energy deposition vs energy
  TH2D *h_dEdx_E_2 = new TH2D("h_dEdx_E_2","",100,4,1e3,100,0.2,10); // Energy deposition vs energy
  TH2D *h_reco_dEdx_E_0 = new TH2D("h_reco_dEdx_E_0","",100,5e-1,8,100,0.2,6); // Energy deposition vs energy
  TH2D *h_reco_dEdx_E_1 = new TH2D("h_reco_dEdx_E_1","",100,5e-1,8,100,0.2,6); // Energy deposition vs energy
  TH2D *h_reco_dEdx_E_2 = new TH2D("h_reco_dEdx_E_2","",100,5e-1,8,100,0.2,6); // Energy deposition vs energy
  TH2D *h_reco_Y_E      = new TH2D("h_reco_Y_E","",100,5e-1,8,100,-600,600); // Energy vs reconstructed Y position 
  TH2D *h_reco_Y_E_zoom = new TH2D("h_reco_Y_E_zoom","",100,5e-1,8,100,596,600); // Energy vs reconstructed Y position 
  TH2D *h_reco_len_E    = new TH2D("h_reco_len_E","",100,5e-1,8,100,2,1800); // Energy vs reconstructed len position 

  SetLogX(h_eDep_E_0);
  SetLogX(h_eDep_E_1);
  SetLogX(h_eDep_E_2);
  SetLogX(h_eDep_E_BP);
  SetLogX(h_qDep_E_0);
  SetLogX(h_qDep_E_1);
  SetLogX(h_qDep_E_2);
  SetLogX(h_qDep_E_BP);
  SetLogX(h_reco_eDep_E_0);
  SetLogX(h_reco_eDep_E_1);
  SetLogX(h_reco_eDep_E_2);
  SetLogX(h_reco_dEdx_E_0);
  SetLogX(h_reco_dEdx_E_1);
  SetLogX(h_reco_dEdx_E_2);
  SetLogX(h_reco_eng);
  SetLogX(h_reco_eng_long);
  SetLogX(h_reco_eng_highy);
  SetLogX(h_reco_eng_long_highy);
  SetLogX(h_energy);
  SetLogX(h_energy_long);
  SetLogX(h_dEdx_E_0);
  SetLogX(h_dEdx_E_1);
  SetLogX(h_dEdx_E_2);
  SetLogX(h_reco_Y_E);
  SetLogX(h_reco_Y_E_zoom);
  SetLogX(h_reco_len_E);

  std::vector<TH1D*> h_eDep{h_eDep_0,h_eDep_1,h_eDep_2};
  std::vector<TH1D*> h_eDepPerL{h_eDepPerL_0,h_eDepPerL_1,h_eDepPerL_2};
  std::vector<TH1D*> h_qDepPerL{h_qDepPerL_0,h_qDepPerL_1,h_qDepPerL_2};
  std::vector<TH2D*> h_eDep_nDaught{h_eDep_nDaught_0,h_eDep_nDaught_1,h_eDep_nDaught_2};
  std::vector<TH2D*> h_dEdx_nDaught{h_dEdx_nDaught_0,h_dEdx_nDaught_1,h_dEdx_nDaught_2};
  std::vector<TH2D*> h_eDep_E{h_eDep_E_0,h_eDep_E_1,h_eDep_E_2};
  std::vector<TH2D*> h_qDep_E{h_qDep_E_0,h_qDep_E_1,h_qDep_E_2};
  std::vector<TH2D*> h_reco_eDep_E{h_reco_eDep_E_0,h_reco_eDep_E_1,h_reco_eDep_E_2};
  std::vector<TH2D*> h_dEdx_E{h_dEdx_E_0,h_dEdx_E_1,h_dEdx_E_2};
  std::vector<TH2D*> h_reco_dEdx_E{h_reco_dEdx_E_0,h_reco_dEdx_E_1,h_reco_dEdx_E_2};
  
  // Setup counters
  
  // Now loop over the events
  unsigned int nEvts = tree->GetEntries();
  unsigned int iIt = 1;

  double total_energy_true = 0.;
  double total_energy_reco = 0.;
  int n_mu_true = 0;
  int n_mu_reco = 0;
  int nHugeE = 0;

  std::cout << " |";
  for(unsigned int iEvt = 0; iEvt < nEvts; ++iEvt){
    tree->GetEntry(iEvt);
    if(!evtProc.SelectEvent(evt)) continue;
    
    // Print the processing rate
    double evtFrac  = iEvt/static_cast<double>(nEvts);

    if(std::abs(0.1*iIt-evtFrac) < std::numeric_limits<double>::epsilon()){
      std::cout << " --- " << evtFrac*100 << " %";
      std::cout.flush();
      iIt++;
    }

    // Geant and hit iterator definitions
    int nGeant = evt->geant_list_size;
    int nHits  = evt->no_hits_stored;
    int nTrks  = evt->ntracks_pandoraTrack;

    int n_true_mus            = 0;
    int n_true_primary_mus    = 0;
    int n_reco_mus            = 0;
    int n_reco_long_mus       = 0;
    int n_reco_long_highy_mus = 0;
    //
    // Truth-level studies
    //
    // Loop over geant tracks to plot things
    for(int iG4 = 0; iG4 < nGeant; ++iG4){

      // Check the particle enters the TPC volume
      if(!evt->inTPCActive[iG4]) continue;

      TVector3 vtxAV(evt->StartPointx_tpcAV[iG4],evt->StartPointy_tpcAV[iG4],evt->StartPointz_tpcAV[iG4]);
      TVector3 endAV(evt->EndPointx_tpcAV[iG4],evt->EndPointy_tpcAV[iG4],evt->EndPointz_tpcAV[iG4]);

      // Determine if the track enters at the top and leaves through the bottom
      float vtxDy = abs(vtxAV.Y()-evt->StartPointy[iG4]);
      float endDy = abs(endAV.Y()-evt->EndPointy[iG4]);

      // If these don't match, the TPC end point and general end point are not same, 
      // therefore the particle goes through the top and bottom faces of the detector
      bool throughGoing = false;
      if(vtxDy+endDy > 1e-10) throughGoing = true;

      int pdg        = evt->pdg[iG4];
      int id         = evt->TrackId[iG4];
      float lengthAV = (endAV-vtxAV).Mag();

      // Make sure we are looking at a primary muon
      if(abs(pdg) != 13) continue;
      n_true_mus++;
      
      if(evt->Mother[iG4] != 0) continue;
      n_true_primary_mus++;

      // Fill the truth-track quantities 
      h_length->Fill(lengthAV);
      h_mom->Fill(evt->P[iG4]);
      h_energy->Fill(evt->Eng[iG4]);
      total_energy_true += evt->Eng[iG4];

      // For the deposition studies, make sure we are looking at a long track (2m)
      if(lengthAV < 300) continue;
      n_mu_true++;
      h_energy_long->Fill(evt->Eng[iG4]);

      h_nDaughters->Fill(evt->NumberDaughters[iG4]);
      h_E_nDaught->Fill(evt->Eng[iG4],evt->NumberDaughters[iG4]);
      
      // Get the best plane
      int bestPlane = 0;
      std::vector<int> hitsOnPlane(3,0);
      for(int iPlane = 0; iPlane < 3; ++iPlane){
        for(int iHit = 0; iHit < nHits; ++iHit){
          // If we are not looking at the current G4 track, continue
          bool currentG4 = false;
          for(int iTrk = 0; iTrk < evt->ntracks_pandoraTrack; ++iTrk){
            int recoId = evt->trkId_pandoraTrack[iTrk];
            int hitId  = evt->hit_trkid[iHit];

            // Check if the current hit is in the reco track
            if(recoId != hitId) continue;

            // If it is, check if the reco track is the G4 track
            if(evt->trkidtruth_pandoraTrack[iTrk][iPlane] == id){
              currentG4 = true;
              break;
            }
          }
          if(!currentG4) continue;
          if(evt->hit_plane[iHit] == iPlane)
            hitsOnPlane.at(iPlane)++;
        } // Hits
      } // Planes
      bestPlane = std::max_element(hitsOnPlane.begin(), hitsOnPlane.end()) - hitsOnPlane.begin();
      h_nHitsPerL_BP->Fill(hitsOnPlane.at(bestPlane)/static_cast<double>(lengthAV));
//      std::cout << " Best plane: " << bestPlane << " with " << hitsOnPlane.at(bestPlane) << " hits and length " << lengthAV << std::endl;

      // Now loop over hits 
      // First loop over wire planes
      for(int iWire = 0; iWire < 3; ++iWire){
        float totalEDep = 0.;
        float totalQDep = 0.;
        // Skip the first and last hits for 'reconstructability' purposes
        for(int iHit = 0; iHit < nHits; ++iHit){
          // Skip the current hit if it wasn't deposited on this plane
          if(evt->hit_plane[iHit] != iWire) continue;

          // If we are not looking at the current G4 track, continue
          bool currentG4 = false;
          for(int iTrk = 0; iTrk < evt->ntracks_pandoraTrack; ++iTrk){
            int recoId = evt->trkId_pandoraTrack[iTrk];
            int hitId  = evt->hit_trkid[iHit];

            // Check if the current hit is in the reco track
            if(recoId != hitId) continue;

            // If it is, check if the reco track is the G4 track
            if(evt->trkidtruth_pandoraTrack[iTrk][iWire] == id){
              currentG4 = true;
              break;
            }
          }
          if(!currentG4) continue;
          if(hitsOnPlane.at(iWire)/static_cast<double>(lengthAV) < 0.2) continue; // If the number of hits on this plane is silly w.r.t the length
          if(thru && !throughGoing) continue; // Check if we should be looking at through-going tracks

          // Then get the parameters of interest for this hit
          float hitX      = evt->hit_trueX[iHit];
          float hitE      = evt->hit_energy[iHit];
          float hitQ      = evt->hit_charge[iHit];

          // Check if x is lower than the APA bound, charge seems to accumulate there
          if(hitX < evtProc.APA_X_POSITIONS[0] || hitX > evtProc.APA_X_POSITIONS[2]) continue;

          int tpc =evtProc.WhichTPC(hitX) + 1;
          float dx = ( -1 + 2*(tpc%2) )*(hitX - evtProc.APA_X_POSITIONS[tpc/2]);
          
          totalEDep += hitE;
          totalQDep += hitQ;

          h_dEdx_E.at(iWire)->Fill(evt->Eng[iG4],hitE);
          h_dEdx_nDaught.at(iWire)->Fill(evt->NumberDaughters[iG4],hitE);

          // Now fill the single energy bins with dEdx on the collection plane
          if(iWire == 2){
            if(evt->Eng[iG4] >= 6 && evt->Eng[iG4] <= 10){
              h_dEdx_E_2_6_10->Fill(hitE);
            } // 6-10 GeV
            else if(evt->Eng[iG4] >= 50 && evt->Eng[iG4] <= 54){
              h_dEdx_E_2_50_54->Fill(hitE);
            } // 50-54 GeV
            else if(evt->Eng[iG4] >= 278 && evt->Eng[iG4] <= 282){
              h_dEdx_E_2_278_282->Fill(hitE);
            } // 278-282 GeV
            else if(evt->Eng[iG4] >= 548 && evt->Eng[iG4] <= 552){
              h_dEdx_E_2_548_552->Fill(hitE);
            } // 548-552 GeV
          } // Wire plane
        }// Hits
        float totalEDepPerLength = totalEDep/static_cast<double>(lengthAV);
        float totalQDepPerLength = totalQDep/static_cast<double>(lengthAV);
        if(totalEDepPerLength < std::numeric_limits<float>::epsilon()) continue;
        if(totalQDepPerLength < std::numeric_limits<float>::epsilon()) continue;
        h_eDep_nDaught.at(iWire)->Fill(evt->NumberDaughters[iG4],totalEDepPerLength);
        h_qDep_E.at(iWire)->Fill(evt->Eng[iG4],totalQDepPerLength);
        h_eDep_E.at(iWire)->Fill(evt->Eng[iG4],totalEDepPerLength);
        h_eDepPerL.at(iWire)->Fill(totalEDepPerLength);
        h_qDepPerL.at(iWire)->Fill(totalQDepPerLength);
        h_eDep.at(iWire)->Fill(totalEDep);
        if(iWire == bestPlane){
          h_qDep_E_BP->Fill(evt->Eng[iG4],totalQDepPerLength);
          h_eDep_E_BP->Fill(evt->Eng[iG4],totalEDepPerLength);
          h_qDepPerL_BP->Fill(totalQDepPerLength);
          h_eDepPerL_BP->Fill(totalEDepPerLength);
        }
        if(iWire == 2){
          if(evt->Eng[iG4] >= 6 && evt->Eng[iG4] <= 10){
            h_eDepPerL_E_2_6_10->Fill(totalEDepPerLength);
          } // 6-10 GeV
          else if(evt->Eng[iG4] >= 50 && evt->Eng[iG4] <= 54){
            h_eDepPerL_E_2_50_54->Fill(totalEDepPerLength);
          } // 50-54 GeV
          else if(evt->Eng[iG4] >= 278 && evt->Eng[iG4] <= 282){
            h_eDepPerL_E_2_278_282->Fill(totalEDepPerLength);
          } // 278-282 GeV
          else if(evt->Eng[iG4] >= 548 && evt->Eng[iG4] <= 552){
            h_eDepPerL_E_2_548_552->Fill(totalEDepPerLength);
          } // 548-552 GeV
        } // Collection plane
      } // Wire plane
    }// Geant
    h_true_mus->Fill(n_true_mus);
    h_true_primary_mus->Fill(n_true_primary_mus);

    //
    // Reco-level studies
    //
    for(int iTrk = 0; iTrk < nTrks; ++iTrk){
      
      // Get the track geometry
      TVector3 startVtx(evt->trkstartx_pandoraTrack[iTrk],
                        evt->trkstarty_pandoraTrack[iTrk],
                        evt->trkstartz_pandoraTrack[iTrk]);
      TVector3 endVtx(evt->trkendx_pandoraTrack[iTrk],
                      evt->trkendy_pandoraTrack[iTrk],
                      evt->trkendz_pandoraTrack[iTrk]);
      
      // Since we are generating downwards-going tracks - if the start y < end y then 
      // assume the reconstruction has got them the wrong way around and flip them
      if(startVtx.Y() < endVtx.Y()){
        TVector3 temp(endVtx);
        endVtx = startVtx;
        startVtx = temp;
      }

      // Get the best plane
      int bestPlane = 0;
      int currHits  = -999;
      for(int iPlane = 0; iPlane < 3; ++iPlane){
        if(evt->ntrkhits_pandoraTrack[iTrk][iPlane] > currHits){
          currHits  = evt->ntrkhits_pandoraTrack[iTrk][iPlane];
          bestPlane = iPlane; 
        } // CurrHits
      } // Planes

      float len = evt->trklen_pandoraTrack[iTrk];
      h_reco_len->Fill(len);
      
      for(int iWire = 0; iWire < 3; ++iWire){
        // Only look at long muons
        if(abs(evt->trkpdgtruth_pandoraTrack[iTrk][iWire]) != 13) continue;

        // Make sure there is some information on the plane
        // Sometimes this is -999 and will break stuff
        int nHitsR = evt->ntrkhits_pandoraTrack[iTrk][iWire];
        if(nHitsR <= 0){
          continue;
        }

        float energy = evt->trkke_pandoraTrack[iTrk][iWire]/1000.; // GeV
        
        // Make sure the energy isn't stupid
        // Make sure it doesn't exceed the maximum energy of a true track: 1e5 GeV
        if(energy > 1e5) {
          nHugeE++;
          continue;
        }

        if(iWire == bestPlane){
          // Reconstructed energy vs start y position
          h_reco_Y_E->Fill(energy,startVtx.Y());
          h_reco_Y_E_zoom->Fill(energy,startVtx.Y());
          h_reco_len_E->Fill(energy,len);
          h_reco_eng->Fill(energy);
          if(startVtx.Y() > 599.5)
            h_reco_eng_highy->Fill(energy);
        
          // If the current reconstructed track is associated to a true muon
          n_reco_mus++;
        }
        
        if(!evtProc.SelectTrack(evt,iTrk)) continue;

        if(iWire == bestPlane){
          // Now access the variables of interest
          h_reco_eng_long->Fill(energy);
          total_energy_reco += energy;
          if(startVtx.Y() > 599.5){
            n_mu_reco++;
            h_reco_eng_long_highy->Fill(energy);
            n_reco_long_highy_mus++;
          }
          
          // For long muon multiplicity plot
          n_reco_long_mus++;
        }

        // Make sure it doesn't exceed the maximum size of the array
        // Count if it does so we can see how often it happens
        if(nHitsR > MAX_TRACK_HITS){
          nHitsR = MAX_TRACK_HITS;
        }

        float *dEdxArr = evt->trkdedx_pandoraTrack[iTrk][iWire];
        std::vector<float> dEdx(dEdxArr, dEdxArr + nHitsR);

        float totalDep = 0.;

        // And fill
        for(int iHit = 0; iHit < nHitsR; ++iHit){
          // Check if x is lower or higher than the APA bounds, charge seems to accumulate there
          float x = evt->trkxyz_pandoraTrack[iTrk][iWire][iHit][0];

          if(x < evtProc.APA_X_POSITIONS[0] || x > evtProc.APA_X_POSITIONS[2]) continue;
          
          int tpc =evtProc.WhichTPC(x) + 1;
          float dx = ( -1 + 2*(tpc%2) )*(x - evtProc.APA_X_POSITIONS[tpc/2]);

          float dEdxVal = dEdx.at(iHit);
          h_reco_dEdx_E.at(iWire)->Fill(energy,dEdxVal);

          totalDep += dEdxVal;//*dx;
        } // iHit

        float eDepLen = totalDep/(len); // MeV/cm
        h_reco_eDep_E.at(iWire)->Fill(energy,eDepLen);

      } // iWire
    } // iTrk
    h_reco_mus->Fill(n_reco_mus);
    h_reco_long_mus->Fill(n_reco_long_mus);
    h_reco_long_highy_mus->Fill(n_reco_long_highy_mus);
  }// Event loop
  std::cout << " --- 100 % --- |" << std::endl;

  TCanvas *c0 = new TCanvas("c0","",900,900);
  SetCanvasStyle(c0, 0.12,0.08,0.06,0.12,0,0,0);
  c0->SetLogy();

  TLegend *l = new TLegend(0.22,0.94,0.92,0.995);
  l->SetNColumns(2);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->SetTextFont(132);

  SetHistogramStyle1D(h_true_mus,"Muon multiplicity", "Rate");
  SetHistogramStyle1D(h_reco_mus,"Muon multiplicity", "Rate");

  // Now sort out the range
  double max_y = 1.1*std::max(h_true_mus->GetMaximum(),h_reco_mus->GetMaximum());

  // Now draw
  h_true_mus->Draw("hist");
  h_reco_mus->Draw("hist same");
  h_true_mus->SetLineWidth(3);
  h_true_mus->SetLineStyle(2);
  h_true_mus->SetLineColor(kTeal-5);
  h_reco_mus->SetLineStyle(2);
  h_reco_mus->SetLineWidth(3);
  h_reco_mus->SetLineColor(kViolet-5);
  h_true_mus->GetYaxis()->SetRangeUser(0.5,max_y);

  l->AddEntry(h_true_mus,"True muons", "l");
  l->AddEntry(h_reco_mus,"Reco muons", "l");
  l->Draw("same");

  c0->SaveAs((location+"/muon_multiplicities"+tag+".png").c_str());
  c0->SaveAs((location+"/muon_multiplicities"+tag+".root").c_str());
  c0->Clear();
  l->Clear();

  SetHistogramStyle1D(h_true_primary_mus,"Muon multiplicity", "Rate");
  SetHistogramStyle1D(h_reco_long_mus,"Muon multiplicity", "Rate");

  // Now sort out the range
  max_y = 1.1*std::max(h_true_primary_mus->GetMaximum(),h_reco_long_mus->GetMaximum());

  // Now draw
  h_true_primary_mus->Draw("hist");
  h_reco_long_mus->Draw("hist same");
  h_true_primary_mus->SetLineWidth(3);
  h_true_primary_mus->SetLineStyle(2);
  h_true_primary_mus->SetLineColor(kTeal-5);
  h_reco_long_mus->SetLineStyle(2);
  h_reco_long_mus->SetLineWidth(3);
  h_reco_long_mus->SetLineColor(kViolet-5);
  h_true_primary_mus->GetYaxis()->SetRangeUser(0.5,max_y);

  l->AddEntry(h_true_primary_mus,"True muons", "l");
  l->AddEntry(h_reco_long_mus,"Reco muons", "l");
  l->Draw("same");

  c0->SaveAs((location+"/long_muon_multiplicities"+tag+".png").c_str());
  c0->SaveAs((location+"/long_muon_multiplicities"+tag+".root").c_str());
  c0->Clear();
  l->Clear();

  SetHistogramStyle1D(h_true_primary_mus,"Muon multiplicity", "Rate");
  SetHistogramStyle1D(h_reco_long_highy_mus,"Muon multiplicity", "Rate");

  // Now sort out the range
  max_y = 1.1*std::max(h_true_primary_mus->GetMaximum(),h_reco_long_highy_mus->GetMaximum());

  // Now draw
  h_true_primary_mus->Draw("hist");
  h_reco_long_highy_mus->Draw("hist same");
  h_true_primary_mus->SetLineWidth(3);
  h_true_primary_mus->SetLineStyle(2);
  h_true_primary_mus->SetLineColor(kTeal-5);
  h_reco_long_highy_mus->SetLineStyle(2);
  h_reco_long_highy_mus->SetLineWidth(3);
  h_reco_long_highy_mus->SetLineColor(kViolet-5);
  h_true_primary_mus->GetYaxis()->SetRangeUser(0.5,max_y);

  l->AddEntry(h_true_primary_mus,"True muons", "l");
  l->AddEntry(h_reco_long_highy_mus,"Reco muons", "l");
  l->Draw("same");

  c0->SaveAs((location+"/long_highy_muon_multiplicities"+tag+".png").c_str());
  c0->SaveAs((location+"/long_highy_muon_multiplicities"+tag+".root").c_str());
  c0->Clear();
  l->Clear();

  c0->SetLogx();
  
  SetHistogramStyle1D(h_energy,"Muon energy [GeV]", "Rate/GeV");
  h_energy->Scale(1,"width");
  h_energy->Draw("hist");
  h_energy->SetLineWidth(3);
  h_energy->SetLineColor(kTeal-5);
  c0->SaveAs((location+"/energy"+tag+".png").c_str());
  c0->SaveAs((location+"/energy"+tag+".root").c_str());
  c0->Clear();

  SetHistogramStyle1D(h_energy_long,"Muon energy [GeV]", "Rate/GeV");
  h_energy_long->Scale(1,"width");
  h_energy_long->Draw("hist");
  h_energy_long->SetLineWidth(3);
  h_energy_long->SetLineColor(kTeal-5);
  c0->SaveAs((location+"/energy_long"+tag+".png").c_str());
  c0->SaveAs((location+"/energy_long"+tag+".root").c_str());
  c0->Clear();

  TCanvas *c1 = new TCanvas("c1","",900,900);
  SetCanvasStyle(c1, 0.12,0.08,0.06,0.12,0,0,0);

  SetHistogramStyle1D(h_length,"Muon length [cm]", "Rate");
  h_length->Draw("hist");
  h_length->SetLineWidth(3);
  h_length->SetLineColor(kTeal-5);
  c1->SaveAs((location+"/length"+tag+".png").c_str());
  c1->SaveAs((location+"/length"+tag+".root").c_str());
  c1->Clear();

  SetHistogramStyle1D(h_mom,"Muon momentum [GeV]", "Rate");
  h_mom->Draw("hist");
  h_mom->SetLineWidth(3);
  h_mom->SetLineColor(kTeal-5);
  c1->SaveAs((location+"/mom"+tag+".png").c_str());
  c1->SaveAs((location+"/mom"+tag+".root").c_str());
  c1->Clear();

  SetHistogramStyle1D(h_nDaughters,"Muon daughters", "Rate");
  h_nDaughters->Draw("hist");
  h_nDaughters->SetLineWidth(3);
  h_nDaughters->SetLineColor(kTeal-5);
  c1->SaveAs((location+"/nDaughters"+tag+".png").c_str());
  c1->SaveAs((location+"/nDaughters"+tag+".root").c_str());
  c1->Clear();
  
  SetHistogramStyle1D(h_reco_len,"Muon length [cm]", "Rate");
  h_reco_len->Draw("hist");
  h_reco_len->SetLineWidth(3);
  h_reco_len->SetLineColor(kTeal-5);
  c1->SaveAs((location+"/reco_length"+tag+".png").c_str());
  c1->SaveAs((location+"/reco_length"+tag+".root").c_str());
  c1->Clear();
 
  l->SetNColumns(1);
  l->SetX1NDC(0.465);
  l->SetY1NDC(0.688);
  l->SetX2NDC(0.894);
  l->SetY2NDC(0.923);
 
  SetHistogramStyle1D(h_dEdx_E_2_6_10,"dE/dx [MeV/cm]", "Normalised rate");
  SetHistogramStyle1D(h_dEdx_E_2_50_54,"dE/dx [MeV/cm]", "Normalised rate");
  SetHistogramStyle1D(h_dEdx_E_2_278_282,"dE/dx [MeV/cm]", "Normalised rate");
  SetHistogramStyle1D(h_dEdx_E_2_548_552,"dE/dx [MeV/cm]", "Normalised rate");

  // Scale the histograms to compare shapes only
  h_dEdx_E_2_6_10->Scale(1/static_cast<double>(h_dEdx_E_2_6_10->Integral()),"width");
  h_dEdx_E_2_50_54->Scale(1/static_cast<double>(h_dEdx_E_2_50_54->Integral()),"width");
  h_dEdx_E_2_278_282->Scale(1/static_cast<double>(h_dEdx_E_2_278_282->Integral()),"width");
  h_dEdx_E_2_548_552->Scale(1/static_cast<double>(h_dEdx_E_2_548_552->Integral()),"width");

  // Now sort out the range
  max_y = 1.1*std::max(h_dEdx_E_2_6_10->GetMaximum(),std::max(h_dEdx_E_2_50_54->GetMaximum(),std::max(h_dEdx_E_2_278_282->GetMaximum(),h_dEdx_E_2_548_552->GetMaximum())));

  // Now draw
  h_dEdx_E_2_6_10->Draw("hist");
  h_dEdx_E_2_50_54->Draw("hist same");
  h_dEdx_E_2_278_282->Draw("hist same");
 // h_dEdx_E_2_548_552->Draw("hist same");
  h_dEdx_E_2_6_10->SetLineWidth(3);
  h_dEdx_E_2_6_10->SetLineColor(kTeal-5);
  h_dEdx_E_2_50_54->SetLineWidth(3);
  h_dEdx_E_2_50_54->SetLineColor(kViolet-5);
  h_dEdx_E_2_278_282->SetLineWidth(3);
  h_dEdx_E_2_278_282->SetLineColor(kOrange+5);
  h_dEdx_E_2_548_552->SetLineWidth(3);
  h_dEdx_E_2_548_552->SetLineColor(kPink-5);
  h_dEdx_E_2_6_10->GetYaxis()->SetRangeUser(0,max_y);

  l->AddEntry(h_dEdx_E_2_6_10,"6 #leq E_{#mu} #leq 10 GeV", "l");
  l->AddEntry(h_dEdx_E_2_50_54,"50 #leq E_{#mu} #leq 54 GeV", "l");
  l->AddEntry(h_dEdx_E_2_278_282,"278 #leq E_{#mu} #leq 282 GeV", "l");
 // l->AddEntry(h_dEdx_E_2_548_552,"548 #leq E_{#mu} #leq 552 GeV", "l");
  l->Draw("same");

  c1->SaveAs((location+"/dEdx_slices"+tag+".png").c_str());
  c1->SaveAs((location+"/dEdx_slices"+tag+".root").c_str());
  c1->Clear();
  l->Clear();

  SetHistogramStyle1D(h_eDepPerL_E_2_6_10,"Total energy deposition per unit length [MeV/cm]", "Normalised rate");
  SetHistogramStyle1D(h_eDepPerL_E_2_50_54,"Total energy deposition per unit length [MeV/cm]", "Normalised rate");
  SetHistogramStyle1D(h_eDepPerL_E_2_278_282,"Total energy deposition per unit length [MeV/cm]", "Normalised rate");
  SetHistogramStyle1D(h_eDepPerL_E_2_548_552,"Total energy deposition per unit length [MeV/cm]", "Normalised rate");

  // Scale the histograms to compare shapes only
  h_eDepPerL_E_2_6_10->Scale(1/static_cast<double>(h_eDepPerL_E_2_6_10->Integral()),"width");
  h_eDepPerL_E_2_50_54->Scale(1/static_cast<double>(h_eDepPerL_E_2_50_54->Integral()),"width");
  h_eDepPerL_E_2_278_282->Scale(1/static_cast<double>(h_eDepPerL_E_2_278_282->Integral()),"width");
  h_eDepPerL_E_2_548_552->Scale(1/static_cast<double>(h_eDepPerL_E_2_548_552->Integral()),"width");

  // Now sort out the range
  max_y = 1.1*std::max(h_eDepPerL_E_2_6_10->GetMaximum(),std::max(h_eDepPerL_E_2_50_54->GetMaximum(),std::max(h_eDepPerL_E_2_278_282->GetMaximum(),h_eDepPerL_E_2_548_552->GetMaximum())));

  // Now draw
  h_eDepPerL_E_2_6_10->Draw("hist");
  h_eDepPerL_E_2_50_54->Draw("hist same");
  h_eDepPerL_E_2_278_282->Draw("hist same");
 // h_eDepPerL_E_2_548_552->Draw("hist same");
  h_eDepPerL_E_2_6_10->SetLineWidth(3);
  h_eDepPerL_E_2_6_10->SetLineColor(kTeal-5);
  h_eDepPerL_E_2_50_54->SetLineWidth(3);
  h_eDepPerL_E_2_50_54->SetLineColor(kViolet-5);
  h_eDepPerL_E_2_278_282->SetLineWidth(3);
  h_eDepPerL_E_2_278_282->SetLineColor(kOrange+5);
  h_eDepPerL_E_2_548_552->SetLineWidth(3);
  h_eDepPerL_E_2_548_552->SetLineColor(kPink-5);
  h_eDepPerL_E_2_6_10->GetYaxis()->SetRangeUser(0,max_y);

  l->AddEntry(h_eDepPerL_E_2_6_10,"6 #leq E_{#mu} #leq 10 GeV", "l");
  l->AddEntry(h_eDepPerL_E_2_50_54,"50 #leq E_{#mu} #leq 54 GeV", "l");
  l->AddEntry(h_eDepPerL_E_2_278_282,"278 #leq E_{#mu} #leq 282 GeV", "l");
 // l->AddEntry(h_eDepPerL_E_2_548_552,"548 #leq E_{#mu} #leq 552 GeV", "l");
  l->Draw("same");

  c1->SaveAs((location+"/eDepPerL_slices"+tag+".png").c_str());
  c1->SaveAs((location+"/eDepPerL_slices"+tag+".root").c_str());
  c1->Clear();
  l->Clear();

  h_dEdx_E_2_6_10->Draw("hist");
  h_dEdx_E_2_6_10->SetLineWidth(3);
  h_dEdx_E_2_6_10->SetLineColor(kTeal-5);

  c1->SaveAs((location+"/dEdx_slices_6-10_GeV"+tag+".png").c_str());
  c1->SaveAs((location+"/dEdx_slices_6-10_GeV"+tag+".root").c_str());
  c1->Clear();

  h_dEdx_E_2_50_54->Draw("hist");
  h_dEdx_E_2_50_54->SetLineWidth(3);
  h_dEdx_E_2_50_54->SetLineColor(kViolet-5);

  c1->SaveAs((location+"/dEdx_slices_50-54_GeV"+tag+".png").c_str());
  c1->SaveAs((location+"/dEdx_slices_50-54_GeV"+tag+".root").c_str());
  c1->Clear();

  h_dEdx_E_2_278_282->Draw("hist");
  h_dEdx_E_2_278_282->SetLineWidth(3);
  h_dEdx_E_2_278_282->SetLineColor(kOrange+5);

  c1->SaveAs((location+"/dEdx_slices_278-282_GeV"+tag+".png").c_str());
  c1->SaveAs((location+"/dEdx_slices_278-282_GeV"+tag+".root").c_str());
  c1->Clear();

  h_dEdx_E_2_548_552->Draw("hist");
  h_dEdx_E_2_548_552->SetLineWidth(3);
  h_dEdx_E_2_548_552->SetLineColor(kPink-5);

  c1->SaveAs((location+"/dEdx_slices_548-552_GeV"+tag+".png").c_str());
  c1->SaveAs((location+"/dEdx_slices_548-552_GeV"+tag+".root").c_str());
  c1->Clear();

  h_eDepPerL_E_2_6_10->Draw("hist");
  h_eDepPerL_E_2_6_10->SetLineWidth(3);
  h_eDepPerL_E_2_6_10->SetLineColor(kTeal-5);

  c1->SaveAs((location+"/eDepPerL_slices_6-10_GeV"+tag+".png").c_str());
  c1->SaveAs((location+"/eDepPerL_slices_6-10_GeV"+tag+".root").c_str());
  c1->Clear();

  h_eDepPerL_E_2_50_54->Draw("hist");
  h_eDepPerL_E_2_50_54->SetLineWidth(3);
  h_eDepPerL_E_2_50_54->SetLineColor(kViolet-5);

  c1->SaveAs((location+"/eDepPerL_slices_50-54_GeV"+tag+".png").c_str());
  c1->SaveAs((location+"/eDepPerL_slices_50-54_GeV"+tag+".root").c_str());
  c1->Clear();

  h_eDepPerL_E_2_278_282->Draw("hist");
  h_eDepPerL_E_2_278_282->SetLineWidth(3);
  h_eDepPerL_E_2_278_282->SetLineColor(kOrange+5);

  c1->SaveAs((location+"/eDepPerL_slices_278-282_GeV"+tag+".png").c_str());
  c1->SaveAs((location+"/eDepPerL_slices_278-282_GeV"+tag+".root").c_str());
  c1->Clear();

  h_eDepPerL_E_2_548_552->Draw("hist");
  h_eDepPerL_E_2_548_552->SetLineWidth(3);
  h_eDepPerL_E_2_548_552->SetLineColor(kPink-5);

  c1->SaveAs((location+"/eDepPerL_slices_548-552_GeV"+tag+".png").c_str());
  c1->SaveAs((location+"/eDepPerL_slices_548-552_GeV"+tag+".root").c_str());
  c1->Clear();

  for(unsigned int iWire = 0; iWire < 3; ++iWire){
    SetHistogramStyle1D(h_eDepPerL.at(iWire),"Energy deposition per unit length [MeV/cm]", "Rate");
    h_eDepPerL.at(iWire)->Draw("hist");
    h_eDepPerL.at(iWire)->SetLineWidth(3);
    h_eDepPerL.at(iWire)->SetLineColor(kTeal-5);
    c1->SaveAs((location+"/eDep_perL"+std::to_string(iWire)+tag+".png").c_str());
    c1->SaveAs((location+"/eDep_perL"+std::to_string(iWire)+tag+".root").c_str());
    c1->Clear();
    
    SetHistogramStyle1D(h_qDepPerL.at(iWire),"Charge deposition per unit length [MeV/cm]", "Rate");
    h_qDepPerL.at(iWire)->Draw("hist");
    h_qDepPerL.at(iWire)->SetLineWidth(3);
    h_qDepPerL.at(iWire)->SetLineColor(kTeal-5);
    c1->SaveAs((location+"/qDep_perL"+std::to_string(iWire)+tag+".png").c_str());
    c1->SaveAs((location+"/qDep_perL"+std::to_string(iWire)+tag+".root").c_str());
    c1->Clear();
    
   SetHistogramStyle1D(h_eDep.at(iWire),"Total energy deposited [MeV]", "Rate");
    h_eDep.at(iWire)->Draw("hist");
    h_eDep.at(iWire)->SetLineWidth(3);
    h_eDep.at(iWire)->SetLineColor(kTeal-5);
    c1->SaveAs((location+"/eDep"+std::to_string(iWire)+tag+".png").c_str());
    c1->SaveAs((location+"/eDep"+std::to_string(iWire)+tag+".root").c_str());
    c1->Clear();
  }
  SetHistogramStyle1D(h_eDepPerL_BP,"Energy deposition per unit length [MeV/cm]", "Rate");
  h_eDepPerL_BP->Draw("hist");
  h_eDepPerL_BP->SetLineWidth(3);
  h_eDepPerL_BP->SetLineColor(kTeal-5);
  c1->SaveAs((location+"/eDep_perL_BP"+tag+".png").c_str());
  c1->SaveAs((location+"/eDep_perL_BP"+tag+".root").c_str());
  c1->Clear();

  SetHistogramStyle1D(h_qDepPerL_BP,"Charge deposition per unit length [MeV/cm]", "Rate");
  h_qDepPerL_BP->Draw("hist");
  h_qDepPerL_BP->SetLineWidth(3);
  h_qDepPerL_BP->SetLineColor(kTeal-5);
  c1->SaveAs((location+"/qDep_perL_BP"+tag+".png").c_str());
  c1->SaveAs((location+"/qDep_perL_BP"+tag+".root").c_str());
  c1->Clear();
  
  SetHistogramStyle1D(h_nHitsPerL_BP,"Hits per unit length [cm^{-1}]", "Rate");
  h_nHitsPerL_BP->Draw("hist");
  h_nHitsPerL_BP->SetLineWidth(3);
  h_nHitsPerL_BP->SetLineColor(kTeal-5);
  c1->SaveAs((location+"/hits_per_L"+tag+".png").c_str());
  c1->SaveAs((location+"/hits_per_L"+tag+".root").c_str());
  c1->Clear();

  c1->SetLogx();
  c1->SetLogy();
  
  SetHistogramStyle1D(h_reco_eng,"Total muon deposition [GeV]", "Rate");
  h_reco_eng->Scale(1,"width");
  h_reco_eng->Draw("hist");
  h_reco_eng->SetLineWidth(3);
  h_reco_eng->SetLineColor(kTeal-5);
  c1->SaveAs((location+"/reco_energy"+tag+".png").c_str());
  c1->SaveAs((location+"/reco_energy"+tag+".root").c_str());
  c1->Clear();
  
  SetHistogramStyle1D(h_reco_eng_long,"Total muon deposition [GeV]", "Rate");
  h_reco_eng_long->Scale(1,"width");
  h_reco_eng_long->Draw("hist");
  h_reco_eng_long->SetLineWidth(3);
  h_reco_eng_long->SetLineColor(kTeal-5);
  c1->SaveAs((location+"/reco_energy_long"+tag+".png").c_str());
  c1->SaveAs((location+"/reco_energy_long"+tag+".root").c_str());
  c1->Clear();
  
  SetHistogramStyle1D(h_reco_eng_highy,"Total muon deposition [GeV]", "Rate");
  h_reco_eng_highy->Scale(1,"width");
  h_reco_eng_highy->Draw("hist");
  h_reco_eng_highy->SetLineWidth(3);
  h_reco_eng_highy->SetLineColor(kTeal-5);
  c1->SaveAs((location+"/reco_energy_highy"+tag+".png").c_str());
  c1->SaveAs((location+"/reco_energy_highy"+tag+".root").c_str());
  c1->Clear();
  
  SetHistogramStyle1D(h_reco_eng_long_highy,"Total muon deposition [GeV]", "Rate");
  h_reco_eng_long_highy->Scale(1,"width");
  h_reco_eng_long_highy->Draw("hist");
  h_reco_eng_long_highy->SetLineWidth(3);
  h_reco_eng_long_highy->SetLineColor(kTeal-5);
  c1->SaveAs((location+"/reco_energy_long_highy"+tag+".png").c_str());
  c1->SaveAs((location+"/reco_energy_long_highy"+tag+".root").c_str());
  c1->Clear();
  
  TCanvas *c2 = new TCanvas("c2","",1000,800);
  SetCanvasStyle(c2, 0.1,0.15,0.05,0.12,0,0,0);

  SetHistogramStyle2D(h_E_nDaught,"Muon energy [GeV]", "Number of daughters",false);
  h_E_nDaught->Draw("colz");
  c2->SaveAs((location+"/nDaught_vs_E"+tag+".png").c_str());
  c2->SaveAs((location+"/nDaught_vs_E"+tag+".root").c_str());
  c2->Clear();

  for(unsigned int iWire = 0; iWire < 3; ++iWire){
    
    SetHistogramStyle2D(h_eDep_nDaught.at(iWire),"Muon daughters", "Energy deposition per unit length [MeV/cm]",false);
    h_eDep_nDaught.at(iWire)->Draw("colz");
    c2->SaveAs((location+"/eDep_vs_nDaughters"+std::to_string(iWire)+tag+".png").c_str());
    c2->SaveAs((location+"/eDep_vs_nDaughters"+std::to_string(iWire)+tag+".root").c_str());
    c2->Clear();

    SetHistogramStyle2D(h_dEdx_nDaught.at(iWire),"Muon daughters", "dE/dx [MeV/cm]",false);
    h_dEdx_nDaught.at(iWire)->Draw("colz");
    c2->SaveAs((location+"/dEdx_vs_nDaughters"+std::to_string(iWire)+tag+".png").c_str());
    c2->SaveAs((location+"/dEdx_vs_nDaughters"+std::to_string(iWire)+tag+".root").c_str());
    c2->Clear();

  } // Wire planes

  TCanvas *c3 = new TCanvas("c3","",1000,800);
  SetCanvasStyle(c3, 0.1,0.15,0.05,0.12,0,0,0);
  c3->SetLogx();
  
  SetHistogramStyle2D(h_reco_Y_E,"Total muon deposition [GeV]", "Muon start Y Position [cm]",false);
  h_reco_Y_E->Scale(1,"width");
  h_reco_Y_E->Draw("colz");
  c3->SaveAs((location+"/reco_Y_vs_E"+tag+".png").c_str());
  c3->SaveAs((location+"/reco_Y_vs_E"+tag+".root").c_str());
  c3->Clear();

  SetHistogramStyle2D(h_reco_Y_E_zoom,"Total muon deposition [GeV]", "Muon start Y Position [cm]",false);
  h_reco_Y_E_zoom->Scale(1,"width");
  h_reco_Y_E_zoom->Draw("colz");
  c3->SaveAs((location+"/reco_Y_vs_E_zoom"+tag+".png").c_str());
  c3->SaveAs((location+"/reco_Y_vs_E_zoom"+tag+".root").c_str());
  c3->Clear();

  SetHistogramStyle2D(h_reco_len_E,"Total muon deposition [GeV]", "Muon track length [cm]",false);
  h_reco_len_E->Scale(1,"width");
  h_reco_len_E->Draw("colz");
  c3->SaveAs((location+"/reco_len_vs_E"+tag+".png").c_str());
  c3->SaveAs((location+"/reco_len_vs_E"+tag+".root").c_str());
  c3->Clear();

  for(unsigned int iWire = 0; iWire < 3; ++iWire){
  
    SetHistogramStyle2D(h_eDep_E.at(iWire),"Muon energy [GeV]", "Energy deposition per unit length [MeV/cm]",false);
    h_eDep_E.at(iWire)->Scale(1,"width");
    h_eDep_E.at(iWire)->Draw("colz");
    c3->SaveAs((location+"/eDep_vs_E"+std::to_string(iWire)+tag+".png").c_str());
    c3->SaveAs((location+"/eDep_vs_E"+std::to_string(iWire)+tag+".root").c_str());
    c3->Clear();

    SetHistogramStyle2D(h_qDep_E.at(iWire),"Muon energy [GeV]", "Charge deposition per unit length [MeV/cm]",false);
    h_qDep_E.at(iWire)->Scale(1,"width");
    h_qDep_E.at(iWire)->Draw("colz");
    c3->SaveAs((location+"/qDep_vs_E"+std::to_string(iWire)+tag+".png").c_str());
    c3->SaveAs((location+"/qDep_vs_E"+std::to_string(iWire)+tag+".root").c_str());
    c3->Clear();

    SetHistogramStyle2D(h_dEdx_E.at(iWire),"Muon energy [GeV]", "dE/dx [MeV/cm]",false);
    h_dEdx_E.at(iWire)->Scale(1,"width");
    h_dEdx_E.at(iWire)->Draw("colz");
    c3->SaveAs((location+"/dEdx_vs_E"+std::to_string(iWire)+tag+".png").c_str());
    c3->SaveAs((location+"/dEdx_vs_E"+std::to_string(iWire)+tag+".root").c_str());
    c3->Clear();
    
    SetHistogramStyle2D(h_reco_eDep_E.at(iWire),"Total muon deposition [GeV]", "Energy deposition per unit length [MeV/cm]",false);
    h_reco_eDep_E.at(iWire)->Scale(1,"width");
    h_reco_eDep_E.at(iWire)->Draw("colz");
    c3->SaveAs((location+"/reco_eDep_vs_E"+std::to_string(iWire)+tag+".png").c_str());
    c3->SaveAs((location+"/reco_eDep_vs_E"+std::to_string(iWire)+tag+".root").c_str());
    c3->Clear();

    SetHistogramStyle2D(h_reco_dEdx_E.at(iWire),"Total muon deposition [GeV]", "dE/dx [MeV/cm]",false);
    h_reco_dEdx_E.at(iWire)->Scale(1,"width");
    h_reco_dEdx_E.at(iWire)->Draw("colz");
    c3->SaveAs((location+"/reco_dEdx_vs_E"+std::to_string(iWire)+tag+".png").c_str());
    c3->SaveAs((location+"/reco_dEdx_vs_E"+std::to_string(iWire)+tag+".root").c_str());
    c3->Clear();

  } // Wire planes
  SetHistogramStyle2D(h_eDep_E_BP,"Muon energy [GeV]", "Energy deposition per unit length [MeV/cm]",false);
  h_eDep_E_BP->Scale(1,"width");
  h_eDep_E_BP->Draw("colz");
  c3->SaveAs((location+"/eDep_vs_E_BP"+tag+".png").c_str());
  c3->SaveAs((location+"/eDep_vs_E_BP"+tag+".root").c_str());
  c3->Clear();

  SetHistogramStyle2D(h_qDep_E_BP,"Muon energy [GeV]", "Charge deposition per unit length [MeV/cm]",false);
  h_qDep_E_BP->Scale(1,"width");
  h_qDep_E_BP->Draw("colz");
  c3->SaveAs((location+"/qDep_vs_E_BP"+tag+".png").c_str());
  c3->SaveAs((location+"/qDep_vs_E_BP"+tag+".root").c_str());
  c3->Clear();
  
  // Now calculate averages
  double average_energy_true = total_energy_true / static_cast<double>(n_mu_true);
  double average_energy_reco = total_energy_reco / static_cast<double>(n_mu_reco);

  // Now find the peak
  double peak_energy_true = GetPeakBinCentre(h_energy_long);
  double peak_energy_reco = GetPeakBinCentre(h_reco_eng_long);

  ofstream oFile;
  oFile.open((location+"/statistics"+tag+".txt").c_str());

  oFile << " Total muons true      : " << n_mu_true << ", reco E < 1e5 GeV : " << n_mu_reco << ", reco E > 1e5 GeV : " << nHugeE << std::endl;
  oFile << " Total energy true     : " << total_energy_true << ", reco: " << total_energy_reco << " GeV " << std::endl;
  oFile << " Average energy true   : " << average_energy_true << ", reco: " << average_energy_reco << " GeV " << std::endl;
  oFile << " Peak energy true      : " << peak_energy_true << ", reco: " << peak_energy_reco << " GeV " << std::endl;

  oFile.close();

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
