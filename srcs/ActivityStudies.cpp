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
  
  // Start of analysis (loop over chain and events
  std::cout << " Running analysis..." << std::endl;

  //TFile *fSake = new TFile((location+"/dummy"+tag+".root").c_str(), "RECREATE");

  // Then setup the histograms, counters and any other variables to add to
  // Setup histograms
  // Truth-level track quantities
  TH1D *h_length     = new TH1D("h_length","",100,0,2.2e3);   // Length of the muons
  TH1D *h_mom        = new TH1D("h_mom","",100,0,2000);       // Momentum of the muons [GeV]
  TH1D *h_energy     = new TH1D("h_energy","",100,0,1000);       // Energy of the muons [GeV]
  TH1D *h_nDaughters = new TH1D("h_nDaughters","",100,0,100); // Number of muon daughters
  TH1D *h_eDepPerL_0 = new TH1D("h_eDepPerL_0","",100,0,5); // Energy deposition per unit length
  TH1D *h_eDepPerL_1 = new TH1D("h_eDepPerL_1","",100,0,5); // Energy deposition per unit length
  TH1D *h_eDepPerL_2 = new TH1D("h_eDepPerL_2","",100,0,5); // Energy deposition per unit length
  TH2D *h_E_nDaught  = new TH2D("h_E_nDaught","",200,0,800,50,0,50); // Number of muon daughters per unit energy

  TH2D *h_eDep_nDaught_0 = new TH2D("h_eDep_nDaught_0","",50,0,50,100,0,5); // Number of muon daughters vs energy depositions per unit length
  TH2D *h_eDep_nDaught_1 = new TH2D("h_eDep_nDaught_1","",50,0,50,100,0,5); // Number of muon daughters vs energy depositions per unit length
  TH2D *h_eDep_nDaught_2 = new TH2D("h_eDep_nDaught_2","",50,0,50,100,0,5); // Number of muon daughters vs energy depositions per unit length
  TH2D *h_eDep_E_0 = new TH2D("h_eDep_E_0","",100,4,1e3,100,0,8); // Number of muon daughters vs energy depositions per unit length
  TH2D *h_eDep_E_1 = new TH2D("h_eDep_E_1","",100,4,1e3,100,0,8); // Number of muon daughters vs energy depositions per unit length
  TH2D *h_eDep_E_2 = new TH2D("h_eDep_E_2","",100,4,1e3,100,0,8); // Number of muon daughters vs energy depositions per unit length
  TH2D *h_log_eDep_E_0 = new TH2D("h_log_eDep_E_0","",100,4,1e3,100,0.01,8); // Number of muon daughters vs energy depositions per unit length
  TH2D *h_log_eDep_E_1 = new TH2D("h_log_eDep_E_1","",100,4,1e3,100,0.01,8); // Number of muon daughters vs energy depositions per unit length
  TH2D *h_log_eDep_E_2 = new TH2D("h_log_eDep_E_2","",100,4,1e3,100,0.01,8); // Number of muon daughters vs energy depositions per unit length
  TH2D *h_reco_eDep_E_0 = new TH2D("h_reco_eDep_E_0","",100,5e-1,10,100,1e-3,1e-2); // Number of muon daughters vs energy depositions per unit length
  TH2D *h_reco_eDep_E_1 = new TH2D("h_reco_eDep_E_1","",100,5e-1,10,100,1e-3,1e-2); // Number of muon daughters vs energy depositions per unit length
  TH2D *h_reco_eDep_E_2 = new TH2D("h_reco_eDep_E_2","",100,5e-1,10,100,1e-3,1e-2); // Number of muon daughters vs energy depositions per unit length

  SetLogX(h_eDep_E_0);
  SetLogX(h_eDep_E_1);
  SetLogX(h_eDep_E_2);
  SetLogX(h_reco_eDep_E_0);
  SetLogX(h_reco_eDep_E_1);
  SetLogX(h_reco_eDep_E_2);
  SetLogY(h_reco_eDep_E_0);
  SetLogY(h_reco_eDep_E_1);
  SetLogY(h_reco_eDep_E_2);
  SetLogX(h_log_eDep_E_0);
  SetLogX(h_log_eDep_E_1);
  SetLogX(h_log_eDep_E_2);
  SetLogY(h_log_eDep_E_0);
  SetLogY(h_log_eDep_E_1);
  SetLogY(h_log_eDep_E_2);

  std::vector<TH1D*> h_eDepPerL{h_eDepPerL_0,h_eDepPerL_1,h_eDepPerL_2};
  std::vector<TH2D*> h_eDep_nDaught{h_eDep_nDaught_0,h_eDep_nDaught_1,h_eDep_nDaught_2};
  std::vector<TH2D*> h_eDep_E{h_eDep_E_0,h_eDep_E_1,h_eDep_E_2};
  std::vector<TH2D*> h_log_eDep_E{h_log_eDep_E_0,h_log_eDep_E_1,h_log_eDep_E_2};
  std::vector<TH2D*> h_reco_eDep_E{h_reco_eDep_E_0,h_reco_eDep_E_1,h_reco_eDep_E_2};
  
  // Setup counters
  
  // Now loop over the events
  unsigned int nEvts = tree->GetEntries();
  unsigned int iIt = 1;

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

    //
    // Truth-level studies
    //
    int iMu = -1;

    // Find the muon
    for(int iG4 = 0; iG4 < nGeant; ++iG4){
      // Check the particle enters the TPC volume
      if(!evt->inTPCActive[iG4]) continue;

      int pdg        = evt->pdg[iG4];
      int id         = evt->TrackId[iG4];

      // Make sure we are looking at a primary muon
      if(evt->Mother[iG4] != 0) continue;
      if(abs(pdg) != 13) continue;
      iMu = id;
      break;
    }
    // Loop over geant tracks to find the delta rays
    int nDelta = 0;
    for(int iG4 = 0; iG4 < nGeant; ++iG4){

      // Check the particle enters the TPC volume
      if(!evt->inTPCActive[iG4]) continue;

      // Make sure we are looking at a primary muon
      if(evt->Mother[iG4] != iMu) continue;

      //int pdg2 = evt->pdg[iG4];
      //if(abs(pdg2) != 11) continue;
      nDelta++;

    } // iGeant
    // Loop over geant tracks to plot things
    for(int iG4 = 0; iG4 < nGeant; ++iG4){

      // Check the particle enters the TPC volume
      if(!evt->inTPCActive[iG4]) continue;

      TVector3 vtxAV(evt->StartPointx_tpcAV[iG4],evt->StartPointy_tpcAV[iG4],evt->StartPointz_tpcAV[iG4]);
      TVector3 endAV(evt->EndPointx_tpcAV[iG4],evt->EndPointy_tpcAV[iG4],evt->EndPointz_tpcAV[iG4]);

      int pdg        = evt->pdg[iG4];
      int id         = evt->TrackId[iG4];
      float lengthAV = (endAV-vtxAV).Mag();

      // Make sure we are looking at a primary muon
      if(evt->Mother[iG4] != 0) continue;
      if(abs(pdg) != 13) continue;

      // Fill the truth-track quantities 
      h_length->Fill(lengthAV);
      h_mom->Fill(evt->P[iG4]);
      h_energy->Fill(evt->Eng[iG4]);

      // For the deposition studies, make sure we are looking at a long track (2m)
      if(lengthAV < 200) continue;

      h_nDaughters->Fill(nDelta);
      h_E_nDaught->Fill(evt->Eng[iG4],nDelta);
      
      // Now loop over hits 
      // First loop over wire planes
      for(int iWire = 0; iWire < 3; ++iWire){
        float totalEDep = 0.;
        float totalQDep = 0.;
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

          // Then get the parameters of interest for this hit
          float hitX      = evt->hit_trueX[iHit];
          float hitE      = evt->hit_energy[iHit];
          float hitQ      = evt->hit_charge[iHit];

          // Check if x is lower than the APA bound, charge seems to accumulate there
          if(hitX < evtProc.APA_X_POSITIONS[0] || hitX > evtProc.APA_X_POSITIONS[2]) continue;

          totalEDep += hitE;
          totalQDep += hitQ;

        }// Hits
        float totalEDepPerLength = totalEDep/lengthAV;
        if(totalEDepPerLength < std::numeric_limits<float>::epsilon()) continue;
        h_eDep_nDaught.at(iWire)->Fill(nDelta,totalEDepPerLength);
        h_eDep_E.at(iWire)->Fill(evt->Eng[iG4],totalEDepPerLength);
        h_log_eDep_E.at(iWire)->Fill(evt->Eng[iG4],totalEDepPerLength);
        h_eDepPerL.at(iWire)->Fill(totalEDepPerLength);
      } // Wire plane
    }// Geant

    //
    // Reco-level studies
    //
    for(int iTrk = 0; iTrk < nTrks; ++iTrk){
      for(int iWire = 0; iWire < 3; ++iWire){
        // Only look at long muons
        if(!(abs(evt->trkpdgtruth_pandoraTrack[iTrk][iWire]) == 13 && evtProc.SelectTrack(evt,iTrk))) continue;

        // Now access the variables of interest
        int nHitsR = evt->ntrkhits_pandoraTrack[iTrk][iWire];
        float energy = evt->trkke_pandoraTrack[iTrk][iWire]/1000.; 

        // Make sure it doesn't exceed the maximum size of the array
        // Count if it does so we can see how often it happens
        if(nHitsR > MAX_TRACK_HITS){
          nHitsR = MAX_TRACK_HITS;
        }

        // Make sure there is some information on the plane
        // Sometimes this is -999 and will break stuff
        if(nHitsR <= 0){
          continue;
        }

        float *dEdxArr = evt->trkdedx_pandoraTrack[iTrk][iWire];
        std::vector<float> dEdx(dEdxArr, dEdxArr + nHitsR);

        // And fill
        for(int iHit = 0; iHit < nHitsR; ++iHit){
          // Check if x is lower or higher than the APA bounds, charge seems to accumulate there
          float x = evt->trkxyz_pandoraTrack[iTrk][iWire][iHit][0];
          if(x < evtProc.APA_X_POSITIONS[0] || x > evtProc.APA_X_POSITIONS[2]) continue;

          float dEdxVal = dEdx.at(iHit);
          float len     = evt->trklen_pandoraTrack[iTrk];
          float eDepLen = dEdxVal/len;
          h_reco_eDep_E.at(iWire)->Fill(energy,eDepLen);
        } // iHit
      } // iWire
    } // iTrk
  }// Event loop
  std::cout << " --- 100 % --- |" << std::endl;

  TCanvas *c1 = new TCanvas("c1","",900,900);
  SetCanvasStyle(c1, 0.12,0.08,0.06,0.12,0,0,0);

  SetHistogramStyle1D(h_length,"Muon length [cm]", "Rate");
  h_length->Draw("hist");
  h_length->SetLineWidth(3);
  h_length->SetLineColor(kTeal-5);
  h_length->GetYaxis()->SetTitleOffset(0.95);
  c1->SaveAs((location+"/truth_tracks_length"+tag+".png").c_str());
  c1->SaveAs((location+"/truth_tracks_length"+tag+".root").c_str());
  c1->Clear();

  SetHistogramStyle1D(h_mom,"Muon momentum [GeV]", "Rate");
  h_mom->Draw("hist");
  h_mom->SetLineWidth(3);
  h_mom->SetLineColor(kTeal-5);
  h_mom->GetYaxis()->SetTitleOffset(0.95);
  c1->SaveAs((location+"/truth_tracks_mom"+tag+".png").c_str());
  c1->SaveAs((location+"/truth_tracks_mom"+tag+".root").c_str());
  c1->Clear();

  SetHistogramStyle1D(h_energy,"Muon energy [GeV]", "Rate");
  h_energy->Draw("hist");
  h_energy->SetLineWidth(3);
  h_energy->SetLineColor(kTeal-5);
  h_energy->GetYaxis()->SetTitleOffset(0.95);
  c1->SaveAs((location+"/truth_tracks_energy"+tag+".png").c_str());
  c1->SaveAs((location+"/truth_tracks_energy"+tag+".root").c_str());
  c1->Clear();

  SetHistogramStyle1D(h_nDaughters,"Muon daughters", "Rate");
  h_nDaughters->Draw("hist");
  h_nDaughters->SetLineWidth(3);
  h_nDaughters->SetLineColor(kTeal-5);
  h_nDaughters->GetYaxis()->SetTitleOffset(0.95);
  c1->SaveAs((location+"/truth_tracks_nDaughters"+tag+".png").c_str());
  c1->SaveAs((location+"/truth_tracks_nDaughters"+tag+".root").c_str());
  c1->Clear();

  for(unsigned int iWire = 0; iWire < 3; ++iWire){
    SetHistogramStyle1D(h_eDepPerL.at(iWire),"Energy deposition per unit length [MeV/cm]", "Rate");
    h_eDepPerL.at(iWire)->Draw("hist");
    h_eDepPerL.at(iWire)->SetLineWidth(3);
    h_eDepPerL.at(iWire)->SetLineColor(kTeal-5);
    h_eDepPerL.at(iWire)->GetYaxis()->SetTitleOffset(0.95);
    c1->SaveAs((location+"/truth_tracks_eDep_perL"+tag+".png").c_str());
    c1->SaveAs((location+"/truth_tracks_eDep_perL"+tag+".root").c_str());
    c1->Clear();
  }
  
  TCanvas *c2 = new TCanvas("c2","",1000,800);
  SetCanvasStyle(c2, 0.1,0.15,0.05,0.12,0,0,0);

  SetHistogramStyle2D(h_E_nDaught,"Muon energy [GeV]", "Number of daughters",false);
  h_E_nDaught->Draw("colz");
  c2->SaveAs((location+"/truth_tracks_nDaught_vs_E"+tag+".png").c_str());
  c2->SaveAs((location+"/truth_tracks_nDaught_vs_E"+tag+".root").c_str());
  c2->Clear();

  for(unsigned int iWire = 0; iWire < 3; ++iWire){
    
    SetHistogramStyle2D(h_eDep_nDaught.at(iWire),"Muon daughters", "Energy deposition per unit length [MeV/cm]",false);
    h_eDep_nDaught.at(iWire)->Draw("colz");
    c2->SaveAs((location+"/truth_tracks_eDep_vs_nDaughters"+std::to_string(iWire)+tag+".png").c_str());
    c2->SaveAs((location+"/truth_tracks_eDep_vs_nDaughters"+std::to_string(iWire)+tag+".root").c_str());
    c2->Clear();

  } // Wire planes

  TCanvas *c3 = new TCanvas("c3","",1000,800);
  SetCanvasStyle(c3, 0.1,0.15,0.05,0.12,0,0,0);
  
  for(unsigned int iWire = 0; iWire < 3; ++iWire){

    c3->SetLogx();
  
    SetHistogramStyle2D(h_eDep_E.at(iWire),"Muon energy [GeV]", "Energy deposition per unit length [MeV/cm]",false);
    h_eDep_E.at(iWire)->Draw("colz");
    c3->SaveAs((location+"/truth_tracks_eDep_vs_E"+std::to_string(iWire)+tag+".png").c_str());
    c3->SaveAs((location+"/truth_tracks_eDep_vs_E"+std::to_string(iWire)+tag+".root").c_str());
    c3->Clear();

  } // Wire planes

  TCanvas *c4 = new TCanvas("c4","",1000,800);
  SetCanvasStyle(c4, 0.1,0.15,0.05,0.12,0,0,0);
  
  for(unsigned int iWire = 0; iWire < 3; ++iWire){

    c4->SetLogx();
    c4->SetLogy();
    
    SetHistogramStyle2D(h_reco_eDep_E.at(iWire),"Muon energy [GeV]", "Energy deposition per unit length [MeV/cm]",false);
    h_reco_eDep_E.at(iWire)->Draw("colz");
    gPad->Update();
    TPaletteAxis *pal = static_cast<TPaletteAxis*>(h_reco_eDep_E.at(iWire)->GetListOfFunctions()->FindObject("palette"));
    pal->SetLabelFont(132);
    pal->SetLabelSize(0.02);
    c4->Modified();
    c4->Update();

    c4->SaveAs((location+"/truth_tracks_reco_eDep_vs_E"+std::to_string(iWire)+tag+".png").c_str());
    c4->SaveAs((location+"/truth_tracks_reco_eDep_vs_E"+std::to_string(iWire)+tag+".root").c_str());
    c4->Clear();

    SetHistogramStyle2D(h_log_eDep_E.at(iWire),"Muon energy [GeV]", "Energy deposition per unit length [MeV/cm]",false);
    h_log_eDep_E.at(iWire)->Draw("colz");
    c4->SaveAs((location+"/truth_tracks_log_eDep_vs_E"+std::to_string(iWire)+tag+".png").c_str());
    c4->SaveAs((location+"/truth_tracks_log_eDep_vs_E"+std::to_string(iWire)+tag+".root").c_str());
    c4->Clear();

  } // Wire planes
  
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
