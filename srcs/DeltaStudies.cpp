/************************************************************************
 *
 * A Macro to plot various distrubtions in order to 
 * understand the contents of the CR sample for 
 * through-going muons
 *
 * Example file list located here:
 *   /home/jones/work/cosmics/LArSoft-v08_50_00/work/files/anafiles.list
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
   "genie_primaries_pdg",
   "inTPCActive",
   "pdg",
   "TrackId",
   "Mother",
   "Eng",
   "process_primary",
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
   "hit_plane",
   "hit_trkid",
   "hit_nelec",
   "nclusters",
   "cluster_NHits",
   "clusterId",
   "process_primary",
   "trkpdgtruth_pandoraTrack",
   "trkg4id_pandoraTrack",
   "trkId_pandoraTrack",
   "trkidtruth_pandoraTrack",
   "ntracks_pandoraTrack",
   "ntrkhits_pandoraTrack",
   "trkxyz_pandoraTrack",
   "trkstartx_pandoraTrack",
   "trkstarty_pandoraTrack",
   "trkstartz_pandoraTrack",
   "trkendx_pandoraTrack",
   "trkendy_pandoraTrack",
   "trkendz_pandoraTrack",
   "trklen_pandoraTrack",
   "nPFParticles",
   "pfp_selfID",
   "pfp_isPrimary",
   "pfp_isTrack",
   "pfp_numClusters",
   "pfp_clusterIDs",
   "pfp_numDaughters",
   "pfp_daughterIDs"
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
     
int deltaStudies(const char *config){

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
  int yCut = 1;
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
  p->getValue("YCut",      yCut);
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
  std::cout << " Number of files: " << n << std::endl;
  
  // Start of analysis (loop over chain and events
  std::cout << " Running analysis..." << std::endl;

  // Then setup the histograms, counters and any other variables to add to
  // Setup histograms
  TH1D *h_muon_length        = new TH1D("h_muon_length","",200,0,2000); // Muon length
  TH1D *h_delta_length       = new TH1D("h_delta_length","",200,0,40); // Secondary length
  TH1D *h_pfp_numDaughters   = new TH1D("h_pfp_numDaughters","",100,0,100); // Number of primary pfparticle daughters
  TH1D *h_pfp_numDauHits     = new TH1D("h_pfp_numDauHits","",200,0,3000); // Number of primary pfparticle daughters
  TH1D *h_hits_muL           = new TH1D("h_hits_muL","",100,0,0.4); // Number of hits per unit length of the primary muon track
  TH1D *h_secondaries_muL    = new TH1D("h_secondaries_muL","",100,0,0.05); // Number of secondaries per unit length of the primary muon track
  
  TH2D *h_muon_length_E        = new TH2D("h_muon_length_E","",100,5,5e3,100,0,2000); // Muon length
  TH2D *h_hits_muL_E           = new TH2D("h_hits_muL_E","",100,5,5e3,100,0,0.12);
  TH2D *h_hits_E               = new TH2D("h_hits_E","",100,5,5e3,80,0,240);
  TH2D *h_secondaries_muL_E    = new TH2D("h_secondaries_muL_E","",100,5,5e3,100,0,0.03);
  TH2D *h_secondaries_E        = new TH2D("h_secondaries_E","",100,5,5e3,50,0,50);
  TH2D *h_pfp_numDaughters_E   = new TH2D("h_pfp_numDaughters_E","",100,5,5e3,80,0,80); // Number of primary pfparticle daughters vs energy
  TH2D *h_pfp_numDau_muL_E     = new TH2D("h_pfp_numDau_muL_E","",100,5,5e3,80,0,0.1); // Number of primary pfparticle daughters vs energy
  TH2D *h_pfp_numDauHits_E     = new TH2D("h_pfp_numDauHits_E","",100,5,5e3,100,0,3000); // Number of primary pfparticle daughters vs energy
  TH2D *h_pfp_numDauHits_muL_E = new TH2D("h_pfp_numDauHits_muL_E","",100,5,5e3,100,0,5); // Number of primary pfparticle daughters vs energy
 
  // Sort out log scales if needed 
  SetLogX(h_muon_length_E);
  SetLogX(h_hits_muL_E);
  SetLogX(h_hits_E);
  SetLogX(h_secondaries_muL_E);
  SetLogX(h_secondaries_E);
  SetLogX(h_pfp_numDaughters_E);
  SetLogX(h_pfp_numDauHits_E);
  SetLogX(h_pfp_numDau_muL_E);
  SetLogX(h_pfp_numDauHits_muL_E);
  
  // Setup counters
  
  // Now loop over the events
  unsigned int nEvts = tree->GetEntries();
  unsigned int iIt = 1;

  std::cout << " |";
  for(unsigned int iEvt = 0; iEvt < nEvts; ++iEvt){
    if(!tree->GetEntry(iEvt)) {
      std::cout << " Entry " << iEvt << " is somehow corrupt, skipping this event" << std::endl;
      continue;
    }
    tree->GetEntry(iEvt);
    if(!evtProc.SelectEvent(evt)) continue;
    unsigned int nTrks = evt->ntracks_pandoraTrack;
    unsigned int nPfps = evt->nPFParticles;
    int nGeant = evt->geant_list_size;
    
    // Print the processing rate
    double evtFrac  = iEvt/static_cast<double>(nEvts);

    if((std::abs(0.1*iIt)-evtFrac) < std::numeric_limits<double>::epsilon()){
      std::cout << " --- " << evtFrac*100 << " %";
      std::cout.flush();
      iIt++;
    }

    // First, find the primary muon and assign the relevant quantities
    //    Length
    //    True energy
    //    Number of hits/length
    float primaryLength     = -99999.;
    float primaryTrueEnergy = -99999.;
    float primaryHitsPerL   = -99999.;
    for(unsigned int iTrk = 0; iTrk < nTrks; ++iTrk){
      // Get the best plane
      int bestPlane = 0;
      int currHits  = -999;
      for(int iPlane = 0; iPlane < 3; ++iPlane){
        if(evt->ntrkhits_pandoraTrack[iTrk][iPlane] > currHits){
          currHits  = evt->ntrkhits_pandoraTrack[iTrk][iPlane];
          bestPlane = iPlane; 
        } // CurrHits
      } // Planes

      int ID = evt->trkId_pandoraTrack[iTrk];

      // Only look at primary muons
      if(abs(evt->trkpdgtruth_pandoraTrack[iTrk][bestPlane]) != 13) continue;

      // Make sure it is associated with the true primary particle
      
      // Length cuts (2m)
      if(!evtProc.SelectTrack(evt,iTrk)) continue;
      
      // Get the track geometry
      TVector3 startVtx(evt->trkstartx_pandoraTrack[iTrk],
                        evt->trkstarty_pandoraTrack[iTrk],
                        evt->trkstartz_pandoraTrack[iTrk]);
      TVector3 endVtx(evt->trkendx_pandoraTrack[iTrk],
                      evt->trkendy_pandoraTrack[iTrk],
                      evt->trkendz_pandoraTrack[iTrk]);

      // Since we are generating downwards-going tracks - if the start y < end y then 
      // assume the reconstruction has got them the wrong way around and flip them
      //
      // Temporarily remove these tracks for energy-deposition's sake - only the longest muon candidate
      if(startVtx.Y() < endVtx.Y()){
        TVector3 temp(endVtx);
        endVtx = startVtx;
        startVtx = temp;
        continue;
      }
      
      // The following studies should be conducted with through-going muons to start with
      // If the number of external planes crossed is >= 2, the track is through-going
      primaryLength = evt->trklen_pandoraTrack[iTrk];
      bool throughGoing = IsThroughGoing(primaryLength,startVtx,endVtx,extPlanes,fidExtPlanes);
      if(thru && !throughGoing) continue;

      h_muon_length->Fill(primaryLength);
        
      for(int iPlane = 0; iPlane < 3; ++iPlane){

        // Use only best plane for now
        if(iPlane != bestPlane) continue;

        // Get the associated true energy of the muon
        // Try to get the true energy
        // Get the list iterator from matching ID's
        // Also make sure this muon is primary
        bool primaryMu = false;
        for(int iG4 = 0; iG4 < nGeant; ++iG4){
          int trueID = evt->TrackId[iG4];

          if(evt->trkidtruth_pandoraTrack[iTrk][bestPlane] == trueID){
            primaryTrueEnergy = evt->Eng[iG4];
            if(evt->process_primary[iG4] == 1) primaryMu = true;
          }
        }
        if(!primaryMu) continue;
        if(primaryTrueEnergy < 0){
          std::cout << " Warning: Energy is below zero, skipping track with energy: " << primaryTrueEnergy << std::endl;
          continue;
        }
        unsigned int nHits = evt->ntrkhits_pandoraTrack[iTrk][iPlane];
        primaryHitsPerL = nHits/primaryLength;

      } // Plane
      h_muon_length_E->Fill(primaryTrueEnergy,primaryLength);
    } // iTrk

    // Make sure we found a primary muon
    if(primaryLength < 0 && primaryTrueEnergy < 0) continue;

    // Now loop over the tracks so we can look at delta rays
    // Setup the counter for the total number of hits per muon
    unsigned int totDeltaHits = 0;
    unsigned int totDelta     = 0;
    for(unsigned int iTrk = 0; iTrk < nTrks; ++iTrk){

      // Get the best plane
      int bestPlane = 0;
      int currHits  = -999;
      for(int iPlane = 0; iPlane < 3; ++iPlane){
        if(evt->ntrkhits_pandoraTrack[iTrk][iPlane] > currHits){
          currHits  = evt->ntrkhits_pandoraTrack[iTrk][iPlane];
          bestPlane = iPlane; 
        } // CurrHits
      } // Planes

      int ID = evt->trkId_pandoraTrack[iTrk];

      // Only look at primary muons
      if(abs(evt->trkpdgtruth_pandoraTrack[iTrk][bestPlane]) != 13) continue;
      
      // Length cuts (3m)
      // If we are looking at a long track associated to the true muon, it is likely to be the muon itself
      // The shorter associated tracks are likely to be the delta rays - skip the muon itself
      if(evtProc.SelectTrack(evt,iTrk)) continue;
      
      // Make sure this is associated to a primary muon
      bool primaryMu = false;
      for(int iG4 = 0; iG4 < nGeant; ++iG4){
        int trueID = evt->TrackId[iG4];

        if(evt->trkidtruth_pandoraTrack[iTrk][bestPlane] == trueID){
          if(evt->process_primary[iG4] == 1) primaryMu = true;
        }
      }
      if(!primaryMu) continue;
      if(evt->ntrkhits_pandoraTrack[iTrk][bestPlane] == 0) continue;
      totDelta++;
      totDeltaHits += evt->ntrkhits_pandoraTrack[iTrk][bestPlane];
      float deltaLength = evt->trklen_pandoraTrack[iTrk];
      h_delta_length->Fill(deltaLength);
    } // Tracks
    float hitsPerMuonL  = totDeltaHits/primaryLength;
    float deltaPerMuonL = totDelta/primaryLength;
    if(totDeltaHits == 0) continue;
    h_hits_muL->Fill(hitsPerMuonL);
    h_secondaries_muL->Fill(deltaPerMuonL);
    h_hits_muL_E->Fill(primaryTrueEnergy,hitsPerMuonL);
    h_secondaries_muL_E->Fill(primaryTrueEnergy,deltaPerMuonL);
    h_hits_E->Fill(primaryTrueEnergy,totDeltaHits);
    h_secondaries_E->Fill(primaryTrueEnergy,totDelta);

    // Now loop over the PFParticles in the event, count how many are 'primary' and then plot daughter information
    // Use cluster associations to access primary muon from most hits, 
    // then access number of daughters and get the clusters/hits from there
    unsigned int nPrimaryPFPs = 0;
    int maxHitsPFP = -999;
    int primaryPFPID = -999;
    for(unsigned int iPfp = 0; iPfp < nPfps; ++iPfp){
      bool isPrimaryPFP = evt->pfp_isPrimary[iPfp];
      if(isPrimaryPFP){
        nPrimaryPFPs++;

        // Now loop over the clusters and corresponding hits to determine the 'dominant' primary
        int nPFPHits = 0;
        unsigned int nClu = evt->pfp_numClusters[iPfp];
        if(nClu > MAX_PFPCLUSTERS){
          nClu = MAX_PFPCLUSTERS;
        }
        for(unsigned int iClu = 0; iClu < nClu; ++iClu){
          // Get the corresponding cluster ID and find it in the cluster list
          int cluIDPFP = evt->pfp_clusterIDs[iPfp][iClu];

          // Now the clusters on their own
          unsigned int nclu = evt->nclusters;
          for(unsigned int iclu = 0; iclu < nclu; ++iclu){
            int cluID = evt->clusterId[iclu];
            if(cluID == cluIDPFP){
              nPFPHits += evt->cluster_NHits[iclu];     
            }
          } // clusters
        } // PFP clusters
        if(nPFPHits > maxHitsPFP){
          maxHitsPFP = nPFPHits;
          primaryPFPID = evt->pfp_selfID[iPfp];
        }
      } // Primary
    } // PFParticles

    // Make sure we found a primary
    if(primaryPFPID < 0) continue;

    // Now we have identified the primary PFP, get its daughters and corresponding hits
    for(unsigned int iPfp = 0; iPfp < nPfps; ++iPfp){
      int PFPID = evt->pfp_selfID[iPfp];
      if(PFPID != primaryPFPID) continue;

      unsigned int nDau = evt->pfp_numDaughters[iPfp];
      if(nDau == 0) continue;
      unsigned int nDauCut = nDau;
      if(nDau > MAX_PFPDAUGHTERS){
        nDauCut = MAX_PFPDAUGHTERS;
      }
      // Now loop over daughters and get the number of hits associated with it (HACK)
      unsigned int nDauHits = 0;
      for(unsigned int iDau = 0; iDau < nDauCut; ++iDau){
        int dauID = evt->pfp_daughterIDs[iPfp][iDau];

        // Now loop back over the PFParticles and get the clusters/hits from this daughter ID
        for(unsigned int jPfp = 0; jPfp < nPfps; ++jPfp){
          if(evt->pfp_selfID[jPfp] != dauID) continue;
          
          // Now loop over the clusters and corresponding hits to determine the 'dominant' primary
          unsigned int nClu = evt->pfp_numClusters[jPfp];
          if(nClu > MAX_PFPCLUSTERS){
            nClu = MAX_PFPCLUSTERS;
          }
          for(unsigned int iClu = 0; iClu < nClu; ++iClu){
            // Get the corresponding cluster ID and find it in the cluster list
            int cluIDPFP = evt->pfp_clusterIDs[jPfp][iClu];

            // Now the clusters on their own
            unsigned int nclu = evt->nclusters;
            for(unsigned int iclu = 0; iclu < nclu; ++iclu){
              int cluID = evt->clusterId[iclu];
              if(cluID == cluIDPFP){
                nDauHits += evt->cluster_NHits[iclu];     
              }
            } // clusters
          } // PFP clusters
        } // jPfp
      } // PFParticles
      if(nDauHits == 0) continue;
      float daughtersPerL = nDau/primaryLength;
      float daughterHitsPerL = nDauHits/primaryLength;
      h_pfp_numDaughters->Fill(nDau);
      h_pfp_numDauHits->Fill(nDauHits);
      h_pfp_numDaughters_E->Fill(primaryTrueEnergy,nDau);
      h_pfp_numDauHits_E->Fill(primaryTrueEnergy,nDauHits);
      h_pfp_numDau_muL_E->Fill(primaryTrueEnergy,daughtersPerL);
      h_pfp_numDauHits_muL_E->Fill(primaryTrueEnergy,daughterHitsPerL);
    } // PFParticles
  }// Event loop
  std::cout << " --- 100 % --- |" << std::endl;

  TCanvas *c0 = new TCanvas("c0","",900,900);
  SetCanvasStyle(c0, 0.12,0.05,0.05,0.12,0,0,0);

  SetHistogramStyle1D(h_muon_length,"Muon length [cm]", " Rate");
  h_muon_length->Draw("hist");
  h_muon_length->SetLineWidth(2);
  h_muon_length->SetLineColor(kViolet-5);
  h_muon_length->GetYaxis()->SetTitleOffset(0.95);
  
  c0->SaveAs((location+"/muon_length"+tag+".png").c_str());
  c0->SaveAs((location+"/muon_length"+tag+".root").c_str());
  c0->Clear();
  
  SetHistogramStyle1D(h_delta_length,"Delta length [cm]", " Rate");
  h_delta_length->Draw("hist");
  h_delta_length->SetLineWidth(2);
  h_delta_length->SetLineColor(kViolet-5);
  h_delta_length->GetYaxis()->SetTitleOffset(0.95);
  
  c0->SaveAs((location+"/delta_length"+tag+".png").c_str());
  c0->SaveAs((location+"/delta_length"+tag+".root").c_str());
  c0->Clear();
  
  SetHistogramStyle1D(h_hits_muL,"N_{Hits}/L_{#mu} [cm^{-1}]", " Rate");
  h_hits_muL->Draw("hist");
  h_hits_muL->SetLineWidth(2);
  h_hits_muL->SetLineColor(kViolet-5);
  h_hits_muL->GetYaxis()->SetTitleOffset(0.95);
  
  c0->SaveAs((location+"/hits_per_Lmu"+tag+".png").c_str());
  c0->SaveAs((location+"/hits_per_Lmu"+tag+".root").c_str());
  c0->Clear();
  
  SetHistogramStyle1D(h_secondaries_muL,"N_{#delta}/L_{#mu} [cm^{-1}]", " Rate");
  h_secondaries_muL->Draw("hist");
  h_secondaries_muL->SetLineWidth(2);
  h_secondaries_muL->SetLineColor(kViolet-5);
  h_secondaries_muL->GetYaxis()->SetTitleOffset(0.95);
  
  c0->SaveAs((location+"/secondaries_per_Lmu"+tag+".png").c_str());
  c0->SaveAs((location+"/secondaries_per_Lmu"+tag+".root").c_str());
  c0->Clear();
  
  SetHistogramStyle1D(h_pfp_numDaughters,"N_{Daughters}", " Rate");
  h_pfp_numDaughters->Draw("hist");
  h_pfp_numDaughters->SetLineWidth(2);
  h_pfp_numDaughters->SetLineColor(kViolet-5);
  h_pfp_numDaughters->GetYaxis()->SetTitleOffset(0.95);
  
  c0->SaveAs((location+"/pfp_numDaughters"+tag+".png").c_str());
  c0->SaveAs((location+"/pfp_numDaughters"+tag+".root").c_str());
  c0->Clear();
  
  SetHistogramStyle1D(h_pfp_numDauHits,"N_{DaughterHits}", " Rate");
  h_pfp_numDauHits->Draw("hist");
  h_pfp_numDauHits->SetLineWidth(2);
  h_pfp_numDauHits->SetLineColor(kViolet-5);
  h_pfp_numDauHits->GetYaxis()->SetTitleOffset(0.95);
  
  c0->SaveAs((location+"/pfp_numDaughterHits"+tag+".png").c_str());
  c0->SaveAs((location+"/pfp_numDaughterHits"+tag+".root").c_str());
  c0->Clear();
  
  TCanvas *c1 = new TCanvas("c1","",1000,800);
  SetCanvasStyle(c1, 0.12,0.12,0.05,0.12,0,0,0);

  c1->SetLogx();
  SetHistogramStyle2D(h_muon_length_E, "True (generated) #mu energy [GeV]", "L_{#mu} [cm]", false);
  h_muon_length_E->GetZaxis()->SetLabelSize(0.03);
  h_muon_length_E->GetZaxis()->SetLabelFont(132);
  h_muon_length_E->Draw("colz");
  
  c1->SaveAs((location+"/length_vs_E"+tag+".png").c_str());
  c1->SaveAs((location+"/length_vs_E"+tag+".root").c_str());
  c1->Clear();
  
  SetHistogramStyle2D(h_hits_muL_E, "True (generated) #mu energy [GeV]", "N_{Hits}/L_{#mu} [cm^{-1}]", false);
  h_hits_muL_E->GetZaxis()->SetLabelSize(0.03);
  h_hits_muL_E->GetZaxis()->SetLabelFont(132);
  h_hits_muL_E->Draw("colz");
  
  c1->SaveAs((location+"/hits_per_Lmu_vs_E"+tag+".png").c_str());
  c1->SaveAs((location+"/hits_per_Lmu_vs_E"+tag+".root").c_str());
  c1->Clear();
  
  SetHistogramStyle2D(h_secondaries_muL_E, "True (generated) #mu energy [GeV]", "N_{#delta}/L_{#mu} [cm^{-1}]", false);
  h_secondaries_muL_E->GetZaxis()->SetLabelSize(0.03);
  h_secondaries_muL_E->GetZaxis()->SetLabelFont(132);
  h_secondaries_muL_E->Draw("colz");
  
  c1->SaveAs((location+"/secondaries_per_Lmu_vs_E"+tag+".png").c_str());
  c1->SaveAs((location+"/secondaries_per_Lmu_vs_E"+tag+".root").c_str());
  c1->Clear();
  
  SetHistogramStyle2D(h_hits_E, "True (generated) #mu energy [GeV]", "N_{Hits}", false);
  h_hits_E->GetZaxis()->SetLabelSize(0.03);
  h_hits_E->GetZaxis()->SetLabelFont(132);
  h_hits_E->Draw("colz");
  
  c1->SaveAs((location+"/hits_vs_E"+tag+".png").c_str());
  c1->SaveAs((location+"/hits_vs_E"+tag+".root").c_str());
  c1->Clear();
  
  SetHistogramStyle2D(h_secondaries_E, "True (generated) #mu energy [GeV]", "N_{#delta}", false);
  h_secondaries_E->GetZaxis()->SetLabelSize(0.03);
  h_secondaries_E->GetZaxis()->SetLabelFont(132);
  h_secondaries_E->Draw("colz");
  
  c1->SaveAs((location+"/secondaries_vs_E"+tag+".png").c_str());
  c1->SaveAs((location+"/secondaries_vs_E"+tag+".root").c_str());
  c1->Clear();
  
  SetHistogramStyle2D(h_pfp_numDaughters_E, "True (generated) #mu energy [GeV]", "N_{Daughters}", false);
  h_pfp_numDaughters_E->GetZaxis()->SetLabelSize(0.03);
  h_pfp_numDaughters_E->GetZaxis()->SetLabelFont(132);
  h_pfp_numDaughters_E->Draw("colz");
  
  c1->SaveAs((location+"/pfp_numDaughters_vs_E"+tag+".png").c_str());
  c1->SaveAs((location+"/pfp_numDaughters_vs_E"+tag+".root").c_str());
  c1->Clear();
  
  SetHistogramStyle2D(h_pfp_numDauHits_E, "True (generated) #mu energy [GeV]", "N_{DaughterHits}", false);
  h_pfp_numDauHits_E->GetZaxis()->SetLabelSize(0.03);
  h_pfp_numDauHits_E->GetZaxis()->SetLabelFont(132);
  h_pfp_numDauHits_E->Draw("colz");
  
  c1->SaveAs((location+"/pfp_numDauHits_vs_E"+tag+".png").c_str());
  c1->SaveAs((location+"/pfp_numDauHits_vs_E"+tag+".root").c_str());
  c1->Clear();
  
  SetHistogramStyle2D(h_pfp_numDau_muL_E, "True (generated) #mu energy [GeV]", "N_{Daughters}/L_{#mu} [cm^{-1}]", false);
  h_pfp_numDau_muL_E->GetZaxis()->SetLabelSize(0.03);
  h_pfp_numDau_muL_E->GetZaxis()->SetLabelFont(132);
  h_pfp_numDau_muL_E->Draw("colz");
  
  c1->SaveAs((location+"/pfp_numDaughters_per_Lmu_vs_E"+tag+".png").c_str());
  c1->SaveAs((location+"/pfp_numDaughters_per_Lmu_vs_E"+tag+".root").c_str());
  c1->Clear();
  
  SetHistogramStyle2D(h_pfp_numDauHits_muL_E, "True (generated) #mu energy [GeV]", "N_{DaughterHits}/L_{#mu} [cm^{-1}]", false);
  h_pfp_numDauHits_muL_E->GetZaxis()->SetLabelSize(0.03);
  h_pfp_numDauHits_muL_E->GetZaxis()->SetLabelFont(132);
  h_pfp_numDauHits_muL_E->Draw("colz");
  
  c1->SaveAs((location+"/pfp_numDauHits_per_Lmu_vs_E"+tag+".png").c_str());
  c1->SaveAs((location+"/pfp_numDauHits_per_Lmu_vs_E"+tag+".root").c_str());
  c1->Clear();
  
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
