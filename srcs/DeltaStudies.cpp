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

#include "Classes/EventProcessor.h"
#include "Classes/ConfigReader.h"

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
   "trkdqdx_pandoraTrack",
   "trkresrg_pandoraTrack",
   "trkdedx_pandoraTrack",
   "nPFParticles",
   "pfp_selfID",
   "pfp_isPrimary",
   "pfp_isTrack",
   "pfp_trackID",
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
  TH1D *h_muon_length               = new TH1D("h_muon_length","",200,0,2000); // Muon length
  TH1D *h_delta_length              = new TH1D("h_delta_length","",200,0,40); // Secondary length
  TH1D *h_pfp_numDaughters          = new TH1D("h_pfp_numDaughters","",100,0,100); // Number of primary pfparticle daughters
  TH1D *h_pfp_numDauHitsTot         = new TH1D("h_pfp_numDauHitsTot","",200,0,3000); // Number of primary pfparticle daughter hits in total for the muon
  TH1D *h_pfp_numDauHits            = new TH1D("h_pfp_numDauHits","",25,0,25); // Number of primary pfparticle daughter hits
  TH1D *h_pfp_numDaughters_thresh   = new TH1D("h_pfp_numDaughters_thresh","",100,0,100); // Number of primary pfparticle daughters
  TH1D *h_pfp_numDauHits_thresh     = new TH1D("h_pfp_numDauHits_thresh","",100,0,100); // Number of primary pfparticle daughter hits
  TH1D *h_hits_muL                  = new TH1D("h_hits_muL","",100,0,0.4); // Number of hits per unit length of the primary muon track
  TH1D *h_secondaries_muL           = new TH1D("h_secondaries_muL","",100,0,0.05); // Number of secondaries per unit length of the primary muon track
  
  TH2D *h_muon_length_E               = new TH2D("h_muon_length_E","",100,5,5e3,100,0,2000); // Muon length
  TH2D *h_hits_muL_E                  = new TH2D("h_hits_muL_E","",100,5,5e3,100,0,0.12);
  TH2D *h_hits_E                      = new TH2D("h_hits_E","",100,5,5e3,80,0,240);
  TH2D *h_secondaries_muL_E           = new TH2D("h_secondaries_muL_E","",100,5,5e3,100,0,0.03);
  TH2D *h_secondaries_E               = new TH2D("h_secondaries_E","",100,5,5e3,50,0,50);
  TH2D *h_pfp_numDaughters_E          = new TH2D("h_pfp_numDaughters_E","",100,5,5e3,80,0,80); // Number of primary pfparticle daughters vs energy
  TH2D *h_pfp_numDau_muL_E            = new TH2D("h_pfp_numDau_muL_E","",100,5,5e3,80,0,0.1); // Number of primary pfparticle daughters vs energy
  TH2D *h_pfp_numDauHitsTot_E         = new TH2D("h_pfp_numDauHitsTot_E","",100,5,5e3,100,0,3000); // Number of primary pfparticle total daughter hits vs energy
  TH2D *h_pfp_numDauHits_E            = new TH2D("h_pfp_numDauHits_E","",100,5,5e3,25,0,25); // Number of primary pfparticle daughter hits vs energy
  TH2D *h_pfp_numDauHits_muL_E        = new TH2D("h_pfp_numDauHits_muL_E","",100,5,5e3,100,0,5); // Number of primary pfparticle daughters vs energy
  TH2D *h_pfp_numDaughters_thresh_E   = new TH2D("h_pfp_numDaughters_thresh_E","",100,5,5e3,40,0,40); // Number of primary pfparticle daughters vs energy
  TH2D *h_pfp_numDau_muL_thresh_E     = new TH2D("h_pfp_numDau_muL_thresh_E","",100,5,5e3,80,0,0.04); // Number of primary pfparticle daughters vs energy
  TH2D *h_pfp_numDauHits_thresh_E     = new TH2D("h_pfp_numDauHits_thresh_E","",100,5,5e3,50,0,50); // Number of primary pfparticle daughter hits vs energy
  TH2D *h_pfp_numDauHits_muL_thresh_E = new TH2D("h_pfp_numDauHits_muL_thresh_E","",100,5,5e3,100,0,5); // Number of primary pfparticle daughters vs energy
  TH2D *h_corr_dqdx_E                 = new TH2D("h_corr_dqdx_E","",300,5,5e3,300,80,600);
  TH2D *h_dedx_rr_lt30GeV             = new TH2D("h_dedx_rr_lt30GeV","",200,0,1500,200,0,7); // dE/dx vs RR in the last 50cm of a track (below 30GeV)
  TH2D *h_dedx_rr_gt30GeV             = new TH2D("h_dedx_rr_gt30GeV","",200,0,1500,200,0,7); // dE/dx vs RR in the last 50cm of a track (above 30GeV)
 
  // Sort out log scales if needed 
  SetLogX(h_muon_length_E);
  SetLogX(h_hits_muL_E);
  SetLogX(h_hits_E);
  SetLogX(h_secondaries_muL_E);
  SetLogX(h_secondaries_E);
  SetLogX(h_pfp_numDaughters_E);
  SetLogX(h_pfp_numDauHitsTot_E);
  SetLogX(h_pfp_numDauHits_E);
  SetLogX(h_pfp_numDau_muL_E);
  SetLogX(h_pfp_numDauHits_muL_E);
  SetLogX(h_pfp_numDaughters_thresh_E);
  SetLogX(h_pfp_numDauHits_thresh_E);
  SetLogX(h_pfp_numDau_muL_thresh_E);
  SetLogX(h_pfp_numDauHits_muL_thresh_E);
  SetLogX(h_corr_dqdx_E);
  
  // Setup counters
  unsigned int nPFParticlesDeltaCut = 0;
  unsigned int nPFParticles = 0;

  unsigned int nBelow30     = 0;
  unsigned int nAbove30     = 0;
  unsigned int nBelow30DCut = 0;
  unsigned int nAbove30DCut = 0;
 
  unsigned int noDeposit    = 0;

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
    int primaryID = -999;
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
      
      // Set the primary ID
      primaryID = ID;

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
        if(primaryTrueEnergy < 30)
          nBelow30++;
        else
          nAbove30++;

        // Now loop over hits and fill dE/dx vs residual range in each energy case (>/< 30 GeV)
        // Make sure it doesn't exceed the maximum size of the array
        // Count if it does so we can see how often it happens
        int nHitsR = evt->ntrkhits_pandoraTrack[iTrk][iPlane];
        if(nHitsR <= 0){
          continue;
        }

        Float_t *RRArr   = evt->trkresrg_pandoraTrack[iTrk][iPlane];
        Float_t *dEdxArr = evt->trkdedx_pandoraTrack[iTrk][iPlane];

        // Convert them to vectors
        std::vector<float> RR(RRArr, RRArr + nHitsR);
        std::vector<float> dEdx(dEdxArr, dEdxArr + nHitsR);
        
        // And fill, removing first and last hit for completeness
        for(int iHit = 1; iHit < nHitsR-1; ++iHit){

          // General geometry of the track
          // Check if x is lower or higher than the APA bounds, charge seems to accumulate there
          float x = evt->trkxyz_pandoraTrack[iTrk][bestPlane][iHit][0];
          if(x < evtProc.APA_X_POSITIONS[0] || x > evtProc.APA_X_POSITIONS[2]) continue;

          // Lifetime correction
          int tpc     = evtProc.WhichTPC(x) + 1;
          float dx    = ( -1 + 2*(tpc%2) )*(x - evtProc.APA_X_POSITIONS[tpc/2]);
          float dt    = dx*evtProc.kXtoT;
          float eCorr = TMath::Exp(-dt/2.88) / TMath::Exp(-dt/3.); // Correct for the already-corrected energy

          // New values
          float RRVal     = RR.at(iHit);
          float dEdxVal   = dEdx.at(iHit);
          float dEdxCorr  = dEdxVal/eCorr;
          if(dEdxCorr < std::numeric_limits<float>::epsilon()){
            noDeposit++;
            continue;
          }

          // Now fill the histograms
          if(primaryTrueEnergy < 30)
            h_dedx_rr_lt30GeV->Fill(RRVal,dEdxCorr);
          else
            h_dedx_rr_gt30GeV->Fill(RRVal,dEdxCorr);

        } // Hits
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

    // Make sure we found a primary
    if(primaryID < 0) continue;
    
    // Now we have identified the primary PFP, get its daughters and corresponding hits
    // Loop over the PFParticles, match to the trackID of the identified muon and get the daughter information
    for(unsigned int iPfp = 0; iPfp < nPfps; ++iPfp){
      int PFPTrackID = evt->pfp_trackID[iPfp];
      int PFPID = evt->pfp_selfID[iPfp];
      if(PFPTrackID != primaryID) continue;

      unsigned int nDau = evt->pfp_numDaughters[iPfp];
      if(nDau == 0) continue;
      unsigned int nDauLim = nDau;
      if(nDau > MAX_PFPDAUGHTERS){
        nDauLim = MAX_PFPDAUGHTERS;
      }
      // Now loop over daughters and get the number of hits associated with it (HACK)
      unsigned int nDauHits   = 0;
      unsigned int nDauHitsCut = 0;
      unsigned int nDauCut     = 0;
      for(unsigned int iDau = 0; iDau < nDauLim; ++iDau){
        int dauID = evt->pfp_daughterIDs[iPfp][iDau];
        bool dauHitThreshold = false;

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
                h_pfp_numDauHits->Fill(evt->cluster_NHits[iclu]);
                h_pfp_numDauHits_E->Fill(primaryTrueEnergy,evt->cluster_NHits[iclu]);

                // Now check if the number of daughter hits is above threshold
                if(evt->cluster_NHits[iclu] > 6){
                  dauHitThreshold = true;
                  nDauHitsCut += evt->cluster_NHits[iclu];
                  h_pfp_numDauHits_thresh->Fill(evt->cluster_NHits[iclu]);
                  h_pfp_numDauHits_thresh_E->Fill(primaryTrueEnergy,evt->cluster_NHits[iclu]);
                }
              }
            } // clusters
          } // PFP clusters
        } // jPfp
        if(dauHitThreshold) nDauCut++;
      } // Daughters
      if(nDauHits == 0) continue;
      nPFParticles++;
      float daughtersPerL    = nDau/primaryLength;
      float daughterHitsPerL = nDauHits/primaryLength;
      h_pfp_numDaughters->Fill(nDau);
      h_pfp_numDauHitsTot->Fill(nDauHits);
      h_pfp_numDaughters_E->Fill(primaryTrueEnergy,nDau);
      h_pfp_numDauHitsTot_E->Fill(primaryTrueEnergy,nDauHits);
      h_pfp_numDau_muL_E->Fill(primaryTrueEnergy,daughtersPerL);
      h_pfp_numDauHits_muL_E->Fill(primaryTrueEnergy,daughterHitsPerL);
      
      // Now for those above the hit threshold
      if(nDauHitsCut == 0) continue;
      float daughtersPerLCut    = nDauCut/primaryLength;
      float daughterHitsPerLCut = nDauHitsCut/primaryLength;
      h_pfp_numDaughters_thresh->Fill(nDauCut);
      h_pfp_numDaughters_thresh_E->Fill(primaryTrueEnergy,nDauCut);
      h_pfp_numDau_muL_thresh_E->Fill(primaryTrueEnergy,daughtersPerLCut);
      h_pfp_numDauHits_muL_thresh_E->Fill(primaryTrueEnergy,daughterHitsPerLCut);

      // Now that the studies have been performed,
      // make the dQ/dx distributions with a cut on Ndaughters/Lmu < 17e-3
      // to try and remove the low energy-dependence
      if(daughtersPerL < 17e-3) continue;
      nPFParticlesDeltaCut++;
      
      // Now loop over the tracks so we can do some stuff!!
      for(unsigned int iTrk = 0; iTrk < nTrks; ++iTrk){
        // Make sure we've got the right track
        int ID = evt->trkId_pandoraTrack[iTrk];
        if(ID != primaryID) continue;

        // Setup list of plane labels the track has crossed
        std::vector<std::string> labelsCrossed;

        // Get the best plane
        int bestPlane = 0;
        int currHits  = -999;
        for(int iPlane = 0; iPlane < 3; ++iPlane){
          if(evt->ntrkhits_pandoraTrack[iTrk][iPlane] > currHits){
            currHits  = evt->ntrkhits_pandoraTrack[iTrk][iPlane];
            bestPlane = iPlane; 
          } // CurrHits
        } // Planes


        // Only look at primary muons
        if(abs(evt->trkpdgtruth_pandoraTrack[iTrk][bestPlane]) != 13) continue;

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
        // Temporarily remove these tracks for energy-deposition's sake
        if(startVtx.Y() < endVtx.Y()){
          TVector3 temp(endVtx);
          endVtx = startVtx;
          startVtx = temp;
        }
        
        float length = evt->trklen_pandoraTrack[iTrk];
        bool thruGoing = IsThroughGoing(length,startVtx,endVtx,extPlanes,fidExtPlanes);
        
        unsigned int nHits = evt->ntrkhits_pandoraTrack[iTrk][bestPlane];

        // Get the associated true energy of the muon
        // Try to get the true energy
        // Get the list iterator from matching ID's
        float eng = -1.;
        for(int iG4 = 0; iG4 < nGeant; ++iG4){
          int trueID = evt->TrackId[iG4];

          if(evt->trkidtruth_pandoraTrack[iTrk][bestPlane] == trueID)
            eng = evt->Eng[iG4];
        }
        if(eng < 0){
          std::cout << " Warning: Energy is below zero, skipping track with energy: " << eng << std::endl;
          continue;
        }

        if(eng < 30)
          nBelow30DCut++;
        else
          nAbove30DCut++;
     
        // Make sure it doesn't exceed the maximum size of the array
        // Count if it does so we can see how often it happens
        if(nHits > MAX_TRACK_HITS){
          nHits = MAX_TRACK_HITS;
        }
        Float_t *dQdxArr = evt->trkdqdx_pandoraTrack[iTrk][bestPlane];
        std::vector<float> dQdx(dQdxArr, dQdxArr + nHits);
        
        for(unsigned int iHit = 0; iHit < nHits; ++iHit){

          // General geometry of the track
          // Check if x is lower or higher than the APA bounds, charge seems to accumulate there
          float x = evt->trkxyz_pandoraTrack[iTrk][bestPlane][iHit][0];
          if(x < evtProc.APA_X_POSITIONS[0] || x > evtProc.APA_X_POSITIONS[2]) continue;

          // Lifetime correction
          int tpc  = evtProc.WhichTPC(x) + 1;
          float dx = ( -1 + 2*(tpc%2) )*(x - evtProc.APA_X_POSITIONS[tpc/2]);
          float dt = dx*evtProc.kXtoT;
          float corr  = TMath::Exp(-dt/2.88);

          // New values
          float dQdxVal   = dQdx.at(iHit);
          float dQdxCorr  = dQdxVal/corr;

          // the following studies should be conducted with top-bottom muons to start with
          if(thru == 1 && !thruGoing) continue;
          h_corr_dqdx_E->Fill(eng,dQdxCorr);
        } // Hits
      } // Tracks
    } // PFParticles
  }// Event loop
  std::cout << " --- 100 % --- |" << std::endl;

  std::cout << " Number of PFParticles:         " << nPFParticles << std::endl;
  std::cout << " PFParticles after delta cut:   " << nPFParticlesDeltaCut << std::endl;

  std::cout << " Number below 30: " << nBelow30 << ", % remaining after delta-cut: " << (nBelow30DCut/static_cast<double>(nBelow30))*100 << std::endl;
  std::cout << " Number above 30: " << nAbove30 << ", % remaining after delta-cut: " << (nAbove30DCut/static_cast<double>(nAbove30))*100 << std::endl;

  std::cout << "Number of hits which don't deposit anything: " << noDeposit << std::endl;

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
  
  SetHistogramStyle1D(h_pfp_numDaughters_thresh,"N_{Daughters}", " Rate");
  h_pfp_numDaughters_thresh->Draw("hist");
  h_pfp_numDaughters_thresh->SetLineWidth(2);
  h_pfp_numDaughters_thresh->SetLineColor(kViolet-5);
  h_pfp_numDaughters_thresh->GetYaxis()->SetTitleOffset(0.95);
  
  c0->SaveAs((location+"/pfp_numDaughters_thresh"+tag+".png").c_str());
  c0->SaveAs((location+"/pfp_numDaughters_thresh"+tag+".root").c_str());
  c0->Clear();
  
  SetHistogramStyle1D(h_pfp_numDauHitsTot,"N_{DaughterHits}", " Rate");
  h_pfp_numDauHitsTot->Draw("hist");
  h_pfp_numDauHitsTot->SetLineWidth(2);
  h_pfp_numDauHitsTot->SetLineColor(kViolet-5);
  h_pfp_numDauHitsTot->GetYaxis()->SetTitleOffset(0.95);
  
  c0->SaveAs((location+"/pfp_numDaughterHitsTot"+tag+".png").c_str());
  c0->SaveAs((location+"/pfp_numDaughterHitsTot"+tag+".root").c_str());
  c0->Clear();
  
  SetHistogramStyle1D(h_pfp_numDauHits,"N_{Hits}/Daughter", " Rate");
  h_pfp_numDauHits->Draw("hist");
  h_pfp_numDauHits->SetLineWidth(2);
  h_pfp_numDauHits->SetLineColor(kViolet-5);
  h_pfp_numDauHits->GetYaxis()->SetTitleOffset(0.95);
  
  c0->SaveAs((location+"/pfp_numDaughterHits"+tag+".png").c_str());
  c0->SaveAs((location+"/pfp_numDaughterHits"+tag+".root").c_str());
  c0->Clear();
  
  SetHistogramStyle1D(h_pfp_numDauHits_thresh,"N_{Hits}/Daughter", " Rate");
  h_pfp_numDauHits_thresh->Draw("hist");
  h_pfp_numDauHits_thresh->SetLineWidth(2);
  h_pfp_numDauHits_thresh->SetLineColor(kViolet-5);
  h_pfp_numDauHits_thresh->GetYaxis()->SetTitleOffset(0.95);
  
  c0->SaveAs((location+"/pfp_numDaughterHits_thresh"+tag+".png").c_str());
  c0->SaveAs((location+"/pfp_numDaughterHits_thresh"+tag+".root").c_str());
  c0->Clear();
  
  TCanvas *c1 = new TCanvas("c1","",1000,800);
  SetCanvasStyle(c1, 0.12,0.12,0.05,0.12,0,0,0);

  SetHistogramStyle2D(h_dedx_rr_gt30GeV,"Residual range [cm]", " Reconstructed dE/dx [MeV/cm]", false);
  h_dedx_rr_gt30GeV->GetZaxis()->SetLabelSize(0.03);
  h_dedx_rr_gt30GeV->GetZaxis()->SetLabelFont(132);
  h_dedx_rr_gt30GeV->Draw("colz");

  c1->SaveAs((location+"/dEdx_vs_RR_gt30GeV"+tag+".png").c_str());
  c1->SaveAs((location+"/dEdx_vs_RR_gt30GeV"+tag+".root").c_str());
  c1->Clear();

  SetHistogramStyle2D(h_dedx_rr_lt30GeV,"Residual range [cm]", " Reconstructed dE/dx [MeV/cm]", false);
  h_dedx_rr_lt30GeV->GetZaxis()->SetLabelSize(0.03);
  h_dedx_rr_lt30GeV->GetZaxis()->SetLabelFont(132);
  h_dedx_rr_lt30GeV->Draw("colz");

  c1->SaveAs((location+"/dEdx_vs_RR_lt30GeV"+tag+".png").c_str());
  c1->SaveAs((location+"/dEdx_vs_RR_lt30GeV"+tag+".root").c_str());
  c1->Clear();
  
  c1->SetLogx();
  SetHistogramStyle2D(h_corr_dqdx_E,"True (generated) #mu energy [GeV]", " dQ/dx [ADC/cm]", false);
  h_corr_dqdx_E->GetZaxis()->SetLabelSize(0.03);
  h_corr_dqdx_E->GetZaxis()->SetLabelFont(132);
  h_corr_dqdx_E->Draw("colz");

  h_corr_dqdx_E->SaveAs((location+"/hist_delta_cut_corr_charge_vs_E"+tag+".root").c_str());
  c1->SaveAs((location+"/delta_cut_corr_charge_vs_E"+tag+".png").c_str());
  c1->SaveAs((location+"/delta_cut_corr_charge_vs_E"+tag+".root").c_str());
  c1->Clear();

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
  
  SetHistogramStyle2D(h_pfp_numDauHitsTot_E, "True (generated) #mu energy [GeV]", "N_{DaughterHits}", false);
  h_pfp_numDauHitsTot_E->GetZaxis()->SetLabelSize(0.03);
  h_pfp_numDauHitsTot_E->GetZaxis()->SetLabelFont(132);
  h_pfp_numDauHitsTot_E->Draw("colz");
  
  c1->SaveAs((location+"/pfp_numDauHitsTot_vs_E"+tag+".png").c_str());
  c1->SaveAs((location+"/pfp_numDauHitsTot_vs_E"+tag+".root").c_str());
  c1->Clear();
  
  SetHistogramStyle2D(h_pfp_numDaughters_E, "True (generated) #mu energy [GeV]", "N_{Daughters}", false);
  h_pfp_numDaughters_E->GetZaxis()->SetLabelSize(0.03);
  h_pfp_numDaughters_E->GetZaxis()->SetLabelFont(132);
  h_pfp_numDaughters_E->Draw("colz");
  
  c1->SaveAs((location+"/pfp_numDaughters_vs_E"+tag+".png").c_str());
  c1->SaveAs((location+"/pfp_numDaughters_vs_E"+tag+".root").c_str());
  c1->Clear();
  
  SetHistogramStyle2D(h_pfp_numDaughters_thresh_E, "True (generated) #mu energy [GeV]", "N_{Daughters}", false);
  h_pfp_numDaughters_thresh_E->GetZaxis()->SetLabelSize(0.03);
  h_pfp_numDaughters_thresh_E->GetZaxis()->SetLabelFont(132);
  h_pfp_numDaughters_thresh_E->Draw("colz");
  
  c1->SaveAs((location+"/pfp_numDaughters_vs_thresh_E"+tag+".png").c_str());
  c1->SaveAs((location+"/pfp_numDaughters_vs_thresh_E"+tag+".root").c_str());
  c1->Clear();
  

  SetHistogramStyle2D(h_pfp_numDauHits_E, "True (generated) #mu energy [GeV]", "N_{Hits}/Daughter", false);
  h_pfp_numDauHits_E->GetZaxis()->SetLabelSize(0.03);
  h_pfp_numDauHits_E->GetZaxis()->SetLabelFont(132);
  h_pfp_numDauHits_E->Draw("colz");
  
  c1->SaveAs((location+"/pfp_numDauHits_vs_E"+tag+".png").c_str());
  c1->SaveAs((location+"/pfp_numDauHits_vs_E"+tag+".root").c_str());
  c1->Clear();
  
  SetHistogramStyle2D(h_pfp_numDauHits_thresh_E, "True (generated) #mu energy [GeV]", "N_{Hits}/Daughter", false);
  h_pfp_numDauHits_thresh_E->GetZaxis()->SetLabelSize(0.03);
  h_pfp_numDauHits_thresh_E->GetZaxis()->SetLabelFont(132);
  h_pfp_numDauHits_thresh_E->Draw("colz");
  
  c1->SaveAs((location+"/pfp_numDauHits_vs_thresh_E"+tag+".png").c_str());
  c1->SaveAs((location+"/pfp_numDauHits_vs_thresh_E"+tag+".root").c_str());
  c1->Clear();
  
  SetHistogramStyle2D(h_pfp_numDau_muL_E, "True (generated) #mu energy [GeV]", "N_{Daughters}/L_{#mu} [cm^{-1}]", false);
  h_pfp_numDau_muL_E->GetZaxis()->SetLabelSize(0.03);
  h_pfp_numDau_muL_E->GetZaxis()->SetLabelFont(132);
  h_pfp_numDau_muL_E->Draw("colz");
  
  c1->SaveAs((location+"/pfp_numDaughters_per_Lmu_vs_E"+tag+".png").c_str());
  c1->SaveAs((location+"/pfp_numDaughters_per_Lmu_vs_E"+tag+".root").c_str());
  c1->Clear();
  
  SetHistogramStyle2D(h_pfp_numDau_muL_thresh_E, "True (generated) #mu energy [GeV]", "N_{Daughters}/L_{#mu} [cm^{-1}]", false);
  h_pfp_numDau_muL_thresh_E->GetZaxis()->SetLabelSize(0.03);
  h_pfp_numDau_muL_thresh_E->GetZaxis()->SetLabelFont(132);
  h_pfp_numDau_muL_thresh_E->Draw("colz");
  
  c1->SaveAs((location+"/pfp_numDaughters_per_Lmu_vs_thresh_E"+tag+".png").c_str());
  c1->SaveAs((location+"/pfp_numDaughters_per_Lmu_vs_thresh_E"+tag+".root").c_str());
  c1->Clear();
  
  SetHistogramStyle2D(h_pfp_numDauHits_muL_E, "True (generated) #mu energy [GeV]", "N_{DaughterHits}/L_{#mu} [cm^{-1}]", false);
  h_pfp_numDauHits_muL_E->GetZaxis()->SetLabelSize(0.03);
  h_pfp_numDauHits_muL_E->GetZaxis()->SetLabelFont(132);
  h_pfp_numDauHits_muL_E->Draw("colz");
  
  c1->SaveAs((location+"/pfp_numDauHits_per_Lmu_vs_E"+tag+".png").c_str());
  c1->SaveAs((location+"/pfp_numDauHits_per_Lmu_vs_E"+tag+".root").c_str());
  c1->Clear();
  
  SetHistogramStyle2D(h_pfp_numDauHits_muL_thresh_E, "True (generated) #mu energy [GeV]", "N_{DaughterHits}/L_{#mu} [cm^{-1}]", false);
  h_pfp_numDauHits_muL_thresh_E->GetZaxis()->SetLabelSize(0.03);
  h_pfp_numDauHits_muL_thresh_E->GetZaxis()->SetLabelFont(132);
  h_pfp_numDauHits_muL_thresh_E->Draw("colz");
  
  c1->SaveAs((location+"/pfp_numDauHits_per_Lmu_vs_thresh_E"+tag+".png").c_str());
  c1->SaveAs((location+"/pfp_numDauHits_per_Lmu_vs_thresh_E"+tag+".root").c_str());
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
