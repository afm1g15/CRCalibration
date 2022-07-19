/************************************************************************
 *
 * A macro to look at the statistical contributions to the CR muon
 * sample
 *
 * Example file list located here:
 *   /home/jones/work/cosmics/LArSoft-v08_50_00/work/files/v09_41_00_02_files.list
 *
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
   "ntracks_pandoraTrack",
   "trkId_pandoraTrack",
   "trkidtruth_pandoraTrack",
   "trkpdgtruth_pandoraTrack",
   "trkg4id_pandoraTrack",
   "ntrkhits_pandoraTrack",
   "trkresrg_pandoraTrack",
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
   "Eng",
   "pathlen"
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
     
int statisticalStudies(const char *config){

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
  
  // Start of analysis (loop over chain and events
  std::cout << " Running analysis..." << std::endl;

  // Then setup the histograms, counters and any other variables to add to
  // Setup histograms
  // dE/dx, all best plane
  TH2D *h_true_eDep_E = new TH2D("h_true_eDep_E",    "",100,4,5e3,100,0.2,10);   // True total energy deposition vs energy
  TH2D *h_true_dEdx_E = new TH2D("h_true_dEdx_E",    "",100,4,5e3,100,0.2,10);   // True hit energy depositions vs energy
  TH2D *h_reco_dEdx_E = new TH2D("h_reco_dEdx_E",    "",100,4,5e3,100,0.2,10);   // Reconstructed dE/dx vs energy
  SetLogX(h_true_eDep_E);
  SetLogX(h_true_dEdx_E);
  SetLogX(h_reco_dEdx_E);

  // dQ/dx, all best plane
  TH2D *h_true_qDep_E        = new TH2D("h_true_qDep_E",       "",100,4,5e3,100,0,1e3);   // Reconstructed dQ/dx vs energy
  TH2D *h_true_dQdx_E        = new TH2D("h_true_dQdx_E",       "",100,4,5e3,100,0,1e3);   // Reconstructed dQ/dx vs energy
  TH2D *h_reco_dQdx_E        = new TH2D("h_reco_dQdx_E",       "",100,4,5e3,100,0,1e3);   // Reconstructed dQ/dx vs energy
  TH2D *h_reco_dQdx_pos      = new TH2D("h_reco_dQdx_pos",     "",100,-800,800,100,0,1e3); // Reconstructed dQ/dx vs x position
  TH2D *h_reco_dQdx_dP       = new TH2D("h_reco_dQdx_dP",      "",100,0.3,1,100,0,1e3);    // Reconstructed dQ/dx vs pitch
  TH2D *h_reco_dQdx_RR       = new TH2D("h_reco_dQdx_RR",      "",100,0,1000,100,0,1e3);   // Reconstructed dQ/dx vs residual range
  TH2D *h_reco_dQdx_width    = new TH2D("h_reco_dQdx_width",   "",100,1,10,100,0,1e3);     // Reconstructed dQ/dx vs hit width
  TH2D *h_reco_dQdx_cosDrift = new TH2D("h_reco_dQdx_cosDrift","",100,-1,1,100,0,1e3);     // Reconstructed dQ/dx vs angle to drift direction (x)
  SetLogX(h_true_qDep_E);
  SetLogX(h_true_dQdx_E);
  SetLogX(h_reco_dQdx_E);
  
  // Setup counters
  unsigned int n_true_tpc_mu                   = 0;
  unsigned int n_true_tpc_mu_L                 = 0;
  unsigned int n_true_tpc_mu_L_E               = 0;
  unsigned int n_true_tpc_mu_L_E_cos           = 0;
  unsigned int n_true_tpc_mu_L_E_cos_thru      = 0;
  unsigned int n_true_tpc_mu_L_E_cos_thru_hitL = 0;
  
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

    // Initialise a vector of 'good G4' ID's 
    // Those which pass the defined cuts
    std::vector<int> goodG4;

    //
    // Truth-level studies
    //
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
      
      int pdg        = evt->pdg[iG4];
      int id         = evt->TrackId[iG4];
      float energy   = evt->Eng[iG4];
      float lengthAV = (endAV-vtxAV).Mag();
      float cosDrift = GetCosDrift(vtxAV,endAV); 

      // Make sure we are looking at a primary muon
      if(abs(pdg) != 13) continue;
      if(evt->Mother[iG4] != 0) continue;
      n_true_tpc_mu++;
      if(lengthAV < 300) continue;
      n_true_tpc_mu_L++;

      // Apply energy and angular requirements
      //  E > 6.504 GeV
      //  -0.567 < cosθDrift < 0.564
      if(energy < 6.504) continue;
      n_true_tpc_mu_L_E++;

      if(cosDrift < -0.567 || cosDrift > 0.564) continue;
      n_true_tpc_mu_L_E_cos++;

      // If we want only through-going muons, apply this here
      if(thru && !throughGoing) continue;
      n_true_tpc_mu_L_E_cos_thru++;
     
      // Now look more into the muons themselves, hits etc
      // Get the best plane
      // Boolean for each hit as to whether it is associated to the current G4 track
      // Create a vector for each plane too
      int bestPlane = 0;
      std::vector<int> hitsOnPlane(3,0);
      std::vector<std::vector<bool>> hitAssocOnPlane(3,std::vector<bool>(nHits,false));
      
      // Now fill those vectors
      GetNHitsOnPlane(id, nHits, evt, hitAssocOnPlane, hitsOnPlane);
      
      // Get the best plane
      bestPlane = std::max_element(hitsOnPlane.begin(), hitsOnPlane.end()) - hitsOnPlane.begin();
      
      // Sort out the minimum number of hits per unit length for a through-going muon
      double hitsPerL = hitsOnPlane.at(bestPlane)/static_cast<double>(lengthAV);
      if(hitsPerL < 0.8) continue;
      n_true_tpc_mu_L_E_cos_thru_hitL++;

      goodG4.push_back(id);
      
      // Now loop over hits 
      // Do everything on the best plane
      // Setup counters for the track-level quantities
      float totalEDep = 0.;
      float totalQDep = 0.;
      // Skip the first and last hits for 'reconstructability' purposes
      for(int iHit = 0; iHit < nHits-1; ++iHit){

        // Check if the current hit on the current plane is associated to the current track
        if(!hitAssocOnPlane.at(bestPlane).at(iHit)) continue;

        // Then get the parameters of interest for this hit
        float hitX = evt->hit_trueX[iHit];
        float hitE = evt->hit_energy[iHit];
        float hitQ = evt->hit_charge[iHit];

        // Check if x is lower than the APA bound, charge seems to accumulate there
        if(hitX < evtProc.APA_X_POSITIONS[0] || hitX > evtProc.APA_X_POSITIONS[2]) continue;

        // Sum up the current depositions
        totalEDep += hitE;
        totalQDep += hitQ;

        // Now fill the hit-based histograms in truth
        h_true_dEdx_E->Fill(energy,hitE);
        h_true_dQdx_E->Fill(energy,hitQ);
      } // Hit

      // Now fill the track-level histograms
      float totalEDepPerLength = totalEDep/static_cast<double>(lengthAV);
      float totalQDepPerLength = totalQDep/static_cast<double>(lengthAV);
      if(totalEDepPerLength < std::numeric_limits<float>::epsilon()) continue;
      if(totalQDepPerLength < std::numeric_limits<float>::epsilon()) continue;
      h_true_qDep_E->Fill(energy,totalQDepPerLength);
      h_true_eDep_E->Fill(energy,totalEDepPerLength);
    }// iG4

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
      CheckAndFlip(startVtx,endVtx);
      
      // Get the reconstructed best plane for this track
      int bestPlane = 0;
      std::vector<int> hitsOnPlane(3,0);
      GetRecoBestPlane(iTrk, evt, bestPlane, hitsOnPlane);

      // Make sure this track matches to a true track which passes the relevant cuts
      // Eventually should make the cuts in reco, but testing the dependence for now
      int trueID = evt->trkidtruth_pandoraTrack[iTrk][bestPlane];
      if(!CheckTrueIDAssoc(trueID,goodG4)) continue;
     
      // Make sure the reconstructed track is also long enough, for sanity
      if(!evtProc.SelectTrack(evt,iTrk)) continue;
      
      // Now get the true energy associated to this reconstructed track
      // and make sure it is physical
      double eng = GetTrueEnergyAssoc(iTrk,nGeant,evt,bestPlane);
      if(eng < 0) continue;

      // Only look at the best plane from now on
      // Make sure there is some information on the plane
      // Sometimes this is -999 and will break stuff
      int nHitsR = evt->ntrkhits_pandoraTrack[iTrk][bestPlane];
      if(nHitsR <= 0) continue;

      // Make sure it doesn't exceed the maximum energy of a true track: 1e5 GeV
      float energy = evt->trkke_pandoraTrack[iTrk][bestPlane]/1000.; // GeV
      if(energy > 1e5) continue;

      // Reset the nHitsR to not exceed the defined maximum
      if(nHitsR > MAX_TRACK_HITS) nHitsR = MAX_TRACK_HITS;

      // Now access the variables of interest
      Float_t *RRArr   = evt->trkresrg_pandoraTrack[iTrk][bestPlane];
      Float_t *dEdxArr = evt->trkdedx_pandoraTrack[iTrk][bestPlane];
      Float_t *dQdxArr = evt->trkdqdx_pandoraTrack[iTrk][bestPlane];

      // Convert them to vectors
      std::vector<float> RR(RRArr, RRArr + nHits);
      std::vector<float> dEdx(dEdxArr, dEdxArr + nHits);
      std::vector<float> dQdx(dQdxArr, dQdxArr + nHits);

      // Loop over the hits
      // And fill, removing first and last hit for completeness
      for(int iHit = 1; iHit < nHitsR-1; ++iHit){

        // Get the location of the current and following hits to determine the pitch
        TVector3 trkXYZ(evt->trkxyz_pandoraTrack[iTrk][bestPlane][iHit][0],
                        evt->trkxyz_pandoraTrack[iTrk][bestPlane][iHit][1],
                        evt->trkxyz_pandoraTrack[iTrk][bestPlane][iHit][2]);
        TVector3 nextXYZ(evt->trkxyz_pandoraTrack[iTrk][bestPlane][iHit+1][0],
                         evt->trkxyz_pandoraTrack[iTrk][bestPlane][iHit+1][1],
                         evt->trkxyz_pandoraTrack[iTrk][bestPlane][iHit+1][2]);

        float x         = trkXYZ.X();
        double dp       = GetHitPitch(bestPlane, trkXYZ, nextXYZ);
        double cosDrift = GetCosDrift(trkXYZ, nextXYZ);
          
        // Check if x is lower or higher than the APA bounds, charge seems to accumulate there
        if(x < evtProc.APA_X_POSITIONS[0] || x > evtProc.APA_X_POSITIONS[2]) continue;

        // Determine the lifetime corrections
        int tpc     = evtProc.WhichTPC(x) + 1;
        float dx    = ( -1 + 2*(tpc%2) )*(x - evtProc.APA_X_POSITIONS[tpc/2]);
        float dt    = dx*evtProc.kXtoT;
        float corr  = TMath::Exp(-dt/2.88);
        float eCorr = TMath::Exp(-dt/2.88) / TMath::Exp(-dt/3.); // Correct for the already-corrected energy

        // Convert with the lifetime correction
        float RRVal    = RR.at(iHit);
        float dEdxVal  = dEdx.at(iHit);
        float dQdxVal  = dQdx.at(iHit);
        float dEdxCorr = dEdxVal/eCorr;
        float dQdxCorr = dQdxVal/corr;
        float hitWidth = evt->hit_endT[iHit] - evt->hit_startT[iHit];

        // And fill the histograms
        h_reco_dEdx_E->Fill(eng,dEdxCorr);
        h_reco_dQdx_E->Fill(eng,dQdxCorr);
        h_reco_dQdx_pos->Fill(x,dQdxCorr);
        h_reco_dQdx_RR->Fill(RRVal,dQdxCorr);
        h_reco_dQdx_dP->Fill(dp,dQdxCorr);
        h_reco_dQdx_width->Fill(hitWidth,dQdxCorr);
        h_reco_dQdx_cosDrift->Fill(cosDrift,dQdxCorr);
      } // iHit
    } // iTrk
  }// Event loop

  // Draw histograms
  // First, setup the logged canvas for energy distributions
  TCanvas *c0 = new TCanvas("c0","",1000,800);
  SetCanvasStyle(c0, 0.1,0.15,0.05,0.12,0,0,0);
  c0->SetLogx();
  
  SetHistogramStyle2D(h_true_dEdx_E,"Muon energy [GeV]", "True dE/dx [MeV/cm]",false);
  h_true_dEdx_E->Scale(1,"width");
  h_true_dEdx_E->Draw("colz");
  h_true_dEdx_E->SaveAs((location+"/hist_true_dEdx_vs_E"+tag+".root").c_str());
  c0->SaveAs((location+"/true_dEdx_vs_E"+tag+".root").c_str());
  c0->SaveAs((location+"/true_dEdx_vs_E"+tag+".root").c_str());
  c0->Clear();
  
  SetHistogramStyle2D(h_true_dQdx_E,"Muon energy [GeV]", "True dQ/dx [ADC/cm]",false);
  h_true_dQdx_E->Scale(1,"width");
  h_true_dQdx_E->Draw("colz");
  h_true_dQdx_E->SaveAs((location+"/hist_true_dQdx_vs_E"+tag+".root").c_str());
  c0->SaveAs((location+"/true_dQdx_vs_E"+tag+".png").c_str());
  c0->SaveAs((location+"/true_dQdx_vs_E"+tag+".root").c_str());
  c0->Clear();
  
  SetHistogramStyle2D(h_true_eDep_E,"Muon energy [GeV]", "Total energy deposition per unit length [MeV/cm]",false);
  h_true_eDep_E->Scale(1,"width");
  h_true_eDep_E->Draw("colz");
  h_true_eDep_E->SaveAs((location+"/hist_true_eDep_vs_E"+tag+".root").c_str());
  c0->SaveAs((location+"/true_eDep_vs_E"+tag+".png").c_str());
  c0->SaveAs((location+"/true_eDep_vs_E"+tag+".root").c_str());
  c0->Clear();
  
  SetHistogramStyle2D(h_true_qDep_E,"Muon energy [GeV]", "Total charge deposition per unit length [ADC/cm]",false);
  h_true_qDep_E->Scale(1,"width");
  h_true_qDep_E->Draw("colz");
  h_true_qDep_E->SaveAs((location+"/hist_true_qDep_vs_E"+tag+".root").c_str());
  c0->SaveAs((location+"/true_qDep_vs_E"+tag+".png").c_str());
  c0->SaveAs((location+"/true_qDep_vs_E"+tag+".root").c_str());
  c0->Clear();
  
  SetHistogramStyle2D(h_reco_dEdx_E,"Muon energy [GeV]", "Reco dE/dx [MeV/cm]",false);
  h_reco_dEdx_E->Scale(1,"width");
  h_reco_dEdx_E->Draw("colz");
  h_reco_dEdx_E->SaveAs((location+"/hist_reco_dEdx_vs_E"+tag+".root").c_str());
  c0->SaveAs((location+"/reco_dEdx_vs_E"+tag+".png").c_str());
  c0->SaveAs((location+"/reco_dEdx_vs_E"+tag+".root").c_str());
  c0->Clear();
  
  SetHistogramStyle2D(h_reco_dQdx_E,"Muon energy [GeV]", "Reco dQ/dx [ADC/cm]",false);
  h_reco_dQdx_E->Draw("colz");
  h_reco_dQdx_E->SaveAs((location+"/hist_reco_dQdx_vs_E_noScale"+tag+".root").c_str());
  c0->SaveAs((location+"/reco_dQdx_vs_E_noScale"+tag+".png").c_str());
  c0->SaveAs((location+"/reco_dQdx_vs_E_noScale"+tag+".root").c_str());
  c0->Clear();
 
  SetHistogramStyle2D(h_reco_dQdx_E,"Muon energy [GeV]", "Reco dQ/dx [ADC/cm]",false);
  h_reco_dQdx_E->Scale(1,"width");
  h_reco_dQdx_E->Draw("colz");
  h_reco_dQdx_E->SaveAs((location+"/hist_reco_dQdx_vs_E"+tag+".root").c_str());
  c0->SaveAs((location+"/reco_dQdx_vs_E"+tag+".png").c_str());
  c0->SaveAs((location+"/reco_dQdx_vs_E"+tag+".root").c_str());
  c0->Clear();
 
  // Next, setup the canvas for the non-logged histograms
  TCanvas *c1 = new TCanvas("c1","",1000,800);
  SetCanvasStyle(c1, 0.1,0.15,0.05,0.12,0,0,0);
  
  SetHistogramStyle2D(h_reco_dQdx_pos,"X-Position [cm]", "Reco dQ/dx [ADC/cm]",false);
  h_reco_dQdx_pos->Scale(1,"width");
  h_reco_dQdx_pos->Draw("colz");
  h_reco_dQdx_pos->SaveAs((location+"/hist_reco_dQdx_vs_pos"+tag+".root").c_str());
  c1->SaveAs((location+"/reco_dQdx_vs_pos"+tag+".png").c_str());
  c1->SaveAs((location+"/reco_dQdx_vs_pos"+tag+".root").c_str());
  c1->Clear();
 
  SetHistogramStyle2D(h_reco_dQdx_RR,"Residual Range [cm]", "Reco dQ/dx [ADC/cm]",false);
  h_reco_dQdx_RR->Scale(1,"width");
  h_reco_dQdx_RR->Draw("colz");
  h_reco_dQdx_RR->SaveAs((location+"/hist_reco_dQdx_vs_RR"+tag+".root").c_str());
  c1->SaveAs((location+"/reco_dQdx_vs_RR"+tag+".png").c_str());
  c1->SaveAs((location+"/reco_dQdx_vs_RR"+tag+".root").c_str());
  c1->Clear();
 
  SetHistogramStyle2D(h_reco_dQdx_dP,"Hit Pitch [cm]", "Reco dQ/dx [ADC/cm]",false);
  h_reco_dQdx_dP->Scale(1,"width");
  h_reco_dQdx_dP->Draw("colz");
  h_reco_dQdx_dP->SaveAs((location+"/hist_reco_dQdx_vs_dP"+tag+".root").c_str());
  c1->SaveAs((location+"/reco_dQdx_vs_dP"+tag+".png").c_str());
  c1->SaveAs((location+"/reco_dQdx_vs_dP"+tag+".root").c_str());
  c1->Clear();
 
  SetHistogramStyle2D(h_reco_dQdx_width,"Hit Width [ticks]", "Reco dQ/dx [ADC/cm]",false);
  h_reco_dQdx_width->Scale(1,"width");
  h_reco_dQdx_width->Draw("colz");
  h_reco_dQdx_width->SaveAs((location+"/hist_reco_dQdx_vs_width"+tag+".root").c_str());
  c1->SaveAs((location+"/reco_dQdx_vs_width"+tag+".png").c_str());
  c1->SaveAs((location+"/reco_dQdx_vs_width"+tag+".root").c_str());
  c1->Clear();
 
  SetHistogramStyle2D(h_reco_dQdx_cosDrift,"cos#theta_{Drift}", "Reco dQ/dx [ADC/cm]",false);
  h_reco_dQdx_cosDrift->Scale(1,"width");
  h_reco_dQdx_cosDrift->Draw("colz");
  h_reco_dQdx_cosDrift->SaveAs((location+"/hist_reco_dQdx_vs_cosDrift"+tag+".root").c_str());
  c1->SaveAs((location+"/reco_dQdx_vs_cosDrift"+tag+".png").c_str());
  c1->SaveAs((location+"/reco_dQdx_vs_cosDrift"+tag+".root").c_str());
  c1->Clear();
 
  // Write the stats out to a file
  ofstream oFile;
  oFile.open((location+"/statistics"+tag+".txt").c_str());

  oFile << " -----------------------------------------------------------------------------------------------------------------" << std::endl;
  oFile << " TRUTH " << std::endl;
  oFile << " -----------------------------------------------------------------------------------------------------------------" << std::endl;
  oFile << "  Total TPC muons true            : " << std::setw(8) << n_true_tpc_mu       
        << " = "<< std::setprecision(3) << std::setw(8) << 100*n_true_tpc_mu/static_cast<double>(n_true_tpc_mu) 
        << " % true tpc muons" << std::endl;
  oFile << "     & L > 3m                     : " << std::setw(8) << n_true_tpc_mu_L 
        << " = "<< std::setprecision(3) << std::setw(8) << 100*n_true_tpc_mu_L/static_cast<double>(n_true_tpc_mu) 
        << " % true tpc muons" << std::endl;
  oFile << "     & E > 6.504 GeV              : " << std::setw(8) << n_true_tpc_mu_L_E     
        << " = "<< std::setprecision(3) << std::setw(8) << 100*n_true_tpc_mu_L_E/static_cast<double>(n_true_tpc_mu) 
        << " % true tpc muons" << std::endl;
  oFile << "     & -0.567 < cosθDrift < 0.564 : " << std::setw(8) << n_true_tpc_mu_L_E_cos 
        << " = "<< std::setprecision(3) << std::setw(8) << 100*n_true_tpc_mu_L_E_cos/static_cast<double>(n_true_tpc_mu) 
        << " % true tpc muons" << std::endl;
  oFile << "     & through                    : " << std::setw(8) << n_true_tpc_mu_L_E_cos_thru
        << " = "<< std::setprecision(3) << std::setw(8) << 100*n_true_tpc_mu_L_E_cos_thru/static_cast<double>(n_true_tpc_mu) 
        << " % true tpc muons" << std::endl;
  oFile << "     & hits/L > 0.8               : " << std::setw(8) << n_true_tpc_mu_L_E_cos_thru_hitL
        << " = "<< std::setprecision(3) << std::setw(8) << 100*n_true_tpc_mu_L_E_cos_thru_hitL/static_cast<double>(n_true_tpc_mu) 
        << " % true tpc muons" << std::endl;
  oFile << " -----------------------------------------------------------------------------------------------------------------" << std::endl;
  oFile << " RECONSTRUCTION " << std::endl;
  oFile << " -----------------------------------------------------------------------------------------------------------------" << std::endl;

  oFile.close();


  std::cout << " --- 100 % --- |" << std::endl;

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
