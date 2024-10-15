/************************************************************************
 * 
 * A macro to make an event selection on true and reco stopping muons for 
 * dE/dx calibration studies
 *
 *
 * Input is a list of ana files.
 * Example file list located here:
 *   /exp/dune/app/users/amoor/duneCalibration/anafiles.list
 *
 *
 *************************************************************************/

#include "EventProcessor.h"
#include "ConfigReader.h"
#include "TTree.h"

using namespace calib;
using namespace cppsecrets;

// Allowed branches to read from the tree (all branches in anatree_core.h)
std::vector<TString> allowed = {
   "run",
   "event",
   "geant_list_size",
   "trkpdgtruth_pandoraTrack",
   "trkg4id_pandoraTrack",
   "trkId_pandoraTrack",
   "trkidtruth_pandoraTrack",
   "ntracks_pandoraTrack",
   "trkId_pandoraTrack",
   "ntrkhits_pandoraTrack",
   "trkdqdx_pandoraTrack",
   "trkdedx_pandoraTrack",
   "trkresrg_pandoraTrack",
   "trkxyz_pandoraTrack",
   "trkstartx_pandoraTrack",
   "trkstarty_pandoraTrack",
   "trkstartz_pandoraTrack",
   "trkendx_pandoraTrack",
   "trkendy_pandoraTrack",
   "trkendz_pandoraTrack",
   "trklen_pandoraTrack",
   "pdg",
   "Mother",
   "EndPointx_tpcAV",
   "EndPointy_tpcAV",
   "EndPointz_tpcAV",
   "EndPointx",
   "EndPointy",
   "EndPointz",
   "TrackId",
   "trkthetaxz_pandoraTrack",
   "trkthetayz_pandoraTrack",
   "trkpurity_pandoraTrack",
   "nvtx_pandora",
   "trkcompleteness_pandoraTrack",
   "ntrkhits_pandoraTrack",
   "trkstartd_pandoraTrack",
   "trkmom_pandoraTrack"
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
     
int stoppingMuonStudy(const char *config){

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

  // Get configuration variables and initiate the relevant ones
  int n = -1;   // How many files from the file list to run. Default: All (-1)
  int thru = 0; // Do we want to select only through-going muons? Default: No (0)
  int stop = 0; // Do we want to select only stopping muons? Default: No (0)
  std::string input_list = "";
  std::string location="";
  std::string tag="";
  std::vector<double> minx_fid, miny_fid, minz_fid;
  std::vector<double> maxx_fid, maxy_fid, maxz_fid;
  std::vector<double> minx_av, miny_av, minz_av;
  std::vector<double> maxx_av, maxy_av, maxz_av;

  // Access corresponding parameter in the configuration file
  p->getValue("InputList", input_list);
  p->getValue("Location",  location);
  p->getValue("Tag",       tag);
  p->getValue("NFiles",    n);
  p->getValue("Thru",      thru);
  p->getValue("Stopping",  stop);
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

  // Sanity check the geometry definitions
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Total number of planes in the active volume of the DUNE SP module: " << allPlanes.size() << std::endl;
  std::cout << " Consisting of " << extPlanes.size() << " external planes and " << intPlanes.size() << " internal planes" << std::endl; 
  std::cout << "-----------------------------------------------------------" << std::endl;
 
  // Sort out the file tag by adding an underscore
  if(tag != "")
    tag = "_"+tag;

  //--------------------------------------------------------------------------------- ---------
  //                                    Initialise
  //--------------------------------------------------------------------------------- ---------

  // Setup TTree from input file list
  std::cout << " Reading files and filling tree..." << std::endl;

  EventProcessor evtProc(allowed, input_list, n);
  evtProc.Initialize();

  // Now setup the tree and event objects to work with
  TChain *tree = evtProc.GetTree();
  anatree *evt = evtProc.GetEvents();

  //Now the trees for the TMVA
  std::unique_ptr<TFile> mySignalFile( TFile::Open("signal.root", "RECREATE") );
  auto sigtree = std::make_unique<TTree>("sigtree", "Signal Tree");
  float trkpuritytree, trkthetaxztree, trkthetayztree, lengthtree, trkcomptree, trkstartdtree;
  float nvtxtree;  // from int
  float distEntertree, distExittree; //from double
  //sigtree->Branch("trkpurity", &trkpuritytree);
  sigtree->Branch("nvtx", &nvtxtree);
  sigtree->Branch("trkthetaxz", &trkthetaxztree);
  sigtree->Branch("trkthetayz", &trkthetayztree);
  sigtree->Branch("length", &lengthtree);
  sigtree->Branch("distEnter", &distEntertree);
  sigtree->Branch("distExit", &distExittree);
  //sigtree->Branch("completeness", &trkcomptree);
  sigtree->Branch("trkstartd", &trkstartdtree);

  std::unique_ptr<TFile> myBkgFile( TFile::Open("background.root", "RECREATE") );
  auto bkgtree = std::make_unique<TTree>("bkgtree", "Background Tree");
  //bkgtree->Branch("trkpurity", &trkpuritytree);
  bkgtree->Branch("nvtx", &nvtxtree);
  bkgtree->Branch("trkthetaxz", &trkthetaxztree);
  bkgtree->Branch("trkthetayz", &trkthetayztree);
  bkgtree->Branch("length", &lengthtree);
  bkgtree->Branch("distEnter", &distEntertree);
  bkgtree->Branch("distExit", &distExittree);
  bkgtree->Branch("trkstartd", &trkstartdtree);
  //bkgtree->Branch("completeness", &trkcomptree);
  
  // Start of analysis (loop over chain and events
  std::cout << " Running analysis..." << std::endl;

  // Then setup the histograms, counters and any other variables to add to
  // Setup histograms
  TH1D *h_muon_len_true   = new TH1D("h_muon_len_true","",100,0,2.2e3);   // Reconstructed length of true selected stopping muon signal
  TH1D *h_muon_len_genpop   = new TH1D("h_muon_len_genpop","",100,0,2.2e3);   // Reconstructed length of all events
  //TH1D *h_muon_resrng_true   = new TH1D("h_muon_resrng_true","",100,0,2.2e3);   // Reconstructed residuial range of true selected stopping muon signal
  //TH1D *h_muon_resrng_genpop   = new TH1D("h_muon_resrng_genpop","",100,0,2.2e3);   // Reconstructed residiucal range of all events
  //TH1D *h_muon_trkpur_true   = new TH1D("h_muon_trkpur_true","",20,0,1);   // Reconstructed track purity of true selected stopping muon signal
  //TH1D *h_muon_trkpur_genpop   = new TH1D("h_muon_trkpur_genpop","",20,0,1);   // Reconstructed track purity of all events
  TH1D *h_muon_nvtx_true   = new TH1D("h_muon_nvtx_true","",20,0,20);   // Reconstructed # vertex of true selected stopping muon signal
  TH1D *h_muon_nvtx_genpop   = new TH1D("h_muon_nvtx_genpop","",20,0,20);   // Reconstructed # vertex of all events
  TH1D *h_muon_txz_true   = new TH1D("h_muon_txz_true","",28,-3.5,3.5);   // Reconstructed thetaxz of true selected stopping muon signal
  TH1D *h_muon_txz_genpop   = new TH1D("h_muon_txz_genpop","",28,-3.5,3.5);   // Reconstructed thetaxz of all events
  TH1D *h_muon_tyz_true   = new TH1D("h_muon_tyz_true","",28,-3.5,3.5);   // Reconstructed thetayz of true selected stopping muon signal
  TH1D *h_muon_tyz_genpop   = new TH1D("h_muon_tyz_genpop","",28,-3.5,3.5);   // Reconstructed thetayz of all events
  //TH1D *h_muon_comp_true   = new TH1D("h_muon_comp_true","",20,0,1);   // Reconstructed trk completeness of true selected stopping muon signal
  //TH1D *h_muon_comp_genpop   = new TH1D("h_muon_comp_genpop","",20,0,1);   // Reconstructed trk completenss of all events
  TH1D *h_muon_hits_true   = new TH1D("h_muon_hits_true","",100,0,1000);   // Reconstructed trk hits of true selected stopping muon signal
  TH1D *h_muon_hits_genpop   = new TH1D("h_muon_hits_genpop","",100,0,1000);   // Reconstructed trk hits of all events
  TH1D *h_muon_startd_true   = new TH1D("h_muon_startd_true","",70,-100,600);   // Reconstructed trk startd of true selected stopping muon signal
  TH1D *h_muon_startd_genpop   = new TH1D("h_muon_startd_genpop","",70,-100,600);   // Reconstructed trk startd of all events
  TH1D *h_muon_mom_true   = new TH1D("h_muon_mom_true","",200,-50,50);   // Reconstructed trk mom of true selected stopping muon signal
  TH1D *h_muon_mom_genpop   = new TH1D("h_muon_mom_genpop","",200,-50,50);   // Reconstructed trk mom of all events
  
  // Setup counters
  unsigned int totalTracksTrue = 0;
  unsigned int trueSignalMuons = 0;
  unsigned int totalTracksReco = 0;
  unsigned int recoSelectedMuons = 0;
  unsigned int recoSelectedSignalMuons = 0;
  
  // Now loop over the events
  unsigned int nEvts = tree->GetEntries();
  unsigned int iIt = 1;
  unsigned int trackRepeats = 0;
  unsigned int eventNum = 0;

  //std::vector<int> recoTrkPassId; //vector of true track IDs that pass reco cuts
  std::cout << " |";
  for(unsigned int iEvt = 0; iEvt < nEvts; ++iEvt){
    tree->GetEntry(iEvt);
    if(!evtProc.SelectEvent(evt)) continue;
    
    // Get the total number of true and reconstructed tracks to loop over
    int nTrks = evt->ntracks_pandoraTrack;   //reco
    int nGeant = evt->geant_list_size;                //true
    //std::cout << "Reco tracks = " << nTrks << std::endl;
    //std::cout << "True tracks = " << nGeant << std::endl;
    
    // Print the processing rate
    double evtFrac  = iEvt/static_cast<double>(nEvts);

    // Prints out how much has been completed so far
    if(std::abs(0.1*iIt-evtFrac) < std::numeric_limits<double>::epsilon()){
      std::cout << " --- " << evtFrac*100 << " %";
      std::cout.flush();
      iIt++;
    }
  
    ///////////////////////////////////
    //          TRUTH                //
    ///////////////////////////////////          
    std::vector<int> trueTrkPassId; //vector of track IDs that pass true signal cuts
    // Now loop over the true tracks
    for(int iTrktru = 0; iTrktru < nGeant; ++iTrktru){

      // Count tracks
      totalTracksTrue++;

      // Look for true pdg
      int trupdg = evt->pdg[iTrktru];
      if (abs(trupdg) != 13)
        continue;

      // Look for true mother, to check is primary
      int truMother = evt->Mother[iTrktru];
      if (truMother != 0)
        continue;

      // Check the the true end coordinates are within the TPC active volume
      // The general start and end points (including cryostat and TPC)
      TVector3 start(evt->StartPointx[iTrktru],evt->StartPointy[iTrktru],evt->StartPointz[iTrktru]);
      TVector3 end(evt->EndPointx[iTrktru],evt->EndPointy[iTrktru],evt->EndPointz[iTrktru]);

      // The tpc AV start and end points
      // me/ph1afm/DUNE/cosmicpandora/muon_startd_genpop_stoppingmuons_v02.png' 
      TVector3 startAV(evt->StartPointx_tpcAV[iTrktru],evt->StartPointy_tpcAV[iTrktru],evt->StartPointz_tpcAV[iTrktru]);
      TVector3 endAV(evt->EndPointx_tpcAV[iTrktru],evt->EndPointy_tpcAV[iTrktru],evt->EndPointz_tpcAV[iTrktru]);

      // Get the differences between the two
      float dx = abs(endAV.X()-end.X())+abs(startAV.X()-start.X());
      float dy = abs(endAV.Y()-end.Y())+abs(startAV.Y()-start.Y());
      float dz = abs(endAV.Z()-end.Z())+abs(startAV.Z()-start.Z());

      // If they don't match, it doesn't stop (i.e. it left the TPC so it's end point will be one of the walls)
      if (dx+dy+dz > 1e-10)
        continue;

      //Fill truth plots (with the reco variables!)
      h_muon_len_true->Fill(evt->trklen_pandoraTrack[iTrktru]);
      //h_muon_resrng_true->Fill(evt->trkresrg_pandoraTrack[iTrktru][2]);
      //h_muon_trkpur_true->Fill(evt->trkpurity_pandoraTrack[iTrktru]);
      h_muon_nvtx_true->Fill(evt->nvtx_pandora);
      h_muon_txz_true->Fill(evt->trkthetaxz_pandoraTrack[iTrktru]);
      h_muon_tyz_true->Fill(evt->trkthetayz_pandoraTrack[iTrktru]);
      //h_muon_comp_true->Fill(evt->trkcompleteness_pandoraTrack[iTrktru]);
      //h_muon_hits_true->Fill(evt->ntrkhits_pandoraTrack[iTrktru]);
      h_muon_startd_true->Fill(evt->trkstartd_pandoraTrack[iTrktru]);
      h_muon_mom_true->Fill(evt->trkmom_pandoraTrack[iTrktru]);

      trueSignalMuons++;
      trueTrkPassId.push_back(evt->TrackId[iTrktru]);


    } // iTrktru, truth loop

    ///////////////////////////////////
    //            RECO               //
    ///////////////////////////////////
    std::vector<int> recoTrkPassId; //vector of true track IDs that pass reco cuts
    for(int iTrk = 0; iTrk < nTrks; ++iTrk){

      // Count tracks
      totalTracksReco++;

      // Get the track verticies points
      TVector3 startVtx(evt->trkstartx_pandoraTrack[iTrk],
                   evt->trkstarty_pandoraTrack[iTrk],
                   evt->trkstartz_pandoraTrack[iTrk]);
      TVector3 endVtx(evt->trkendx_pandoraTrack[iTrk],
                   evt->trkendy_pandoraTrack[iTrk],
                   evt->trkendz_pandoraTrack[iTrk]);

      CheckAndFlip(startVtx,endVtx);

      // Get the reconstructed best plane for this track (the one with most hits)
      int bestPlane = 0;
      std::vector<int> hitsOnPlane(3,0);
      GetRecoBestPlane(iTrk, evt, bestPlane, hitsOnPlane);
      
      //Check the track only crosses one external plane
      float length = evt->trklen_pandoraTrack[iTrk];       //cm

      Plane enteringPlane = GetClosestPlane(extPlanes, startVtx, endVtx);
      double distFromEntrance = GetDistanceToPlane(enteringPlane, startVtx, endVtx);

      Plane exitingPlane = GetClosestPlane(extPlanes, endVtx, startVtx);
      double distFromExit = GetDistanceToPlane(exitingPlane, endVtx, startVtx);

      unsigned int nExtCrossed    = 0;
      for(const Plane &pl : extPlanes){
        if(enteringPlane.GetLabel() == pl.GetLabel()){
          if(distFromEntrance < 1){
            nExtCrossed++;
          }
        } // Intersects
        else if(exitingPlane.GetLabel() == pl.GetLabel()){
          if(distFromExit < 1){
            nExtCrossed++;
          }
        } // Intersects
        else if(CheckIfIntersectsPlane(pl,startVtx,endVtx,length)){
          nExtCrossed++;
        } // Intersects
      } // Planes

      //if it crosses more then one external plane, it doesn't stop
      //and if it crosses less than one external plane it's not a primary cosmic muon
      //if (nExtCrossed != 1)
      //  continue;

      //Also need to ensure the track is not a fragment, so set a minimum length
      if (length < 15)
        continue;

      //Now apply angular conditions
      float thetaYZ = evt->trkthetayz_pandoraTrack[iTrk];

      if ((thetaYZ > -0.5))
        continue;

     //consider the number of reco verticies in the event
     int nvtx = evt->nvtx_pandora;
     if (nvtx > 25)
        continue;

     //float purity = evt->trkpurity_pandoraTrack[iTrk];
     //if (purity < 0.8)
     //     continue;
     
     //float completeness = evt->trkcompleteness_pandoraTrack[iTrk];
     //if (completeness < 0.45)
     //     continue;
     
     //consider the track's start direction
     float trkstartd = evt->trkstartd_pandoraTrack[iTrk];
     if ((trkstartd > 50))
       continue;

     h_muon_len_genpop->Fill(evt->trklen_pandoraTrack[iTrk]);
     //h_muon_resrng_genpop->Fill(evt->trkresrg_pandoraTrack[iTrk][2]);
     //h_muon_trkpur_genpop->Fill(evt->trkpurity_pandoraTrack[iTrk]);
     h_muon_nvtx_genpop->Fill(evt->nvtx_pandora);
     h_muon_txz_genpop->Fill(evt->trkthetaxz_pandoraTrack[iTrk]);
     h_muon_tyz_genpop->Fill(evt->trkthetayz_pandoraTrack[iTrk]);
     //h_muon_comp_genpop->Fill(evt->trkcompleteness_pandoraTrack[iTrk]);
     //h_muon_hits_genpop->Fill(evt->ntrkhits_pandoraTrack[iTrk]);
     h_muon_startd_genpop->Fill(evt->trkstartd_pandoraTrack[iTrk]);
     h_muon_mom_genpop->Fill(evt->trkmom_pandoraTrack[iTrk]);


     recoSelectedMuons++;

     //Now need to check how many of the selected tracks are also true signal
     int trueID = evt->trkidtruth_pandoraTrack[iTrk][bestPlane];
     recoTrkPassId.push_back(trueID);
     if(CheckTrueIDAssoc(trueID,trueTrkPassId)) {
       recoSelectedSignalMuons++; 
       //trkpuritytree = purity;
       nvtxtree = static_cast<float>(nvtx);
       trkthetayztree = thetaYZ;
       trkthetaxztree = evt->trkthetaxz_pandoraTrack[iTrk];
       lengthtree = length;
       distEntertree = static_cast<float>(distFromEntrance);
       distExittree = static_cast<float>(distFromExit);
       trkstartdtree = trkstartd;
       //trkcomptree = completeness;
       sigtree->Fill();
     } //if selected signal
     else {
       //trkpuritytree = purity;
       nvtxtree = static_cast<float>(nvtx);
       trkthetayztree = thetaYZ;
       trkthetaxztree = evt->trkthetaxz_pandoraTrack[iTrk];
       lengthtree = length;
       distEntertree = static_cast<float>(distFromEntrance);
       distExittree = static_cast<float>(distFromExit);
       trkstartdtree = trkstartd;
       //trkcomptree = completeness;
       bkgtree->Fill();
     } //else

     //For the TMVA to work, it will need a tree of true signal and a tree of true background
     //But these will be in reconstructed space, so it can learn how they look
     //The signal tree will be true signal events which pass the signal cuts (reco selected signal muons)
     //The background tree will be events which pass the precuts, but are not true signal events

    } // iTrk, reco loop

    if (recoTrkPassId.size() > 1) {
      std::sort(recoTrkPassId.begin(), recoTrkPassId.end());
      auto i1 = std::adjacent_find(recoTrkPassId.begin(), recoTrkPassId.end());
      bool isUnique = (i1 == recoTrkPassId.end());
      if (isUnique == 0) {
        trackRepeats++;
      }  //is unique
    } //recopassID
  eventNum++;
  }// Event loop


  std::cout << " --- 100 % --- |" << std::endl;
  //sigtree->Scan();
  //bkgtree->Scan();
  mySignalFile->Write();
  myBkgFile->Write();


  //Calculate the efficiency and purity of the selection
  float purity = 0;
  float efficiency = 100;
  if (recoSelectedMuons != 0)
    purity = ((float)recoSelectedSignalMuons/(float)recoSelectedMuons)*100;
  if (trueSignalMuons != 0)
    efficiency = ((float)recoSelectedSignalMuons/(float)trueSignalMuons)*100;

  // Print Stats
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Results..." << std::endl;
  std::cout << " True Track #           = " << totalTracksTrue << std::endl;
  std::cout << " True Signal #          = " << trueSignalMuons << std::endl;
  std::cout << " Reco Track #           = " << totalTracksReco << std::endl;
  std::cout << " Reco Selected #        = " << recoSelectedMuons << std::endl;
  std::cout << " Reco Selected Signal # = " << recoSelectedSignalMuons << std::endl;
  std::cout << " Selection Purity       = " << purity << std::endl;
  std::cout << " Selection Efficiency   = " << efficiency << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;

  // Now write the histograms
  TCanvas *c1 = new TCanvas("c1","",900,900);
  SetCanvasStyle(c1, 0.12,0.08,0.06,0.12,0,0,0);

  TLegend *l = new TLegend(0.22,0.94,0.98,0.995);
  l->SetNColumns(3);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->SetTextFont(132);

  //length
  SetHistogramStyle1D(h_muon_len_true,"Track length [cm]", "Rate");
  h_muon_len_true->Draw("hist");
  h_muon_len_true->SetLineWidth(3);
  h_muon_len_true->SetLineColor(kTeal-5);
  h_muon_len_true->GetYaxis()->SetTitleOffset(0.95);
  c1->SaveAs((location+"/muon_len_true"+tag+".png").c_str());
  c1->SaveAs((location+"/muon_len_true"+tag+".root").c_str());
  c1->Clear();

  SetHistogramStyle1D(h_muon_len_genpop,"Track length [cm]", "Rate");
  h_muon_len_genpop->Draw("hist");
  h_muon_len_genpop->SetLineWidth(3);
  h_muon_len_genpop->SetLineColor(kTeal-5);
  h_muon_len_genpop->GetYaxis()->SetTitleOffset(0.95);
  c1->SaveAs((location+"/muon_len_genpop"+tag+".png").c_str());
  c1->SaveAs((location+"/muon_len_genpop"+tag+".root").c_str());
  c1->Clear();
 
  //n vertex
  SetHistogramStyle1D(h_muon_nvtx_true,"N Vertex", "Rate");
  h_muon_nvtx_true->Draw("hist");
  h_muon_nvtx_true->SetLineWidth(3);
  h_muon_nvtx_true->SetLineColor(kTeal-5);
  h_muon_nvtx_true->GetYaxis()->SetTitleOffset(0.95);
  c1->SaveAs((location+"/muon_nvtx_true"+tag+".png").c_str());
  c1->SaveAs((location+"/muon_nvtx_true"+tag+".root").c_str());
  c1->Clear();

  SetHistogramStyle1D(h_muon_nvtx_genpop,"N Vertex", "Rate");
  h_muon_nvtx_genpop->Draw("hist");
  h_muon_nvtx_genpop->SetLineWidth(3);
  h_muon_nvtx_genpop->SetLineColor(kTeal-5);
  h_muon_nvtx_genpop->GetYaxis()->SetTitleOffset(0.95);
  c1->SaveAs((location+"/muon_nvtx_genpop"+tag+".png").c_str());
  c1->SaveAs((location+"/muon_nvtx_genpop"+tag+".root").c_str());
  c1->Clear();

  //Theta xz
  SetHistogramStyle1D(h_muon_txz_true,"Theta XZ", "Rate");
  h_muon_txz_true->Draw("hist");
  h_muon_txz_true->SetLineWidth(3);
  h_muon_txz_true->SetLineColor(kTeal-5);
  h_muon_txz_true->GetYaxis()->SetTitleOffset(0.95);
  c1->SaveAs((location+"/muon_txz_true"+tag+".png").c_str());
  c1->SaveAs((location+"/muon_txz_true"+tag+".root").c_str());
  c1->Clear();

  SetHistogramStyle1D(h_muon_txz_genpop,"Theta XZ", "Rate");
  h_muon_txz_genpop->Draw("hist");
  h_muon_txz_genpop->SetLineWidth(3);
  h_muon_txz_genpop->SetLineColor(kTeal-5);
  h_muon_txz_genpop->GetYaxis()->SetTitleOffset(0.95);
  c1->SaveAs((location+"/muon_txz_genpop"+tag+".png").c_str());
  c1->SaveAs((location+"/muon_txz_genpop"+tag+".root").c_str());
  c1->Clear();

  //Theta yz
  SetHistogramStyle1D(h_muon_tyz_true,"Theta YZ", "Rate");
  h_muon_tyz_true->Draw("hist");
  h_muon_tyz_true->SetLineWidth(3);
  h_muon_tyz_true->SetLineColor(kTeal-5);
  h_muon_tyz_true->GetYaxis()->SetTitleOffset(0.95);
  c1->SaveAs((location+"/muon_tyz_true"+tag+".png").c_str());
  c1->SaveAs((location+"/muon_tyz_true"+tag+".root").c_str());
  c1->Clear();

  SetHistogramStyle1D(h_muon_tyz_genpop,"Theta YZ", "Rate");
  h_muon_tyz_genpop->Draw("hist");
  h_muon_tyz_genpop->SetLineWidth(3);
  h_muon_tyz_genpop->SetLineColor(kTeal-5);
  h_muon_tyz_genpop->GetYaxis()->SetTitleOffset(0.95);
  c1->SaveAs((location+"/muon_tyz_genpop"+tag+".png").c_str());
  c1->SaveAs((location+"/muon_tyz_genpop"+tag+".root").c_str());
  c1->Clear();

  //ntrkhits
  SetHistogramStyle1D(h_muon_hits_true,"Hits", "Rate");
  h_muon_hits_true->Draw("hist");
  h_muon_hits_true->SetLineWidth(3);
  h_muon_hits_true->SetLineColor(kTeal-5);
  h_muon_hits_true->GetYaxis()->SetTitleOffset(0.95);
  c1->SaveAs((location+"/muon_hits_true"+tag+".png").c_str());
  c1->SaveAs((location+"/muon_hits_true"+tag+".root").c_str());
  c1->Clear();

  SetHistogramStyle1D(h_muon_hits_genpop,"Hits", "Rate");
  h_muon_hits_genpop->Draw("hist");
  h_muon_hits_genpop->SetLineWidth(3);
  h_muon_hits_genpop->SetLineColor(kTeal-5);
  h_muon_hits_genpop->GetYaxis()->SetTitleOffset(0.95);
  c1->SaveAs((location+"/muon_hits_genpop"+tag+".png").c_str());
  c1->SaveAs((location+"/muon_hits_genpop"+tag+".root").c_str());
  c1->Clear();

  //startd
  SetHistogramStyle1D(h_muon_startd_true,"Startd", "Rate");
  h_muon_startd_true->Draw("hist");
  h_muon_startd_true->SetLineWidth(3);
  h_muon_startd_true->SetLineColor(kTeal-5);
  h_muon_startd_true->GetYaxis()->SetTitleOffset(0.95);
  c1->SaveAs((location+"/muon_startd_true"+tag+".png").c_str());
  c1->SaveAs((location+"/muon_startd_true"+tag+".root").c_str());
  c1->Clear();

  SetHistogramStyle1D(h_muon_startd_genpop,"Startd", "Rate");
  h_muon_startd_genpop->Draw("hist");
  h_muon_startd_genpop->SetLineWidth(3);
  h_muon_startd_genpop->SetLineColor(kTeal-5);
  h_muon_startd_genpop->GetYaxis()->SetTitleOffset(0.95);
  c1->SaveAs((location+"/muon_startd_genpop"+tag+".png").c_str());
  c1->SaveAs((location+"/muon_startd_genpop"+tag+".root").c_str());
  c1->Clear();

  //mom
  SetHistogramStyle1D(h_muon_mom_true,"Mom", "Rate");
  h_muon_mom_true->Draw("hist");
  h_muon_mom_true->SetLineWidth(3);
  h_muon_mom_true->SetLineColor(kTeal-5);
  h_muon_mom_true->GetYaxis()->SetTitleOffset(0.95);
  c1->SaveAs((location+"/muon_mom_true"+tag+".png").c_str());
  c1->SaveAs((location+"/muon_mom_true"+tag+".root").c_str());
  c1->Clear();

  SetHistogramStyle1D(h_muon_mom_genpop,"Mom", "Rate");
  h_muon_mom_genpop->Draw("hist");
  h_muon_mom_genpop->SetLineWidth(3);
  h_muon_mom_genpop->SetLineColor(kTeal-5);
  h_muon_mom_genpop->GetYaxis()->SetTitleOffset(0.95);
  c1->SaveAs((location+"/muon_mom_genpop"+tag+".png").c_str());
  c1->SaveAs((location+"/muon_mom_genpop"+tag+".root").c_str());
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
