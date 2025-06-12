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
   "trkstartd_pandoraTrack",
   "trkendx_pandoraTrack",
   "trkendy_pandoraTrack",
   "trkendz_pandoraTrack",
   "trklen_pandoraTrack",
   "pdg",
   "Mother",
   "EndPointx_tpcAV",
   "EndPointy_tpcAV",
   "EndPointz_tpcAV",
   "StartPointx_tpcAV",
   "StartPointy_tpcAV",
   "StartPointz_tpcAV",
   "EndPointx",
   "EndPointx_drifted",
   "EndPointy",
   "EndPointz",
   "StartPointx",
   "StartPointx_drifted",
   "StartPointy",
   "StartPointz",
   "TrackId",
   "trkthetaxz_pandoraTrack",
   "trkthetayz_pandoraTrack",
   "trkpurity_pandoraTrack",
   "nvtx_pandora",
   "trkpidpdg_pandoraTrack",
   "trkcompleteness_pandoraTrack",
   "trkorig_pandoraTrack",
   "trkflashT0_pandoraTrack",
   "trktrueT0_pandoraTrack",
   "pathlen",
   "evttime",
   "beamtime",
   "triggertime",
   "StartT",
   "EndT",
   "StartT_tpcAV",
   "EndT_tpcAV",
   "trkke_pandoraTrack",
   "trkmom_pandoraTrack",
   "P",
   "trkrange_pandoraTrack",
   "trkpitchc_pandoraTrack"
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
     
int stoppingMuonStudy10kt(const char *config){

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
  float trkpuritytree, trkthetaxztree, trkthetayztree, trkcomptree, trkstartdtree, lengthtree, ketree, rangetree;
  float nvtxtree;  // from int
  float distExittree, startytree, endytree, startxtree, endxtree, startztree, distEntertree, endztree; //from double

  sigtree->Branch("nvtx", &nvtxtree);
  sigtree->Branch("trkthetaxz", &trkthetaxztree);
  sigtree->Branch("trkthetayz", &trkthetayztree);
  sigtree->Branch("length", &lengthtree);
  sigtree->Branch("distEnter", &distEntertree);
  sigtree->Branch("distExit", &distExittree);
  sigtree->Branch("trkstartd", &trkstartdtree);
  sigtree->Branch("starty", &startytree);
  sigtree->Branch("endy", &endytree);
  sigtree->Branch("startx", &startxtree);
  sigtree->Branch("endx", &endxtree);
  sigtree->Branch("startz", &startztree);
  sigtree->Branch("endz", &endztree);
  sigtree->Branch("ke", &ketree);
  sigtree->Branch("range", &rangetree);

  std::unique_ptr<TFile> myBkgFile( TFile::Open("background.root", "RECREATE") );
  auto bkgtree = std::make_unique<TTree>("bkgtree", "Background Tree");

  bkgtree->Branch("nvtx", &nvtxtree);
  bkgtree->Branch("trkthetaxz", &trkthetaxztree);
  bkgtree->Branch("trkthetayz", &trkthetayztree);
  bkgtree->Branch("length", &lengthtree);
  bkgtree->Branch("distEnter", &distEntertree);
  bkgtree->Branch("distExit", &distExittree);
  bkgtree->Branch("trkstartd", &trkstartdtree);
  bkgtree->Branch("starty", &startytree);
  bkgtree->Branch("endy", &endytree);
  bkgtree->Branch("startx", &startxtree);
  bkgtree->Branch("endx", &endxtree);
  bkgtree->Branch("startz", &startztree);
  bkgtree->Branch("endz", &endztree);
  bkgtree->Branch("ke", &ketree);
  bkgtree->Branch("range", &rangetree);
  
  // Start of analysis (loop over chain and events
  std::cout << " Running analysis..." << std::endl;

  // Then setup the histograms, counters and any other variables to add to
  // Setup histograms if wanted
  
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

  std::cout << " |";
  for(unsigned int iEvt = 0; iEvt < nEvts; ++iEvt){
    tree->GetEntry(iEvt);
    if(!evtProc.SelectEvent(evt)) continue;
    
    // Get the total number of true and reconstructed tracks to loop over
    int nTrks = evt->ntracks_pandoraTrack;            //reco
    int nGeant = evt->geant_list_size;                //true
    
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
      TVector3 startAV(evt->StartPointx_tpcAV[iTrktru],evt->StartPointy_tpcAV[iTrktru],evt->StartPointz_tpcAV[iTrktru]);
      TVector3 endAV(evt->EndPointx_tpcAV[iTrktru],evt->EndPointy_tpcAV[iTrktru],evt->EndPointz_tpcAV[iTrktru]);

      // Get the differences between the two
      float dx = abs(endAV.X()-end.X());
      float dy = abs(endAV.Y()-end.Y());
      float dz = abs(endAV.Z()-end.Z());

      // If they don't match, it doesn't stop (i.e. it left the TPC so it's end point will be one of the walls)
      if (dx+dy+dz > 1e-10)
        continue;

      int bestPlane = 0;
      std::vector<int> hitsOnPlane(3,0);
      GetRecoBestPlane(iTrktru, evt, bestPlane, hitsOnPlane);

      //Fill truth plots (with the reco variables!) if wanted

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


      //Also need to ensure the track is not a fragment, so set a minimum length
      if (length < 50)
        continue;

      //Now apply angular conditions
      float thetaYZ = evt->trkthetayz_pandoraTrack[iTrk];
      if ((thetaYZ > 0.0))
        continue;

     //consider the number of reco verticies in the event
     int nvtx = evt->nvtx_pandora;
     if (nvtx > 20)
        continue;
     
     //consider the tracks start
     float trkstartd = evt->trkstartd_pandoraTrack[iTrk];
     if ((trkstartd > 20))
       continue;

     if ((startVtx.Y() < -50))
       continue;

     if ((endVtx.Y() < -500))
       continue;

     if (( endVtx.X() < -725 || endVtx.X() > 725 ))
       continue;

     if ( (endVtx.Z() < 100 || endVtx.Z() > 5600))
       continue;

      float ke = evt->trkke_pandoraTrack[iTrk][bestPlane];
      if ( ke < 150 ) {
         continue;
      }

      float range = evt->trkrange_pandoraTrack[iTrk][bestPlane];
      if ( range < 60 ) {
         continue;
      }

     //Fill general plots here if you want

     recoSelectedMuons++;

     //Now need to check how many of the selected tracks are also true signal
     int trueID = evt->trkidtruth_pandoraTrack[iTrk][bestPlane];
     recoTrkPassId.push_back(trueID);

     //Fill trees for BDT
     if(CheckTrueIDAssoc(trueID,trueTrkPassId)) {
       recoSelectedSignalMuons++; 
       nvtxtree = static_cast<float>(nvtx);
       trkthetayztree = thetaYZ;
       trkthetaxztree = evt->trkthetaxz_pandoraTrack[iTrk];
       lengthtree = length;
       distEntertree = static_cast<float>(distFromEntrance);
       distExittree = static_cast<float>(distFromExit);
       trkstartdtree = trkstartd;
       startytree = startVtx.Y();
       endytree = endVtx.Y();
       startxtree = startVtx.X();
       endxtree = endVtx.X();
       startztree = startVtx.Z();
       endztree = endVtx.Z();
       ketree = ke;
       rangetree = range;
       sigtree->Fill();
     } //if selected signal
     else {
       nvtxtree = static_cast<float>(nvtx);
       trkthetayztree = thetaYZ;
       trkthetaxztree = evt->trkthetaxz_pandoraTrack[iTrk];
       lengthtree = length;
       distEntertree = static_cast<float>(distFromEntrance);
       distExittree = static_cast<float>(distFromExit);
       trkstartdtree = trkstartd;
       startytree = startVtx.Y();
       endytree = endVtx.Y();
       startxtree = startVtx.X();
       endxtree = endVtx.X();
       startztree = startVtx.Z();
       endztree = endVtx.Z();
       ketree = ke;
       rangetree = range;
       bkgtree->Fill();
     } //else

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

  // Now write the histograms if wanted

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
