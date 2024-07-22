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
#include "TMVA/GeneticAlgorithm.h"
#include "TMVA/GeneticFitter.h"
#include "TMVA/IFitterTarget.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Reader.h"
#include "TMVA/ClassifierFactory.h"
#include "TMVA/RReader.hxx"

//using namespace TMVA::Experimental;
using namespace calib;
using namespace cppsecrets;
using namespace TMVA;

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
   "pdg",  //<---pdg for true tracks (so only need to input track id)
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
   "trkpidpdg_pandoraTrack",
   "trkcompleteness_pandoraTrack",
   "trkorig_pandoraTrack",
   "trkflashT0_pandoraTrack"
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
     
int stoppingMuonStudyBDTApplied(const char *config){

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
  
  // Start of analysis (loop over chain and events
  std::cout << " Running analysis..." << std::endl;

  // Then setup the histograms, counters and any other variables to add to
  // Setup histograms
  TH1D *h_true_background_pdg   = new TH1D("h_true_background_pdg","",100,0,100);   // Reconstructed selected background true pdg codes
  
  // Setup counters
  unsigned int totalTracksTrue = 0;
  unsigned int trueSignalMuons = 0;
  unsigned int totalTracksReco = 0;
  unsigned int recoSelectedMuons = 0;
  unsigned int recoSelectedSignalMuons = 0;
  unsigned int trueSecondaryMuons = 0;
  
  // Now loop over the events
  unsigned int nEvts = tree->GetEntries();
  unsigned int iIt = 1;
  unsigned int trackRepeats = 0;
  unsigned int eventNum = 0;
  unsigned int dups_tot = 0;

  std::cout << " |";
  for(unsigned int iEvt = 0; iEvt < nEvts; ++iEvt){
    tree->GetEntry(iEvt);
    if(!evtProc.SelectEvent(evt)) continue;
    
    // Get the total number of true and reconstructed tracks to loop over
    int nTrks = evt->ntracks_pandoraTrack;   //reco
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
    //std::cout << "Looping over true tracks..." << std::endl;
    for(int iTrktru = 0; iTrktru < nGeant; ++iTrktru){

      // Count tracks
      totalTracksTrue++;

      // Look for true pdg
      int trupdg = evt->pdg[iTrktru];
     // std::cout << "true pdg = " << trupdg << std::endl;
      if (abs(trupdg) != 13)
        continue;

      // Look for true mother, to check is primary
      int truMother = evt->Mother[iTrktru];
      if (truMother != 0) {
        trueSecondaryMuons++;
        continue;
      }

      // Check the the true end coordinates are within the TPC active volume
      // The general start and end points (including cryostat and TPC)
      TVector3 start(evt->StartPointx[iTrktru],evt->StartPointy[iTrktru],evt->StartPointz[iTrktru]);
      TVector3 end(evt->EndPointx[iTrktru],evt->EndPointy[iTrktru],evt->EndPointz[iTrktru]);

      // The tpc AV start and end points
      TVector3 startAV(evt->StartPointx_tpcAV[iTrktru],evt->StartPointy_tpcAV[iTrktru],evt->StartPointz_tpcAV[iTrktru]);
      TVector3 endAV(evt->EndPointx_tpcAV[iTrktru],evt->EndPointy_tpcAV[iTrktru],evt->EndPointz_tpcAV[iTrktru]);

      // Get the differences between the two
      float dx = abs(endAV.X()-end.X())+abs(startAV.X()-start.X());
      float dy = abs(endAV.Y()-end.Y())+abs(startAV.Y()-start.Y());
      float dz = abs(endAV.Z()-end.Z())+abs(startAV.Z()-start.Z());

      // If they don't match, it doesn't stop (i.e. it left the TPC so it's end point will be one of the walls)
      if (dx+dy+dz > 1e-10)
        continue;

      trueSignalMuons++;
      //std::cout << "evt->TrackId[iTrktru] = " << evt->TrackId[iTrktru] << std::endl;
      trueTrkPassId.push_back(evt->TrackId[iTrktru]);

      //std::sort(trueTrkPassId.begin(), trueTrkPassId.end());

      //for(auto it = std::cbegin(trueTrkPassId); it != std::cend(trueTrkPassId); ) {
        //std:;cout << "it " << *it <<std::endl;
        //int dups = std::count(it, std::cend(trueTrkPassId), *it);
        //if ( dups > 1 )
        //   cout << *it << " is a true duplicated number, times: " << dups << endl;
        //for(auto last = *it;*++it == last;);
     // }


    } // iTrktru, truth loop

    ///////////////////////////////////
    //            RECO               //
    ///////////////////////////////////
    std::vector<int> recoTrkPassId; //vector of true track IDs that pass reco cuts
    //std::cout << "Looping over reco tracks..." << std::endl;
    for(int iTrk = 0; iTrk < nTrks; ++iTrk){

      // Count tracks
      //totalTracksReco++;

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

      //for only true signal
      int trueID = evt->trkidtruth_pandoraTrack[iTrk][bestPlane];
      if (!CheckTrueIDAssoc(trueID,trueTrkPassId))
         continue;

      totalTracksReco++;
      
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
      if (nExtCrossed != 1)
        continue;

      //Also need to ensure the track is not a fragment, so set a minimum length
     // if (length < 50)
     //   continue;

     //Now apply angular conditions
     float thetaYZ = evt->trkthetayz_pandoraTrack[iTrk];

     //if ((thetaYZ < -2.5) || (thetaYZ > -0.5 && thetaYZ < 0.5) || ( thetaYZ > 2.5))
     //  continue;

     //consider the number of reco verticies in the event
     int nvtx = evt->nvtx_pandora;
     //if (nvtx > 7)
     //    continue;

     float purity = evt->trkpurity_pandoraTrack[iTrk];
     //if (purity < 0.8)
     //     continue;

     float completeness = evt->trkcompleteness_pandoraTrack[iTrk];
     //if (completeness < 0.45)
     //     continue;
     
     //apply the BDT here
     //Load in the model from the TMMA xml file
     TMVA::Experimental::RReader model("datasetBkg0/weights/TMVAMultiBkg0_BDTG.weights.xml");

     float thetaXZ = evt->trkthetaxz_pandoraTrack[iTrk];

     //Apply model
     auto prediction = model.Compute({purity, static_cast<float>(nvtx), thetaXZ, thetaYZ, length, static_cast<float>(distFromEntrance), static_cast<float>(distFromExit), completeness});
     //std::cout << "Single-event inference: " << prediction[0] << std::endl;;

     //if (prediction[0] < -0.09161)
	//continue;       

     recoSelectedMuons++;

     //Now need to check how many of the selected tracks are also true signal
     //int trueID = evt->trkidtruth_pandoraTrack[iTrk][bestPlane];
     recoTrkPassId.push_back(trueID);
     if(CheckTrueIDAssoc(trueID,trueTrkPassId)) {
       recoSelectedSignalMuons++;
     } //if selected signal
     else {
      //add to background pdg plot
      int truetrkpdg = evt->trkpdgtruth_pandoraTrack[iTrk][bestPlane];
      //std::cout << "abs true pdg of selected background = " << abs(truetrkpdg)  << std::endl;
      TVector3 vtxAV(evt->StartPointx_tpcAV[trueID],evt->StartPointy_tpcAV[trueID],evt->StartPointz_tpcAV[trueID]);
      TVector3 endAV(evt->EndPointx_tpcAV[trueID],evt->EndPointy_tpcAV[trueID],evt->EndPointz_tpcAV[trueID]);
      float lengthAV = (endAV-vtxAV).Mag();
      h_true_background_pdg->Fill(abs(truetrkpdg));
     }


    } // iTrk, reco loop

    if (recoTrkPassId.size() > 1) {
      std::sort(recoTrkPassId.begin(), recoTrkPassId.end());

      //for(auto it = std::cbegin(recoTrkPassId); it != std::cend(recoTrkPassId); ) {  

        //int dups = std::count(it, std::cend(recoTrkPassId), *it);
        //if ( dups > 1 )
        //   dups_tot = dups_tot+dups-1;   //-1 so that original not included
        //   cout << *it << " is a duplicated number, times: " << dups << endl;
        //for(auto last = *it;*++it == last;);
      //}

      auto i1 = std::adjacent_find(recoTrkPassId.begin(), recoTrkPassId.end());
      bool isUnique = (i1 == recoTrkPassId.end());
      if (isUnique == 0) {
       // std::cout << "Repeats in the list of true IDs from selected reco tracks!" << std::endl;
        trackRepeats++;
      }
      //else if (isUnique == 1) {
        //std::cout << "NO repeats in the list of true IDs from selected reco tracks!" << std::endl;
     // }
    } //if
  eventNum++;
  }// Event loop
  std::cout << " --- 100 % --- |" << std::endl;

  //std::cout << "number of tracks in event = " << recoTrkPassId.size() <<std::endl;

  //Calculate the efficiency and purity of the selection
  float purity = 0;
  float efficiency = 100;
  //float puritySignalOnly = 0;
  //float efficiencySignalOnly = 100;
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
  std::cout << " True 2ndary Muon #     = " << trueSecondaryMuons << std::endl;
  std::cout << trackRepeats << " events with track repeats of " << eventNum << " events total" << std::endl;
  std::cout << " Number of duplicate reco tracks present = " << dups_tot << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;


  // Now write the histograms
  TCanvas *c1 = new TCanvas("c1","",900,900);
  SetCanvasStyle(c1, 0.12,0.08,0.06,0.12,0,0,0);
  TLegend *l = new TLegend(0.22,0.94,0.98,0.995);
  l->SetNColumns(3);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->SetTextFont(132);

  SetHistogramStyle1D(h_true_background_pdg,"PDG Code", "Rate");
  h_true_background_pdg->Draw("hist");
  h_true_background_pdg->SetLineWidth(3);
  h_true_background_pdg->SetLineColor(kTeal-5);
  h_true_background_pdg->GetYaxis()->SetTitleOffset(0.95);
  c1->SaveAs((location+"/true_background_pdg"+tag+".png").c_str());
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
