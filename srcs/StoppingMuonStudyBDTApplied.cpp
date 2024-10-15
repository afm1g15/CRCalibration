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
#include "TMinuit.h"
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include "TProfile.h"

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
   "trkstartd_pandoraTrack",
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

double betap = 0.212;      //(kV/cm)(g/cm^2)/MeV // taken from ArgoNeuT experiment
double Rho = 1.3936;       // g/cm^3 (liquid argon density at temp 87.596 K 18.0 psia)
double Wion = 23.6e-6;     // work function of argon // parameter from ArgoNeuT experiment at 0.481kV/cm
double alpha = 0.93;       // parameter from ArgoNeuT experiment at 0.481kV/cm
double Efield = 0.50;      // kV/cm protoDUNE electric filed

double mass_muon = 105.658; // [MeV]
TSpline3 * muon_sp_range_to_KE = Get_sp_range_KE(mass_muon);

Float_t Dedx(float dqdx, float Ef)
{
  return (exp(dqdx * (betap / (Rho * Ef) * Wion)) - alpha) / (betap / (Rho * Ef));
}
     
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

  TH2F *h_reco_dQdx_RR = new TH2F("h_reco_dQdx_RR", "; Residual range [cm];dQ/dx [ADC/cm]", 200, 0, 200, 100, 0, 1000);
  TH2D *h_RR_bin_MPVs   = new TH2D("h_RR_bin_MPVs","",40,0,200, 100, 0, 500);   // RR bin MP dqdx values
  TH2D *h_RR_bin_dEdx_MPVs   = new TH2D("h_RR_bin_dEdx_MPVs","",40,0,200, 600, 0, 6);   // RR bin MP dEdx values
  TH2D *h_pitch_vs_RR   = new TH2D("h_RR_vs_pitch","",200,0,200, 100, 0, 10);   // pitch vs RR
  TH2D *h_RR_bin_MPVs_Th   = new TH2D("h_RR_bin_MPVs_Th","",40,0,200, 600, 0, 6);   // RR bin MP dEdx values (theory)
  TH2D *h_RR_vs_Ratio   = new TH2D("h_RR_vs_Ratio","",40,0,200, 100, 0, 0.01);   // RR binned ratio
  TH2F *h_reco_dEdx_RR = new TH2F("h_reco_dEdx_RR", "; Residual range [cm];dE/dx [MeV/cm]", 200, 0, 200, 100, 0, 10);

  int nbin = 40;
  int binsize = 5;

  TH1D *dqdx[nbin];
  for (int i = 0; i < nbin; ++i)
  {
   //std::cout << "i = " << i << std::endl;
    if (i == 0)
      dqdx[i] = new TH1D(Form("dqdx_%d", i), "; dQ/dx [ADC/cm]; Number of entries", 100, 0.0, 2000);

    if (i != 0)
      dqdx[i] = new TH1D(Form("dqdx_%d", i), "; dQ/dx [ADC/cm]; Number of entries", 50, 0.0, 1000);

    dqdx[i]->SetLineColor(kBlack);
    dqdx[i]->Sumw2(); // also store errors
  }
 
  // Setup counters
  unsigned int totalTracksTrue = 0;
  unsigned int trueSignalMuons = 0;
  unsigned int totalTracksReco = 0;
  unsigned int recoSelectedMuons = 0;
  unsigned int recoSelectedSignalMuons = 0;
  unsigned int trueSecondaryMuons = 0;
  unsigned int failAtPlaneCross = 0;
  unsigned int failAtLength = 0;
  unsigned int failAtAngle = 0;
  unsigned int failAtVertex = 0;
  unsigned int failAtBoundDist = 0;
  unsigned int failAtBDT = 0;
  unsigned int secondPass = 0;
  unsigned int firstPass = 0;
  unsigned int recoTrkIdLength = 0;


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
    std::vector<int> recoTrkId;
    std::vector<int> alreadyPassed;
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

      //for only true signal - why?
      int trueID = evt->trkidtruth_pandoraTrack[iTrk][bestPlane];
      //if (!CheckTrueIDAssoc(trueID,trueTrkPassId))
      //   continue;

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
          if(distFromEntrance < 1){  //1
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
          //std::cout << "intersects plane" << std::endl;
        } // Intersects
      } // Planes

      //if it crosses more then one external plane, it doesn't stop
      //and if it crosses less than one external plane it's not a primary cosmic muon
      //if (nExtCrossed != 1) {
      //  if(CheckTrueIDAssoc(trueID,trueTrkPassId)) {
      //    if(!CheckTrueIDAssoc(trueID,recoTrkId)) {
      //      failAtPlaneCross++;
      //      recoTrkId.push_back(trueID);
      //      std::cout << "# planes crossed = " << nExtCrossed << std::endl;
      //    }
      //  }
      //  continue;
     // }

      //Also need to ensure the track is not a fragment, so set a minimum length
      if (length < 15) {
        if(CheckTrueIDAssoc(trueID,trueTrkPassId)) {
          if(!CheckTrueIDAssoc(trueID,recoTrkId)) {
            failAtLength++;
            std::cout << "event number (length) " << eventNum << " and reco trk id " << iTrk << std::endl;
            recoTrkId.push_back(trueID);
          }
        }
        continue;
      }

      //Now apply angular conditions
      float thetaYZ = evt->trkthetayz_pandoraTrack[iTrk];

      if ((thetaYZ > -0.5)) {
        if(CheckTrueIDAssoc(trueID,trueTrkPassId)) {
          if(!CheckTrueIDAssoc(trueID,recoTrkId)) {
            failAtAngle++;
            std::cout << "event number (angle) " << eventNum << " and reco trk id " << iTrk << std::endl;
            recoTrkId.push_back(trueID);
          }
        }
        continue;
      }

      //consider the number of reco verticies in the event
      int nvtx = evt->nvtx_pandora;
      if (nvtx > 25) {
        if(CheckTrueIDAssoc(trueID,trueTrkPassId)) {
          if(!CheckTrueIDAssoc(trueID,recoTrkId)) {
            failAtVertex++;
            std::cout << "event number (vtx) " << eventNum << " and reco trk id " << iTrk << std::endl;
            recoTrkId.push_back(trueID);
          }
        }
         continue;
      }

      float trkstartd = evt->trkstartd_pandoraTrack[iTrk];
      if (trkstartd > 50) {
        if(CheckTrueIDAssoc(trueID,trueTrkPassId)) {
          if(!CheckTrueIDAssoc(trueID,recoTrkId)) {
            failAtBoundDist++;
            std::cout << "event number (trkstartd) " << eventNum << " and reco trk id " << iTrk << std::endl;
            recoTrkId.push_back(trueID);
          }
        }
         continue;
      }

      //float purity = evt->trkpurity_pandoraTrack[iTrk];
      //if (purity < 0.8)
      //    continue;

      //float completeness = evt->trkcompleteness_pandoraTrack[iTrk];
      //if (completeness < 0.45)
      //    continue;
     
      //apply the BDT here
      //Load in the model from the TMMA xml file
      TMVA::Experimental::RReader model("datasetBkg0/weights/TMVAMultiBkg0_BDTG.weights.xml");

      float thetaXZ = evt->trkthetaxz_pandoraTrack[iTrk];

      //Apply model
      auto prediction = model.Compute({static_cast<float>(nvtx), thetaXZ, thetaYZ, length, static_cast<float>(distFromEntrance), static_cast<float>(distFromExit), trkstartd});
      //std::cout << "Single-event inference: " << prediction[0] << std::endl;;

      if (prediction[0] < -0.0111086) {
        if(CheckTrueIDAssoc(trueID,trueTrkPassId)) {
          if(!CheckTrueIDAssoc(trueID,recoTrkId)) {
            failAtBDT++;
            std::cout << "event number (BDT) " << eventNum << " and reco trk id " << iTrk << std::endl;
            recoTrkId.push_back(trueID);
          }
        }
       // std::cout << "Single-event inference (signal failure): " << prediction[0] << std::endl;
	continue;       
      }

      recoSelectedMuons++;

      //Now need to check how many of the selected tracks are also true signal
      //int trueID = evt->trkidtruth_pandoraTrack[iTrk][bestPlane];
      recoTrkPassId.push_back(trueID);
      if(CheckTrueIDAssoc(trueID,recoTrkId) && CheckTrueIDAssoc(trueID,trueTrkPassId)) {
        int splitnum = std::count(recoTrkId.begin(), recoTrkId.end(), trueID);
        //std::cout << " track split passed at number " << splitnum  << std::endl;
        secondPass++;
        if (!CheckTrueIDAssoc(trueID,alreadyPassed)) {
          alreadyPassed.push_back(trueID);
        }
        else {
          dups_tot++;
        }
      }
      else if (!CheckTrueIDAssoc(trueID,recoTrkId) && CheckTrueIDAssoc(trueID,trueTrkPassId)) {
        //std::cout << "passed first time" << std::endl;
        firstPass++;
        recoTrkId.push_back(trueID);
        alreadyPassed.push_back(trueID);        
      }

      if(CheckTrueIDAssoc(trueID,trueTrkPassId)) {
        recoSelectedSignalMuons++;

        //Now take these selected tracks and use them in the calibration part
        int nHitsR = evt->ntrkhits_pandoraTrack[iTrk][bestPlane];
        for (int iHit = 1; iHit < nHitsR - 1; ++iHit)
        {
          if(evt->hit_plane[iHit] != bestPlane) continue;

          //Get the location of the current hit to determine the pitch
          TVector3 trkXYZ(evt->trkxyz_pandoraTrack[iTrk][bestPlane][iHit][0],
                          evt->trkxyz_pandoraTrack[iTrk][bestPlane][iHit][1],
                          evt->trkxyz_pandoraTrack[iTrk][bestPlane][iHit][2]);
          TVector3 nextXYZ(evt->trkxyz_pandoraTrack[iTrk][bestPlane][iHit+1][0],
                           evt->trkxyz_pandoraTrack[iTrk][bestPlane][iHit+1][1],
                           evt->trkxyz_pandoraTrack[iTrk][bestPlane][iHit+1][2]);
          
          //get the x value of the hit and convert it to a time
          float x = trkXYZ.X();
          float t = x * evtProc.kXtoT;
          double dp  = GetHitPitch(bestPlane, trkXYZ, nextXYZ);

          // Check if x is lower than the APA bound, charge seems to accumulate there
          if(x < evtProc.APA_X_POSITIONS[0] || x > evtProc.APA_X_POSITIONS[2]) continue;
          
          //Which tpc is the hit in
          int tpc = evtProc.WhichTPC(x) + 1;

          //make lifetime correction
          float dx = (-1 + 2 * (tpc % 2)) * (x - evtProc.APA_X_POSITIONS[tpc / 2]); // calculate the distance(positive) between hitx and apa plane
          float dt = dx * evtProc.kXtoT;
          float corr = TMath::Exp(-dt / 2.88); //2.88 appears to be 2.88 ms corrected lifetime value from through going muons

          //want to use the absolute energy scale, not bethe-bloch
          float corrected_dq_dx = evt->trkdqdx_pandoraTrack[iTrk][bestPlane][iHit] / corr;

          int bin = int(evt->trkresrg_pandoraTrack[iTrk][bestPlane][iHit]) / binsize;
         // if (bin < 2) {
         //   std::cout << "RR = " << int(evt->trkresrg_pandoraTrack[iTrk][bestPlane][iHit]) << " and binsize = " << binsize << " then bin = " << bin << std::endl;
         // }

          double hit_RR = evt->trkresrg_pandoraTrack[iTrk][bestPlane][iHit];

          h_reco_dQdx_RR->Fill(hit_RR, corrected_dq_dx);
          h_pitch_vs_RR->Fill(hit_RR, dp);

          double dEdx_corr = (0.0044741 + (0.00667017*(1/hit_RR)) + (2.19884e-06*hit_RR))*corrected_dq_dx;
          h_reco_dEdx_RR->Fill(hit_RR, dEdx_corr);
          
          if (bin < nbin)
          {
             dqdx[bin]->Fill(corrected_dq_dx);
             if (bin < 2) {
               //std::cout << "Filling bin " << bin << " with " << corrected_dq_dx << std::endl;
               //std::cout << "   Entries: " << dqdx[bin]->GetEntries() << std::endl;            
             }
          } 

        } //iHit

      } //if selected signal
      else {
       //add to background pdg plot
       int truetrkpdg = evt->trkpdgtruth_pandoraTrack[iTrk][bestPlane];
       TVector3 vtxAV(evt->StartPointx_tpcAV[trueID],evt->StartPointy_tpcAV[trueID],evt->StartPointz_tpcAV[trueID]);
       TVector3 endAV(evt->EndPointx_tpcAV[trueID],evt->EndPointy_tpcAV[trueID],evt->EndPointz_tpcAV[trueID]);
       float lengthAV = (endAV-vtxAV).Mag();
       h_true_background_pdg->Fill(abs(truetrkpdg));
      // if(CheckTrueIDAssoc(trueID,trueTrkPassId)) {
       // Have a look at true signal that ends up as background

      // }
      }


    } // iTrk, reco loop
    recoTrkIdLength = recoTrkIdLength + recoTrkId.size();

    if (recoTrkPassId.size() > 1) {
      std::sort(recoTrkPassId.begin(), recoTrkPassId.end());

      auto i1 = std::adjacent_find(recoTrkPassId.begin(), recoTrkPassId.end());
      bool isUnique = (i1 == recoTrkPassId.end());
      if (isUnique == 0) {
        trackRepeats++;
      }
      //else if (isUnique == 1) {
        //std::cout << "NO repeats in the list of true IDs from selected reco tracks!" << std::endl;
     // }
    } //if
  eventNum++;
  }// Event loop
  std::cout << " --- 100 % --- |" << std::endl;


  //Calculate the efficiency and purity of the selection
  float purity = 0;
  float efficiency = 100;
  float dup_adj_eff = 100;
  if (recoSelectedMuons != 0)
    purity = ((float)recoSelectedSignalMuons/(float)recoSelectedMuons)*100;
  if (trueSignalMuons != 0)
    efficiency = ((float)recoSelectedSignalMuons/(float)trueSignalMuons)*100;   
    dup_adj_eff = (((float)recoSelectedSignalMuons - (float)dups_tot)/(float)trueSignalMuons)*100;

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
  std::cout << " Selection Efficiency (adj for dups)   = " << dup_adj_eff << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " True 2ndary Muon #     = " << trueSecondaryMuons << std::endl;
  std::cout << " " << trackRepeats << " events with reco track repeats passing selection of " << eventNum << " events total" << std::endl;
  std::cout << " Number of duplicate reco tracks (of individual true tracks) selected = " << dups_tot << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Signal Failures At # Planes Crossed  = " << failAtPlaneCross << std::endl;
  std::cout << " Signal Failures At Length            = " << failAtLength << std::endl;
  std::cout << " Signal Failures At Angle             = " << failAtAngle << std::endl;
  std::cout << " Signal Failures At # Vertex in event = " << failAtVertex << std::endl;
  std::cout << " Signal Failures At Boundry Distance  = " << failAtBoundDist << std::endl;
  std::cout << " Signal Failures At BDT               = " << failAtBDT << std::endl;
  std::cout << "  " << std::endl;
  std::cout << " Signal Pass at 2nd (or more) split   = " << secondPass << std::endl;
  std::cout << " Signal Pass first time               = " << firstPass << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Signal Failures Sum                  = " << failAtPlaneCross + failAtLength + failAtAngle + failAtVertex + failAtBoundDist + failAtBDT  << std::endl;
  std::cout << " Signal TrueIds seen in Reco          = " << recoTrkIdLength << std::endl;

  //Now move on to the fitting
  std::cout << "Creating outfiles..." << std::endl;

  //for (int i = 0; i < nbin; i++)
 // {
 //   std::cout << "In bin " << i << "  Entries: " << dqdx[i]->GetEntries() << std::endl;
 // }

  ofstream myfile1;

  std::vector<double> mostProbValues;
  mostProbValues.resize(40, 0.0);
  std::cout << "Number of bins = " << nbin << std::endl;
  //for (int i = 0; i < nbin; i++)
 // {
 //   std::cout << "In bin " << i << "  Entries: " << dqdx[i]->GetEntries() << std::endl;
 // }


  for (int i = 0; i < nbin; i++)
  {

    std::cout << "Fitting bin ************************************ " << i << std::endl;
    std::cout << "In bin " << i << "  Entries: " << dqdx[i]->GetEntries() << std::endl;

    //set up a plot of the dqdx values in this bin
    TCanvas *d[i];
    d[i] = new TCanvas(Form("d_%d", i), Form("d_%d", i));

    //set up variables for fitting
    Double_t fp[4], fpe[4];
    Double_t chisqr;
    Int_t ndf;
    Int_t i2;
    Char_t FunName[100];

    sprintf(FunName, "Fitfcn_%s", dqdx[i]->GetName());
    TF1 *ffitold = (TF1 *)gROOT->GetListOfFunctions()->FindObject(FunName);
    if (ffitold)
      delete ffitold;

    //calculate some variables for fitting
    double norm = dqdx[i]->GetEntries() * dqdx[i]->GetBinWidth(1);
    int maxbin          = dqdx[i]->GetMaximumBin();
    double maxloc       = dqdx[i]->GetBinCenter(maxbin);
    std::cout << "-----------------------------" << std::endl;
    std::cout << "Starting parameter values..." << std::endl;
    std::cout << "   Landau scale: " << 0.2*maxloc << std::endl;
    std::cout << "   Landau MPV: " << maxloc << std::endl;
    std::cout << "   Norm: " << norm << std::endl;
    std::cout << "   Gauss Sigma: " << 0.2*maxloc << std::endl;
    double sv[4] = {0.2*maxloc, maxloc, norm, 0.2*maxloc}; // starting values for parameters: Landau scale, Landau MPV, Norm, Gauss sigma
    double minR = 0;
    double maxR = 0;
    if (i == 0) {
      minR = 500; //dedx[i]->GetXaxis()->GetXmin();
      maxR = 1500; //dedx[i]->GetXaxis()->GetXmax();
    }
    else {
      minR = 200; //dedx[i]->GetXaxis()->GetXmin();
      maxR = 800; //dedx[i]->GetXaxis()->GetXmax();
    }
    //double nBinsFromPeak = 30;
    //if(maxbin-nBinsFromPeak > 1)
    //  minR = dqdx[i]->GetBinCenter(maxbin-nBinsFromPeak);
    //if(maxbin+nBinsFromPeak < dqdx[i]->GetNbinsX())
    //  maxR = dqdx[i]->GetBinCenter(maxbin+nBinsFromPeak);
    std::cout << "   Min range: " << minR << std::endl;
    std::cout << "   Max range: " << maxR << std::endl;
    std::cout << "   Entries: " << dqdx[i]->GetEntries() << std::endl;
    std::cout << "-----------------------------" << std::endl;

    //if there's no entries in the bin, don't fit it
    if (dqdx[i]->GetEntries() == 0)
       continue;

    std::cout << "Entries present." << std::endl;

    //set up a new fit
    TF1 *ffit;
    ffit = new TF1(FunName, langaufun, minR, maxR, 4);
    ffit->SetParameters(sv);
    ffit->SetParNames("Width", "MPV", "Area", "GSigma"); 

    //Calculate the fit for this bin
    auto result = dqdx[i]->Fit(ffit, "QSMR", "");
    dqdx[i]->Print();

    ffit->GetParameters(fp); // obtain fit parameters
    for (i2 = 0; i2 < 4; i2++)
    {
      fpe[i2] = ffit->GetParError(i2); // obtain fit parameter errors
    }
    chisqr = ffit->GetChisquare(); // obtain chi^2
    ndf = ffit->GetNDF();          // obt
    std::cout << "chisqr = " << chisqr << std::endl;
    std::cout << "ndf = " << ndf << std::endl;

    double mpv = result->Parameter(1);

    std::cout << "------------------------------------" << std::endl;
    for(unsigned int p = 0; p < result->NPar(); ++p){
      std::cout << " " << result->ParName(p) << " : " << result->Parameter(p) << std::endl;
    }

    //Put the MPVs in plots
    std::cout << " Peak (MPV): " << mpv << std::endl;
    mostProbValues[i] = mpv;
    double range_for_bin = ((i+1)*binsize)-(binsize/2.0);
    double dEdx_corr_MPV = (0.0044741 + (0.00667017*(1/range_for_bin)) + (2.19884e-06*range_for_bin))*mpv;

    h_RR_bin_MPVs->Fill(range_for_bin, mpv);
    h_RR_bin_dEdx_MPVs->Fill(range_for_bin, dEdx_corr_MPV, 4);
    std::cout << "------------------------------------" << std::endl;

    //----------------------------------------------------------------------------------- 

    //------------------------------------------------------------------------------------
    //add the fit to the plots as a red line
    ffit->SetLineColor(kRed);

    dqdx[i]->SetStats(1);
    dqdx[i]->Draw("hist");
    ffit->Draw("same");

    gStyle->SetOptFit(1111);
    TPaveStats *ps = (TPaveStats *)dqdx[i]->FindObject("stats");
    ps->SetX1NDC(0.6);
    ps->SetX2NDC(0.85);
    ps->SetY1NDC(0.65);
    ps->SetY2NDC(0.9);

    d[i]->Write();
    d[i]->SaveAs(Form("d_%d.pdf", i));
    d[i]->Close();

    //If a fit passes the conditions, then this adds the range and energy measurements to the plots
    std::cout << "Checking results of fit..." << std::endl;
    if (gMinuit && dqdx[i]->GetEntries() > 100)
    {
      std::cout << "More than 100 entries in bin..." << std::endl;
      TString test = gMinuit->fCstatu.Data();
      if (ffit->GetNDF() != 0)
      {
        std::cout << "More than 0 DOF in bin..." << std::endl;
        std::cout << "ffit->GetParError(1) " << ffit->GetParError(1) << std::endl;
        std::cout << "ffit->GetChisquare() / ffit->GetNDF() " << ffit->GetChisquare() / ffit->GetNDF() << std::endl;
        std::cout << "test " << test << std::endl;
        if ((test.EqualTo("CONVERGED ") or test.EqualTo("OK ") ) && (ffit->GetParError(1) < 1000) && ((ffit->GetChisquare() / ffit->GetNDF() < 10)))
        {
         std::cout << "Adding results of this bin to main plots..." << std::endl;
          //range.push_back(r1);
          //erange.push_back(0);
          //range_measured[i] = r1;
          //energy_measured[i] = ffit->GetParameter(1);
          //energy.push_back(ffit->GetParameter(1));
          //eenergy.push_back(ffit->GetParError(1));
	  std::cout << "energy measured in bin " <<ffit->GetParameter(1) <<std::endl;
          //std::cout << "range measured in bin " << r1 <<std::endl;
          std::cout << "------------------------------------------" <<std::endl;
        }
      }
    }
  }


  //-------------------------------------------------------------------
  //calculate theoretical values
  
  //std::vector<double> theory_RR = [];
  //std::vector<double> theory_MPV_dEdx = [];
  std::cout << "" << std::endl;
  std::cout << "Number of MPVs = " << mostProbValues.size() << std::endl;
  for (int i = 0; i < nbin; i++)
  {
    //calculate theoretical KE for RR bin
    double range_for_bin = ((i+1)*binsize)-(binsize/2.0);
    double this_KE= muon_sp_range_to_KE -> Eval(range_for_bin); // == from rr

    std::cout << "Bin center RR = " << range_for_bin << std::endl;
    std::cout << "KE for this bin (from RR) = " << this_KE << std::endl;
    //std::cout << "-----------------------------" << std::endl;

    double pitch = 0.55; //Average per bin, compared to Rhiannon
    double this_MPV_dEdx = MPVdEdx(this_KE, pitch, mass_muon);
    std::cout << "dEdx MPV for this bin = " << this_MPV_dEdx << std::endl;
    // std::cout << "-----------------------------" << std::endl;
  
    //theory_RR.push_back(range_for_bin);
    //theory_MPV_dEdx.push_back(this_MPV_dEdx);

    h_RR_bin_MPVs_Th->Fill(range_for_bin, this_MPV_dEdx, 4);

    double ratio;
    if (mostProbValues[i] != 0.0)
    {    
      ratio = this_MPV_dEdx/mostProbValues[i];
      std::cout << "ratio for this ^ bin = " << ratio << std::endl;
      //std::cout << "-----------------------------" << std::endl;
    }
    else {
      ratio = 0.0;
      std::cout << "ratio for this ^ bin (failed) = " << ratio << std::endl;
      //std::cout << "-----------------------------" << std::endl;
    }

    h_RR_vs_Ratio->Fill(range_for_bin, ratio);
    std::cout << "filled ratio bin i = " << i << std::endl;
    std::cout << "-----------------------------" << std::endl;
  }
  //-------------------------------------------------------------------

  double sum = 0;
  std::cout << "Make some plots..." << std::endl;

  //-----------------------------------------------

  std::cout << "************************* Calibration.C has ended ***************************" << std::endl;


  //--------------------------------------------
  TCanvas *c3 = new TCanvas();
  h_reco_dQdx_RR->Draw("COLZ");
  h_reco_dQdx_RR->Write(" h_reco_dQdx_RR");
  c3->SaveAs("reco_dqdx_rr.png");
  c3->Write("reco_dQdx_RR");
  c3->Clear();

  //--------------------------------------------
  h_reco_dEdx_RR->Draw("COLZ");
  h_reco_dEdx_RR->Write(" h_reco_dEdx_RR");
  c3->SaveAs("reco_dEdx_rr.png");
  c3->Write("reco_dEdx_RR");
  c3->Clear();
  //------------------------------------------------

  c3->Close();

  gStyle->SetPalette(1, 0);
  gStyle->SetNumberContours(64);


  // Now write the histograms
  TCanvas *ca = new TCanvas("ca","",900,900);
  SetCanvasStyle(ca, 0.12,0.08,0.06,0.12,0,0,0);
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
  ca->SaveAs((location+"/true_background_pdg"+tag+".png").c_str());
  ca->Clear();

  //
  SetHistogramStyle2D(h_RR_bin_MPVs,"RR bin (5cm)", "dQdx MPV");
  h_RR_bin_MPVs->Draw("hist");
  h_RR_bin_MPVs->SetLineWidth(3);
  h_RR_bin_MPVs->SetLineColor(kTeal-5);
  h_RR_bin_MPVs->SetMarkerStyle(2);
  h_RR_bin_MPVs->SetMarkerSize(2);
  h_RR_bin_MPVs->GetYaxis()->SetTitleOffset(0.95);
  ca->SaveAs((location+"/RR_MPV_in_bin"+tag+".png").c_str());
  ca->Clear();

  SetHistogramStyle2D(h_RR_bin_dEdx_MPVs,"RR bin (5cm)", "dEdx MPV");
  h_RR_bin_dEdx_MPVs->Draw("box");
  h_RR_bin_dEdx_MPVs->SetLineWidth(3);
  h_RR_bin_dEdx_MPVs->SetLineColor(kTeal-5);
  h_RR_bin_dEdx_MPVs->SetMarkerStyle(2);
  h_RR_bin_dEdx_MPVs->SetMarkerSize(2);
  h_RR_bin_dEdx_MPVs->GetYaxis()->SetTitleOffset(0.95);
  ca->SaveAs((location+"/RR_MPV_dEdx_in_bin"+tag+".png").c_str());
  ca->Clear();

  //
  //SetHistogramStyle2D(h_pitch_vs_RR,"RR", "Pitch");
  //h_pitch_vs_RR->Draw("colz");
  //h_pitch_vs_RR->SetLineWidth(3);
  //h_pitch_vs_RR->SetLineColor(kTeal-5);
  //h_pitch_vs_RR->SetMarkerSize(5);
  //h_pitch_vs_RR->GetYaxis()->SetTitleOffset(0.95);
  TCanvas *cb = new TCanvas();
  h_pitch_vs_RR->Draw("COLZ");
  cb->SaveAs((location+"/RR_vs_pitch"+tag+".png").c_str());
  cb->Clear();

  //
  SetHistogramStyle2D(h_RR_bin_MPVs_Th,"RR bin (5cm)", "dEdx MPV (Theory)");
  h_RR_bin_MPVs_Th->Draw("box");
  h_RR_bin_MPVs_Th->SetLineWidth(3);
  h_RR_bin_MPVs_Th->SetLineColor(kTeal-5);
  h_RR_bin_MPVs_Th->SetMarkerSize(5);
  h_RR_vs_Ratio->SetMarkerStyle(2);
  h_RR_bin_MPVs_Th->GetYaxis()->SetTitleOffset(0.95);
  ca->SaveAs((location+"/RR_MPV_dEdx_Th_in_bin"+tag+".png").c_str());
  ca->Clear();

  //
  //TCanvas *c1 = new TCanvas("c1","show profile",600,900);
  //c1->Divide(1,2);
  //c1->cd(1);
  //TF1 *f1 = new TF1("f1", "expo", 5, 200);
  TF1 *f1 = new TF1("f1", "[0] + ([1]*(1/x)) + [2]*x", 0, 200);
  f1->SetParameters(0.,200.);
  f1->SetLineColor(kRed);
  h_RR_vs_Ratio->Fit(f1, "MR");

  SetHistogramStyle2D(h_RR_vs_Ratio,"RR bin (5cm)", "dEdx/dqdx");
  h_RR_vs_Ratio->Draw("box");
  h_RR_vs_Ratio->SetLineWidth(3);
  h_RR_vs_Ratio->SetLineColor(kTeal-5);
  h_RR_vs_Ratio->SetMarkerStyle(2);
  h_RR_vs_Ratio->SetMarkerSize(2);
  h_RR_vs_Ratio->GetYaxis()->SetTitleOffset(0.95);
  //c1->cd(2);
  std::cout << "------------------------------------------Fit--" <<std::endl;
  //TProfile *prof = h_RR_vs_Ratio->ProfileX();
  //prof->Fit("expo", "WW");
  //prof->Print();
  f1->Draw("same");

  double par0 = f1->GetParameter(0);
  double par1 = f1->GetParameter(1);
  double par2 = f1->GetParameter(2);

  std::cout << "Found parameter values..." << std::endl;
  std::cout << "    Par 0 = " << par0 << std::endl;
  std::cout << "    Par 1 = " << par1 << std::endl;
  std::cout << "    Par 2 = " << par2 << std::endl;

  std::cout << " " << std::endl;

  ca->SaveAs((location+"/RR_vs_Ratio"+tag+".png").c_str());
  ca->Clear();

  //------
  TCanvas *c4 = new TCanvas();
  h_reco_dEdx_RR->Draw("COLZ");
  h_RR_bin_MPVs_Th->Draw("box same");
  h_RR_bin_MPVs_Th->SetMarkerStyle(4);
  h_RR_bin_MPVs_Th->SetMarkerSize(2);
  h_RR_bin_MPVs_Th->SetMarkerColor(1);
  h_RR_bin_MPVs_Th->SetFillColor(1);
  h_RR_bin_dEdx_MPVs->Draw("box same");
  h_RR_bin_dEdx_MPVs->SetMarkerStyle(4);
  h_RR_bin_dEdx_MPVs->SetMarkerSize(2);
  h_RR_bin_dEdx_MPVs->SetMarkerColor(7);
  h_RR_bin_dEdx_MPVs->SetFillColor(7);

  TLegend *legend = new TLegend(0.1,0.7,0.48,0.9);
  legend->AddEntry(h_RR_bin_MPVs_Th,"Theory","l");
  legend->AddEntry(h_RR_bin_dEdx_MPVs,"Reco","l");

  c4->SaveAs((location+"/RR_dEdx_Th_Overlay"+tag+".png").c_str());
  c4-> Clear();

  h_RR_bin_MPVs_Th->Draw("box");
  h_RR_bin_MPVs_Th->SetMarkerStyle(4);
  h_RR_bin_MPVs_Th->SetMarkerSize(2);
  h_RR_bin_MPVs_Th->SetMarkerColor(1);
  h_RR_bin_MPVs_Th->SetFillColor(1);     //black
  h_RR_bin_dEdx_MPVs->Draw("box same");
  h_RR_bin_dEdx_MPVs->SetMarkerStyle(4);
  h_RR_bin_dEdx_MPVs->SetMarkerSize(2);
  h_RR_bin_dEdx_MPVs->SetMarkerColor(7);
  h_RR_bin_dEdx_MPVs->SetFillColor(7);   //light blue

  TLegend *legend2 = new TLegend(0.1,0.7,0.48,0.9);
  legend2->AddEntry(h_RR_bin_MPVs_Th,"Theory","l");
  legend2->AddEntry(h_RR_bin_dEdx_MPVs,"Reco","l");
  //legend2->Draw("same");

  c4->SaveAs((location+"/RR_dEdx_Th_Minus_Overlay"+tag+".png").c_str());
  c4-> Clear();

  h_RR_bin_dEdx_MPVs->Draw("box");
  h_RR_bin_dEdx_MPVs->SetMarkerStyle(4);
  h_RR_bin_dEdx_MPVs->SetMarkerSize(2);
  h_RR_bin_dEdx_MPVs->SetMarkerColor(7);
  h_RR_bin_dEdx_MPVs->SetFillColor(7);   //light blue
  h_RR_bin_MPVs_Th->Draw("box same");
  h_RR_bin_MPVs_Th->SetMarkerStyle(4);
  h_RR_bin_MPVs_Th->SetMarkerSize(2);
  h_RR_bin_MPVs_Th->SetMarkerColor(1);
  h_RR_bin_MPVs_Th->SetFillColor(1);     //black

  c4->SaveAs((location+"/RR_dEdx_Th_Minus_Overlay_Flipped"+tag+".png").c_str());
  c4-> Clear();

  c4->Close();

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
