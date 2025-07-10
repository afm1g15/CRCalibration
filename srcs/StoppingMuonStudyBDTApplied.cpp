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
#include "TView3D.h"
#include "TAxis3D.h"
#include "TLine.h"

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
   "trkpitchc_pandoraTrack",
   "hit_energy",
   "hit_startT",
   "hit_endT"
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
  TH1D *h_true_background_mother   = new TH1D("h_true_background_mother","",100,0,100);

  TH2F *h_reco_dQdx_RR = new TH2F("h_reco_dQdx_RR", ";Residual Range [cm];dQ/dx [ADC/cm]", 200, 0, 200, 100, 0, 1000);
  TH2D *h_RR_bin_MPVs   = new TH2D("h_RR_bin_MPVs","",40,0,200, 100, 0, 500);   // RR bin MP dqdx values
  TH2D *h_RR_bin_dEdx_MPVs   = new TH2D("h_RR_bin_dEdx_MPVs","",40,0,200, 600, 0, 6);   // RR bin MP dEdx values
  TH2D *h_pitch_vs_RR   = new TH2D("h_RR_vs_pitch",";Residual Range [cm];Pitch [cm]",200,0,200, 100, 0, 10);   // pitch vs RR
  TH2D *h_RR_bin_MPVs_Th   = new TH2D("h_RR_bin_MPVs_Th","",40,0,200, 600, 0, 6);   // RR bin MP dEdx values (theory)
  TH2D *h_RR_vs_Ratio   = new TH2D("h_RR_vs_Ratio","",40,0,200, 100, 0, 0.01);   // RR binned ratio
  TH2F *h_reco_dEdx_RR = new TH2F("h_reco_dEdx_RR", ";Residual Range [cm];dE/dx [MeV/cm]", 200, 0, 200, 100, 0, 10);


  int nbin = 40;
  int binsize = 5;
  //These are the latest fit values for the translation ratio in the AE methos
  double fitVar1 = 0.00462586;
  double fitVar2 = 0.00663633;
  double fitVar3 = 9.9079e-07;

  //Create the plots for RR bins to get the MPV values
  TH1D *dqdx[nbin];
  for (int i = 0; i < nbin; ++i)
  {
    if (i == 0)
      dqdx[i] = new TH1D(Form("dqdx_%d", i), "; dQ/dx [ADC/cm]; Number of entries", 100, 0.0, 2000);

    if (i != 0)
      dqdx[i] = new TH1D(Form("dqdx_%d", i), "; dQ/dx [ADC/cm]; Number of entries", 50, 0.0, 1000);

    dqdx[i]->SetLineColor(kBlack);
    dqdx[i]->GetXaxis()->SetTitleSize(0.04);
    dqdx[i]->GetYaxis()->SetTitleSize(0.04);
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
  unsigned int failAtStartY = 0;
  unsigned int failAtEndY = 0;
  unsigned int failAtEndYX = 0;
  unsigned int failAtEndYZ = 0;
  unsigned int failAtEndX = 0;
  unsigned int failAtEndZ = 0;
  unsigned int failAtKE = 0;
  unsigned int failAtRange = 0;
  unsigned int secondPass = 0;
  unsigned int firstPass = 0;
  unsigned int recoTrkIdLength = 0;
  unsigned int wrongByPDG = 0;
  unsigned int wrongByMother = 0;
  unsigned int wrongByWall = 0;
  unsigned int negTrueID = 0;
  unsigned int negTrueIDSignal = 0;
  unsigned int flippedRecoTracks = 0;
  unsigned int flippedSelectedRecoTracks = 0;


  // Now loop over the events
  unsigned int nEvts = tree->GetEntries();
  unsigned int iIt = 1;
  unsigned int trackRepeats = 0;
  unsigned int eventNum = 0;
  unsigned int dups_tot = 0;

  std::cout << "Total number of events = " << nEvts << std::endl;
  //Events numbers to be graphically printed
  std::vector<int> evtsToPrint = {128853};
  //Print all events out in one graphic?
  bool printSelectedEvts = true;

  TCanvas *c2 = new TCanvas("c2","",1000,1000);
    if (printSelectedEvts == true) {
      c2->cd();
      //Draw the Fid and Active Volumes
      for (int i=0; i < 3; i++) {
        double rmin[3] = {minx_fid[i],
                          miny_fid[i],
                          minz_fid[i]};
        double rmax[3] = {maxx_fid[i],
                          maxy_fid[i],
                          maxz_fid[i]};
        DrawCube(c2, rmin, rmax, 1, 3.);
        double rmin2[3] = {minx_av[i],
                           miny_av[i],
                           minz_av[i]};
        double rmax2[3] = {maxx_av[i],
                           maxy_av[i],
                           maxz_av[i]};
        DrawCube(c2, rmin2, rmax2, 1, 3.);
      }
    }


  std::cout << " |";
  for(unsigned int iEvt = 0; iEvt < nEvts; ++iEvt){
    tree->GetEntry(iEvt);
    if(!evtProc.SelectEvent(evt)) continue;
    bool printThisEvt = false;
    
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

    //Create a canvas and box to work with if it's an event number to print
    if (std::count(evtsToPrint.begin(), evtsToPrint.end(), iEvt)) {
      printThisEvt = true;
    }
 
    TCanvas *c1 = new TCanvas("c1","",1000,1000);
    if (printThisEvt == true) {
      c1->cd();
      //Draw the Fid and Active Volumes
      for (int i=0; i < 3; i++) {
        double rmin[3] = {minx_fid[i],
                          miny_fid[i],
                          minz_fid[i]};
        double rmax[3] = {maxx_fid[i],
                          maxy_fid[i],
                          maxz_fid[i]};
        DrawCube(c1, rmin, rmax, 1, 3.);
        double rmin2[3] = {minx_av[i],
                           miny_av[i],
                           minz_av[i]};
        double rmax2[3] = {maxx_av[i],
                           maxy_av[i],
                           maxz_av[i]};
        DrawCube(c1, rmin2, rmax2, 1, 3.);
      }
    }
  
    ///////////////////////////////////
    //          TRUTH                //
    ///////////////////////////////////          
    std::vector<int> trueTrkPassId; //vector of track IDs that pass true signal cuts
    unsigned int primaryInEvent = 0;
    // Now loop over the true tracks
    for(int iTrktru = 0; iTrktru < nGeant; ++iTrktru){

      // Count tracks
      totalTracksTrue++;

      TVector3 start(evt->StartPointx[iTrktru],evt->StartPointy[iTrktru],evt->StartPointz[iTrktru]);
      TVector3 end(evt->EndPointx[iTrktru],evt->EndPointy[iTrktru],evt->EndPointz[iTrktru]);

      if (printThisEvt) {
        c1->cd();
        TPolyLine3D *line = new TPolyLine3D(2);
        line->SetPoint(0, start.X(), start.Y(), start.Z());
        line->SetPoint(1, end.X(), end.Y(), end.Z());

        line->SetLineColor(6.);   //pink
        line->SetLineWidth(2.);
        line->Draw();
      }

      // Look for true pdg
      int trupdg = evt->pdg[iTrktru];
      if (abs(trupdg) != 13)
        continue;

      // Look for true mother, to check is primary
      int truMother = evt->Mother[iTrktru];
      if (truMother != 0) {
        trueSecondaryMuons++;
        continue;
      }
      else {
        primaryInEvent++;
      }

      // Check the the true end coordinates are within the TPC active volume
      // The general start and end points (including cryostat and TPC)

      // The tpc AV start and end points
      TVector3 endAV(evt->EndPointx_tpcAV[iTrktru],evt->EndPointy_tpcAV[iTrktru],evt->EndPointz_tpcAV[iTrktru]);

      // Get the differences between the two
      float dx = abs(endAV.X()-end.X());
      float dy = abs(endAV.Y()-end.Y());
      float dz = abs(endAV.Z()-end.Z());


      // If they don't match, it doesn't stop (i.e. it left the TPC so it's end point will be one of the walls)
      if (dx+dy+dz > 1e-10)
        continue;

      trueSignalMuons++;

      trueTrkPassId.push_back(evt->TrackId[iTrktru]);
      if (evt->TrackId[iTrktru] < 0) {
        negTrueIDSignal++;
      }

    } // iTrktru, truth loop


    ///////////////////////////////////
    //            RECO               //
    ///////////////////////////////////
    std::vector<int> recoTrkPassId; //vector of true track IDs that pass reco cuts
    std::vector<int> recoTrkId;
    std::vector<int> alreadyPassed;
    std::vector<int> allRecoTrkId;
    //Loop over all reco tracks
    for(int iTrk = 0; iTrk < nTrks; ++iTrk){

      // Get the track verticies points
      TVector3 startVtx(evt->trkstartx_pandoraTrack[iTrk],
                   evt->trkstarty_pandoraTrack[iTrk],
                   evt->trkstartz_pandoraTrack[iTrk]);
      TVector3 endVtx(evt->trkendx_pandoraTrack[iTrk],
                   evt->trkendy_pandoraTrack[iTrk],
                   evt->trkendz_pandoraTrack[iTrk]);

      if(startVtx.Y() < endVtx.Y()) {
        flippedRecoTracks++;
      }


      if (printThisEvt) {
        c1->cd();
        std::cout << "# tracks in event " << iEvt << " = " << nTrks << std::endl;
        TPolyLine3D *line = new TPolyLine3D(2);
        line->SetPoint(0, startVtx.X(), startVtx.Y(), startVtx.Z());
        line->SetPoint(1, endVtx.X(), endVtx.Y(), endVtx.Z());

        line->SetLineColor(4.);  //blue
        line->SetLineWidth(2.);
        line->Draw();
      }

      // Find the best plane (with most hits)
      int bestPlane = 0;
      std::vector<int> hitsOnPlane(3,0);
      GetRecoBestPlane(iTrk, evt, bestPlane, hitsOnPlane);

      //look out for delta rays with negative true IDs
      int trueID = evt->trkidtruth_pandoraTrack[iTrk][bestPlane];
      if (trueID < 0) {
        negTrueID++;
      }
      
      allRecoTrkId.push_back(trueID);
      totalTracksReco++;
      
      //Check the track only crosses one external plane
      float length = evt->trklen_pandoraTrack[iTrk];       //cm

      Plane enteringPlane = GetClosestPlane(extPlanes, startVtx, endVtx);
      double distFromEntrance = GetDistanceToPlane(enteringPlane, startVtx, endVtx);

      Plane exitingPlane = GetClosestPlane(extPlanes, endVtx, startVtx);
      double distFromExit = GetDistanceToPlane(exitingPlane, endVtx, startVtx);
      
      //Also need to ensure the track is not a fragment, so set a minimum length
      if (length < 20) {
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
      if ((thetaYZ > 0.0)) {
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
      if (nvtx > 20) {
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
      if (trkstartd > 20) {
        if(CheckTrueIDAssoc(trueID,trueTrkPassId)) {
          if(!CheckTrueIDAssoc(trueID,recoTrkId)) {
            failAtBoundDist++;
            std::cout << "event number (trkstartd) " << eventNum << " and reco trk id " << iTrk << std::endl;
            recoTrkId.push_back(trueID);
          }
        }
         continue;
      }

      
      if (startVtx.Y() < -100) {
        if(CheckTrueIDAssoc(trueID,trueTrkPassId)) {
          if(!CheckTrueIDAssoc(trueID,recoTrkId)) {
            failAtStartY++;
            std::cout << "event number (startVtxY) " << eventNum << " and reco trk id " << iTrk << std::endl;
            recoTrkId.push_back(trueID);
          }
        }
         continue;
      }


      if (endVtx.Y() < -550) {
        if(CheckTrueIDAssoc(trueID,trueTrkPassId)) {
          if(!CheckTrueIDAssoc(trueID,recoTrkId)) {
            failAtEndY++;
            std::cout << "event number (endVtxY) " << eventNum << " and reco trk id " << iTrk << std::endl;
            recoTrkId.push_back(trueID);
          }
        }
         continue;
      }
      
      if (( endVtx.X() < -355 || endVtx.X() > 355 )) {
        if(CheckTrueIDAssoc(trueID,trueTrkPassId)) {
          if(!CheckTrueIDAssoc(trueID,recoTrkId)) {
            failAtEndX++;
            std::cout << "event number (endVtxX) " << eventNum << " and reco trk id " << iTrk << std::endl;
            recoTrkId.push_back(trueID);
          }
        }
         continue;
      }

      if ( (endVtx.Z() < 50 || endVtx.Z() > 1350)) {
        if(CheckTrueIDAssoc(trueID,trueTrkPassId)) {
          if(!CheckTrueIDAssoc(trueID,recoTrkId)) {
            failAtEndZ++;
            std::cout << "event number (endVtxZ) " << eventNum << " and reco trk id " << iTrk << std::endl;
            recoTrkId.push_back(trueID);
          }
        }
         continue;
      }


      float recoKE = evt->trkke_pandoraTrack[iTrk][bestPlane];
      if ( recoKE < 150 ) {
        if(CheckTrueIDAssoc(trueID,trueTrkPassId)) {
          if(!CheckTrueIDAssoc(trueID,recoTrkId)) {
            failAtKE++;
            std::cout << "event number (KE) " << eventNum << " and reco trk id " << iTrk << std::endl;
            recoTrkId.push_back(trueID);
          }
        }
         continue;
      }
    
      float recoRange = evt->trkrange_pandoraTrack[iTrk][bestPlane];
      if ( recoRange < 60 ) {
        if(CheckTrueIDAssoc(trueID,trueTrkPassId)) {
          if(!CheckTrueIDAssoc(trueID,recoTrkId)) {
            failAtRange++;
            std::cout << "event number (range) " << eventNum << " and reco trk id " << iTrk << std::endl;
            recoTrkId.push_back(trueID);
          }
        }
         continue;
      }

      //apply the BDT here
      //Load in the model from the TMMA xml file
      TMVA::Experimental::RReader model("datasetBkg0/weights/TMVAMultiBkg0_BDTG.weights.xml");

      float thetaXZ = evt->trkthetaxz_pandoraTrack[iTrk];

      
      //Apply model
      auto prediction = model.Compute({static_cast<float>(nvtx), thetaXZ, thetaYZ, length, static_cast<float>(distFromEntrance), static_cast<float>(distFromExit), trkstartd, static_cast<float>(startVtx.Y()), static_cast<float>(endVtx.Y()), static_cast<float>(startVtx.X()), static_cast<float>(endVtx.X()), static_cast<float>(startVtx.Z()), static_cast<float>(endVtx.Z()), recoKE, recoRange});

      if (prediction[0] < -0.963129) {
        if(CheckTrueIDAssoc(trueID,trueTrkPassId)) {
          if(!CheckTrueIDAssoc(trueID,recoTrkId)) {
            failAtBDT++;
            std::cout << "event number (BDT) " << eventNum << " and reco trk id " << iTrk << std::endl;
            recoTrkId.push_back(trueID);
          }
        }
	continue;       
     }

      recoSelectedMuons++;
    

      //Now need to check how many of the selected tracks are also true signal
      recoTrkPassId.push_back(trueID);
      if(CheckTrueIDAssoc(trueID,recoTrkId) && CheckTrueIDAssoc(trueID,trueTrkPassId)) {
        int splitnum = std::count(recoTrkId.begin(), recoTrkId.end(), trueID);
        secondPass++;
        if (!CheckTrueIDAssoc(trueID,alreadyPassed)) {
          alreadyPassed.push_back(trueID);
        }
        else {
          dups_tot++;
        }
      }
      else if (!CheckTrueIDAssoc(trueID,recoTrkId) && CheckTrueIDAssoc(trueID,trueTrkPassId)) {
        firstPass++;
        recoTrkId.push_back(trueID);
        alreadyPassed.push_back(trueID);        
      }


      if ((printSelectedEvts)) {
        c2->cd();
        TPolyLine3D *line = new TPolyLine3D(2);
        line->SetPoint(0, startVtx.X(), startVtx.Y(), startVtx.Z());
        line->SetPoint(1, endVtx.X(), endVtx.Y(), endVtx.Z());

        if(CheckTrueIDAssoc(trueID,trueTrkPassId)) {
          line->SetLineColor(3.);
        }
        else {
          line->SetLineColor(2.);
        }
        line->SetLineWidth(2.);
        line->Draw();
      }

      /////////////////////////
      // TRUE SIG OR NOT?  ////
      // //////////////////////

      //if(CheckTrueIDAssoc(trueID,trueTrkPassId)) {  //if you want to use only true signal

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

          double hit_RR = evt->trkresrg_pandoraTrack[iTrk][bestPlane][iHit];
          if (!(hit_RR < 0.5 && corrected_dq_dx < 5)) {
            h_reco_dQdx_RR->Fill(hit_RR, corrected_dq_dx);
          }

          h_pitch_vs_RR->Fill(hit_RR, dp);

          double dEdx_corr = (fitVar1 + (fitVar2*(1/hit_RR)) + (fitVar3*hit_RR))*corrected_dq_dx;
          h_reco_dEdx_RR->Fill(hit_RR, dEdx_corr);
          
          //Put the hits in the correct RR bin for fitting
          if (bin < nbin)
          {
             dqdx[bin]->Fill(corrected_dq_dx);
          } 

        } //iHit

       ///////////////////////////
       //  TRUE OR NOT PT 2 //////
       //  ///////////////////////


      if(CheckTrueIDAssoc(trueID,trueTrkPassId)) {
        recoSelectedSignalMuons++;
      } //if selected signal
      else {
        //add to background pdg plot
        int truetrkpdg = evt->trkpdgtruth_pandoraTrack[iTrk][bestPlane];

        //this section will contain background incorrectly selected as signal
        if (abs(truetrkpdg) != 13) {
          wrongByPDG++;
        }
       
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
    } //if
     
    //see if any true tracks are missing a reco partner
    for(int iT = 0; iT < trueTrkPassId.size(); ++iT){
      int current = trueTrkPassId[iT];
      std::cout << "checking ID... " << current << std::endl;
      if (!CheckTrueIDAssoc(current,allRecoTrkId)) {
        std::cout << "This true signal track with trueID " << current << " in event " << iEvt << " doesn't seem to have a reco partner (selected or otherwise) " << std::endl;
      }
    }
    

    //Now save the plots if needed
    if (printThisEvt) {  
      c1->cd();
      TView3D *view = (TView3D*) TView::CreateView(1);

      double c[3] = { 0., 0., 3000. };   //centre
      double s[3] = { 1200., -900., 6100. };     //scale

      view->SetRange(-800, -650, -10, 800, 650, 6000);

      view->ToggleRulers();


      TAxis3D *axis = TAxis3D::GetPadAxis(gPad);
      axis->GetXaxis()->SetTitle("X (W)");
      axis->GetYaxis()->SetTitle("Y (Up)");
      axis->GetZaxis()->SetTitle("Z (N)");
      axis->GetXaxis()->SetAxisColor(kBlack);
      axis->GetYaxis()->SetAxisColor(kBlack);
      axis->GetZaxis()->SetAxisColor(kBlack);
      axis->GetXaxis()->SetLabelColor(kBlack);
      axis->GetYaxis()->SetLabelColor(kBlack);
      axis->GetZaxis()->SetLabelColor(kBlack);
      axis->GetXaxis()->SetLabelSize(0.024);
      axis->GetYaxis()->SetLabelSize(0.024);
      axis->GetZaxis()->SetLabelSize(0.024);


      view->DefineViewDirection(s, c,
                                0, 1,
                                1, 0,
                                1, 0,
                                view->GetTnorm(),
                                view->GetTback());


      axis->GetXaxis()->SetTitleOffset(2);
      axis->GetYaxis()->SetTitleOffset(-1.7);
      axis->GetXaxis()->SetLabelOffset(-0.065);
      axis->GetYaxis()->SetLabelOffset(-0.2);
      c1->SaveAs(Form("%u_front.png", iEvt));



      view->DefineViewDirection(s, c,
                                0, 1,
                                0, 1,
                                1, 0,
                                view->GetTnorm(),
                                view->GetTback());
      axis->GetXaxis()->SetTitleOffset(2);
      axis->GetZaxis()->SetTitleOffset(-1.7);
      axis->GetXaxis()->SetLabelOffset(-0.065);
      axis->GetZaxis()->SetLabelOffset(-0.2);
      c1->SaveAs(Form("%u_top.png", iEvt));

      view->DefineViewDirection(s, c,
                                1, 0,
                                0, 1,
                                0, 1,
                                view->GetTnorm(),
                                view->GetTback());

      axis->GetYaxis()->SetTitleOffset(-2);
      axis->GetZaxis()->SetTitleOffset(-1.7);
      axis->GetYaxis()->SetLabelOffset(0.005);
      axis->GetZaxis()->SetLabelOffset(0.005);
      c1->SaveAs(Form("%u_side.png", iEvt));
    }

    delete c1;

    eventNum++;
  }// Event loop
  std::cout << " --- 100 % --- |" << std::endl;

  //Now plot whole events
  if (printSelectedEvts) {
      c2->cd();
      TView3D *view = (TView3D*) TView::CreateView(1);

      double c[3] = { 0., 0., 700. };   //centre, 3000, 700
      double s[3] = { 1200., -900., 1500 };     //scale, 6100, 1500

      view->SetRange(-800, -650, -10, 800, 650, 1500);  //6000, 1500

      view->ToggleRulers();


      TAxis3D *axis = TAxis3D::GetPadAxis(gPad);
      axis->GetXaxis()->SetTitle("X [cm]");
      axis->GetYaxis()->SetTitle("Y [cm]");
      axis->GetZaxis()->SetTitle("Z [cm]");
      axis->GetXaxis()->SetAxisColor(kBlack);
      axis->GetYaxis()->SetAxisColor(kBlack);
      axis->GetZaxis()->SetAxisColor(kBlack);
      axis->GetXaxis()->SetLabelColor(kBlack);
      axis->GetYaxis()->SetLabelColor(kBlack);
      axis->GetZaxis()->SetLabelColor(kBlack);
      axis->GetXaxis()->SetLabelSize(0.024);
      axis->GetYaxis()->SetLabelSize(0.024);
      axis->GetZaxis()->SetLabelSize(0.024);


      view->DefineViewDirection(s, c,
                                0, 1,
                                1, 0,
                                1, 0,
                                view->GetTnorm(),
                                view->GetTback());


      axis->GetXaxis()->SetTitleOffset(2);
      axis->GetYaxis()->SetTitleOffset(-1.7);
      axis->GetXaxis()->SetLabelOffset(-0.065);
      axis->GetYaxis()->SetLabelOffset(-0.2);
      c2->SaveAs("full_front.png");

      view->DefineViewDirection(s, c,
                                0, 1,
                                0, 1,
                                1, 0,
                                view->GetTnorm(),
                                view->GetTback());
      axis->GetXaxis()->SetTitleOffset(2);
      axis->GetZaxis()->SetTitleOffset(-1.7);
      axis->GetXaxis()->SetLabelOffset(-0.065);
      axis->GetZaxis()->SetLabelOffset(-0.2);
      c2->SaveAs("full_top.png");

      view->DefineViewDirection(s, c,
                                1, 0,
                                0, 1,
                                0, 1,
                                view->GetTnorm(),
                                view->GetTback());

      axis->GetYaxis()->SetTitleOffset(-2);
      axis->GetZaxis()->SetTitleOffset(-1.7);
      axis->GetYaxis()->SetLabelOffset(0.005);
      axis->GetZaxis()->SetLabelOffset(0.005);
      c2->SaveAs("full_side.png");
  }

  delete c2;


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
  std::cout << " Signal Failures At Length            = " << failAtLength << std::endl;
  std::cout << " Signal Failures At Angle             = " << failAtAngle << std::endl;
  std::cout << " Signal Failures At # Vertex in event = " << failAtVertex << std::endl;
  std::cout << " Signal Failures At Boundry Distance  = " << failAtBoundDist << std::endl;
  std::cout << " Signal Failures At Start Y           = " << failAtStartY << std::endl;
  std::cout << " Signal Failures At End Y             = " << failAtEndY << std::endl;
  std::cout << " Signal Failures At End X             = " << failAtEndX << std::endl;
  std::cout << " Signal Failures At End Z             = " << failAtEndZ << std::endl;
  std::cout << " Signal Failures At KE                = " << failAtKE << std::endl;
  std::cout << " Signal Failures At Range             = " << failAtRange << std::endl;
  std::cout << " Signal Failures At BDT               = " << failAtBDT << std::endl;
  std::cout << "  " << std::endl;
  std::cout << " Signal Pass at 2nd (or more) split   = " << secondPass << std::endl;
  std::cout << " Signal Pass first time               = " << firstPass << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Signal Failures Sum                  = " << failAtLength + failAtAngle + failAtVertex + failAtBoundDist + failAtStartY + failAtEndY + failAtEndX + failAtEndZ + failAtKE + failAtRange + failAtBDT  << std::endl;
  std::cout << " Signal TrueIds seen in Reco          = " << recoTrkIdLength << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Reco selected background fails being signal by..." << std::endl;
  std::cout << "   PDG not 13         = " << wrongByPDG << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << "Reco tracks total with TrueID < 0               = " << negTrueID << std::endl;
  std::cout << "Signal tracks total with TrueID < 0             = " << negTrueIDSignal << std::endl;
  std::cout << "Reco tracks that would be flipped               = " << flippedRecoTracks << std::endl;
  std::cout << "Selected Reco tracks that would be flipped      = " << flippedSelectedRecoTracks << std::endl;

  //Now move on to the fitting
  std::cout << "Creating outfiles..." << std::endl;

  std::vector<double> mostProbValues;
  mostProbValues.resize(40, 0.0);
  std::vector<double> mostProbdEdxValues;
  mostProbdEdxValues.resize(40, 0.0);
  std::cout << "Number of bins = " << nbin << std::endl;


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
      minR = 400;
      maxR = 1200;
    }
    else {
      minR = 200; //dedx[i]->GetXaxis()->GetXmin();
      maxR = 800; //dedx[i]->GetXaxis()->GetXmax();
    }
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
    double dEdx_corr_MPV = (fitVar1 + (fitVar2*(1/range_for_bin)) + (fitVar3*range_for_bin))*mpv;
    mostProbdEdxValues[i] = dEdx_corr_MPV;

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
  }


  //-------------------------------------------------------------------
  //calculate theoretical values
 
  std::cout << "" << std::endl;
  std::cout << "Number of MPVs = " << mostProbValues.size() << std::endl;
  std::vector<double> thVals;
  thVals.resize(40, 0.0);

  for (int i = 0; i < nbin; i++)
  {
    //calculate theoretical KE for RR bin
    double range_for_bin = ((i+1)*binsize)-(binsize/2.0);
    double this_KE= muon_sp_range_to_KE -> Eval(range_for_bin); // == from rr

    std::cout << "Bin center RR = " << range_for_bin << std::endl;
    std::cout << "KE for this bin (from RR) = " << this_KE << std::endl;

    double pitch = 0.55; //Average per bin, compared to Rhiannon
    double this_MPV_dEdx = MPVdEdx(this_KE, pitch, mass_muon);
    std::cout << "dEdx MPV for this bin = " << this_MPV_dEdx << std::endl;

    h_RR_bin_MPVs_Th->Fill(range_for_bin, this_MPV_dEdx, 4);

    double ratio;
    if (mostProbValues[i] != 0.0)
    {    
      ratio = this_MPV_dEdx/mostProbValues[i];
      std::cout << "ratio for this ^ bin = " << ratio << std::endl;
    }
    else {
      ratio = 0.0;
      std::cout << "ratio for this ^ bin (failed) = " << ratio << std::endl;
    }

    //need to compare dEdx theory to dEdx reco, NOT dQdx, if wanted
    thVals[i] = this_MPV_dEdx;
    h_RR_vs_Ratio->Fill(range_for_bin, ratio);
    std::cout << "filled ratio bin i = " << i << std::endl;
    std::cout << "-----------------------------" << std::endl;
  }
  //-------------------------------------------------------------------
  std::cout << "Print most probable dEdx values in each bin..." << std::endl;
  for (int i = 0; i < nbin; i++)
  {
    std::cout << mostProbdEdxValues[i] << ", ";
  }

  std::cout << "Print most theoretical dEdx values in each bin..." << std::endl;
  for (int i = 0; i < nbin; i++)
  {
    std::cout << thVals[i] << ", ";
  }

  double sum = 0;
  std::cout << "Make some plots..." << std::endl;

  //-----------------------------------------------

  std::cout << "************************* Calibration.C has ended ***************************" << std::endl;


  //--------------------------------------------
  TCanvas *c3 = new TCanvas();
  h_reco_dQdx_RR->SetStats(0);
  h_reco_dQdx_RR->GetXaxis()->SetTitleSize(0.04);
  h_reco_dQdx_RR->GetYaxis()->SetTitleSize(0.04);
  h_reco_dQdx_RR->Draw("COLZ");
  h_reco_dQdx_RR->Write(" h_reco_dQdx_RR");
  c3->SaveAs("reco_dqdx_rr.png");
  c3->Write("reco_dQdx_RR");
  c3->Clear();

  //--------------------------------------------
  h_reco_dEdx_RR->SetStats(0);
  h_reco_dEdx_RR->GetXaxis()->SetTitleSize(0.04);
  h_reco_dEdx_RR->GetYaxis()->SetTitleSize(0.04);
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
  SetHistogramStyle2D(h_RR_bin_MPVs_Th,"Residual Range [cm]", "dE/dx MPV_{Theory} [MeV/cm]");
  h_RR_bin_MPVs_Th->Draw("box");
  h_RR_bin_MPVs_Th->SetLineWidth(3);
  h_RR_bin_MPVs_Th->SetLineColor(kTeal-5);
  h_RR_bin_MPVs_Th->SetMarkerSize(5);
  h_RR_bin_MPVs_Th->SetMarkerStyle(2);
  h_RR_bin_MPVs_Th->GetYaxis()->SetTitleOffset(0.95);
  ca->SaveAs((location+"/RR_MPV_dEdx_Th_in_bin"+tag+".png").c_str());
  ca->Clear();

  TCanvas *cb = new TCanvas();
  h_pitch_vs_RR->SetStats(0);
  h_pitch_vs_RR->GetXaxis()->SetTitleSize(0.04);
  h_pitch_vs_RR->GetYaxis()->SetTitleSize(0.04);
  h_pitch_vs_RR->Draw("COLZ");
  cb->SaveAs((location+"/RR_vs_pitch"+tag+".png").c_str());
  cb->Clear();

  //
  TCanvas *cc1 = new TCanvas("cc1","",900,900);
  SetCanvasStyle(cc1, 0.12,0.08,0.06,0.15,0,0,0);
  TF1 *f1 = new TF1("f1", "[0] + ([1]*(1/x)) + [2]*x", 0, 200);
  f1->SetParameters(0.,200.);
  f1->SetLineColor(kRed);
  h_RR_vs_Ratio->Fit(f1, "MR");

  SetHistogramStyle2D(h_RR_vs_Ratio,"Residual Range [cm]", "(dE/dx_{Theory})/(dQ/dx_{Reconstructed}) [MeV/ADC]");
  h_RR_vs_Ratio->Draw("box");
  h_RR_vs_Ratio->SetLineWidth(3);
  h_RR_vs_Ratio->SetLineColor(kTeal-5);
  h_RR_vs_Ratio->SetMarkerStyle(2);
  h_RR_vs_Ratio->SetMarkerSize(2);
  h_RR_vs_Ratio->GetYaxis()->SetTitleOffset(0.95);
  std::cout << "------------------------------------------Fit--" <<std::endl;
  f1->Draw("same");

  double par0 = f1->GetParameter(0);
  double par1 = f1->GetParameter(1);
  double par2 = f1->GetParameter(2);

  std::cout << "Found parameter values..." << std::endl;
  std::cout << "    Par 0 = " << par0 << std::endl;
  std::cout << "    Par 1 = " << par1 << std::endl;
  std::cout << "    Par 2 = " << par2 << std::endl;

  std::cout << " " << std::endl;

  cc1->SaveAs((location+"/RR_vs_Ratio"+tag+".png").c_str());
  cc1->SaveAs((location+"/RR_vs_Ratio"+tag+".root").c_str());
  cc1->Clear();

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
