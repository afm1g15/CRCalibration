/************************************************************************
 * 
 * A macro to plot various quantities, mostly at the truth-level in the 
 * CR muon calibration studies
 *
 * Plots for the calibration procedure are made using ActivityStudies.cpp
 *
 * Example file list located here:
 *   /home/jones/work/cosmics/LArSoft-v08_50_00/work/files/v09_41_00_02_files.list
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
   "process_primary",
   "ntracks_pandoraTrack",
   "trkId_pandoraTrack",
   "trkidtruth_pandoraTrack",
   "trkg4id_pandoraTrack",
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
   "theta_xz",
   "theta_yz",
   "theta",
   "phi",
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
     
int calibrationStudies(const char *config){

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
  n = evtProc.GetFiles();
  
  // Start of analysis (loop over chain and events
  std::cout << " Running analysis..." << std::endl;

  // Then setup the histograms, counters and any other variables to add to
  // Setup histograms
  // Truth-level track quantities
  TH1D *h_startX     = new TH1D("h_startX","",100,-770,770);  // Start X position of the muons
  TH1D *h_startY     = new TH1D("h_startY","",100,-630,630);  // Start Y position of the muons
  TH1D *h_startZ     = new TH1D("h_startZ","",100,-80,5900);  // Start Z position of the muons
  TH1D *h_endX       = new TH1D("h_endX","",100,-770,770);    // end X position of the muons
  TH1D *h_endY       = new TH1D("h_endY","",100,-630,630);    // end Y position of the muons
  TH1D *h_endZ       = new TH1D("h_endZ","",100,-80,5900);    // end Z position of the muons
  TH1D *h_length     = new TH1D("h_length","",100,0,2.2e3);   // Length of the muons
  TH1D *h_nDaughters = new TH1D("h_nDaughters","",100,0,100); // Number of muon daughters
  TH1D *h_mom        = new TH1D("h_mom","",100,3,1e5);        // Momentum of the muons [GeV]
  TH1D *h_energy     = new TH1D("h_energy","",100,3,1e5);     // Energy of the muons [GeV]
  SetLogX(h_mom);
  SetLogX(h_energy);
  
  TH1D *h_thetaxz    = new TH1D("h_thetaxz","",100,-3.2,3.2); // ThetaXZ
  TH1D *h_thetayz    = new TH1D("h_thetayz","",100,-3.2,0);   // ThetaYZ
  TH1D *h_theta      = new TH1D("h_theta","",100,0,3.2);      // Theta
  TH1D *h_phi        = new TH1D("h_phi","",100,-3.2,0);     // Phi
  TH1D *h_eDepPerL_0 = new TH1D("h_eDepPerL_0","",100,0,5); // Energy deposition per unit length
  TH1D *h_eDepPerL_1 = new TH1D("h_eDepPerL_1","",100,0,5); // Energy deposition per unit length
  TH1D *h_eDepPerL_2 = new TH1D("h_eDepPerL_2","",100,0,5); // Energy deposition per unit length
  TH2D *h_E_nDaught  = new TH2D("h_E_nDaught","",100,0,800,50,0,50); // Number of muon daughters per unit energy

  // Truth-level plane quantities
  TH1D *h_plane_cross  = new TH1D("h_plane_cross","",9,0,9); // Number of tracks crossing each plane
  TH1D *h_n_crossed    = new TH1D("h_n_crossed","",9,0,9); // Number of planes crossed by each track
  TH1D *h_plane_enter  = new TH1D("h_plane_enter","",9,0,9); // Number of tracks entering from each external plane
  TH1D *h_plane_exit   = new TH1D("h_plane_exit","",9,0,9); // Number of tracks exiting from each external plane
  TH1D *h_enter_dist   = new TH1D("h_enter_dist","",100,0,0.1); // Distance from candidate entrance plane 
  TH1D *h_exit_dist    = new TH1D("h_exit_dist","",100,0,0.1); // Distance from candidate exit plane

  // Hit-level depositions
  TH2D *h_hit_energy_x_0 = new TH2D("h_hit_energy_x_plane0","",100,-750,750,100,0,10); // Hit deposition energy vs X primary TPC muons, plane0
  TH2D *h_hit_charge_x_0 = new TH2D("h_hit_charge_x_plane0","",100,-750,750,100,0,500); // Hit deposition energy vs X primary TPC muons, plane0
  TH2D *h_hit_nelec_x_0  = new TH2D("h_hit_nelec_x_plane0","",100,-750,750,100,0,1.5e5); // Hit deposition energy vs X primary TPC muons, plane0
  TH2D *h_hit_energy_x_1 = new TH2D("h_hit_energy_x_plane1","",100,-750,750,100,0,10); // Hit deposition energy vs X primary TPC muons, plane1
  TH2D *h_hit_charge_x_1 = new TH2D("h_hit_charge_x_plane1","",100,-750,750,100,0,500); // Hit deposition energy vs X primary TPC muons, plane1
  TH2D *h_hit_nelec_x_1  = new TH2D("h_hit_nelec_x_plane1","",100,-750,750,100,0,1.5e5); // Hit deposition energy vs X primary TPC muons, plane1
  TH2D *h_hit_energy_x_2 = new TH2D("h_hit_energy_x_plane2","",100,-750,750,100,0,10); // Hit deposition energy vs X primary TPC muons, plane2
  TH2D *h_hit_charge_x_2 = new TH2D("h_hit_charge_x_plane2","",100,-750,750,100,0,1000); // Hit deposition energy vs X primary TPC muons, plane2
  TH2D *h_hit_nelec_x_2  = new TH2D("h_hit_nelec_x_plane2","",100,-750,750,100,0,2.5e5); // Hit deposition energy vs X primary TPC muons, plane2
  TH2D *h_corr_hit_energy_x_0 = new TH2D("h_corr_hit_energy_x_plane0","",100,-750,750,100,0,10); // Hit deposition energy vs X primary TPC muons, plane0
  TH2D *h_corr_hit_charge_x_0 = new TH2D("h_corr_hit_charge_x_plane0","",100,-750,750,100,0,500); // Hit deposition energy vs X primary TPC muons, plane0
  TH2D *h_corr_hit_nelec_x_0  = new TH2D("h_corr_hit_nelec_x_plane0","",100,-750,750,100,0,1.5e5); // Hit deposition energy vs X primary TPC muons, plane0
  TH2D *h_corr_hit_energy_x_1 = new TH2D("h_corr_hit_energy_x_plane1","",100,-750,750,100,0,10); // Hit deposition energy vs X primary TPC muons, plane1
  TH2D *h_corr_hit_charge_x_1 = new TH2D("h_corr_hit_charge_x_plane1","",100,-750,750,100,0,500); // Hit deposition energy vs X primary TPC muons, plane1
  TH2D *h_corr_hit_nelec_x_1  = new TH2D("h_corr_hit_nelec_x_plane1","",100,-750,750,100,0,1.5e5); // Hit deposition energy vs X primary TPC muons, plane1
  TH2D *h_corr_hit_energy_x_2 = new TH2D("h_corr_hit_energy_x_plane2","",100,-750,750,100,0,10); // Hit deposition energy vs X primary TPC muons, plane2
  TH2D *h_corr_hit_charge_x_2 = new TH2D("h_corr_hit_charge_x_plane2","",100,-750,750,100,0,1000); // Hit deposition energy vs X primary TPC muons, plane2
  TH2D *h_corr_hit_nelec_x_2  = new TH2D("h_corr_hit_nelec_x_plane2","",100,-750,750,100,0,2.5e5); // Hit deposition energy vs X primary TPC muons, plane2
  TH2D *h_corr_hit_nelec_energy_x_0  = new TH2D("h_corr_hit_nelec_energy_x_plane0","",100,-750,750,100,30e-6,43e-6); // Energy per electron [MeV/el] plane 0 
  TH2D *h_corr_hit_nelec_energy_x_1  = new TH2D("h_corr_hit_nelec_energy_x_plane1","",100,-750,750,100,30e-6,43e-6); // Energy per electron [MeV/el] plane 1
  TH2D *h_corr_hit_nelec_energy_x_2  = new TH2D("h_corr_hit_nelec_energy_x_plane2","",100,-750,750,100,30e-6,43e-6); // Energy per electron [MeV/el] plane 2
  TH2D *h_corr_hit_charge_energy_x_0  = new TH2D("h_corr_hit_charge_energy_x_plane0","",100,-750,750,100,5e-3,9e-3); // Energy per charge [MeV/ADC] plane 0 
  TH2D *h_corr_hit_charge_energy_x_1  = new TH2D("h_corr_hit_charge_energy_x_plane1","",100,-750,750,100,5e-3,9e-3); // Energy per charge [MeV/ADC] plane 1
  TH2D *h_corr_hit_charge_energy_x_2  = new TH2D("h_corr_hit_charge_energy_x_plane2","",100,-750,750,100,5e-3,9e-3); // Energy per charge [MeV/ADC] plane 2
  TH2D *h_corr_hit_nelec_corr_energy_x_0  = new TH2D("h_corr_hit_nelec_corr_energy_x_plane0","",100,-750,750,100,30e-6,43e-6); // Energy per electron [MeV/el] plane 0 
  TH2D *h_corr_hit_nelec_corr_energy_x_1  = new TH2D("h_corr_hit_nelec_corr_energy_x_plane1","",100,-750,750,100,30e-6,43e-6); // Energy per electron [MeV/el] plane 1
  TH2D *h_corr_hit_nelec_corr_energy_x_2  = new TH2D("h_corr_hit_nelec_corr_energy_x_plane2","",100,-750,750,100,30e-6,43e-6); // Energy per electron [MeV/el] plane 2
  TH2D *h_corr_hit_charge_corr_energy_x_0  = new TH2D("h_corr_hit_charge_corr_energy_x_plane0","",100,-750,750,100,5e-3,9e-3); // Energy per charge [MeV/ADC] plane 0 
  TH2D *h_corr_hit_charge_corr_energy_x_1  = new TH2D("h_corr_hit_charge_corr_energy_x_plane1","",100,-750,750,100,5e-3,9e-3); // Energy per charge [MeV/ADC] plane 1
  TH2D *h_corr_hit_charge_corr_energy_x_2  = new TH2D("h_corr_hit_charge_corr_energy_x_plane2","",100,-750,750,100,5e-3,9e-3); // Energy per charge [MeV/ADC] plane 2
  TH2D *h_corr_hit_nelec_corr_charge_x_0  = new TH2D("h_corr_hit_nelec_corr_charge_x_plane0","",100,-750,750,100,160,300); // charge per electron [MeV/el] plane 0 
  TH2D *h_corr_hit_nelec_corr_charge_x_1  = new TH2D("h_corr_hit_nelec_corr_charge_x_plane1","",100,-750,750,100,160,300); // charge per electron [MeV/el] plane 1
  TH2D *h_corr_hit_nelec_corr_charge_x_2  = new TH2D("h_corr_hit_nelec_corr_charge_x_plane2","",100,-750,750,100,160,300); // charge per electron [MeV/el] plane 2
  TH2D *h_eDep_nDaught_0 = new TH2D("h_eDep_nDaught_0","",50,0,50,100,0,5); // Number of muon daughters vs energy depositions per unit length
  TH2D *h_eDep_nDaught_1 = new TH2D("h_eDep_nDaught_1","",50,0,50,100,0,5); // Number of muon daughters vs energy depositions per unit length
  TH2D *h_eDep_nDaught_2 = new TH2D("h_eDep_nDaught_2","",50,0,50,100,0,5); // Number of muon daughters vs energy depositions per unit length
  TH2D *h_eDep_E_0 = new TH2D("h_eDep_E_0","",100,0,800,100,0,5); // Number of muon daughters vs energy depositions per unit length
  TH2D *h_eDep_E_1 = new TH2D("h_eDep_E_1","",100,0,800,100,0,5); // Number of muon daughters vs energy depositions per unit length
  TH2D *h_eDep_E_2 = new TH2D("h_eDep_E_2","",100,0,800,100,0,5); // Number of muon daughters vs energy depositions per unit length

  std::vector<TH2D*> h_energies{h_hit_energy_x_0,h_hit_energy_x_1,h_hit_energy_x_2};
  std::vector<TH2D*> h_charges{h_hit_charge_x_0,h_hit_charge_x_1,h_hit_charge_x_2};
  std::vector<TH2D*> h_nelecs{h_hit_nelec_x_0,h_hit_nelec_x_1,h_hit_nelec_x_2};
  
  std::vector<TH2D*> h_corr_energies{h_corr_hit_energy_x_0,h_corr_hit_energy_x_1,h_corr_hit_energy_x_2};
  std::vector<TH2D*> h_corr_charges{h_corr_hit_charge_x_0,h_corr_hit_charge_x_1,h_corr_hit_charge_x_2};
  std::vector<TH2D*> h_corr_nelecs{h_corr_hit_nelec_x_0,h_corr_hit_nelec_x_1,h_corr_hit_nelec_x_2};
  
  std::vector<TH2D*> h_corr_nelec_energies{h_corr_hit_nelec_energy_x_0,h_corr_hit_nelec_energy_x_1,h_corr_hit_nelec_energy_x_2};
  std::vector<TH2D*> h_corr_nelec_corr_energies{h_corr_hit_nelec_corr_energy_x_0,h_corr_hit_nelec_corr_energy_x_1,h_corr_hit_nelec_corr_energy_x_2};
  std::vector<TH2D*> h_corr_charge_energies{h_corr_hit_charge_energy_x_0,h_corr_hit_charge_energy_x_1,h_corr_hit_charge_energy_x_2};
  std::vector<TH2D*> h_corr_charge_corr_energies{h_corr_hit_charge_corr_energy_x_0,h_corr_hit_charge_corr_energy_x_1,h_corr_hit_charge_corr_energy_x_2};
  std::vector<TH2D*> h_corr_nelec_corr_charges{h_corr_hit_nelec_corr_charge_x_0,h_corr_hit_nelec_corr_charge_x_1,h_corr_hit_nelec_corr_charge_x_2};

  std::vector<TH1D*> h_eDepPerL{h_eDepPerL_0,h_eDepPerL_1,h_eDepPerL_2};
  std::vector<TH2D*> h_eDep_nDaught{h_eDep_nDaught_0,h_eDep_nDaught_1,h_eDep_nDaught_2};
  std::vector<TH2D*> h_eDep_E{h_eDep_E_0,h_eDep_E_1,h_eDep_E_2};
  
  // Setup counters
  unsigned int nMu        = 0;
  unsigned int nPrimaryMu = 0;
  unsigned int nLongMu    = 0;
  
  unsigned int noPlanes     = 0;
  unsigned int topBottom    = 0;
  unsigned int topOrBottom  = 0;
  unsigned int min2APACPA   = 0;
  unsigned int min1APACPA   = 0;
  unsigned int stopping     = 0;
  unsigned int exiting      = 0;
  unsigned int extCrossed2  = 0;
  
  // Setup counters
  int allHits      = 0;
  int hitId        = 0;
  int pandoraHitId = 0;

  // Now loop over the events
  unsigned int nEvts = tree->GetEntries();
  unsigned int iIt = 1;

  std::cout << " |";
  for(unsigned int iEvt = 0; iEvt < nEvts; ++iEvt){
    tree->GetEntry(iEvt);
    if(!evtProc.SelectEvent(evt)) continue;
    
    // Print the processing rate
    double evtFrac  = iEvt/static_cast<double>(nEvts);

    if((std::abs(0.1*iIt)-evtFrac) < std::numeric_limits<double>::epsilon()){
      std::cout << " --- " << evtFrac*100 << " %";
      std::cout.flush();
      iIt++;
    }

    // Geant and hit iterator definitions
    int nGeant = evt->geant_list_size;
    int nHits  = evt->no_hits_stored;

    // Loop over geant tracks
    for(int iG4 = 0; iG4 < nGeant; ++iG4){

      // Check the particle enters the TPC volume
      if(!evt->inTPCActive[iG4]) continue;

      // For the entrance tests
      // If the start or end locations are outside the detector, set them to be at the edge of the detector
      // Only do this for 1 coordinate
      TVector3 vtx(evt->StartPointx[iG4],evt->StartPointy[iG4],evt->StartPointz[iG4]);
      TVector3 end(evt->EndPointx[iG4],evt->EndPointy[iG4],evt->EndPointz[iG4]);
      TVector3 vtxAV(evt->StartPointx_tpcAV[iG4],evt->StartPointy_tpcAV[iG4],evt->StartPointz_tpcAV[iG4]);
      TVector3 endAV(evt->EndPointx_tpcAV[iG4],evt->EndPointy_tpcAV[iG4],evt->EndPointz_tpcAV[iG4]);

      // Determine if the track enters at the top and leaves through the bottom
      float vtxDy = abs(vtxAV.Y()-evt->StartPointy[iG4]);
      float endDy = abs(endAV.Y()-evt->EndPointy[iG4]);

      // If these don't match, the TPC end point and general end point are not same, 
      // and the locations are close to the top and bottom face, 
      // the particle goes through the top and bottom faces of the detector
      bool topAndBottom = false;
      if(vtxDy > std::numeric_limits<float>::epsilon() && endDy > std::numeric_limits<float>::epsilon()) {
        if(vtxAV.Y() > 599.5 && endAV.Y() < -599.5){
          topAndBottom = true;
        }
      }

      // Check for number of through-going and stopping muons in truth
      bool throughGoing = IsTrueThroughGoing(vtx,end,vtxAV,endAV);
      bool isStopping   = IsTrueStopping(vtx,end,vtxAV,endAV);

      int pdg        = evt->pdg[iG4];
      int id         = evt->TrackId[iG4];
      float lengthAV = (endAV-vtxAV).Mag();

      // Make sure we are looking at a primary muon
      if(abs(pdg) != 13) continue;
      nMu++;

      if(evt->Mother[iG4] != 0) continue;
      nPrimaryMu++;

      // Fill the truth-track quantities 
      h_startX->Fill(vtxAV.X());
      h_startY->Fill(vtxAV.Y());
      h_startZ->Fill(vtxAV.Z());
      h_endX->Fill(endAV.X());
      h_endY->Fill(endAV.Y());
      h_endZ->Fill(endAV.Z());
      h_length->Fill(lengthAV);
      h_mom->Fill(evt->P[iG4]);
      h_energy->Fill(evt->Eng[iG4]);
      h_nDaughters->Fill(evt->NumberDaughters[iG4]);
      h_phi->Fill(evt->phi[iG4]);
      h_theta->Fill(evt->theta[iG4]);
      h_thetaxz->Fill(evt->theta_xz[iG4]);
      h_thetayz->Fill(evt->theta_yz[iG4]);
      h_E_nDaught->Fill(evt->Eng[iG4],evt->NumberDaughters[iG4]);

      // Find the closest plane to the start vertex and count it as a crossing plane
      Plane enteringPlane = GetClosestPlane(extPlanes, vtxAV, endAV);
      double distFromEntrance = GetDistanceToPlane(enteringPlane, vtxAV, endAV);
      h_enter_dist->Fill(distFromEntrance);

      // Find the closest plane to the end vertex and count it as a crossing plane
      Plane exitingPlane = GetClosestPlane(extPlanes, endAV, vtxAV);
      double distFromExit = GetDistanceToPlane(exitingPlane, endAV, vtxAV);
      h_exit_dist->Fill(distFromExit);

      // Now check which planes are crossed
      unsigned int planeN    = 0;
      unsigned int extPlaneN = 0;
      
      // Counter for the number of planes this track has crossed
      unsigned int nPlanesCrossed = 0;
      unsigned int nExtCrossed    = 0;

      // Setup list of plane labels the track has crossed
      std::vector<std::string> labelsCrossed;

      // Loop over planes
      for(const Plane &pl : allPlanes){
        if(planeN > allPlanes.size()){
          std::cerr << " Error: Somehow the current plane iterator exceeds the number of planes in the detector: " << std::endl;
          std::cerr << " Iterator: " << planeN << " of " << allPlanes.size() << " total possible planes " << std::endl;
          std::exit(1);
        } // Debug
        // Check if this is the plane it (likely) entered the detector through 
        // Determine a maximum allowed distance from the plane to count it as the entrance
        if(enteringPlane.GetLabel() == pl.GetLabel()){
          if(distFromEntrance < 0.05){
            h_plane_cross->Fill(planeN);
            h_plane_enter->Fill(planeN);
            nPlanesCrossed++;
            labelsCrossed.push_back(pl.GetLabel());
          }
        }
        else if(exitingPlane.GetLabel() == pl.GetLabel()){
          if(!isStopping && distFromExit < 0.05){
            h_plane_cross->Fill(planeN);
            h_plane_exit->Fill(planeN);
            nPlanesCrossed++;
            labelsCrossed.push_back(pl.GetLabel());
          }
        }
        // Otherwise check if the track intersects the current plane
        else if(CheckIfIntersectsPlane(pl,vtxAV,endAV,lengthAV)){
          h_plane_cross->Fill(planeN);
          nPlanesCrossed++;
          labelsCrossed.push_back(pl.GetLabel());
        } // Intersects
        
        // Sort out the bin label
        h_plane_enter->GetXaxis()->SetBinLabel(planeN+1,planeLabels.find(pl.GetLabel())->second.c_str());
        h_plane_exit->GetXaxis()->SetBinLabel(planeN+1,planeLabels.find(pl.GetLabel())->second.c_str());
        h_plane_cross->GetXaxis()->SetBinLabel(planeN+1,planeLabels.find(pl.GetLabel())->second.c_str());
        ++planeN;
      } // Planes
      for(const Plane &pl : extPlanes){
        if(enteringPlane.GetLabel() == pl.GetLabel()){
          if(distFromEntrance < 0.05){
            nExtCrossed++;
          }
        } // Entrance
        else if(exitingPlane.GetLabel() == pl.GetLabel()){
          if(!isStopping && distFromExit < 0.05){
            nExtCrossed++;
          }
        } // Exit
        else if(CheckIfIntersectsPlane(pl,vtxAV,endAV,lengthAV)){
          nExtCrossed++;
        } // Intersects
        ++extPlaneN;
      } // ExtPlanes
     
      // Now fill the number of planes crossed histogram
      h_n_crossed->Fill(nPlanesCrossed);
      if(nPlanesCrossed == 0){
        noPlanes++;
      }
      // For the deposition studies, make sure we are looking at a long track (2m)
      if(lengthAV < 300) continue;
      nLongMu++;

      if(nExtCrossed >= 2) extCrossed2++;
      if(throughGoing) exiting++;
      if(isStopping) stopping++;
      
      // Loop over the labels crossed by this track and fill the appropriate counters
      bool APA = false;
      bool CPA = false;
      bool top = false;
      bool bot = false;

      for(std::string &str : labelsCrossed){
        // Get the converted label
        std::string longLab = planeLabels.find(str)->second;
        TString lab(longLab);
        if(lab.Contains("APA"))
          APA = true;
        else if(lab.Contains("CPA"))
          CPA = true;
        else if(lab.Contains("Top"))
          top = true;
        else if(lab.Contains("Bot"))
          bot = true;
      }
      if(APA || CPA)
        min1APACPA++;
      if(APA && CPA)
        min2APACPA++;
      if(top || bot)
        topOrBottom++;
      if(topAndBottom){
        topBottom++;
      }

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
          float hitEl     = evt->hit_nelec[iHit];

          totalEDep += hitE;
          totalQDep += hitQ;

          // Check if x is lower than the APA bound, charge seems to accumulate there
          if(hitX < evtProc.APA_X_POSITIONS[0] || hitX > evtProc.APA_X_POSITIONS[2]) continue;

          // Lifetime correction
          int tpc =evtProc.WhichTPC(hitX) + 1;
          float dx = ( -1 + 2*(tpc%2) )*(hitX - evtProc.APA_X_POSITIONS[tpc/2]);
          float dt = dx*evtProc.kXtoT;

          float corr      = TMath::Exp(-dt/2.88);
          float nCorr     = TMath::Exp(-dt/3);
          float eCorr     = TMath::Exp(-dt/2.88) / TMath::Exp(-dt/3); // Correct for the already-corrected energy
          float hitECorr  = hitE/eCorr;
          float hitQCorr  = hitQ/corr;
          float hitElCorr = hitEl/nCorr;///corr; // This may not be correct

          float hitENEl     = hitE/hitElCorr;
          float hitEQ       = hitE/hitQCorr;
          float hitENElCorr = hitECorr/hitElCorr;
          float hitEQCorr   = hitECorr/hitQCorr;
          float hitNElQCorr = hitElCorr/hitQCorr;

          // Now fill some histograms
          h_energies.at(iWire)->Fill(hitX,hitE);
          h_charges.at(iWire)->Fill(hitX,hitQ);
          h_nelecs.at(iWire)->Fill(hitX,hitEl);

          // Corrected
          h_corr_energies.at(iWire)->Fill(hitX,hitECorr);
          h_corr_charges.at(iWire)->Fill(hitX,hitQCorr);
          h_corr_nelecs.at(iWire)->Fill(hitX,hitElCorr);

          // Ratio
          h_corr_nelec_energies.at(iWire)->Fill(hitX,hitENEl);
          h_corr_nelec_corr_energies.at(iWire)->Fill(hitX,hitENElCorr);
          h_corr_charge_energies.at(iWire)->Fill(hitX,hitEQ);
          h_corr_charge_corr_energies.at(iWire)->Fill(hitX,hitEQCorr);
          h_corr_nelec_corr_charges.at(iWire)->Fill(hitX,hitNElQCorr);

        }// Hits
        float totalEDepPerLength = totalEDep/lengthAV;
        if(totalEDepPerLength < std::numeric_limits<float>::epsilon()) continue;
        h_eDep_nDaught.at(iWire)->Fill(evt->NumberDaughters[iG4],totalEDepPerLength);
        h_eDep_E.at(iWire)->Fill(evt->Eng[iG4],totalEDepPerLength);
        h_eDepPerL.at(iWire)->Fill(totalEDepPerLength);
      } // Wire plane
    }// Geant

  }// Event loop
  std::cout << " --- 100 % --- |" << std::endl;

  // Sort out the TeX file
  std::vector<std::string> contents{
    "Events",
    "TPC $\\mu$",
    "Primary TPC $\\mu$",
    "Long, primary TPC $\\mu$",
    "Crosses top or bottom",
    "Crosses top and bottom",
    "Crosses $ \\geq $ 1 APA/CPA",
    "Crosses $ \\geq $ 2 APA/CPA",
    "Stopping",
    "Exiting"
  };
  std::vector<unsigned int> rates{
    nEvts,
    nMu,
    nPrimaryMu,
    nLongMu,
    topOrBottom,
    topBottom,
    min1APACPA,
    min2APACPA,
    stopping,
    exiting
  };

  ofstream texFile;
  texFile.open(location+"/truth_contents"+tag+".tex");
  WriteStatsToTeX(texFile, n, contents, rates, static_cast<double>(nLongMu), "Long, primary TPC $\\mu$");

  // Plane crossing
  TCanvas *c1 = new TCanvas("c1","",900,900);
  SetCanvasStyle(c1, 0.12,0.05,0.06,0.15,0,0,0);

  TLegend *l = new TLegend(0.22,0.94,0.98,0.995);
  l->SetNColumns(3);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->SetTextFont(132);

  SetHistogramStyle1D(h_plane_cross,"Plane label", " Number of tracks crossing plane");
  h_plane_cross->Draw("hist");
  h_plane_enter->Draw("same");
  h_plane_exit->Draw("same");
  h_plane_cross->SetLineWidth(3);
  h_plane_enter->SetLineWidth(3);
  h_plane_exit->SetLineWidth(3);
  h_plane_cross->SetLineStyle(2);
  h_plane_enter->SetLineStyle(3);
  h_plane_exit->SetLineStyle(4);
  h_plane_enter->SetLineColor(pal.at(0));
  h_plane_exit->SetLineColor(pal.at(1));
  h_plane_cross->SetLineColor(pal.at(2));
  h_plane_cross->LabelsOption("v");
  h_plane_cross->GetXaxis()->SetTitleOffset(1.4);
  h_plane_cross->GetYaxis()->SetTitleOffset(0.95);

  l->AddEntry(h_plane_cross, "Crossed planes", "l");
  l->AddEntry(h_plane_enter, "Entry planes", "l");
  l->AddEntry(h_plane_exit,  "Exit planes", "l");
  l->Draw("same");

  c1->SaveAs((location+"/truth_planes_crossed_entered_exited"+tag+".png").c_str());
  c1->SaveAs((location+"/truth_planes_crossed_entered_exited"+tag+".root").c_str());
  c1->Clear();
  l->Clear();
  
  TCanvas *c2 = new TCanvas("c2","",900,900);
  SetCanvasStyle(c2, 0.12,0.08,0.06,0.12,0,1,0);

  l->SetNColumns(2);
  l->SetX1NDC(0.47);

  SetHistogramStyle1D(h_startX,"Muon X [cm]", "Rate");
  h_startX->Draw("hist");
  h_endX->Draw("same");
  h_startX->SetLineWidth(3);
  h_startX->SetLineStyle(2);
  h_endX->SetLineWidth(3);
  h_endX->SetLineStyle(2);
  h_startX->SetLineColor(pal.at(0));
  h_endX->SetLineColor(pal.at(1));
  h_startX->GetYaxis()->SetTitleOffset(0.95);

  l->AddEntry(h_startX, "Start", "l");
  l->AddEntry(h_endX, "End", "l");
  l->Draw();
  
  c2->SaveAs((location+"/truth_tracks_start_endX"+tag+".png").c_str());
  c2->SaveAs((location+"/truth_tracks_start_endX"+tag+".root").c_str());
  c2->Clear();
  l->Clear();

  SetHistogramStyle1D(h_startY,"Muon Y [cm]", "Rate");
  h_startY->Draw("hist");
  h_endY->Draw("same");
  h_startY->SetLineWidth(3);
  h_startY->SetLineStyle(2);
  h_endY->SetLineWidth(3);
  h_endY->SetLineStyle(2);
  h_startY->SetLineColor(pal.at(0));
  h_endY->SetLineColor(pal.at(1));
  h_startY->GetYaxis()->SetTitleOffset(0.95);

  l->AddEntry(h_startY, "Start", "l");
  l->AddEntry(h_endY, "End", "l");
  l->Draw();
  
  c2->SaveAs((location+"/truth_tracks_start_endY"+tag+".png").c_str());
  c2->SaveAs((location+"/truth_tracks_start_endY"+tag+".root").c_str());
  c2->Clear();
  l->Clear();

  SetHistogramStyle1D(h_startZ,"Muon Z [cm]", "Rate");
  h_startZ->Draw("hist");
  h_endZ->Draw("same");
  h_startZ->SetLineWidth(3);
  h_startZ->SetLineStyle(2);
  h_endZ->SetLineWidth(3);
  h_endZ->SetLineStyle(2);
  h_startZ->SetLineColor(pal.at(0));
  h_endZ->SetLineColor(pal.at(1));
  h_startZ->GetYaxis()->SetTitleOffset(0.95);

  l->AddEntry(h_startZ, "Start", "l");
  l->AddEntry(h_endZ, "End", "l");
  l->Draw();
  
  c2->SaveAs((location+"/truth_tracks_start_endZ"+tag+".png").c_str());
  c2->SaveAs((location+"/truth_tracks_start_endZ"+tag+".root").c_str());
  c2->Clear();
  l->Clear();

  TCanvas *c3 = new TCanvas("c3","",900,900);
  SetCanvasStyle(c3, 0.12,0.08,0.06,0.12,0,0,0);

  SetHistogramStyle1D(h_n_crossed,"Number of planes crossed [P]", " Number of tracks crossing P planes");
  h_n_crossed->Draw("hist");
  h_n_crossed->SetLineWidth(3);
  h_n_crossed->SetLineColor(pal.at(0));
  h_n_crossed->GetYaxis()->SetTitleOffset(0.95);
  c3->SaveAs((location+"/truth_tracks_crossed_nplanes"+tag+".png").c_str());
  c3->SaveAs((location+"/truth_tracks_crossed_nplanes"+tag+".root").c_str());
  c3->Clear();

  SetHistogramStyle1D(h_length,"Muon length [cm]", "Rate");
  h_length->Draw("hist");
  h_length->SetLineWidth(3);
  h_length->SetLineColor(pal.at(0));
  h_length->GetYaxis()->SetTitleOffset(0.95);
  c3->SaveAs((location+"/truth_tracks_length"+tag+".png").c_str());
  c3->SaveAs((location+"/truth_tracks_length"+tag+".root").c_str());
  c3->Clear();

  SetHistogramStyle1D(h_nDaughters,"Muon daughters", "Rate");
  h_nDaughters->Draw("hist");
  h_nDaughters->SetLineWidth(3);
  h_nDaughters->SetLineColor(pal.at(0));
  h_nDaughters->GetYaxis()->SetTitleOffset(0.95);
  c3->SaveAs((location+"/truth_tracks_nDaughters"+tag+".png").c_str());
  c3->SaveAs((location+"/truth_tracks_nDaughters"+tag+".root").c_str());
  c3->Clear();

  SetHistogramStyle1D(h_thetaxz,"Muon #theta_{XZ} [rad]", "Rate");
  h_thetaxz->Draw("hist");
  h_thetaxz->SetLineWidth(3);
  h_thetaxz->SetLineColor(pal.at(0));
  h_thetaxz->GetYaxis()->SetTitleOffset(0.95);
  c3->SaveAs((location+"/truth_tracks_thetaxz"+tag+".png").c_str());
  c3->SaveAs((location+"/truth_tracks_thetaxz"+tag+".root").c_str());
  c3->Clear();

  SetHistogramStyle1D(h_thetayz,"Muon #theta_{YZ} [rad]", "Rate");
  h_thetayz->Draw("hist");
  h_thetayz->SetLineWidth(3);
  h_thetayz->SetLineColor(pal.at(0));
  h_thetayz->GetYaxis()->SetTitleOffset(0.95);
  c3->SaveAs((location+"/truth_tracks_thetayz"+tag+".png").c_str());
  c3->SaveAs((location+"/truth_tracks_thetayz"+tag+".root").c_str());
  c3->Clear();

  SetHistogramStyle1D(h_theta,"Muon #theta [rad]", "Rate");
  h_theta->Draw("hist");
  h_theta->SetLineWidth(3);
  h_theta->SetLineColor(pal.at(0));
  h_theta->GetYaxis()->SetTitleOffset(0.95);
  c3->SaveAs((location+"/truth_tracks_theta"+tag+".png").c_str());
  c3->SaveAs((location+"/truth_tracks_theta"+tag+".root").c_str());
  c3->Clear();

  SetHistogramStyle1D(h_phi,"Muon #phi [rad]", "Rate");
  h_phi->Draw("hist");
  h_phi->SetLineWidth(3);
  h_phi->SetLineColor(pal.at(0));
  h_phi->GetYaxis()->SetTitleOffset(0.95);
  c3->SaveAs((location+"/truth_tracks_phi"+tag+".png").c_str());
  c3->SaveAs((location+"/truth_tracks_phi"+tag+".root").c_str());
  c3->Clear();

  SetHistogramStyle1D(h_enter_dist,"Distance from candidate entrance/exit [cm]", " Rate");
  h_enter_dist->Draw("hist");
  h_exit_dist->Draw("same");
  h_enter_dist->SetLineWidth(3);
  h_exit_dist->SetLineWidth(3);
  h_enter_dist->SetLineStyle(2);
  h_exit_dist->SetLineStyle(3);
  h_enter_dist->SetLineColor(pal.at(1));
  h_exit_dist->SetLineColor(pal.at(0));
  h_enter_dist->GetYaxis()->SetTitleOffset(0.95);

  l->AddEntry(h_enter_dist, "Entrance", "l");
  l->AddEntry(h_exit_dist, "Exit", "l");
  l->Draw();

  c3->SaveAs((location+"/truth_distance_to_entrance_exit_planes"+tag+".png").c_str());
  c3->SaveAs((location+"/truth_distance_to_entrance_exit_planes"+tag+".root").c_str());
  c3->Clear();
  l->Clear();
  
  for(unsigned int iWire = 0; iWire < 3; ++iWire){
    SetHistogramStyle1D(h_eDepPerL.at(iWire),"Energy deposition per unit length [MeV/cm]", "Rate");
    h_eDepPerL.at(iWire)->Draw("hist");
    h_eDepPerL.at(iWire)->SetLineWidth(3);
    h_eDepPerL.at(iWire)->SetLineColor(pal.at(0));
    h_eDepPerL.at(iWire)->GetYaxis()->SetTitleOffset(0.95);
    c3->SaveAs((location+"/truth_tracks_eDep_perL"+tag+".png").c_str());
    c3->SaveAs((location+"/truth_tracks_eDep_perL"+tag+".root").c_str());
    c3->Clear();
  }

  c3->SetLogx();
  c3->SetLogy();
  SetHistogramStyle1D(h_mom,"Muon momentum [GeV]", "Rate");
  h_mom->Draw("hist");
  h_mom->SetLineWidth(3);
  h_mom->SetLineColor(pal.at(0));
  h_mom->GetYaxis()->SetTitleOffset(0.95);
  c3->SaveAs((location+"/truth_tracks_mom"+tag+".png").c_str());
  c3->SaveAs((location+"/truth_tracks_mom"+tag+".root").c_str());
  c3->Clear();

  SetHistogramStyle1D(h_energy,"Muon energy [GeV]", "Rate");
  h_energy->Draw("hist");
  h_energy->SetLineWidth(3);
  h_energy->SetLineColor(pal.at(0));
  h_energy->GetYaxis()->SetTitleOffset(0.95);
  c3->SaveAs((location+"/truth_tracks_energy"+tag+".png").c_str());
  c3->SaveAs((location+"/truth_tracks_energy"+tag+".root").c_str());
  c3->Clear();

  
  TCanvas *c4 = new TCanvas("c4","",1000,800);
  SetCanvasStyle(c4, 0.1,0.12,0.05,0.12,0,0,0);

  std::vector<TLine*> APACPALines;
  // Now draw lines and labels where the APA and CPAs are
  for(unsigned int iA = 0; iA < 3; ++iA){
    TLine *l = new TLine(evtProc.APA_X_POSITIONS[iA], 0, evtProc.APA_X_POSITIONS[iA], 10);
    l->SetLineColor(kWhite);
    l->SetLineWidth(3);
    l->SetLineStyle(2);
    APACPALines.push_back(l);
  }
  for(unsigned int iC = 0; iC < 2; ++iC){
    TLine *l = new TLine(evtProc.CPA_X_POSITIONS[iC], 0, evtProc.CPA_X_POSITIONS[iC], 10);
    l->SetLineColor(kWhite);
    l->SetLineWidth(3);
    l->SetLineStyle(2);
    APACPALines.push_back(l);
  }

  for(unsigned int iWire = 0; iWire < 3; ++iWire){
    // Energies
    SetHistogramStyle2D(h_energies.at(iWire),"x [cm]", " Energy deposition [MeV/cm]",false);
    h_energies.at(iWire)->GetZaxis()->SetLabelSize(0.03);
    h_energies.at(iWire)->GetZaxis()->SetLabelFont(132);
    h_energies.at(iWire)->Draw("colz");

    // Draw the APA and CPA lines and labels
    for(unsigned int iLine = 0; iLine < APACPALines.size(); ++iLine){
      APACPALines.at(iLine)->SetY1(0);
      APACPALines.at(iLine)->SetY2(10);
      APACPALines.at(iLine)->Draw();
    }

    FormatLatex(evtProc.APA_X_POSITIONS[0]+10,9, "#color[0]{APA}");
    FormatLatex(evtProc.CPA_X_POSITIONS[0]+10,9, "#color[0]{CPA}");
    FormatLatex(evtProc.APA_X_POSITIONS[1]+10,9, "#color[0]{APA}");
    FormatLatex(evtProc.CPA_X_POSITIONS[1]+10,9, "#color[0]{CPA}");
    FormatLatex(evtProc.APA_X_POSITIONS[2]+10,9, "#color[0]{APA}");

    c4->SaveAs((location+"/energy_vs_X_plane"+std::to_string(iWire)+tag+".png").c_str());
    c4->SaveAs((location+"/energy_vs_X_plane"+std::to_string(iWire)+tag+".root").c_str());
    c4->Clear();
    
    SetHistogramStyle2D(h_corr_energies.at(iWire),"x [cm]", " Energy deposition [MeV/cm]",false);
    h_corr_energies.at(iWire)->GetZaxis()->SetLabelSize(0.03);
    h_corr_energies.at(iWire)->GetZaxis()->SetLabelFont(132);
    h_corr_energies.at(iWire)->Draw("colz");

    // Draw the APA and CPA lines and labels
    for(unsigned int iLine = 0; iLine < APACPALines.size(); ++iLine){
      APACPALines.at(iLine)->Draw();
    }

    FormatLatex(evtProc.APA_X_POSITIONS[0]+10,9, "#color[0]{APA}");
    FormatLatex(evtProc.CPA_X_POSITIONS[0]+10,9, "#color[0]{CPA}");
    FormatLatex(evtProc.APA_X_POSITIONS[1]+10,9, "#color[0]{APA}");
    FormatLatex(evtProc.CPA_X_POSITIONS[1]+10,9, "#color[0]{CPA}");
    FormatLatex(evtProc.APA_X_POSITIONS[2]+10,9, "#color[0]{APA}");

    c4->SaveAs((location+"/corr_energy_vs_X_plane"+std::to_string(iWire)+tag+".png").c_str());
    c4->SaveAs((location+"/corr_energy_vs_X_plane"+std::to_string(iWire)+tag+".root").c_str());
    c4->Clear();
    
    SetHistogramStyle2D(h_corr_nelec_energies.at(iWire),"x [cm]", " Energy deposition per e^{-} [MeV/e^{-}]",false);
    h_corr_nelec_energies.at(iWire)->GetZaxis()->SetLabelSize(0.03);
    h_corr_nelec_energies.at(iWire)->GetZaxis()->SetLabelFont(132);
    h_corr_nelec_energies.at(iWire)->Draw("colz");

    // Draw the APA and CPA lines and labels
    for(unsigned int iLine = 0; iLine < APACPALines.size(); ++iLine){
      APACPALines.at(iLine)->SetY1(30e-6);
      APACPALines.at(iLine)->SetY2(43e-6);
      APACPALines.at(iLine)->Draw();
    }

    FormatLatex(evtProc.APA_X_POSITIONS[0]+10,41.5e-6, "#color[0]{APA}");
    FormatLatex(evtProc.CPA_X_POSITIONS[0]+10,41.5e-6, "#color[0]{CPA}");
    FormatLatex(evtProc.APA_X_POSITIONS[1]+10,41.5e-6, "#color[0]{APA}");
    FormatLatex(evtProc.CPA_X_POSITIONS[1]+10,41.5e-6, "#color[0]{CPA}");
    FormatLatex(evtProc.APA_X_POSITIONS[2]+10,41.5e-6, "#color[0]{APA}");

    c4->SaveAs((location+"/energy_nelec_vs_X_plane"+std::to_string(iWire)+tag+".png").c_str());
    c4->SaveAs((location+"/energy_nelec_vs_X_plane"+std::to_string(iWire)+tag+".root").c_str());
    c4->Clear();
    
    SetHistogramStyle2D(h_corr_nelec_corr_energies.at(iWire),"x [cm]", " Energy deposition per e^{-} [MeV/e^{-}]",false);
    h_corr_nelec_corr_energies.at(iWire)->GetZaxis()->SetLabelSize(0.03);
    h_corr_nelec_corr_energies.at(iWire)->GetZaxis()->SetLabelFont(132);
    h_corr_nelec_corr_energies.at(iWire)->Draw("colz");

    // Draw the APA and CPA lines and labels
    for(unsigned int iLine = 0; iLine < APACPALines.size(); ++iLine){
      APACPALines.at(iLine)->Draw();
    }

    FormatLatex(evtProc.APA_X_POSITIONS[0]+10,41.5e-6, "#color[0]{APA}");
    FormatLatex(evtProc.CPA_X_POSITIONS[0]+10,41.5e-6, "#color[0]{CPA}");
    FormatLatex(evtProc.APA_X_POSITIONS[1]+10,41.5e-6, "#color[0]{APA}");
    FormatLatex(evtProc.CPA_X_POSITIONS[1]+10,41.5e-6, "#color[0]{CPA}");
    FormatLatex(evtProc.APA_X_POSITIONS[2]+10,41.5e-6, "#color[0]{APA}");

    c4->SaveAs((location+"/energy_nelec_corr_vs_X_plane"+std::to_string(iWire)+tag+".png").c_str());
    c4->SaveAs((location+"/energy_nelec_corr_vs_X_plane"+std::to_string(iWire)+tag+".root").c_str());
    c4->Clear();
    
    SetHistogramStyle2D(h_corr_charge_energies.at(iWire),"x [cm]", " dE/dQ [MeV/ADC]",false);
    h_corr_charge_energies.at(iWire)->GetZaxis()->SetLabelSize(0.03);
    h_corr_charge_energies.at(iWire)->GetZaxis()->SetLabelFont(132);
    h_corr_charge_energies.at(iWire)->Draw("colz");

    // Draw the APA and CPA lines and labels
    for(unsigned int iLine = 0; iLine < APACPALines.size(); ++iLine){
      APACPALines.at(iLine)->SetY1(5e-3);
      APACPALines.at(iLine)->SetY2(9e-3);
      APACPALines.at(iLine)->Draw();
    }

    FormatLatex(evtProc.APA_X_POSITIONS[0]+10,8.5e-3, "#color[0]{APA}");
    FormatLatex(evtProc.CPA_X_POSITIONS[0]+10,8.5e-3, "#color[0]{CPA}");
    FormatLatex(evtProc.APA_X_POSITIONS[1]+10,8.5e-3, "#color[0]{APA}");
    FormatLatex(evtProc.CPA_X_POSITIONS[1]+10,8.5e-3, "#color[0]{CPA}");
    FormatLatex(evtProc.APA_X_POSITIONS[2]+10,8.5e-3, "#color[0]{APA}");

    c4->SaveAs((location+"/energy_charge_vs_X_plane"+std::to_string(iWire)+tag+".png").c_str());
    c4->SaveAs((location+"/energy_charge_vs_X_plane"+std::to_string(iWire)+tag+".root").c_str());
    c4->Clear();
    
    SetHistogramStyle2D(h_corr_charge_corr_energies.at(iWire),"x [cm]", " dE/dQ [MeV/ADC]",false);
    h_corr_charge_corr_energies.at(iWire)->GetZaxis()->SetLabelSize(0.03);
    h_corr_charge_corr_energies.at(iWire)->GetZaxis()->SetLabelFont(132);
    h_corr_charge_corr_energies.at(iWire)->Draw("colz");

    // Draw the APA and CPA lines and labels
    for(unsigned int iLine = 0; iLine < APACPALines.size(); ++iLine){
      APACPALines.at(iLine)->Draw();
    }

    FormatLatex(evtProc.APA_X_POSITIONS[0]+10,8.5e-3, "#color[0]{APA}");
    FormatLatex(evtProc.CPA_X_POSITIONS[0]+10,8.5e-3, "#color[0]{CPA}");
    FormatLatex(evtProc.APA_X_POSITIONS[1]+10,8.5e-3, "#color[0]{APA}");
    FormatLatex(evtProc.CPA_X_POSITIONS[1]+10,8.5e-3, "#color[0]{CPA}");
    FormatLatex(evtProc.APA_X_POSITIONS[2]+10,8.5e-3, "#color[0]{APA}");

    c4->SaveAs((location+"/energy_charge_corr_vs_X_plane"+std::to_string(iWire)+tag+".png").c_str());
    c4->SaveAs((location+"/energy_charge_corr_vs_X_plane"+std::to_string(iWire)+tag+".root").c_str());
    c4->Clear();
    
    // Charges
    SetHistogramStyle2D(h_charges.at(iWire),"x [cm]", " Charge deposition [ADC/cm]",false);
    h_charges.at(iWire)->GetZaxis()->SetLabelSize(0.03);
    h_charges.at(iWire)->GetZaxis()->SetLabelFont(132);
    h_charges.at(iWire)->Draw("colz");

    // Draw the APA and CPA lines and labels
    for(unsigned int iLine = 0; iLine < APACPALines.size(); ++iLine){
      APACPALines.at(iLine)->SetY1(0);
      if(iWire == 2)
        APACPALines.at(iLine)->SetY2(1000);
      else
        APACPALines.at(iLine)->SetY2(500);
      APACPALines.at(iLine)->Draw();
    }

    if(iWire == 2){
      FormatLatex(evtProc.APA_X_POSITIONS[0]+10,900, "#color[0]{APA}");
      FormatLatex(evtProc.CPA_X_POSITIONS[0]+10,900, "#color[0]{CPA}");
      FormatLatex(evtProc.APA_X_POSITIONS[1]+10,900, "#color[0]{APA}");
      FormatLatex(evtProc.CPA_X_POSITIONS[1]+10,900, "#color[0]{CPA}");
      FormatLatex(evtProc.APA_X_POSITIONS[2]+10,900, "#color[0]{APA}");
    }
    else{
      FormatLatex(evtProc.APA_X_POSITIONS[0]+10,450, "#color[0]{APA}");
      FormatLatex(evtProc.CPA_X_POSITIONS[0]+10,450, "#color[0]{CPA}");
      FormatLatex(evtProc.APA_X_POSITIONS[1]+10,450, "#color[0]{APA}");
      FormatLatex(evtProc.CPA_X_POSITIONS[1]+10,450, "#color[0]{CPA}");
      FormatLatex(evtProc.APA_X_POSITIONS[2]+10,450, "#color[0]{APA}");
    }

    c4->SaveAs((location+"/charge_vs_X_plane"+std::to_string(iWire)+tag+".png").c_str());
    c4->SaveAs((location+"/charge_vs_X_plane"+std::to_string(iWire)+tag+".root").c_str());
    c4->Clear();
    
    SetHistogramStyle2D(h_corr_charges.at(iWire),"x [cm]", " Charge deposition [ADC/cm]",false);
    h_corr_charges.at(iWire)->GetZaxis()->SetLabelSize(0.03);
    h_corr_charges.at(iWire)->GetZaxis()->SetLabelFont(132);
    h_corr_charges.at(iWire)->Draw("colz");

    // Draw the APA and CPA lines and labels
    for(unsigned int iLine = 0; iLine < APACPALines.size(); ++iLine){
      APACPALines.at(iLine)->Draw();
    }

    if(iWire == 2){
      FormatLatex(evtProc.APA_X_POSITIONS[0]+10,900, "#color[0]{APA}");
      FormatLatex(evtProc.CPA_X_POSITIONS[0]+10,900, "#color[0]{CPA}");
      FormatLatex(evtProc.APA_X_POSITIONS[1]+10,900, "#color[0]{APA}");
      FormatLatex(evtProc.CPA_X_POSITIONS[1]+10,900, "#color[0]{CPA}");
      FormatLatex(evtProc.APA_X_POSITIONS[2]+10,900, "#color[0]{APA}");
    }
    else{
      FormatLatex(evtProc.APA_X_POSITIONS[0]+10,450, "#color[0]{APA}");
      FormatLatex(evtProc.CPA_X_POSITIONS[0]+10,450, "#color[0]{CPA}");
      FormatLatex(evtProc.APA_X_POSITIONS[1]+10,450, "#color[0]{APA}");
      FormatLatex(evtProc.CPA_X_POSITIONS[1]+10,450, "#color[0]{CPA}");
      FormatLatex(evtProc.APA_X_POSITIONS[2]+10,450, "#color[0]{APA}");
    }

    c4->SaveAs((location+"/corr_charge_vs_X_plane"+std::to_string(iWire)+tag+".png").c_str());
    c4->SaveAs((location+"/corr_charge_vs_X_plane"+std::to_string(iWire)+tag+".root").c_str());
    c4->Clear();
    
    // # Electrons
    SetHistogramStyle2D(h_nelecs.at(iWire),"x [cm]", " Number of electrons",false);
    h_nelecs.at(iWire)->GetZaxis()->SetLabelSize(0.03);
    h_nelecs.at(iWire)->GetZaxis()->SetLabelFont(132);
    h_nelecs.at(iWire)->Draw("colz");

    // Draw the APA and CPA lines and labels
    for(unsigned int iLine = 0; iLine < APACPALines.size(); ++iLine){
      if(iWire == 2)
        APACPALines.at(iLine)->SetY2(2.5e5);
      else
        APACPALines.at(iLine)->SetY2(1.5e5);
      APACPALines.at(iLine)->Draw();
    }

    if(iWire == 2){
      FormatLatex(evtProc.APA_X_POSITIONS[0]+10,2.25e5, "#color[0]{APA}");
      FormatLatex(evtProc.CPA_X_POSITIONS[0]+10,2.25e5, "#color[0]{CPA}");
      FormatLatex(evtProc.APA_X_POSITIONS[1]+10,2.25e5, "#color[0]{APA}");
      FormatLatex(evtProc.CPA_X_POSITIONS[1]+10,2.25e5, "#color[0]{CPA}");
      FormatLatex(evtProc.APA_X_POSITIONS[2]+10,2.25e5, "#color[0]{APA}");
    }
    else{
      FormatLatex(evtProc.APA_X_POSITIONS[0]+10,1.35e5, "#color[0]{APA}");
      FormatLatex(evtProc.CPA_X_POSITIONS[0]+10,1.35e5, "#color[0]{CPA}");
      FormatLatex(evtProc.APA_X_POSITIONS[1]+10,1.35e5, "#color[0]{APA}");
      FormatLatex(evtProc.CPA_X_POSITIONS[1]+10,1.35e5, "#color[0]{CPA}");
      FormatLatex(evtProc.APA_X_POSITIONS[2]+10,1.35e5, "#color[0]{APA}");
    }

    c4->SaveAs((location+"/nelecs_vs_X_plane"+std::to_string(iWire)+tag+".png").c_str());
    c4->SaveAs((location+"/nelecs_vs_X_plane"+std::to_string(iWire)+tag+".root").c_str());
    c4->Clear();
    
    SetHistogramStyle2D(h_corr_nelecs.at(iWire),"x [cm]", " Number of electrons",false);
    h_corr_nelecs.at(iWire)->GetZaxis()->SetLabelSize(0.03);
    h_corr_nelecs.at(iWire)->GetZaxis()->SetLabelFont(132);
    h_corr_nelecs.at(iWire)->Draw("colz");

    // Draw the APA and CPA lines and labels
    for(unsigned int iLine = 0; iLine < APACPALines.size(); ++iLine){
      APACPALines.at(iLine)->Draw();
    }

    if(iWire == 2){
      FormatLatex(evtProc.APA_X_POSITIONS[0]+10,2.25e5, "#color[0]{APA}");
      FormatLatex(evtProc.CPA_X_POSITIONS[0]+10,2.25e5, "#color[0]{CPA}");
      FormatLatex(evtProc.APA_X_POSITIONS[1]+10,2.25e5, "#color[0]{APA}");
      FormatLatex(evtProc.CPA_X_POSITIONS[1]+10,2.25e5, "#color[0]{CPA}");
      FormatLatex(evtProc.APA_X_POSITIONS[2]+10,2.25e5, "#color[0]{APA}");
    }
    else{
      FormatLatex(evtProc.APA_X_POSITIONS[0]+10,1.35e5, "#color[0]{APA}");
      FormatLatex(evtProc.CPA_X_POSITIONS[0]+10,1.35e5, "#color[0]{CPA}");
      FormatLatex(evtProc.APA_X_POSITIONS[1]+10,1.35e5, "#color[0]{APA}");
      FormatLatex(evtProc.CPA_X_POSITIONS[1]+10,1.35e5, "#color[0]{CPA}");
      FormatLatex(evtProc.APA_X_POSITIONS[2]+10,1.35e5, "#color[0]{APA}");
    }

    c4->SaveAs((location+"/corr_nelecs_vs_X_plane"+std::to_string(iWire)+tag+".png").c_str());
    c4->SaveAs((location+"/corr_nelecs_vs_X_plane"+std::to_string(iWire)+tag+".root").c_str());
    c4->Clear();
    
    SetHistogramStyle2D(h_corr_nelec_corr_charges.at(iWire),"x [cm]", " e^{-} per charge [N_{e^{-}}/ADC]",false);
    h_corr_nelec_corr_charges.at(iWire)->GetZaxis()->SetLabelSize(0.03);
    h_corr_nelec_corr_charges.at(iWire)->GetZaxis()->SetLabelFont(132);
    h_corr_nelec_corr_charges.at(iWire)->Draw("colz");

    // Draw the APA and CPA lines and labels
    for(unsigned int iLine = 0; iLine < APACPALines.size(); ++iLine){
      APACPALines.at(iLine)->SetY1(160);
      APACPALines.at(iLine)->SetY2(300);
      APACPALines.at(iLine)->Draw();
    }

    FormatLatex(evtProc.APA_X_POSITIONS[0]+10,285, "#color[0]{APA}");
    FormatLatex(evtProc.CPA_X_POSITIONS[0]+10,285, "#color[0]{CPA}");
    FormatLatex(evtProc.APA_X_POSITIONS[1]+10,285, "#color[0]{APA}");
    FormatLatex(evtProc.CPA_X_POSITIONS[1]+10,285, "#color[0]{CPA}");
    FormatLatex(evtProc.APA_X_POSITIONS[2]+10,285, "#color[0]{APA}");

    c4->SaveAs((location+"/corr_nelecs_corr_charge_vs_X_plane"+std::to_string(iWire)+tag+".png").c_str());
    c4->SaveAs((location+"/corr_nelecs_corr_charge_vs_X_plane"+std::to_string(iWire)+tag+".root").c_str());
    c4->Clear();
    
    SetHistogramStyle2D(h_eDep_nDaught.at(iWire),"Muon daughters", "Energy deposition per unit length [MeV/cm]",false);
    h_eDep_nDaught.at(iWire)->GetZaxis()->SetLabelSize(0.03);
    h_eDep_nDaught.at(iWire)->GetZaxis()->SetLabelFont(132);
    h_eDep_nDaught.at(iWire)->Draw("colz");
    c4->SaveAs((location+"/truth_tracks_eDep_vs_nDaughters"+std::to_string(iWire)+tag+".png").c_str());
    c4->SaveAs((location+"/truth_tracks_eDep_vs_nDaughters"+std::to_string(iWire)+tag+".root").c_str());
    c4->Clear();

    SetHistogramStyle2D(h_eDep_E.at(iWire),"Muon energy [GeV]", "Energy deposition per unit length [MeV/cm]",false);
    h_eDep_E.at(iWire)->GetZaxis()->SetLabelSize(0.03);
    h_eDep_E.at(iWire)->GetZaxis()->SetLabelFont(132);
    h_eDep_E.at(iWire)->Draw("colz");
    c4->SaveAs((location+"/truth_tracks_eDep_vs_E"+std::to_string(iWire)+tag+".png").c_str());
    c4->SaveAs((location+"/truth_tracks_eDep_vs_E"+std::to_string(iWire)+tag+".root").c_str());
    c4->Clear();

  } // Wire planes
  
  SetHistogramStyle2D(h_E_nDaught,"Muon energy [GeV]", "Number of daughters",false);
  h_E_nDaught->GetZaxis()->SetLabelSize(0.03);
  h_E_nDaught->GetZaxis()->SetLabelFont(132);
  h_E_nDaught->Draw("colz");
  c4->SaveAs((location+"/truth_tracks_nDaught_vs_E"+tag+".png").c_str());
  c4->SaveAs((location+"/truth_tracks_nDaught_vs_E"+tag+".root").c_str());
  c4->Clear();

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
