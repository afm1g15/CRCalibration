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
   "taulife",
   "geant_list_size",
   "inTPCActive",
   "TrackId",
   "pdg",
   "origin",
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
   "StartE_tpcAV",
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
  int n   = -1;
  int cry = 0; // Whether or not to save the cosmic (cry) information (trkorigin = 2) or the beam information (trkorigin = 4)
  int drawMCParticles = 0; // Should we draw the particles or not?

  double tau = 3.5; // measured electron lifetime

  std::string input_list = "";
  std::string location="";
  std::string tag="";
  std::vector<double> minx_fid, miny_fid, minz_fid;
  std::vector<double> maxx_fid, maxy_fid, maxz_fid;
  std::vector<double> minx_av, miny_av, minz_av;
  std::vector<double> maxx_av, maxy_av, maxz_av;

  p->getValue("InputList",       input_list);
  p->getValue("Location",        location);
  p->getValue("DrawMCParticles", drawMCParticles);
  p->getValue("Tag",             tag);
  p->getValue("NFiles",          n);
  p->getValue("Cry",             cry);
  p->getValue("Tau",             tau);
  p->getValue("MinXFid",         minx_fid);
  p->getValue("MinYFid",         miny_fid);
  p->getValue("MinZFid",         minz_fid);
  p->getValue("MaxXFid",         maxx_fid);
  p->getValue("MaxYFid",         maxy_fid);
  p->getValue("MaxZFid",         maxz_fid);
  p->getValue("MinXAV",          minx_av);
  p->getValue("MinYAV",          miny_av);
  p->getValue("MinZAV",          minz_av);
  p->getValue("MaxXAV",          maxx_av);
  p->getValue("MaxYAV",          maxy_av);
  p->getValue("MaxZAV",          maxz_av);

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
  TH1D *h_startX     = new TH1D("h_startX","",80,-400,400);  // Start X position of the muons
  TH1D *h_startY     = new TH1D("h_startY","",61,0,610);  // Start Y position of the muons
  TH1D *h_startZ     = new TH1D("h_startZ","",70,0,700);  // Start Z position of the muons
  TH1D *h_endX       = new TH1D("h_endX","",60,-400,400);    // end X position of the muons
  TH1D *h_endY       = new TH1D("h_endY","",60,0,610);    // end Y position of the muons
  TH1D *h_endZ       = new TH1D("h_endZ","",60,-1,700);    // end Z position of the muons
  TH1D *h_length     = new TH1D("h_length","",50,0,800);   // Length of the muons
  TH1D *h_nDaughters = new TH1D("h_nDaughters","",50,0,200); // Number of muon daughters
  TH1D *h_mom        = new TH1D("h_mom","",50,0,40);        // Momentum of the muons [GeV]
  TH1D *h_energy     = new TH1D("h_energy","",50,0,40);     // Start TPC Energy of the muons [GeV]
  TH1D *h_energyLog  = new TH1D("h_energyLog","",50,0.8,40);     // Start TPC Energy of the muons [GeV]
  TH1D *h_energyPlug = new TH1D("h_energyPlug","",50,0.8,40);     // Start TPC Energy of the muons [GeV]
  TH1D *h_energyHalo = new TH1D("h_energyHalo","",50,0.8,40);     // Start TPC Energy of the muons [GeV]
  TH1D *h_eng        = new TH1D("h_eng","",50,0,40);     // Generated Energy of the muons [GeV]
  SetLogX(h_energyLog);
  SetLogX(h_energyHalo);
  SetLogX(h_energyPlug);
 
  // 2D truth studies
  TH2D *h_startX_energy = new TH2D("h_startX_energy","",60,-200,400,47,0.8,40); // Energy vs start position
  TH2D *h_startY_energy = new TH2D("h_startY_energy","",60,10,610,47,0.8,40); // Energy vs start position
  SetLogY(h_startX_energy);
  SetLogY(h_startY_energy);

  // Angular distributions
  TH1D *h_thetaxz    = new TH1D("h_thetaxz","",50,-180,180); // ThetaXZ
  TH1D *h_thetayz    = new TH1D("h_thetayz","",50,-180,0);   // ThetaYZ
  TH1D *h_theta      = new TH1D("h_theta","",50,0,180);      // Theta
  TH1D *h_phi        = new TH1D("h_phi","",50,-180,0);     // Phi
  TH1D *h_eDepPerL_0 = new TH1D("h_eDepPerL_0","",50,0,40); // Energy deposition per unit length
  TH1D *h_eDepPerL_1 = new TH1D("h_eDepPerL_1","",50,0,40); // Energy deposition per unit length
  TH1D *h_eDepPerL_2 = new TH1D("h_eDepPerL_2","",50,0,40); // Energy deposition per unit length

  // Truth-level plane quantities
  TH1D *h_plane_cross  = new TH1D("h_plane_cross","",7,0,7); // Number of tracks crossing each plane
  TH1D *h_n_crossed    = new TH1D("h_n_crossed","",7,0,7); // Number of planes crossed by each track
  TH1D *h_plane_enter  = new TH1D("h_plane_enter","",7,0,7); // Number of tracks entering from each external plane
  TH1D *h_plane_exit   = new TH1D("h_plane_exit","",7,0,7); // Number of tracks exiting from each external plane
  TH1D *h_enter_dist   = new TH1D("h_enter_dist","",50,0,1); // Distance from candidate entrance plane 
  TH1D *h_exit_dist    = new TH1D("h_exit_dist","",50,0,1); // Distance from candidate exit plane

  // Hit-level depositions
  TH2D *h_hit_energy_x_0 = new TH2D("h_hit_energy_x_plane0","",50,-400,400,50,0,50); // Hit deposition energy vs X primary TPC muons, plane0
  TH2D *h_hit_charge_x_0 = new TH2D("h_hit_charge_x_plane0","",50,-400,400,50,0,80); // Hit deposition energy vs X primary TPC muons, plane0
  TH2D *h_hit_nelec_x_0  = new TH2D("h_hit_nelec_x_plane0","",50,-400,400,50,0,1e6); // Hit deposition energy vs X primary TPC muons, plane0
  TH2D *h_hit_energy_x_1 = new TH2D("h_hit_energy_x_plane1","",50,-400,400,50,0,50); // Hit deposition energy vs X primary TPC muons, plane1
  TH2D *h_hit_charge_x_1 = new TH2D("h_hit_charge_x_plane1","",50,-400,400,50,0,80); // Hit deposition energy vs X primary TPC muons, plane1
  TH2D *h_hit_nelec_x_1  = new TH2D("h_hit_nelec_x_plane1","",50,-400,400,50,0,1e6); // Hit deposition energy vs X primary TPC muons, plane1
  TH2D *h_hit_energy_x_2 = new TH2D("h_hit_energy_x_plane2","",50,-400,400,50,0,50); // Hit deposition energy vs X primary TPC muons, plane2
  TH2D *h_hit_charge_x_2 = new TH2D("h_hit_charge_x_plane2","",50,-400,400,50,0,80); // Hit deposition energy vs X primary TPC muons, plane2
  TH2D *h_hit_nelec_x_2  = new TH2D("h_hit_nelec_x_plane2","",50,-400,400,50,0,1e6); // Hit deposition energy vs X primary TPC muons, plane2
  TH2D *h_corr_hit_energy_x_0 = new TH2D("h_corr_hit_energy_x_plane0","",50,-400,400,50,0,50); // Hit deposition energy vs X primary TPC muons, plane0
  TH2D *h_corr_hit_charge_x_0 = new TH2D("h_corr_hit_charge_x_plane0","",50,-400,400,50,0,750); // Hit deposition energy vs X primary TPC muons, plane0
  TH2D *h_corr_hit_nelec_x_0  = new TH2D("h_corr_hit_nelec_x_plane0","",50,-400,400,50,0,1e6); // Hit deposition energy vs X primary TPC muons, plane0
  TH2D *h_corr_hit_energy_x_1 = new TH2D("h_corr_hit_energy_x_plane1","",50,-400,400,50,0,50); // Hit deposition energy vs X primary TPC muons, plane1
  TH2D *h_corr_hit_charge_x_1 = new TH2D("h_corr_hit_charge_x_plane1","",50,-400,400,50,0,750); // Hit deposition energy vs X primary TPC muons, plane1
  TH2D *h_corr_hit_nelec_x_1  = new TH2D("h_corr_hit_nelec_x_plane1","",50,-400,400,50,0,1e6); // Hit deposition energy vs X primary TPC muons, plane1
  TH2D *h_corr_hit_energy_x_2 = new TH2D("h_corr_hit_energy_x_plane2","",50,-400,400,50,0,50); // Hit deposition energy vs X primary TPC muons, plane2
  TH2D *h_corr_hit_charge_x_2 = new TH2D("h_corr_hit_charge_x_plane2","",50,-400,400,50,0,1000); // Hit deposition energy vs X primary TPC muons, plane2
  TH2D *h_corr_hit_nelec_x_2  = new TH2D("h_corr_hit_nelec_x_plane2","",50,-400,400,50,0,1e6); // Hit deposition energy vs X primary TPC muons, plane2
  TH2D *h_corr_hit_nelec_energy_x_0  = new TH2D("h_corr_hit_nelec_energy_x_plane0","",50,-400,400,50,0,1e-4); // Energy per electron [MeV/el] plane 0 
  TH2D *h_corr_hit_nelec_energy_x_1  = new TH2D("h_corr_hit_nelec_energy_x_plane1","",50,-400,400,50,0,1e-4); // Energy per electron [MeV/el] plane 1
  TH2D *h_corr_hit_nelec_energy_x_2  = new TH2D("h_corr_hit_nelec_energy_x_plane2","",50,-400,400,50,0,1e-4); // Energy per electron [MeV/el] plane 2
  TH2D *h_corr_hit_charge_energy_x_0  = new TH2D("h_corr_hit_charge_energy_x_plane0","",50,-400,400,50,0,1); // Energy per charge [MeV/ADC] plane 0 
  TH2D *h_corr_hit_charge_energy_x_1  = new TH2D("h_corr_hit_charge_energy_x_plane1","",50,-400,400,50,0,1); // Energy per charge [MeV/ADC] plane 1
  TH2D *h_corr_hit_charge_energy_x_2  = new TH2D("h_corr_hit_charge_energy_x_plane2","",50,-400,400,50,0,1); // Energy per charge [MeV/ADC] plane 2
  TH2D *h_corr_hit_nelec_corr_energy_x_0  = new TH2D("h_corr_hit_nelec_corr_energy_x_plane0","",50,-400,400,50,0,1e-4); // Energy per electron [MeV/el] plane 0 
  TH2D *h_corr_hit_nelec_corr_energy_x_1  = new TH2D("h_corr_hit_nelec_corr_energy_x_plane1","",50,-400,400,50,0,1e-4); // Energy per electron [MeV/el] plane 1
  TH2D *h_corr_hit_nelec_corr_energy_x_2  = new TH2D("h_corr_hit_nelec_corr_energy_x_plane2","",50,-400,400,50,0,1e-4); // Energy per electron [MeV/el] plane 2
  TH2D *h_corr_hit_charge_corr_energy_x_0  = new TH2D("h_corr_hit_charge_corr_energy_x_plane0","",50,-400,400,50,0,1); // Energy per charge [MeV/ADC] plane 0 
  TH2D *h_corr_hit_charge_corr_energy_x_1  = new TH2D("h_corr_hit_charge_corr_energy_x_plane1","",50,-400,400,50,0,1); // Energy per charge [MeV/ADC] plane 1
  TH2D *h_corr_hit_charge_corr_energy_x_2  = new TH2D("h_corr_hit_charge_corr_energy_x_plane2","",50,-400,400,50,0,1); // Energy per charge [MeV/ADC] plane 2
  TH2D *h_corr_hit_nelec_corr_charge_x_0  = new TH2D("h_corr_hit_nelec_corr_charge_x_plane0","",50,-400,400,50,0,30000); // charge per electron [MeV/el] plane 0 
  TH2D *h_corr_hit_nelec_corr_charge_x_1  = new TH2D("h_corr_hit_nelec_corr_charge_x_plane1","",50,-400,400,50,0,20000); // charge per electron [MeV/el] plane 1
  TH2D *h_corr_hit_nelec_corr_charge_x_2  = new TH2D("h_corr_hit_nelec_corr_charge_x_plane2","",50,-400,400,50,0,20000); // charge per electron [MeV/el] plane 2
  TH2D *h_eDep_E_0 = new TH2D("h_eDep_E_0","",50,0,40,50,0,40); // Number of muon daughters vs energy depositions per unit length
  TH2D *h_eDep_E_1 = new TH2D("h_eDep_E_1","",50,0,40,50,0,40); // Number of muon daughters vs energy depositions per unit length
  TH2D *h_eDep_E_2 = new TH2D("h_eDep_E_2","",50,0,40,50,0,40); // Number of muon daughters vs energy depositions per unit length

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
  std::vector<TH2D*> h_eDep_E{h_eDep_E_0,h_eDep_E_1,h_eDep_E_2};
 
  // Setup vectors of MCParticles
  std::vector<TVector3> startPoints, endPoints;

  // Setup counters
  unsigned int nMu        = 0;
  unsigned int nPrimaryMu = 0;
  unsigned int nLongMu    = 0;
  
  unsigned int noPlanes     = 0;
  unsigned int topBottom    = 0;
  unsigned int topOrBottom  = 0;
  unsigned int min2APACPA   = 0;
  unsigned int eq1APACPA    = 0;
  unsigned int frontOrBack  = 0;
  unsigned int frontBack    = 0;
  unsigned int stopping     = 0;
  unsigned int exiting      = 0;
  unsigned int extCrossed2  = 0;

  unsigned int nThru        = 0;
  
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

      // Check that the origin of the particle matches our choice
      // origin = 2 (Cry, cosmic) origin = 4 (NuWro, beam)
      if((evt->origin[iG4] != 2 && cry) || (evt->origin[iG4] != 4 && !cry)) continue;

      int pdg        = evt->pdg[iG4];
      int id         = evt->TrackId[iG4];

      // Check the particle enters the TPC volume
      if(!evt->inTPCActive[iG4] || evt->StartE_tpcAV[iG4] < 0) continue;

      // Make sure we are looking at a muon
      if(abs(pdg) != 13) continue;
      nMu++;

      // Make sure we are looking at the daughter of a primary particle
      if(evt->Mother[iG4] != 0 || cry) continue;
      nPrimaryMu++;

      // For the entrance tests
      // If the start or end locations are outside the detector, set them to be at the edge of the detector
      // Only do this for 1 coordinate
      TVector3 vtx(evt->StartPointx[iG4],evt->StartPointy[iG4],evt->StartPointz[iG4]);
      TVector3 end(evt->EndPointx[iG4],evt->EndPointy[iG4],evt->EndPointz[iG4]);
      TVector3 vtxAV(evt->StartPointx_tpcAV[iG4],evt->StartPointy_tpcAV[iG4],evt->StartPointz_tpcAV[iG4]);
      TVector3 endAV(evt->EndPointx_tpcAV[iG4],evt->EndPointy_tpcAV[iG4],evt->EndPointz_tpcAV[iG4]);
      
      float lengthAV = (endAV-vtxAV).Mag();
      float energy   = evt->StartE_tpcAV[iG4];
      float eng      = evt->Eng[iG4];

      // Determine if the track enters at the top and leaves through the bottom
      float vtxDy = abs(vtxAV.Y()-evt->StartPointy[iG4]);
      float endDy = abs(endAV.Y()-evt->EndPointy[iG4]);

      // Check for number of through-going and stopping muons in truth
      bool throughGoing = IsTrueThroughGoing(vtx,end,vtxAV,endAV);
      bool isStopping   = IsTrueStopping(vtx,end,vtxAV,endAV);

      // Fill the truth-track quantities 
      h_startX->Fill(vtxAV.X());
      h_startY->Fill(vtxAV.Y());
      h_startZ->Fill(vtxAV.Z());
      h_endX->Fill(endAV.X());
      h_endY->Fill(endAV.Y());
      h_endZ->Fill(endAV.Z());
      h_length->Fill(lengthAV);
      h_mom->Fill(evt->P[iG4]);
      h_energy->Fill(energy);
      h_energyLog->Fill(energy);
      h_eng->Fill(eng);
      h_nDaughters->Fill(evt->NumberDaughters[iG4]);
      h_phi->Fill(evt->phi[iG4]*(180/TMath::Pi()));
      h_theta->Fill(evt->theta[iG4]*(180/TMath::Pi()));
      h_thetaxz->Fill(evt->theta_xz[iG4]*(180/TMath::Pi()));
      h_thetayz->Fill(evt->theta_yz[iG4]*(180/TMath::Pi()));

      // Fill energies for the plug location and the halo location
      if(vtxAV.X() > -100 && vtxAV.X() < 50 && vtxAV.Y() > 380 && vtxAV.Y() < 480)
        h_energyPlug->Fill(energy);
      if(vtxAV.X() > 100 && vtxAV.Y() > 500)
        h_energyHalo->Fill(energy);

      // Check for through going and fill 
      if(throughGoing){
        nThru++;
      }
      h_startX_energy->Fill(vtxAV.X(),energy);
      h_startY_energy->Fill(vtxAV.Y(),energy);
      
      // Fill the start and end point vectors for drawing the mcparticles
      startPoints.push_back(vtxAV);
      endPoints.push_back(endAV);

      // Find the closest plane to the start vertex and count it as a crossing plane
      Plane enteringPlane = GetClosestPlane(extPlanes, vtxAV, endAV);
      double distFromEntrance = GetDistanceToPlane(enteringPlane, vtxAV, endAV);
      h_enter_dist->Fill(distFromEntrance*1000.); // mm

      // Find the closest plane to the end vertex and count it as a crossing plane
      Plane exitingPlane = GetClosestPlane(extPlanes, endAV, vtxAV);
      double distFromExit = GetDistanceToPlane(exitingPlane, endAV, vtxAV);
      h_exit_dist->Fill(distFromExit*1000.); // mm

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
      if(lengthAV < 200) continue;
      nLongMu++;

      if(nExtCrossed >= 2) extCrossed2++;
      if(throughGoing) exiting++;
      if(isStopping) stopping++;
      
      // Loop over the labels crossed by this track and fill the appropriate counters
      bool APA = false;
      bool CPA = false;
      bool top = false;
      bool bot = false;
      bool fro = false;
      bool bac = false;

      unsigned int apaCPA = 0;
      for(std::string &str : labelsCrossed){
        // Get the converted label
        std::string longLab = planeLabels.find(str)->second;
        TString lab(longLab);
        if(lab.Contains("APA") || lab.Contains("CPA")){
          apaCPA++;
        }
        else if(lab.Contains("Top"))
          top = true;
        else if(lab.Contains("Bot"))
          bot = true;
        else if(lab.Contains("Fro"))
          fro = true;
        else if(lab.Contains("Bac"))
          bac = true;
      }
      if(apaCPA == 1)
        eq1APACPA++;
      if(apaCPA > 1)
        min2APACPA++;
      if(top || bot)
        topOrBottom++;
      if(top && bot)
        topBottom++;
      if(fro || bac)
        frontOrBack++;
      if(fro && bac)
        frontBack++;

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
          float hitQ      = evt->hit_charge[iHit]*5.; // Scale to FD ADC definition: 1 ADC = 200 e vs PD: 1 ADC = 1000 e
          float hitEl     = evt->hit_nelec[iHit];

          totalEDep += hitE;
          totalQDep += hitQ;

          // Check if x is lower than the APA bound, charge seems to accumulate there
          if(hitX < evtProc.PD_APA_X_POSITIONS[0] || hitX > evtProc.PD_APA_X_POSITIONS[1]) continue;

          // Lifetime correction
          int tpc      = evtProc.WhichTPC(hitX) + 1;
          float dx     = ( -1 + 2*(tpc%2) )*(hitX - evtProc.PD_APA_X_POSITIONS[tpc/2]);
          float dt     = dx*evtProc.kXtoT;
          float simTau = evt->taulife/1.0e4;

          float corr      = TMath::Exp(-dt/tau);
          float nCorr     = TMath::Exp(-dt/simTau);
          float eCorr     = TMath::Exp(-dt/tau) / TMath::Exp(-dt/simTau); // Correct for the already-corrected energy
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
        float totalEDepPerLength = totalEDep/(lengthAV);
        if(totalEDepPerLength < std::numeric_limits<float>::epsilon()) continue;
        h_eDep_E.at(iWire)->Fill(evt->StartE_tpcAV[iG4],totalEDepPerLength);
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
    "Long TPC $\\mu$",
    "Crosses top or bottom",
    "Crosses top and bottom",
    "Crosses front or back",
    "Crosses front and back",
    "Crosses 1 APA or CPA",
    "Crosses $ \\geq $ 2 APA or CPA",
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
    frontOrBack,
    frontBack,
    eq1APACPA,
    min2APACPA,
    stopping,
    exiting
  };

  ofstream texFile;
  texFile.open(location+"/truth_contents"+tag+".tex");
  WriteStatsToTeX(texFile, n, contents, rates, static_cast<double>(nLongMu), "Long TPC $\\mu$");

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
  SetCanvasStyle(c2, 0.12,0.08,0.06,0.12,0,0,0);

  l->SetNColumns(2);
  l->SetX1NDC(0.47);

  SetHistogramStyle1D(h_startX,"Muon X [cm]", "Rate");
  h_startX->Draw("hist");
  h_endX->Draw("same");
  h_startX->SetLineWidth(3);
  h_startX->SetLineStyle(2);
  h_endX->SetLineWidth(3);
  h_endX->SetLineStyle(2);
  h_startX->SetLineColor(pal.at(1));
  h_endX->SetLineColor(pal.at(0));
  h_startX->GetYaxis()->SetTitleOffset(0.95);

  l->AddEntry(h_startX, "Start", "l");
  l->AddEntry(h_endX, "End", "l");
  l->Draw();
  
  // Draw the location of the beam pipe
  double lineMinX = std::min(h_startX->GetMinimum(),h_endX->GetMinimum());
  double lineMaxX = std::max(h_startX->GetMaximum(),h_endX->GetMaximum());

  TLine *lX = new TLine(beamStart.X(), lineMinX, beamStart.X(), lineMaxX*1.05);
  lX->SetLineColor(pal.at(1));
  lX->SetLineWidth(3);
  lX->SetLineStyle(7);
  lX->Draw();

  TLine *lXE = new TLine(beamEnd.X(), lineMinX, beamEnd.X(), lineMaxX*1.05);
  lXE->SetLineColor(pal.at(0));
  lXE->SetLineWidth(3);
  lXE->SetLineStyle(7);
  lXE->Draw();

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
  h_startY->SetLineColor(pal.at(1));
  h_endY->SetLineColor(pal.at(0));
  h_startY->GetYaxis()->SetTitleOffset(0.95);

  l->AddEntry(h_startY, "Start", "l");
  l->AddEntry(h_endY, "End", "l");
  l->Draw();
  
  // Draw the location of the beam pipe
  double lineMinY = std::min(h_startY->GetMinimum(),h_endY->GetMinimum());
  double lineMaxY = std::max(h_startY->GetMaximum(),h_endY->GetMaximum());

  TLine *lY = new TLine(beamStart.Y(), lineMinY, beamStart.Y(), lineMaxY*1.05);
  lY->SetLineColor(pal.at(1));
  lY->SetLineWidth(3);
  lY->SetLineStyle(7);
  lY->Draw();

  TLine *lYE = new TLine(beamEnd.Y(), lineMinY, beamEnd.Y(), lineMaxY*1.05);
  lYE->SetLineColor(pal.at(0));
  lYE->SetLineWidth(3);
  lYE->SetLineStyle(7);
  lYE->Draw();

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
  h_startZ->SetLineColor(pal.at(1));
  h_endZ->SetLineColor(pal.at(0));
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
  h_n_crossed->SetLineColor(pal.at(1));
  h_n_crossed->GetYaxis()->SetTitleOffset(0.95);
  c3->SaveAs((location+"/truth_tracks_crossed_nplanes"+tag+".png").c_str());
  c3->SaveAs((location+"/truth_tracks_crossed_nplanes"+tag+".root").c_str());
  c3->Clear();

  SetHistogramStyle1D(h_length,"Muon length [cm]", "Rate");
  h_length->Draw("hist");
  h_length->SetLineWidth(3);
  h_length->SetLineColor(pal.at(1));
  h_length->GetYaxis()->SetTitleOffset(0.95);
  c3->SaveAs((location+"/truth_tracks_length"+tag+".png").c_str());
  c3->SaveAs((location+"/truth_tracks_length"+tag+".root").c_str());
  c3->Clear();

  SetHistogramStyle1D(h_nDaughters,"Muon daughters", "Rate");
  h_nDaughters->Draw("hist");
  h_nDaughters->SetLineWidth(3);
  h_nDaughters->SetLineColor(pal.at(1));
  h_nDaughters->GetYaxis()->SetTitleOffset(0.95);
  c3->SaveAs((location+"/truth_tracks_nDaughters"+tag+".png").c_str());
  c3->SaveAs((location+"/truth_tracks_nDaughters"+tag+".root").c_str());
  c3->Clear();

  SetHistogramStyle1D(h_thetaxz,"Muon #theta_{XZ} [^{#circ}]", "Rate");
  h_thetaxz->Draw("hist");
  h_thetaxz->SetLineWidth(3);
  h_thetaxz->SetLineColor(pal.at(1));
  h_thetaxz->GetYaxis()->SetTitleOffset(0.95);
  c3->SaveAs((location+"/truth_tracks_thetaxz"+tag+".png").c_str());
  c3->SaveAs((location+"/truth_tracks_thetaxz"+tag+".root").c_str());
  c3->Clear();

  SetHistogramStyle1D(h_thetayz,"Muon #theta_{YZ} [^{#circ}]", "Rate");
  h_thetayz->Draw("hist");
  h_thetayz->SetLineWidth(3);
  h_thetayz->SetLineColor(pal.at(1));
  h_thetayz->GetYaxis()->SetTitleOffset(0.95);
  c3->SaveAs((location+"/truth_tracks_thetayz"+tag+".png").c_str());
  c3->SaveAs((location+"/truth_tracks_thetayz"+tag+".root").c_str());
  c3->Clear();

  SetHistogramStyle1D(h_theta,"Muon #theta [^{#circ}]", "Rate");
  h_theta->Draw("hist");
  h_theta->SetLineWidth(3);
  h_theta->SetLineColor(pal.at(1));
  h_theta->GetYaxis()->SetTitleOffset(0.95);
  c3->SaveAs((location+"/truth_tracks_theta"+tag+".png").c_str());
  c3->SaveAs((location+"/truth_tracks_theta"+tag+".root").c_str());
  c3->Clear();

  SetHistogramStyle1D(h_phi,"Muon #phi [^{#circ}]", "Rate");
  h_phi->Draw("hist");
  h_phi->SetLineWidth(3);
  h_phi->SetLineColor(pal.at(1));
  h_phi->GetYaxis()->SetTitleOffset(0.95);
  c3->SaveAs((location+"/truth_tracks_phi"+tag+".png").c_str());
  c3->SaveAs((location+"/truth_tracks_phi"+tag+".root").c_str());
  c3->Clear();

  for(unsigned int iWire = 0; iWire < 3; ++iWire){
    SetHistogramStyle1D(h_eDepPerL.at(iWire),"Energy deposition per unit length [MeV/cm]", "Rate");
    h_eDepPerL.at(iWire)->Draw("hist");
    h_eDepPerL.at(iWire)->SetLineWidth(3);
    h_eDepPerL.at(iWire)->SetLineColor(pal.at(1));
    h_eDepPerL.at(iWire)->GetYaxis()->SetTitleOffset(0.95);
    c3->SaveAs((location+"/truth_tracks_eDep_perL"+tag+".png").c_str());
    c3->SaveAs((location+"/truth_tracks_eDep_perL"+tag+".root").c_str());
    c3->Clear();
  }

  //c3->SetLogx();
  SetHistogramStyle1D(h_mom,"Muon momentum [GeV]", "Rate");
  h_mom->Draw("hist");
  h_mom->SetLineWidth(3);
  h_mom->SetLineColor(pal.at(1));
  h_mom->GetYaxis()->SetTitleOffset(0.95);
  c3->SaveAs((location+"/truth_tracks_mom"+tag+".png").c_str());
  c3->SaveAs((location+"/truth_tracks_mom"+tag+".root").c_str());
  c3->Clear();

  SetHistogramStyle1D(h_energy,"Muon energy [GeV]", "Rate");
  h_energy->Draw("hist");
  h_energy->SetLineWidth(3);
  h_energy->SetLineStyle(2);
  h_energy->SetLineColor(pal.at(1));
  h_energy->GetYaxis()->SetTitleOffset(0.95);
  c3->SaveAs((location+"/truth_tracks_energy"+tag+".png").c_str());
  c3->SaveAs((location+"/truth_tracks_energy"+tag+".root").c_str());
  c3->Clear();

  // Overlay energy definitions: Eng & StartE_tpcAV
  SetHistogramStyle1D(h_eng,"Muon energy [GeV]", "Rate");
  l->SetX1NDC(0.27);
  l->SetX2NDC(0.96);

  float maxy = std::max(h_energy->GetMaximum(),h_eng->GetMaximum());

  h_eng->Draw("hist");
  h_energy->Draw("hist same");
  h_eng->SetLineWidth(3);
  h_eng->SetLineStyle(2);
  h_eng->SetLineColor(pal.at(0));
  h_eng->GetYaxis()->SetTitleOffset(0.95);
  h_eng->GetYaxis()->SetRangeUser(0,maxy*1.1);
  
  l->AddEntry(h_energy, "TPC Energy", "l");
  l->AddEntry(h_eng, "Generated Energy", "l");
  l->Draw();

  c3->SaveAs((location+"/truth_tracks_energy_eng"+tag+".png").c_str());
  c3->SaveAs((location+"/truth_tracks_energy_eng"+tag+".root").c_str());
  c3->Clear();
  l->Clear();


  // Overlay scaled and unscaled energy
  c3->SetLogx();
  TH1D *h_scaled = static_cast<TH1D*>(h_energyLog->Clone("h_scaled"));
  SetHistogramStyle1D(h_energyLog,"Muon energy [GeV]", "Rate");
  SetHistogramStyle1D(h_scaled,"Muon energy [GeV]", "Rate");
  h_scaled->Scale(1,"width");

  maxy = std::max(h_energy->GetMaximum(),h_scaled->GetMaximum());
  
  h_energyLog->Draw("hist");
  h_scaled->Draw("hist same");
  h_energyLog->SetLineWidth(3);
  h_energyLog->SetLineStyle(2);
  h_energyLog->SetLineColor(pal.at(1));
  h_scaled->SetLineWidth(3);
  h_scaled->SetLineStyle(2);
  h_scaled->SetLineColor(pal.at(0));
  h_energyLog->GetYaxis()->SetTitleOffset(0.95);
  h_energyLog->GetYaxis()->SetRangeUser(0,maxy*1.1);
  
  l->AddEntry(h_energyLog, "No scaling", "l");
  l->AddEntry(h_scaled, "Scaled by bin width", "l");
  l->Draw();

  c3->SaveAs((location+"/truth_tracks_energy_scale_noScale"+tag+".png").c_str());
  c3->SaveAs((location+"/truth_tracks_energy_scale_noScale"+tag+".root").c_str());
  c3->Clear();
  l->Clear();

  // Separating the beam plug and halo muons
  SetHistogramStyle1D(h_energyPlug,"Muon energy [GeV]", "Rate");
  SetHistogramStyle1D(h_energyHalo,"Muon energy [GeV]", "Rate");

  maxy = std::max(h_energyPlug->GetMaximum(),h_energyHalo->GetMaximum());

  h_energyPlug->Draw("hist");
  h_energyHalo->Draw("hist same");
  h_energyPlug->SetLineWidth(3);
  h_energyHalo->SetLineWidth(3);
  h_energyPlug->SetLineStyle(2);
  h_energyHalo->SetLineStyle(2);
  h_energyPlug->SetLineColor(pal.at(1));
  h_energyHalo->SetLineColor(pal.at(0));
  h_energyPlug->GetYaxis()->SetTitleOffset(0.95);
  h_energyPlug->GetYaxis()->SetRangeUser(0,maxy*1.1);

  l->AddEntry(h_energyPlug, "Beam plug #mu", "l");
  l->AddEntry(h_energyHalo, "Halo #mu", "l");
  l->Draw();

  c3->SaveAs((location+"/truth_tracks_energy_plug_halo"+tag+".png").c_str());
  c3->SaveAs((location+"/truth_tracks_energy_plug_halo"+tag+".root").c_str());
  c3->Clear();
  l->Clear();

  TCanvas *c4 = new TCanvas("c4","",900,900);
  SetCanvasStyle(c4, 0.12,0.08,0.06,0.12,0,0,0);

  c4->SetLogy();
  SetHistogramStyle1D(h_enter_dist,"Distance from candidate entrance/exit [m]", " Rate");
  h_enter_dist->Draw("hist");
  h_exit_dist->Draw("same");
  h_enter_dist->SetLineWidth(3);
  h_exit_dist->SetLineWidth(3);
  h_enter_dist->SetLineStyle(2);
  h_exit_dist->SetLineStyle(2);
  h_enter_dist->SetLineColor(pal.at(0));
  h_exit_dist->SetLineColor(pal.at(1));
  h_enter_dist->GetYaxis()->SetTitleOffset(0.95);

  l->AddEntry(h_enter_dist, "Entrance", "l");
  l->AddEntry(h_exit_dist, "Exit", "l");
  l->Draw();

  c4->SaveAs((location+"/truth_distance_to_entrance_exit_planes"+tag+".png").c_str());
  c4->SaveAs((location+"/truth_distance_to_entrance_exit_planes"+tag+".root").c_str());
  c4->Clear();
  l->Clear();
  
  TCanvas *c5 = new TCanvas("c5","",1000,800);
  SetCanvasStyle(c5, 0.1,0.12,0.05,0.12,0,0,0);

  for(unsigned int iWire = 0; iWire < 3; ++iWire){
    // Energies
    SetHistogramStyle2D(h_energies.at(iWire),"x [cm]", " Energy deposition [MeV/cm]",false);
    h_energies.at(iWire)->GetZaxis()->SetLabelSize(0.03);
    h_energies.at(iWire)->GetZaxis()->SetLabelFont(132);
    h_energies.at(iWire)->Draw("colz");

    c5->SaveAs((location+"/energy_vs_X_plane"+std::to_string(iWire)+tag+".png").c_str());
    c5->SaveAs((location+"/energy_vs_X_plane"+std::to_string(iWire)+tag+".root").c_str());
    c5->Clear();
    
    SetHistogramStyle2D(h_corr_energies.at(iWire),"x [cm]", " Energy deposition [MeV/cm]",false);
    h_corr_energies.at(iWire)->GetZaxis()->SetLabelSize(0.03);
    h_corr_energies.at(iWire)->GetZaxis()->SetLabelFont(132);
    h_corr_energies.at(iWire)->Draw("colz");

    c5->SaveAs((location+"/corr_energy_vs_X_plane"+std::to_string(iWire)+tag+".png").c_str());
    c5->SaveAs((location+"/corr_energy_vs_X_plane"+std::to_string(iWire)+tag+".root").c_str());
    c5->Clear();
    
    SetHistogramStyle2D(h_corr_nelec_energies.at(iWire),"x [cm]", " Energy deposition per e^{-} [MeV/e^{-}]",false);
    h_corr_nelec_energies.at(iWire)->GetZaxis()->SetLabelSize(0.03);
    h_corr_nelec_energies.at(iWire)->GetZaxis()->SetLabelFont(132);
    h_corr_nelec_energies.at(iWire)->Draw("colz");

    c5->SaveAs((location+"/energy_nelec_vs_X_plane"+std::to_string(iWire)+tag+".png").c_str());
    c5->SaveAs((location+"/energy_nelec_vs_X_plane"+std::to_string(iWire)+tag+".root").c_str());
    c5->Clear();
    
    SetHistogramStyle2D(h_corr_nelec_corr_energies.at(iWire),"x [cm]", " Energy deposition per e^{-} [MeV/e^{-}]",false);
    h_corr_nelec_corr_energies.at(iWire)->GetZaxis()->SetLabelSize(0.03);
    h_corr_nelec_corr_energies.at(iWire)->GetZaxis()->SetLabelFont(132);
    h_corr_nelec_corr_energies.at(iWire)->Draw("colz");

    c5->SaveAs((location+"/energy_nelec_corr_vs_X_plane"+std::to_string(iWire)+tag+".png").c_str());
    c5->SaveAs((location+"/energy_nelec_corr_vs_X_plane"+std::to_string(iWire)+tag+".root").c_str());
    c5->Clear();
    
    SetHistogramStyle2D(h_corr_charge_energies.at(iWire),"x [cm]", " dE/dQ [MeV/ADC]",false);
    h_corr_charge_energies.at(iWire)->GetZaxis()->SetLabelSize(0.03);
    h_corr_charge_energies.at(iWire)->GetZaxis()->SetLabelFont(132);
    h_corr_charge_energies.at(iWire)->Draw("colz");

    c5->SaveAs((location+"/energy_charge_vs_X_plane"+std::to_string(iWire)+tag+".png").c_str());
    c5->SaveAs((location+"/energy_charge_vs_X_plane"+std::to_string(iWire)+tag+".root").c_str());
    c5->Clear();
    
    SetHistogramStyle2D(h_corr_charge_corr_energies.at(iWire),"x [cm]", " dE/dQ [MeV/ADC]",false);
    h_corr_charge_corr_energies.at(iWire)->GetZaxis()->SetLabelSize(0.03);
    h_corr_charge_corr_energies.at(iWire)->GetZaxis()->SetLabelFont(132);
    h_corr_charge_corr_energies.at(iWire)->Draw("colz");

    c5->SaveAs((location+"/energy_charge_corr_vs_X_plane"+std::to_string(iWire)+tag+".png").c_str());
    c5->SaveAs((location+"/energy_charge_corr_vs_X_plane"+std::to_string(iWire)+tag+".root").c_str());
    c5->Clear();
    
    // Charges
    SetHistogramStyle2D(h_charges.at(iWire),"x [cm]", " Charge deposition [ADC/cm]",false);
    h_charges.at(iWire)->GetZaxis()->SetLabelSize(0.03);
    h_charges.at(iWire)->GetZaxis()->SetLabelFont(132);
    h_charges.at(iWire)->Draw("colz");

    c5->SaveAs((location+"/charge_vs_X_plane"+std::to_string(iWire)+tag+".png").c_str());
    c5->SaveAs((location+"/charge_vs_X_plane"+std::to_string(iWire)+tag+".root").c_str());
    c5->Clear();
    
    SetHistogramStyle2D(h_corr_charges.at(iWire),"x [cm]", " Charge deposition [ADC/cm]",false);
    h_corr_charges.at(iWire)->GetZaxis()->SetLabelSize(0.03);
    h_corr_charges.at(iWire)->GetZaxis()->SetLabelFont(132);
    h_corr_charges.at(iWire)->Draw("colz");

    c5->SaveAs((location+"/corr_charge_vs_X_plane"+std::to_string(iWire)+tag+".png").c_str());
    c5->SaveAs((location+"/corr_charge_vs_X_plane"+std::to_string(iWire)+tag+".root").c_str());
    c5->Clear();
    
    // # Electrons
    SetHistogramStyle2D(h_nelecs.at(iWire),"x [cm]", " Number of electrons",false);
    h_nelecs.at(iWire)->GetZaxis()->SetLabelSize(0.03);
    h_nelecs.at(iWire)->GetZaxis()->SetLabelFont(132);
    h_nelecs.at(iWire)->Draw("colz");

    c5->SaveAs((location+"/nelecs_vs_X_plane"+std::to_string(iWire)+tag+".png").c_str());
    c5->SaveAs((location+"/nelecs_vs_X_plane"+std::to_string(iWire)+tag+".root").c_str());
    c5->Clear();
    
    SetHistogramStyle2D(h_corr_nelecs.at(iWire),"x [cm]", " Number of electrons",false);
    h_corr_nelecs.at(iWire)->GetZaxis()->SetLabelSize(0.03);
    h_corr_nelecs.at(iWire)->GetZaxis()->SetLabelFont(132);
    h_corr_nelecs.at(iWire)->Draw("colz");

    c5->SaveAs((location+"/corr_nelecs_vs_X_plane"+std::to_string(iWire)+tag+".png").c_str());
    c5->SaveAs((location+"/corr_nelecs_vs_X_plane"+std::to_string(iWire)+tag+".root").c_str());
    c5->Clear();
    
    SetHistogramStyle2D(h_corr_nelec_corr_charges.at(iWire),"x [cm]", " e^{-} per charge [N_{e^{-}}/ADC]",false);
    h_corr_nelec_corr_charges.at(iWire)->GetZaxis()->SetLabelSize(0.03);
    h_corr_nelec_corr_charges.at(iWire)->GetZaxis()->SetLabelFont(132);
    h_corr_nelec_corr_charges.at(iWire)->Draw("colz");

    c5->SaveAs((location+"/corr_nelecs_corr_charge_vs_X_plane"+std::to_string(iWire)+tag+".png").c_str());
    c5->SaveAs((location+"/corr_nelecs_corr_charge_vs_X_plane"+std::to_string(iWire)+tag+".root").c_str());
    c5->Clear();
    
    SetHistogramStyle2D(h_eDep_E.at(iWire),"Muon energy [GeV]", "Energy deposition per unit length [MeV/cm]",false);
    h_eDep_E.at(iWire)->GetZaxis()->SetLabelSize(0.03);
    h_eDep_E.at(iWire)->GetZaxis()->SetLabelFont(132);
    h_eDep_E.at(iWire)->Draw("colz");
    c5->SaveAs((location+"/truth_tracks_eDep_vs_E"+std::to_string(iWire)+tag+".png").c_str());
    c5->SaveAs((location+"/truth_tracks_eDep_vs_E"+std::to_string(iWire)+tag+".root").c_str());
    c5->Clear();

  } // Wire planes
  
  TCanvas *c6 = new TCanvas("c6","",1100,800);
  SetCanvasStyle(c6, 0.1,0.132,0.05,0.12,0,0,0);

  SetHistogramStyle2D(h_startX_energy, "Start X [cm]", " TPC Energy [GeV]",false);
  h_startX_energy->GetZaxis()->SetLabelSize(0.03);
  h_startX_energy->GetZaxis()->SetLabelFont(132);
  h_startX_energy->Draw("colz");
 
  c6->SetLogy();
  c6->SaveAs((location+"/energy_vs_startX"+tag+".png").c_str());
  c6->SaveAs((location+"/energy_vs_startX"+tag+".root").c_str());
  c6->Clear();

  // Now with scaling
  TH2D *h_startX_energy_scaled = static_cast<TH2D*>(h_startX_energy->Clone("h_startX_energy_scaled"));
  SetHistogramStyle2D(h_startX_energy_scaled, "Start X [cm]", " TPC Energy [GeV]",false);
  h_startX_energy_scaled->Scale(1,"width");
  h_startX_energy_scaled->GetZaxis()->SetLabelSize(0.03);
  h_startX_energy_scaled->GetZaxis()->SetLabelFont(132);
  h_startX_energy_scaled->Draw("colz");
 
  c6->SaveAs((location+"/energy_vs_startX_scaled"+tag+".png").c_str());
  c6->SaveAs((location+"/energy_vs_startX_scaled"+tag+".root").c_str());
  c6->Clear();
  
  // Now try overlaying the 1D x distribution
  h_startX_energy->GetZaxis()->SetLabelSize(0.03);
  h_startX_energy->GetZaxis()->SetLabelFont(132);
  h_startX_energy->Draw("colz");
 
  TPaletteAxis *paletteX = static_cast<TPaletteAxis*>(h_startX_energy->GetListOfFunctions()->FindObject("palette"));

  // the following lines move the palette. Choose the values you need for the position.
  paletteX->SetX1NDC(0.915);
  paletteX->SetX2NDC(0.965);
  gPad->Modified();
  gPad->Update();

  // Try and overlay the rate of events
  // Scale startX to the pad coordinates
  TPad *pX = new TPad("pX", "", 0, 0, 1, 1);
  pX->SetFillStyle(0);
  SetPadStyle(pX, 0.1,0.132,0.05,0.12);
  pX->Draw();
  pX->cd();

  h_startX->SetLineColor(pal.at(2));
  h_startX->SetLineStyle(2);
  h_startX->SetLineWidth(3);
  h_startX->Draw("hist");
  
  //draw an axis on the right side
  float rightmax = 1.05*h_startX->GetMaximum();
  TGaxis *axisX = new TGaxis(h_startX->GetXaxis()->GetXmax(),
                             h_startX->GetYaxis()->GetXmin(),
                             h_startX->GetXaxis()->GetXmax(),
                             rightmax,
                             0, rightmax, 510, "+L");
  axisX->SetLineColor(pal.at(1));
  axisX->SetLabelColor(pal.at(1));
  axisX->SetLabelFont(132);
  axisX->SetMoreLogLabels();
  axisX->Draw();

  h_startX->GetXaxis()->SetRangeUser(h_startX_energy->GetXaxis()->GetXmin(),h_startX_energy->GetXaxis()->GetXmax());
  h_startX->GetXaxis()->SetTickSize(0);
  h_startX->GetXaxis()->SetTitleSize(0);
  h_startX->GetXaxis()->SetLabelSize(0);
  h_startX->GetYaxis()->SetTickSize(0);
  h_startX->GetYaxis()->SetTitleSize(0);
  h_startX->GetYaxis()->SetLabelSize(0);
  pX->Draw();

  TFrame *fX = pX->GetFrame();
  fX->SetFillStyle(0);
  fX->SetLineStyle(0);
  fX->SetBorderMode(0);
  fX->SetBorderSize(0);

  pX->Modified();
  pX->Update();
  pX->Draw();
  c6->SaveAs((location+"/energy_vs_startX_overlay"+tag+".png").c_str());
  c6->SaveAs((location+"/energy_vs_startX_overlay"+tag+".root").c_str());
  c6->Clear();
 
  TCanvas *c7 = new TCanvas("c7","",1100,800);
  SetCanvasStyle(c7, 0.1,0.132,0.05,0.12,0,0,0);

  SetHistogramStyle2D(h_startY_energy, "Start Y [cm]", " TPC Energy [GeV]",false);
  h_startY_energy->GetZaxis()->SetLabelSize(0.03);
  h_startY_energy->GetZaxis()->SetLabelFont(132);
  h_startY_energy->Draw("colz");

  c7->SetLogy();
  c7->SaveAs((location+"/energy_vs_startY"+tag+".png").c_str());
  c7->SaveAs((location+"/energy_vs_startY"+tag+".root").c_str());
  c7->Clear();
    
  // Now with scaling
  TH2D *h_startY_energy_scaled = static_cast<TH2D*>(h_startY_energy->Clone("h_startY_energy_scaled"));
  SetHistogramStyle2D(h_startY_energy_scaled, "Start Y [cm]", " TPC Energy [GeV]",false);
  h_startY_energy_scaled->Scale(1,"width");
  h_startY_energy_scaled->GetZaxis()->SetLabelSize(0.03);
  h_startY_energy_scaled->GetZaxis()->SetLabelFont(132);
  h_startY_energy_scaled->Draw("colz");
 
  c7->SetLogy();
  c7->SaveAs((location+"/energy_vs_startY_scaled"+tag+".png").c_str());
  c7->SaveAs((location+"/energy_vs_startY_scaled"+tag+".root").c_str());
  c7->Clear();
  
  // Try and overlay the rate of events
  // Now try overlaying the 1D x distribution
  h_startY_energy->GetZaxis()->SetLabelSize(0.03);
  h_startY_energy->GetZaxis()->SetLabelFont(132);
  h_startY_energy->Draw("colz");
 
  TPaletteAxis *paletteY = static_cast<TPaletteAxis*>(h_startY_energy->GetListOfFunctions()->FindObject("palette"));

  // the following lines move the palette. Choose the values you need for the position.
  paletteY->SetX1NDC(0.915);
  paletteY->SetX2NDC(0.965);
  gPad->Modified();
  gPad->Update();

  // Scale startY to the pad coordinates
  TPad *pY = new TPad("pY", "", 0, 0, 1, 1);
  pY->SetFillStyle(0);
  SetPadStyle(pY, 0.1,0.132,0.05,0.12);
  pY->Draw();
  pY->cd();

  h_startY->SetLineColor(pal.at(2));
  h_startY->SetLineStyle(2);
  h_startY->SetLineWidth(3);
  h_startY->Draw("hist");
  
  //draw an axis on the right side
  rightmax = 1.05*h_startY->GetMaximum();
  TGaxis *axisY = new TGaxis(h_startY->GetXaxis()->GetXmax(),
                             h_startY->GetYaxis()->GetXmin(),
                             h_startY->GetXaxis()->GetXmax(),
                             rightmax,
                             0, rightmax, 510, "+L");
  axisY->SetLineColor(pal.at(1));
  axisY->SetLabelColor(pal.at(1));
  axisY->SetLabelFont(132);
  axisY->SetMoreLogLabels();
  axisY->Draw();

  h_startY->GetXaxis()->SetRangeUser(h_startY_energy->GetXaxis()->GetXmin(),h_startY_energy->GetXaxis()->GetXmax());
  h_startY->GetXaxis()->SetTickSize(0);
  h_startY->GetXaxis()->SetTitleSize(0);
  h_startY->GetXaxis()->SetLabelSize(0);
  h_startY->GetYaxis()->SetTickSize(0);
  h_startY->GetYaxis()->SetTitleSize(0);
  h_startY->GetYaxis()->SetLabelSize(0);
  pY->Draw();

  TFrame *fY = pY->GetFrame();
  fY->SetFillStyle(0);
  fY->SetLineStyle(0);
  fY->SetBorderMode(0);
  fY->SetBorderSize(0);

  pY->Modified();
  pY->Update();
  pY->Draw();
  c7->SaveAs((location+"/energy_vs_startY_overlay"+tag+".png").c_str());
  c7->SaveAs((location+"/energy_vs_startY_overlay"+tag+".root").c_str());
  c7->Clear();
 
  // Now try making an event display, if we like
  if(drawMCParticles){
    gStyle->SetTextFont(132);
    TCanvas *cEvd = new TCanvas("cEvd","",1200,900);
    SetCanvasStyle(cEvd, 0.06,0.06,0.06,0.06,0,0,0);
    
    // Creating a view
    TView3D *view = (TView3D*) TView::CreateView(1);
    view->SetRange(-450.,-100.,-50.,450.,720.,620.);
    view->SetOutlineToCube();
    view->RotateView(110,135);

    // Loop over the number of planes we have
    unsigned int p = 0;
    TPolyMarker3D *pl3detmids = new TPolyMarker3D(6);
    for(const Plane &pl : allPlanes){
      std::vector<TVector3> boundaryCoordinates;
      GetPlaneBoundaryCoordinates(pl, boundaryCoordinates);
      // Draw the detector boundaries
      TPolyLine3D *pl3det = new TPolyLine3D(5);
      for(unsigned int bc = 0; bc < boundaryCoordinates.size(); ++bc){
        pl3det->SetPoint(bc, boundaryCoordinates.at(bc).X(), boundaryCoordinates.at(bc).Z(), boundaryCoordinates.at(bc).Y());
      }
      pl3det->SetPoint(4, boundaryCoordinates.at(0).X(), boundaryCoordinates.at(0).Z(), boundaryCoordinates.at(0).Y());
      pl3det->SetLineColor(kGray+3);
      pl3det->SetLineWidth(1);
      pl3det->SetLineStyle(2);
      pl3det->Draw();
      pl3detmids->SetPoint(p, pl.GetV().X(), pl.GetV().Z(), pl.GetV().Y());
      pl3detmids->SetMarkerStyle(43);
      pl3detmids->SetMarkerColor(kGray+3);
      p++;
    } // planes

    // Now draw the mcparticles
    unsigned int palIt = 0;
    for(unsigned int i = 0; i < startPoints.size(); ++i){
      if(palIt >= pal.size())
        palIt = 0;
    
      TPolyLine3D *pl3part = new TPolyLine3D(2);
      pl3part->SetPoint(0, startPoints.at(i).X(), startPoints.at(i).Z(), startPoints.at(i).Y());
      pl3part->SetPoint(1, endPoints.at(i).X(), endPoints.at(i).Z(), endPoints.at(i).Y());
      pl3part->SetLineColor(pal.at(palIt));
      pl3part->SetLineWidth(1);
      pl3part->SetLineStyle(1);
      pl3part->Draw();
      
      TPolyMarker3D *pl3partdir = new TPolyMarker3D(1);
      pl3partdir->SetPoint(0, endPoints.at(i).X(), endPoints.at(i).Z(), endPoints.at(i).Y());
      pl3partdir->SetMarkerStyle(33);
      pl3partdir->SetMarkerColor(pal.at(palIt));
      pl3partdir->Draw();
      
      ++palIt;
    }
    pl3detmids->Draw();
    
    TPolyLine3D *pl3beam = new TPolyLine3D(2);
    pl3beam->SetPoint(0, beamStart.X(), beamStart.Z(), beamStart.Y());
    pl3beam->SetPoint(1, beamEnd.X(), beamEnd.Z(), beamEnd.Y());
    pl3beam->SetLineColor(kGray+3);
    pl3beam->SetLineWidth(4);
    pl3beam->SetLineStyle(1);
    pl3beam->Draw();
      
    TPolyMarker3D *pl3beamends = new TPolyMarker3D(2);
    pl3beamends->SetPoint(0, beamStart.X(), beamStart.Z(), beamStart.Y());
    pl3beamends->SetPoint(1, beamEnd.X(), beamEnd.Z(), beamEnd.Y());
    pl3beamends->SetMarkerStyle(49);
    pl3beamends->SetMarkerColor(kGray+3);
    pl3beamends->Draw();
    
    view->ToggleRulers();

    std::cout << " There are " << startPoints.size() << " muons " << std::endl;

    cEvd->SaveAs((location+"/mcparticle_display"+tag+".png").c_str());
    cEvd->SaveAs((location+"/mcparticle_display"+tag+".root").c_str());
  
  } // Draw MCParticles
  std::cout << " There are " << nThru << " through-going muons " << std::endl;

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
