/************************************************************************
 * 
 * A macro to understand the behaviour of delta-ray activity surrounding 
 * cosmic-ray muons and to produce plots for the absolute dE/dx 
 * calibration procedure
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
   "origin",
   "taulife",
   "inTPCActive",
   "TrackId",
   "pdg",
   "Mother",
   "ntracks_pandoraTrack",
   "trkId_pandoraTrack",
   "trkidtruth_pandoraTrack",
   "trkpdgtruth_pandoraTrack",
   "trkg4id_pandoraTrack",
   "trkorigin_pandoraTrack",
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
   "StartE_tpcAV",
   "EndPointx",
   "EndPointy",
   "EndPointz",
   "EndE_tpcAV",
   "StartPointx_tpcAV",
   "StartPointy_tpcAV",
   "StartPointz_tpcAV",
   "EndPointx_tpcAV",
   "EndPointy_tpcAV",
   "EndPointz_tpcAV",
   "NumberDaughters",
   "P",
   "Eng",
   "EndE",
 };

// A translation list from plane labels to longer labels for plotting
std::map<std::string, std::string> planeLabels = {
  {"h0", "APA 1"},
  {"h1", "CPA 1"},
  {"h2", "APA 2"},
  {"t",  "Top"},
  {"bo", "Bot."},
  {"f",  "Fro."},
  {"ba", "Back"},
};

typedef std::vector<Plane> PlaneList;
     
int activityStudies(const char *config){

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
  int cry = 0; // Whether or not to save the cosmic (cry) information (trkorigin = 2) or the beam information (trkorigin = 4)
  int n = -1;
  int thru = 0;
  std::vector<double> etau; // measured electron lifetime, one per TPC if desired
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
  p->getValue("ETau",      etau);
  p->getValue("NFiles",    n);
  p->getValue("Cry",       cry);
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

  // Global variables
  // Normal to the APA plane
  TVector3 apaNorm={0,0,0};
  for(const Plane &pl: allPlanes){
    if(pl.GetLabel() != "h0") continue;
    apaNorm = pl.GetUnitN();
    break;
  }

  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Total number of planes in the active volume of the DUNE SP module: " << allPlanes.size() << std::endl;
  std::cout << " Consisting of " << extPlanes.size() << " external planes and " << intPlanes.size() << " internal planes" << std::endl; 
  std::cout << "-----------------------------------------------------------" << std::endl;
 
  // Sort out the file tag
  if(tag != "")
    tag = "_"+tag;

  // Get the lifetime measured in each TPC
  if(etau.size() == 0) etau.push_back(35);
  bool LifetimePerTPC = etau.size() > 1 ? 1 : 0;

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

  //TFile *fSake = new TFile((location+"/dummy"+tag+".root").c_str(), "RECREATE");

  // Then setup the histograms, counters and any other variables to add to
  // Setup histograms
  // Truth-level track quantities
  TH1D *h_hit_pitch            = new TH1D("h_hit_pitch","",80,0,2);   // Pitch of the muon hits
  TH1D *h_length               = new TH1D("h_length","",80,0,12);   // Length of the muons
  TH1D *h_mom                  = new TH1D("h_mom","",80,0,40);       // Momentum of the muons [GeV]
  TH1D *h_energy               = new TH1D("h_energy","",80,0.8,40);       // Energy of the muons [GeV]
  TH1D *h_energy_long          = new TH1D("h_energy_long","",80,0.8,40);       // Energy of the muons [GeV]
  TH1D *h_energy_nolog         = new TH1D("h_energy_nolog","",80,0,40);       // Energy of the muons [GeV]
  TH1D *h_energy_long_nolog    = new TH1D("h_energy_long_nolog","",80,0,40);       // Energy of the muons [GeV]
  TH1D *h_hitE_long_BP         = new TH1D("h_hitE_long_BP","",80,0.8,40);       // Total true energy depositions of the muons [GeV]
  TH1D *h_nDaughters           = new TH1D("h_nDaughters","",50,0,200); // Number of muon daughters
  TH1D *h_reco_eng             = new TH1D("h_reco_eng","",80,0.8,40);       // Energy of the muons [GeV]
  TH1D *h_reco_eng_long        = new TH1D("h_reco_eng_long","",80,0.8,40);       // Energy of the muons [GeV]
  TH1D *h_reco_eng_highy       = new TH1D("h_reco_eng_highy","",80,0.8,40);       // Energy of the muons [GeV]
  TH1D *h_reco_eng_long_highy  = new TH1D("h_reco_eng_long_highy","",80,0.8,40);       // Energy of the muons [GeV]
  TH1D *h_reco_len             = new TH1D("h_reco_len","",80,0,12);   // Length of the muons [m]
  TH1D *h_EoL_BP               = new TH1D("h_EoL_BP","",80,0,50); // True energy over trajectory length of the muon
  TH1D *h_dEdx_BP              = new TH1D("h_dEdx_BP","",80,0,50);  // dE/dx of the muons [GeV]
  TH1D *h_hit_widthX_BP        = new TH1D("h_hit_widthX_BP","",80,6,7.5);  // hit width [cm]
  TH1D *h_hit_widthTicks_BP    = new TH1D("h_hit_widthTicks_BP","",80,4,10);  // hit width [cm]
  TH1D *h_dEdx_hitCut_BP       = new TH1D("h_dEdx_hitCut_BP","",80,0,7);  // dE/dx of the muons [GeV]
  TH1D *h_true_mus             = new TH1D("h_true_mus","",101,-0.5,100.5); // True muon multiplicity
  TH1D *h_true_primary_mus     = new TH1D("h_true_primary_mus","",12,-0.5,11.5); // True muon multiplicity
  TH1D *h_reco_mus             = new TH1D("h_reco_mus","",101,-0.5,100.5); // Reco muon multiplicity
  TH1D *h_reco_long_mus        = new TH1D("h_reco_long_mus","",12,-0.5,11.5); // Reco muon multiplicity
  TH1D *h_reco_long_highy_mus  = new TH1D("h_reco_long_highy_mus","",12,-0.5,11.5); // Reco muon multiplicity
  TH1D *h_nHitsPerL_BP         = new TH1D("h_nHitsPerL_BP","",80,0,2); // Number of hits per unit length
  TH1D *h_costheta_BP_0        = new TH1D("h_costheta_BP_0","",80,-1,1); // Angle to the best plane if the best plane is 0
  TH1D *h_costheta_BP_1        = new TH1D("h_costheta_BP_1","",80,-1,1); // Angle to the best plane if the best plane is 1
  TH1D *h_costheta_BP_2        = new TH1D("h_costheta_BP_2","",80,-1,1); // Angle to the best plane if the best plane is 2
  TH1D *h_costheta_BP_0_hitCut = new TH1D("h_costheta_BP_0_hitCut","",80,-1,1); // Angle to the best plane if the best plane is 0 with the hits/length cut
  TH1D *h_costheta_BP_1_hitCut = new TH1D("h_costheta_BP_1_hitCut","",80,-1,1); // Angle to the best plane if the best plane is 1 with the hits/length cut
  TH1D *h_costheta_BP_2_hitCut = new TH1D("h_costheta_BP_2_hitCut","",80,-1,1); // Angle to the best plane if the best plane is 2 with the hits/length cut
  TH1D *h_xPos_BP              = new TH1D("h_xPos_BP","",80,-400,400); // Distribution of deposition x positions

  TH2D *h_E_nDaught                 = new TH2D("h_E_nDaught","",200,0,800,51,-0.5,50.5); // Number of muon daughters per unit energy
  TH2D *h_dEdx_nDaught_0            = new TH2D("h_dEdx_nDaught_0","",51,-0.5,50.5,50,0.2,4); // Number of muon daughters vs energy depositions per unit length
  TH2D *h_dEdx_nDaught_1            = new TH2D("h_dEdx_nDaught_1","",51,-0.5,50.5,50,0.2,4); // Number of muon daughters vs energy depositions per unit length
  TH2D *h_dEdx_nDaught_2            = new TH2D("h_dEdx_nDaught_2","",51,-0.5,50.5,50,0.2,5); // Number of muon daughters vs energy depositions per unit length
  TH2D *h_dEdx_E_0                  = new TH2D("h_dEdx_E_0","",50,0.8,40,50,0.2,10); // Energy deposition vs energy
  TH2D *h_dEdx_E_1                  = new TH2D("h_dEdx_E_1","",50,0.8,40,50,0.2,10); // Energy deposition vs energy
  TH2D *h_dEdx_E_2                  = new TH2D("h_dEdx_E_2","",50,0.8,40,50,0.2,10); // Energy deposition vs energy
  TH2D *h_dEdx_E_BP                 = new TH2D("h_dEdx_E_BP","",50,0.8,40,50,0.2,10); // Energy deposition vs energy
  TH2D *h_dEdx_hitCut_E_0           = new TH2D("h_dEdx_hitCut_E_0","",50,0.8,40,50,0.2,10); // Energy deposition vs energy
  TH2D *h_dEdx_hitCut_E_1           = new TH2D("h_dEdx_hitCut_E_1","",50,0.8,40,50,0.2,10); // Energy deposition vs energy
  TH2D *h_dEdx_hitCut_E_2           = new TH2D("h_dEdx_hitCut_E_2","",50,0.8,40,50,0.2,10); // Energy deposition vs energy
  TH2D *h_dEdx_hitCut_E_BP          = new TH2D("h_dEdx_hitCut_E_BP","",50,0.8,40,50,0.2,10); // Energy deposition vs energy
  TH2D *h_dQdx_hitWidth_BP          = new TH2D("h_dQdx_hitWidth_BP","",50,6,7.5,50,80,600); // dQ/dx vs hit width
  TH2D *h_reco_dQdx_simCorr_E       = new TH2D("h_reco_dQdx_simCorr_E","",50,0.8,40,50,80,600); // dQ/dx vs energy with simulated lifetime correction
  TH2D *h_reco_dQdx_E               = new TH2D("h_reco_dQdx_E","",50,0.8,40,50,80,600); // dQ/dx vs energy
  TH2D *h_reco_dQdx_pos             = new TH2D("h_reco_dQdx_pos","",50,-400,400,50,80,600); // dQ/dx vs x position
  TH2D *h_reco_dQdx_dP              = new TH2D("h_reco_dQdx_dP","",50,0.45,0.6,50,80,600); // dQ/dx vs pitch
  TH2D *h_reco_dQdx_RR              = new TH2D("h_reco_dQdx_RR","",50,0,8,50,80,600); // dQ/dx vs residual range
  TH2D *h_reco_dQdx_width           = new TH2D("h_reco_dQdx_width","",50,6,7.5,50,80,600); // dQ/dx vs hit width
  TH2D *h_reco_dQdx_cosDrift        = new TH2D("h_reco_dQdx_cosDrift","",50,-1,1,50,80,600); // dQ/dx vs angle to drift direction (x)
  TH2D *h_reco_dEdx_RR_BP           = new TH2D("h_reco_dEdx_RR_BP","",50,0,8,50,0,7); // dE/dx vs energy, best plane
  TH2D *h_reco_dEdx_dP_BP           = new TH2D("h_reco_dEdx_dP_BP","",50,0.45,0.6,50,0,7); // dE/dx vs hit pitch, best plane
  TH2D *h_reco_dEdx_dQdx_BP         = new TH2D("h_reco_dEdx_dQdx_BP","",50,0,7,50,80,600); // dE/dx vs dQ/dx, best plane
  TH2D *h_reco_dEdx_noCorr_E_BP     = new TH2D("h_reco_dEdx_noCorr_E_BP","",50,0.8,40,50,0,7); // dE/dx vs energy, best plane, no lifetime correction
  TH2D *h_reco_dEdx_E_BP            = new TH2D("h_reco_dEdx_E_BP","",50,0.8,40,50,0,7); // dE/dx vs energy, best plane
  TH2D *h_reco_dEdx_pos_BP          = new TH2D("h_reco_dEdx_pos_BP","",50,-400,400,50,0,7); // dE/dx vs position, best plane
  TH2D *h_reco_dEdx_E_hitCut_0      = new TH2D("h_reco_dEdx_E_hitCut_0","",50,0.8,40,50,0,7); // dE/dx vs energy, best plane, hit/length cut
  TH2D *h_reco_dEdx_E_hitCut_1      = new TH2D("h_reco_dEdx_E_hitCut_1","",50,0.8,40,50,0,7); // dE/dx vs energy, best plane, hit/length cut
  TH2D *h_reco_dEdx_E_hitCut_2      = new TH2D("h_reco_dEdx_E_hitCut_2","",50,0.8,40,50,0,7); // dE/dx vs energy, best plane, hit/length cut
  TH2D *h_reco_dEdx_E_hitCut_BP     = new TH2D("h_reco_dEdx_E_hitCut_BP","",50,0.8,40,50,0,7); // dE/dx vs energy, best plane, hit/length cut
  TH2D *h_reco_dEdx_E_0             = new TH2D("h_reco_dEdx_E_0","",50,0.8,40,50,0.2,6); // Energy deposition vs true energy
  TH2D *h_reco_dEdx_E_1             = new TH2D("h_reco_dEdx_E_1","",50,0.8,40,50,0.2,6); // Energy deposition vs true energy
  TH2D *h_reco_dEdx_E_2             = new TH2D("h_reco_dEdx_E_2","",50,0.8,40,50,0.2,6); // Energy deposition vs energy
  TH2D *h_reco_dQdx_thetaY_0        = new TH2D("h_reco_dQdx_thetaY_0","",50,0.5,1,50,80,600); // Energy deposition vs thetaY
  TH2D *h_reco_dQdx_thetaY_1        = new TH2D("h_reco_dQdx_thetaY_1","",50,0.5,1,50,80,600); // Energy deposition vs thetaY
  TH2D *h_reco_dQdx_thetaY_2        = new TH2D("h_reco_dQdx_thetaY_2","",50,0.5,1,50,80,600); // Energy deposition vs thetaY
  TH2D *h_reco_dEdx_thetaY_0        = new TH2D("h_reco_dEdx_thetaY_0","",50,0.5,1,50,0,5); // Energy deposition vs thetaY
  TH2D *h_reco_dEdx_thetaY_1        = new TH2D("h_reco_dEdx_thetaY_1","",50,0.5,1,50,0,5); // Energy deposition vs thetaY
  TH2D *h_reco_dEdx_thetaY_2        = new TH2D("h_reco_dEdx_thetaY_2","",50,0.5,1,50,0,5); // Energy deposition vs thetaY
  TH2D *h_reco_dEdx_thetaY_hitCut_0 = new TH2D("h_reco_dEdx_thetaY_hitCut_0","",50,0.5,1,50,0,5); // Energy deposition vs thetaY
  TH2D *h_reco_dEdx_thetaY_hitCut_1 = new TH2D("h_reco_dEdx_thetaY_hitCut_1","",50,0.5,1,50,0,5); // Energy deposition vs thetaY
  TH2D *h_reco_dEdx_thetaY_hitCut_2 = new TH2D("h_reco_dEdx_thetaY_hitCut_2","",50,0.5,1,50,0,5); // Energy deposition vs thetaY
  TH2D *h_reco_Y_E                  = new TH2D("h_reco_Y_E","",50,5e-1,8,50,-600,600); // Energy vs reconstructed Y position 
  TH2D *h_reco_Y_E_zoom             = new TH2D("h_reco_Y_E_zoom","",50,5e-1,8,50,596,600); // Energy vs reconstructed Y position 
  TH2D *h_reco_len_E                = new TH2D("h_reco_len_E","",50,5e-1,8,50,2,12); // Energy vs reconstructed len position 
  TH2D *h_nHitsPerL_costheta_BP     = new TH2D("h_nHitsPerL_costheta_BP","",50,-1,1,50,0,2); // Number of hits per unit length vs best plane
  TH2D *h_nHitsPerL_costheta_CP     = new TH2D("h_nHitsPerL_costheta_CP","",50,-0.5,1,50,0.5,2.2); // Number of hits per unit length vs collection plane
  TH2D *h_nHitsPerL_cos_to_plane    = new TH2D("h_nHitsPerL_cos_to_plane","",50,-1,1,50,0,2); // Number of hits per unit length
  TH2D *h_nHitsPerL_recoLength_BP   = new TH2D("h_nHitsPerL_recoLength_BP","",50,0,12,50,0,1); // Number of hits per unit length

  SetLogX(h_reco_dEdx_E_0);
  SetLogX(h_reco_dEdx_E_1);
  SetLogX(h_reco_dEdx_E_2);
  SetLogX(h_reco_eng);
  SetLogX(h_reco_eng_long);
  SetLogX(h_reco_eng_highy);
  SetLogX(h_reco_eng_long_highy);
  SetLogX(h_energy);
  SetLogX(h_energy_long);
  SetLogX(h_dEdx_E_0);
  SetLogX(h_dEdx_E_1);
  SetLogX(h_dEdx_E_2);
  SetLogX(h_dEdx_E_BP);
  SetLogX(h_dEdx_E_BP);
  SetLogX(h_dEdx_hitCut_E_0);
  SetLogX(h_dEdx_hitCut_E_1);
  SetLogX(h_dEdx_hitCut_E_2);
  SetLogX(h_reco_dQdx_E);
  SetLogX(h_reco_dQdx_simCorr_E);
  SetLogX(h_reco_dEdx_noCorr_E_BP);
  SetLogX(h_reco_dEdx_E_BP);
  SetLogX(h_reco_dEdx_E_hitCut_0);
  SetLogX(h_reco_dEdx_E_hitCut_1);
  SetLogX(h_reco_dEdx_E_hitCut_2);
  SetLogX(h_reco_dEdx_E_hitCut_BP);
  SetLogX(h_reco_Y_E);
  SetLogX(h_reco_Y_E_zoom);
  SetLogX(h_reco_len_E);

  std::vector<TH2D*> h_dEdx_nDaught{h_dEdx_nDaught_0,h_dEdx_nDaught_1,h_dEdx_nDaught_2};
  std::vector<TH2D*> h_dEdx_E{h_dEdx_E_0,h_dEdx_E_1,h_dEdx_E_2};
  std::vector<TH2D*> h_dEdx_hitCut_E{h_dEdx_hitCut_E_0,h_dEdx_hitCut_E_1,h_dEdx_hitCut_E_2};
  std::vector<TH2D*> h_reco_dEdx_E_hitCut{h_reco_dEdx_E_hitCut_0,h_reco_dEdx_E_hitCut_1,h_reco_dEdx_E_hitCut_2};
  std::vector<TH2D*> h_reco_dEdx_E{h_reco_dEdx_E_0,h_reco_dEdx_E_1,h_reco_dEdx_E_2};
  std::vector<TH2D*> h_reco_dQdx_thetaY{h_reco_dQdx_thetaY_0,h_reco_dQdx_thetaY_1,h_reco_dQdx_thetaY_2};
  std::vector<TH2D*> h_reco_dEdx_thetaY{h_reco_dEdx_thetaY_0,h_reco_dEdx_thetaY_1,h_reco_dEdx_thetaY_2};
  std::vector<TH2D*> h_reco_dEdx_thetaY_hitCut{h_reco_dEdx_thetaY_hitCut_0,h_reco_dEdx_thetaY_hitCut_1,h_reco_dEdx_thetaY_hitCut_2};
  
  // Setup counters
  double total_energy_true = 0.;
  double total_energy_reco = 0.;
  int n_mu_true = 0;
  int n_mu_reco = 0;
  int n_mu_thru = 0;
  int nHugeE = 0;
  int nBP_0 = 0;
  int nBP_1 = 0;
  int nBP_2 = 0;

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
    unsigned int nTrks  = evt->ntracks_pandoraTrack;
    unsigned int nHits  = evt->no_hits_stored;
    unsigned int nGeant = evt->geant_list_size;

    int n_true_mus            = 0;
    int n_true_primary_mus    = 0;
    int n_reco_mus            = 0;
    int n_reco_long_mus       = 0;
    int n_reco_long_highy_mus = 0;
    
    //
    // Truth-level studies
    //
    // Loop over geant tracks to plot things
    for(unsigned int iG4 = 0; iG4 < nGeant; ++iG4){

      // Check that the origin of the particle matches our choice
      // origin = 2 (Cry, cosmic) origin = 4 (NuWro, beam)
      if((evt->origin[iG4] != 2 && cry) || (evt->origin[iG4] != 4 && !cry)) continue;

      int pdg        = evt->pdg[iG4];
      int id         = evt->TrackId[iG4];

      // Check the particle enters the TPC volume
      if(!evt->inTPCActive[iG4] || evt->StartE_tpcAV[iG4] < 0) continue;

      // Make sure we are looking at a muon
      if(std::abs(pdg) != 13) continue;
      n_true_mus++;
      
      h_energy->Fill(evt->StartE_tpcAV[iG4]);

      // Make sure we are looking at the daughter of a primary particle
      if(evt->Mother[iG4] != 0 && cry) continue;
      n_true_primary_mus++;

      // For the entrance tests
      // If the start or end locations are outside the detector, set them to be at the edge of the detector
      // Only do this for 1 coordinate
      TVector3 vtx(evt->StartPointx[iG4],evt->StartPointy[iG4],evt->StartPointz[iG4]);
      TVector3 end(evt->EndPointx[iG4],evt->EndPointy[iG4],evt->EndPointz[iG4]);
      TVector3 vtxAV(evt->StartPointx_tpcAV[iG4],evt->StartPointy_tpcAV[iG4],evt->StartPointz_tpcAV[iG4]);
      TVector3 endAV(evt->EndPointx_tpcAV[iG4],evt->EndPointy_tpcAV[iG4],evt->EndPointz_tpcAV[iG4]);
      
      float lengthAV = (endAV-vtxAV).Mag()/100.;

      // Check for number of through-going and stopping muons in truth
      bool throughGoing = IsTrueThroughGoing(vtx,end,vtxAV,endAV);
      bool isStopping   = IsTrueStopping(vtx,end,vtxAV,endAV);

      // Fill the truth-track quantities 
      h_length->Fill(lengthAV);
      h_mom->Fill(evt->P[iG4]);
      h_energy_nolog->Fill(evt->StartE_tpcAV[iG4]);

      // For the deposition studies, make sure we are looking at a long track (3m)
      if(lengthAV < 2) continue;
      total_energy_true += evt->Eng[iG4];
      n_mu_true++;
      h_energy_long->Fill(evt->StartE_tpcAV[iG4]);
      h_energy_long_nolog->Fill(evt->StartE_tpcAV[iG4]);

      h_nDaughters->Fill(evt->NumberDaughters[iG4]);
      h_E_nDaught->Fill(evt->StartE_tpcAV[iG4],evt->NumberDaughters[iG4]);
      
      // Get the best plane
      int bestPlane = 0;
      std::vector<int> hitsOnPlane(3,0);
      for(int iPlane = 0; iPlane < 3; ++iPlane){
        for(unsigned int iHit = 0; iHit < nHits; ++iHit){
      
          // Get the ID of the hit
          int hitId  = evt->hit_trkid[iHit];

          // If we are not looking at the current G4 track, continue
          bool currentG4 = false;
          for(int iTrk = 0; iTrk < evt->ntracks_pandoraTrack; ++iTrk){
            int recoId = evt->trkId_pandoraTrack[iTrk];

            // Check if the current hit is in the reco track
            if(recoId != hitId) continue;

            // If it is, check if the reco track is the G4 track
            if(evt->trkidtruth_pandoraTrack[iTrk][iPlane] == id){
              currentG4 = true;
              break;
            }
          }
          if(!currentG4) continue;
          
          if(evt->hit_plane[iHit] == iPlane){
            hitsOnPlane.at(iPlane)++;
          }
        } // Hits
      } // Planes
      if(thru){
        if(!throughGoing) continue;
        else n_mu_thru++;
      }
     
      // Get the best plane
      bestPlane = std::max_element(hitsOnPlane.begin(), hitsOnPlane.end()) - hitsOnPlane.begin();
     
      // Number of hits per unit length
      double hitsPerL = hitsOnPlane.at(bestPlane)/static_cast<double>(lengthAV*100.);
      // Get the number of hits per unit length on the best plane
      h_nHitsPerL_BP->Fill(hitsPerL);

      // Now get the angle of the track to the wires in the best plane
      double costheta      = GetCosTheta(bestPlane,vtxAV,endAV);
      double cosCollection = GetCosTheta(2,vtxAV,endAV);
      h_nHitsPerL_costheta_BP->Fill(costheta,hitsPerL);
      h_nHitsPerL_costheta_CP->Fill(std::abs(cosCollection),hitsPerL);

      // Count up and find out what we're dealing with
      if(bestPlane == 0){
        nBP_0++;
        h_costheta_BP_0->Fill(costheta);
        if(hitsPerL < 0.2){
          h_costheta_BP_0_hitCut->Fill(costheta);
        }
      }
      else if(bestPlane == 1){
        nBP_1++;
        h_costheta_BP_1->Fill(costheta);
        if(hitsPerL < 0.2){
          h_costheta_BP_1_hitCut->Fill(costheta);
        }
      }
      else if(bestPlane == 2){
        nBP_2++;
        h_costheta_BP_2->Fill(costheta);
        if(hitsPerL < 0.2){
          h_costheta_BP_2_hitCut->Fill(costheta);
        }
      }
      
      // Now get the angle of the track to the wire plane
      double costoplane = GetAngleToAPAs(apaNorm,vtxAV,endAV);
      h_nHitsPerL_cos_to_plane->Fill(costoplane,hitsPerL);

      // Now check the reconstructed track length
      for(int iTrk = 0; iTrk < evt->ntracks_pandoraTrack; ++iTrk){
        if(evt->trkidtruth_pandoraTrack[iTrk][bestPlane] != id) continue;
        h_nHitsPerL_recoLength_BP->Fill(evt->trklen_pandoraTrack[iTrk],hitsPerL);
      }

      // Now loop over hits 
      // First loop over wire planes
      for(int iPlane = 0; iPlane < 3; ++iPlane){
        // Skip the first and last hits for 'reconstructability' purposes
        for(unsigned int iHit = 0; iHit < nHits-1; ++iHit){
          // Skip the current hit if it wasn't deposited on this plane
          if(evt->hit_plane[iHit] != iPlane) continue;

          // If we are not looking at the current G4 track, continue
          int hitId      = evt->hit_trkid[iHit];
          bool currentG4 = false;
          for(int iTrk = 0; iTrk < evt->ntracks_pandoraTrack; ++iTrk){
            int recoId = evt->trkId_pandoraTrack[iTrk];

            // Check if the current hit is in the reco track
            if(recoId != hitId) continue;

            // If it is, check if the reco track is the G4 track
            if(evt->trkidtruth_pandoraTrack[iTrk][iPlane] == id){
              currentG4 = true;
              break;
            }
          }
          if(!currentG4) continue;

          // Then get the parameters of interest for this hit
          float hitX      = evt->hit_trueX[iHit];
          float hitE      = evt->hit_energy[iHit];
          float hitQ      = evt->hit_charge[iHit];
          
          // Check if x is lower than the APA bound, charge seems to accumulate there
          // Put a 5 cm fiducial border around the APAs and CPAs too
          if(hitX < evtProc.PD_APA_X_POSITIONS[0]+5 || hitX > evtProc.PD_APA_X_POSITIONS[1]-5) continue;
          if(hitX > evtProc.PD_CPA_X_POSITIONS[0]-5 && hitX < evtProc.PD_CPA_X_POSITIONS[0]+5) continue;
          
          int tpc         = evtProc.WhichTPC(hitX) + 1;
          float dx        = ( -1 + 2*(tpc%2) )*(hitX - evtProc.PD_APA_X_POSITIONS[tpc/2]);
          float hit_width = evt->hit_endT[iHit] - evt->hit_startT[iHit]; // In ticks
          float widthT    = hit_width*0.5e-6; // To s
          float widthX    = widthT/static_cast<float>(evtProc.kXtoT); // To cm
          float pitch     = 0.53; // 5.3 mm peak pitch from activity studies

          h_dEdx_E.at(iPlane)->Fill(evt->StartE_tpcAV[iG4],hitE/pitch);
          if(iPlane == bestPlane){
            h_dEdx_BP->Fill(hitE/pitch);
            h_hit_widthX_BP->Fill(widthX);
            h_hit_widthTicks_BP->Fill(hit_width);
            h_dEdx_E_BP->Fill(evt->StartE_tpcAV[iG4],hitE/pitch);
          }
          
          if(hitsOnPlane.at(iPlane)/static_cast<double>(100*lengthAV) < 0.8) continue; // If the number of hits on this plane is silly w.r.t the length

          // Truth-level energy
          h_dEdx_hitCut_E.at(iPlane)->Fill(evt->StartE_tpcAV[iG4],hitE/pitch);
          h_dEdx_nDaught.at(iPlane)->Fill(evt->NumberDaughters[iG4],hitE/pitch);

          // Reco charge versus hit width, to compare with Gray's ICARUS studies
          if(iPlane == bestPlane){
            h_dQdx_hitWidth_BP->Fill(hit_width, hitQ/pitch);
            h_dEdx_hitCut_E_BP->Fill(evt->StartE_tpcAV[iG4],hitE/pitch);
            //h_dEdx_hitCut_E_BP->Fill(evt->StartE_tpcAV[iG4],hitdEdx);
            h_dEdx_hitCut_BP->Fill(hitE/pitch);
          }
        }// Hits
      } // Plane plane
    }// Geant
    h_true_mus->Fill(n_true_mus);
    h_true_primary_mus->Fill(n_true_primary_mus);

    //
    // Reco-level studies
    //
    for(unsigned int iTrk = 0; iTrk < nTrks; ++iTrk){
      
      // Check if we should be looking at cosmics (cry && origin(2)) or beam (NuWro && origin(4))
      if((*evt->trkorigin_pandoraTrack[iTrk] == 2 && !cry) || (*evt->trkorigin_pandoraTrack[iTrk] == 4 && cry) || *evt->trkorigin_pandoraTrack[iTrk] == -1) continue;

      // Get the track geometry
      TVector3 startVtx(evt->trkstartx_pandoraTrack[iTrk],
                        evt->trkstarty_pandoraTrack[iTrk],
                        evt->trkstartz_pandoraTrack[iTrk]);
      TVector3 endVtx(evt->trkendx_pandoraTrack[iTrk],
                      evt->trkendy_pandoraTrack[iTrk],
                      evt->trkendz_pandoraTrack[iTrk]);
      
      // Since we are generating downwards-going tracks - if the start y < end y then 
      // assume the reconstruction has got them the wrong way around and flip them
      if(startVtx.Y() < endVtx.Y()){
        TVector3 temp(endVtx);
        endVtx = startVtx;
        startVtx = temp;
      }

      // Get the best plane
      int bestPlane = 0;
      int currHits  = -999;
      std::vector<int> hitsOnPlane(3,0);
      for(int iPlane = 0; iPlane < 3; ++iPlane){
        hitsOnPlane.at(iPlane) = evt->ntrkhits_pandoraTrack[iTrk][iPlane];
        if(evt->ntrkhits_pandoraTrack[iTrk][iPlane] > currHits){
          currHits  = evt->ntrkhits_pandoraTrack[iTrk][iPlane];
          bestPlane = iPlane; 
        } // CurrHits
      } // Planes

      // Check for muon pdg code
      if(std::abs(evt->trkpdgtruth_pandoraTrack[iTrk][bestPlane]) != 13) continue;
    
      // Get the length to assess the geometry of the muons
      float length = evt->trklen_pandoraTrack[iTrk];

      // The following studies should be conducted with through-going muons to start with
      // If the number of external planes crossed is >= 2, the track is through-going
      bool throughGoing = IsThroughGoing(length,startVtx,endVtx,extPlanes,fidExtPlanes);
      
      // Fill the length histogram
      if(evtProc.SelectTrack(evt,iTrk)){
        h_reco_len->Fill(length/100.); // [m]
      }

      // Try to get the true energy
      // Get the list iterator from matching ID's
      float eng = -1.;
      for(unsigned int iG4 = 0; iG4 < nGeant; ++iG4){
        int trueID = evt->TrackId[iG4];

        if(evt->trkidtruth_pandoraTrack[iTrk][bestPlane] == trueID){
          eng = evt->StartE_tpcAV[iG4];
          break;
        }
      }
      if(eng < 0) continue;

      for(int iPlane = 0; iPlane < 3; ++iPlane){
        // Only look at long muons
        if(std::abs(evt->trkpdgtruth_pandoraTrack[iTrk][iPlane]) != 13) continue;

        // Make sure there is some information on the plane
        // Sometimes this is -999 and will break stuff
        unsigned int nHitsR = evt->ntrkhits_pandoraTrack[iTrk][iPlane];
        // We need at least 2 hits for our studies
        if(nHitsR < 2) continue;

        float energy = evt->trkke_pandoraTrack[iTrk][iPlane]/1000.; // GeV
        
        // Make sure the energy isn't stupid
        // Make sure it doesn't exceed the maximum energy of a true track: 1e5 GeV
        if(energy > 1e5) {
          nHugeE++;
          continue;
        }

        if(iPlane == bestPlane){
          // Reconstructed energy vs start y position
          h_reco_Y_E->Fill(energy,startVtx.Y());
          h_reco_Y_E_zoom->Fill(energy,startVtx.Y());
          h_reco_len_E->Fill(energy,length);
          h_reco_eng->Fill(energy);
          if(startVtx.Y() > 599.5)
            h_reco_eng_highy->Fill(energy);
        
          // If the current reconstructed track is associated to a true muon
          n_reco_mus++;
        } // if bestPlane
        
        if(!evtProc.SelectTrack(evt,iTrk)) continue;

        if(iPlane == bestPlane){
          // Now access the variables of interest
          h_reco_eng_long->Fill(energy);
          if(startVtx.Y() > 599.5){
            n_mu_reco++;
            total_energy_reco += energy;
            h_reco_eng_long_highy->Fill(energy);
            n_reco_long_highy_mus++;
          }
          
          // For long muon multiplicity plot
          n_reco_long_mus++;
        }
        
        if(thru != throughGoing) continue;

        // Make sure it doesn't exceed the maximum size of the array
        // Count if it does so we can see how often it happens
        if(nHitsR > MAX_TRACK_HITS){
          nHitsR = MAX_TRACK_HITS;
        }

        // Now access the variables of interest
        Float_t *RRArr   = evt->trkresrg_pandoraTrack[iTrk][iPlane];
        Float_t *dEdxArr = evt->trkdedx_pandoraTrack[iTrk][iPlane];
        Float_t *dQdxArr = evt->trkdqdx_pandoraTrack[iTrk][iPlane];

        // Convert them to vectors
        std::vector<float> RR(RRArr, RRArr + nHitsR);
        std::vector<float> dEdx(dEdxArr, dEdxArr + nHitsR);
        std::vector<float> dQdx(dQdxArr, dQdxArr + nHitsR);

        // And fill, removing first and last hit for completeness
        /* Not sure this is necessary
        // Set the hit counter to 0 for this plane
        unsigned int iHit = 0;
        for(unsigned int itHit = 0; itHit < nHits; ++itHit){
          // Once we have reached the number of reconstructed hits, break out of the loop
          if(iHit >= nHitsR) break;
          if(evt->hit_plane[itHit] != bestPlane) continue;
          */

        // Now loop over hits so we can work our calo magic
        for(unsigned int iHit = 0; iHit < nHitsR-1; ++iHit){

          
          // Get the location of the current and following hits to determine the pitch
          TVector3 trkXYZ(evt->trkxyz_pandoraTrack[iTrk][iPlane][iHit][0],
                          evt->trkxyz_pandoraTrack[iTrk][iPlane][iHit][1],
                          evt->trkxyz_pandoraTrack[iTrk][iPlane][iHit][2]);
          TVector3 nextXYZ(evt->trkxyz_pandoraTrack[iTrk][iPlane][iHit+1][0],
                           evt->trkxyz_pandoraTrack[iTrk][iPlane][iHit+1][1],
                           evt->trkxyz_pandoraTrack[iTrk][iPlane][iHit+1][2]);
          TVector3 endXYZ(evt->trkxyz_pandoraTrack[iTrk][iPlane][nHitsR-2][0],
                          evt->trkxyz_pandoraTrack[iTrk][iPlane][nHitsR-2][1],
                          evt->trkxyz_pandoraTrack[iTrk][iPlane][nHitsR-2][2]);

          float x = trkXYZ.X();
          float t = x * evtProc.kXtoT;
          double dp = GetHitPitch(iPlane, trkXYZ, nextXYZ);
          double cosDrift = GetCosDrift(trkXYZ, nextXYZ);

          // Look at the x position distribution on the best plane
          if(iPlane == bestPlane)
            h_xPos_BP->Fill(x);

          // Get the angle of the track to the y-direction
          TVector3 yDir(0,1,0);
          double thetaY = GetAngleToPlane(trkXYZ,nextXYZ,yDir);

          // Check if x is lower or higher than the APA bounds, charge seems to accumulate there
          // Put a 5 cm fiducial border around the APAs and CPAs too
          if(x < evtProc.PD_APA_X_POSITIONS[0]+5 || x > evtProc.PD_APA_X_POSITIONS[1]-5) continue;
          if(x > evtProc.PD_CPA_X_POSITIONS[0]-5 && x < evtProc.PD_CPA_X_POSITIONS[0]+5) continue;
          
          int tpc =evtProc.WhichTPC(x, true) + 1;
          float dx      = ( -1 + 2*(tpc%2) )*(x - evtProc.PD_APA_X_POSITIONS[tpc/2]);
          float dt      = dx*evtProc.kXtoT;
          float simtau  = evt->taulife/1.0e3; // convert from us to ms
          float simCorr = TMath::Exp(-dt/simtau); // Correction using the simulation value of the lifetime
          float corr    = TMath::Exp(-dt/etau.at(0));
          float eCorr   = TMath::Exp(-dt/etau.at(0)) / TMath::Exp(-dt/simtau); // Correct for the already-corrected energy
          if(LifetimePerTPC){
            corr  = TMath::Exp(-dt/etau.at(tpc-1));
            eCorr = TMath::Exp(-dt/etau.at(tpc-1)) / TMath::Exp(-dt/simtau);
          }

          // New values
          float RRVal       = RR.at(iHit)/100.;
          float dEdxVal     = dEdx.at(iHit);
          float dEdxCorr    = dEdxVal/eCorr;
          float dQdxVal     = dQdx.at(iHit)*5.; // Scale to 1 ADC = 200 e from 1 ADC = 1000 e
          float dQdxSimCorr = dQdxVal/simCorr;
          float dQdxCorr    = dQdxVal/corr;
          float hitWidth = evt->hit_endT[iHit] - evt->hit_startT[iHit];

          h_reco_dQdx_thetaY.at(iPlane)->Fill(thetaY,dQdxCorr);
          h_reco_dEdx_thetaY.at(iPlane)->Fill(thetaY,dEdxCorr);
          h_reco_dEdx_E.at(iPlane)->Fill(eng,dEdxCorr);
          if(hitsOnPlane.at(iPlane)/length > 0.2){
            h_reco_dEdx_E_hitCut.at(iPlane)->Fill(eng,dEdxCorr);
            h_reco_dEdx_thetaY_hitCut.at(iPlane)->Fill(thetaY,dEdxCorr);
          } 
          
          if(iPlane == bestPlane){
            h_reco_dQdx_simCorr_E->Fill(eng,dQdxSimCorr);
            h_reco_dQdx_pos->Fill(x,dQdxCorr);
            h_reco_dQdx_E->Fill(eng,dQdxCorr);
            h_reco_dQdx_RR->Fill(RRVal,dQdxCorr);
            h_reco_dQdx_dP->Fill(dp,dQdxCorr);
            h_reco_dQdx_width->Fill(hitWidth,dQdxCorr);
            h_reco_dQdx_cosDrift->Fill(cosDrift,dQdxCorr);
            h_reco_dEdx_noCorr_E_BP->Fill(eng,dEdxVal);
            h_reco_dEdx_pos_BP->Fill(x,dEdxCorr);
            h_reco_dEdx_E_BP->Fill(eng,dEdxCorr);
            h_reco_dEdx_RR_BP->Fill(RRVal,dEdxCorr);
            h_reco_dEdx_dP_BP->Fill(dp,dEdxCorr);
            h_reco_dEdx_dQdx_BP->Fill(dEdxCorr,dQdxCorr);
            h_hit_pitch->Fill(dp);
            // Now apply the minimum hits/length requirement
            if(hitsOnPlane.at(iPlane)/length > 0.8){
              h_reco_dEdx_E_hitCut_BP->Fill(eng,dEdxCorr);
            } 
          }
          // Increment the hit counter for this plane
          //iHit++;
        } // iHit
      } // iPlane
    } // iTrk

    h_reco_mus->Fill(n_reco_mus);
    h_reco_long_mus->Fill(n_reco_long_mus);
    h_reco_long_highy_mus->Fill(n_reco_long_highy_mus);
  }// Event loop
  std::cout << " --- 100 % --- |" << std::endl;

  std::vector<TH1D*> h_costheta_BPs{h_costheta_BP_0,h_costheta_BP_1,h_costheta_BP_2};
  std::vector<TH1D*> h_costheta_BP_hitCuts{h_costheta_BP_0_hitCut,h_costheta_BP_1_hitCut,h_costheta_BP_2_hitCut};

  TCanvas *c0 = new TCanvas("c0","",900,900);
  SetCanvasStyle(c0, 0.12,0.08,0.06,0.12,0,0,0);
  c0->SetLogy();

  TLegend *l = new TLegend(0.22,0.94,0.92,0.995);
  l->SetNColumns(2);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->SetTextFont(132);

  SetHistogramStyle1D(h_true_mus,"Muon multiplicity", "Rate");
  SetHistogramStyle1D(h_reco_mus,"Muon multiplicity", "Rate");

  // Now sort out the range
  double max_y = 1.1*std::max(h_true_mus->GetMaximum(),h_reco_mus->GetMaximum());

  // Now draw
  h_true_mus->Draw("hist");
  h_reco_mus->Draw("hist same");
  h_true_mus->SetLineWidth(3);
  h_true_mus->SetLineStyle(2);
  h_true_mus->SetLineColor(kTeal-5);
  h_reco_mus->SetLineStyle(2);
  h_reco_mus->SetLineWidth(3);
  h_reco_mus->SetLineColor(kViolet-5);
  h_true_mus->GetYaxis()->SetRangeUser(0.5,max_y);

  l->AddEntry(h_true_mus,"True muons", "l");
  l->AddEntry(h_reco_mus,"Reco muons", "l");
  l->Draw("same");

  c0->SaveAs((location+"/muon_multiplicities"+tag+".png").c_str());
  c0->SaveAs((location+"/muon_multiplicities"+tag+".root").c_str());
  c0->Clear();
  l->Clear();

  SetHistogramStyle1D(h_true_primary_mus,"Muon multiplicity", "Rate");
  SetHistogramStyle1D(h_reco_long_mus,"Muon multiplicity", "Rate");

  // Now sort out the range
  max_y = 1.1*std::max(h_true_primary_mus->GetMaximum(),h_reco_long_mus->GetMaximum());

  // Now draw
  h_true_primary_mus->Draw("hist");
  h_reco_long_mus->Draw("hist same");
  h_true_primary_mus->SetLineWidth(3);
  h_true_primary_mus->SetLineStyle(2);
  h_true_primary_mus->SetLineColor(kTeal-5);
  h_reco_long_mus->SetLineStyle(2);
  h_reco_long_mus->SetLineWidth(3);
  h_reco_long_mus->SetLineColor(kViolet-5);
  h_true_primary_mus->GetYaxis()->SetRangeUser(0.5,max_y);

  l->AddEntry(h_true_primary_mus,"True muons", "l");
  l->AddEntry(h_reco_long_mus,"Reco muons", "l");
  l->Draw("same");

  c0->SaveAs((location+"/long_muon_multiplicities"+tag+".png").c_str());
  c0->SaveAs((location+"/long_muon_multiplicities"+tag+".root").c_str());
  c0->Clear();
  l->Clear();

  SetHistogramStyle1D(h_true_primary_mus,"Muon multiplicity", "Rate");
  SetHistogramStyle1D(h_reco_long_highy_mus,"Muon multiplicity", "Rate");

  // Now sort out the range
  max_y = 1.1*std::max(h_true_primary_mus->GetMaximum(),h_reco_long_highy_mus->GetMaximum());

  // Now draw
  h_true_primary_mus->Draw("hist");
  h_reco_long_highy_mus->Draw("hist same");
  h_true_primary_mus->SetLineWidth(3);
  h_true_primary_mus->SetLineStyle(2);
  h_true_primary_mus->SetLineColor(kTeal-5);
  h_reco_long_highy_mus->SetLineStyle(2);
  h_reco_long_highy_mus->SetLineWidth(3);
  h_reco_long_highy_mus->SetLineColor(kViolet-5);
  h_true_primary_mus->GetYaxis()->SetRangeUser(0.5,max_y);

  l->AddEntry(h_true_primary_mus,"True muons", "l");
  l->AddEntry(h_reco_long_highy_mus,"Reco muons", "l");
  l->Draw("same");

  c0->SaveAs((location+"/long_highy_muon_multiplicities"+tag+".png").c_str());
  c0->SaveAs((location+"/long_highy_muon_multiplicities"+tag+".root").c_str());
  c0->Clear();
  l->Clear();

  c0->SetLogx();
  
  SetHistogramStyle1D(h_energy,"Muon energy [GeV]", "Rate/GeV");
  h_energy->Scale(1,"width");
  h_energy->Draw("hist");
  h_energy->SetLineWidth(3);
  h_energy->SetLineColor(kTeal-5);
  c0->SaveAs((location+"/energy"+tag+".png").c_str());
  c0->SaveAs((location+"/energy"+tag+".root").c_str());
  c0->Clear();

  SetHistogramStyle1D(h_energy_long,"Muon energy [GeV]", "Rate/GeV");
  h_energy_long->Scale(1,"width");
  h_energy_long->Draw("hist");
  h_energy_long->SetLineWidth(3);
  h_energy_long->SetLineColor(kTeal-5);
  c0->SaveAs((location+"/energy_long"+tag+".png").c_str());
  c0->SaveAs((location+"/energy_long"+tag+".root").c_str());
  c0->Clear();

  SetHistogramStyle1D(h_hitE_long_BP,"Total true muon energy depositions [GeV]", "Rate/GeV");
  h_hitE_long_BP->Scale(1,"width");
  h_hitE_long_BP->Draw("hist");
  h_hitE_long_BP->SetLineWidth(3);
  h_hitE_long_BP->SetLineColor(kTeal-5);
  c0->SaveAs((location+"/hitE_long_BP"+tag+".png").c_str());
  c0->SaveAs((location+"/hitE_long_BP"+tag+".root").c_str());
  c0->Clear();

  SetHistogramStyle1D(h_hitE_long_BP,"True muon energies [GeV]", "Rate/GeV");
  h_energy_long->Draw("hist");
  h_hitE_long_BP->Draw("hist same");
  h_hitE_long_BP->SetLineWidth(3);
  h_hitE_long_BP->SetLineColor(kPink+5);
  c0->SaveAs((location+"/energy_hitE_long_BP"+tag+".png").c_str());
  c0->SaveAs((location+"/energy_hitE_long_BP"+tag+".root").c_str());
  c0->Clear();

  TCanvas *c1 = new TCanvas("c1","",900,900);
  SetCanvasStyle(c1, 0.12,0.03,0.05,0.12,0,0,0);

  SetHistogramStyle1D(h_xPos_BP,"Deposition x position [cm]", "Rate");
  h_xPos_BP->Draw("hist");
  h_xPos_BP->SetLineWidth(3);
  h_xPos_BP->SetLineColor(kTeal-5);
  c1->SaveAs((location+"/xPos_BP"+tag+".png").c_str());
  c1->SaveAs((location+"/xPos_BP"+tag+".root").c_str());
  c1->Clear();

  SetHistogramStyle1D(h_hit_pitch,"Muon hit pitch [cm]", "Rate");
  h_hit_pitch->Draw("hist");
  h_hit_pitch->SetLineWidth(3);
  h_hit_pitch->SetLineColor(kTeal-5);
  c1->SaveAs((location+"/hit_pitch"+tag+".png").c_str());
  c1->SaveAs((location+"/hit_pitch"+tag+".root").c_str());
  c1->Clear();

  SetHistogramStyle1D(h_length,"Muon length [m]", "Rate");
  h_length->Draw("hist");
  h_length->SetLineWidth(3);
  h_length->SetLineColor(kTeal-5);
  c1->SaveAs((location+"/length"+tag+".png").c_str());
  c1->SaveAs((location+"/length"+tag+".root").c_str());
  c1->Clear();

  SetHistogramStyle1D(h_mom,"Muon momentum [GeV]", "Rate");
  h_mom->Draw("hist");
  h_mom->SetLineWidth(3);
  h_mom->SetLineColor(kTeal-5);
  c1->SaveAs((location+"/mom"+tag+".png").c_str());
  c1->SaveAs((location+"/mom"+tag+".root").c_str());
  c1->Clear();

  SetHistogramStyle1D(h_nDaughters,"Muon daughters", "Rate");
  h_nDaughters->Draw("hist");
  h_nDaughters->SetLineWidth(3);
  h_nDaughters->SetLineColor(kTeal-5);
  c1->SaveAs((location+"/nDaughters"+tag+".png").c_str());
  c1->SaveAs((location+"/nDaughters"+tag+".root").c_str());
  c1->Clear();
  
  SetHistogramStyle1D(h_reco_len,"Muon length [m]", "Rate");
  h_reco_len->Draw("hist");
  h_reco_len->SetLineWidth(3);
  h_reco_len->SetLineColor(kTeal-5);
  c1->SaveAs((location+"/reco_length"+tag+".png").c_str());
  c1->SaveAs((location+"/reco_length"+tag+".root").c_str());
  c1->Clear();
 
  SetHistogramStyle1D(h_energy_nolog,"Muon energy [GeV]", "Rate");
  h_energy_nolog->Draw("hist");
  h_energy_nolog->SetLineWidth(3);
  h_energy_nolog->SetLineColor(kTeal-5);
  c1->SaveAs((location+"/energy_nolog"+tag+".png").c_str());
  c1->SaveAs((location+"/energy_nolog"+tag+".root").c_str());
  c1->Clear();

  SetHistogramStyle1D(h_energy_long_nolog,"Muon energy [GeV]", "Rate");
  h_energy_long_nolog->Draw("hist");
  h_energy_long_nolog->SetLineWidth(3);
  h_energy_long_nolog->SetLineColor(kTeal-5);
  c1->SaveAs((location+"/energy_nolog_long"+tag+".png").c_str());
  c1->SaveAs((location+"/energy_nolog_long"+tag+".root").c_str());
  c1->Clear();

  for(unsigned int iPlane = 0; iPlane < 3; ++iPlane){
    SetHistogramStyle1D(h_costheta_BPs.at(iPlane),("cos#theta_{Plane,"+std::to_string(iPlane)+"}").c_str(),"Rate");
    h_costheta_BPs.at(iPlane)->Draw("hist");
    h_costheta_BPs.at(iPlane)->SetLineWidth(3);
    h_costheta_BPs.at(iPlane)->SetLineColor(pal.at(iPlane));
    c1->SaveAs((location+"/cosTheta_BP_"+std::to_string(iPlane)+tag+".png").c_str());
    c1->SaveAs((location+"/cosTheta_BP_"+std::to_string(iPlane)+tag+".root").c_str());
    c1->Clear();

    SetHistogramStyle1D(h_costheta_BP_hitCuts.at(iPlane),("cos#theta_{Plane,"+std::to_string(iPlane)+"}").c_str(),"Rate");
    h_costheta_BP_hitCuts.at(iPlane)->Draw("hist");
    h_costheta_BP_hitCuts.at(iPlane)->SetLineWidth(3);
    h_costheta_BP_hitCuts.at(iPlane)->SetLineColor(pal.at(iPlane));
    c1->SaveAs((location+"/cosTheta_BP_hitCut_"+std::to_string(iPlane)+tag+".png").c_str());
    c1->SaveAs((location+"/cosTheta_BP_hitCut_"+std::to_string(iPlane)+tag+".root").c_str());
    c1->Clear();

  }
  SetHistogramStyle1D(h_dEdx_BP,"True dE/dx [MeV/cm]", "Rate");
  h_dEdx_BP->Draw("hist");
  h_dEdx_BP->SetLineWidth(3);
  h_dEdx_BP->SetLineColor(kTeal-5);
  c1->SaveAs((location+"/dEdx_BP"+tag+".png").c_str());
  c1->SaveAs((location+"/dEdx_BP"+tag+".root").c_str());
  c1->Clear();

  SetHistogramStyle1D(h_EoL_BP,"True #Delta E / L [MeV/cm]", "Rate");
  h_EoL_BP->Draw("hist");
  h_EoL_BP->SetLineWidth(3);
  h_EoL_BP->SetLineColor(kTeal-5);
  c1->SaveAs((location+"/EoL_BP"+tag+".png").c_str());
  c1->SaveAs((location+"/EoL_BP"+tag+".root").c_str());
  c1->Clear();

  SetHistogramStyle1D(h_hit_widthX_BP,"Hit width [cm]", "Rate");
  h_hit_widthX_BP->Draw("hist");
  h_hit_widthX_BP->SetLineWidth(3);
  h_hit_widthX_BP->SetLineColor(kTeal-5);
  c1->SaveAs((location+"/hit_widthX_BP"+tag+".png").c_str());
  c1->SaveAs((location+"/hit_widthX_BP"+tag+".root").c_str());
  c1->Clear();

  SetHistogramStyle1D(h_hit_widthTicks_BP,"Hit width [ticks]", "Rate");
  h_hit_widthTicks_BP->Draw("hist");
  h_hit_widthTicks_BP->SetLineWidth(3);
  h_hit_widthTicks_BP->SetLineColor(kTeal-5);
  c1->SaveAs((location+"/hit_widthTicks_BP"+tag+".png").c_str());
  c1->SaveAs((location+"/hit_widthTicks_BP"+tag+".root").c_str());
  c1->Clear();

  SetHistogramStyle1D(h_dEdx_hitCut_BP,"True dE/dx [MeV/cm]", "Rate");
  h_dEdx_hitCut_BP->Draw("hist");
  h_dEdx_hitCut_BP->SetLineWidth(3);
  h_dEdx_hitCut_BP->SetLineColor(kTeal-5);
  c1->SaveAs((location+"/dEdx_hitCut_BP"+tag+".png").c_str());
  c1->SaveAs((location+"/dEdx_hitCut_BP"+tag+".root").c_str());
  c1->Clear();

  l->SetNColumns(1);
  l->SetX1NDC(0.347);
  l->SetY1NDC(0.710);
  l->SetX2NDC(0.817);
  l->SetY2NDC(0.957);
 
  // Now plot the overlay and legend
  h_costheta_BPs.at(0)->Draw("hist");
  double max_cos_y = -999.;
  for(unsigned int iPlane = 0; iPlane < 3; ++iPlane){
    h_costheta_BPs.at(iPlane)->Draw("hist same");
    h_costheta_BP_hitCuts.at(iPlane)->Draw("hist same");
    h_costheta_BPs.at(iPlane)->SetLineWidth(2);
    h_costheta_BPs.at(iPlane)->SetLineStyle(1);
    h_costheta_BP_hitCuts.at(iPlane)->SetLineWidth(3);
    h_costheta_BP_hitCuts.at(iPlane)->SetLineStyle(2);
    l->AddEntry(h_costheta_BPs.at(iPlane),("All events, plane "+std::to_string(iPlane)).c_str(), "l");
    l->AddEntry(h_costheta_BP_hitCuts.at(iPlane),("Events with < 0.2 hits/cm, plane "+std::to_string(iPlane)).c_str(), "l");

    double max_wire = std::max(h_costheta_BPs.at(iPlane)->GetMaximum(),h_costheta_BP_hitCuts.at(iPlane)->GetMaximum());
    if(max_wire > max_cos_y)
      max_cos_y = max_wire;
  }
  h_costheta_BPs.at(0)->GetYaxis()->SetRangeUser(0,1.1*max_cos_y);
  l->Draw("same");
  c1->SaveAs((location+"/cosTheta_BP_hitCut_overlay"+tag+".png").c_str());
  c1->SaveAs((location+"/cosTheta_BP_hitCut_overlay"+tag+".root").c_str());
  c1->Clear();
  l->Clear();

  SetHistogramStyle1D(h_nHitsPerL_BP,"Hits per unit length [cm^{-1}]", "Rate");
  h_nHitsPerL_BP->Draw("hist");
  h_nHitsPerL_BP->SetLineWidth(3);
  h_nHitsPerL_BP->SetLineColor(kTeal-5);
  c1->SaveAs((location+"/hits_per_L"+tag+".png").c_str());
  c1->SaveAs((location+"/hits_per_L"+tag+".root").c_str());
  c1->Clear();

  c1->SetLogx();
  c1->SetLogy();
  
  SetHistogramStyle1D(h_reco_eng,"Total muon deposition [GeV]", "Rate/GeV");
  h_reco_eng->Scale(1,"width");
  h_reco_eng->Draw("hist");
  h_reco_eng->SetLineWidth(3);
  h_reco_eng->SetLineColor(kTeal-5);
  c1->SaveAs((location+"/reco_energy"+tag+".png").c_str());
  c1->SaveAs((location+"/reco_energy"+tag+".root").c_str());
  c1->Clear();
  
  SetHistogramStyle1D(h_reco_eng_long,"Total muon deposition [GeV]", "Rate/GeV");
  h_reco_eng_long->Scale(1,"width");
  h_reco_eng_long->Draw("hist");
  h_reco_eng_long->SetLineWidth(3);
  h_reco_eng_long->SetLineColor(kTeal-5);
  c1->SaveAs((location+"/reco_energy_long"+tag+".png").c_str());
  c1->SaveAs((location+"/reco_energy_long"+tag+".root").c_str());
  c1->Clear();
  
  SetHistogramStyle1D(h_reco_eng_highy,"Total muon deposition [GeV]", "Rate/GeV");
  h_reco_eng_highy->Scale(1,"width");
  h_reco_eng_highy->Draw("hist");
  h_reco_eng_highy->SetLineWidth(3);
  h_reco_eng_highy->SetLineColor(kTeal-5);
  c1->SaveAs((location+"/reco_energy_highy"+tag+".png").c_str());
  c1->SaveAs((location+"/reco_energy_highy"+tag+".root").c_str());
  c1->Clear();
  
  SetHistogramStyle1D(h_reco_eng_long_highy,"Total muon deposition [GeV]", "Rate/GeV");
  h_reco_eng_long_highy->Scale(1,"width");
  h_reco_eng_long_highy->Draw("hist");
  h_reco_eng_long_highy->SetLineWidth(3);
  h_reco_eng_long_highy->SetLineColor(kTeal-5);
  c1->SaveAs((location+"/reco_energy_long_highy"+tag+".png").c_str());
  c1->SaveAs((location+"/reco_energy_long_highy"+tag+".root").c_str());
  c1->Clear();
  
  TCanvas *c2 = new TCanvas("c2","",1000,800);
  SetCanvasStyle(c2, 0.1,0.15,0.05,0.12,0,0,0);

  SetHistogramStyle2D(h_nHitsPerL_costheta_BP,"cos#theta_{Plane}", "Hits per unit length [cm^{-1}]",false);
  h_nHitsPerL_costheta_BP->Draw("colz");
  c2->SaveAs((location+"/cosTheta_hitsPerLength"+tag+".png").c_str());
  c2->SaveAs((location+"/cosTheta_hitsPerLength"+tag+".root").c_str());
  c2->Clear();

  SetHistogramStyle2D(h_nHitsPerL_costheta_CP,"cos#theta_{Collection}", "Hits per unit length [cm^{-1}]",false);
  h_nHitsPerL_costheta_CP->Draw("colz");
  c2->SaveAs((location+"/cosThetaCollection_hitsPerLength"+tag+".png").c_str());
  c2->SaveAs((location+"/cosThetaCollection_hitsPerLength"+tag+".root").c_str());
  c2->Clear();

  SetHistogramStyle2D(h_nHitsPerL_cos_to_plane,"cos#theta_{APA}", "Hits per unit length [cm^{-1}]",false);
  h_nHitsPerL_cos_to_plane->Draw("colz");
  c2->SaveAs((location+"/cosToPlane_hitsPerLength"+tag+".png").c_str());
  c2->SaveAs((location+"/cosToPlane_hitsPerLength"+tag+".root").c_str());
  c2->Clear();

  SetHistogramStyle2D(h_nHitsPerL_recoLength_BP,"Reconstructed length", "Hits per unit length [cm^{-1}]",false);
  h_nHitsPerL_recoLength_BP->Draw("colz");
  c2->SaveAs((location+"/recoLength_hitsPerLength"+tag+".png").c_str());
  c2->SaveAs((location+"/recoLength_hitsPerLength"+tag+".root").c_str());
  c2->Clear();

  SetHistogramStyle2D(h_E_nDaught,"Muon energy [GeV]", "Number of daughters",false);
  h_E_nDaught->Draw("colz");
  c2->SaveAs((location+"/nDaught_vs_E"+tag+".png").c_str());
  c2->SaveAs((location+"/nDaught_vs_E"+tag+".root").c_str());
  c2->Clear();

  SetHistogramStyle2D(h_dQdx_hitWidth_BP,"Hit width [ticks]"," Reconstructed dQ/dx [ADC/cm]",false);
  h_dQdx_hitWidth_BP->Draw("colz");
  c2->SaveAs((location+"/dQdx_vs_hitWidth"+tag+".png").c_str());
  c2->SaveAs((location+"/dQdx_vs_hitWidth"+tag+".root").c_str());
  c2->Clear();

  SetHistogramStyle2D(h_reco_dEdx_pos_BP,"x Position [cm]", "Reconstructed dE/dx [MeV/cm]",false);
  h_reco_dEdx_pos_BP->Draw("colz");
  c2->SaveAs((location+"/reco_dEdx_vs_pos_BP"+tag+".png").c_str());
  c2->SaveAs((location+"/reco_dEdx_vs_pos_BP"+tag+".root").c_str());
  c2->Clear();

  SetHistogramStyle2D(h_reco_dEdx_RR_BP,"Residual range [cm]", "Reconstructed dE/dx [MeV/cm]",false);
  h_reco_dEdx_RR_BP->Draw("colz");
  c2->SaveAs((location+"/reco_dEdx_vs_RR_BP"+tag+".png").c_str());
  c2->SaveAs((location+"/reco_dEdx_vs_RR_BP"+tag+".root").c_str());
  c2->Clear();

  SetHistogramStyle2D(h_reco_dEdx_dQdx_BP,"dE/dx [MeV/cm]","Reconstructed dQ/dx [ADC/cm]",false);
  h_reco_dEdx_dQdx_BP->Draw("colz");
  c2->SaveAs((location+"/reco_dEdx_vs_dQdx_BP"+tag+".png").c_str());
  c2->SaveAs((location+"/reco_dEdx_vs_dQdx_BP"+tag+".root").c_str());
  c2->Clear();

  SetHistogramStyle2D(h_reco_dQdx_pos,"x Position [cm]","Reconstructed dQ/dx [ADC/cm]",false);
  h_reco_dQdx_pos->Draw("colz");
  c2->SaveAs((location+"/reco_dQdx_vs_pos"+tag+".png").c_str());
  c2->SaveAs((location+"/reco_dQdx_vs_pos"+tag+".root").c_str());
  h_reco_dQdx_pos->SaveAs((location+"/hist_reco_dQdx_vs_pos"+tag+".root").c_str());
  c2->Clear();

  SetHistogramStyle2D(h_reco_dQdx_RR,"Residual range [m]","Reconstructed dQ/dx [ADC/cm]",false);
  h_reco_dQdx_RR->Draw("colz");
  c2->SaveAs((location+"/reco_dQdx_vs_RR"+tag+".png").c_str());
  c2->SaveAs((location+"/reco_dQdx_vs_RR"+tag+".root").c_str());
  h_reco_dQdx_RR->SaveAs((location+"/hist_reco_dQdx_vs_RR"+tag+".root").c_str());
  c2->Clear();

  SetHistogramStyle2D(h_reco_dQdx_dP,"Hit pitch [cm]","Reconstructed dQ/dx [ADC/cm]",false);
  h_reco_dQdx_dP->Draw("colz");
  c2->SaveAs((location+"/reco_dQdx_vs_dP"+tag+".png").c_str());
  c2->SaveAs((location+"/reco_dQdx_vs_dP"+tag+".root").c_str());
  h_reco_dQdx_dP->SaveAs((location+"/hist_reco_dQdx_vs_dP"+tag+".root").c_str());
  c2->Clear();

  SetHistogramStyle2D(h_reco_dQdx_width,"Hit width [ticks]","Reconstructed dQ/dx [ADC/cm]",false);
  h_reco_dQdx_width->Draw("colz");
  c2->SaveAs((location+"/reco_dQdx_vs_hitWidth"+tag+".png").c_str());
  c2->SaveAs((location+"/reco_dQdx_vs_hitWidth"+tag+".root").c_str());
  h_reco_dQdx_width->SaveAs((location+"/hist_reco_dQdx_vs_hitWidth"+tag+".root").c_str());
  c2->Clear();

  SetHistogramStyle2D(h_reco_dQdx_cosDrift,"cos #theta_{Drift}","Reconstructed dQ/dx [ADC/cm]",false);
  h_reco_dQdx_cosDrift->Draw("colz");
  c2->SaveAs((location+"/reco_dQdx_vs_cosDrift"+tag+".png").c_str());
  c2->SaveAs((location+"/reco_dQdx_vs_cosDrift"+tag+".root").c_str());
  h_reco_dQdx_cosDrift->SaveAs((location+"/hist_reco_dQdx_vs_cosDrift"+tag+".root").c_str());
  c2->Clear();

  for(unsigned int iPlane = 0; iPlane < 3; ++iPlane){
    
    SetHistogramStyle2D(h_dEdx_nDaught.at(iPlane),"Muon daughters", "dE/dx [MeV/cm]",false);
    h_dEdx_nDaught.at(iPlane)->Draw("colz");
    c2->SaveAs((location+"/dEdx_vs_nDaughters"+std::to_string(iPlane)+tag+".png").c_str());
    c2->SaveAs((location+"/dEdx_vs_nDaughters"+std::to_string(iPlane)+tag+".root").c_str());
    c2->Clear();
    
    SetHistogramStyle2D(h_reco_dQdx_thetaY.at(iPlane),"cos#theta_{Y}","dQ/dx [ADC/cm]", false);
    h_reco_dQdx_thetaY.at(iPlane)->Draw("colz");
    c2->SaveAs((location+"/reco_dQdx_costhetaY_"+std::to_string(iPlane)+tag+".png").c_str());
    c2->SaveAs((location+"/reco_dQdx_costhetaY_"+std::to_string(iPlane)+tag+".root").c_str());
    c2->Clear();

    SetHistogramStyle2D(h_reco_dEdx_thetaY.at(iPlane),"cos#theta_{Y}","dE/dx [MeV/cm]", false);
    h_reco_dEdx_thetaY.at(iPlane)->Draw("colz");
    c2->SaveAs((location+"/reco_dEdx_costhetaY_"+std::to_string(iPlane)+tag+".png").c_str());
    c2->SaveAs((location+"/reco_dEdx_costhetaY_"+std::to_string(iPlane)+tag+".root").c_str());
    c2->Clear();

    SetHistogramStyle2D(h_reco_dEdx_thetaY_hitCut.at(iPlane),"cos#theta_{Y}","dE/dx [MeV/cm]", false);
    h_reco_dEdx_thetaY_hitCut.at(iPlane)->Draw("colz");
    c2->SaveAs((location+"/reco_dEdx_costhetaY_hitCut_"+std::to_string(iPlane)+tag+".png").c_str());
    c2->SaveAs((location+"/reco_dEdx_costhetaY_hitCut_"+std::to_string(iPlane)+tag+".root").c_str());
    c2->Clear();


  } // Plane planes

  TCanvas *c3 = new TCanvas("c3","",1000,800);
  SetCanvasStyle(c3, 0.1,0.15,0.05,0.12,0,0,0);
  c3->SetLogx();
  
  SetHistogramStyle2D(h_reco_Y_E,"Total muon deposition [GeV]", "Muon start Y Position [cm]",false);
  h_reco_Y_E->Scale(1,"width");
  h_reco_Y_E->Draw("colz");
  c3->SaveAs((location+"/reco_Y_vs_E"+tag+".png").c_str());
  c3->SaveAs((location+"/reco_Y_vs_E"+tag+".root").c_str());
  c3->Clear();

  SetHistogramStyle2D(h_reco_Y_E_zoom,"Total muon deposition [GeV]", "Muon start Y Position [cm]",false);
  h_reco_Y_E_zoom->Scale(1,"width");
  h_reco_Y_E_zoom->Draw("colz");
  c3->SaveAs((location+"/reco_Y_vs_E_zoom"+tag+".png").c_str());
  c3->SaveAs((location+"/reco_Y_vs_E_zoom"+tag+".root").c_str());
  c3->Clear();

  SetHistogramStyle2D(h_reco_len_E,"Total muon deposition [GeV]", "Muon track length [m]",false);
  h_reco_len_E->Scale(1,"width");
  h_reco_len_E->Draw("colz");
  c3->SaveAs((location+"/reco_len_vs_E"+tag+".png").c_str());
  c3->SaveAs((location+"/reco_len_vs_E"+tag+".root").c_str());
  c3->Clear();

  for(unsigned int iPlane = 0; iPlane < 3; ++iPlane){
  
    SetHistogramStyle2D(h_dEdx_E.at(iPlane),"Muon energy [GeV]", "dE/dx [MeV/cm]",false);
    h_dEdx_E.at(iPlane)->Scale(1,"width");
    h_dEdx_E.at(iPlane)->Draw("colz");
    c3->SaveAs((location+"/dEdx_vs_E"+std::to_string(iPlane)+tag+".png").c_str());
    c3->SaveAs((location+"/dEdx_vs_E"+std::to_string(iPlane)+tag+".root").c_str());
    c3->Clear();
    
    SetHistogramStyle2D(h_dEdx_hitCut_E.at(iPlane),"Muon energy [GeV]", "dE/dx [MeV/cm]",false);
    h_dEdx_hitCut_E.at(iPlane)->Scale(1,"width");
    h_dEdx_hitCut_E.at(iPlane)->Draw("colz");
    c3->SaveAs((location+"/dEdx_hitCut_vs_E"+std::to_string(iPlane)+tag+".png").c_str());
    c3->SaveAs((location+"/dEdx_hitCut_vs_E"+std::to_string(iPlane)+tag+".root").c_str());
    c3->Clear();
    
    SetHistogramStyle2D(h_reco_dEdx_E.at(iPlane),"Total muon deposition [GeV]", "Reconstructed dE/dx [MeV/cm]",false);
    h_reco_dEdx_E.at(iPlane)->Scale(1,"width");
    h_reco_dEdx_E.at(iPlane)->Draw("colz");
    c3->SaveAs((location+"/reco_dEdx_vs_E"+std::to_string(iPlane)+tag+".png").c_str());
    c3->SaveAs((location+"/reco_dEdx_vs_E"+std::to_string(iPlane)+tag+".root").c_str());
    c3->Clear();

    SetHistogramStyle2D(h_reco_dEdx_E_hitCut.at(iPlane),"Total muon deposition [GeV]", "Reconstructed dE/dx [MeV/cm]",false);
    h_reco_dEdx_E_hitCut.at(iPlane)->Scale(1,"width");
    h_reco_dEdx_E_hitCut.at(iPlane)->Draw("colz");
    c3->SaveAs((location+"/reco_dEdx_vs_E_hitCut"+std::to_string(iPlane)+tag+".png").c_str());
    c3->SaveAs((location+"/reco_dEdx_vs_E_hitCut"+std::to_string(iPlane)+tag+".root").c_str());
    c3->Clear();

  } // Plane planes
  SetHistogramStyle2D(h_dEdx_E_BP,"Muon energy [GeV]", "dE/dx [MeV/cm]",false);
  h_dEdx_E_BP->Scale(1,"width");
  h_dEdx_E_BP->Draw("colz");
  c3->SaveAs((location+"/dEdx_vs_E_BP"+tag+".png").c_str());
  c3->SaveAs((location+"/dEdx_vs_E_BP"+tag+".root").c_str());
  c3->Clear();

  SetHistogramStyle2D(h_dEdx_hitCut_E_BP,"Muon energy [GeV]", "dE/dx [MeV/cm]",false);
  h_dEdx_hitCut_E_BP->Scale(1,"width");
  h_dEdx_hitCut_E_BP->Draw("colz");
  c3->SaveAs((location+"/dEdx_hitCut_vs_E_BP"+tag+".png").c_str());
  c3->SaveAs((location+"/dEdx_hitCut_vs_E_BP"+tag+".root").c_str());
  c3->Clear();
    
  SetHistogramStyle2D(h_reco_dQdx_simCorr_E,"True muon energy [GeV]", "Reconstructed dQ/dx [ADC/cm]",false);
  h_reco_dQdx_simCorr_E->Draw("colz");
  c3->SaveAs((location+"/reco_dQdx_simCorr_vs_E_BP_noScale"+tag+".png").c_str());
  c3->SaveAs((location+"/reco_dQdx_simCorr_vs_E_BP_noScale"+tag+".root").c_str());
  h_reco_dQdx_simCorr_E->SaveAs((location+"/hist_reco_dQdx_simCorr_vs_E_BP_noScale"+tag+".root").c_str());
  c3->Clear();

  SetHistogramStyle2D(h_reco_dQdx_E,"True muon energy [GeV]", "Reconstructed dQ/dx [ADC/cm]",false);
  h_reco_dQdx_E->Draw("colz");
  c3->SaveAs((location+"/reco_dQdx_vs_E_BP_noScale"+tag+".png").c_str());
  c3->SaveAs((location+"/reco_dQdx_vs_E_BP_noScale"+tag+".root").c_str());
  h_reco_dQdx_E->SaveAs((location+"/hist_reco_dQdx_vs_E_BP_noScale"+tag+".root").c_str());
  c3->Clear();

  SetHistogramStyle2D(h_reco_dQdx_simCorr_E,"True muon energy [GeV]", "Reconstructed dQ/dx [ADC/cm]",false);
  h_reco_dQdx_simCorr_E->Scale(1,"width");
  h_reco_dQdx_simCorr_E->Draw("colz");
  c3->SaveAs((location+"/reco_dQdx_simCorr_vs_E_BP"+tag+".png").c_str());
  c3->SaveAs((location+"/reco_dQdx_simCorr_vs_E_BP"+tag+".root").c_str());
  h_reco_dQdx_simCorr_E->SaveAs((location+"/hist_reco_dQdx_simCorr_vs_E_BP"+tag+".root").c_str());
  c3->Clear();

  SetHistogramStyle2D(h_reco_dQdx_E,"True muon energy [GeV]", "Reconstructed dQ/dx [ADC/cm]",false);
  h_reco_dQdx_E->Scale(1,"width");
  h_reco_dQdx_E->Draw("colz");
  c3->SaveAs((location+"/reco_dQdx_vs_E_BP"+tag+".png").c_str());
  c3->SaveAs((location+"/reco_dQdx_vs_E_BP"+tag+".root").c_str());
  h_reco_dQdx_E->SaveAs((location+"/hist_reco_dQdx_vs_E_BP"+tag+".root").c_str());
  c3->Clear();

  SetHistogramStyle2D(h_reco_dEdx_noCorr_E_BP,"True muon energy [GeV]", "Reconstructed dE/dx [MeV/cm]",false);
  h_reco_dEdx_noCorr_E_BP->Scale(1,"width");
  h_reco_dEdx_noCorr_E_BP->Draw("colz");
  c3->SaveAs((location+"/reco_dEdx_noCorr_vs_E_BP"+tag+".png").c_str());
  c3->SaveAs((location+"/reco_dEdx_noCorr_vs_E_BP"+tag+".root").c_str());
  c3->Clear();

  SetHistogramStyle2D(h_reco_dEdx_E_BP,"True muon energy [GeV]", "Reconstructed dE/dx [MeV/cm]",false);
  h_reco_dEdx_E_BP->Scale(1,"width");
  h_reco_dEdx_E_BP->Draw("colz");
  c3->SaveAs((location+"/reco_dEdx_vs_E_BP"+tag+".png").c_str());
  c3->SaveAs((location+"/reco_dEdx_vs_E_BP"+tag+".root").c_str());
  c3->Clear();

  SetHistogramStyle2D(h_reco_dEdx_E_hitCut_BP,"True muon energy [GeV]", "Reconstructed dE/dx [MeV/cm]",false);
  h_reco_dEdx_E_hitCut_BP->Scale(1,"width");
  h_reco_dEdx_E_hitCut_BP->Draw("colz");
  c3->SaveAs((location+"/reco_dEdx_vs_E_hitCut_BP"+tag+".png").c_str());
  c3->SaveAs((location+"/reco_dEdx_vs_E_hitCut_BP"+tag+".root").c_str());
  c3->Clear();

  // Now calculate averages from muons passing all criteria
  double average_energy_true = total_energy_true / static_cast<double>(n_mu_true);
  double average_energy_reco = total_energy_reco / static_cast<double>(n_mu_reco);

  // Now find the peak
  double peak_energy_true = GetPeakBinCentre(h_energy_long);
  double peak_energy_reco = GetPeakBinCentre(h_reco_eng_long_highy);

  ofstream oFile;
  oFile.open((location+"/statistics"+tag+".txt").c_str());

  oFile << " Total muons true              : " << n_mu_true << ", reco E < 1e5 GeV : " << n_mu_reco << ", reco E > 1e5 GeV : " << nHugeE << std::endl;
  oFile << " Total energy true             : " << total_energy_true << ", reco: " << total_energy_reco << " GeV " << std::endl;
  oFile << " Average energy true           : " << average_energy_true << ", reco: " << average_energy_reco << " GeV " << std::endl;
  oFile << " Peak energy true              : " << peak_energy_true << ", reco: " << peak_energy_reco << " GeV " << std::endl;
  oFile << " Total through-going muons     : " << n_mu_thru << std::endl;
  oFile << " Frequency of each plane as BP : " << std::setw(10) << " 0" << std::setw(10) << " 1" << std::setw(10) << " 2" << std::endl;
  oFile << "                               : " << std::setw(10) << nBP_0 << std::setw(10) << nBP_1 << std::setw(10) << nBP_2 << std::endl;

  oFile.close();

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
