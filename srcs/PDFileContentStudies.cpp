/************************************************************************
 *  
 * A macro to plot various reconstructed quantities in the CR muon 
 * calibration studies
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
   "genie_primaries_pdg",
   "inTPCActive",
   "pdg",
   "origin",
   "TrackId",
   "Mother",
   "Eng",
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
   "process_primary",
   "trkpdgtruth_pandoraTrack",
   "trkg4id_pandoraTrack",
   "trkorigin_pandoraTrack",
   "trkId_pandoraTrack",
   "trkidtruth_pandoraTrack",
   "ntracks_pandoraTrack",
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
   "trklen_pandoraTrack"
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
     
int pdFileContentStudies(const char *config){

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
  int n               = -1;
  int yCut            = 0;
  int thru            = 0;
  int cry             = 0; // Whether or not to save the cosmic (cry) information (trkorigin = 2) or the beam information (trkorigin = 4)
  int drawPFParticles = 0; // Whether or not to draw the PFParticles
  
  int nBinsFromPeak   = -1; // How many bins to traverse either side of the peak in the fit
  int nBinsFromPeakL  = -1; // How many bins to traverse left of the peak in the fit
  int nBinsFromPeakR  = -1; // How many bins to traverse right of the peak in the fit

  std::vector<double> etau; // measured electron lifetime, one per TPC if desired

  std::string input_list = "";
  std::string location="";
  std::string tag="";
  std::vector<double> minx_fid, miny_fid, minz_fid;
  std::vector<double> maxx_fid, maxy_fid, maxz_fid;
  std::vector<double> minx_av, miny_av, minz_av;
  std::vector<double> maxx_av, maxy_av, maxz_av;

  p->getValue("InputList",       input_list);
  p->getValue("Location",        location);
  p->getValue("Tag",             tag);
  p->getValue("NFiles",          n);
  p->getValue("YCut",            yCut);
  p->getValue("Thru",            thru);
  p->getValue("Cry",             cry);
  p->getValue("ETau",            etau);
  p->getValue("DrawPFParticles", drawPFParticles);
  p->getValue("NBinsFromPeak",   nBinsFromPeak);
  p->getValue("NBinsFromPeakL",  nBinsFromPeakL);
  p->getValue("NBinsFromPeakR",  nBinsFromPeakR);
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

  SortBinsFromPeak(nBinsFromPeak,nBinsFromPeakL,nBinsFromPeakR);

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
  std::cout << " Number of files: " << n << std::endl;
  
  // Start of analysis (loop over chain and events
  std::cout << " Running analysis..." << std::endl;

  // Then setup the histograms, counters and any other variables to add to
  // Setup histograms
  TH2D *h_dedx_x            = new TH2D("h_dedx_x","",50,-400,400,50,0.2,6);
  TH2D *h_dqdx_x            = new TH2D("h_dqdx_x","",50,-400,400,50,80,600);
  TH2D *h_dqdx_t            = new TH2D("h_dqdx_t","",50,-2.49,2.49,50,80,600);
  TH2D *h_dqdx_RR_stop      = new TH2D("h_dqdx_RR_stop","",50,0,3,50,80,600);
  TH2D *h_dqdx_RR           = new TH2D("h_dqdx_RR","",50,0,6,50,80,600);
  TH2D *h_dqdx_E            = new TH2D("h_dqdx_E","",50,0.8,40,50,80,600);
  TH2D *h_corr_dedx_x       = new TH2D("h_corr_dedx_x","",50,-400,400,50,0.2,6);
  TH2D *h_corr_dqdx_x       = new TH2D("h_corr_dqdx_x","",50,-400,400,50,80,600);
  TH2D *h_corr_dqdx_t       = new TH2D("h_corr_dqdx_t","",50,-2.49,2.49,50,80,600);
  TH2D *h_corr_dqdx_E       = new TH2D("h_corr_dqdx_E","",50,0.8,40,50,80,600);
  TH2D *h_corr_dqdx_RR_stop = new TH2D("h_corr_dqdx_RR_stop","",50,0,3,50,80,600);
  TH2D *h_corr_dqdx_RR      = new TH2D("h_corr_dqdx_RR","",50,0,6,50,80,600);
  TH2D *h_corr_dedq_x       = new TH2D("h_corr_dedq_x","",50,-400,400,50,5.8e-3,7.5e-3);
  TH2D *h_corr2_dedq_x      = new TH2D("h_corr2_dedq_x","",50,-400,400,50,5.8e-3,7.5e-3);
  TH2D *h_hits_xy           = new TH2D("h_hits_xy","",50,-400,400,50,0,610);
  TH2D *h_hits_xz           = new TH2D("h_hits_xz","",50,-400,400,300,-1,700);
  TH2D *h_hits_yz           = new TH2D("h_hits_yz","",50,0,610,300,-1,700);
  TH3D *h_hits_xyz          = new TH3D("h_hits_xyz","",50,-400,400,50,0,610,50,-1,700);
  TH2D *h_dqdx_xy           = new TH2D("h_dqdx_xy","",50,-400,400,65,0,610);
  TH2D *h_dqdx_xz           = new TH2D("h_dqdx_xz","",50,-400,400,310,-1,700);
  TH2D *h_dqdx_yz           = new TH2D("h_dqdx_yz","",65,0,610,310,-1,700);
  TH2D *h_statErr_xy        = new TH2D("h_statErr_xy","",50,-400,400,65,0,610);
  TH2D *h_statErr_xz        = new TH2D("h_statErr_xz","",50,-400,400,310,-1,700);
  TH2D *h_statErr_yz        = new TH2D("h_statErr_yz","",65,0,610,310,-1,700);
  TH1D *h_plane_cross       = new TH1D("h_plane_cross","",7,0,7); // Number of tracks crossing each plane
  TH1D *h_plane_enter       = new TH1D("h_plane_enter","",7,0,7); // Number of tracks entering from each external plane
  TH1D *h_plane_exit        = new TH1D("h_plane_exit","",7,0,7); // Number of tracks exiting from each external plane
  TH1D *h_enter_dist        = new TH1D("h_enter_dist","",50,0,100); // Number of tracks entering from each external plane
  TH1D *h_exit_dist         = new TH1D("h_exit_dist","",50,0,100); // Number of tracks entering from each external plane
  TH1D *h_muon_length       = new TH1D("h_muon_length","",50,0,12); // Muon length
  TH1D *h_n_crossed         = new TH1D("h_n_crossed","",4,0,4); // Number of planes crossed by each track
  TH1D *h_median_dqdx_xy    = new TH1D("h_median_dqdx_xy","XY Plane",40,80,600);
  TH1D *h_median_dqdx_xz    = new TH1D("h_median_dqdx_xz","XZ Plane",40,80,600);
  TH1D *h_median_dqdx_yz    = new TH1D("h_median_dqdx_yz","YZ Plane",40,80,600);
  TH1D *h_startX            = new TH1D("h_startX","",50,-400,400);  // Start X position of the muons
  TH1D *h_startY            = new TH1D("h_startY","",50,0,610);  // Start Y position of the muons
  TH1D *h_startZ            = new TH1D("h_startZ","",50,-1,700);  // Start Z position of the muons
  TH1D *h_endX              = new TH1D("h_endX","",50,-400,400);  // end X position of the muons
  TH1D *h_endY              = new TH1D("h_endY","",50,0,610);  // end Y position of the muons
  TH1D *h_endZ              = new TH1D("h_endZ","",50,-1,700);  // end Z position of the muons
  TH1D *h_depX              = new TH1D("h_depX","",50,-400,400);  // Start X position of the muons

  // Map of global bin numbers and vector of dQ/dx entries for each dimensional histogram
  // Use these to calculate the median dQ/dx in each bin and fill the median maps
  // Can then assess the distribution of medians and establish the variation in depositions in space
  std::map<int,std::vector<float>> bin_dqdx_xy, bin_dqdx_xz, bin_dqdx_yz;
  std::map<int,float> median_dqdx_xy, median_dqdx_xz, median_dqdx_yz;

  // Setup a vector of histograms and maps for filling
  std::vector<std::map<int,std::vector<float>>> bin_dqdx_maps{bin_dqdx_xy, bin_dqdx_xz, bin_dqdx_yz};
  std::vector<std::map<int,float>> median_dqdx_maps{median_dqdx_xy, median_dqdx_xz, median_dqdx_yz};
  std::vector<TH2D*> dqdx_hists{h_dqdx_xy,h_dqdx_xz,h_dqdx_yz};
  std::vector<TH2D*> error_hists{h_statErr_xy,h_statErr_xz,h_statErr_yz};
  std::vector<TH1D*> median_hists{h_median_dqdx_xy,h_median_dqdx_xz,h_median_dqdx_yz};

  // Setup vectors of PFParticles
  std::vector<TVector3> startPoints, endPoints;

  // Setup the maps
  for(unsigned int h = 0; h < dqdx_hists.size(); ++h){
    TH2D *hist = dqdx_hists.at(h);
    for(int xBin = 0; xBin <= hist->GetNbinsX(); ++xBin){
      for(int yBin = 0; yBin <= hist->GetNbinsY(); ++yBin){
        std::vector<float> dqdx;
        int global = hist->GetBin(xBin,yBin);
        bin_dqdx_maps.at(h).emplace(global,dqdx);
      } // ybin
    } // Xbin
  } // Hists

  // Sort out log scales if needed 
  SetLogX(h_dqdx_E);
  SetLogX(h_corr_dqdx_E);
  
  // Setup counters
  unsigned int maxHitsLimit     = 0;
  unsigned int wrongWay         = 0;
  unsigned int noPlanes         = 0;
  unsigned int totalTracks      = 0;
  unsigned int nMu              = 0;
  unsigned int nLongTracks      = 0;
  unsigned int topBottom        = 0;
  unsigned int topOrBottom      = 0;
  unsigned int frontBack        = 0;
  unsigned int frontOrBack      = 0;
  unsigned int min2APACPA       = 0;
  unsigned int eq1APACPA        = 0;
  unsigned int stopping         = 0;
  unsigned int exiting          = 0;
  unsigned int outOfRange       = 0;

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
    unsigned int nHits  = evt->no_hits_stored;
    unsigned int nGeant = evt->geant_list_size;
    
    // Print the processing rate
    double evtFrac  = iEvt/static_cast<double>(nEvts);

    if((std::abs(0.1*iIt)-evtFrac) < std::numeric_limits<double>::epsilon()){
      std::cout << " --- " << evtFrac*100 << " %";
      std::cout.flush();
      iIt++;
    }

    // Now loop over the tracks so we can do some stuff!!
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
      //
      // Temporarily remove these tracks for energy-deposition's sake
      if(startVtx.Y() < endVtx.Y()){
        wrongWay++;
        TVector3 temp(endVtx);
        endVtx = startVtx;
        startVtx = temp;
      }

      // Setup list of plane labels the track has crossed
      std::vector<std::string> labelsCrossed;

      // Get the best plane
      int bestPlane = -1;
      int currHits  = -999;
      for(int iPlane = 0; iPlane < 3; ++iPlane){
        if(evt->ntrkhits_pandoraTrack[iTrk][iPlane] > currHits){
          currHits  = evt->ntrkhits_pandoraTrack[iTrk][iPlane];
          bestPlane = iPlane; 
        } // CurrHits
      } // Planes

      // Count tracks
      totalTracks++;

      int ID = evt->trkId_pandoraTrack[iTrk];

      // Only look at primary muons
      if(std::abs(evt->trkpdgtruth_pandoraTrack[iTrk][bestPlane]) != 13) continue;
      nMu++;
     
      // Length cuts (2m)
      if(!evtProc.SelectTrack(evt,iTrk)) continue;
      nLongTracks++;
      
      float length = evt->trklen_pandoraTrack[iTrk];
      h_muon_length->Fill(length/100.);

      // Fill the start positions to make sure everything looks sensible
      h_startX->Fill(startVtx.X());
      h_startY->Fill(startVtx.Y());
      h_startZ->Fill(startVtx.Z());
      h_endX->Fill(endVtx.X());
      h_endY->Fill(endVtx.Y());
      h_endZ->Fill(endVtx.Z());

      // Fill the start and end point vectors for drawing the pfparticles
      startPoints.push_back(startVtx);
      endPoints.push_back(endVtx);

      // Find the closest plane to the start vertex and count it as a crossing plane
      Plane enteringPlane = GetClosestPlane(extPlanes, startVtx, endVtx);
      double distFromEntrance = GetDistanceToPlane(enteringPlane, startVtx, endVtx);
      h_enter_dist->Fill(distFromEntrance);

      Plane exitingPlane = GetClosestPlane(extPlanes, endVtx, startVtx);
      double distFromExit = GetDistanceToPlane(exitingPlane, endVtx, startVtx);
      h_exit_dist->Fill(distFromExit);

      // Now determine if the current track crossed each detector plane individually
      unsigned int planeN    = 0;
      unsigned int extPlaneN = 0;
      
      // Counter for the number of planes this track has crossed
      unsigned int nPlanesCrossed = 0;
      unsigned int nExtCrossed    = 0;

      // Loop over planes
      for(const Plane &pl : fidAllPlanes){
        if(planeN > fidAllPlanes.size()){
          std::cerr << " Error: Somehow the current plane iterator exceeds the number of planes in the detector: " << std::endl;
          std::cerr << " Iterator: " << planeN << " of " << allPlanes.size() << " total possible planes " << std::endl;
          std::exit(1);
        } // Debug
        // Check if this is the plane it (likely) entered the detector through 
        // Determine a maximum allowed distance from the plane to count it as the entrance
        if(enteringPlane.GetLabel() == pl.GetLabel()){
          if(distFromEntrance < 1){
            h_plane_cross->Fill(planeN);
            h_plane_enter->Fill(planeN);
            nPlanesCrossed++;
            labelsCrossed.push_back(pl.GetLabel());
          }
        }
        else if(exitingPlane.GetLabel() == pl.GetLabel()){
          if(distFromExit < 1){
            h_plane_cross->Fill(planeN);
            h_plane_exit->Fill(planeN);
            nPlanesCrossed++;
            labelsCrossed.push_back(pl.GetLabel());
          }
        }
        // Otherwise, check if the track intersects the current plane
        else if(CheckIfIntersectsPlane(pl,startVtx,endVtx,length)){
          h_plane_cross->Fill(planeN);
          nPlanesCrossed++;
          labelsCrossed.push_back(pl.GetLabel());
        } // Intersects
        // Sort out the bin label
        h_plane_cross->GetXaxis()->SetBinLabel(planeN+1,planeLabels.find(pl.GetLabel())->second.c_str());
        h_plane_enter->GetXaxis()->SetBinLabel(planeN+1,planeLabels.find(pl.GetLabel())->second.c_str());
        h_plane_exit->GetXaxis()->SetBinLabel(planeN+1,planeLabels.find(pl.GetLabel())->second.c_str());
        
        planeN++;
      } // Planes
      
      // Now fill the number of planes crossed histogram
      h_n_crossed->Fill(nPlanesCrossed);
      if(nPlanesCrossed == 0){
        noPlanes++;
      }
      
      // Now count the planes crossed for the stats table
      bool isStopping = IsStopping(length,startVtx,endVtx,extPlanes,fidExtPlanes);
      if(isStopping) stopping++;

      bool thruGoing = IsThroughGoing(length,startVtx,endVtx,extPlanes,fidExtPlanes);
      if(thruGoing) exiting++;

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

      // Now fill dQ/dx and dE/dx and hit histograms for each of the three wire planes
      // Somehow flag the best wire plane histogram
      // Look at everything on the best plane
      unsigned int nHitsR = evt->ntrkhits_pandoraTrack[iTrk][bestPlane];

      // Get the associated true energy of the muon
      // Try to get the true energy
      // Get the list iterator from matching ID's
      float eng = -1.;
      bool isTrueStopping = false;
      for(unsigned int iG4 = 0; iG4 < nGeant; ++iG4){
        int trueID = evt->TrackId[iG4];

        if(evt->trkidtruth_pandoraTrack[iTrk][bestPlane] == trueID){
          eng = evt->Eng[iG4];

          // Check for stopping muons in truth
          TVector3 endAV(evt->EndPointx_tpcAV[iG4],evt->EndPointy_tpcAV[iG4],evt->EndPointz_tpcAV[iG4]);
          float dx = abs(endAV.X()-evt->EndPointx[iG4]);
          float dy = abs(endAV.Y()-evt->EndPointy[iG4]);
          float dz = abs(endAV.Z()-evt->EndPointz[iG4]);

          // If these match, the TPC end point and general end point are the same, therefore the particle stops
          if(dx+dy+dz < 1e-10) isTrueStopping = true;

          break;
        }
      }
      if(eng < 0){
        std::cout << " Warning: Energy is below zero, skipping track with energy: " << eng << std::endl;
        continue;
      }

      // Make sure it doesn't exceed the maximum size of the array
      // Count if it does so we can see how often it happens
      if(nHitsR > MAX_TRACK_HITS){
        maxHitsLimit++;
        nHitsR = MAX_TRACK_HITS;
      }

      // Now access the variables of interest
      Float_t *dEdxArr = evt->trkdedx_pandoraTrack[iTrk][bestPlane];
      Float_t *RRArr   = evt->trkresrg_pandoraTrack[iTrk][bestPlane];
      Float_t *dQdxArr = evt->trkdqdx_pandoraTrack[iTrk][bestPlane];

      // Convert them to vectors
      std::vector<float> dEdx(dEdxArr, dEdxArr + nHitsR);
      std::vector<float> ResRg(RRArr, RRArr + nHitsR);
      std::vector<float> dQdx(dQdxArr, dQdxArr + nHitsR);

      /* Not sure this is necessary
      // Set the hit counter to 0 for this plane
      unsigned int iHit = 0;
      for(unsigned int itHit = 0; itHit < nHits; ++itHit){
        // Once we have reached the number of reconstructed hits, break out of the loop
        if(iHit >= nHitsR) break;
        if(evt->hit_plane[itHit] != bestPlane) continue;
        */

      // Now loop over hits so we can work our calo magic
      for(unsigned int iHit = 0; iHit < nHitsR; ++iHit){

        // General geometry of the track
        float x = evt->trkxyz_pandoraTrack[iTrk][bestPlane][iHit][0];
        float y = evt->trkxyz_pandoraTrack[iTrk][bestPlane][iHit][1];
        float z = evt->trkxyz_pandoraTrack[iTrk][bestPlane][iHit][2];
        float t = x * evtProc.kXtoT;

        h_depX->Fill(x);

        // Check if x is lower or higher than the APA bounds, charge seems to accumulate there
        // Put a 5 cm fiducial border around the APAs and CPAs too
        if(x < evtProc.PD_APA_X_POSITIONS[0]+5 || x > evtProc.PD_APA_X_POSITIONS[1]-5) continue;
        if(x > evtProc.PD_CPA_X_POSITIONS[0]-5 && x < evtProc.PD_CPA_X_POSITIONS[0]+5) continue;

        // Lifetime correction
        int tpc  = evtProc.WhichTPC(x,true) + 1;
        float simtau = evt->taulife/1.0e3; // convert from us to ms
        float dx = ( -1 + 2*(tpc%2) )*(x - evtProc.PD_APA_X_POSITIONS[tpc/2]);
        float dt = dx*evtProc.kXtoT;
        float corr  = TMath::Exp(-dt/etau.at(0));
        float eCorr = TMath::Exp(-dt/etau.at(0)) / TMath::Exp(-dt/simtau); // Correct for the already-corrected energy
        if(LifetimePerTPC){
          corr  = TMath::Exp(-dt/etau.at(tpc-1));
          eCorr = TMath::Exp(-dt/etau.at(tpc-1)) / TMath::Exp(-dt/simtau);
        }

        // New values
        float dEdxVal   = dEdx.at(iHit);
        float dQdxVal   = dQdx.at(iHit)*5.; // Translate to FD ADC definition: 1 ADC = 200 e vs PD: 1 ADC = 1000 e (from Tingjun, need to find in code)
        float RRVal     = ResRg.at(iHit)/100.;
        float dQdxCorr  = dQdxVal/corr;
        float dEdxCorr  = dEdxVal/eCorr;
        float dEdQVal   = dEdxVal/dQdxCorr;
        float dEdQCorr  = dEdxCorr/dQdxCorr;

        if(isTrueStopping && !thruGoing){
          h_dqdx_RR_stop->Fill(RRVal,dQdxVal);
          h_corr_dqdx_RR_stop->Fill(RRVal,dQdxCorr);
        }

        // the following studies should be conducted with top-bottom muons to start with
        if(thru != thruGoing) continue;

        h_dedx_x->Fill(x,dEdxVal);
        h_dqdx_x->Fill(x,dQdxVal);
        h_dqdx_t->Fill(t,dQdxVal);
        h_dqdx_E->Fill(eng,dQdxVal);
        h_dqdx_RR->Fill(RRVal,dQdxVal);
        h_corr_dedx_x->Fill(x,dEdxCorr);
        h_corr_dqdx_x->Fill(x,dQdxCorr);
        h_corr_dqdx_t->Fill(t,dQdxCorr);
        h_corr_dqdx_E->Fill(eng,dQdxCorr);
        h_corr_dqdx_RR->Fill(RRVal,dQdxCorr);
        h_corr_dedq_x->Fill(x,dEdQVal);
        h_corr2_dedq_x->Fill(x,dEdQCorr);

        h_hits_xy->Fill(x,y);
        h_hits_xz->Fill(x,z);
        h_hits_yz->Fill(y,z);
        h_hits_xyz->Fill(x,y,z);

        // Fill the maps to sort the median
        // First, XY
        int globalBinXY = h_dqdx_xy->FindBin(x,y);
        int globalBinXZ = h_dqdx_xz->FindBin(x,z);
        int globalBinYZ = h_dqdx_yz->FindBin(y,z);
        std::vector<int> globalBins{globalBinXY, globalBinXZ, globalBinYZ};

        for(unsigned int h = 0; h < dqdx_hists.size(); ++h){
          // Now fill the vectors
          int globalBin = globalBins.at(h);
          if(bin_dqdx_maps.at(h).find(globalBin) != bin_dqdx_maps.at(h).end()){
            bin_dqdx_maps.at(h).at(globalBin).push_back(dQdxCorr);
          }
          else
            outOfRange++;
        } // Histogram loop
        // Increment the hit counter for this plane
        //iHit++;
      } // Hits
    } // Tracks
  }// Event loop
  std::cout << " --- 100 % --- |" << std::endl;

  // Now calculate the median and set the entry in the map
  for(unsigned int h = 0; h < dqdx_hists.size(); ++h){
    TH2D *hist    = dqdx_hists.at(h);
    TH2D *errHist = error_hists.at(h);
    // Fill the histogram as normal
    for(int xBin = 0; xBin <= hist->GetNbinsX(); ++xBin){
      for(int yBin = 0; yBin <= hist->GetNbinsY(); ++yBin){
        int globalBin = hist->GetBin(xBin,yBin);
        const float *dqdxVals = bin_dqdx_maps.at(h).at(globalBin).data();
        median_dqdx_maps.at(h).emplace(globalBin,TMath::Median(bin_dqdx_maps.at(h).at(globalBin).size(),dqdxVals));
        if(median_dqdx_maps.at(h).at(globalBin) > 0)
          hist->SetBinContent(xBin,yBin,median_dqdx_maps.at(h).at(globalBin));

        double statErr = TMath::Sqrt(bin_dqdx_maps.at(h).at(globalBin).size())/static_cast<double>(bin_dqdx_maps.at(h).at(globalBin).size());
        errHist->SetBinContent(xBin,yBin,statErr);
      } // ybin
    } // Xbin
    // Now fill the median 1D histograms so we can evaluate the spread
    for(std::map<int,float>::const_iterator it = median_dqdx_maps.at(h).begin(); it != median_dqdx_maps.at(h).end(); ++it){
      median_hists.at(h)->Fill(it->second);
    } // 1D histogram loop
  } // Histogram loop

  std::cout << " Number of times max hits exceeds limit: " << maxHitsLimit << std::endl;
  std::cout << " Wrong way tracks:                       " << wrongWay << std::endl;

  // Sort out the TeX file
  std::vector<std::string> contents{
    "Events",
    "Tracks",
    "Muons",
    "$\\mu > 3$~m",
    "Crosses top or bottom",
    "Crosses top and bottom",
    "Crosses front or back",
    "Crosses front and back",
    "Crosses 1 APA or CPA",
    "Crosses $ \\geq $ 2 APA or CPA",
    "Stopping",
    "Exiting",
  };
  std::vector<unsigned int> rates{
    nEvts,
    totalTracks,
    nMu,
    nLongTracks,
    topOrBottom,
    topBottom,
    frontOrBack,
    frontBack,
    eq1APACPA,
    min2APACPA,
    stopping,
    exiting,
  };

  ofstream texFile;
  texFile.open(location+"/sample_contents_events"+tag+".tex");
  WriteStatsToTeX(texFile, n, contents, rates, static_cast<double>(nLongTracks), "$\\mu > 3$~m");
  
  TCanvas *c0 = new TCanvas("c0","",900,900);
  SetCanvasStyle(c0, 0.12,0.08,0.06,0.12,0,0,0);

  TLegend *l = new TLegend(0.47,0.94,0.98,0.995);
  l->SetNColumns(2);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->SetTextFont(132);

  // Draw the location of the beam pipe
  double lineMinX = std::min(h_startX->GetMinimum(),h_endX->GetMinimum());
  double lineMaxX = std::max(h_startX->GetMaximum(),h_endX->GetMaximum());

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
  h_startX->GetYaxis()->SetRangeUser(0,lineMaxX*1.05);

  l->AddEntry(h_startX, "Start", "l");
  l->AddEntry(h_endX, "End", "l");
  l->Draw();
  
  TLine *lX = new TLine(beamStart.X(), 0, beamStart.X(), lineMaxX*1.05);
  lX->SetLineColor(pal.at(0));
  lX->SetLineWidth(3);
  lX->SetLineStyle(7);
  lX->Draw();

  TLine *lXE = new TLine(beamEnd.X(), 0, beamEnd.X(), lineMaxX*1.05);
  lXE->SetLineColor(pal.at(1));
  lXE->SetLineWidth(3);
  lXE->SetLineStyle(7);
  lXE->Draw();

  c0->SaveAs((location+"/reco_tracks_start_endX"+tag+".png").c_str());
  c0->SaveAs((location+"/reco_tracks_start_endX"+tag+".root").c_str());
  c0->Clear();
  l->Clear();

  // Draw the location of the beam pipe
  double lineMinY = std::min(h_startY->GetMinimum(),h_endY->GetMinimum());
  double lineMaxY = std::max(h_startY->GetMaximum(),h_endY->GetMaximum());

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
  h_startY->GetYaxis()->SetRangeUser(0,lineMaxY*1.05);

  l->AddEntry(h_startY, "Start", "l");
  l->AddEntry(h_endY, "End", "l");
  l->Draw();
  
  TLine *lY = new TLine(beamStart.Y(), 0, beamStart.Y(), lineMaxY*1.05);
  lY->SetLineColor(pal.at(0));
  lY->SetLineWidth(3);
  lY->SetLineStyle(7);
  lY->Draw();

  TLine *lYE = new TLine(beamEnd.Y(), 0, beamEnd.Y(), lineMaxY*1.05);
  lYE->SetLineColor(pal.at(1));
  lYE->SetLineWidth(3);
  lYE->SetLineStyle(7);
  lYE->Draw();

  c0->SaveAs((location+"/reco_tracks_start_endY"+tag+".png").c_str());
  c0->SaveAs((location+"/reco_tracks_start_endY"+tag+".root").c_str());
  c0->Clear();
  l->Clear();

  // Draw the location of the beam pipe
  double lineMinZ = std::min(h_startZ->GetMinimum(),h_endZ->GetMinimum());
  double lineMaxZ = std::max(h_startZ->GetMaximum(),h_endZ->GetMaximum());

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
  h_startZ->GetYaxis()->SetRangeUser(0,lineMaxZ*1.05);

  l->AddEntry(h_startZ, "Start", "l");
  l->AddEntry(h_endZ, "End", "l");
  l->Draw();
  
  c0->SaveAs((location+"/reco_tracks_start_endZ"+tag+".png").c_str());
  c0->SaveAs((location+"/reco_tracks_start_endZ"+tag+".root").c_str());
  c0->Clear();
  l->Clear();

  SetHistogramStyle1D(h_depX,"Deposition x position [cm]", "Rate");
  h_depX->Draw("hist");
  h_depX->SetLineWidth(3);
  h_depX->SetLineColor(pal.at(0));
  c0->SaveAs((location+"/depX"+tag+".png").c_str());
  c0->SaveAs((location+"/depX"+tag+".root").c_str());
  c0->Clear();

  TCanvas *c1 = new TCanvas("c1","",1000,800);
  SetCanvasStyle(c1, 0.1,0.12,0.05,0.12,0,0,0);

  std::vector<TLine*> APACPALines;
  // Now draw lines and labels where the APA and CPAs are
  for(unsigned int iA = 0; iA < 2; ++iA){
    TLine *l = new TLine(evtProc.PD_APA_X_POSITIONS[iA], 0, evtProc.PD_APA_X_POSITIONS[iA], 10);
    l->SetLineColor(kWhite);
    l->SetLineWidth(3);
    l->SetLineStyle(2);
    APACPALines.push_back(l);
  }
  for(unsigned int iC = 0; iC < 1; ++iC){
    TLine *l = new TLine(evtProc.PD_CPA_X_POSITIONS[iC], 0, evtProc.PD_CPA_X_POSITIONS[iC], 10);
    l->SetLineColor(kWhite);
    l->SetLineWidth(3);
    l->SetLineStyle(2);
    APACPALines.push_back(l);
  }

  // dEdx vs x
  SetHistogramStyle2D(h_dedx_x,"Reconstructed drift coordinate [cm]", " dE/dx [MeV/cm]",false);
  h_dedx_x->GetZaxis()->SetLabelSize(0.03);
  h_dedx_x->GetZaxis()->SetLabelFont(132);
  h_dedx_x->Draw("colz");

  // Draw the APA and CPA lines and labels
  for(unsigned int iLine = 0; iLine < APACPALines.size(); ++iLine){
    APACPALines.at(iLine)->SetY2(10);
    APACPALines.at(iLine)->Draw();
  }

  FormatLatex(evtProc.PD_APA_X_POSITIONS[0]+10,5.2, "#color[0]{APA}");
  FormatLatex(evtProc.PD_CPA_X_POSITIONS[0]+10,5.2, "#color[0]{CPA}");
  FormatLatex(evtProc.PD_APA_X_POSITIONS[1]+10,5.2, "#color[0]{APA}");

  c1->SaveAs((location+"/dEdx_vs_X"+tag+".png").c_str());
  c1->SaveAs((location+"/dEdx_vs_X"+tag+".root").c_str());
  c1->Clear();

  // corrected dEdx vs x
  SetHistogramStyle2D(h_corr_dedx_x,"Reconstructed drift coordinate [cm]", " dE/dx [MeV/cm]",false);
  h_corr_dedx_x->GetZaxis()->SetLabelSize(0.03);
  h_corr_dedx_x->GetZaxis()->SetLabelFont(132);
  h_corr_dedx_x->Draw("colz");

  // Draw the APA and CPA lines and labels
  for(unsigned int iLine = 0; iLine < APACPALines.size(); ++iLine){
    APACPALines.at(iLine)->Draw();
  }

  FormatLatex(evtProc.PD_APA_X_POSITIONS[0]+10,5.2, "#color[0]{APA}");
  FormatLatex(evtProc.PD_CPA_X_POSITIONS[0]+10,5.2, "#color[0]{CPA}");
  FormatLatex(evtProc.PD_APA_X_POSITIONS[1]+10,5.2, "#color[0]{APA}");

  c1->SaveAs((location+"/corr_dEdx_vs_X"+tag+".png").c_str());
  c1->SaveAs((location+"/corr_dEdx_vs_X"+tag+".root").c_str());
  c1->Clear();

  // Charge
  SetHistogramStyle2D(h_dqdx_x,"Reconstructed drift coordinate [cm]", " dQ/dx [ADC/cm]",false);
  h_dqdx_x->GetZaxis()->SetLabelSize(0.03);
  h_dqdx_x->GetZaxis()->SetLabelFont(132);
  h_dqdx_x->Draw("colz");

  // Draw the APA and CPA lines and labels
  for(unsigned int iLine = 0; iLine < APACPALines.size(); ++iLine){
    APACPALines.at(iLine)->SetY2(200);
    APACPALines.at(iLine)->Draw();
  }

  FormatLatex(evtProc.PD_APA_X_POSITIONS[0]+10,500, "#color[0]{APA}");
  FormatLatex(evtProc.PD_CPA_X_POSITIONS[0]+10,500, "#color[0]{CPA}");
  FormatLatex(evtProc.PD_APA_X_POSITIONS[1]+10,500, "#color[0]{APA}");

  c1->SaveAs((location+"/charge_vs_X"+tag+".png").c_str());
  c1->SaveAs((location+"/charge_vs_X"+tag+".root").c_str());
  c1->Clear();

  SetHistogramStyle2D(h_corr_dqdx_x,"Reconstructed drift coordinate [cm]", " dQ/dx [ADC/cm]", false);
  h_corr_dqdx_x->GetZaxis()->SetLabelSize(0.03);
  h_corr_dqdx_x->GetZaxis()->SetLabelFont(132);
  h_corr_dqdx_x->Draw("colz");

  // Draw the APA and CPA lines and labels
  for(unsigned int iLine = 0; iLine < APACPALines.size(); ++iLine){
    APACPALines.at(iLine)->SetY2(200);
    APACPALines.at(iLine)->Draw();
  }

  FormatLatex(evtProc.PD_APA_X_POSITIONS[0]+10,500, "#color[0]{APA}");
  FormatLatex(evtProc.PD_CPA_X_POSITIONS[0]+10,500, "#color[0]{CPA}");
  FormatLatex(evtProc.PD_APA_X_POSITIONS[1]+10,500, "#color[0]{APA}");

  c1->SaveAs((location+"/corr_charge_vs_X"+tag+".png").c_str());
  c1->SaveAs((location+"/corr_charge_vs_X"+tag+".root").c_str());
  c1->Clear();

  // Charge vs time
  SetHistogramStyle2D(h_dqdx_t,"Reconstructed drift time [ms]", " dQ/dx [ADC/cm]",false);
  h_dqdx_t->GetZaxis()->SetLabelSize(0.03);
  h_dqdx_t->GetZaxis()->SetLabelFont(132);
  h_dqdx_t->Draw("colz");

  // Draw the APA and CPA lines and labels
  for(unsigned int iLine = 0; iLine < APACPALines.size(); ++iLine){
    APACPALines.at(iLine)->SetY2(200);
    APACPALines.at(iLine)->Draw();
  }

  FormatLatex((evtProc.PD_APA_X_POSITIONS[0]+10)*evtProc.kXtoT,500, "#color[0]{APA}");
  FormatLatex((evtProc.PD_CPA_X_POSITIONS[0]+10)*evtProc.kXtoT,500, "#color[0]{CPA}");
  FormatLatex((evtProc.PD_APA_X_POSITIONS[1]+10)*evtProc.kXtoT,500, "#color[0]{APA}");

  c1->SaveAs((location+"/charge_vs_T"+tag+".png").c_str());
  c1->SaveAs((location+"/charge_vs_T"+tag+".root").c_str());
  c1->Clear();

  SetHistogramStyle2D(h_corr_dqdx_t,"Reconstructed drift time [ms]", " dQ/dx [ADC/cm]", false);
  h_corr_dqdx_t->GetZaxis()->SetLabelSize(0.03);
  h_corr_dqdx_t->GetZaxis()->SetLabelFont(132);
  h_corr_dqdx_t->Draw("colz");

  // Draw the APA and CPA lines and labels
  for(unsigned int iLine = 0; iLine < APACPALines.size(); ++iLine){
    APACPALines.at(iLine)->SetY2(200);
    APACPALines.at(iLine)->Draw();
  }

  FormatLatex((evtProc.PD_APA_X_POSITIONS[0]+10)*evtProc.kXtoT,500, "#color[0]{APA}");
  FormatLatex((evtProc.PD_CPA_X_POSITIONS[0]+10)*evtProc.kXtoT,500, "#color[0]{CPA}");
  FormatLatex((evtProc.PD_APA_X_POSITIONS[1]+10)*evtProc.kXtoT,500, "#color[0]{APA}");

  c1->SaveAs((location+"/corr_charge_vs_T"+tag+".png").c_str());
  c1->SaveAs((location+"/corr_charge_vs_T"+tag+".root").c_str());
  c1->Clear();

  SetHistogramStyle2D(h_corr_dedq_x,"Reconstructed drift coordinate [cm]", " Energy per charge deposition [MeV/ADC]", false);
  h_corr_dedq_x->GetZaxis()->SetLabelSize(0.03);
  h_corr_dedq_x->GetZaxis()->SetLabelFont(132);
  h_corr_dedq_x->Draw("colz");

  // Draw the APA and CPA lines and labels
  for(unsigned int iLine = 0; iLine < APACPALines.size(); ++iLine){
    APACPALines.at(iLine)->SetY1(6.6e-3);
    APACPALines.at(iLine)->SetY2(7.2e-3);
    APACPALines.at(iLine)->Draw();
  }

  FormatLatex(evtProc.PD_APA_X_POSITIONS[0]+10,7.1e-3, "#color[0]{APA}");
  FormatLatex(evtProc.PD_CPA_X_POSITIONS[0]+10,7.1e-3, "#color[0]{CPA}");
  FormatLatex(evtProc.PD_APA_X_POSITIONS[1]+10,7.1e-3, "#color[0]{APA}");

  c1->SaveAs((location+"/corr_energy_charge_vs_X"+tag+".png").c_str());
  c1->SaveAs((location+"/corr_energy_charge_vs_X"+tag+".root").c_str());
  c1->Clear();

  SetHistogramStyle2D(h_corr2_dedq_x,"Reconstructed drift coordinate [cm]", " Energy per charge deposition [MeV/ADC]", false);
  h_corr2_dedq_x->GetZaxis()->SetLabelSize(0.03);
  h_corr2_dedq_x->GetZaxis()->SetLabelFont(132);
  h_corr2_dedq_x->Draw("colz");

  // Draw the APA and CPA lines and labels
  for(unsigned int iLine = 0; iLine < APACPALines.size(); ++iLine){
    APACPALines.at(iLine)->SetY1(6.6e-3);
    APACPALines.at(iLine)->SetY2(8e-3);
    APACPALines.at(iLine)->Draw();
  }

  FormatLatex(evtProc.PD_APA_X_POSITIONS[0]+10,7.9e-3, "#color[0]{APA}");
  FormatLatex(evtProc.PD_CPA_X_POSITIONS[0]+10,7.9e-3, "#color[0]{CPA}");
  FormatLatex(evtProc.PD_APA_X_POSITIONS[1]+10,7.9e-3, "#color[0]{APA}");

  c1->SaveAs((location+"/corr2_energy_charge_vs_X"+tag+".png").c_str());
  c1->SaveAs((location+"/corr2_energy_charge_vs_X"+tag+".root").c_str());
  c1->Clear();

  SetHistogramStyle2D(h_dqdx_RR_stop,"Reconstructed residual range [m]", " dQ/dx [ADC/cm]", false);
  h_dqdx_RR_stop->GetZaxis()->SetLabelSize(0.03);
  h_dqdx_RR_stop->GetZaxis()->SetLabelFont(132);
  h_dqdx_RR_stop->Draw("colz");

  c1->SaveAs((location+"/charge_vs_RR_stop"+tag+".png").c_str());
  c1->SaveAs((location+"/charge_vs_RR_stop"+tag+".root").c_str());
  c1->Clear();
  
  SetHistogramStyle2D(h_corr_dqdx_RR_stop,"Reconstructed residual range [m]", " dQ/dx [ADC/cm]", false);
  h_corr_dqdx_RR_stop->GetZaxis()->SetLabelSize(0.03);
  h_corr_dqdx_RR_stop->GetZaxis()->SetLabelFont(132);
  h_corr_dqdx_RR_stop->Draw("colz");

  c1->SaveAs((location+"/corr_charge_vs_RR_stop"+tag+".png").c_str());
  c1->SaveAs((location+"/corr_charge_vs_RR_stop"+tag+".root").c_str());
  c1->Clear();

  SetHistogramStyle2D(h_dqdx_RR,"Reconstructed residual range [m]", " dQ/dx [ADC/cm]", false);
  h_dqdx_RR->GetZaxis()->SetLabelSize(0.03);
  h_dqdx_RR->GetZaxis()->SetLabelFont(132);
  h_dqdx_RR->Draw("colz");

  c1->SaveAs((location+"/charge_vs_RR"+tag+".png").c_str());
  c1->SaveAs((location+"/charge_vs_RR"+tag+".root").c_str());
  c1->Clear();

  SetHistogramStyle2D(h_corr_dqdx_RR,"Reconstructed residual range [m]", " dQ/dx [ADC/cm]", false);
  h_corr_dqdx_RR->GetZaxis()->SetLabelSize(0.03);
  h_corr_dqdx_RR->GetZaxis()->SetLabelFont(132);
  h_corr_dqdx_RR->Draw("colz");

  c1->SaveAs((location+"/corr_charge_vs_RR"+tag+".png").c_str());
  c1->SaveAs((location+"/corr_charge_vs_RR"+tag+".root").c_str());
  c1->Clear();

  c1->SetLogx();
  SetHistogramStyle2D(h_dqdx_E,"True (generated) #mu energy [GeV]", " dQ/dx [ADC/cm]", false);
  h_dqdx_E->GetZaxis()->SetLabelSize(0.03);
  h_dqdx_E->GetZaxis()->SetLabelFont(132);
  h_dqdx_E->Draw("colz");

  c1->SaveAs((location+"/charge_vs_E"+tag+".png").c_str());
  c1->SaveAs((location+"/charge_vs_E"+tag+".root").c_str());
  c1->Clear();

  SetHistogramStyle2D(h_corr_dqdx_E,"True (generated) #mu energy [GeV]", " dQ/dx [ADC/cm]", false);
  h_corr_dqdx_E->GetZaxis()->SetLabelSize(0.03);
  h_corr_dqdx_E->GetZaxis()->SetLabelFont(132);
  h_corr_dqdx_E->Draw("colz");

  c1->SaveAs((location+"/corr_charge_vs_E"+tag+".png").c_str());
  c1->SaveAs((location+"/corr_charge_vs_E"+tag+".root").c_str());
  c1->Clear();

  // Number of hits XY
  TCanvas *c2 = new TCanvas("c2","",1000,800);
  SetCanvasStyle(c2, 0.1,0.12,0.05,0.12,0,0,0);
  SetHistogramStyle2D(h_hits_xy,"x [cm]", " y [cm]", false);
  h_hits_xy->GetZaxis()->SetLabelSize(0.03);
  h_hits_xy->GetZaxis()->SetLabelFont(132);
  h_hits_xy->Draw("colz");
  c2->SaveAs((location+"/xy_hits"+tag+".png").c_str());
  c2->SaveAs((location+"/xy_hits"+tag+".root").c_str());
  c2->Clear();

  // Number of hits XY
  SetHistogramStyle2D(h_hits_xz,"x [cm]"," z [cm]", false);
  h_hits_xz->GetZaxis()->SetLabelSize(0.03);
  h_hits_xz->GetZaxis()->SetLabelFont(132);
  h_hits_xz->Draw("colz");
  c2->SaveAs((location+"/xz_hits"+tag+".png").c_str());
  c2->SaveAs((location+"/xz_hits"+tag+".root").c_str());
  c2->Clear();

  // Number of hits YZ
  SetHistogramStyle2D(h_hits_yz,"y [cm]"," z [cm]", false);
  h_hits_yz->GetZaxis()->SetLabelSize(0.03);
  h_hits_yz->GetZaxis()->SetLabelFont(132);
  h_hits_yz->Draw("colz");
  c2->SaveAs((location+"/yz_hits"+tag+".png").c_str());
  c2->SaveAs((location+"/yz_hits"+tag+".root").c_str());
  c2->Clear();

  // Charge depositions in space
  // dQ/dx in XY space
  SetHistogramStyle2D(h_dqdx_xy,"x [cm]", " y [cm]", false);
  h_dqdx_xy->GetZaxis()->SetLabelSize(0.03);
  h_dqdx_xy->GetZaxis()->SetLabelFont(132);
  h_dqdx_xy->GetZaxis()->SetTitle("Reconstructed dQ/dx [ADC/cm]");
  h_dqdx_xy->GetZaxis()->SetTitleSize(0.03);
  h_dqdx_xy->GetZaxis()->SetTitleFont(132);
  h_dqdx_xy->Draw("colz");
  c2->SaveAs((location+"/xy_dqdx"+tag+".png").c_str());
  c2->SaveAs((location+"/xy_dqdx"+tag+".root").c_str());
  c2->Clear();

  // dQ/dx in XZ space
  SetHistogramStyle2D(h_dqdx_xz,"x [cm]"," z [cm]", false);
  h_dqdx_xz->GetZaxis()->SetLabelSize(0.03);
  h_dqdx_xz->GetZaxis()->SetLabelFont(132);
  h_dqdx_xz->GetZaxis()->SetTitle("Reconstructed dQ/dx [ADC/cm]");
  h_dqdx_xz->GetZaxis()->SetTitleSize(0.03);
  h_dqdx_xz->GetZaxis()->SetTitleFont(132);
  h_dqdx_xz->Draw("colz");
  c2->SaveAs((location+"/xz_dqdx"+tag+".png").c_str());
  c2->SaveAs((location+"/xz_dqdx"+tag+".root").c_str());
  c2->Clear();

  // dQ/dx in YZ space
  SetHistogramStyle2D(h_dqdx_yz,"y [cm]"," z [cm]", false);
  h_dqdx_yz->GetZaxis()->SetLabelSize(0.03);
  h_dqdx_yz->GetZaxis()->SetLabelFont(132);
  h_dqdx_yz->GetZaxis()->SetTitle("Reconstructed dQ/dx [ADC/cm]");
  h_dqdx_yz->GetZaxis()->SetTitleSize(0.03);
  h_dqdx_yz->GetZaxis()->SetTitleFont(132);
  h_dqdx_yz->Draw("colz");
  c2->SaveAs((location+"/yz_dqdx"+tag+".png").c_str());
  c2->SaveAs((location+"/yz_dqdx"+tag+".root").c_str());
  c2->Clear();

  // dQ/dx in XY space
  SetHistogramStyle2D(h_statErr_xy,"x [cm]", " y [cm]", false);
  h_statErr_xy->GetZaxis()->SetLabelSize(0.03);
  h_statErr_xy->GetZaxis()->SetLabelFont(132);
  h_statErr_xy->GetZaxis()->SetTitle("Fractional statistical uncertainty");
  h_statErr_xy->GetZaxis()->SetTitleSize(0.03);
  h_statErr_xy->GetZaxis()->SetTitleFont(132);
  h_statErr_xy->Draw("colz");
  c2->SaveAs((location+"/xy_statErr"+tag+".png").c_str());
  c2->SaveAs((location+"/xy_statErr"+tag+".root").c_str());
  c2->Clear();

  // dQ/dx in XZ space
  SetHistogramStyle2D(h_statErr_xz,"x [cm]"," z [cm]", false);
  h_statErr_xz->GetZaxis()->SetLabelSize(0.03);
  h_statErr_xz->GetZaxis()->SetLabelFont(132);
  h_statErr_xz->GetZaxis()->SetTitle("Fractional statistical uncertainty");
  h_statErr_xz->GetZaxis()->SetTitleSize(0.03);
  h_statErr_xz->GetZaxis()->SetTitleFont(132);
  h_statErr_xz->Draw("colz");
  c2->SaveAs((location+"/xz_statErr"+tag+".png").c_str());
  c2->SaveAs((location+"/xz_statErr"+tag+".root").c_str());
  c2->Clear();

  // dQ/dx in YZ space
  SetHistogramStyle2D(h_statErr_yz,"y [cm]"," z [cm]", false);
  h_statErr_yz->GetZaxis()->SetLabelSize(0.03);
  h_statErr_yz->GetZaxis()->SetLabelFont(132);
  h_statErr_yz->GetZaxis()->SetTitle("Fractional statistical uncertainty");
  h_statErr_yz->GetZaxis()->SetTitleSize(0.03);
  h_statErr_yz->GetZaxis()->SetTitleFont(132);
  h_statErr_yz->Draw("colz");
  c2->SaveAs((location+"/yz_statErr"+tag+".png").c_str());
  c2->SaveAs((location+"/yz_statErr"+tag+".root").c_str());
  c2->Clear();

  TCanvas *c3 = new TCanvas("c3","",1000,800);
  SetCanvasStyle(c3, 0.1,0.12,0.05,0.12,0,0,0);
  SetHistogramStyle3D(h_hits_xyz,"x [cm]","y [cm]","z [cm]");
  h_hits_xyz->SetMarkerStyle(33);
  h_hits_xyz->SetMarkerColor(pal.at(1));
  h_hits_xyz->Draw();
  c3->SaveAs((location+"/xyz_hits"+tag+".png").c_str());
  c3->SaveAs((location+"/xyz_hits"+tag+".root").c_str());
  
  // Plane crossing
  TCanvas *c4 = new TCanvas("c4","",900,900);
  SetCanvasStyle(c4, 0.12,0.05,0.06,0.15,0,0,0);

  l->SetNColumns(3);
  l->SetX1NDC(0.22);

  SetHistogramStyle1D(h_plane_cross,"Plane label", " Number of tracks crossing plane");
  h_plane_cross->Draw("hist");
  h_plane_enter->Draw("same");
  h_plane_exit->Draw("same");
  h_plane_cross->SetLineWidth(2);
  h_plane_enter->SetLineWidth(2);
  h_plane_exit->SetLineWidth(2);
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

  c4->SaveAs((location+"/planes_crossed_entered_exited"+tag+".png").c_str());
  c4->SaveAs((location+"/planes_crossed_entered_exited"+tag+".root").c_str());
  c4->Clear();
  l->Clear();
  
  TCanvas *c5 = new TCanvas("c5","",900,900);
  SetCanvasStyle(c5, 0.12,0.05,0.05,0.12,0,0,0);

  SetHistogramStyle1D(h_n_crossed,"Number of planes crossed [P]", " Number of tracks crossing P planes");
  h_n_crossed->Draw("hist");
  h_n_crossed->SetLineWidth(2);
  h_n_crossed->SetLineColor(pal.at(0));
  h_n_crossed->GetYaxis()->SetTitleOffset(0.95);
  c5->SaveAs((location+"/tracks_crossed_nplanes"+tag+".png").c_str());
  c5->SaveAs((location+"/tracks_crossed_nplanes"+tag+".root").c_str());
  c5->Clear();
  
  TCanvas *c6 = new TCanvas("c6","",900,900);
  SetCanvasStyle(c6, 0.12,0.05,0.06,0.15,0,0,0);

  SetHistogramStyle1D(h_enter_dist,"Distance from candidate entrance/exit [cm]", " Rate");
  h_enter_dist->Draw("hist");
  h_exit_dist->Draw("same");
  h_enter_dist->SetLineWidth(2);
  h_exit_dist->SetLineWidth(2);
  h_enter_dist->SetLineStyle(2);
  h_exit_dist->SetLineStyle(3);
  h_enter_dist->SetLineColor(pal.at(1));
  h_exit_dist->SetLineColor(pal.at(0));
  h_enter_dist->GetYaxis()->SetTitleOffset(0.95);

  l->SetNColumns(2);
  l->SetX1NDC(0.47);
  l->AddEntry(h_enter_dist, "Entrance", "l");
  l->AddEntry(h_exit_dist, "Exit", "l");
  l->Draw();

  c6->SaveAs((location+"/distance_to_entrance_exit_planes"+tag+".png").c_str());
  c6->SaveAs((location+"/distance_to_entrance_exit_planes"+tag+".root").c_str());
  c6->Clear();
  l->Clear();
  
  SetHistogramStyle1D(h_muon_length,"Muon length [m]", " Rate");
  h_muon_length->Draw("hist");
  h_muon_length->SetLineWidth(2);
  h_muon_length->SetLineColor(pal.at(1));
  h_muon_length->GetYaxis()->SetTitleOffset(0.95);
  c6->SaveAs((location+"/muon_length"+tag+".png").c_str());
  c6->SaveAs((location+"/muon_length"+tag+".root").c_str());
  c6->Clear();
 
  ofstream oFile;
  oFile.open((location+"/statistics"+tag+".txt").c_str());

  oFile << " Fits to the variation of the median dQ/dx values in 20cm bins across the detector" << std::endl;

  // Fit to each of the median distributions
  double maxy = -99999.; 
  for(unsigned int h = 0; h < median_hists.size(); ++h){
    TH1D *hist = median_hists.at(h);
    oFile << " " << hist->GetTitle() << std::endl;

    SetHistogramStyle1D(hist,"Median dQ/dx [ADC/cm]", " Rate");
    l->AddEntry(hist,hist->GetTitle(),"L");
    hist->Scale(1/static_cast<double>(hist->Integral()));
    if(hist->GetMaximum() > maxy)
      maxy = hist->GetMaximum();

    hist->SetTitle("");

    // Fit to the distributions to extract the width
    // First, get the range to fit
    int maxbin  = hist->GetMaximumBin();
    double minX = hist->GetXaxis()->GetXmin();
    double maxX = hist->GetXaxis()->GetXmax();
    if(maxbin-nBinsFromPeak > 1)
      minX = hist->GetBinCenter(maxbin-nBinsFromPeakL);
    if(maxbin+nBinsFromPeak < hist->GetNbinsX())
      maxX = hist->GetBinCenter(maxbin+nBinsFromPeakR);
  
    if(h == 0){
      hist->Draw("hist");
    }
    else{
      hist->Draw("hist same");
    }
    
    // Now define the fit function and do the fit
    TF1 *fit    = new TF1(("fit_"+std::to_string(h)).c_str(),"gaus",minX,maxX);
    auto result = hist->Fit(fit, "LEQSMR", "");
    //fit->Draw("hist same");

    // Now write the results
    oFile << " Chi2: " << result->Chi2() << " / " << result->Ndf() << std::endl;
    for(unsigned int p = 0; p < result->NPar(); ++p){
      oFile << " " << result->ParName(p) << " : " << result->Parameter(p) << " +/- " << result->Error(p) << " [ADC/cm] "  << std::endl;
    }

    //FormatStats(hist,1110,101);
    hist->SetStats(kFALSE);
    hist->GetYaxis()->SetTitleOffset(0.95);
    hist->SetLineWidth(2);
    hist->SetLineColor(pal.at(h));
    fit->SetLineStyle(7);
    fit->SetLineWidth(3);
    fit->SetLineColor(pal.at(h));
  }
  median_hists.at(0)->GetYaxis()->SetRangeUser(0,maxy*1.1);
  l->SetNColumns(3);
  l->SetX1NDC(0.22);
  l->Draw();
  c6->SaveAs((location+"/median_dqdx"+tag+".png").c_str());
  c6->SaveAs((location+"/median_dqdx"+tag+".root").c_str());
  c6->Clear();
  l->Clear();
  
  oFile.close();
  
  // Now try making an event display, if we like
  if(drawPFParticles){
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

    cEvd->SaveAs((location+"/pfparticle_display"+tag+".png").c_str());
    cEvd->SaveAs((location+"/pfparticle_display"+tag+".root").c_str());
  
  } // Draw PFParticles

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
