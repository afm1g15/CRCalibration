/************************************************************************
 *
 * A Macro to plot various distrubtions in order to 
 * understand the contents of the CR sample for 
 * through-going muons
 *
 * Example file list located here:
 *   /home/jones/work/cosmics/LArSoft-v08_50_00/work/files/anafiles.list
 *
 * Parameters to look into:
 *   - Muon momentum
 *   - ThetaXZ
 *   - ThetaYZ
 *   - Y vs Z position
 *   - X vs Y position
 *   - Entry plane
 *   - Exit plane
 *   - Hit distribution in YZ
 *   - Hist distribution in XY
 *
 * Ultimately, compare with Viktor and Praveen's analogous studies
 *
 *************************************************************************/

#include "EventProcessor.h"
#include "ConfigReader.h"

using namespace calib;
using namespace cppsecrets;

// Allowed branches to read from the tree
std::vector<TString> allowed = {
   "run",
   "event",
   "genie_primaries_pdg",
   "inTPCActive",
   "pdg",
   "Mother",
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
   "trkidtruth_pandoraTrack",
   "ntracks_pandoraTrack",
   "trkId_pandoraTrack",
   "ntrkhits_pandoraTrack",
   "trkdqdx_pandoraTrack",
   "trkdedx_pandoraTrack",
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
  {"h3", "CPA 2"},
  {"h4", "APA 3"},
  {"t",  "Top"},
  {"bo", "Bot."},
  {"f",  "Fro."},
  {"ba", "Back"},
};

typedef std::vector<Plane> PlaneList;
     
int fileContentStudies(const char *config){

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

  // Get configuration variables
  int n = -1;
  int yCut = 1;
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
  p->getValue("YCut",      yCut);
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
  n = evtProc.GetFiles();
  std::cout << " Number of files: " << n << std::endl;
  
  // Start of analysis (loop over chain and events
  std::cout << " Running analysis..." << std::endl;

  // Then setup the histograms, counters and any other variables to add to
  // Setup histograms
  TH2D *h_dedx_x       = new TH2D("h_dedx_x","",100,-800,800,100,0,10);
  TH2D *h_dqdx_x       = new TH2D("h_dqdx_x","",100,-800,800,100,0,1000);
  TH2D *h_corr_dqdx_x  = new TH2D("h_corr_dqdx_x","",100,-800,800,100,0,1000);
  TH2D *h_corr_dedq_x  = new TH2D("h_corr_dedq_x","",100,-800,800,100,6.6e-3,7.2e-3);
  TH2D *h_corr2_dedq_x = new TH2D("h_corr2_dedq_x","",100,-800,800,100,6.5e-3,8e-3);
  TH2D *h_hits_xy      = new TH2D("h_hits_xy","",100,-800,800,100,-650,650);
  TH2D *h_hits_xz      = new TH2D("h_hits_xz","",100,-800,800,300,-200,6000);
  TH2D *h_hits_yz      = new TH2D("h_hits_yz","",100,-700,700,300,-200,6000);
  TH3D *h_hits_xyz     = new TH3D("h_hits_xyz","",100,-800,800,100,-700,700,300,-200,6000);
  TH1D *h_plane_cross  = new TH1D("h_plane_cross","",9,0,9); // Number of tracks crossing each plane
  TH1D *h_plane_enter  = new TH1D("h_plane_enter","",9,0,9); // Number of tracks entering from each external plane
  TH1D *h_plane_exit   = new TH1D("h_plane_exit","",9,0,9); // Number of tracks exiting from each external plane
  TH1D *h_enter_dist   = new TH1D("h_enter_dist","",200,0,10); // Number of tracks entering from each external plane
  TH1D *h_exit_dist    = new TH1D("h_exit_dist","",200,0,10); // Number of tracks entering from each external plane
  TH1D *h_muon_length  = new TH1D("h_muon_length","",200,0,2000); // Muon length
  TH1D *h_n_crossed    = new TH1D("h_n_crossed","",9,0,9); // Number of planes crossed by each track
  
  // Setup counters
  unsigned int maxHitsLimit     = 0;
  unsigned int wrongWay         = 0;
  unsigned int noPlanes         = 0;
  unsigned int totalTracks      = 0;
  unsigned int nMu              = 0;
  unsigned int nPrimaryMu       = 0;
  unsigned int nLongTracks      = 0;
  unsigned int nLongHighYTracks = 0;
  unsigned int topBottom        = 0;
  unsigned int topOrBottom      = 0;
  unsigned int min2APACPA       = 0;
  unsigned int min1APACPA       = 0;
  unsigned int stopping         = 0;

  // Now loop over the events
  unsigned int nEvts = tree->GetEntries();
  unsigned int iIt = 1;

  std::cout << " |";
  for(unsigned int iEvt = 0; iEvt < nEvts; ++iEvt){
    tree->GetEntry(iEvt);
    if(!evtProc.SelectEvent(evt)) continue;
    unsigned int nTrks = evt->ntracks_pandoraTrack;
    
    // Print the processing rate
    double evtFrac  = iEvt/static_cast<double>(nEvts);

    if(std::abs(0.1*iIt-evtFrac) < std::numeric_limits<double>::epsilon()){
      std::cout << " --- " << evtFrac*100 << " %";
      std::cout.flush();
      iIt++;
    }

    // Now loop over the tracks so we can do some stuff!!
    for(unsigned int iTrk = 0; iTrk < nTrks; ++iTrk){

      // Setup list of plane labels the track has crossed
      std::vector<std::string> labelsCrossed;

      // Get the best plane
      int bestPlane = 0;
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
      if(abs(evt->trkpdgtruth_pandoraTrack[iTrk][bestPlane]) != 13) continue;
      nMu++;
      
      // Length cuts (2m)
      if(!evtProc.SelectTrack(evt,iTrk)) continue;
      nLongTracks++;
      
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
        wrongWay++;
        TVector3 temp(endVtx);
        endVtx = startVtx;
        startVtx = temp;
      }

      float length = evt->trklen_pandoraTrack[iTrk];
      h_muon_length->Fill(length);

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
      
      for(const Plane &pl : fidExtPlanes){
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
      
      // Now fill the number of planes crossed histogram
      h_n_crossed->Fill(nPlanesCrossed);
      if(nPlanesCrossed == 0){
        noPlanes++;
      }
      if(nExtCrossed == 1){
        stopping++;
      }

      // Loop over the labels crossed by this track and fill the appropriate counters
      bool APA = false;
      bool CPA = false;
      bool top = false;
      bool bot = false;

      bool thruGoing = false;
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
      if(top && bot){
        thruGoing = true;
        topBottom++;
      }

      // the following studies should be conducted with top-bottom muons to start with
      if(thru == 1 && !thruGoing) continue;

      if(yCut == 1 && startVtx.Y() < 599.5) continue; // Make sure the tracks start at the top of the detector
      nLongHighYTracks++;

      // Now fill dQ/dx and dE/dx and hit histograms for each of the three wire planes
      // Somehow flag the best wire plane histogram
      for(int iPlane = 0; iPlane < 3; ++iPlane){

        // Use only best plane for now
        if(iPlane != bestPlane) continue;

        unsigned int nHits = evt->ntrkhits_pandoraTrack[iTrk][iPlane];
     
        // Make sure it doesn't exceed the maximum size of the array
        // Count if it does so we can see how often it happens
        if(nHits > MAX_TRACK_HITS){
          maxHitsLimit++;
          nHits = MAX_TRACK_HITS;
        }

        // Now access the variables of interest
        Float_t *dEdxArr = evt->trkdedx_pandoraTrack[iTrk][iPlane];
        Float_t *dQdxArr = evt->trkdqdx_pandoraTrack[iTrk][iPlane];

        // Convert them to vectors
        std::vector<float> dEdx(dEdxArr, dEdxArr + nHits);
        std::vector<float> dQdx(dQdxArr, dQdxArr + nHits);

        // Now loop over hits so we can work our calo magic
        for(unsigned int iHit = 0; iHit < nHits; ++iHit){

          // General geometry of the track
          float x = evt->trkxyz_pandoraTrack[iTrk][iPlane][iHit][0];
          float y = evt->trkxyz_pandoraTrack[iTrk][iPlane][iHit][1];
          float z = evt->trkxyz_pandoraTrack[iTrk][iPlane][iHit][2];
          float t = x * evtProc.kXtoT;
          
          // Check if x is lower or higher than the APA bounds, charge seems to accumulate there
          if(x < evtProc.APA_X_POSITIONS[0] || x > evtProc.APA_X_POSITIONS[2]) continue;

          // Lifetime correction
          int tpc =evtProc.WhichTPC(x) + 1;
          float dx = ( -1 + 2*(tpc%2) )*(x - evtProc.APA_X_POSITIONS[tpc/2]);
          float dt = dx*evtProc.kXtoT;
          float corr      = TMath::Exp(-dt/2.88);
          float eCorr     = TMath::Exp(-dt/2.88) / TMath::Exp(-dt/3.); // Correct for the already-corrected energy

          // New values
          float dEdxVal   = dEdx.at(iHit);
          float dQdxVal   = dQdx.at(iHit);
          float dQdxCorr  = dQdxVal/corr;
          float dEdxCorr  = dEdxVal/eCorr;
          float dEdQVal   = dEdxVal/dQdxCorr;
          float dEdQCorr  = dEdxCorr/dQdxCorr;

          h_dedx_x->Fill(x,dEdxVal);
          h_dqdx_x->Fill(x,dQdxVal);
          h_corr_dqdx_x->Fill(x,dQdxCorr);
          h_corr_dedq_x->Fill(x,dEdQVal);
          h_corr2_dedq_x->Fill(x,dEdQCorr);
          
          h_hits_xy->Fill(x,y);
          h_hits_xz->Fill(x,z);
          h_hits_yz->Fill(y,z);
          h_hits_xyz->Fill(x,y,z);
        } // Hits
      } // Planes
    } // Tracks
  }// Event loop
  std::cout << " --- 100 % --- |" << std::endl;

  std::cout << " Number of times max hits exceeds limit: " << maxHitsLimit << std::endl;

  // Sort out the TeX file
  std::vector<std::string> contents{
    "Events",
    "Tracks",
    "Muons",
    "$\\mu > 3$m",
    "$\\mu > 3$m \\& $y_{i}$ > 599.5 cm",
    "Crosses top or bottom",
    "Crosses top and bottom",
    "Crosses $ \\geq $ 1 APA/CPA",
    "Crosses $ \\geq $ 2 APA/CPA",
    "Stopping"
  };
  std::vector<unsigned int> rates{
    nEvts,
    totalTracks,
    nMu,
    nLongTracks,
    nLongHighYTracks,
    topOrBottom,
    topBottom,
    min1APACPA,
    min2APACPA,
    stopping
  };

  ofstream texFile;
  texFile.open(location+"sample_contents_events"+tag+".tex");
  WriteStatsToTeX(texFile, n, contents, rates, static_cast<double>(nEvts), "All Events");

  TCanvas *c1 = new TCanvas("c1","",1000,800);
  SetCanvasStyle(c1, 0.1,0.12,0.05,0.12,0,0,0);

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

  // dEdx vs x
  SetHistogramStyle2D(h_dedx_x,"x [cm]", " dE/dx [MeV/cm]",false);
  h_dedx_x->Draw("colz");

  // Draw the APA and CPA lines and labels
  for(unsigned int iLine = 0; iLine < APACPALines.size(); ++iLine){
    APACPALines.at(iLine)->SetY2(10);
    APACPALines.at(iLine)->Draw();
  }

  FormatLatex(evtProc.APA_X_POSITIONS[0]+10,9.3, "#color[0]{APA}");
  FormatLatex(evtProc.CPA_X_POSITIONS[0]+10,9.3, "#color[0]{CPA}");
  FormatLatex(evtProc.APA_X_POSITIONS[1]+10,9.3, "#color[0]{APA}");
  FormatLatex(evtProc.CPA_X_POSITIONS[1]+10,9.3, "#color[0]{CPA}");
  FormatLatex(evtProc.APA_X_POSITIONS[2]+10,9.3, "#color[0]{APA}");

  c1->SaveAs((location+"/dEdx_vs_X"+tag+".png").c_str());
  c1->SaveAs((location+"/dEdx_vs_X"+tag+".root").c_str());
  c1->Clear();

  // Charge
  SetHistogramStyle2D(h_dqdx_x,"x [cm]", " Charge deposition [ADC/cm]",false);
  h_dqdx_x->Draw("colz");

  // Draw the APA and CPA lines and labels
  for(unsigned int iLine = 0; iLine < APACPALines.size(); ++iLine){
    APACPALines.at(iLine)->SetY2(1000);
    APACPALines.at(iLine)->Draw();
  }

  FormatLatex(evtProc.APA_X_POSITIONS[0]+10,900, "#color[0]{APA}");
  FormatLatex(evtProc.CPA_X_POSITIONS[0]+10,900, "#color[0]{CPA}");
  FormatLatex(evtProc.APA_X_POSITIONS[1]+10,900, "#color[0]{APA}");
  FormatLatex(evtProc.CPA_X_POSITIONS[1]+10,900, "#color[0]{CPA}");
  FormatLatex(evtProc.APA_X_POSITIONS[2]+10,900, "#color[0]{APA}");

  c1->SaveAs((location+"/charge_vs_X"+tag+".png").c_str());
  c1->SaveAs((location+"/charge_vs_X"+tag+".root").c_str());
  c1->Clear();

  SetHistogramStyle2D(h_corr_dqdx_x,"x [cm]", " Charge deposition [ADC/cm]", false);
  h_corr_dqdx_x->Draw("colz");

  // Draw the APA and CPA lines and labels
  for(unsigned int iLine = 0; iLine < APACPALines.size(); ++iLine){
    APACPALines.at(iLine)->SetY2(1000);
    APACPALines.at(iLine)->Draw();
  }

  FormatLatex(evtProc.APA_X_POSITIONS[0]+10,900, "#color[0]{APA}");
  FormatLatex(evtProc.CPA_X_POSITIONS[0]+10,900, "#color[0]{CPA}");
  FormatLatex(evtProc.APA_X_POSITIONS[1]+10,900, "#color[0]{APA}");
  FormatLatex(evtProc.CPA_X_POSITIONS[1]+10,900, "#color[0]{CPA}");
  FormatLatex(evtProc.APA_X_POSITIONS[2]+10,900, "#color[0]{APA}");

  c1->SaveAs((location+"/corr_charge_vs_X"+tag+".png").c_str());
  c1->SaveAs((location+"/corr_charge_vs_X"+tag+".root").c_str());
  c1->Clear();

  SetHistogramStyle2D(h_corr_dedq_x,"x [cm]", " Energy per charge deposition [MeV/ADC]", false);
  h_corr_dedq_x->Draw("colz");

  // Draw the APA and CPA lines and labels
  for(unsigned int iLine = 0; iLine < APACPALines.size(); ++iLine){
    APACPALines.at(iLine)->SetY1(6.6e-3);
    APACPALines.at(iLine)->SetY2(7.2e-3);
    APACPALines.at(iLine)->Draw();
  }

  FormatLatex(evtProc.APA_X_POSITIONS[0]+10,7.1e-3, "#color[0]{APA}");
  FormatLatex(evtProc.CPA_X_POSITIONS[0]+10,7.1e-3, "#color[0]{CPA}");
  FormatLatex(evtProc.APA_X_POSITIONS[1]+10,7.1e-3, "#color[0]{APA}");
  FormatLatex(evtProc.CPA_X_POSITIONS[1]+10,7.1e-3, "#color[0]{CPA}");
  FormatLatex(evtProc.APA_X_POSITIONS[2]+10,7.1e-3, "#color[0]{APA}");

  c1->SaveAs((location+"/corr_energy_charge_vs_X"+tag+".png").c_str());
  c1->SaveAs((location+"/corr_energy_charge_vs_X"+tag+".root").c_str());
  c1->Clear();

  SetHistogramStyle2D(h_corr2_dedq_x,"x [cm]", " Energy per charge deposition [MeV/ADC]", false);
  h_corr2_dedq_x->Draw("colz");

  // Draw the APA and CPA lines and labels
  for(unsigned int iLine = 0; iLine < APACPALines.size(); ++iLine){
    APACPALines.at(iLine)->SetY1(6.6e-3);
    APACPALines.at(iLine)->SetY2(8e-3);
    APACPALines.at(iLine)->Draw();
  }

  FormatLatex(evtProc.APA_X_POSITIONS[0]+10,7.9e-3, "#color[0]{APA}");
  FormatLatex(evtProc.CPA_X_POSITIONS[0]+10,7.9e-3, "#color[0]{CPA}");
  FormatLatex(evtProc.APA_X_POSITIONS[1]+10,7.9e-3, "#color[0]{APA}");
  FormatLatex(evtProc.CPA_X_POSITIONS[1]+10,7.9e-3, "#color[0]{CPA}");
  FormatLatex(evtProc.APA_X_POSITIONS[2]+10,7.9e-3, "#color[0]{APA}");

  c1->SaveAs((location+"/corr2_energy_charge_vs_X"+tag+".png").c_str());
  c1->SaveAs((location+"/corr2_energy_charge_vs_X"+tag+".root").c_str());
  c1->Clear();

  // Number of hits XY
  TCanvas *c2 = new TCanvas("c2","",1000,800);
  SetCanvasStyle(c2, 0.1,0.12,0.05,0.12,0,0,0);
  SetHistogramStyle2D(h_hits_xy,"x [cm]", " y [cm]", false);
  h_hits_xy->Draw("colz");
  c2->SaveAs((location+"/xy_hits"+tag+".png").c_str());
  c2->SaveAs((location+"/xy_hits"+tag+".root").c_str());
  c2->Clear();

  // Number of hits XY
  SetHistogramStyle2D(h_hits_xz,"x [cm]"," z [cm]", false);
  h_hits_xz->Draw("colz");
  c2->SaveAs((location+"/xz_hits"+tag+".png").c_str());
  c2->SaveAs((location+"/xz_hits"+tag+".root").c_str());
  c2->Clear();

  // Number of hits YZ
  SetHistogramStyle2D(h_hits_yz,"y [cm]"," z [cm]", false);
  h_hits_yz->Draw("colz");
  c2->SaveAs((location+"/yz_hits"+tag+".png").c_str());
  c2->SaveAs((location+"/yz_hits"+tag+".root").c_str());
  c2->Clear();

  TCanvas *c3 = new TCanvas("c3","",1000,800);
  SetCanvasStyle(c3, 0.1,0.12,0.05,0.12,0,0,0);
  SetHistogramStyle3D(h_hits_xyz,"x [cm]","y [cm]","z [cm]");
  h_hits_xyz->SetMarkerStyle(33);
  h_hits_xyz->SetMarkerColor(kViolet-5);
  h_hits_xyz->Draw();
  c3->SaveAs((location+"/xyz_hits"+tag+".png").c_str());
  c3->SaveAs((location+"/xyz_hits"+tag+".root").c_str());
  
  // Plane crossing
  TCanvas *c4 = new TCanvas("c4","",900,900);
  SetCanvasStyle(c4, 0.12,0.05,0.06,0.15,0,0,0);

  TLegend *l = new TLegend(0.22,0.94,0.98,0.995);
  l->SetNColumns(3);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->SetTextFont(132);

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
  h_plane_enter->SetLineColor(kViolet-5);
  h_plane_exit->SetLineColor(kTeal-5);
  h_plane_cross->SetLineColor(kOrange+5);
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
  h_n_crossed->SetLineColor(kTeal-5);
  h_n_crossed->GetYaxis()->SetTitleOffset(0.95);
  c5->SaveAs((location+"/tracks_crossed_nplanes"+tag+".png").c_str());
  c5->SaveAs((location+"/tracks_crossed_nplanes"+tag+".root").c_str());
  c5->Clear();
  
  TCanvas *c6 = new TCanvas("c6","",900,900);
  SetCanvasStyle(c6, 0.12,0.05,0.06,0.12,0,0,0);

  SetHistogramStyle1D(h_enter_dist,"Distance from candidate entrance/exit [cm]", " Rate");
  h_enter_dist->Draw("hist");
  h_exit_dist->Draw("same");
  h_enter_dist->SetLineWidth(2);
  h_exit_dist->SetLineWidth(2);
  h_enter_dist->SetLineStyle(2);
  h_exit_dist->SetLineStyle(3);
  h_enter_dist->SetLineColor(kViolet-5);
  h_exit_dist->SetLineColor(kTeal-5);
  h_enter_dist->GetYaxis()->SetTitleOffset(0.95);

  l->SetNColumns(2);
  l->SetX1NDC(0.47);
  l->AddEntry(h_enter_dist, "Entrance", "l");
  l->AddEntry(h_exit_dist, "Exit", "l");
  l->Draw();

  c6->SaveAs((location+"/distance_to_entrance_exit_planes"+tag+".png").c_str());
  c6->SaveAs((location+"/distance_to_entrance_exit_planes"+tag+".root").c_str());
  c6->Clear();
  
  SetHistogramStyle1D(h_muon_length,"Muon length [cm]", " Rate");
  h_muon_length->Draw("hist");
  h_muon_length->SetLineWidth(2);
  h_muon_length->SetLineColor(kViolet-5);
  h_muon_length->GetYaxis()->SetTitleOffset(0.95);
  c6->SaveAs((location+"/muon_length"+tag+".png").c_str());
  c6->SaveAs((location+"/muon_length"+tag+".root").c_str());
  c6->Clear();
  
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
