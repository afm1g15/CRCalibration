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
   "TrackId",
   "Mother",
   "Eng",
   "pathlen",
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
   "trkId_pandoraTrack",
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
   "trklen_pandoraTrack",
   "trkke_pandoraTrack"
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
     
int pandoraStudies(const char *config){

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
  // 1D
  TH1D *h_nHits_bp     = new TH1D("h_nHits_bp","",100,0,5000);
  TH1D *h_reco_length  = new TH1D("h_reco_length","",100,0,2000);
  TH1D *h_true_length  = new TH1D("h_true_length","",100,0,2000);
  TH1D *h_reco_energy  = new TH1D("h_reco_energy","",100,0.001,3000);
  TH1D *h_true_energy  = new TH1D("h_true_energy","",100,0.001,3000);

  // Log scale for energies
  SetLogX(h_reco_energy);
  SetLogX(h_true_energy);

  // 2D
  TH2D *h_hits_xy      = new TH2D("h_hits_xy","",200,-800,800,200,-650,650);
  TH2D *h_hits_xz      = new TH2D("h_hits_xz","",200,-800,800,200,0,6000);
  TH2D *h_hits_yz      = new TH2D("h_hits_yz","",200,-700,700,200,0,6000);
  
  // Sort out log scales if needed 
  
  // Setup counters
  unsigned int nTracks          = 0;
  unsigned int nMu              = 0;
  unsigned int nLongTracks      = 0;
  unsigned int nLongHighYTracks = 0;
  unsigned int topBottom        = 0;
  unsigned int topOrBottom      = 0;
  unsigned int min2APACPA       = 0;
  unsigned int min1APACPA       = 0;
  unsigned int stopping         = 0;
  unsigned int exiting          = 0;

  // Now loop over the events
  unsigned int nEvts = tree->GetEntries();
  unsigned int iIt = 1;

  std::cout << " |";
  for(unsigned int iEvt = 0; iEvt < nEvts; ++iEvt){
    tree->GetEntry(iEvt);
    if(!evtProc.SelectEvent(evt)) continue;
    unsigned int nTrks = evt->ntracks_pandoraTrack;
    int nGeant = evt->geant_list_size;
    
    // Print the processing rate
    double evtFrac  = iEvt/static_cast<double>(nEvts);

    if(std::abs(0.1*iIt-evtFrac) < std::numeric_limits<double>::epsilon()){
      std::cout << " --- " << evtFrac*100 << " %";
      std::cout.flush();
      iIt++;
    }

    // Now loop over the tracks so we can do some stuff!!
    for(unsigned int iTrk = 0; iTrk < nTrks; ++iTrk){
      
      // Increment and initialise counters
      unsigned int nHits_BP = 0;
      nTracks++;

      // Get the best plane
      int bestPlane = 0;
      int currHits  = -999;
      for(int iPlane = 0; iPlane < 3; ++iPlane){
        if(evt->ntrkhits_pandoraTrack[iTrk][iPlane] > currHits){
          currHits  = evt->ntrkhits_pandoraTrack[iTrk][iPlane];
          bestPlane = iPlane; 
        } // CurrHits
      } // Planes

      // Only look at primary muons
      if(abs(evt->trkpdgtruth_pandoraTrack[iTrk][bestPlane]) != 13) continue;
      nMu++;
      
      // Length cuts (3m)
      if(evtProc.SelectTrack(evt,iTrk)) nLongTracks++;
      
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

      if(evtProc.SelectTrack(evt,iTrk) && startVtx.Y() < 599.5) nLongHighYTracks++;

      // Parameters of interest
      int ID         = evt->trkId_pandoraTrack[iTrk];
      float reco_len = evt->trklen_pandoraTrack[iTrk];
      float reco_eng = evt->trkke_pandoraTrack[bestPlane][iTrk]/1000.;

      // Find the closest plane to the start vertex and count it as a crossing plane
      Plane enteringPlane = GetClosestPlane(extPlanes, startVtx, endVtx);
      double distFromEntrance = GetDistanceToPlane(enteringPlane, startVtx, endVtx);

      Plane exitingPlane = GetClosestPlane(extPlanes, endVtx, startVtx);
      double distFromExit = GetDistanceToPlane(exitingPlane, endVtx, startVtx);

      // Now determine if the current track crossed each detector plane individually
      unsigned int planeN    = 0;
      unsigned int extPlaneN = 0;
      
      // Counter for the number of planes this track has crossed
      unsigned int nPlanesCrossed = 0;
      unsigned int nExtCrossed    = 0;

      // Setup list of plane labels the track has crossed
      std::vector<std::string> labelsCrossed;
      
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
            nPlanesCrossed++;
            labelsCrossed.push_back(pl.GetLabel());
          }
        }
        else if(exitingPlane.GetLabel() == pl.GetLabel()){
          if(distFromExit < 1){
            nPlanesCrossed++;
            labelsCrossed.push_back(pl.GetLabel());
          }
        }
        // Otherwise, check if the track intersects the current plane
        else if(CheckIfIntersectsPlane(pl,startVtx,endVtx,reco_len)){
          nPlanesCrossed++;
          labelsCrossed.push_back(pl.GetLabel());
        } // Intersects
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
        else if(CheckIfIntersectsPlane(pl,startVtx,endVtx,reco_len)){
          nExtCrossed++;
        } // Intersects
      } // Planes
      
      if(yCut == 1 && startVtx.Y() < 599.5) continue; // Make sure the tracks start at the top of the detector

      // Now count the planes crossed for the stats table
      bool thruGoing = false;
      if(nExtCrossed == 1){
        stopping++;
      }
      if(nExtCrossed >= 2){
        thruGoing = true;
        exiting++;
      }

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
      if(top && bot){
        topBottom++;
      }

      // the following studies should be conducted with top-bottom muons to start with
      if(thru == 1 && !thruGoing) continue;

      // Now fill dQ/dx and dE/dx and hit histograms for each of the three wire planes
      // Somehow flag the best wire plane histogram
      for(int iPlane = 0; iPlane < 3; ++iPlane){

        unsigned int nHits = evt->ntrkhits_pandoraTrack[iTrk][iPlane];
        
        // Fill counter for best plane number of hits
        if(iPlane == bestPlane)
          nHits_BP = nHits;

        // Get the associated true energy of the muon
        // Try to get the true energy
        // Get the list iterator from matching ID's
        float true_eng = -1.;
        float true_len = -1.;
        for(int iG4 = 0; iG4 < nGeant; ++iG4){
          int trueID = evt->TrackId[iG4];

          if(evt->trkidtruth_pandoraTrack[iTrk][bestPlane] == trueID){
            true_eng = evt->Eng[iG4];
            true_len = evt->pathlen[iG4];

            break;
          }
        }
        if(true_eng < 0){
          std::cout << " Warning: Energy is below zero, skipping track with energy: " << true_eng << std::endl;
          continue;
        }
     
        // Make sure it doesn't exceed the maximum size of the array
        // Count if it does so we can see how often it happens
        if(nHits > MAX_TRACK_HITS){
          nHits = MAX_TRACK_HITS;
        }
      
        // Fill track-level histograms
        h_reco_length->Fill(reco_len);
        h_reco_energy->Fill(reco_eng);
        h_true_length->Fill(true_len);
        h_true_energy->Fill(true_eng);

        // Now loop over hits so we can work our calo magic
        for(unsigned int iHit = 0; iHit < nHits; ++iHit){

          // General geometry of the track
          float x = evt->trkxyz_pandoraTrack[iTrk][iPlane][iHit][0];
          float y = evt->trkxyz_pandoraTrack[iTrk][iPlane][iHit][1];
          float z = evt->trkxyz_pandoraTrack[iTrk][iPlane][iHit][2];

          h_hits_xy->Fill(x,y);
          h_hits_xz->Fill(x,z);
          h_hits_yz->Fill(y,z);
          
        } // Hits
      } // Planes
      h_nHits_bp->Fill(nHits_BP);
    } // Tracks
  }// Event loop
  std::cout << " --- 100 % --- |" << std::endl;

  // Sort out the TeX file
  std::vector<std::string> contents{
    "Events",
    "Tracks",
    "Muons",
    "$\\mu > 3$~m",
    "$\\mu > 3$~m \\& $y_{i} > 599.5$~cm",
    "Crosses top or bottom",
    "Crosses top and bottom",
    "Crosses $ \\geq $ 1 APA/CPA",
    "Crosses $ \\geq $ 2 APA/CPA",
    "Stopping",
    "Exiting"
  };
  std::vector<unsigned int> rates{
    nEvts,
    nTracks,
    nMu,
    nLongTracks,
    nLongHighYTracks,
    topOrBottom,
    topBottom,
    min1APACPA,
    min2APACPA,
    stopping,
    exiting
  };

  ofstream texFile;
  texFile.open(location+"sample_contents_events"+tag+".tex");
  WriteStatsToTeX(texFile, n, contents, rates, static_cast<double>(nEvts), "All Events");

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

  // -------------------------------------------------------------------------
  //                           2D CANVAS
  // -------------------------------------------------------------------------
  // Number of hits XY
  TCanvas *c1 = new TCanvas("c1","",1000,800);
  c1->SetLogz();
  SetCanvasStyle(c1, 0.1,0.12,0.05,0.12,0,0,0);
  SetHistogramStyle2D(h_hits_xy,"x [cm]", " y [cm]", false);
  h_hits_xy->GetZaxis()->SetLabelSize(0.03);
  h_hits_xy->GetZaxis()->SetLabelFont(132);
  h_hits_xy->Draw("colz");
  
  // Draw the APA and CPA lines and labels
  for(unsigned int iLine = 0; iLine < APACPALines.size(); ++iLine){
    APACPALines.at(iLine)->SetY1(-650);
    APACPALines.at(iLine)->SetY2(650);
    APACPALines.at(iLine)->Draw();
  }

  FormatLatex(evtProc.APA_X_POSITIONS[0]+10,500, "#color[0]{APA}");
  FormatLatex(evtProc.CPA_X_POSITIONS[0]+10,500, "#color[0]{CPA}");
  FormatLatex(evtProc.APA_X_POSITIONS[1]+10,500, "#color[0]{APA}");
  FormatLatex(evtProc.CPA_X_POSITIONS[1]+10,500, "#color[0]{CPA}");
  FormatLatex(evtProc.APA_X_POSITIONS[2]-180,500, "#color[0]{APA}");

  c1->SaveAs((location+"/xy_hits"+tag+".png").c_str());
  c1->SaveAs((location+"/xy_hits"+tag+".root").c_str());
  c1->Clear();

  // Number of hits XY
  SetHistogramStyle2D(h_hits_xz,"x [cm]"," z [cm]", false);
  h_hits_xz->GetZaxis()->SetLabelSize(0.03);
  h_hits_xz->GetZaxis()->SetLabelFont(132);
  h_hits_xz->Draw("colz");
  
  // Draw the APA and CPA lines and labels
  for(unsigned int iLine = 0; iLine < APACPALines.size(); ++iLine){
    APACPALines.at(iLine)->SetY1(0);
    APACPALines.at(iLine)->SetY2(6000);
    APACPALines.at(iLine)->Draw();
  }

  FormatLatex(evtProc.APA_X_POSITIONS[0]+10,5500, "#color[0]{APA}");
  FormatLatex(evtProc.CPA_X_POSITIONS[0]+10,5500, "#color[0]{CPA}");
  FormatLatex(evtProc.APA_X_POSITIONS[1]+10,5500, "#color[0]{APA}");
  FormatLatex(evtProc.CPA_X_POSITIONS[1]+10,5500, "#color[0]{CPA}");
  FormatLatex(evtProc.APA_X_POSITIONS[2]-180,5500, "#color[0]{APA}");

  c1->SaveAs((location+"/xz_hits"+tag+".png").c_str());
  c1->SaveAs((location+"/xz_hits"+tag+".root").c_str());
  c1->Clear();

  // Number of hits YZ
  SetHistogramStyle2D(h_hits_yz,"y [cm]"," z [cm]", false);
  h_hits_yz->GetZaxis()->SetLabelSize(0.03);
  h_hits_yz->GetZaxis()->SetLabelFont(132);
  h_hits_yz->Draw("colz");
  
  // Draw the APA and CPA lines and labels
  for(unsigned int iLine = 0; iLine < APACPALines.size(); ++iLine){
    APACPALines.at(iLine)->Draw();
  }

  FormatLatex(evtProc.APA_X_POSITIONS[0]+10,5500, "#color[0]{APA}");
  FormatLatex(evtProc.CPA_X_POSITIONS[0]+10,5500, "#color[0]{CPA}");
  FormatLatex(evtProc.APA_X_POSITIONS[1]+10,5500, "#color[0]{APA}");
  FormatLatex(evtProc.CPA_X_POSITIONS[1]+10,5500, "#color[0]{CPA}");
  FormatLatex(evtProc.APA_X_POSITIONS[2]-180,5500, "#color[0]{APA}");

  c1->SaveAs((location+"/yz_hits"+tag+".png").c_str());
  c1->SaveAs((location+"/yz_hits"+tag+".root").c_str());
  c1->Clear();
  
  // -------------------------------------------------------------------------
  //                           1D NO LOG CANVAS
  // -------------------------------------------------------------------------
  TCanvas *c2 = new TCanvas("c2","",900,900);
  SetCanvasStyle(c2, 0.12,0.08,0.06,0.12,0,0,0);
  
  SetHistogramStyle1D(h_nHits_bp,"Number of hits on the best plane", "Rate");
  h_nHits_bp->Draw("hist");
  h_nHits_bp->SetLineWidth(3);
  h_nHits_bp->SetLineColor(kTeal-5);
  c2->SaveAs((location+"/nHits_bp"+tag+".png").c_str());
  c2->SaveAs((location+"/nHits_bp"+tag+".root").c_str());
  c2->Clear();

  h_nHits_bp->Scale(1./static_cast<double>(h_nHits_bp->Integral()));
  h_nHits_bp->Draw("hist");
  h_nHits_bp->SetLineWidth(3);
  h_nHits_bp->SetLineColor(kTeal-5);
  c2->SaveAs((location+"/nHits_bp"+tag+".png").c_str());
  c2->SaveAs((location+"/nHits_bp"+tag+".root").c_str());
  c2->Clear();
 
  SetHistogramStyle1D(h_reco_length,"Track length [cm]", "Rate");
  SetHistogramStyle1D(h_true_length,"Track length [cm]", "Rate");
  
  h_reco_length->Draw("hist");
  h_reco_length->SetLineWidth(3);
  h_reco_length->SetLineColor(kViolet-5);

  c2->SaveAs((location+"/track_reco_length"+tag+".png").c_str());
  c2->SaveAs((location+"/track_reco_length"+tag+".root").c_str());
  c2->Clear();
  
  h_true_length->Draw("hist");
  h_true_length->SetLineWidth(3);
  h_true_length->SetLineColor(kTeal-5);
  
  c2->SaveAs((location+"/track_true_length"+tag+".png").c_str());
  c2->SaveAs((location+"/track_true_length"+tag+".root").c_str());
  c2->Clear();
  
  h_reco_length->Scale(1./static_cast<double>(h_reco_length->Integral()));
  h_true_length->Scale(1./static_cast<double>(h_true_length->Integral()));
  double maxy = std::max(h_reco_length->GetMaximum(),h_true_length->GetMaximum())*1.1;

  h_reco_length->Draw("hist");
  h_true_length->Draw("hist same");
  h_reco_length->SetLineWidth(3);
  h_reco_length->SetLineColor(kViolet-5);
  h_true_length->SetLineWidth(3);
  h_true_length->SetLineColor(kTeal-5);
  h_reco_length->GetYaxis()->SetRangeUser(0,maxy);

  TLegend *l = new TLegend(0.52,0.94,0.98,0.995);
  l->SetNColumns(2);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->SetTextFont(132);
  l->SetTextSize(0.035);
  
  l->AddEntry(h_true_length, "Truth", "l");
  l->AddEntry(h_reco_length, "Reconstruction", "l");
  l->Draw();
  
  c2->SaveAs((location+"/track_length"+tag+".png").c_str());
  c2->SaveAs((location+"/track_length"+tag+".root").c_str());
  c2->Clear();
  l->Clear();
  
  // -------------------------------------------------------------------------
  //                           1D LOG CANVAS
  // -------------------------------------------------------------------------
  TCanvas *c3 = new TCanvas("c3","",900,900);
  SetCanvasStyle(c3, 0.12,0.08,0.06,0.12,0,0,0);
  
  SetHistogramStyle1D(h_reco_energy,"Track energy [GeV]", "Rate");
  
  h_reco_energy->Draw("hist");
  h_reco_energy->SetLineWidth(3);
  h_reco_energy->SetLineColor(kViolet-5);

  c3->SaveAs((location+"/track_reco_energy"+tag+".png").c_str());
  c3->SaveAs((location+"/track_reco_energy"+tag+".root").c_str());
  c3->Clear();
  
  SetHistogramStyle1D(h_true_energy,"Track energy [GeV]", "Rate");
  
  h_true_energy->Draw("hist");
  h_true_energy->SetLineWidth(3);
  h_true_energy->SetLineColor(kTeal-5);

  c3->SaveAs((location+"/track_true_energy"+tag+".png").c_str());
  c3->SaveAs((location+"/track_true_energy"+tag+".root").c_str());
  c3->Clear();
  
  h_reco_energy->Scale(1./static_cast<double>(h_reco_energy->Integral()), "width");
  h_true_energy->Scale(1./static_cast<double>(h_true_energy->Integral()), "width");
  maxy = std::max(h_reco_energy->GetMaximum(),h_true_energy->GetMaximum())*1.1;
  
  h_reco_energy->Draw("hist");
  h_true_energy->Draw("hist same");
  h_reco_energy->GetYaxis()->SetRangeUser(1,maxy);

  l->AddEntry(h_true_energy, "Truth", "l");
  l->AddEntry(h_reco_energy, "Reconstruction", "l");
  l->Draw();
  
  c3->SetTopMargin(0.06);
  c3->SetLogx();
  c3->SetLogy();
  c3->SaveAs((location+"/track_energy"+tag+".png").c_str());
  c3->SaveAs((location+"/track_energy"+tag+".root").c_str());
  c3->Clear();
  l->Clear();

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
