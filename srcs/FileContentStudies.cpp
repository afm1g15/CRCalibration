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
   "trkpdgtruth_pandoraTrack",
   "trkg4id_pandoraTrack",
   "process_primary",
   "ntracks_pandoraTrack",
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
  {"t0", "Top 1"},
  {"t1", "Top 2"},
  {"t2", "Top 3"},
  {"t3", "Top 4"},
  {"bo0", "Bot. 1"},
  {"bo1", "Bot. 2"},
  {"bo2", "Bot. 3"},
  {"bo3", "Bot. 4"},
  {"f0", "Fro. 1"},
  {"f1", "Fro. 2"},
  {"f2", "Fro. 3"},
  {"f3", "Fro. 4"},
  {"ba0", "Bck 1"},
  {"ba1", "Bck 2"},
  {"ba2", "Bck 3"},
  {"ba3", "Bck 4"}
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

  std::cout << " Total number of planes in the active volume of the DUNE SP module: " << allPlanes.size() << std::endl;
  std::cout << " Consisting of " << extPlanes.size() << " external planes and " << intPlanes.size() << " internal planes with labels:" << std::endl; 
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
  
  // Start of analysis (loop over chain and events
  std::cout << " Running analysis..." << std::endl;

  // Then setup the histograms, counters and any other variables to add to
  // Setup histograms
  TH2D *h_dedx_x      = new TH2D("h_dedx_x","",100,-800,800,100,0,10);
  TH2D *h_hits_xy     = new TH2D("h_hits_xy","",100,-800,800,100,-650,650);
  TH2D *h_hits_xz     = new TH2D("h_hits_xz","",100,-800,800,300,-200,6000);
  TH2D *h_hits_yz     = new TH2D("h_hits_yz","",100,-700,700,300,-200,6000);
  TH3D *h_hits_xyz    = new TH3D("h_hits_xyz","",100,-800,800,100,-700,700,300,-200,6000);
  TH1D *h_plane_cross = new TH1D("h_plane_cross","",21,0,21); // Number of tracks crossing each plane
  TH1D *h_plane_enter = new TH1D("h_plane_enter","",21,0,21); // Number of tracks entering from each external plane
  TH1D *h_enter_dist  = new TH1D("h_enter_dist","",100,0,200); // Number of tracks entering from each external plane
  TH1D *h_n_crossed   = new TH1D("h_n_crossed","",21,0,21); // Number of planes crossed by each track
  
  // Setup counters
  unsigned int maxHitsLimit = 0;
  unsigned int wrongWay     = 0;
  unsigned int noPlanes     = 0;
  unsigned int totalTracks  = 0;
  unsigned int nLongTracks  = 0;
  unsigned int nPrimaryMu   = 0;
  unsigned int topBottom    = 0;
  unsigned int topOrBottom  = 0;
  unsigned int min2APACPA   = 0;
  unsigned int min1APACPA   = 0;
  unsigned int stopping     = 0;

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

      // Count tracks
      totalTracks++;
      
      // Length cuts (2m)
      if(!evtProc.SelectTrack(evt,iTrk)) continue;
      nLongTracks++;
      
      // Get the best plane
      unsigned int bestPlane = 0;
      int currHits  = -999;
      for(unsigned int iPlane = 0; iPlane < 3; ++iPlane){
        if(evt->ntrkhits_pandoraTrack[iTrk][iPlane] > currHits){
          currHits  = evt->ntrkhits_pandoraTrack[iTrk][iPlane];
          bestPlane = iPlane; 
        } // CurrHits
      } // Planes

      // Only look at primary muons
      if(abs(evt->trkpdgtruth_pandoraTrack[iTrk][bestPlane]) != 13) continue;
      nPrimaryMu++;

      // Get the track geometry
      TVector3 startVtx(evt->trkstartx_pandoraTrack[iTrk],
                        evt->trkstarty_pandoraTrack[iTrk],
                        evt->trkstartz_pandoraTrack[iTrk]);
      TVector3 endVtx(evt->trkendx_pandoraTrack[iTrk],
                      evt->trkendy_pandoraTrack[iTrk],
                      evt->trkendz_pandoraTrack[iTrk]);

      // Since we are generating downwards-going tracks - if the start y < end y then 
      // assume the reconstruction has got them the wrong way around and flip them
      /*
      if(startVtx.Y() < endVtx.Y()){
        wrongWay++;
        TVector3 temp(endVtx);
        endVtx = startVtx;
        startVtx = temp;
      }*/

      float length = evt->trklen_pandoraTrack[iTrk];

      // Find the closest plane to the start vertex and count it as a crossing plane
      Plane enteringPlane = GetClosestPlane(extPlanes, startVtx, endVtx);
      double distFromEntrance = GetDistanceToPlane(enteringPlane, startVtx, endVtx);
      h_enter_dist->Fill(distFromEntrance);

      // Now determine if the current track crossed each detector plane individually
      unsigned int planeN    = 0;
      unsigned int extPlaneN = 0;
      
      // Counter for the number of planes this track has crossed
      unsigned int nPlanesCrossed = 0;
      unsigned int nExtCrossed    = 0;

      // Loop over planes
      for(const Plane &pl : allPlanes){
        if(planeN > allPlanes.size()){
          std::cerr << " Error: Somehow the current plane iterator exceeds the number of planes in the detector: " << std::endl;
          std::cerr << " Iterator: " << planeN << " of " << allPlanes.size() << " total possible planes " << std::endl;
          std::exit(1);
        } // Debug
        // Check if the track intersects the current plane
        if(CheckIfIntersectsPlane(pl,startVtx,endVtx,length)){
          h_plane_cross->Fill(planeN);
          nPlanesCrossed++;
          labelsCrossed.push_back(pl.GetLabel());
        } // Intersects
        // Otherwise check if this is the plane it (likely) entered the detector through 
        // Determine a maximum allowed distance from the plane to count it as the entrance
        else if(enteringPlane.GetLabel() == pl.GetLabel() && distFromEntrance < 100){
          h_plane_cross->Fill(planeN);
          nPlanesCrossed++;
          labelsCrossed.push_back(pl.GetLabel());
        }
        // Now just fill the entrance distribution
        if(enteringPlane.GetLabel() == pl.GetLabel() && distFromEntrance < 100){
          h_plane_enter->Fill(planeN);
        }
        // Sort out the bin label
        h_plane_cross->GetXaxis()->SetBinLabel(planeN+1,planeLabels.find(pl.GetLabel())->second.c_str());
        h_plane_enter->GetXaxis()->SetBinLabel(planeN+1,planeLabels.find(pl.GetLabel())->second.c_str());
        
        planeN++;
      } // Planes
      
      for(const Plane &pl : extPlanes){
        if(CheckIfIntersectsPlane(pl,startVtx,endVtx,length)){
          nExtCrossed++;
        } // Intersects
        else if(enteringPlane.GetLabel() == pl.GetLabel() && distFromEntrance < 100){
          nExtCrossed++;
        } // Intersects
      } // Planes
      
      // Now fill the number of planes crossed histogram
      h_n_crossed->Fill(nPlanesCrossed);
      if(nPlanesCrossed == 0){
        noPlanes++;
        /*
        std::cout << " Start : ( " << startVtx.X() << ", " << startVtx.Y() << ", " << startVtx.Z() << " ) " << std::endl;
        std::cout << " End   : ( " << endVtx.X() << ", " << endVtx.Y() << ", " << endVtx.Z() << " ) " << std::endl;
        std::cout << " Length:   " << length << std::endl;
        std::cin.get();
        */
      }
      if(nExtCrossed == 1)
        stopping++;

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
      if(top && bot)
        topBottom++;

      // Now fill dQ/dx and dE/dx and hit histograms for each of the three wire planes
      // Somehow flag the best wire plane histogram
      for(unsigned int iPlane = 0; iPlane < 3; ++iPlane){

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
        Float_t *rRArr   = evt->trkresrg_pandoraTrack[iTrk][iPlane];

        // Convert them to vectors
        std::vector<float> dEdx(dEdxArr, dEdxArr + nHits);
        std::vector<float> dQdx(dQdxArr, dQdxArr + nHits);

        // Now loop over hits so we can work our calo magic
        for(unsigned int iHit = 0; iHit < nHits; ++iHit){
          double x = evt->trkxyz_pandoraTrack[iTrk][iPlane][iHit][0];
          double y = evt->trkxyz_pandoraTrack[iTrk][iPlane][iHit][1];
          double z = evt->trkxyz_pandoraTrack[iTrk][iPlane][iHit][2];
          double t = x * evtProc.kXtoT;

          h_dedx_x->Fill(x,dEdx.at(iHit));
          h_hits_xy->Fill(x,y);
          h_hits_xz->Fill(x,z);
          h_hits_yz->Fill(y,z);
          h_hits_xyz->Fill(x,y,z);
        } // Hits
      } // Planes
    } // Tracks
  }// Event loop
  std::cout << " --- 100 % --- |" << std::endl;

  // Sort out the TeX file
  std::vector<std::string> contents{
    "Events",
    "Tracks",
    "Tracks $>$ 2m",
    "Long muons",
    "Crosses top or bottom",
    "Crosses top and bottom",
    "Crosses $ \\geq $ 1 APA/CPA",
    "Crosses $ \\geq $ 2 APA/CPA",
    "Stopping"
  };
  std::vector<unsigned int> rates{
    nEvts,
    totalTracks,
    nLongTracks,
    nPrimaryMu,
    topOrBottom,
    topBottom,
    min1APACPA,
    min2APACPA,
    stopping
  };

  ofstream texFile;
  texFile.open(location+"sample_contents"+tag+".tex");
  WriteStatsToTeX(texFile, n, contents, rates, static_cast<double>(totalTracks));

  TCanvas *c1 = new TCanvas("c1","",1000,800);
  SetCanvasStyle(c1, 0.1,0.12,0.05,0.12,0,0,1);
  
  // dEdx vs x
  SetHistogramStyle2D(h_dedx_x,"x [cm]", " dE/dx [MeV/cm]");
  h_dedx_x->Draw("colz");
  
  // Now draw lines and labels where the APA and CPAs are
  for(unsigned int iA = 0; iA < 3; ++iA){
    TLine *l = new TLine(evtProc.APA_X_POSITIONS[iA], 0, evtProc.APA_X_POSITIONS[iA], 10);
    l->SetLineColor(kBlack);
    l->SetLineWidth(3);
    l->SetLineStyle(2);
    l->Draw();
  }
  for(unsigned int iC = 0; iC < 2; ++iC){
    TLine *l = new TLine(evtProc.CPA_X_POSITIONS[iC], 0, evtProc.CPA_X_POSITIONS[iC], 10);
    l->SetLineColor(kWhite);
    l->SetLineWidth(3);
    l->SetLineStyle(2);
    l->Draw();
  }

  FormatLatex(evtProc.APA_X_POSITIONS[0]+10,9.3, "#color[1]{APA}");
  FormatLatex(evtProc.CPA_X_POSITIONS[0]+10,9.3, "#color[0]{CPA}");

  c1->SaveAs((location+"/dEdx_vs_X.png").c_str());
  c1->SaveAs((location+"/dEdx_vs_X.root").c_str());
  c1->Clear();

  // Number of hits XY
  TCanvas *c2 = new TCanvas("c2","",1000,800);
  SetCanvasStyle(c2, 0.1,0.12,0.05,0.12,0,0,0);
  SetHistogramStyle2D(h_hits_xy,"x [cm]", " y [cm]");
  h_hits_xy->Draw("colz");
  c2->SaveAs((location+"/xy_hits.png").c_str());
  c2->SaveAs((location+"/xy_hits.root").c_str());
  c2->Clear();

  // Number of hits XY
  SetHistogramStyle2D(h_hits_xz,"x [cm]"," z [cm]");
  h_hits_xz->Draw("colz");
  c2->SaveAs((location+"/xz_hits.png").c_str());
  c2->SaveAs((location+"/xz_hits.root").c_str());
  c2->Clear();

  // Number of hits YZ
  SetHistogramStyle2D(h_hits_yz,"y [cm]"," z [cm]");
  h_hits_yz->Draw("colz");
  c2->SaveAs((location+"/yz_hits.png").c_str());
  c2->SaveAs((location+"/yz_hits.root").c_str());
  c2->Clear();

  TCanvas *c3 = new TCanvas("c3","",1000,800);
  SetCanvasStyle(c3, 0.1,0.12,0.05,0.12,0,0,0);
  SetHistogramStyle3D(h_hits_xyz,"x [cm]","y [cm]","z [cm]");
  h_hits_xyz->SetMarkerStyle(33);
  h_hits_xyz->SetMarkerColor(kViolet-5);
  h_hits_xyz->Draw();
  c3->SaveAs((location+"/xyz_hits.png").c_str());
  c3->SaveAs((location+"/xyz_hits.root").c_str());
  
  // Plane crossing
  TCanvas *c4 = new TCanvas("c4","",900,900);
  SetCanvasStyle(c4, 0.12,0.05,0.08,0.15,0,0,0);

  TLegend *l = new TLegend(0.22,0.92,0.98,0.98);
  l->SetNColumns(2);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->SetTextFont(132);

  SetHistogramStyle1D(h_plane_cross,"Plane label", " Number of tracks crossing plane");
  h_plane_cross->Draw("hist");
  h_plane_enter->Draw("same");
  h_plane_cross->SetLineWidth(2);
  h_plane_enter->SetLineWidth(2);
  h_plane_cross->SetLineStyle(2);
  h_plane_enter->SetLineStyle(2);
  h_plane_cross->SetLineColor(kViolet-5);
  h_plane_enter->SetLineColor(kTeal-5);
  h_plane_cross->LabelsOption("v");
  h_plane_cross->GetXaxis()->SetTitleOffset(1.4);
  h_plane_cross->GetYaxis()->SetTitleOffset(0.95);

  l->AddEntry(h_plane_cross, "Crossed planes", "l");
  l->AddEntry(h_plane_enter, "Entry planes", "l");
  l->Draw("same");

  c4->SaveAs((location+"/planes_crossed_entered"+tag+".png").c_str());
  c4->SaveAs((location+"/planes_crossed_entered"+tag+".root").c_str());
  c4->Clear();
  
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
  SetCanvasStyle(c6, 0.12,0.05,0.05,0.12,0,0,0);

  SetHistogramStyle1D(h_enter_dist,"Distance from candidate entrance", " Rate");
  h_enter_dist->Draw("hist");
  h_enter_dist->SetLineWidth(2);
  h_enter_dist->SetLineColor(kViolet-5);
  h_enter_dist->GetYaxis()->SetTitleOffset(0.95);
  c6->SaveAs((location+"/distance_to_entrance_plane"+tag+".png").c_str());
  c6->SaveAs((location+"/distance_to_entrance_plane"+tag+".root").c_str());
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
