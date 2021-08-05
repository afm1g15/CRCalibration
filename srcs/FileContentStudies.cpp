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

// Location to save plots from this macro:
std::string location="/home/jones/work/cosmics/LArSoft-v08_50_00/work/plots/v08_50_00/dEdxCalib/";

// Allowed branches to read from the tree
std::vector<TString> allowed = {
   "run",
   "event",
   "ntracks_pandoraTrack",
   "ntrkhits_pandoraTrack",
   "trkdqdx_pandoraTrack",
   "trkdedx_pandoraTrack",
   "trkxyz_pandoraTrack",
   "trklen_pandoraTrack",
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
  std::vector<double> minx_fid, miny_fid, minz_fid;
  std::vector<double> maxx_fid, maxy_fid, maxz_fid;
  std::vector<double> minx_av, miny_av, minz_av;
  std::vector<double> maxx_av, maxy_av, maxz_av;

  p->getValue("InputList", input_list);
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

  std::cout << " Total number of planes in the active volume of the DUNE SP module: " << allPlanes.size() << std::endl;
  std::cout << " Corresponds to " << active.GetNTPCs() << " TPC's in the DUNE SP module" << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;
  
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
  TH2D *h_dedx_x   = new TH2D("h_dedx_x","",100,-800,800,100,0,10);
  TH2D *h_hits_xy  = new TH2D("h_hits_xy","",100,-800,800,100,-650,650);
  TH2D *h_hits_xz  = new TH2D("h_hits_xz","",100,-800,800,300,-200,6000);
  TH2D *h_hits_yz  = new TH2D("h_hits_yz","",100,-700,700,300,-200,6000);
  TH3D *h_hits_xyz = new TH3D("h_hits_xyz","",100,-800,800,100,-700,700,300,-200,6000);
  
  // Setup counters
  unsigned int maxHitsLimit = 0;

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
      // Length cuts (2m)
      if(!evtProc.SelectTrack(evt,iTrk)) continue;

      // Get the best plane
      unsigned int bestPlane = 0;
      int currHits  = -999;
      for(unsigned int iPlane = 0; iPlane < 3; ++iPlane){
        if(evt->ntrkhits_pandoraTrack[iTrk][iPlane] > currHits){
          currHits  = evt->ntrkhits_pandoraTrack[iTrk][iPlane];
          bestPlane = iPlane; 
        } // CurrHits
      } // Planes

      // Now fill dQ/dx and dE/dx and hit histograms for each of the three wire planes
      // Somehow flag the best plane histogram
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
