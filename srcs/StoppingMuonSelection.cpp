/************************************************************************
 * 
 * A macro to understand the truth-level hit and track distributions for 
 * dE/dx calibration studies
 *
 *
 * Example file list located here:
 *   /home/jones/work/cosmics/LArSoft-v08_50_00/work/files/anafiles.list
 *
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
     
int stoppingMuonSelection(const char *config){

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
  TH1D *h_plane_cross   = new TH1D("h_plane_cross","",9,0,9);   // Number of tracks crossing each plane
  TH1D *h_plane_enter   = new TH1D("h_plane_enter","",9,0,9);   // Number of tracks entering from each external plane
  TH1D *h_plane_exit    = new TH1D("h_plane_exit","",9,0,9);    // Number of tracks exiting from each external plane
  TH1D *h_enter_dist    = new TH1D("h_enter_dist","",200,0,10); // Number of tracks entering from each external plane
  TH1D *h_exit_dist     = new TH1D("h_exit_dist","",200,0,10);  // Number of tracks entering from each external plane
  TH1D *h_n_crossed     = new TH1D("h_n_crossed","",9,0,9);     // Number of planes crossed by each track
  TH1D *h_n_ext_crossed = new TH1D("h_n_ext_crossed","",6,0,6); // Number of external planes crossed by each track
  
  // Setup counters
  unsigned int totalTracks = 0;
  unsigned int nMu         = 0;
  unsigned int nLongMu     = 0;
  unsigned int stopping    = 0;
  unsigned int exiting     = 0;
  
  // Now loop over the events
  unsigned int nEvts = tree->GetEntries();
  unsigned int iIt = 1;

  std::cout << " |";
  for(unsigned int iEvt = 0; iEvt < nEvts; ++iEvt){
    tree->GetEntry(iEvt);
    if(!evtProc.SelectEvent(evt)) continue;
    
    // Get the total number of true and reconstructed tracks to loop over
    unsigned int nTrks = evt->ntracks_pandoraTrack;
    int nGeant = evt->geant_list_size;
    
    // Print the processing rate
    double evtFrac  = iEvt/static_cast<double>(nEvts);

    // Prints out how much has been completed so far
    if(std::abs(0.1*iIt-evtFrac) < std::numeric_limits<double>::epsilon()){
      std::cout << " --- " << evtFrac*100 << " %";
      std::cout.flush();
      iIt++;
    }
    
    // Now loop over the reconstructed tracks so we can do some stuff!!
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

      // Only look at primary muons
      if(abs(evt->trkpdgtruth_pandoraTrack[iTrk][bestPlane]) != 13) continue;
      nMu++;
      
      // Length cuts to make sure your muons aren't fragments (3m)
      if(!evtProc.SelectTrack(evt,iTrk)) continue;
      nLongMu++;

      // Get the track length and geometry
      float length = evt->trklen_pandoraTrack[iTrk];
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
      for(const Plane &pl : allPlanes){
        if(planeN > allPlanes.size()){
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
     
      // Now look only at the external planes to determine if the track 
      // was likely to have left each plane
      for(const Plane &pl : extPlanes){
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
      
      // Now fill the number of planes crossed histograms
      h_n_crossed->Fill(nPlanesCrossed);
      h_n_ext_crossed->Fill(nExtCrossed);

      // Now count the planes crossed for the stats table
      if(nExtCrossed == 1){
        stopping++;
      }
      if(nExtCrossed >= 2){
        exiting++;
      }

      // Now you can do some analysis on your stopping muons
      // First, make sure you're only looking at the stopping muons
      if(nExtCrossed != 1) continue;

      // .
      // .
      // .
      
    } // iTrk
  }// Event loop
  std::cout << " --- 100 % --- |" << std::endl;

  // Sort out the TeX file
  std::vector<std::string> contents{
    "Tracks",
    "Muons",
    "Muons > 3m",
    "Stopping",
    "Exiting"
  };
  std::vector<unsigned int> rates{
    totalTracks,
    nMu,
    nLongMu,
    stopping,
    exiting
  };

  // Write the TeX file and print the numbers as a percentage of the total reconstructed tracks
  ofstream texFile;
  texFile.open(location+"sample_contents"+tag+".tex");
  WriteStatsToTeX(texFile, n, contents, rates, static_cast<double>(totalTracks), "All Tracks");

  // Now write the histograms
  // Plane crossing
  TCanvas *c1 = new TCanvas("c1","",900,900);
  SetCanvasStyle(c1, 0.12,0.06,0.06,0.12,0,0,0);

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

  l->AddEntry(h_plane_cross, "Crossed planes", "l");
  l->AddEntry(h_plane_enter, "Entry planes", "l");
  l->AddEntry(h_plane_exit,  "Exit planes", "l");
  l->Draw("same");

  c1->SaveAs((location+"/planes_crossed_entered_exited"+tag+".png").c_str());
  c1->SaveAs((location+"/planes_crossed_entered_exited"+tag+".root").c_str());
  c1->Clear();
  l->Clear();
  
  SetHistogramStyle1D(h_n_crossed,"Number of planes crossed [P]", " Number of tracks crossing P planes");
  // Get the maximum in y to define the axis range
  double maxy = std::max(h_n_crossed->GetMaximum(),h_n_ext_crossed->GetMaximum());
  h_n_crossed->Draw("hist");
  h_n_ext_crossed->Draw("hist same");
  h_n_crossed->SetLineWidth(2);
  h_n_ext_crossed->SetLineWidth(2);
  h_n_crossed->SetLineColor(kTeal-5);
  h_n_ext_crossed->SetLineColor(kViolet-5);
  h_n_crossed->SetLineStyle(2);
  h_n_ext_crossed->SetLineStyle(3);
  h_n_crossed->GetYaxis()->SetRangeUser(0,maxy*1.1);
  
  l->AddEntry(h_n_crossed,     "Total planes crossed", "l");
  l->AddEntry(h_n_ext_crossed, "External planes crossed", "l");
  l->Draw("same");
  
  c1->SaveAs((location+"/tracks_crossed_nplanes"+tag+".png").c_str());
  c1->SaveAs((location+"/tracks_crossed_nplanes"+tag+".root").c_str());
  c1->Clear();
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
