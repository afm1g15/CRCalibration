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
     
int dEdxComparison(const char *config){

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
  int yCut = 0;
  int thru = 0;
  int trimmed = 0; // Whether to apply angle and hit cuts to the sample
  int tunedScale = 0; // Whether or not to use the tuned value of Cscale from the calibration studies
  double chargeScale = 300; // Default, update from fits in the configuration
  double energyScale = 1.8002; // dE/dx MPV for a 292.264 GeV muon in 5.3 mm thickness (both values from truth MPVs) from here: https://lar.bnl.gov/properties/
  std::string input_list = "";
  std::string location="";
  std::string tag="";
  std::vector<double> minx_fid, miny_fid, minz_fid;
  std::vector<double> maxx_fid, maxy_fid, maxz_fid;
  std::vector<double> minx_av, miny_av, minz_av;
  std::vector<double> maxx_av, maxy_av, maxz_av;

  p->getValue("InputList",   input_list);
  p->getValue("Location",    location);
  p->getValue("Tag",         tag);
  p->getValue("NFiles",      n);
  p->getValue("YCut",        yCut);
  p->getValue("Thru",        thru);
  p->getValue("Trimmed",     trimmed);
  p->getValue("TunedScale",  tunedScale);
  p->getValue("ChargeScale", chargeScale);
  p->getValue("EnergyScale", energyScale);
  p->getValue("MinXFid",     minx_fid);
  p->getValue("MinYFid",     miny_fid);
  p->getValue("MinZFid",     minz_fid);
  p->getValue("MaxXFid",     maxx_fid);
  p->getValue("MaxYFid",     maxy_fid);
  p->getValue("MaxZFid",     maxz_fid);
  p->getValue("MinXAV",      minx_av);
  p->getValue("MinYAV",      miny_av);
  p->getValue("MinZAV",      minz_av);
  p->getValue("MaxXAV",      maxx_av);
  p->getValue("MaxYAV",      maxy_av);
  p->getValue("MaxZAV",      maxz_av);

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
  if(trimmed)
    tag = tag+"_trimmed";

  std::cout << "File tag: " << tag << std::endl;

  // Calculate the energy/charge scale factor
  float energyChargeScale = energyScale/static_cast<double>(chargeScale);
  std::cout << " Energy scale: " << energyScale << ", chargeScale: " << chargeScale << ", energyChargeScale: " << energyChargeScale << std::endl;
  
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
  TH2D *h_reco_dedx_dqdx              = new TH2D("h_reco_dedx_dqdx","",200,0,6,200,0,1012);
  TH2D *h_reco_dedx_conv_dedx         = new TH2D("h_reco_dedx_conv_dedx","",200,0,6,200,0,6);
  TH1D *h_modbox_dedx_conv_dedx       = new TH1D("h_modbox_dedx_conv_dedx","",200,0,6);
  TH1D *h_modbox_dedx_conv_dedx_tuned = new TH1D("h_modbox_dedx_conv_dedx_tuned","",200,0,6);
  TH1D *h_modbox_dedx_dqdx            = new TH1D("h_modbox_dedx_dqdx","",200,0,6);
  TH1D *h_modbox_dedx_dqdx_tuned      = new TH1D("h_modbox_dedx_dqdx_tuned","",200,0,6);
  h_modbox_dedx_conv_dedx->GetYaxis()->SetRangeUser(0,8);
  h_modbox_dedx_conv_dedx_tuned->GetYaxis()->SetRangeUser(0,8);
  h_modbox_dedx_dqdx->GetYaxis()->SetRangeUser(0,1012);
  h_modbox_dedx_dqdx_tuned->GetYaxis()->SetRangeUser(0,1012);
  
  // Setup counters

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

      // Setup list of plane labels the track has crossed
      std::vector<std::string> labelsCrossed;

      // Get the best plane
      std::vector<int> nHitsPerPlane{evt->ntrkhits_pandoraTrack[iTrk][0],
                                     evt->ntrkhits_pandoraTrack[iTrk][1],
                                     evt->ntrkhits_pandoraTrack[iTrk][2]};

      int bestPlane = GetBestPlane(nHitsPerPlane);
      int ID        = evt->trkId_pandoraTrack[iTrk];
      double length = evt->trklen_pandoraTrack[iTrk];

      // Only look at primary muons
      if(abs(evt->trkpdgtruth_pandoraTrack[iTrk][bestPlane]) != 13) continue;
      
      // Length cuts (2m)
      if(!evtProc.SelectTrack(evt,iTrk)) continue;
      
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

      // Make sure the tracks start at the top of the detector
      if(yCut && startVtx.Y() < 599.5) continue; 

      // The following studies should be conducted with through-going muons to start with
      // If the number of external planes crossed is >= 2, the track is through-going
      bool throughGoing = IsThroughGoing(length,startVtx,endVtx,extPlanes,fidExtPlanes);
      if(thru && !throughGoing) continue;

      // If we want the 'trimmed' sample, find out the appropriate parameters and apply the cut
      float cosDrift = GetCosDrift(startVtx,endVtx); 
      double hitsPerL = nHitsPerPlane.at(bestPlane)/static_cast<double>(length);
      if(trimmed){
        if(cosDrift < -0.567 || cosDrift > 0.564) continue;
        if(hitsPerL < 0.8) continue;
      }

      // Now fill dQ/dx and dE/dx and hit histograms for each of the three wire planes
      // Somehow flag the best wire plane histogram
      for(int iPlane = 0; iPlane < 3; ++iPlane){

        // Use only best plane for now
        if(iPlane != bestPlane) continue;

        unsigned int nHits = evt->ntrkhits_pandoraTrack[iTrk][iPlane];

        // Get the associated true energy of the muon
        // Try to get the true energy
        // Get the list iterator from matching ID's
        float eng = -1.;
        for(int iG4 = 0; iG4 < nGeant; ++iG4){
          int trueID = evt->TrackId[iG4];

          if(evt->trkidtruth_pandoraTrack[iTrk][bestPlane] == trueID){
            eng = evt->Eng[iG4];
            break;
          }
        }
        if(eng < 0){
          std::cout << " Warning: Energy is below zero, skipping track with energy: " << eng << std::endl;
          continue;
        }
     
        // Make sure it doesn't exceed the maximum size of the array
        // Count if it does so we can see how often it happens
        if(nHits > MAX_TRACK_HITS){
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
          int tpc  = evtProc.WhichTPC(x) + 1;
          float dx = ( -1 + 2*(tpc%2) )*(x - evtProc.APA_X_POSITIONS[tpc/2]);
          float dt = dx*evtProc.kXtoT;
          float corr  = TMath::Exp(-dt/2.88);
          float eCorr = TMath::Exp(-dt/2.88) / TMath::Exp(-dt/3.); // Correct for the already-corrected energy

          // New values
          float dEdxVal   = dEdx.at(iHit);
          float dQdxVal   = dQdx.at(iHit);
          float dQdxCorr  = dQdxVal/corr;
          float dEdxCorr  = dEdxVal/eCorr;

          // Fill dE/dx vs dQ/dx
          h_reco_dedx_dqdx->Fill(dEdxCorr,dQdxCorr);

          // Now get the corresponding ModifiedBox value of dEdx for the current, lifetime-corrected dQdx
          float dEdxModBox      = ModBoxCorrection(dQdxCorr);
          float dEdxModBoxTuned = ModBoxCorrection(dQdxCorr,tunedScale);

          // Now define the energy-charge conversion and plot the 2D comparison
          // scale = energy/charge and should be applied to the reconstructed charge
          float dEdxScaled = dQdxCorr*energyChargeScale;
          h_reco_dedx_conv_dedx->Fill(dEdxCorr,dEdxScaled);

          // Find the relevant dEdx bin
          int dEdxBin      = h_modbox_dedx_conv_dedx->FindBin(dEdxModBox);
          int dEdxBinTuned = h_modbox_dedx_conv_dedx->FindBin(dEdxModBoxTuned);

          h_modbox_dedx_conv_dedx->SetBinContent(dEdxBin,dEdxScaled);
          h_modbox_dedx_conv_dedx_tuned->SetBinContent(dEdxBinTuned,dEdxScaled);
          h_modbox_dedx_dqdx->SetBinContent(dEdxBin,dQdxCorr);
          h_modbox_dedx_dqdx_tuned->SetBinContent(dEdxBinTuned,dQdxCorr);
          
        } // Hits
      } // Planes
    } // Tracks
  }// Event loop
  std::cout << " --- 100 % --- |" << std::endl;

  TCanvas *c1 = new TCanvas("c1","",1000,800);
  SetCanvasStyle(c1, 0.1,0.12,0.05,0.12,0,0,0);

  SetHistogramStyle2D(h_reco_dedx_conv_dedx,"dEdx (Reconstructed) [MeV/cm]", " dE/dx (Converted) [MeV/cm]",false);
  h_reco_dedx_conv_dedx->GetZaxis()->SetLabelSize(0.03);
  h_reco_dedx_conv_dedx->GetZaxis()->SetLabelFont(132);
  h_reco_dedx_conv_dedx->Draw("colz");
  h_modbox_dedx_conv_dedx->SetLineColor(kCyan-5);
  h_modbox_dedx_conv_dedx_tuned->SetLineColor(kPink+4);
  h_modbox_dedx_conv_dedx->SetLineStyle(1);
  h_modbox_dedx_conv_dedx_tuned->SetLineStyle(1);
  h_modbox_dedx_conv_dedx->SetLineWidth(3);
  h_modbox_dedx_conv_dedx_tuned->SetLineWidth(3);
  h_modbox_dedx_conv_dedx->Draw("hist same l");
  h_modbox_dedx_conv_dedx_tuned->Draw("hist same l");
  h_modbox_dedx_conv_dedx->Smooth();
  h_modbox_dedx_conv_dedx_tuned->Smooth();

  // Overlay a line at x == y
  TLine *l = new TLine(0,0,6,6);
  l->SetLineWidth(3);
  l->SetLineColor(kViolet+5);
  l->SetLineStyle(2);
  l->Draw();

  TLegend *leg = new TLegend(0.15,0.75,0.60,0.90);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(132);
  leg->SetTextSize(0.036);

  // Double to string with precision
  // Create an output string stream
  std::ostringstream ssCScale;
  // Set Fixed -Point Notation
  ssCScale << std::fixed;
  // Set precision to 2 digits
  ssCScale << std::setprecision(2);
  //Add double to stream
  ssCScale << energyChargeScale*1000;
  // Get string from output string stream
  std::string eCScale = ssCScale.str();
  

  leg->AddEntry(h_modbox_dedx_conv_dedx,"ModBox model, C_{Scale} = 5 #times 10^{-3}","l");
  leg->AddEntry(h_modbox_dedx_conv_dedx_tuned,("ModBox model, C_{Scale} = "+eCScale+" #times 10^{-3}").c_str(),"l");
  leg->Draw("same");

  h_reco_dedx_conv_dedx->SaveAs((location+"/hist_reco_vs_conv_dEdx"+tag+".root").c_str());
  c1->SaveAs((location+"/reco_vs_conv_dEdx"+tag+".png").c_str());
  c1->SaveAs((location+"/reco_vs_conv_dEdx"+tag+".root").c_str());
  c1->Clear();
  leg->Clear();

  SetHistogramStyle2D(h_reco_dedx_dqdx,"dEdx [MeV/cm]", " dQ/dx [ADC/cm]",false);
  h_reco_dedx_dqdx->GetZaxis()->SetLabelSize(0.03);
  h_reco_dedx_dqdx->GetZaxis()->SetLabelFont(132);
  h_reco_dedx_dqdx->Draw("colz");
  h_modbox_dedx_dqdx->SetLineColor(kCyan-5);
  h_modbox_dedx_dqdx_tuned->SetLineColor(kPink+4);
  h_modbox_dedx_dqdx->SetLineStyle(1);
  h_modbox_dedx_dqdx_tuned->SetLineStyle(1);
  h_modbox_dedx_dqdx->SetLineWidth(3);
  h_modbox_dedx_dqdx_tuned->SetLineWidth(3);
  h_modbox_dedx_dqdx->Draw("hist same l");
  h_modbox_dedx_dqdx_tuned->Draw("hist same l");
  h_modbox_dedx_dqdx->Smooth();
  h_modbox_dedx_dqdx_tuned->Smooth();
  l->SetY2(1e3);
  l->Draw();

  leg->AddEntry(h_modbox_dedx_dqdx,"ModBox model, C_{Scale} = 5 #times 10^{-3}","l");
  leg->AddEntry(h_modbox_dedx_dqdx_tuned,("ModBox model, C_{Scale} = "+eCScale+" #times 10^{-3}").c_str(),"l");
  leg->Draw("same");
  
  h_reco_dedx_dqdx->SaveAs((location+"/hist_reco_dEdx_vs_dQdx"+tag+".root").c_str());
  c1->SaveAs((location+"/reco_dEdx_vs_dQdx"+tag+".png").c_str());
  c1->SaveAs((location+"/reco_dEdx_vs_dQdx"+tag+".root").c_str());
  c1->Clear();
  leg->Clear();

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
