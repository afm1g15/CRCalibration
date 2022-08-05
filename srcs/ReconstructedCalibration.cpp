/************************************************************************
 *
 * Script to run the absolute energy scale calibration procedure using 
 * reconstructed quantites only
 *
 * First defined: August 2022
 * Author:        R. S. Jones
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
   "trklen_pandoraTrack",
   "trkke_pandoraTrack",
   "nPFParticles",
   "pfp_selfID",
   "pfp_isPrimary",
   "pfp_isTrack",
   "pfp_trackID",
   "pfp_numClusters",
   "pfp_clusterIDs",
   "pfp_numDaughters",
   "pfp_daughterIDs"
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

//Useful functions
TF1 *FitToHistogram(TCanvas *canvas, 
                    TH1D *hist, 
                    ofstream &out, 
                    const double minX, 
                    const double maxX, 
                    double &mpv,
                    double &mpv_error, 
                    double &stat_error);
     
int reconstructedCalibration(const char *config){

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
  int thru = 0;
 
  double deltaCut     = 17e-3; // Calculated by eye from the sliced delta studies
  double measuredTau  = 2.88; // [ms] measured by V. Pek
  double simulatedTau = 3; // [ms]

  double eDepError    = 0.451; // [ADC/cm] from 050822 delta-cut analysis
  double tauError     = 3.561; // [ADC/cm] from 050822 lifetimeBuffer analysis

  double nominalMPV   = 1.802; // dE/dx MPV [MeV/cm] for a 292.264 GeV muon in 5.3 mm thickness (both values from truth MPVs) from here: https://lar.bnl.gov/properties/
 
  int nBinsFromPeak   = -1; // How many bins to traverse either side of the peak in the fit
  int nBinsFromPeakL  = -1; // How many bins to traverse left of the peak in the fit
  int nBinsFromPeakR  = -1; // How many bins to traverse right of the peak in the fit
  
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
  p->getValue("Thru",            thru);
  p->getValue("DeltaCut",        deltaCut);
  p->getValue("MeasuredTau",     measuredTau);
  p->getValue("SimulatedTau",    simulatedTau);
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

  // Open a text file to print the outputs to
  ofstream ofile;
  ofile.open((location+"/statistics"+tag+".txt").c_str());
  ofile << " --------------------------------------------------------- " << std::endl;
  ofile << std::endl;

  if((nBinsFromPeak+nBinsFromPeakL+nBinsFromPeakR) < 0){
    std::cerr << " Error: Fitting a number of bins from the peak but haven't defined by how many bins" << std::endl;
    std::exit(1);
  }
  else if(nBinsFromPeak > 0 && (nBinsFromPeakL+nBinsFromPeakR) > 0){
    std::cerr << " Error: Fitting a number of bins from the peak but have set both the number (symmetric) and the L&R values (asymmetric)" << std::endl;
    std::exit(1);
  }
  else if(nBinsFromPeak < 0 && (nBinsFromPeakL+nBinsFromPeakR) > 0){
    if(nBinsFromPeakL < 0 || nBinsFromPeakR < 0){
      std::cerr << " Error: Setting the asymmetric binning from the peak but have only set one value (L or R)" << std::endl;
      std::exit(1);
    }
  }
  else{
    nBinsFromPeakL = nBinsFromPeak;
    nBinsFromPeakR = nBinsFromPeak;
  }
  std::cout << " Number of bins to fit from the peak L : " << nBinsFromPeakL << ", R : " << nBinsFromPeakR << std::endl;
  
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
  TH1D *h_dqdx              = new TH1D("h_dqdx","",400,80,600);
  TH1D *h_dedx              = new TH1D("h_dedx","",400,0.2,6);
  TH1D *h_corr_dqdx         = new TH1D("h_corr_dqdx","",400,80,600);
  TH1D *h_corr_dedx         = new TH1D("h_corr_dedx","",400,0.2,6);
  
  TH2D *h_dqdx_RR           = new TH2D("h_dqdx_RR","",400,0,15,400,20,1200);
  TH2D *h_dedx_RR           = new TH2D("h_dedx_RR","",400,0,15,400,20,1200);
  TH2D *h_corr_dqdx_RR      = new TH2D("h_corr_dqdx_RR","",400,0,15,400,100,1000);
  TH2D *h_corr_dedx_RR      = new TH2D("h_corr_dedx_RR","",400,0,15,400,100,1000);
  
  // Setup counters
  
  // Now loop over the events
  unsigned int nEvts = tree->GetEntries();
  unsigned int iIt = 1;

  std::cout << " |";
  for(unsigned int iEvt = 0; iEvt < nEvts; ++iEvt){
    tree->GetEntry(iEvt);
    if(!evtProc.SelectEvent(evt)) continue;
    unsigned int nTrks = evt->ntracks_pandoraTrack;
    unsigned int nPfps = evt->nPFParticles;
    
    // Print the processing rate
    double evtFrac  = iEvt/static_cast<double>(nEvts);

    if(std::abs(0.1*iIt-evtFrac) < std::numeric_limits<double>::epsilon()){
      std::cout << " --- " << evtFrac*100 << " %";
      std::cout.flush();
      iIt++;
    }

    // Start by looping over the PFParticles to access the primary and number-of-daughters information
    unsigned int nMultipleTrackAssoc = 0;
    for(unsigned int iPfp = 0; iPfp < nPfps; ++iPfp){
      // PFParticle quantities
      int PFPTrackID      = evt->pfp_trackID[iPfp];
      int PFPID           = evt->pfp_selfID[iPfp];
      bool PFPIsPrimary   = evt->pfp_isPrimary[iPfp];
      int PFPNumDaughters = evt->pfp_numDaughters[iPfp];

      // Make sure the PFParticle is primary
      if(!PFPIsPrimary) continue;

      // Now loop over the tracks to get length and charge deposition information
      unsigned int nTracksAssociated = 0;
      for(unsigned int iTrk = 0; iTrk < nTrks; ++iTrk){

        // Get the best plane
        int bestPlane = GetBestPlane(std::vector<int>(std::begin(evt->ntrkhits_pandoraTrack[iTrk]),std::end(evt->ntrkhits_pandoraTrack[iTrk])));

        // Make sure we're looking at reconstructed tracks associated to the current PFParticle
        if(evt->trkId_pandoraTrack[iTrk] != PFPTrackID) continue;
        
        // Make sure the reconstructed track is long enough
        // And make sure there is only 1 long track associated, 
        // skip if there is more than 1 and find out how often this happens
        if(!evtProc.SelectTrack(evt,iTrk)) continue;
        nTracksAssociated++;

        // Break out of the track loop if there is more than 1 track associated and don't fill any quantities
        if(nTracksAssociated > 1) {
          nMultipleTrackAssoc++;
          break;
        }

        // Get the reconstructed length and apply the delta-cut at 17e-3 [cm^-1]
        // This corresponds to a cut at ~30 [GeV] true muon energy
        float reco_length = evt->trklen_pandoraTrack[iTrk]; // [cm]
        float delta_perL  = PFPNumDaughters/reco_length;
        if(delta_perL < deltaCut) continue;

        // Now that we have applied the cuts, access the quantities and fill the relevant histograms
        // Get the number of hits
        unsigned int nHits = evt->ntrkhits_pandoraTrack[iTrk][bestPlane];
        
        // Make sure it doesn't exceed the maximum size of the array
        // Count if it does so we can see how often it happens
        if(nHits > MAX_TRACK_HITS){
          nHits = MAX_TRACK_HITS;
        }
        
        // Now access the hit-level variables of interest
        Float_t *dEdxArr = evt->trkdedx_pandoraTrack[iTrk][bestPlane];
        Float_t *RRArr   = evt->trkresrg_pandoraTrack[iTrk][bestPlane];
        Float_t *dQdxArr = evt->trkdqdx_pandoraTrack[iTrk][bestPlane];

        // Convert them to vectors
        std::vector<float> dEdx(dEdxArr, dEdxArr + nHits);
        std::vector<float> ResRg(RRArr, RRArr + nHits);
        std::vector<float> dQdx(dQdxArr, dQdxArr + nHits);

        // Now loop over hits so we can work our calo magic
        for(unsigned int iHit = 0; iHit < nHits; ++iHit){

          // General geometry of the track
          float x = evt->trkxyz_pandoraTrack[iTrk][bestPlane][iHit][0];
          float t = x * evtProc.kXtoT;
          
          // Check if x is lower or higher than the APA bounds, charge seems to accumulate there
          if(x < evtProc.APA_X_POSITIONS[0] || x > evtProc.APA_X_POSITIONS[2]) continue;

          // Lifetime correction
          int tpc  = evtProc.WhichTPC(x) + 1;
          float dx = ( -1 + 2*(tpc%2) )*(x - evtProc.APA_X_POSITIONS[tpc/2]);
          float dt = dx*evtProc.kXtoT;
          float corr  = TMath::Exp(-dt/measuredTau);
          float eCorr = TMath::Exp(-dt/measuredTau) / TMath::Exp(-dt/simulatedTau); // Correct for the already-corrected energy

          // New values
          float dEdxVal   = dEdx.at(iHit);
          float dQdxVal   = dQdx.at(iHit);
          float RRVal     = ResRg.at(iHit)/100.; // [m]
          float dQdxCorr  = dQdxVal/corr;
          float dEdxCorr  = dEdxVal/eCorr;
          float dEdQVal   = dEdxVal/dQdxCorr;
          
          // 1D
          h_dedx->Fill(dEdxVal);
          h_dqdx->Fill(dQdxVal);
          h_corr_dedx->Fill(dEdxCorr);
          h_corr_dqdx->Fill(dQdxCorr);
          
          // 2D
          h_dqdx_RR->Fill(RRVal,dQdxVal);
          h_dedx_RR->Fill(RRVal,dEdxVal);
          h_corr_dqdx_RR->Fill(RRVal,dQdxCorr);
          h_corr_dedx_RR->Fill(RRVal,dEdxCorr);

        } // iHit
      } // iTrk
    }// iPfp
  }// Event loop
  std::cout << " --- 100 % --- |" << std::endl;

  // Setup things for the analysis 
  TFile *f = new TFile((location+"/calibration_histograms"+tag+".root").c_str(), "RECREATE");
  TCanvas *c = new TCanvas("c","",900,900);
  SetCanvasStyle(c, 0.12,0.08,0.05,0.12,0,0,0);
  f->cd();

  // Now fit a langaus function to the corrected dQ/dx distribution and extract the MPV
  // Get the min and max range values for the fit
  double minR = h_corr_dqdx->GetXaxis()->GetXmin();
  double maxR = h_corr_dqdx->GetXaxis()->GetXmax();
  int maxbin  = h_corr_dqdx->GetMaximumBin();
  if(maxbin-nBinsFromPeak > 1)
    minR = h_corr_dqdx->GetBinCenter(maxbin-nBinsFromPeakL);
  if(maxbin+nBinsFromPeak < h_corr_dqdx->GetNbinsX())
    maxR = h_corr_dqdx->GetBinCenter(maxbin+nBinsFromPeakR);

  // Define the TF1
  double MPV = 0.;
  double mpvError = 0.;
  double statError = 0.;
  TF1 *fitFunc = FitToHistogram(c, h_corr_dqdx, ofile, minR, maxR, MPV, mpvError, statError);
  double totalError = sqrt(pow(mpvError,2)+pow(eDepError,2)+pow(tauError,2)+pow(statError,2)); // [ADC/cm]

  if(fitFunc == nullptr)
    std::cout << " Error: The returned function is null" << std::endl;

  // Now write the results
  ofile << " --------------------------------------------------------- " << std::endl;
  ofile << std::endl;
  ofile << " MPV  =  " << MPV << " +/-" << mpvError << " (syst,MPV) +/- ";
  ofile << eDepError << " (syst,Edep) +/- " << tauError << " (syst,Lifetime) +/- ";
  ofile << statError << " (stat) [ADC/cm]" << std::endl;
  ofile << std::endl;
  ofile << " MPV  =  " << MPV << " +/- " << totalError << " (stat+syst) [ADC/cm] (" << (totalError/MPV)*100. << "%)" << std::endl;
  ofile << std::endl;
  ofile << " --------------------------------------------------------- " << std::endl;
  ofile << std::endl;

  // Scale factor and uncertainty
  double fracError   = totalError/MPV;
  double scaleFactor = nominalMPV/MPV; // [MeV/ADC]
  double scaleError  = scaleFactor*fracError; // [MeV/ADC]

  ofile << " Nominal dE/dx MPV = " << nominalMPV << " [MeV/cm] " << std::endl;
  ofile << " Scale factor      =  " << scaleFactor << " +/- " << scaleError << " (stat+syst) [MeV/ADC] (" << (scaleError/scaleFactor)*100. << "%)" << std::endl;
  ofile << std::endl;
  ofile << " --------------------------------------------------------- " << std::endl;
  ofile << std::endl;
  
  // Draw the result
  // Rename the canvas
  c->SetName("c_hist");

  // Draw and save
  SetHistogramStyle1D(h_corr_dqdx, "Reconstructed dQ/dx [ADC/cm]", "Rate");
  h_corr_dqdx->SetLineWidth(2);
  h_corr_dqdx->SetLineColor(pal.at(0));
  h_corr_dqdx->SetLineStyle(1);
  h_corr_dqdx->Draw("hist");
  h_corr_dqdx->Write();
  c->Write();
  c->SaveAs((location+"/corr_dQdx"+tag+".png").c_str());
  c->SaveAs((location+"/corr_dQdx"+tag+".root").c_str());
  c->Clear();

  c->SetName("c_fit_hist");
  fitFunc->SetLineWidth(3);
  fitFunc->SetLineColor(pal.at(1));
  fitFunc->SetLineStyle(7);
  h_corr_dqdx->Draw("hist");
  fitFunc->Draw("same");
  FormatStats(h_corr_dqdx,1110,10001);
    
  // Add on the calculated information
  // Calculate the spacing in the stat box for a single entry
  double margin  = 0.05;
  double spacing = ((0.93-(0.58+margin))/9.)*2.;

  // Now construct the strings with a chosen precision
  std::stringstream ss, ss_err;
  ss << std::fixed << std::setprecision(2) << MPV;
  ss_err << std::fixed << std::setprecision(2) << mpvError;
  std::string str = ss.str();
  std::string str_err = ss_err.str();
  FormatLatexNDC(0.55+(margin/2.),0.58-spacing,("MPV  = "+str+" #pm "+str_err).c_str(),0.035,11);
  ss.str(std::string());
  ss_err.str(std::string());

  c->Write();
  c->SaveAs((location+"/fit_corr_dQdx"+tag+".png").c_str());
  c->SaveAs((location+"/fit_corr_dQdx"+tag+".root").c_str());
  c->Clear();

  // Now convert the dQ/dx to dE/dx and draw
  c->SetName("h_converted");
  double minBin = h_corr_dqdx->GetXaxis()->GetXmin()*scaleFactor;
  double maxBin = h_corr_dqdx->GetXaxis()->GetXmax()*scaleFactor;

  // If we simply want to scale the y axis by the constant from the fit, do that
  TH1D *h_conv = new TH1D("h_dedx_from_dqdx","",h_corr_dqdx->GetNbinsX(),minBin,maxBin);
  SetHistogramStyle1D(h_conv,"dE/dx [MeV/cm]", "Rate");

  // Loop over the y axis and scale the bins
  for(int nX = 1; nX <= h_conv->GetNbinsX(); ++nX){
    double chargeContent = h_corr_dqdx->GetBinContent(nX);
    h_conv->SetBinContent(nX,chargeContent);
  } // NBinsX
  h_conv->Draw("hist");
  c->SaveAs((location+"/converted_hist_dEdx"+tag+".root").c_str());
  c->SaveAs((location+"/converted_hist_dEdx"+tag+".png").c_str());
  c->Write();
  c->Clear();

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

// Function definitions
TF1 *FitToHistogram(TCanvas *canvas, 
                    TH1D *hist, 
                    ofstream &out, 
                    const double minX, 
                    const double maxX, 
                    double &mpv,
                    double &mpv_error, 
                    double &stat_error){
  
  TF1 *fit;
  // Get the FWHM and maximum values from langaupro
  int maxbin    = hist->GetMaximumBin();
  double maxloc = hist->GetBinCenter(maxbin);
  double proMPV = 0.; 
  double norm   = hist->GetEntries() * hist->GetBinWidth(1); //hist->Integral("width"); //
  double sv[4]  = {0.2*maxloc, maxloc, norm, 0.2*maxloc}; // starting values for parameters: Landau scale, Landau MPV, Norm, Gauss sigma
  fit = new TF1("fit",langaufun,minX,maxX,4);
  fit->SetParNames("Width","MP","Area","GSigma");
  fit->SetParameters(sv);

  // Now get the results
  auto result = hist->Fit(fit, "LEQSMR", "");
  mpv  = fit->GetMaximumX(result->Parameter(1), result->Parameter(1) + result->Parameter(3));
    
  for(unsigned int p = 0; p < result->NPar(); ++p){
    out << " " << result->ParName(p) << " : " << result->Parameter(p) << " +/- " << result->Error(p) << std::endl;
  }
  out << " MPV: " << mpv << std::endl;
  out << std::endl;

  // Determine the functional form of MPV in terms of MP and Gsigma for the uncertainty calculation
  // MPV = MP + psi*Gsigma
  // Setup temporary canvas for the psi value and draw
  canvas->SetName("c_psi");
  double psi = 1.;
  GetCoefficient(mpv, result->Parameter(1), result->Parameter(3), psi);
  
  TH1D *h_psi = new TH1D("h_psi","",200,0.3,0.5);
  h_psi->Fill(psi);
  SetHistogramStyle1D(h_psi,"G#sigma scaling factor","Rate");
  h_psi->Draw("hist");
  
  canvas->Write();
  canvas->Clear();
  
  // The get the uncertainty on MPV based on this formula, using the uncertainties on mp and Gsigma
  mpv_error = 0.;
  GetMPVUncertainty(mpv, result->Parameter(1), result->Error(1), result->Parameter(3), result->Error(3), psi, mpv_error);
  stat_error = mpv*(sqrt(hist->GetEntries())/static_cast<double>(hist->GetEntries()));

  out << " MPV  = MP + Psi*gsig " << std::endl; 
  out << " MP   =  " << result->Parameter(1) << " +/- " << result->Error(1) << std::endl;
  out << " GSig =  " << result->Parameter(3) << " +/- " << result->Error(3) << std::endl;
  out << " Psi  =  " << psi << std::endl << std::endl;
  out << " --------------------------------------------------------- " << std::endl;
  out << std::endl;
  out << " Uncertainties on the MPV: " << std::endl;
  out << " Total error on fit parameters: " << mpv_error << std::endl;
  out << " Error using width: " << result->Parameter(0)/static_cast<double>(sqrt(hist->GetNbinsX())) << std::endl;
  out << " Error using gsigma: " << result->Parameter(3)/static_cast<double>(sqrt(hist->GetNbinsX())) << std::endl;
  out << std::endl;
  out << " --------------------------------------------------------- " << std::endl;
  out << std::endl;

  return fit;

} // FitHistogram
