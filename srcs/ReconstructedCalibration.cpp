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
                    double &mp,
                    double &mpv,
                    double &mpv_error, 
                    double &stat_error);

void GetThrownUncertainty(const std::vector<double> &dqdx, 
                          TCanvas *c_throw,
                          ofstream &out, 
                          const double &scale, 
                          const double &err, 
                          const int &nThrows,
                          const int &nBins,
                          const double &min,
                          const double &max,
                          const int &nBinsFromPeak,
                          const std::string &loc,
                          const std::string &tag,
                          std::vector<double> &binErrors);
     
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

  int nDivDays        = 1; // Number of days to divide the statistical uncertainty histogram into
  int nBins           = 400; // Number of bins for the dQ/dx and dE/dx histogram

  int nToys           = 1000; // Number of times to throw the scale-factor within it's 1sigma uncertainty to calculate dE/dx 1sigma

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
  p->getValue("NDivDays",        nDivDays);
  p->getValue("NBins",           nBins);
  p->getValue("NToys",           nToys);
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

  // If we've set L & R set the average to be the general nBins
  if(nBinsFromPeak < 0){
    nBinsFromPeak = (nBinsFromPeakL+nBinsFromPeakR)/2.;
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
  TH1D *h_dqdx              = new TH1D("h_dqdx","",nBins,80,600);
  TH1D *h_dedx              = new TH1D("h_dedx","",nBins,0.2,6);
  TH1D *h_corr_dqdx         = new TH1D("h_corr_dqdx","",nBins,80,600);
  TH1D *h_corr_dedx         = new TH1D("h_corr_dedx","",nBins,0.48,3.6);
  
  TH2D *h_dqdx_RR           = new TH2D("h_dqdx_RR","",nBins,0,15,nBins,20,1200);
  TH2D *h_dedx_RR           = new TH2D("h_dedx_RR","",nBins,0,15,nBins,20,1200);
  TH2D *h_corr_dqdx_RR      = new TH2D("h_corr_dqdx_RR","",nBins,0,15,nBins,100,1000);
  TH2D *h_corr_dedx_RR      = new TH2D("h_corr_dedx_RR","",nBins,0,15,nBins,100,1000);

  // Vectors of relevant quantities for the scaling part of the analysis 
  std::vector<double> v_dqdx_corr, v_dedx_corr; 

  // Setup counters
  
  // Total number of events
  unsigned int nEvts = tree->GetEntries();
  
  // Sussing out the statistical uncertainty with time
  // First, estimate the total number of days in the sample
  // Calculate the approximate number of days from the number of files
  // nFiles * 500 events per file / 14114 events per day (0.16356 s^{-1})
  // Setup the histogram to calculate a statistical uncertainty entry every 10 days
  unsigned int nDays = nEvts/14114.;
  double dNDays      = nDays/static_cast<double>(nDivDays);
  int nNDays         = std::floor(dNDays);
  int nPerNDays      = 14114*nDivDays;
  std::cout << " Total number of events: " << nEvts
            << ", corresponding to " << nDays 
            << " days and " << nNDays 
            << " lots of " << nDivDays
            << " day(s) at " << nPerNDays 
            << " events per " << nDivDays << " day(s) " << std::endl;

  TH1D *h_stat_err_time = new TH1D("h_stat_err_time","",nNDays,0,nNDays);

  // Now loop over the events
  unsigned int iIt    = 1;
  unsigned int itDays = 1;
  double frac_stat    = 0.;
  std::cout << " |";
  for(unsigned int iEvt = 0; iEvt < nEvts; ++iEvt){
    tree->GetEntry(iEvt);
    if(!evtProc.SelectEvent(evt)) continue;
    unsigned int nTrks = evt->ntracks_pandoraTrack;
    unsigned int nPfps = evt->nPFParticles;
    
    // Print the processing rate
    double evtFrac  = iEvt/static_cast<double>(nEvts);
    if(std::abs(0.1*iIt-evtFrac) < 1e-5){
      std::cout << " --- " << std::setprecision(2) << evtFrac*100 << " %";
      std::cout.flush();
      iIt++;
    }

    // Fill the fractional uncertainty histogram
    int currEv = iEvt-(itDays*nPerNDays);
    if(std::abs(currEv) == 0){
      // Fractional uncertainty scaled for the number of bins in the dQ/dx distribution
      // sqrt(N/bins)/(N/bins) = sqrt(bins/N)
      frac_stat = sqrt((nBins)/static_cast<double>(iEvt));
      //std::cout << " Event: " << iEvt << ", " << itDays << " lots of " << nDivDays << " days have passed and the statistical uncertainty is: " << frac_stat*100. << " % " << std::endl;
      h_stat_err_time->Fill(itDays-1,frac_stat*100.);
      itDays++;
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
          h_dqdx->Fill(dQdxVal);
          h_dedx->Fill(dEdxVal);
          h_corr_dqdx->Fill(dQdxCorr);
          h_corr_dedx->Fill(dEdxCorr);
          
          // 2D
          h_dqdx_RR->Fill(RRVal,dQdxVal);
          h_dedx_RR->Fill(RRVal,dEdxVal);
          h_corr_dqdx_RR->Fill(RRVal,dQdxCorr);
          h_corr_dedx_RR->Fill(RRVal,dEdxCorr);

          // Vectors
          v_dqdx_corr.push_back(dQdxCorr);
          v_dedx_corr.push_back(dEdxCorr);

        } // iHit
      } // iTrk
    }// iPfp
  }// Event loop
  std::cout << " --- 100 % --- |" << std::endl;

  // Setup things for the analysis 
  TFile *f = new TFile((location+"/calibration_histograms"+tag+".root").c_str(), "RECREATE");
  TCanvas *c = new TCanvas("c","",900,900);
  SetCanvasStyle(c, 0.12,0.06,0.06,0.12,0,0,0);
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
  double MP = 0.;
  double MPV = 0.;
  double mpvError = 0.;
  double statError = 0.;
  TF1 *fitFunc = FitToHistogram(c, h_corr_dqdx, ofile, minR, maxR, MP, MPV, mpvError, statError);
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
  //double fracError   = totalError/MPV;
  //double scaleFactor = nominalMPV/MPV; // [MeV/ADC]
  double fracError   = totalError/MP;
  double scaleFactor = nominalMPV/MP; // [MeV/ADC]
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
  FormatStats(h_corr_dqdx,1110,101);
    
  // Add on the calculated information
  // Calculate the spacing in the stat box for a single entry
  double margin  = 0.05;
  double spacing = ((0.93-(0.58+margin))/9.)*2.;

  // Now construct the strings with a chosen precision
  std::stringstream ss, ss_err;
  ss << std::fixed << std::setprecision(2) << MPV;
  ss_err << std::fixed << std::setprecision(2) << totalError;
  std::string str = ss.str();
  std::string str_err = ss_err.str();
  FormatLatexNDC(0.58+(margin/2.),0.58-spacing,("MPV  = "+str+" #pm "+str_err).c_str(),0.035,11);
  ss.str(std::string());
  ss_err.str(std::string());

  c->Write();
  c->SaveAs((location+"/fit_corr_dQdx"+tag+".png").c_str());
  c->SaveAs((location+"/fit_corr_dQdx"+tag+".root").c_str());
  c->Clear();


  // -------------------------------------------------------------------------------------------
  //                                        APPLICATION
  // -------------------------------------------------------------------------------------------
  //
  // Now convert the dQ/dx to dE/dx and draw
  const double minBin = 0.48; 
  const double maxBin = 3.6;

  // If we simply want to scale the y axis by the constant from the fit, do that
  c->SetName("h_converted");
  TH1D *h_conv = new TH1D("h_dedx_from_dqdx","",nBins,minBin,maxBin);
  SetHistogramStyle1D(h_conv,"dE/dx [MeV/cm]", "Rate [cm/MeV]");
  SetHistogramStyle1D(h_corr_dedx,"dE/dx [MeV/cm]", "Rate [cm/MeV]");

  // Loop over the y axis and scale the bins
  // Average dE/dx
  double avgdEdx = 0.;
  for(const double &dq: v_dqdx_corr){
    h_conv->Fill(dq*scaleFactor);
    avgdEdx += dq*scaleFactor;
  } // dQdx
  
  avgdEdx /= static_cast<double>(v_dqdx_corr.size());

  ofile << " Average converted dE/dx = " << avgdEdx << " [MeV/cm] " << std::endl;

  // Now extract the dE/dx uncertainty from the scale-factor uncertainty 
  // by throwing the scale-factor NToys times within its own uncertainty
  // First, setup vectors of 1sigma errors in NBins
  std::vector<double> binErrors;
  GetThrownUncertainty(v_dqdx_corr, c, ofile, scaleFactor, scaleError, nToys, nBins, minBin, maxBin, nBinsFromPeak, location, tag, binErrors);

  for(int nX = 1; nX <= h_conv->GetNbinsX(); ++nX){
    h_conv->SetBinError(nX,binErrors.at(nX-1)*h_conv->GetBinContent(nX));
  }

  h_conv->Scale(1,"width");
  h_conv->Draw("E1 X0");
  h_conv->Draw("hist same");
  h_corr_dedx->Scale(1,"width");
  h_corr_dedx->Draw("hist same");
  h_conv->SetLineWidth(1);
  h_corr_dedx->SetLineWidth(2);
  h_conv->SetLineColor(pal.at(1));
  h_conv->SetMarkerColor(pal.at(1));
  h_corr_dedx->SetLineColor(pal.at(0));
  h_conv->SetLineStyle(1);
  h_corr_dedx->SetLineStyle(1);

  double max = std::max(h_conv->GetMaximum(),h_corr_dedx->GetMaximum());
  h_conv->GetYaxis()->SetRangeUser(0,max*1.1);

  TLegend *l = new TLegend(0.24,0.95,0.90,0.99);
  l->SetNColumns(2);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->SetTextFont(132);
  l->SetTextSize(0.028);
  l->AddEntry(h_corr_dedx,"Mod-Box dE/dx","l");
  l->AddEntry(h_conv,"Scaled dE/dx","l");

  l->Draw();

  c->SaveAs((location+"/converted_hist_dEdx"+tag+".root").c_str());
  c->SaveAs((location+"/converted_hist_dEdx"+tag+".png").c_str());
  c->Write();
  c->Clear();

  // Fractional statistical uncertainty against time
  c->SetName("frac_stat");
  SetHistogramStyle1D(h_stat_err_time,"Number of days","Statistical uncertainty/N_{Bins} [%]");
  h_stat_err_time->GetYaxis()->SetRangeUser(-0.2,20);
  h_stat_err_time->Draw("hist");
  h_stat_err_time->SetLineWidth(2);
  h_stat_err_time->SetLineColor(pal.at(0));
  h_stat_err_time->SetLineStyle(1);

  // Draw a line at 0%
  TLine *line = new TLine(0,0,nNDays,0);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  line->SetLineColor(kGray+2);
  line->Draw();

  c->SaveAs((location+"/fractional_statistical_error"+tag+".root").c_str());
  c->SaveAs((location+"/fractional_statistical_error"+tag+".png").c_str());
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
                    double &mp,
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
  mp   = result->Parameter(1);
    
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
  out << " Chi^2/NDOF = " << result->Chi2() << " / " << result->Ndf() << std::endl; 
  out << std::endl;
  out << " --------------------------------------------------------- " << std::endl;
  out << std::endl;

  return fit;

} // FitHistogram

void GetThrownUncertainty(const std::vector<double> &dqdx, 
                          TCanvas *c_throw,
                          ofstream &out, 
                          const double &scale, 
                          const double &err, 
                          const int &nThrows,
                          const int &nBins,
                          const double &min,
                          const double &max,
                          const int &nBinsFromPeak,
                          const std::string &loc,
                          const std::string &tag,
                          std::vector<double> &binErrors){

  // If we simply want to scale the y axis by the constant from the fit, do that
  TH1D *hist_conv = new TH1D("hist_conv","",nBins,min,max);
  SetHistogramStyle1D(hist_conv,"dE/dx [MeV/cm]", "Rate [cm/MeV]");

  // Loop over the y axis and scale the bins
  for(const double &dq: dqdx){
    hist_conv->Fill(dq*scale);
  } // dQdx

  hist_conv->Scale(1,"width");
  // Calculate relative number of bins 
  // Define vector of histograms for the throws
  int nSliceBins  = nThrows/15.;
  double binRatio = (nSliceBins)/static_cast<double>(nBins);
  int nFromPeak   = std::round(2*nBinsFromPeak*binRatio);
  std::vector<TH1D*> h_throws;

  c_throw->SetName("c_throws");
  // Loop over the number of throws
  for(int t = 0; t < nThrows; ++t){
    // Extract a random scale factor between +/- 1sigma of the scale factor
    TRandom3 rnd;
    rnd.SetSeed(t);
    //double randScale = rnd.Uniform(scale-err,scale+err);
    double randScale = rnd.Gaus(scale, err);
  //  std::cout << " Scale: " << std::setprecision(6) <<  scale << " error: " << std::setprecision(6) << err << " random: " << std::setprecision(6) <<  randScale << std::endl;

    // Now construct a TH1D for the current throw's dE/dx
    TH1D *h_curr_dedx = new TH1D(("h_dedx_throw_"+std::to_string(t)).c_str(),"",nBins,min,max);
  
    // Loop over the y axis and scale the bins for this throw
    for(const double &dq: dqdx){
      h_curr_dedx->Fill(dq*randScale);
    } // dQdx

    SetHistogramStyle1D(h_curr_dedx,"Thrown dE/dx [MeV/cm]", "Rate [cm/MeV]");
    h_curr_dedx->Scale(1,"width");
    if(t == 0)
      h_curr_dedx->Draw("hist");
    else
      h_curr_dedx->Draw("hist same");
    h_curr_dedx->SetLineWidth(1);
    h_curr_dedx->SetLineColor(kGray+1);
    h_curr_dedx->SetLineStyle(7);

    h_throws.push_back(h_curr_dedx);
  }
  hist_conv->Draw("hist same");
  hist_conv->SetLineWidth(3);
  hist_conv->SetLineColor(kSpring-5);
  c_throw->SaveAs((loc+"/dEdx_throws"+tag+".root").c_str());
  c_throw->SaveAs((loc+"/dEdx_throws"+tag+".png").c_str());
  c_throw->Write();
  c_throw->Clear();

  // Now loop over the nominal histogram and get each bin centre, 
  // then access the bin content for that bin for every throw
  for(int nX = 1; nX <= nBins; ++nX){
    double binCentre  = hist_conv->GetBinCenter(nX);
    double binContent = hist_conv->GetBinContent(nX);
    double minY = binContent-binContent*0.6;
    double maxY = binContent+binContent*0.6;

    // Setup 1D histogram for the current bin
    TH1D *h_bin = new TH1D(("h_rate_bin_"+std::to_string(nX)).c_str(),"",nSliceBins,minY,maxY);

    // Now get the bin content for the current bin in each throw
    for(int t = 0; t < nThrows; ++t){
      TH1D *h_curr = h_throws.at(t);
      h_bin->Fill(h_curr->GetBinContent(nX));
    }

    double minR = h_bin->GetXaxis()->GetXmin();
    double maxR = h_bin->GetXaxis()->GetXmax();
    int maxbin  = h_bin->GetMaximumBin();
    double maxContent = h_bin->GetBinCenter(maxbin);
    if(maxbin-nFromPeak > 1)
      minR = h_bin->GetBinCenter(maxbin-nFromPeak);
    if(maxbin+nFromPeak < h_bin->GetNbinsX())
      maxR = h_bin->GetBinCenter(maxbin+nFromPeak);

    // Now fit to the current TH1 to extract the 1sigma uncertainty
    TF1 *fitGaus = new TF1("fitGaus","gaus(0)",minR,maxR);
    auto gausResult = h_bin->Fit(fitGaus, "LEQSMR", "");
    for(unsigned int p = 0; p < gausResult->NPar(); ++p){
      out << " " << gausResult->ParName(p) << " : " << gausResult->Parameter(p) << " +/- " << gausResult->Error(p) << std::endl;
    }
    
    c_throw->SetName(("c_bin_"+std::to_string(nX)).c_str());
    h_bin->Draw("hist");
    fitGaus->Draw("hist same");
    fitGaus->SetLineStyle(2);
    c_throw->Write();
    c_throw->Clear();
 
    // Parameter #2 is the 1sigma uncertainty, apply this to the dE/dx distribution as a fraction of the mean, Parameter(1)
    binErrors.push_back(gausResult->Parameter(2)/binContent);
    
  } // nX

  return;

} // Get 1sigma dE/dx uncertainty
