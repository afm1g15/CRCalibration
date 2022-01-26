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
     
int langausFit(const char *config){

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
  int fitRange           = 1; // Whether or not to fit the full range of the function
  double fitMin          = 99999.;
  double fitMax          = -99999.;
  int fitExp             = 1; // Whether to fit the MPV's to an exponential
  int fitPol2            = 0; // Whether to fit the MPV's to a second order polynomial
  int fitFromPeak        = 0; // Whether to define the fit from the peak outwards
  int nBinsFromPeak      = 1; // How many bins to traverse either side of the peak in the fit
  int logSpace           = 0; // Whether to spread the bins in log space
  int logX               = 0; // Whether the x dimension is in log space 
  double buffer          = 0;
  int lineAtMPV          = 1; // Whether to draw a line at the MPV or the max, 0 == max, 1 == mpv
  std::string fitFunc    = "langaus";
  std::string inputFile  = "";
  std::string inputHist  = "";
  std::string location   = "";
  std::string tag        = "";

  p->getValue("FitRange",        fitRange);
  p->getValue("FitMin",          fitMin);
  p->getValue("FitMax",          fitMax);
  p->getValue("FitExp",          fitExp);
  p->getValue("FitPol2",         fitPol2);
  p->getValue("FitFromPeak",     fitFromPeak);
  p->getValue("NBinsFromPeak",   nBinsFromPeak);
  p->getValue("LogSpace",        logSpace);
  p->getValue("LogX",            logX);
  p->getValue("Buffer",          buffer);
  p->getValue("LineAtMPV",       lineAtMPV);
  p->getValue("FitFunction",     fitFunc);
  p->getValue("InputFile",       inputFile);
  p->getValue("InputHist",       inputHist);
  p->getValue("Location",        location);
  p->getValue("Tag",             tag);

  // Read in the input histogram from the file
  TFile *histFile = new TFile(inputFile.c_str(),"READ");
  TH1D *h = static_cast<TH1D*>(histFile->Get(inputHist.c_str()));

  // Sort out the file tag
  if(tag != "")
    tag = "_"+tag;

  // Open a text file to print the outputs to
  ofstream ofile;
  ofile.open((location+"/statistics"+tag+".txt").c_str());

  //--------------------------------------------------------------------------------- ---------
  //                                    Initialise
  //--------------------------------------------------------------------------------- ---------

  // Setup TTree from input file list
  std::cout << " Reading histogram file" << std::endl;

  // Start of analysis (loop over chain and events
  std::cout << " Running analysis..." << std::endl;

  // Now define the relevant histograms and fill them
  // Get the y range from the 2D histogram first
  double minX = h->GetXaxis()->GetXmin();
  double maxX = h->GetXaxis()->GetXmax();
  double minY = h->GetMinimum();
  double maxY = h->GetMaximum();

  // Get the x-axis label for writing later
  std::string xLabel = static_cast<std::string>(h->GetXaxis()->GetTitle());

  ofile << " X: (" << minX << ", " << maxX << ")" << std::endl;
  ofile << " Y: (" << minY << ", " << maxY << ")" << std::endl;

  // Sort out the input histogram
  SetHistogramStyle1D(h, xLabel.c_str(), "Rate");
  
  // Now draw and save
  TFile *f = new TFile((location+"/langaus_fit"+tag+".root").c_str(), "RECREATE");
  TCanvas *c = new TCanvas("c","",900,900);
  SetCanvasStyle(c, 0.12,0.08,0.05,0.12,0,0,0);
  f->cd();

  // Also find the minimum and maximum MPV, the average MPV,
  // and get the min-max difference wrt the average
  double minMPV       = 999.;
  double maxMPV       = -999.;
  double avgMPV       = 0.;
  double totFracSigma = 0.;
  int maxbin          = h->GetMaximumBin();
  double maxloc       = h->GetBinCenter(maxbin);

  // Now define the TF1
  TF1 *fit;
  double minR = minX;
  double maxR = maxX;

  // If we're not fitting a defined range, get the limits from the histogram
  if(!fitRange && !fitFromPeak){
    // Check that the limits have been set, otherwise use the histogram values
    if(fitMin < 99999)
      minR = fitMin;
    if(fitMax > -99999.)
      maxR = fitMax;
  }

  // If we want the fit range from the peak, get the range
  if(fitFromPeak){
    // Make sure the desired width doesn't exceed the limits
    if(maxbin-nBinsFromPeak > 1)
      minR = h->GetBinCenter(maxbin-nBinsFromPeak);
    if(maxbin+nBinsFromPeak < h->GetNbinsX())
      maxR = h->GetBinCenter(maxbin+nBinsFromPeak);
  }

  if(fitFunc == "langaus"){
    // Set some approximate start parameters
    fit = new TF1("fit",langaufun,minR,maxR,4);
    fit->SetParNames("Width","MP","Area","GSigma");
    double norm = h->GetEntries() * h->GetBinWidth(1); // h->Integral("width"); //
    double sv[4] = {0.2*maxloc, maxloc, norm, 0.2*maxloc}; // starting values for parameters: Landau scale, Landau MPV, Norm, Gauss sigma
    fit->SetParameters(sv);

  }
  else{
    // Set some approximate start parameters
    fit = new TF1("fit",fitFunc.c_str(),minR,maxR);
  }

  // Now get the results
  auto result = h->Fit(fit, "QSMR", "");
  double mpv = result->Parameter(1);

  // Find the x value corresponding to the maximum of the function in the interval MPV, MPV+GSig
  if(fitFunc == "langaus")
    mpv = fit->GetMaximumX(result->Parameter(1), result->Parameter(1) + result->Parameter(3));

  ofile << " ----------------------------------------" << std::endl;
  for(unsigned int p = 0; p < result->NPar(); ++p){
    ofile << " " << result->ParName(p) << " : " << result->Parameter(p) << std::endl;
  }
  ofile << " ----------------------------------------" << std::endl;
  ofile << " Peak pitch: " << maxloc << std::endl;
  ofile << " ----------------------------------------" << std::endl;

  // Get the error on the mpv from both the fit scaled with the error due to the statistics
  double tot_error = result->Parameter(3)/static_cast<double>(sqrt(h->GetNbinsX()));

  if(mpv > maxMPV)
    maxMPV = mpv;
  if(mpv < minMPV)
    minMPV = mpv;
  totFracSigma += pow(tot_error/static_cast<double>(mpv),2);

  // Rename the canvas
  c->SetName("c_fit");

  // Draw and save
  FormatStats(h,10,10);
  fit->SetLineWidth(1);
  fit->SetLineColor(pal.at(0));
  fit->SetLineStyle(7);
  h->Draw("hist");
  fit->Draw("same");
  c->Write();
  c->SaveAs((location+"/fit"+tag+".png").c_str());
  c->SaveAs((location+"/fit"+tag+".root").c_str());
  c->Clear();
  
  // Rename the canvas
  c->SetName("h_orig");
 
  h->SetStats(0); 
  h->Draw("hist");
  TLine *l = new TLine(maxloc, minY, maxloc, maxY*1.05);
  if(lineAtMPV){
    l->SetX1(mpv);
    l->SetX2(mpv);
  }
  l->SetLineColor(pal.at(0));
  l->SetLineWidth(3);
  l->SetLineStyle(2);
  l->Draw();
  c->Write();
  c->SaveAs((location+"/line"+tag+".png").c_str());
  c->SaveAs((location+"/line"+tag+".root").c_str());
  c->Clear();

  std::cout << " Writing all slices to file: " << (location+"/fit_histograms"+tag+".root").c_str() << std::endl;

  f->Write();
  ofile.close();

  // End of script
  std::cout << " ...finished analysis" << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;
  time_t rawtime_end;
  GetTime(rawtime_end);
  std::cout << "-----------------------------------------------------------" << std::endl;
  GetTotalTime(rawtime, rawtime_end);
  std::cout << "-----------------------------------------------------------" << std::endl;
  
  f->Close();
 
  return 0;
}
