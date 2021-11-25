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
     
int depositionOverlay(const char *config){

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
  int fitFromPeak   = 0; // Whether to fit a certain number of bins around the peak, rather than a particular range
  int nBinsFromPeak = 1; // How many bins to traverse either side of the peak in the fit
  int overlayFit    = 1; // Whether to overlay the fit result on the output histogram
  std::vector< std::string > inputs, hists, labels, titles;
  std::string location  ="";
  std::string tag="";
  std::string xAxis="dE/dx";
  std::string units="MeV/cm";
  std::vector<double> minimums, maximums;

  p->getValue("Inputs",        inputs);
  p->getValue("Hists",         hists);
  p->getValue("Labels",        labels);
  p->getValue("Titles",        titles);
  p->getValue("Minimums",      minimums);
  p->getValue("Maximums",      maximums);
  p->getValue("FitFromPeak",   fitFromPeak);
  p->getValue("NBinsFromPeak", nBinsFromPeak);
  p->getValue("OverlayFit",    overlayFit);
  p->getValue("Location",      location);
  p->getValue("Tag",           tag);
  p->getValue("XAxis",         xAxis);
  p->getValue("Units",         units);

  unsigned int nHists = inputs.size();
  bool fitRange = false;

  // Find out if we should fit the full range -
  // i.e. we haven't specified fit from peak and we haven't given maximums and minumums
  if((minimums.size()+maximums.size() == 0) && !fitFromPeak) fitRange = true;

  // Check things are the same size
  if((inputs.size() + hists.size() + labels.size() + titles.size())/4. != nHists){
    std::cerr << " Error: Vector inputs must all be the same size, but they aren't." << std::endl;
    std::exit(1);
  }
  if((inputs.size() + hists.size() + labels.size() + titles.size()) == 0){
    std::cerr << " Error: Vector inputs must not be empty, but they are" << std::endl;
    std::exit(1);
  }

  // Sort out the file tag
  if(tag != "")
    tag = "_"+tag;

  //--------------------------------------------------------------------------------- ---------
  //                                    Initialise
  //--------------------------------------------------------------------------------- ---------

  // Setup TTree from input file list
  std::cout << " Reading histogram file" << std::endl;

  // Start of analysis (loop over chain and events
  std::cout << " Running analysis..." << std::endl;
  
  // Open a text file to print the outputs to
  ofstream ofile;
  ofile.open((location+"/statistics"+tag+".txt").c_str());
  
  // Read in the input histograms from the files
  // Project out the x-dimension, whatever that might be
  std::vector< TH2D* > histograms;
  std::vector< TH1D* > projections;
  for(unsigned int n = 0; n < nHists; ++n){
    // Access the 2D histogram
    TFile *f = new TFile(inputs.at(n).c_str(),"READ");
    TH2D  *h = static_cast<TH2D*>(f->Get(hists.at(n).c_str()));
    histograms.push_back(h);

    // Project out the x dimension
    std::string name = "h_proj_"+labels.at(n);
    TH1D *hProj = static_cast<TH1D*>(h->ProjectionY(name.c_str(),1,h->GetNbinsX()));

    std::cout << " Number of entries in histogram " << n << " : " << hProj->Integral() << std::endl;
    projections.push_back(hProj);
  }

  std::vector<TF1*> fits;
  std::vector<double> mpvs;

  unsigned int n = 0;
  for(TH1D* h : projections){

    // Scale the histograms by area
    h->Scale(1/static_cast<double>(h->Integral()));

    // Get the y limits
    double minX   = h->GetXaxis()->GetXmin();
    double maxX   = h->GetXaxis()->GetXmax();
    int maxbin    = h->GetMaximumBin();
    double maxloc = h->GetBinCenter(maxbin);
    std::string name(labels.at(n)); 
   
    // Setup the range to fit in this histogram
    double minR = minX;
    double maxR = maxX;

    // If we're not fitting a defined range, get the limits from the histogram
    if(!fitRange && !fitFromPeak){
      minR = minimums.at(n);
      maxR = maximums.at(n);
    }

    // If we want the fit range from the peak, get the range
    if(fitFromPeak){
      // Make sure the desired width doesn't exceed the limits
      if(maxbin-nBinsFromPeak > 1)
        minR = h->GetBinCenter(maxbin-nBinsFromPeak);
      if(maxbin+nBinsFromPeak < h->GetNbinsX())
        maxR = h->GetBinCenter(maxbin+nBinsFromPeak);
    }

    // Now fit a langaus to each histogram and extract the MPV
    TF1 *fit = new TF1(("fit_"+name).c_str(),langaufun,minR,maxR,4);
    fit->SetParNames("Width","MP","Area","GSigma");
    double norm = h->Integral("width"); // h->GetEntries() * h->GetBinWidth(1);
    double sv[4] = {0.2*maxloc, maxloc, norm, 0.2*maxloc}; // starting values for parameters: Landau scale, Landau MPV, Norm, Gauss sigma
    fit->SetParameters(sv);

    auto result = h->Fit(fit, "QSLR");
    double mpv = fit->GetMaximumX(result->Parameter(1), result->Parameter(1) + result->Parameter(3));

    ofile << name << " MPV: " << mpv << std::endl;

    mpvs.push_back(mpv);
    fits.push_back(fit);
    ++n;
  }

  TCanvas *c = new TCanvas("c","",900,900);
  SetCanvasStyle(c, 0.12,0.05,0.05,0.12,0,0,0);

  TLegend *l = new TLegend(0.516,0.687,0.998,0.922);
  l->SetTextFont(132);
  l->SetTextSize(0.025);
  l->SetBorderSize(0);
  l->SetFillStyle(0);

  double maxy = -999;
  for(unsigned int i = 0; i < projections.size(); ++i){
    TH1D *h = projections.at(i);
    TF1  *f = fits.at(i);

    SetHistogramStyle1D(h, (xAxis+" ["+units+"]").c_str(), "Rate");

    if(h->GetMaximum() > maxy)
      maxy = h->GetMaximum();

    h->SetLineColor(pal.at(i));
    f->SetLineColor(pal.at(i));
    h->SetLineWidth(2);
    f->SetLineWidth(1);
    h->SetLineStyle(1);
    f->SetLineStyle(7);

    // Set the precision of the mpv to print
    std::stringstream stream;
    stream << std::fixed << std::setprecision(2) << mpvs.at(i);
    std::string s = stream.str();

    l->AddEntry(h,(titles.at(i)+" MPV: "+s+" "+units).c_str(), "l");
    if(overlayFit)
      l->AddEntry(f,(titles.at(i)+" Fit").c_str(), "l");

    if(i == 0)
      h->Draw("hist");
    else
      h->Draw("hist same");
    if(overlayFit)
      f->Draw("hist same");
  }
  projections.at(0)->GetYaxis()->SetRangeUser(0,maxy*1.1);
  l->Draw();
  c->SaveAs((location+"/deposition_overlay"+tag+".root").c_str());
  c->SaveAs((location+"/deposition_overlay"+tag+".png").c_str());
  
  ofile.close();
  
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
