/************************************************************************
 * 
 * A macro to overlay and compare dE/dx and dQ/dx distributions 
 * constructed in various ways
 *    - Reconstructed
 *    - Converted (using calibration procedure)
 *
 * Note (July 2022):
 *    This has not been updated recently to reflect sample and calibration
 *    changes.
 *
 * Example file list located here:
 *   /home/jones/work/cosmics/LArSoft-v08_50_00/work/files/v09_41_00_02_files.list
 *
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
  // Create object of the classConfigReader
  // Parse the configuration file
  // Dump map on the console after parsing it
  ConfigReader* p = ConfigReader::getInstance();
  p->parseFile(config);
  std::cout << " Variables from configuration file: " << std::endl;
  p->dumpFileValues();
  std::cout << "-----------------------------------------------------------" << std::endl;

  // Get configuration variables
  int plotVariance  = 0; // Whether to plot Converted-Reco/Reco
  int fitFromPeak   = 0; // Whether to fit a certain number of bins around the peak, rather than a particular range
  int nBinsFromPeak = 1; // How many bins to traverse either side of the peak in the fit
  int overlayFit    = 1; // Whether to overlay the fit result on the output histogram
  int zoom          = 0; // Whether to zoom in when plotting the ratio and variance
  double presMinX   = -99999; // Minimum X limit for the output
  double presMaxX   =  99999; // Maximum X limit for the output
  double zoomMinX   = -99999; // See above
  double zoomMaxX   =  99999;
  std::vector< std::string > inputs, hists, labels, titles;
  std::string location  ="";
  std::string tag="";
  std::string xAxis="dE/dx";
  std::string units="MeV/cm";
  std::string denomFile = "";
  std::string denomHist = "h_reco_dEdx_E_BP";
  std::vector<double> minimums, maximums;

  p->getValue("Inputs",        inputs);
  p->getValue("Hists",         hists);
  p->getValue("DenomFile",     denomFile);
  p->getValue("DenomHist",     denomHist);
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
  p->getValue("Zoom",          zoom);
  p->getValue("ZoomMinX",      zoomMinX);
  p->getValue("ZoomMaxX",      zoomMaxX);
  p->getValue("PresMinX",      presMinX);
  p->getValue("PresMaxX",      presMaxX);
  p->getValue("PlotVariance",  plotVariance);

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

  // Check if we're plotting the ratio
  bool plotRatio = false;
  if(denomFile != "")
    plotRatio = true;
  else if(denomFile == "" && plotVariance){
    std::cerr << " Error: Cannot plot variance if denominator file/histogram is not given" << std::endl;
    std::exit(1);
  }

  // Check if we should read the presentation limits from the configuration
  bool presLimits = false;
  if(presMinX > -99999){
    if(presMaxX < 99999)
      presLimits = true;
    else{
      std::cerr << " Error: Set the minimum X limit but not the maximum " << std::endl;
      std::exit(1);
    }
  }
  else if(presMaxX < 99999){
    if(presMinX > -99999)
      presLimits = true;
    else{
      std::cerr << " Error: Set the maximum X limit but not the minimum " << std::endl;
      std::exit(1);
    }
  }
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

    projections.push_back(hProj);
  }

  // If we are plotting the ratio, get the histogram
  TH2D *hDenom = nullptr;
  TH1D *hDenomProj = nullptr;
  if(plotRatio){
    TFile *f = new TFile(denomFile.c_str(),"READ");
    hDenom = static_cast<TH2D*>(f->Get(denomHist.c_str()));
  
    // Project out the x dimension
    std::string name = "h_proj_denom";
    hDenomProj = static_cast<TH1D*>(hDenom->ProjectionY(name.c_str(),1,hDenom->GetNbinsX()));
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

  TLegend *l = new TLegend(0.5,0.675,0.92,0.930);
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
    h->SetLineWidth(3);
    f->SetLineWidth(2);
    h->SetLineStyle(2);
    f->SetLineStyle(5);

    TString tLabel(labels.at(i).data());
    if(tLabel.Contains("true"))
      h->SetLineStyle(7);
    else if(tLabel.Contains("reco"))
      h->SetLineStyle(9);

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
  if(presLimits)
    projections.at(0)->GetXaxis()->SetRangeUser(presMinX,presMaxX);
  projections.at(0)->GetYaxis()->SetRangeUser(0,maxy*1.1);
  l->Draw();
  c->SaveAs((location+"/deposition_overlay"+tag+".root").c_str());
  c->SaveAs((location+"/deposition_overlay"+tag+".png").c_str());
  c->Clear();
  l->Clear();

  // Now do the ratios, if needed
  if(plotRatio){
    for(unsigned int i = 0; i < projections.size(); ++i){
      // Get the current histogram and clone it for the ratio
      TH1D *hProj = projections.at(i);
      std::string name = "h_ratio_"+labels.at(i);
      TH1D *hRatio = static_cast<TH1D*>(hProj->Clone(name.c_str()));
      
      // Divide!
      hDenomProj->Scale(1/static_cast<double>(hDenomProj->Integral()));
      hRatio->Divide(hDenomProj);

      SetHistogramStyle1D(hRatio, (xAxis+" ["+units+"]").c_str(), "Converted/Reconstructed");

      if(presLimits && !zoom)
        hRatio->GetXaxis()->SetRangeUser(presMinX,presMaxX);
      else if(zoom)
        hRatio->GetXaxis()->SetRangeUser(zoomMinX,zoomMaxX);
      hRatio->GetYaxis()->SetRangeUser(0,2);

      hRatio->SetLineColor(pal.at(i));
      hRatio->SetLineWidth(3);
      hRatio->SetLineStyle(2);

      // Set the precision of the mpv to print
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2) << mpvs.at(i);
      std::string s = stream.str();

      TString tLabel(labels.at(i).data());
      if(tLabel.Contains("true"))
        hRatio->SetLineStyle(7);
      else if(tLabel.Contains("reco"))
        hRatio->SetLineStyle(9);

      if(inputs.at(i) == denomFile){
        hRatio->SetLineWidth(2);
        hRatio->SetLineStyle(1);
        l->AddEntry(hRatio,("Nominal MPV: "+s+" "+units).c_str(), "l");
      }
      else
        l->AddEntry(hRatio,(titles.at(i)+" MPV: "+s+" "+units).c_str(), "l");

      if(i == 0)
        hRatio->Draw("hist");
      else
        hRatio->Draw("hist same");
    }
    l->Draw();
    if(zoom){
      c->SaveAs((location+"/deposition_ratio_overlay"+tag+"_zoom.root").c_str());
      c->SaveAs((location+"/deposition_ratio_overlay"+tag+"_zoom.png").c_str());
    }
    else{
      c->SaveAs((location+"/deposition_ratio_overlay"+tag+".root").c_str());
      c->SaveAs((location+"/deposition_ratio_overlay"+tag+".png").c_str());
    } 
    c->Clear();
    l->Clear();
  }
  
  // Now do the variance, if needed
  if(plotVariance){
    for(unsigned int i = 0; i < projections.size(); ++i){
      // Get the current histogram and clone it for the variance
      TH1D *hProj = projections.at(i);
      
      std::string name = "h_variance_"+labels.at(i);
      TH1D *hVariance = static_cast<TH1D*>(hProj->Clone(name.c_str()));
      
      // Calculate the variance
      hVariance->Add(hDenomProj,-1);
      hVariance->Divide(hDenomProj);

      SetHistogramStyle1D(hVariance, (xAxis+" ["+units+"]").c_str(), "(Converted-Reconstructed)/Reconstructed");

      if(presLimits && !zoom)
        hVariance->GetXaxis()->SetRangeUser(presMinX,presMaxX);
      else if(zoom)
        hVariance->GetXaxis()->SetRangeUser(zoomMinX,zoomMaxX);
      hVariance->GetYaxis()->SetRangeUser(-1,1);

      hVariance->SetLineColor(pal.at(i));
      hVariance->SetLineWidth(3);
      hVariance->SetLineStyle(2);

      // Set the precision of the mpv to print
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2) << mpvs.at(i);
      std::string s = stream.str();

      TString tLabel(labels.at(i).data());
      if(tLabel.Contains("true"))
        hVariance->SetLineStyle(7);
      else if(tLabel.Contains("reco"))
        hVariance->SetLineStyle(9);

      if(denomHist == hists.at(i)){
        hVariance->SetLineWidth(2);
        hVariance->SetLineStyle(1);
        l->AddEntry(hVariance,"Nominal", "l");
      }
      else
        l->AddEntry(hVariance,(titles.at(i)+" MPV: "+s+" "+units).c_str(), "l");

      if(i == 0)
        hVariance->Draw("hist");
      else
        hVariance->Draw("hist same");
    }
    l->Draw();
    if(zoom){
      c->SaveAs((location+"/deposition_variance_overlay"+tag+"_zoom.root").c_str());
      c->SaveAs((location+"/deposition_variance_overlay"+tag+"_zoom.png").c_str());
    }
    else{
      c->SaveAs((location+"/deposition_variance_overlay"+tag+".root").c_str());
      c->SaveAs((location+"/deposition_variance_overlay"+tag+".png").c_str());
    }
  }
  
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
