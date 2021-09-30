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
  std::string inputTrue = "";
  std::string inputReco = "";
  std::string inputConv = "";
  std::string trueHist  = "";
  std::string recoHist  = "";
  std::string convHist  = "";
  std::string location  ="";
  std::string tag="";
  std::vector<double> trueRange, recoRange, convRange;

  p->getValue("InputTrue", inputTrue);
  p->getValue("InputReco", inputReco);
  p->getValue("InputConv", inputConv);
  p->getValue("TrueHist",  trueHist);
  p->getValue("RecoHist",  recoHist);
  p->getValue("ConvHist",  convHist);
  p->getValue("Location",  location);
  p->getValue("Tag",       tag);
  p->getValue("TrueRange", trueRange);
  p->getValue("RecoRange", recoRange);
  p->getValue("ConvRange", convRange);

  // Read in the input histograms from the files
  TFile *trueFile = new TFile(inputTrue.c_str(),"READ");
  TH2D  *hTrue    = static_cast<TH2D*>(trueFile->Get(trueHist.c_str()));
  TFile *recoFile = new TFile(inputReco.c_str(),"READ");
  TH2D  *hReco    = static_cast<TH2D*>(recoFile->Get(recoHist.c_str()));
  TFile *convFile = new TFile(inputConv.c_str(),"READ");
  TH2D  *hConv    = static_cast<TH2D*>(convFile->Get(convHist.c_str()));

  // Sort out the file tag
  if(tag != "")
    tag = "_"+tag;

  std::vector< std::vector<double> > ranges = {trueRange, recoRange, convRange};
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
  
  // Now project each on along the entire x dimension
  TH1D *hProjTrue = static_cast<TH1D*>(hTrue->ProjectionY("h_true_dEdx",1,hTrue->GetNbinsX()));
  TH1D *hProjReco = static_cast<TH1D*>(hReco->ProjectionY("h_reco_dEdx",1,hReco->GetNbinsX()));
  TH1D *hProjConv = static_cast<TH1D*>(hConv->ProjectionY("h_conv_dEdx",1,hConv->GetNbinsX()));

  std::vector<TH1D*>      hProj = {hProjTrue,hProjReco,hProjConv};
  std::vector<std::string> hLab = {"True dE/dx", "Reconstructed dE/dx", "Converted dE/dx"};
  std::vector<TF1*> hFits;
  std::vector<std::string> fLab = {"True fit", "Reconstructed fit", "Converted fit"};
  std::vector<double> mpvs;

  unsigned int n = 0;
  for(TH1D* h : hProj){
    // Get the y limits
    double minX   = h->GetXaxis()->GetXmin();
    double maxX   = h->GetXaxis()->GetXmax();
    int maxbin    = h->GetMaximumBin();
    double maxloc = h->GetBinCenter(maxbin);
    std::string name(hLab.at(n)); 
    
    // Now fit a langaus to each histogram and extract the MPV
    TF1 *fit = new TF1(("fit_"+name).c_str(),langaufun,ranges.at(n).at(0),ranges.at(n).at(1),4);
    fit->SetParNames("Width","MP","Area","GSigma");
    double norm = h->GetEntries() * h->GetBinWidth(1);
    double sv[4] = {10., maxloc, norm, 10.}; // starting values for parameters: Landau scale, Landau MPV, Norm, Gauss sigma
    fit->SetParameters(sv);

    auto result = h->Fit(fit, "QSLR");
    double mpv = fit->GetMaximumX(result->Parameter(1), result->Parameter(1) + result->Parameter(3));

    ofile << name << " MPV: " << mpv << std::endl;

    mpvs.push_back(mpv);
    hFits.push_back(fit);
    ++n;
  }

  TCanvas *c = new TCanvas("c","",900,900);
  SetCanvasStyle(c, 0.12,0.05,0.05,0.12,0,0,0);

  TLegend *l = new TLegend(0.556,0.637,0.998,0.922);
  l->SetTextFont(132);
  l->SetBorderSize(0);
  l->SetFillStyle(0);

  double maxy = -999;
  for(unsigned int i = 1; i < hProj.size(); ++i){
    TH1D *h = hProj.at(i);
    TF1  *f = hFits.at(i);

    SetHistogramStyle1D(h,"dE/dx [MeV/cm]", "Rate");

    if(h->GetMaximum() > maxy)
      maxy = h->GetMaximum();

    h->SetLineColor(pal.at(i));
    f->SetLineColor(pal.at(i));
    h->SetLineWidth(2);
    f->SetLineWidth(1);
    h->SetLineStyle(1);
    f->SetLineStyle(7);

    l->AddEntry(h,hLab.at(i).c_str(), "l");
    l->AddEntry(f,fLab.at(i).c_str(), "l");

    if(i == 1)
      h->Draw("hist");
    else
      h->Draw("hist same");
    f->Draw("hist same");
  }
  hProj.at(1)->GetYaxis()->SetRangeUser(0,maxy*1.1);
  l->Draw();
  c->SaveAs((location+"/energy_deposition_overlay"+tag+".root").c_str());
  c->SaveAs((location+"/energy_deposition_overlay"+tag+".png").c_str());
  
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
