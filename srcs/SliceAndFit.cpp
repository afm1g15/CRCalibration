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
     
int sliceAndFit(const char *config){

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
  int nSlices = -1;
  double binWidths = -1;
  std::string inputFile = "";
  std::string inputHist = "";
  std::string location="";
  std::string tag="";
  std::vector<double> sliceMinX, sliceMaxX;
  std::vector<std::string> sliceHistLabel;

  p->getValue("InputFile",      inputFile);
  p->getValue("InputHist",      inputHist);
  p->getValue("NSlices",        nSlices);
  p->getValue("BinWidths",      binWidths);
  p->getValue("Location",       location);
  p->getValue("Tag",            tag);
  p->getValue("SliceMinX",      sliceMinX);
  p->getValue("SliceMaxX",      sliceMaxX);

  // Make sure at least the vectors have been filled or the number of bins and bin widths have been filled
  if((sliceMinX.size() + sliceMaxX.size()) == 0 && (nSlices == -1 || binWidths < 0)){
    std::cerr << " Error: You must either fill the vectors of bins to slice, or give BOTH the number of slices and the bin widths. " << std::endl;
    std::cerr << " MinX size: " << sliceMinX.size() << ", MaxX size: " << sliceMaxX.size() << ", nSlices: " << nSlices << ", binWidths: " << binWidths << std::endl;
    std::exit(1);
  }

  // Make sure the vectors are all the same length
  if(((sliceMinX.size() + sliceMaxX.size())/2) != sliceMaxX.size()){
    std::cerr << " Error: Vector sizes don't match, but they should." << std::endl;
    std::cerr << "    SliceMinX size: " << sliceMinX.size() << ", SliceMaxX size: " << sliceMaxX.size() << std::endl;
    std::exit(1);
  }

  // Read in the input histogram from the file
  TFile *histFile = new TFile(inputFile.c_str(),"READ");
  TH2D *h = static_cast<TH2D*>(histFile->Get(inputHist.c_str()));

  // If no specific bins have been passed, determine them from the number of bins and bin widths in the configuration
  if((sliceMinX.size() + sliceMaxX.size()) == 0){
    std::cout << " Slice vectors not given, therefore calculating bins from widths." << std::endl;
    FillSliceVectors(h,nSlices,binWidths,sliceMinX,sliceMaxX);
  }
  else
    std::cout << " Slice vectors given, therefore not calculating bins from widths." << std::endl;
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

  GetSliceLabels(sliceMinX,sliceMaxX,sliceHistLabel);

  // Now define the relevant histograms and fill them
  // Get the y range from the 2D histogram first
  double minY = h->GetYaxis()->GetXmin();
  double maxY = h->GetYaxis()->GetXmax();

  // Now define the histograms
  std::vector<TH1D*> sliceHists;
  //DefineHistograms(sliceHistLabel,minY,maxY,sliceHists);

  // And fill them
  FillHistograms(h,sliceMinX,sliceMaxX,sliceHistLabel,sliceHists);

  // Now draw and save
  TFile *f = new TFile((location+"/slice_histograms"+tag+".root").c_str(), "RECREATE");
  TCanvas *c = new TCanvas("c","",900,900);
  SetCanvasStyle(c, 0.12,0.1,0.05,0.12,0,0,0);
  f->cd();

  for(unsigned int i = 0; i < sliceHists.size(); ++i){
    // Rename the canvas
    c->SetName(("c_"+sliceHistLabel.at(i)).c_str());

    // Sort out the histogram
    SetHistogramStyle1D(sliceHists.at(i),"Energy per charge deposition [MeV/ADC]", ("Rate, slice "+sliceHistLabel.at(i)).c_str());

    // Draw and save
    sliceHists.at(i)->SetLineWidth(2);
    sliceHists.at(i)->SetLineColor(pal.at(i));
    sliceHists.at(i)->SetLineStyle(styles.at(i));
    sliceHists.at(i)->Draw("hist");
    sliceHists.at(i)->Write();
    c->Write();
    c->Clear();
  } // Histograms and canvases

  TLegend *l = new TLegend(0.52,0.22,0.94,0.94);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->SetTextFont(132);
 
  // Sort out limits
  double maxy = -999.;
  c->SetName("c_overlay");
  for(unsigned int i = 0; i < sliceHists.size(); ++i){
    // Scale the histograms
    sliceHists.at(i)->Scale(1/static_cast<double>(sliceHists.at(i)->Integral()));
    if(sliceHists.at(i)->GetMaximum() > maxy)
      maxy = sliceHists.at(i)->GetMaximum();
  
    if(i == 0)
      sliceHists.at(i)->Draw("hist");
    else
      sliceHists.at(i)->Draw("hist same");

    // Add to the legend
    l->AddEntry(sliceHists.at(i),sliceHistLabel.at(i).c_str(),"l");
  } // For the overlay
  sliceHists.at(0)->GetYaxis()->SetRangeUser(0,maxy*1.1);
  l->Draw("same");
  c->Write();
  c->Clear();
  
  // Now fit the slices to Landau distributions and calculate the MPVs
  TH1D *mpv_x = new TH1D("MPV_vs_X","",100,-800,800);
  SetHistogramStyle1D(mpv_x,"X [cm]","Energy per charge deposition, MPV [MeV/ADC]");
  c->SetName("mpv_x");
  for(unsigned int i = 0; i < sliceHists.size(); ++i){
    double sliceCentre = sliceMinX.at(i)+((sliceMaxX.at(i)-sliceMinX.at(i))/2.);

    // Now define the TF1
    TF1 *fit = new TF1("fit","landau", 6.7e-3,7.8e-3);
    if(sliceHists.at(i)->Fit("fit", "Q L M R")){
      double mpv = fit->GetParameter(1);
      std::cout << " Bin centre: " << sliceCentre << ", MPV: " << mpv << std::endl;
      mpv_x->Fill(sliceCentre,mpv);
    }
  } // Loop for fits
  mpv_x->SetMarkerColor(pal.at(0));
  mpv_x->SetMarkerStyle(33);
  mpv_x->GetYaxis()->SetRangeUser(minY,maxY);
  mpv_x->Draw("P hist");
  mpv_x->Write();
  c->Write();

  std::cout << " Writing all slices to file: " << (location+"/slice_histograms"+tag+".root").c_str() << std::endl;

  f->Write();
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
