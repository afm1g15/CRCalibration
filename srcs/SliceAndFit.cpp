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

  GetSliceLabels(sliceMinX,sliceMaxX,sliceHistLabel);

  // Now define the relevant histograms and fill them
  // Get the y range from the 2D histogram first
  double minX = h->GetXaxis()->GetXmin();
  double maxX = h->GetXaxis()->GetXmax();
  double minY = h->GetYaxis()->GetXmin();
  double maxY = h->GetYaxis()->GetXmax();

  ofile << " X: (" << minX << ", " << maxX << ")" << std::endl;
  ofile << " Y: (" << minY << ", " << maxY << ")" << std::endl;

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
    SetHistogramStyle1D(sliceHists.at(i),"Charge deposition [ADC/cm]", ("Rate, "+sliceHistLabel.at(i)).c_str());

    // Draw and save
    sliceHists.at(i)->SetLineWidth(2);
    sliceHists.at(i)->SetLineColor(pal.at(i));
    sliceHists.at(i)->SetLineStyle(styles.at(i));
    sliceHists.at(i)->Draw("hist");
    sliceHists.at(i)->Write();
    c->Write();
    c->SaveAs((location+"/slice_"+sliceHistLabel.at(i)+"_"+tag+".png").c_str());
    c->SaveAs((location+"/slice_"+sliceHistLabel.at(i)+"_"+tag+".root").c_str());
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
  c->SaveAs((location+"/slice_overlay"+tag+".root").c_str());
  c->SaveAs((location+"/slice_overlay"+tag+".png").c_str());
  c->Clear();
  
  // Now fit the slices to Landau distributions and calculate the MPVs
  TH1D *mpv_x = new TH1D("MPV_vs_X","",100,minX,maxX);
  SetHistogramStyle1D(mpv_x,"X [cm]","Energy per charge deposition, MPV [MeV/ADC]");
  c->SetName("mpv_x");

  // Also find the minimum and maximum MPV, the average MPV,
  // and get the min-max difference wrt the average
  double minMPV = 999.;
  double maxMPV = -999.;
  double avgMPV = 0.;
  for(unsigned int i = 0; i < sliceHists.size(); ++i){
    double sliceCentre = sliceMinX.at(i)+((sliceMaxX.at(i)-sliceMinX.at(i))/2.);
    int maxbin         = sliceHists.at(i)->GetMaximumBin();
    double maxloc      = sliceHists.at(i)->GetBinCenter(maxbin);

    // Now define the TF1
    // Set some approximate start parameters
    TF1 *fit = new TF1("fit",langaufun,minY,maxY,4);
    fit->SetParNames("Width","MP","Area","GSigma");
    double norm = sliceHists.at(i)->GetEntries() * sliceHists.at(i)->GetBinWidth(1);
    double sv[4] = {10., maxloc, norm, 10.}; // starting values for parameters: Landau scale, Landau MPV, Norm, Gauss sigma
    sv[1] = maxloc;
    sv[2] = norm;
    fit->SetParameters(sv);

    auto result = sliceHists.at(i)->Fit(fit, "QSMR", "");
    double mpv = fit->GetMaximumX(result->Parameter(1), result->Parameter(1) + result->Parameter(3));
    ofile << " Bin centre: " << sliceCentre << ", MPV: " << mpv << std::endl;
    mpv_x->Fill(sliceCentre,mpv);
    if(mpv > maxMPV)
      maxMPV = mpv;
    if(mpv < minMPV)
      minMPV = mpv;
    avgMPV += mpv;

    // Rename the canvas
    c->SetName(("c_fit_"+sliceHistLabel.at(i)).c_str());

    // Draw and save
    sliceHists.at(i)->SetLineWidth(2);
    sliceHists.at(i)->SetLineColor(pal.at(i));
    sliceHists.at(i)->SetLineStyle(styles.at(i));
    fit->SetLineWidth(1);
    fit->SetLineColor(pal.at(i));
    fit->SetLineStyle(7);
    sliceHists.at(i)->Draw("hist");
    fit->Draw("same");
    c->Write();
    c->SaveAs((location+"/fit_slice_"+sliceHistLabel.at(i)+"_"+tag+".png").c_str());
    c->SaveAs((location+"/fit_slice_"+sliceHistLabel.at(i)+"_"+tag+".root").c_str());
    c->Clear();
  } // Loop for fits
  // Now calculcate the fractional MPV difference 
  avgMPV /= static_cast<double>(sliceHists.size());
  double mpvDiff = maxMPV - minMPV;
  double fracMPVDiff = mpvDiff/avgMPV;

  ofile << " The minimum MPV is: " << minMPV << ", the maximum MPV is: " << maxMPV << std::endl;
  ofile << " The difference between the max and min MPV is           : " << mpvDiff << std::endl;
  ofile << " The average MPV is                                      : " << avgMPV << std::endl;
  ofile << " The fractional difference between the max and min MPV is: " << fracMPVDiff << std::endl;

  c->SetRightMargin(0.05);
  mpv_x->SetMarkerColor(pal.at(0));
  mpv_x->SetMarkerStyle(33);
  mpv_x->GetYaxis()->SetRangeUser(minY,maxY);
  mpv_x->Draw("P hist");
  mpv_x->Write();
  c->Write();
  c->SaveAs((location+"/mpv_x"+tag+".root").c_str());
  c->SaveAs((location+"/mpv_x"+tag+".png").c_str());
  c->Clear();

  // Now fit a straight line to the MPV vs x distribution
  TF1 *fitLine = new TF1("fitLine","[0]+[1]*x",-800,800);
  // Start the fit at the average MPV
  mpv_x->Fit("fitLine", "Q R");

  double scaleFactor = 2.12/fitLine->GetParameter(0);
  ofile << " Constant : " << fitLine->GetParameter(0) << std::endl;
  ofile << " Gradient : " << fitLine->GetParameter(1) << std::endl;
  ofile << " ChiSquare: " << fitLine->GetChisquare() << std::endl;
  ofile << " ----------------------------------------" << std::endl;
  ofile << " Conversion factor: 2.12 [MeV/cm]/MPV [ADC/cm] = " << scaleFactor << " [MeV/ADC] " << std::endl;

  // Draw the mpv_x and the fit distribution
  c->SetName("mpv_x_fit_line");
  mpv_x->Draw("P hist");
  fitLine->SetLineColor(pal.at(0));
  fitLine->SetLineStyle(7);
  fitLine->SetLineWidth(1);
  fitLine->Draw("same");
  c->Write();
  c->SaveAs((location+"/mpv_x_fit"+tag+".root").c_str());
  c->SaveAs((location+"/mpv_x_fit"+tag+".png").c_str());
  c->Clear();

  // Now overlay the mpvs with the 2D
  TCanvas *c1 = new TCanvas("c1","",1000,800);
  SetCanvasStyle(c1, 0.1,0.12,0.05,0.12,0,0,0);

  c1->SetName("mpv_x_2D_overlay");
  h->Draw("colz");
  mpv_x->Draw("P hist same");
  c1->SaveAs((location+"/mpv_x_2D_overlay"+tag+".root").c_str());
  c1->SaveAs((location+"/mpv_x_2D_overlay"+tag+".png").c_str());
  c1->Write();
  c1->Clear();

  // Now copy the input histogram and scale with the conversion factor
  c1->SetName("h_converted");
  double minBin = minY*scaleFactor;
  double maxBin = maxY*scaleFactor;
  TH2D *h_conv = new TH2D("h_dedx_from_dqdx_x","",h->GetNbinsX(),minX,maxX,h->GetNbinsY(),minBin,maxBin);
  SetHistogramStyle2D(h_conv,"x [cm]", " dE/dx [MeV/cm]", false);

  // Loop over the y axis and scale the bins
  for(int nX = 1; nX <= h_conv->GetNbinsX(); ++nX){
    for(int nY = 1; nY <= h_conv->GetNbinsY(); ++nY){
      double chargeContent = h->GetBinContent(nX,nY);
      h_conv->SetBinContent(nX,nY,chargeContent);
    } // NBinsY
  } // NBinsX
  h_conv->Draw("colz");
  c1->SaveAs((location+"/converted_hist_2D"+tag+".root").c_str());
  c1->SaveAs((location+"/converted_hist_2D"+tag+".png").c_str());
  c1->Write();

  std::cout << " Writing all slices to file: " << (location+"/slice_histograms"+tag+".root").c_str() << std::endl;

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
