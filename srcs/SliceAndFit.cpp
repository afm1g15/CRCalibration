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
  int nSlices            = -1;
  int fitRange           = 1; // Whether or not to fit the full range of the function
  int fitExp             = 1; // Whether to fit the MPV's to an exponential
  int fitPol2            = 0; // Whether to fit the MPV's to a second order polynomial
  int logSpace           = 0; // Whether to spread the bins in log space
  int logX               = 0; // Whether the x dimension is in log space 
  int convert            = 1; // Whether to convert the histogram from charge to energy space
  int fitFromPeak        = 0; // Whether to define the fit from the peak outwards
  int nBinsFromPeak      = 1; // How many bins to traverse either side of the peak in the fit
  int fitBelow           = 0; // Whether or not to fit below to specified minimum for the chosen function
  int drawSliceLines     = 0; // Whether to draw the slice lines on the 2D orig histogram to indicate their size and location
  int rebin              = 0; // Should we rebin the slices
  int constScale         = 0; // Should we convert the histogram with a constant (1) or a function (0)
  int nOverlay           = -1; // Number of slices to overlay, -1 = all
  double buffer          = 0;
  double binWidths       = -1; // I think this should be the percentage/fraction of the space
  double fitMin          = 99999.;
  double fitMax          = -99999.;
  double mpvFitMin       = 99999.;
  double mpvFitMax       = -99999.;
  double nominalMPV      = 1.761; // dE/dx MPV for a 292.264 GeV muon in 3.53 mm thickness (both values from truth MPVs) from here: https://lar.bnl.gov/properties/
  std::string fitFunc    = "langaus";
  std::string lowFitFunc = "lin";
  std::string units      = "";
  std::string inputFile  = "";
  std::string inputHist  = "";
  std::string location   = "";
  std::string tag        = "";
  std::vector<double> sliceMinX, sliceMaxX;
  std::vector<std::string> sliceHistLabel, sliceHistLabelTeX;

  p->getValue("InputFile",       inputFile);
  p->getValue("InputHist",       inputHist);
  p->getValue("NSlices",         nSlices);
  p->getValue("Rebin",           rebin);
  p->getValue("BinWidths",       binWidths);
  p->getValue("DrawSliceLines",  drawSliceLines);
  p->getValue("Location",        location);
  p->getValue("Tag",             tag);
  p->getValue("Buffer",          buffer);
  p->getValue("FitMin",          fitMin);
  p->getValue("FitMax",          fitMax);
  p->getValue("FitBelow",        fitBelow);
  p->getValue("MPVFitBufferMin", mpvFitMin);
  p->getValue("MPVFitBufferMax", mpvFitMax);
  p->getValue("FitRange",        fitRange);
  p->getValue("FitFromPeak",     fitFromPeak);
  p->getValue("NBinsFromPeak",   nBinsFromPeak);
  p->getValue("FitFunction",     fitFunc);
  p->getValue("BelowFitFunc",    lowFitFunc);
  p->getValue("LogSpace",        logSpace);
  p->getValue("ConstScale",      constScale);
  p->getValue("LogX",            logX);
  p->getValue("FitExp",          fitExp);
  p->getValue("FitPol2",         fitPol2);
  p->getValue("Convert",         convert);
  p->getValue("SliceMinX",       sliceMinX);
  p->getValue("SliceMaxX",       sliceMaxX);
  p->getValue("NOverlay",        nOverlay);
  p->getValue("Units",           units);
  p->getValue("NominalMPV",      nominalMPV);

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
    FillSliceVectors(h,nSlices,binWidths,logSpace,sliceMinX,sliceMaxX,buffer);
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
  GetSliceLabelsTeX(sliceMinX,sliceMaxX,sliceHistLabelTeX,units);

  // Now define the relevant histograms and fill them
  // Get the y range from the 2D histogram first
  double minX = h->GetXaxis()->GetXmin();
  double maxX = h->GetXaxis()->GetXmax();
  double minY = h->GetYaxis()->GetXmin();
  double maxY = h->GetYaxis()->GetXmax();

  std::string xLabel = static_cast<std::string>(h->GetXaxis()->GetTitle());

  ofile << " X: (" << minX << ", " << maxX << ")" << std::endl;
  ofile << " Y: (" << minY << ", " << maxY << ")" << std::endl;

  // Now define the histograms
  std::vector<TH1D*> sliceHists;

  // And fill them
  FillHistograms(h,sliceMinX,sliceMaxX,sliceHistLabel,sliceHists);

  // Now draw and save
  TFile *f = new TFile((location+"/slice_histograms"+tag+".root").c_str(), "RECREATE");
  TCanvas *c = new TCanvas("c","",900,900);
  SetCanvasStyle(c, 0.12,0.08,0.05,0.12,0,0,0);
  f->cd();

  for(unsigned int i = 0; i < sliceHists.size(); ++i){
    // Rename the canvas
    c->SetName(("c_"+sliceHistLabel.at(i)).c_str());

    // Sort out the histogram
    SetHistogramStyle1D(sliceHists.at(i),h->GetYaxis()->GetTitle(), ("Rate, "+sliceHistLabelTeX.at(i)).c_str());

    // Draw and save
    // Rebin to account for low stats
    if(rebin)
      sliceHists.at(i)->Rebin(2);
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

  // Define the Y1 position according to how many slices we're overlaying
  double Y1 = 0.08;
  double Y2 = 0.92;
  if(nOverlay > 0){
    double yDiff = Y2-Y1;
    double dPerSlice = yDiff/sliceHists.size();
    Y1 = Y2-(dPerSlice*nOverlay);
  }
  TLegend *l = new TLegend(0.58,Y1,0.94,Y2);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->SetTextFont(132);
  l->SetTextSize(0.028);
 
  // Sort out limits
  double maxy = -999.;
  if(nOverlay != 0){
    c->SetName("c_overlay");
    // Checker for the number to overlay
    unsigned int nAdded = 0;
    unsigned int nSkip = 0;
    // If the number to overlay is not 'all (-1)' then calculate how many to skip
    if(nOverlay > 0){
      nSkip = std::round(sliceHists.size()/static_cast<double>(nOverlay-1)); 
    }
    for(unsigned int i = 0; i < sliceHists.size(); ++i){
      // Determine whether to include this in the overlay
      if(nOverlay > 0){
        unsigned int multipleCheck = i+nSkip;
        
        // If it's not a multiple of nSkip, continue
        if(multipleCheck % nSkip != 0 && i != sliceHists.size()-1){
          continue;
        }
      }
      // Scale the histograms
      sliceHists.at(i)->Scale(1/static_cast<double>(sliceHists.at(i)->Integral()));
      if(sliceHists.at(i)->GetMaximum() > maxy)
        maxy = sliceHists.at(i)->GetMaximum();

      if(i == 0)
        sliceHists.at(i)->Draw("hist");
      else
        sliceHists.at(i)->Draw("hist same");

      // Add to the legend
      l->AddEntry(sliceHists.at(i),sliceHistLabelTeX.at(i).c_str(),"l");
    } // For the overlay
    sliceHists.at(0)->GetYaxis()->SetTitle("Normalised rate");
    sliceHists.at(0)->GetYaxis()->SetRangeUser(0,maxy*1.1);
    l->Draw("same");
    c->Write();
    c->SaveAs((location+"/slice_overlay"+tag+".root").c_str());
    c->SaveAs((location+"/slice_overlay"+tag+".png").c_str());
    c->Clear();
  }
  
  // Now fit the slices to Landau distributions and calculate the MPVs
  TH1D *mpv_x = new TH1D("MPV_vs_X","",100,minX,maxX);
  SetHistogramStyle1D(mpv_x,h->GetXaxis()->GetTitle(),h->GetYaxis()->GetTitle());
  if(logX)
    SetLogX(mpv_x);
  c->SetName("mpv_x");

  // Also find the minimum and maximum MPV, the average MPV,
  // and get the min-max difference wrt the average
  double minMPV       = 999.;
  double maxMPV       = -999.;
  double avgMPV       = 0.;
  double totFracSigma = 0.;
  for(unsigned int i = 0; i < sliceHists.size(); ++i){
    double sliceCentre = sliceMinX.at(i)+((sliceMaxX.at(i)-sliceMinX.at(i))/2.);
    int maxbin         = sliceHists.at(i)->GetMaximumBin();
    double maxloc      = sliceHists.at(i)->GetBinCenter(maxbin);

    // Now define the TF1
    TF1 *fit;
    double minR = minY;
    double maxR = maxY;

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
        minR = sliceHists.at(i)->GetBinCenter(maxbin-nBinsFromPeak);
      if(maxbin+nBinsFromPeak < sliceHists.at(i)->GetNbinsX())
        maxR = sliceHists.at(i)->GetBinCenter(maxbin+nBinsFromPeak);
    }

    if(fitFunc == "langaus"){
      // Set some approximate start parameters
      fit = new TF1("fit",langaufun,minR,maxR,4);
      fit->SetParNames("Width","MP","Area","GSigma");
      double norm = sliceHists.at(i)->Integral("width"); //sliceHists.at(i)->GetEntries() * sliceHists.at(i)->GetBinWidth(1); //
      double sv[4] = {0.2*maxloc, maxloc, norm, 0.2*maxloc}; // starting values for parameters: Landau scale, Landau MPV, Norm, Gauss sigma
      fit->SetParameters(sv);

    }
    else{
      // Set some approximate start parameters
      fit = new TF1("fit",fitFunc.c_str(),minR,maxR);
    }

    // Now get the results
    auto result = sliceHists.at(i)->Fit(fit, "QSMR", "");
    double mpv = result->Parameter(1);

    // Find the x value corresponding to the maximum of the function in the interval MPV, MPV+GSig
    if(fitFunc == "langaus")
      mpv = fit->GetMaximumX(result->Parameter(1), result->Parameter(1) + result->Parameter(3));

    for(unsigned int p = 0; p < result->NPar(); ++p){
      ofile << " " << result->ParName(p) << " : " << result->Parameter(p) << std::endl;
    }

    // Get the error on the mpv from both the fit scaled with the error due to the statistics
    double tot_error = result->Parameter(3)/static_cast<double>(sqrt(sliceHists.at(i)->GetNbinsX()));
    mpv_x->SetBinContent(mpv_x->FindBin(sliceCentre),mpv);
    mpv_x->SetBinError(mpv_x->FindBin(sliceCentre),tot_error);

    if(mpv > maxMPV)
      maxMPV = mpv;
    if(mpv < minMPV)
      minMPV = mpv;
    avgMPV += mpv;
    totFracSigma += pow(tot_error/static_cast<double>(mpv),2);

    // Rename the canvas
    c->SetName(("c_fit_"+sliceHistLabel.at(i)).c_str());

    // Draw and save
    FormatStats(sliceHists.at(i),10,10);
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
  double avgSigma = avgMPV*sqrt(totFracSigma);
  
  double mpvDiff = maxMPV - minMPV;
  double fracMPVDiff = mpvDiff/avgMPV;

  ofile << " The minimum MPV is: " << minMPV << ", the maximum MPV is: " << maxMPV << std::endl;
  ofile << " The difference between the max and min MPV is           : " << mpvDiff << std::endl;
  ofile << " The average MPV is                                      : " << avgMPV << " +/- " << avgSigma << std::endl;
  ofile << " The fractional difference between the max and min MPV is: " << fracMPVDiff << std::endl;

  // Now fit a straight line to the MPV vs x parameter distribution
  double mpvMin = minX; 
  double mpvMax = maxX;
  if(mpvFitMin < 99999)
    mpvMin = mpvFitMin;
  if(mpvFitMax > -99999)
    mpvMax = mpvFitMax;

  TF1 *fitLine;
  if(logSpace && fitExp){
    fitLine = new TF1("fitLine","[0]-[1]*exp(-[2]*x)",mpvMin,mpvMax);
  }
  else if(fitPol2){
    fitLine = new TF1("fitLine","[0]+[1]*x+[2]*x*x",mpvMin,mpvMax);
  }
  else{
    fitLine = new TF1("fitLine","[0]+[1]*x",mpvMin,mpvMax);
  }
  auto mpv_result = mpv_x->Fit(fitLine, "QSMR0");
  for(unsigned int p = 0; p < mpv_result->NPar(); ++p){
    ofile << " " << mpv_result->ParName(p) << " : " << mpv_result->Parameter(p) << std::endl;
  }

  // If we also want to fit below the minimum value, do that with a straight line
  TF1 *lowFit;
  if(fitBelow){
    if(lowFitFunc == "lin")
      lowFit = new TF1("lowFit","[0]+[1]*x",sliceMinX.at(0),mpvMin);
    else if(lowFitFunc == "exp")
      lowFit = new TF1("lowFit","[0]-[1]*exp(-[2]*x)",sliceMinX.at(0),mpvMin);
    else
      lowFit = new TF1("fitLine","[0]",mpvMin,mpvMin);
  }
  else{
    lowFit = new TF1("fitLine","[0]",mpvMin,mpvMin);
  }
  auto mpv_low_result = mpv_x->Fit(lowFit, "QSMR0");
  
  c->SetName("mpv_x");
  c->SetRightMargin(0.05);
  if(logSpace)
    c->SetLogx();
  mpv_x->SetMarkerColor(pal.at(0));
  mpv_x->SetLineColor(pal.at(0));
  mpv_x->SetMarkerStyle(1);
  mpv_x->GetYaxis()->SetRangeUser(minY,maxY);
  mpv_x->Draw("P E1 X0");
  c->Write();
  c->SaveAs((location+"/mpv_x"+tag+".root").c_str());
  c->SaveAs((location+"/mpv_x"+tag+".png").c_str());
  c->Clear();

  double scaleFactor = nominalMPV/avgMPV;
  ofile << " ----------------------------------------" << std::endl;
  ofile << " Conversion factor: " << nominalMPV << " [MeV/cm]/MPV [ADC/cm] = " << scaleFactor << " [MeV/ADC] " << std::endl;
  ofile << " ----------------------------------------" << std::endl;
  ofile << " Fit between : " << mpvMin << ", " << mpvMax << std::endl;
  ofile << " Constant    : " << fitLine->GetParameter(0) << std::endl;
  ofile << " Gradient    : " << fitLine->GetParameter(1) << std::endl;
  ofile << " ChiSquare   : " << fitLine->GetChisquare() << std::endl;
  ofile << " NDOF        : " << fitLine->GetNDF() << std::endl;
  ofile << " ----------------------------------------" << std::endl;
  if(fitBelow){
    ofile << " Fit between : " << mpv_x->GetXaxis()->GetXmin() << ", " << mpvMin << std::endl;
    ofile << " Constant    : " << lowFit->GetParameter(0) << std::endl;
    ofile << " Gradient    : " << lowFit->GetParameter(1) << std::endl;
    ofile << " ChiSquare   : " << lowFit->GetChisquare() << std::endl;
    ofile << " NDOF        : " << lowFit->GetNDF() << std::endl;
    ofile << " ----------------------------------------" << std::endl;
  }

  // If we are fitting a single linear function, set constScale if the gradient is < 1e-2
  if(!fitBelow && !fitPol2 && !fitExp){
    if(abs(fitLine->GetParameter(1)) < 1e-2){
      std::cout << " Fit a line and the gradient is " << fitLine->GetParameter(1) << ", converting with a constant scale factor" << std::endl;
      constScale = true;
    }
  }
  // Draw the mpv_x and the fit distribution
  c->SetName("mpv_x_fit_line");
  mpv_x->Draw("P E1 X0");
  fitLine->SetLineColor(pal.at(0));
  fitLine->SetLineStyle(7);
  fitLine->SetLineWidth(1);
  fitLine->Draw("same");
  if(fitBelow){
    lowFit->SetLineColor(pal.at(1));
    lowFit->SetLineStyle(7);
    lowFit->SetLineWidth(1);
    lowFit->Draw("same");
  }
  c->Write();
  c->SaveAs((location+"/mpv_x_fit"+tag+".root").c_str());
  c->SaveAs((location+"/mpv_x_fit"+tag+".png").c_str());
  c->Clear();

  // Now overlay the mpvs with the 2D
  TCanvas *c1 = new TCanvas("c1","",1000,800);
  SetCanvasStyle(c1, 0.1,0.15,0.05,0.12,0,0,0);

  if(logSpace)
    c1->SetLogx();
  c1->SetName("mpv_x_2D_overlay");
  SetHistogramStyle2D(h,h->GetXaxis()->GetTitle(), h->GetYaxis()->GetTitle(), false);
  h->Draw("colz");
  mpv_x->SetMarkerColor(0);
  mpv_x->SetLineColor(0);
  mpv_x->Draw("P E2 X0 same");
  c1->SaveAs((location+"/mpv_x_2D_overlay"+tag+".root").c_str());
  c1->SaveAs((location+"/mpv_x_2D_overlay"+tag+".png").c_str());
  c1->Write();
  c1->Clear();

  if(convert){
    // Now copy the input histogram and scale with the conversion factor
    c1->SetName("h_converted");
    double k = fitLine->GetParameter(0);
    double minBin = minY*scaleFactor;
    double maxBin = maxY*scaleFactor;

    // If we simply want to scale the y axis by the constant from the fit, do that
    if(constScale){
      TH2D *h_conv = new TH2D("h_dedx_from_dqdx","",h->GetNbinsX(),minX,maxX,h->GetNbinsY(),minBin,maxBin);
      if(logX)
        SetLogX(h_conv);
      SetHistogramStyle2D(h_conv, xLabel.c_str(), " dE/dx [MeV/cm]", false);

      // Loop over the y axis and scale the bins
      for(int nX = 1; nX <= h_conv->GetNbinsX(); ++nX){
        for(int nY = 1; nY <= h_conv->GetNbinsY(); ++nY){
          double chargeContent = h->GetBinContent(nX,nY);
          h_conv->SetBinContent(nX,nY,chargeContent);
        } // NBinsY
      } // NBinsX
      h_conv->Draw("colz");
      h_conv->SaveAs((location+"/hist_converted_hist_2D"+tag+".root").c_str());
      h_conv->SaveAs((location+"/hist_converted_hist_2D"+tag+".png").c_str());
    }
    else{ // Otherwise apply the scale factor in its functional form
      double m = fitLine->GetParameter(1);
      double p = 0;
      // If second order, get the 2nd coefficient
      if(fitPol2)
        p = fitLine->GetParameter(2);

      double kLow  = k;
      double mLow  = m;
      double phase = 0; 
      if(fitBelow){
        kLow = lowFit->GetParameter(0);
        mLow = lowFit->GetParameter(1);
        if(lowFitFunc == "exp")
          phase = lowFit->GetParameter(2);
      }

      TH2D *h_conv = new TH2D("h_dedx_from_dqdx","",h->GetNbinsX(),minX,maxX,h->GetNbinsY(),minBin,maxBin);
      if(logX)
        SetLogX(h_conv);
      SetHistogramStyle2D(h_conv, xLabel.c_str(), " dE/dx [MeV/cm]", false);

      // Loop over the y axis and scale the bins
      for(int nX = 0; nX <= h_conv->GetNbinsX(); ++nX){
        // Check if we are looking above or below 20 GeV
        bool above20 = true;
        double currVal = h->GetXaxis()->GetBinCenter(nX);
        if(currVal < 20) above20 = false;

        // Now loop over the y bins and sort those out
        for(int nY = 0; nY <= h_conv->GetNbinsY(); ++nY){
          double content    = h->GetBinContent(nX,nY);
          double yCentre    = h->GetYaxis()->GetBinCenter(nY);

          // Apply the more involved scaling from the fit above
          // depending on where abouts we are in the parameter space
          //
          //   Q/L   = k + m.E
          //   C     = nominalMPV / (Q/L) =  nominalMPV / (k + m.E)
          //   dE/dx = C.dQ/dx      = (nominalMPV / (k + m.E))*dQ/dx
          //
          double newYCentre  = yCentre*nominalMPV;
          if(fitBelow && !above20){
            newYCentre = yCentre * ( nominalMPV / ( kLow + mLow*currVal ) );
            if(lowFitFunc == "exp")
              newYCentre = yCentre * ( nominalMPV / ( kLow - mLow*exp(-phase*currVal ) ) );
          }
          else
            newYCentre = yCentre * ( nominalMPV / ( k + m*currVal + p*pow(currVal,2)) );

          int nEntries = std::ceil(content);
          for(int n = 0; n < nEntries; ++n){
            h_conv->Fill(currVal, newYCentre);
          }
        } // NBinsY
      } // NBinsX
      h_conv->Draw("colz");
      h_conv->SaveAs((location+"/hist_converted_hist_2D"+tag+".root").c_str());
      h_conv->SaveAs((location+"/hist_converted_hist_2D"+tag+".png").c_str());
    }
    c1->SaveAs((location+"/converted_hist_2D"+tag+".root").c_str());
    c1->SaveAs((location+"/converted_hist_2D"+tag+".png").c_str());
    c1->Write();
    c1->Clear();
  }

  // Now just draw the original histogram for completeness
  SetHistogramStyle2D(h,h->GetXaxis()->GetTitle(), h->GetYaxis()->GetTitle(), false);
  c1->SetName("h_orig");
  h->Draw("colz");
  c1->SaveAs((location+"/orig_hist_2D"+tag+".root").c_str());
  c1->SaveAs((location+"/orig_hist_2D"+tag+".png").c_str());
  c1->Write();
  c1->Clear();

  // If we want to draw the slices lines, do so
  if(drawSliceLines){
    c1->SetName("h_slice_lines");
    h->Draw("colz");

    for(double &min : sliceMinX){
      TLine *l = new TLine(min, 0, min, h->GetYaxis()->GetXmax());
      l->SetLineColor(kMagenta-4);
      l->SetLineWidth(3);
      l->SetLineStyle(2);
      l->Draw();
    }
    for(double &max : sliceMaxX){
      TLine *l = new TLine(max, 0, max, h->GetYaxis()->GetXmax());
      l->SetLineColor(kRed-4);
      l->SetLineWidth(3);
      l->SetLineStyle(7);
      l->Draw();
    }
    FormatLatex(h->GetXaxis()->GetXmin()*1.2, h->GetYaxis()->GetXmax()*1.01, "#color[612]{Slice minimum}", 0.04);
    FormatLatex((h->GetXaxis()->GetXmax()-h->GetXaxis()->GetXmin())/5, h->GetYaxis()->GetXmax()*1.01, "#color[628]{Slice maximum}", 0.04);

    c1->SaveAs((location+"/slice_lines_hist"+tag+".root").c_str());
    c1->SaveAs((location+"/slice_lines_hist"+tag+".png").c_str());
    c1->Write();
    //c1->Clear();
  
  }

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
