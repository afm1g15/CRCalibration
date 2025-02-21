#include "Utilities.h"
#include "TSpline.h"

namespace calib{
  
  //------------------------------------------------------------------------------------------ 
 
  float GetDistanceToPlane(const Plane &plane, const TVector3 &vtx, const TVector3 &end){
    // Get the value of the unit vector of the particle dotted with the normal to the plane
    // If this is zero, the particle is parallel so throw an exception and catch it in the main
    // then continue looping through the list of planes
    TVector3 track_direction   = (end - vtx).Unit();
    float direction_from_plane = track_direction.Dot(plane.GetUnitN());

    if(std::abs(direction_from_plane) <= std::numeric_limits<float>::epsilon()) return std::numeric_limits<float>::max();

    /*
    std::cout << " Plane: " << plane.GetLabel() << std::endl;
    std::cout << " A          : (" << plane.GetA().X() << ", " << plane.GetA().Y() << ", " << plane.GetA().Z() << ") " << std::endl;
    std::cout << " B          : (" << plane.GetB().X() << ", " << plane.GetB().Y() << ", " << plane.GetB().Z() << ") " << std::endl;
    std::cout << " V          : (" << plane.GetV().X() << ", " << plane.GetV().Y() << ", " << plane.GetV().Z() << ") " << std::endl;
    std::cout << " vtx        : (" << vtx.X() << ", " << vtx.Y() << ", " << vtx.Z() << ") " << std::endl;
    std::cout << " end        : (" << end.X() << ", " << end.Y() << ", " << end.Z() << ") " << std::endl;
    std::cout << " V - vtx    : (" << (plane.GetV() - vtx).X() << ", " << (plane.GetV() - vtx).Y() << ", " << (plane.GetV() - vtx).Z() << ") " << std::endl;
    std::cout << " 1/Direction:" << 1/direction_from_plane << std::endl;
    std::cout << " N          : (" << plane.GetUnitN().X() << ", " << plane.GetUnitN().Y() << ", " << plane.GetUnitN().Z() << ") " << std::endl;
    */

    return (1/direction_from_plane)*((plane.GetV() - vtx).Dot(plane.GetUnitN()));
  } // Get distance to plane
  
  //------------------------------------------------------------------------------------------ 

  Plane GetClosestPlane(const PlaneList &planes, const TVector3 &vtx, const TVector3 &end){
    float minDist = 999999.;
    unsigned int planeID = 0;
    unsigned int minID = 0;
    for(const Plane &pl : planes){
      float dist = GetDistanceToPlane(pl, vtx, end);
      //std::cout << " Plane: " << pl.GetLabel() << ", dist: " << dist << std::endl;
      if(abs(dist) < minDist){
        minDist = abs(dist);
        minID   = planeID;
      }
      planeID++;
    }
    return planes.at(minID);
  }
  
  //------------------------------------------------------------------------------------------ 

  bool CheckIfIntersectsPlane(const Plane &plane, const TVector3 &vtx, const TVector3 &end, const float &length){
    float d = GetDistanceToPlane(plane, vtx, end);

    if(abs(d - std::numeric_limits<float>::max()) <= std::numeric_limits<float>::epsilon()){
      return false;
    } 

    if(d < 0 || d > length) {
   //   std::cout << " Doesn't intersect, d: " << d << ", length: " << length << std::endl;
      return false;
    }
 //   if(d > length) return false;
    TVector3 track_direction    = (end - vtx).Unit();
    TVector3 intersection_point = vtx + d * track_direction;

    bool intersects = IsProjectedPointInPlaneBounds(intersection_point, plane);
    /*
    if(!intersects){
      std::cout << "Intersection point: (" << intersection_point.X() << ", " << intersection_point.Y() << ", " << intersection_point.Z() << ") " << std::endl;
    }*/
    return intersects;
  }

  //------------------------------------------------------------------------------------------ 

  bool IsProjectedPointInPlaneBounds(const TVector3 &point, const Plane &plane){
    // Check if the point lies within the bound plane
    return (std::abs((point-plane.GetV()).Dot(plane.GetUnitA())) <= plane.GetAlpha() && std::abs((point-plane.GetV()).Dot(plane.GetUnitB())) <= plane.GetBeta());
  }
  
  //------------------------------------------------------------------------------------------ 
 
  bool CheckExternal(const Geometry &geom, const Plane &pl){
    PlaneList extPlanes = geom.GetExternalPlaneList();
    if(extPlanes.size() == 0){
      std::cerr << " Error: No external planes defined" << std::endl;
      return false;
    }
    for(const Plane &ext : extPlanes){
      if(pl.GetLabel() == ext.GetLabel())
        return true;
    }
    return false;
  }

  //------------------------------------------------------------------------------------------ 
  
  void SetCanvasStyle(TCanvas *c, const double &l, const double &r, const double &t, const double &b, const bool logX, const bool logY, const bool logZ){
    c->SetLeftMargin(l);
    c->SetRightMargin(r);
    c->SetTopMargin(t);
    c->SetBottomMargin(b);

    if(logX)
      c->SetLogx();
    if(logY)
      c->SetLogy();
    if(logZ)
      c->SetLogz();
  } // Canvas Style
  
  // --------------------------------------------------------------------------------------------------------------------------------------------------
  
  void SetUserPalette(){
    const int Number = 4;
    double Red[Number]    = {223/255., 137/255., 61/255., 26/255.};
    double Green[Number]  = {214/255., 119/255., 48/255., 21/255.};
    double Blue[Number]   = {234/255., 187/255., 95/255., 41/255.};
    double Length[Number] = {0, .45, .9, 1};
    /*
    double Red[Number]    = {26/255., 61/255., 137/255., 223/255.};
    double Green[Number]  = {21/255., 48/255., 119/255., 214/255.};
    double Blue[Number]   = {41/255., 95/255., 187/255., 234/255.};
    double Length[Number] = {0, .1, .55, 1};
    */
    int nb = 99;
    TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);
  } // SetUserPalette
  
  // --------------------------------------------------------------------------------------------------------------------------------------------------
  
  void SetHistogramStyle1D(TH1D *h, const char *xLabel, const char *yLabel){
    h->GetXaxis()->SetTitle(xLabel);
    h->GetYaxis()->SetTitle(yLabel);
    h->GetXaxis()->SetTitleFont(132);
    h->GetXaxis()->SetLabelFont(132);
    h->GetYaxis()->SetTitleFont(132);
    h->GetYaxis()->SetLabelFont(132);
    h->GetXaxis()->SetTitleSize(0.055);
    h->GetXaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetTitleSize(0.055);
    h->GetYaxis()->SetLabelSize(0.045);
    h->GetXaxis()->SetMaxDigits(3);
    h->GetYaxis()->SetMaxDigits(3);
    h->GetYaxis()->SetTitleOffset(0.9);
    h->SetStats(0);

  } // 1D Histogram Style

  // --------------------------------------------------------------------------------------------------------------------------------------------------
  
  void SetHistogramStyle2D(TH2D *h, const char *xLabel, const char *yLabel, const bool &palDefault){
    h->GetXaxis()->SetTitle(xLabel);
    h->GetYaxis()->SetTitle(yLabel);
    h->GetXaxis()->SetTitleFont(132);
    h->GetXaxis()->SetLabelFont(132);
    h->GetYaxis()->SetTitleFont(132);
    h->GetYaxis()->SetLabelFont(132);
    h->GetZaxis()->SetTitleFont(132);
    h->GetZaxis()->SetLabelFont(132);
    h->GetXaxis()->SetTitleSize(0.055);
    h->GetXaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetTitleSize(0.055);
    h->GetYaxis()->SetLabelSize(0.045);
    h->GetZaxis()->SetLabelSize(0.045);
    h->GetXaxis()->SetMaxDigits(3);
    h->GetYaxis()->SetMaxDigits(2);
    h->GetZaxis()->SetMaxDigits(3);
    h->GetYaxis()->SetTitleOffset(0.9);
    h->SetContour(99);
    h->SetStats(0);
    if(!palDefault){
      SetUserPalette();
    }
  } // 2D Histogram Style
  
  // --------------------------------------------------------------------------------------------------------------------------------------------------
  
  void SetHistogramStyle3D(TH3D *h, const char *xLabel, const char *yLabel, const char *zLabel){
    h->GetXaxis()->SetTitle(xLabel);
    h->GetYaxis()->SetTitle(yLabel);
    h->GetZaxis()->SetTitle(zLabel);
    h->GetXaxis()->SetTitleFont(132);
    h->GetXaxis()->SetLabelFont(132);
    h->GetYaxis()->SetTitleFont(132);
    h->GetYaxis()->SetLabelFont(132);
    h->GetZaxis()->SetTitleFont(132);
    h->GetZaxis()->SetLabelFont(132);
    h->GetXaxis()->SetTitleSize(0.045);
    h->GetXaxis()->SetLabelSize(0.035);
    h->GetYaxis()->SetTitleSize(0.045);
    h->GetYaxis()->SetLabelSize(0.035);
    h->GetZaxis()->SetTitleSize(0.045);
    h->GetZaxis()->SetLabelSize(0.035);
    h->GetXaxis()->SetMaxDigits(3);
    h->GetYaxis()->SetMaxDigits(3);
    h->GetZaxis()->SetMaxDigits(3);
    h->SetContour(99);
    h->SetStats(0);
  } // 2D Histogram Style
  
  // --------------------------------------------------------------------------------------------------------------------------------------------------
  
  void FormatLatex(const double &x, const double &y, const char *s, double t, int a){
    // Setup the Latex object
    TLatex l;
    l.SetTextAlign(a); // Align at bottom
    l.SetTextSize(t); 
    l.SetTextFont(132);
    l.DrawLatex(x,y,s);
  }

  // --------------------------------------------------------------------------------------------------------------------------------------------------

  void FormatStats(TH1D *h, int o, int f, int t){

    // Firstly, turn the stats on for this histogram
    h->SetStats(1);
    gStyle->SetOptStat(o);
    gStyle->SetOptFit(f);

    TPaveStats *st = static_cast<TPaveStats*>(h->FindObject("stats"));
    st->SetBorderSize(0);
    st->SetFillStyle(0);
    st->SetTextFont(t);
    st->SetTextSize(0.035);
    st->SetX1NDC(0.52);
    st->SetY1NDC(0.63);
    st->SetX2NDC(0.92);
    st->SetY2NDC(0.93);
  }
  
  // --------------------------------------------------------------------------------------------------------------------------------------------------

  void WriteStatsToTeX(ofstream   &file,
                       const int  &nFiles,
                       const std::vector<std::string> &contents,
                       const std::vector<unsigned int> &rates,
                       const double &denom,
                       const std::string &denLab){
 
    // Calculate the approximate number of days from the number of files
    // nFiles * 500 events per file / 14113 events per day //This is for full size
    double nDays = (9976)/5065.;   //this is for 1x2x6 from global intensity ( x s in day = N / day)
    // Start by writing the first few lines of the tex file
    file << "\\begin{document} " << std::endl;
    file << "  \\thispagestyle{empty}" << std::endl;
    file << "  \\renewcommand{\\arraystretch}{1.2}" << std::endl;

    // Setup the table
    file << "  \\begin{table}[h!]" << std::endl;
    file << "    \\centering" << std::endl;
    file << "    \\begin{tabular}{ m{4cm} * {2}{ >{\\centering\\arraybackslash}m{4cm} } }" << std::endl;
    file << "      \\toprule" << std::endl;
    file << "      Statistic & Rate / " << std::setprecision(4) << nDays << " Days & \\% "+denLab+" \\\\" << std::endl;
    file << "      \\midrule" << std::endl;

    for(unsigned int nRate = 0; nRate < rates.size(); ++nRate){
      std::string str = contents.at(nRate);
      unsigned int val = rates.at(nRate);
      if(nRate == 0){
        file << "      " << str << " & " << "\\num{" << std::setprecision(4) << val << "} & - \\\\" << std::endl;
        file << "      \\midrule" << std::endl;
      }
      else
        file << "      " <<  str << " & " << "\\num{" << std::setprecision(4) << val << "} & " << std::setprecision(5) << (val/denom)*100 << "~\\%  \\\\" << std::endl;
    }
    file << "      \\bottomrule" << std::endl;
    file << "    \\end{tabular}" << std::endl;
    file << "  \\end{table}" << std::endl;
    file << "\\end{document}" << std::endl;
  }

  // --------------------------------------------------------------------------------------------------------------------------------------------------

  void SetLogX(TH1* h){
    TAxis* axis = h->GetXaxis();

    double start = TMath::Log10(axis->GetXmin());
    double stop = TMath::Log10(axis->GetXmax());
    double range = stop - start;
    int nbins = axis->GetNbins();
    double binwidth = range / nbins;

    double *bins = new double[nbins+1];
    for (int i = 0; i < (nbins+1); ++i) {
      bins[i] = TMath::Power(10, start + i*binwidth);
    }
    axis->Set(nbins, bins);
    delete[] bins;
  } // Set Log X
  
  // --------------------------------------------------------------------------------------------------------------------------------------------------

  void SetLogY(TH2* h){
    TAxis* axis = h->GetYaxis();

    double start = TMath::Log10(axis->GetXmin());
    double stop = TMath::Log10(axis->GetXmax());
    double range = stop - start;
    int nbins = axis->GetNbins();
    double binwidth = range / nbins;

    double *bins = new double[nbins+1];
    for (int i = 0; i < (nbins+1); ++i) {
      bins[i] = TMath::Power(10, start + i*binwidth);
    }
    axis->Set(nbins, bins);
    delete[] bins;
  } // Set Log X
  
  // --------------------------------------------------------------------------------------------------------------------------------------------------

  double GetPeakBinCentre(TH1 *h){
    // Get the bin containing the maximum content
    int bin = h->GetMaximumBin();
    
    // Calculated this from the low edge and the width in case of log scales
    double low   = h->GetBinLowEdge(bin);
    double width = h->GetBinWidth(bin);

    // Now calculate and return the centre of this bin
    return low+(width/2.);
  }
  
  // --------------------------------------------------------------------------------------------------------------------------------------------------

  void FillSliceVectors(const TH2D *h,
                        const int &nSlices, 
                        const double &binWidths, 
                        const bool &log,
                        std::vector<double> &minX, 
                        std::vector<double> &maxX,
                        double buffer){
    
    // Get the full range and the 'buffered' range
    double buffMin   = h->GetXaxis()->GetXmin()+buffer;
    double buffMax   = h->GetXaxis()->GetXmax()-buffer;
    double buffRange = buffMax-buffMin;
    double width     = binWidths*buffRange;
    double buffStep  = (buffRange-width)/static_cast<double>(nSlices-1);
    
    if(log){
      // If we are dealing with an input histogram with logX values, need to translate before defining the bins
      buffMin   = TMath::Log10(buffMin);
      buffMax   = TMath::Log10(buffMax);
      buffRange = buffMax-buffMin;
      width     = buffRange*binWidths;
      buffStep  = (buffRange-width)/static_cast<double>(nSlices-1);

      // Determine the central location for each new bin
      // Set the buffer range min and max to be the bin and max of the first and last new bin
      minX.push_back(pow(10,buffMin));
      maxX.push_back(pow(10,buffMin+width));
      for(int n = 1; n <= nSlices-2; ++n){
        // Find the min and max bins by adding the step to the first min and max
        double binCentre = buffMin+n*buffStep + width*0.5; 
        minX.push_back(pow(10,binCentre-0.5*width));
        maxX.push_back(pow(10,binCentre+0.5*width));
      }
      minX.push_back(pow(10,buffMax-width));
      maxX.push_back(pow(10,buffMax));
    }
    else{
      // Determine the central location for each new bin
      // Set the buffer range min and max to be the bin and max of the first and last new bin
      minX.push_back(buffMin);
      maxX.push_back(buffMin+width);
      for(int n = 1; n <= nSlices-2; ++n){
        // Find the min and max bins by adding the step to the first min and max
        double binCentre = buffMin+n*buffStep + width*0.5; 
        minX.push_back(binCentre-0.5*width);
        maxX.push_back(binCentre+0.5*width);
      }
      minX.push_back(buffMax-width);
      maxX.push_back(buffMax);
    }

    if((minX.size() + maxX.size())*0.5 != nSlices){
      std::cerr << " Error: The number of slices created does not match the number desired, " << std::endl;
      std::cerr << " Desired: " << nSlices << ", created: " << (minX.size() + maxX.size())*0.5 << std::endl;
    }
      
    std::cout << " Buff Range: " << buffRange << ", bin widths: " << width << ", buffStep: " << buffStep << std::endl;
    std::cout << std::setw(15) << "Defined bins: " << std::setw(10) << " Min" << std::setw(5) << " | " << std::setw(10) << " Max" << std::endl;
    for(int n = 0; n < nSlices; ++n){
    std::cout << std::setw(15) << " " << std::setw(10) << minX.at(n) << std::setw(5) << " | " << std::setw(10) << maxX.at(n) << std::endl;
    }
    return;
  }
  
  // --------------------------------------------------------------------------------------------------------------------------------------------------

  void GetSliceLabels(const std::vector<double> &minX,
                      const std::vector<double> &maxX,
                      std::vector<std::string> &labels){
    unsigned int n = minX.size();

    std::vector< std::vector<double> > minMax{minX,maxX};

    // Loop over the vectors
    for(unsigned int i = 0; i < n; ++i){
      std::vector<std::string> minMaxStr;
      // Now do the same thing for min and max
      for(unsigned int j = 0; j < minMax.size(); ++j){

        // Convert doubles to strings with single-point precision
        // Create an output string stream
        std::ostringstream ss;
        // Set Fixed -Point Notation
        ss << std::fixed;
        // Set precision to 2 digits
        ss << std::setprecision(2);
        // Add double to stream and get string
        ss << minMax.at(j).at(i);
        std::string str = ss.str();

        // Now translate '-' to 'm' and '.' to 'p'
        TString tStr(str);
        tStr.ReplaceAll("-","m");
        tStr.ReplaceAll(".","p");
        minMaxStr.push_back(tStr.Data());

      } // Loop over minmax
      std::string label = "slice_"+minMaxStr.at(0)+"_to_"+minMaxStr.at(1);
      labels.push_back(label);
    } // Loop over slices
    return;
  } // GetSliceLabels
  
  // --------------------------------------------------------------------------------------------------------------------------------------------------

  void GetSliceLabelsTeX(const std::vector<double> &minX,
                         const std::vector<double> &maxX,
                         std::vector<std::string> &labels,
                         const std::string &units){
    unsigned int n = minX.size();
  
    std::vector< std::vector<double> > minMax{minX,maxX};

    // Loop over the vectors
    for(unsigned int i = 0; i < n; ++i){
      std::vector<std::string> minMaxStr;
      // Now do the same thing for min and max
      for(unsigned int j = 0; j < minMax.size(); ++j){

        // Convert doubles to strings with single-point precision
        // Create an output string stream
        std::ostringstream ss;
        // Set Fixed -Point Notation
        ss << std::fixed;
        // Set precision to 2 digits
        ss << std::setprecision(2);
        // Add double to stream and get string
        ss << minMax.at(j).at(i);
        std::string str = ss.str();

        // Now translate '-' to 'm' and '.' to 'p'
        TString tStr(str);
        minMaxStr.push_back(tStr.Data());

      } // Loop over minmax
      std::string label = minMaxStr.at(0)+" to "+minMaxStr.at(1)+" "+units;
      labels.push_back(label);
    } // Loop over slices
    return;
  }
  
  // --------------------------------------------------------------------------------------------------------------------------------------------------

  void DefineHistograms(const std::vector<std::string> &labels, 
                        const double &minY,
                        const double &maxY,
                        std::vector<TH1D*> hists){

    for(unsigned int i = 0; i < labels.size(); ++i){
      hists.push_back(new TH1D(labels.at(i).c_str(),"",100,minY,maxY));
    } // Loop
    return;
  } // DefineHistograms
  
  // --------------------------------------------------------------------------------------------------------------------------------------------------
  
  void FillHistograms(TH2D* h, 
                      const std::vector<double> &min, 
                      const std::vector<double> &max, 
                      const std::vector<std::string> &labels, 
                      std::vector<TH1D*> &hists){
  
    for(unsigned int i = 0; i < min.size(); ++i){
      double minX = min.at(i);
      double maxX = max.at(i);
      int minBin = 999;
      int maxBin = -999;
      // Find the bins to merge based on the min and max requirements and the existing bin centres
      GetBinsToMerge(h,minX,maxX,minBin,maxBin);

      // Now get the projection
      hists.push_back(h->ProjectionY(labels.at(i).c_str(),minBin,maxBin));

    }
    return;
  } // FillHistograms
  
  // --------------------------------------------------------------------------------------------------------------------------------------------------
  
  void GetBinsToMerge(TH2D *h, 
                      const double &minX, 
                      const double &maxX,
                      int &minBin,
                      int &maxBin){

    // Using the 2D histogram and the chosen bin edges, find the id of the min and max bins
    for(int n = 1; n <= h->GetNbinsX(); ++n){
      double centre = h->GetXaxis()->GetBinCenter(n);
      if(centre >= minX && centre <= maxX){
        // Reset the max and min bin integers if needed
        if(n < minBin)
          minBin = n;
        if(n > maxBin)
          maxBin = n;
      }
    }
    return;
  }

  // --------------------------------------------------------------------------------------------------------------------------------------------------
  
  double langaufun(double *x, double *par) {

    //Fit parameters:
    //par[0]=Width (scale) parameter of Landau density
    //par[1]=Most Probable (MP, location) parameter of Landau density
    //par[2]=Total area (integral -inf to inf, normalization constant)
    //par[3]=Width (sigma) of convoluted Gaussian function
    //
    //In the Landau distribution (represented by the CERNLIB approximation), 
    //the maximum is located at x=-0.22278298 with the location parameter=0.
    //This shift is corrected within this function, so that the actual
    //maximum is identical to the MP parameter.

    // Numeric constants
    double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
    double mpshift  = -0.22278298;       // Landau maximum location

    // Control constants
    double np = 100.0;      // number of convolution steps
    double sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

    // Variables
    double xx;
    double mpc;
    double fland;
    double sum = 0.0;
    double xlow,xupp;
    double step;
    double i;

    // MP shift correction
    mpc = par[1] - mpshift * par[0]; 

    // Range of convolution integral
    xlow = x[0] - sc * par[3];
    xupp = x[0] + sc * par[3];

    step = (xupp-xlow) / np;

    // Convolution integral of Landau and Gaussian by sum
    for(i=1.0; i<=np/2; i++) {
      xx = xlow + (i-.5) * step;
      fland = TMath::Landau(xx,mpc,par[0]) / par[0];
      sum += fland * TMath::Gaus(x[0],xx,par[3]);

      xx = xupp - (i-.5) * step;
      fland = TMath::Landau(xx,mpc,par[0]) / par[0];
      sum += fland * TMath::Gaus(x[0],xx,par[3]);
    }

    return (par[2] * step * sum * invsq2pi / par[3]);
  }

  // --------------------------------------------------------------------------------------------------------------------------------------------------
  
  double GetCosTheta(const int &i, const TVector3 &vtx, const TVector3 &end){
  
    // First, get the direction of the wires on the plane
    TVector3 wireDir = wireDirections.at(i);

    // Now get the direction from start and end 
    TVector3 dir = (end-vtx);
    dir *= 1/static_cast<double>(end.Mag()*vtx.Mag());

    // Now get the angle between them
    double costheta = wireDir.Dot(dir)/static_cast<double>(wireDir.Mag()*dir.Mag());

    return costheta;
  }

  // --------------------------------------------------------------------------------------------------------------------------------------------------
  
  double GetSinTheta(const int &i, const TVector3 &vtx, const TVector3 &end){
  
    // First, get the direction of the wires on the plane
    TVector3 wireDir = wireDirections.at(i);

    // Now get the direction from start and end 
    TVector3 dir = (end-vtx);
    dir *= 1/static_cast<double>(end.Mag()*vtx.Mag());

    // Now get the angle between them
    double sintheta = (wireDir.Cross(dir)).Mag()/static_cast<double>(wireDir.Mag()*dir.Mag());

    return sintheta;
  }

  // --------------------------------------------------------------------------------------------------------------------------------------------------
  
  double GetAngleToAPAs(const TVector3 &norm, const TVector3 &vtx, const TVector3 &end){

    // Get the unit direction of the track 
    TVector3 trackDir = (end-vtx);
    trackDir *= 1/static_cast<double>(end.Mag()*vtx.Mag());

    // Get the angle of the track to the wire planes
    double costoplane = norm.Dot(trackDir)/static_cast<double>(norm.Mag()*trackDir.Mag());

    return costoplane;

  }

  // --------------------------------------------------------------------------------------------------------------------------------------------------
  
  double GetHitPitch(const int &plane, const TVector3 &currXYZ, const TVector3 &nextXYZ){

    // pitch = 0.48/sintheta
    double sinTheta = GetSinTheta(plane,currXYZ,nextXYZ);

    // Get the pitch
    double hitPitch = 0.48/sinTheta;

    return hitPitch;
  }
  
  // --------------------------------------------------------------------------------------------------------------------------------------------------
  
  double GetCosDrift(const TVector3 &currXYZ, const TVector3 &nextXYZ){
    
    // Get the angle of the hit to the x axis
    TVector3 xDir(1,0,0);
    TVector3 hitVec = (nextXYZ-currXYZ);

    double cosDrift = hitVec.Dot(xDir)/static_cast<double>(hitVec.Mag()*xDir.Mag());
    return cosDrift;
  }
  
  // --------------------------------------------------------------------------------------------------------------------------------------------------
  
  void CheckAndFlip(TVector3 &start, TVector3 &end){
      if(start.Y() < end.Y()){
        TVector3 temp(end);
        end   = start;
        start = temp;
      }
  }
  
  // --------------------------------------------------------------------------------------------------------------------------------------------------
  
  void GetRecoBestPlane(const int &iTrk, const anatree *evt, int &bP, std::vector<int> &hits){
    // Get the best plane
    int currHits  = -999;
    for(int iPlane = 0; iPlane < 3; ++iPlane){
      hits.at(iPlane) = evt->ntrkhits_pandoraTrack[iTrk][iPlane];
      if(evt->ntrkhits_pandoraTrack[iTrk][iPlane] > currHits){
        currHits = evt->ntrkhits_pandoraTrack[iTrk][iPlane];
        bP       = iPlane; 
      } // CurrHits
    } // Planes
  }
  
  // --------------------------------------------------------------------------------------------------------------------------------------------------
 
  double GetTrueEnergyAssoc(const int &iTrk, const int &nGeant, const anatree *evt, const int &bP){
    double eng = -1;
    for(int iG4 = 0; iG4 < nGeant; ++iG4){
      int trueID = evt->TrackId[iG4];

      if(evt->trkidtruth_pandoraTrack[iTrk][bP] == trueID){
        eng = evt->Eng[iG4];
        break;
      } // ID if
    } // G4
    return eng;
  }
  
  // --------------------------------------------------------------------------------------------------------------------------------------------------
 
  bool CheckTrueIDAssoc(const int &trueID, const std::vector<int> &goodG4){
    for(const int &id: goodG4){
      if(trueID == id)
        return true;
    }
    return false;
  }
  
  // --------------------------------------------------------------------------------------------------------------------------------------------------
 
  void GetNHitsOnPlane(const int &id, const int &nHits, const anatree *evt, std::vector<std::vector<bool>> &hitAssoc, std::vector<int> &hitsOnPlane){
    for(int iPlane = 0; iPlane < 3; ++iPlane){
      for(int iHit = 0; iHit < nHits; ++iHit){

        // Get the ID of the hit
        int hitId  = evt->hit_trkid[iHit];

        // If we are not looking at the current G4 track, continue
        bool currentG4 = false;
        for(int iTrk = 0; iTrk < evt->ntracks_pandoraTrack; ++iTrk){
          int recoId = evt->trkId_pandoraTrack[iTrk];

          // Check if the current hit is in the reco track
          if(recoId != hitId) continue;

          // If it is, check if the reco track is the G4 track
          if(evt->trkidtruth_pandoraTrack[iTrk][iPlane] == id){
            currentG4 = true;
            break;
          } // ID check
        } // iTrk
        if(!currentG4) continue;

        // Now fill the vectors
        if(evt->hit_plane[iHit] == iPlane){
          hitsOnPlane.at(iPlane)++;
          hitAssoc.at(iPlane).at(iHit) = true;
        } // Check current plane
      } // Hits
    } // Planes
  }
  
  // --------------------------------------------------------------------------------------------------------------------------------------------------

  // == Bethe-Bloch parameters, https://indico.fnal.gov/event/14933/contributions/28526/attachments/17961/22583/Final_SIST_Paper.pdf

  const double rho = 1.39; // [g/cm3], density of LAr
  const double K = 0.307075; // [MeV cm2 / mol]
  const double Z = 18.; // atomic number of Ar
  const double A = 39.948; // [g / mol], atomic mass of Ar
  const double I = 188.0e-6; // [MeV], mean excitation energy
  const double me = 0.511; // [Mev], mass of electron
  // == Parameters for the density correction
  const double density_C = 5.2146;
  const double density_y0 = 0.2;
  const double density_y1 = 3.0;
  const double density_a = 0.19559;
  const double density_k = 3.0;

  double densityEffect(double beta, double gamma, double mass){

    double density_y = TMath::Log10(beta * gamma);
    double ln10 = TMath::Log(10);
    double this_delta = 0.;
    if(density_y > density_y1){
      this_delta = 2.0 * ln10 * density_y - density_C;
    }
    else if (density_y < density_y0){
      this_delta = 0.;
    }
    else{
      this_delta = 2.0 * ln10 * density_y - density_C + density_a * pow(density_y1 - density_y, density_k);
    }
  
    return this_delta;
  }

  double betaGamma(double KE, double mass){

    double gamma, beta;
    gamma = (KE + mass) / mass;
    beta = sqrt( 1 - 1/pow(gamma,2));

    return beta*gamma;
  
  }

  double Landau_xi(double KE, double pitch, double mass){
    double gamma = (KE/mass)+1.0;
    double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
    double xi = rho * pitch * 0.5 * K * (Z / A) * pow(1. / beta, 2);
    return xi;
  }  

  double Get_Wmax(double KE, double mass){
    double gamma = (KE/mass)+1.0;
    double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
    double Wmax = (2.0 * me * pow(beta * gamma, 2)) / (1.0 + 2.0 * me * (gamma / mass) + pow((me / mass),2));

    return Wmax;
  }

  double meandEdx(double KE, double mass){

    double gamma = (KE + mass) / mass;
    double beta = sqrt( 1 - 1/pow(gamma,2));
    double wmax = Get_Wmax(KE, mass);
    double dEdX = (rho*K*Z*pow(1.,2))/(A*pow(beta,2))*(0.5*log(2*me*pow(gamma,2)*pow(beta,2)*wmax/pow(I,2)) - pow(beta,2) - densityEffect( beta, gamma, mass )/2 ); // charge = 1

    return dEdX;
  }

  double MPVdEdx(double KE, double pitch, double mass){

    double gamma = (KE + mass) / mass;
    double beta = sqrt( 1 - 1/pow(gamma,2));

    double xi = Landau_xi(KE, pitch, mass);

    double eloss_mpv = xi*(log( 2*me*pow(gamma,2)*pow(beta,2) / I ) + log( xi / I ) + 0.2 - pow(beta,2) - densityEffect( beta, gamma, mass ) )/pitch;

    return eloss_mpv;
  }

  double IntegratedEdx(double mass, double KE0, double KE1, int n){

    if (KE0>KE1) swap(KE0, KE1);

    double step = (KE1-KE0)/n;

    double area = 0;

    for (int i = 0; i<n; ++i){
      double dEdx = meandEdx(KE0 + (i+0.5)*step, mass);
      if (dEdx)
        area += 1/dEdx*step;
    }
    return area;
  }


  double RangeFromKE(double KE, double mass){

    return IntegratedEdx(mass, 0, KE);
  }


  TSpline3 * Get_sp_range_KE(double mass, int np, double minke, double maxke){

    double *KE = new double[np];
    double *Range = new double[np];

    for (int i = 0; i<np; ++i){
      double ke = pow(10, log10(minke)+i*log10(maxke/minke)/np);
      KE[i] = ke;
      Range[i] = RangeFromKE(ke, mass);
    }

    TSpline3 *out = new TSpline3("sp_range_KE", Range, KE, np, "b2e2", 0, 0);
    delete[] KE;
    delete[] Range;

    return out;
  }

  //--------------------------------------------------

  void DrawCube(TCanvas *c1, double *rmin, double *rmax, int colour, int lineWidth)
  {
    c1->cd();
    TList *outline = new TList;
    TPolyLine3D *p1 = new TPolyLine3D(4);
    TPolyLine3D *p2 = new TPolyLine3D(4);
    TPolyLine3D *p3 = new TPolyLine3D(4);
    TPolyLine3D *p4 = new TPolyLine3D(4);
    p1->SetLineColor(colour);

    if(lineWidth == -1)
      p1->SetLineWidth(3.);
    else
      p1->SetLineWidth(3.);

    p1->Copy(*p2);
    p1->Copy(*p3);
    p1->Copy(*p4);
    outline->Add(p1);
    outline->Add(p2);
    outline->Add(p3);
    outline->Add(p4); 
    TPolyLine3D::DrawOutlineCube(outline, rmin, rmax);
    p1->Draw();
    p2->Draw();
    p3->Draw();
    p4->Draw();
  }




} // calib
