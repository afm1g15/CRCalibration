#include "SliceHelpers.h"

namespace calib{
  
  //------------------------------------------------------------------------------------------ 
 
  void FillSliceVectorsEqualRates(TH2D *h,
                                  int &nSlices, 
                                  const bool &log,
                                  std::vector<double> &minX, 
                                  std::vector<double> &maxX,
                                  double buffer){
  
    // First, get the total number of entries in the histogram and estimate the number/slice
    // Make sure to check within the buffer range
    double buffMin   = h->GetXaxis()->FindBin(h->GetXaxis()->GetXmin()+buffer);
    double buffMax   = h->GetXaxis()->FindBin(h->GetXaxis()->GetXmax()-buffer);
    TH1D *hProjBuff  = h->ProjectionY("buffRange",buffMin,buffMax);
    std::cout << " Buffer bin range: " << buffMin << ", " << buffMax << std::endl;

    double nEntries = hProjBuff->GetEntries();
    double nPerSlice = nEntries/static_cast<double>(nSlices);
    std::cout << " Number of entries: " << nEntries << ", number of slices: " << nSlices << ", n/slice: " << nPerSlice << std::endl;

    // Now loop over the x bins (from any buffer) and accumulate the events
    int nAccum = 0;
    bool maxReached = false;
    for(int i = 1; i <= h->GetNbinsX(); ++i){
      double buffCenterMin = h->GetXaxis()->GetBinCenter(buffMin);
      double buffCenterMax = h->GetXaxis()->GetBinCenter(buffMax);
      std::cout << " Bin: " << i << " at " << h->GetXaxis()->GetBinCenter(i);
      if(h->GetXaxis()->GetBinCenter(i) < buffCenterMin){
        std::cout << " Below the buffer" << std::endl;
        continue;
      }
      else if(h->GetXaxis()->GetBinCenter(i) > buffCenterMax) {
        std::cout << " Reached the last bin in the buffer range, pushing back " << h->GetXaxis()->GetBinLowEdge(i) << " to maxX list" << std::endl;
        maxX.push_back(h->GetXaxis()->GetBinLowEdge(i));
        break;
      }
      else{
        TH1D *hCurrProj = h->ProjectionY(("currProj"+std::to_string(i)).c_str(),i,i);

        if(hCurrProj->GetEntries() == 0) continue;
        if(nAccum == 0){
          std::cout << " Pushing back " << h->GetXaxis()->GetBinLowEdge(i) << " to minX list" << std::endl;
          minX.push_back(h->GetXaxis()->GetBinLowEdge(i));
        }
        nAccum += hCurrProj->GetEntries();
        std::cout << " has " << hCurrProj->GetEntries() << " entries, and the current total is: " << nAccum << std::endl;
        if(nAccum > nPerSlice){ // If the accumulated number of entries is within the desired range, set the edges
          std::cout << " Pushing back " << h->GetXaxis()->GetBinUpEdge(i) << " to maxX list" << std::endl;
          maxX.push_back(h->GetXaxis()->GetBinUpEdge(i));
          nAccum = 0;
        } // If nAccum is within the range
        else if(nAccum < nPerSlice && i == h->GetNbinsX()){
          std::cout << " Reached the last bin pushing back " << h->GetXaxis()->GetBinUpEdge(i) << " to maxX list" << std::endl;
          maxX.push_back(h->GetXaxis()->GetBinUpEdge(i));
          maxReached = true;
          break;
        }
      } // Between buffMin and buffMax
    } // X bins
    if((minX.size() + maxX.size())*0.5 != nSlices && !maxReached){
      std::cerr << " The number of slices created does not match the number desired, " << std::endl;
      std::cerr << " Desired: " << nSlices << ", created: " << (minX.size() + maxX.size())*0.5 << std::endl;
      nSlices = (minX.size() + maxX.size())*0.5;
    }
    else if(maxReached && (minX.size() + maxX.size())*0.5 != nSlices){
      nSlices = (minX.size() + maxX.size())*0.5;
    }
     
    std::cout << " Final number of slices: " << nSlices << std::endl; 
    std::cout << std::setw(15) << "Defined bins: " << std::setw(10) << " Min" << std::setw(5) << " | " << std::setw(10) << " Max" << std::endl;
    for(int n = 0; n < nSlices; ++n){
    std::cout << std::setw(15) << " " << std::setw(10) << minX.at(n) << std::setw(5) << " | " << std::setw(10) << maxX.at(n) << std::endl;
    }
    return;
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
  
} // calib
