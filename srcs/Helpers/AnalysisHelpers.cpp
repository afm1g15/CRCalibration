#include "AnalysisHelpers.h"

namespace calib{
  
  // --------------------------------------------------------------------------------------------------------------------------------------------------
  
  void CheckAndFlip(TVector3 &start, TVector3 &end){
      if(start.Y() < end.Y()){
        TVector3 temp(end);
        end   = start;
        start = temp;
      }
  }
 
  //------------------------------------------------------------------------------------------ 
 
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
  
  void GetFWHMFromTF1(TH1 *h, TF1 *f, double &fwhm, double &fitMaxX){

    // Get the x value of the maximum to return and get the maximum on the y axis from the histogram itself to half it
    fitMaxX = f->GetMaximumX();
    double fitMax  = h->GetMaximum();
    double halfMax = fitMax/2.;

    // Now find the full width half max
    // First, take the range of the xaxis to determine an appropriate step size (0.1%)
    double xMin     = f->GetXmin();
    double xMax     = f->GetXmax();
    double x0       = fitMaxX;
    double stepSize = std::abs(xMax-xMin)*0.001;

    // Then, starting from the maximum value (x0), evaluate the function by stepping left and right 
    // to find the point where it approaches halfMax
    //
    // Do this by calculating diff = (hcurr - h/2) for each step
    // When xL < x < x0 and x0 < x < xR, diff will be positive
    // Find the point when either diff approaches 0 or when diff becomes negative
    // Calculate the halway point between xprev and xcurr to extract the left and right x positions
    std::vector<double> hPos(2, f->Eval(x0));
    std::vector<double> xPos(2,x0);
    std::vector<double> xHalf(2,x0);
    std::vector<double> steps{-stepSize,stepSize};
    double &xL = xHalf.at(0);
    double &xR = xHalf.at(1);

    // Loop for checking left and right
    for(unsigned int i = 0; i < 2; ++i){ 
      // Access the left or right value
      double hCurr = hPos.at(i);
      double xCurr = xPos.at(i);
      double xPrev = xCurr;
      double step  = steps.at(i);
      double diff  = hCurr - halfMax;
     
      // Keep an eye on the number of iterations
      unsigned int it = 0;
      // Loop over the number of steps
      while(diff > 0){
    
        // Access and update the current values
        xPrev  = xCurr;
        xCurr += step;

        hCurr  = f->Eval(xCurr);
        diff   = hCurr - halfMax;
       
        /*
        // Print test
        std::cout << " It: " << it
                  << ", hHalf: " << halfMax
                  << ", x0: " << x0
                  << ", hCurr: " << hCurr
                  << ", xCurr: " << xCurr
                  << ", xPrev: " << xPrev
                  << ", step: " << step
                  << ", diff: " << diff << std::endl;
        //std::cin.get();
        */
        // If diff is negative, find the halfway position between the current and previous x positions
        if(diff < 0){ // If negative
          xHalf.at(i) = (xPrev+xCurr)/2.;
          //std::cout << " xHalf at " << i << ": " << xHalf.at(i) << std::endl;
          break;
        }
        // If diff is 0, set xCurr to be the halfway position 
        else if(diff < std::numeric_limits<double>::epsilon()){
          xHalf.at(i) = xCurr;
          //std::cout << " xHalf at " << i << ": " << xHalf.at(i) << std::endl;
          break;
        }
        
        ++it;

      } // While diff > 0
    } // i
    fwhm = std::abs(xR - xL);
  } // GetFWHM
  
  // --------------------------------------------------------------------------------------------------------------------------------------------------
  
  void GetCoefficient(const double &mpv, const double &mp, const double &gsig, double &psi){

    // MPV = MP + psi*GSigma
    // ->
    // psi = (MPV-MP)/GSigma
    psi = (mpv-mp)/gsig;
  }
  
  // --------------------------------------------------------------------------------------------------------------------------------------------------
  
  bool GetMPVUncertainty(const double &mpv, 
                         const double &mp, 
                         const double &mp_error, 
                         const double &gsig, 
                         const double &gsig_error, 
                         const double &psi,
                         double &mpv_error){

    // Make sure the coefficient has been calculated correctly
    double check_mpv = mp + psi*gsig;
    if(abs(check_mpv - mpv) > std::numeric_limits<double>::epsilon()){
      std::cerr << " Error: Psi has not been calculated correctly since mpv: " << mpv << ", does not equal mp + psi*gsigma: " << check_mpv << std::endl;
      return false;
    }

    // Now calculate the uncertainty on MPV using the differential error propagation formula
    //
    // For this linear combination of paramters, this becomes:
    //
    // MPV = phi*MP + psi*GSigma
    // sigMPV^2 = phi^2*sigMP^2 + psi^2*sigGSig^2
    //
    // since phi = 1, it's even more simple:
    //
    double sigMPVSq = pow(mp_error,2) + 
                      pow(psi,2)*pow(gsig_error,2);
    mpv_error       = sqrt(sigMPVSq);    

    return true;
  }
  
} // calib
