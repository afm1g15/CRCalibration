#ifndef UTILITIES_H
#define UTILITIES_H

#include "Setup.h"
#include "TSpline.h"
#include "TPolyLine3D.h"

namespace calib{

  typedef std::vector<Plane> PlaneList;

  std::vector<TVector3> wireDirections{
    TVector3(0, 0.812, 0.584),
    TVector3(0,-0.812, 0.584),
    TVector3(0,     1,     0)
  };

  std::vector<int> pal{
    kViolet-5,
    kSpring-5,
    kOrange+5,
    kAzure-5,
    kPink+5,
    kTeal-5,
    kViolet-5,
    kSpring-5,
    kAzure-5,
    kOrange+5,
    kTeal-5,
    kPink+5,
    kViolet-5,
    kSpring-5,
    kAzure-5,
    kOrange+5,
    kTeal-5,
    kPink+5,
    kViolet-5,
    kSpring-5,
    kAzure-5,
    kOrange+5,
    kTeal-5,
    kPink+5,
    kViolet-5,
    kSpring-5,
    kAzure-5,
    kOrange+5,
    kTeal-5,
    kPink+5
  };

  std::vector<int> styles{
    1,
    1,
    1,
    1,
    1,
    1,
    2,
    2,
    2,
    2,
    2,
    2,
    3,
    3,
    3,
    3,
    3,
    3,
    4,
    4,
    4,
    4,
    4,
    4,
    5,
    5,
    5,
    5,
    5,
    5
  };
  /**
   * @brief  The the distance to a defined plane
   *
   * @param  plane 
   * @param  vtx
   * @param  end
   * 
   * @return the distance from the neutrino vertex to the plane
   */
  float GetDistanceToPlane(const Plane &plane, const TVector3 &vtx, const TVector3 &end);

  /**
   * @brief Find the closest plane to the vertex to determine where the track entered
   *
   * @param planes List of planes to check
   * @param vtx Start position of current track
   * @param end End position of current track
   *
   * @return plane Closest plane to vertex
   */
  Plane GetClosestPlane(const PlaneList &planes, const TVector3 &vtx, const TVector3 &end);

  /**
   * @brief  Check if it intersects a given plane
   *
   * @param  plane
   * @param  vtx
   * @param  end
   * @param  length
   *
   * @return true or false
   */
  bool CheckIfIntersectsPlane(const Plane &plane, const TVector3 &vtx, const TVector3 &end, const float &length);

  /**
   * @brief  Check if the projected point is within the bounds of the given plane
   *
   * @param  point
   * @param  plane
   *
   * @return true or false
   */
  bool IsProjectedPointInPlaneBounds(const TVector3 &point, const Plane &plane);

  /**
   * @brief Check if the current plane is in the external list
   *
   * @param geom The geometry to check
   * @param pl Current plane to check
   *
   * @return if external
   */
  bool CheckExternal(const Geometry &geom, const Plane &pl);
      
  /**
   * @brief Get user palette 
   *
   */
  void SetUserPalette();
  
  /**
   * @brief Function to set general style options for 1D histograms
   *
   * @param h       The histogram
   * @param xLabel  X axis label
   * @param yLAbel  Y axis label
   *
   */
  void SetHistogramStyle1D(TH1D *h, const char *xLabel, const char *yLabel);
  
  /**
   * @brief Function to set general style options for 2D histograms
   *
   * @param h       The histogram
   * @param xLabel  X axis label
   * @param yLAbel  Y axis label
   * @param palDefault Whether to use the default palette or not
   *
   */
  void SetHistogramStyle2D(TH2D *h, const char *xLabel, const char *yLabel, const bool &palDefault=true);
  
  /**
   * @brief Function to set general style options for 3D histograms
   *
   * @param h       The histogram
   * @param xLabel  X axis label
   * @param yLAbel  Y axis label
   * @param zLAbel  Z axis label
   *
   */
  void SetHistogramStyle3D(TH3D *h, const char *xLabel, const char *yLabel, const char *zLabel);
  /**
   * @brief Function to set the general style options for a canvas
   *
   * @param c     The canvas
   * @param l     Left margin size
   * @param r     Right margin size
   * @param t     Top margin size
   * @param b     Bottom margin size
   * @param logX  Whether to set the x-axis to a log scale 
   * @param logY  Whether to set the y-axis to a log scale 
   * @param logZ  Whether to set the z-axis to a log scale 
   *
   */
  void SetCanvasStyle(TCanvas *c, const double &l, const double &r, const double &t, const double &b, const bool logX, const bool logY, const bool logZ);

  /**
   * @brief Function to set the format of the latex to be printed
   *
   * @param x  X position
   * @param y  Y position
   * @param s  String to print
   * @param t  Text size, default 0.05
   * @param a  Text alignment, default 11
   *
   */
  void FormatLatex(const double &x, const double &y, const char *s, double t = 0.05, int a = 11);

  /*
   * @brief Function to set the format of the stats for the specified histogram
   *
   * @param h  Histogram to format
   * @param o  OptStat integer
   * @param f  OptFit integer, default 1 
   * @param t  Font, default 132
   *
   */
  void FormatStats(TH1D *h, int o, int f = 1, int t = 132);
  /**
   * @brief Function to write given statistics to a TeX output file
   *
   * @param file Output file
   * @param nFiles Number of files run over
   * @param contents Vector of strings to label the rates
   * @param rates Vector of absolute value rates
   * @param denom Value to divide by for percentages
   *
   */
  void WriteStatsToTeX(ofstream   &file,
                       const int  &nFiles,
                       const std::vector<std::string> &contents,
                       const std::vector<unsigned int> &rates,
                       const double &denom,
                       const std::string &denLab);

  /**
   * @brief Set the x axis of a histogram to be logged
   *
   * @param h Histogram to fix
   */
  void SetLogX(TH1* h);
  
  /**
   * @brief Set the Y axis of a histogram to be logged
   *
   * @param h Histogram to fix
   */
  void SetLogY(TH2* h);

  /**
   * @brief Get the x bin centre which corresponds to the peak of the 1D histogram
   *
   * @param h Histogram to assess
   *
   * @return x value of the peak
   */
  double GetPeakBinCentre(TH1 *h);

  /**
   * @brief Fill the vectors containing information on how to slice the 2D histograms
   *
   * @param h Histogram to slice
   * @param nSlices Number of slices to define
   * @param binWidths Bin widths of the slices to define as a fraction of the total width
   * @param whether to define slices in log space
   * @param minX Vector of minimum x slice limits to fill
   * @param maxX Vector of maximum x slice limits to fill
   * @buffer buffer to apply to either end so slices don't get too close
   */
  void FillSliceVectors(const TH2D *h,
                        const int &nSlices, 
                        const double &binWidths, 
                        const bool &log,
                        std::vector<double> &minX, 
                        std::vector<double> &maxX,
                        double buffer = 80.);

  /**
   * @brief Get labels for the slices that have been defined in TeX format
   *
   * @param minX Vector of Min X values
   * @param maxX Vector of Max X values
   * @param labels Vector of Labels to fill
   * @param units of the slice, if any
   *
   */
  void GetSliceLabelsTeX(const std::vector<double> &minX,
                         const std::vector<double> &maxX,
                         std::vector<std::string> &labels,
                         const std::string &units = "");
  
  /**
   * @brief Get labels for the slices that have been defined
   *
   * @param minX Vector of Min X values
   * @param maxX Vector of Max X values
   * @param labels Vector of Labels to fill
   *
   */
  void GetSliceLabels(const std::vector<double> &minX,
                      const std::vector<double> &maxX,
                      std::vector<std::string> &labels);
  /**
   * @brief Define the histograms for every pre-defined slice
   *
   * @param labels  Histogram label vector
   * @param minY  Vector of minimum Y values
   * @param maxY  Vector of maximum Y values
   * @param hists  Vector of histograms to fill
   *
   */
  void DefineHistograms(const std::vector<std::string> &labels, 
                        const double &minY,
                        const double &maxY,
                        std::vector<TH1D*> hists);
  
  /*
   * @brief Fill the slice histograms from the 2D and bin definitions
   *
   * @param h  2D histogram to slice
   * @param min  Vector of min x values to define the histograms with
   * @param max  Vector of max x values to define the histograms with
   * @param labels  Histogram label vector
   * @param hists  1D histograms to fill
   *
   */
  void FillHistograms(TH2D* h, 
                      const std::vector<double> &min, 
                      const std::vector<double> &max, 
                      const std::vector<std::string> &labels,
                      std::vector<TH1D*> &hists);
  
  /*
   * @brief Determine which bins to merge in the slices
   *
   * @param h  2D histogram to slice
   * @param minX  min x value to define the bins with
   * @param maxX  max x value to define the bins with
   * @param minBin return value of the min bin
   * @param maxBin return calue of the max bin
   *
   */
  void GetBinsToMerge(TH2D *h, 
                      const double &minX, 
                      const double &maxX,
                      int &minBin,
                      int &maxBin);

  /**
   * @brief Landau*Gaussian convolution function from ROOT: http://merlot.ijs.si/~matevz/docs/RootHtmlDoc-5.22-00/tutorials/fit/langaus.C.html
   *
   * @param x
   * @param par
   *
   * @return function at x and par
   */
  double langaufun(double *x, double *par);

  /**
   * @brief Get the angle between two points and the wires in the current plane
   *
   * @param i current plane index
   * @param vtx Vertex of the track
   * @param end End of the track
   *
   * @return costheta
   */
  double GetCosTheta(const int &i, const TVector3 &vtx, const TVector3 &end);
  
  /**
   * @brief Get the angle between two points and the wires in the current plane
   *
   * @param i current plane index
   * @param vtx Vertex of the track
   * @param end End of the track
   *
   * @return sintheta
   */
  double GetSinTheta(const int &i, const TVector3 &vtx, const TVector3 &end);
  
  /**
   * @brief Get the angle between the track and the wire plane
   *
   * @param norm Normal to the wire planes
   * @param vtx Vertex of the track
   * @param end End of the track
   *
   * @return costheta
   */
  double GetAngleToAPAs(const TVector3 &norm, const TVector3 &vtx, const TVector3 &end);

  /**
   * @brief Get the pitch of the current hit from the current and following hit location
   *
   * @param plane the wire plane 
   * @param currHitXYZ location of the current hit
   * @param nextHitXYZ location of the adjacent hit
   *
   * @return hit pitch
   */
  double GetHitPitch(const int &plane, const TVector3 &currXYZ, const TVector3 &nextXYZ);
  
  /**
   * @brief Get the angle of the current segment to the drift direction
   *
   * @param currHitXYZ location of the current hit
   * @param nextHitXYZ location of another hit hit
   *
   * @return cosDrift
   */
  double GetCosDrift(const TVector3 &currXYZ, const TVector3 &nextXYZ);

  /**
   * @brief Check if the start and end y positions should be flipped
   *
   * @param startY  start y position
   * @param endY   end y position
   *
   */
  void CheckAndFlip(TVector3 &start, TVector3 &end);


  /**
   * @brief Get the identity of the best plane for the current reco track
   *
   * @param iTrk  Iterator of the current reco track
   * @param evt   Anatree event object
   * @param bP    BestPlane to allocate
   *
   */
  void GetRecoBestPlane(const int &iTrk, const anatree *evt, int &bP, std::vector<int> &hits);

  /**
   * @brief Get the true energy associated to the current reconstructed track
   *
   * @param iTrk  Iterator of the current reco track
   * @param nGeant  number of G4 entries
   * @param evt  anatree event objects
   * @param bP  bestPlane iterator
   *
   * @return true energy
   */
  double GetTrueEnergyAssoc(const int &iTrk, const int &nGeant, const anatree *evt, const int &bP);
      
  /**
   * @brief Check if the current track ID is in the good GEANT list
   *
   * @param trueID true ID of the current reco track
   * @param goodG4 vector of G4 IDs which pass the chosen cuts
   *
   * @return true if the IDs match
   */
  bool CheckTrueIDAssoc(const int &trueID, const std::vector<int> &goodG4);

  /**
   * @brief Fill the vectors regarding the number of hits on each plane
   *
   * @param id  Current G4 track id
   * @param nHits  Number of hits to check
   * @param evt  anatree event objects
   * @param hitAssoc Hit association to plane vector to fill
   * @param hitsOnPlane Number of hits on each plane vector to fill
   *
   */
  void GetNHitsOnPlane(const int &id, const int &nHits, const anatree *evt, std::vector<std::vector<bool>> &hitAssoc, std::vector<int> &hitsOnPlane);

  double densityEffect(double beta, double gamma, double mass);
  double betaGamma(double KE, double mass);
  double Landau_xi(double KE, double pitch, double mass);
  double Get_Wmax(double KE, double mass);
  double meandEdx(double KE, double mass);
  double MPVdEdx(double KE, double pitch, double mass);
  double IntegratedEdx(double mass, double KE0, double KE1, int n = 10000);
  double RangeFromKE(double KE, double mass);
  TSpline3 * Get_sp_range_KE(double mass, int np = 1000, double minke = .01, double maxke = 2e5);

  void DrawCube(TCanvas *c1, double *rmin, double *rmax, int colour, int lineWidth); 
} // calib
#endif
