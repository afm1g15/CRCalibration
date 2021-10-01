#ifndef UTILITIES_H
#define UTILITIES_H

#include "Setup.h"

namespace calib{

  typedef std::vector<Plane> PlaneList;

  std::vector<int> pal{
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
   *
   */
  void FormatLatex(const double &x, const double &y, const char *s);

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
   * @param binWidths Bin widths of the slices to define
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
} // calib
#endif
