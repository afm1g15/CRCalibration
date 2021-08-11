#ifndef UTILITIES_H
#define UTILITIES_H

#include "Setup.h"

namespace calib{

  typedef std::vector<Plane> PlaneList;
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
   *
   */
  void SetHistogramStyle2D(TH2D *h, const char *xLabel, const char *yLabel);
  
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
} // calib
#endif
