/*
 * @brief Helper functions for plotting
 *
 */
#ifndef PLOTTINGHELPERS_H
#define PLOTTINGHELPERS_H

#include "../Setup/Setup.h"

namespace calib{

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
   * @brief Get user palette 
   *
   */
  void SetUserPalette();
  
  /**
   * @brief Function to set statistical uncertainties in 1D histogram bins
   *
   * @param h           The histogram
   * @param widthScale  Whether or not we are scaling this histogram by bin width, default [false]
   *
   */
  void SetHistogramStatErrors(TH1D *h, const bool &widthScale = false);

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
  
  /**
   * @brief Function to set the format of the latex to be printed using NDC coordinates
   *
   * @param x  X position
   * @param y  Y position
   * @param s  String to print
   * @param t  Text size, default 0.05
   * @param a  Text alignment, default 11
   *
   */
  void FormatLatexNDC(const double &x, const double &y, const char *s, double t, int a);

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
   * @brief Function to write given statistics to a TeX output file from multiple productions for comparison
   *
   * @param file Output file
   * @param nFiles Number of events run over
   * @param contents Map containing the parameters to print and the vector of counters for each production
   * @param rates Vector of labels for each production
   * @param verbose Boolean for print level
   *
   */
  void WriteStatsToTeXMultiProd(ofstream &file,
                                const int  &nEvents,
                                std::map<std::string,std::vector<unsigned int>> &contentMap,
                                std::vector<std::string> &counterLabels,
                                std::vector<std::string> prodLabels,
                                const bool &verbose);
  
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

} // calib
#endif
