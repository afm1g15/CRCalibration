#ifndef SETUP
#define SETUP

#include "TROOT.h"
#include "TLine.h"
#include "TLatex.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TVector3.h"
#include "TString.h"

#include "anatree.h"

#include <numeric>
#include <time.h>
#include <ctime>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <string>

namespace calib{
  /**
   * @brief  Get the time to monitor the efficiency of the programme
   *
   * @param  rawtime time object to read
   *
   */
  void GetTime(time_t &rawtime);

  /**
   * @brief  Get the total time taken and print nicely
   *
   * @param  starttime start time object to read
   * @param  endtime end time object to read
   *
   */
  void GetTotalTime(time_t &starttime, time_t &endtime);

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
   * @brief Function to set the general style options for a 2D canvas
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
  void SetCanvasStyle2D(TCanvas *c, const double &l, const double &r, const double &t, const double &b, const bool logX, const bool logY, const bool logZ);

  /**
   * @brief Function to set the format of the latex to be printed
   *
   * @param x  X position
   * @param y  Y position
   * @param s  String to print
   *
   */
  void FormatLatex(const double &x, const double &y, const char *s);
} // calib

#endif
