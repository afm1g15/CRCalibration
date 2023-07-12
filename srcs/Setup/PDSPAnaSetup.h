#ifndef PDSPANASETUP_H
#define PDSPANASETUP_H

#include "TROOT.h"
#include "TFile.h"
#include "TFrame.h"
#include "TLine.h"
#include "TLatex.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TColor.h"
#include "TLegend.h"
#include "TMath.h"
#include "TVector3.h"
#include "TString.h"
#include "TPaletteAxis.h"
#include "TPaveStats.h"
#include "TRandom3.h"
#include "TView3D.h"
#include "TPolyLine3D.h"
#include "TPolyMarker3D.h"

#include "anatree_pdspana.h"
#include "../Classes/Plane.h"
#include "../Classes/Geometry.h"
#include "../Helpers/Constants.h"

#include <cstdlib>
#include <numeric>
#include <time.h>
#include <ctime>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <map>
#include <vector>
#include <string>

namespace calib{

  /*
   * @brief Global typedefs
   *
   */
  typedef std::vector<Plane> PlaneList;

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
   * @brief Find and replace function for c++ strings
   *
   * @param s            The main string to find and replace within
   * @param toReplace    The substring to find
   * @param replaceWith  The substring to replace the initial substring with
   * @param verbose      Whether or not to print in the verbose format, for debugging
   *
   */
  void FindReplace(std::string& s,
                   std::string const& toReplace,
                   std::string const& replaceWith,
                   bool verbose = false);

} // calib

#endif
