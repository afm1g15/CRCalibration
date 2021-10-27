#ifndef SETUP_H
#define SETUP_H

#include "TROOT.h"
#include "TFile.h"
#include "TLine.h"
#include "TLatex.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TVector3.h"
#include "TString.h"
#include "TPaletteAxis.h"
#include "TPaveStats.h"

#include "anatree.h"
#include "Plane.h"
#include "Geometry.h"

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

} // calib

#endif
