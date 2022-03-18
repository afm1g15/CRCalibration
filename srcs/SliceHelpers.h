#ifndef SLICEHELPERS_H
#define SLICEHELPERS_H

#include "Setup.h"

namespace calib{

  /**
   * @brief Fill the vectors containing information on how to slice the 2D histograms
   *
   * @param h Histogram to slice
   * @param nSlices Number of slices to define
   * @param binWidths Width of the bin slices
   * @param log whether to define slices in log space
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
   * @brief Fill the vectors containing information on how to slice the 2D histograms
   *        distribute the total number of entries as evenly as possible
   *
   * @param h Histogram to slice
   * @param nSlices Number of slices to define
   * @param log whether to define slices in log space
   * @param minX Vector of minimum x slice limits to fill
   * @param maxX Vector of maximum x slice limits to fill
   * @buffer buffer to apply to either end so slices don't get too close
   */
  void FillSliceVectorsEqualRates(TH2D *h,
                                  int &nSlices, 
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

} // calib
#endif
