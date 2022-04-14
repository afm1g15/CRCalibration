#ifndef EXTERNALFUNCTIONHELPERS_H
#define EXTERNALFUNCTIONHELPERS_H

#include "Setup.h"

namespace calib{

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
   * @brief Extract the MPV and FWHM from the langau fit output from ROOT: http://merlot.ijs.si/~matevz/docs/RootHtmlDoc-5.22-00/tutorials/fit/langaus.C.html
   *
   * @param params Langau fit paramters
   * @param maxx  Langau MPV
   * @param FWHM  Langau full width half max
   *
   * @return function at x and par
   */
  int langaupro(double *params, double &maxx, double &FWHM);

  /*
   * @brief ModBox correction for calculating the reconstructed dE/dx from dQ/dx from $LARDATAALG_DIR/.../DetectorInfo/DetectorPropertiesStandard.cxx
   *
   * @param dQdx    The reconstructed dQ/dx we wish to convert
   * @param tuned   Whether or not to use the tuned value for Cscale from my analysis
   * @param eField  The electric field across the detector, default 500 V/cm
   *
   * @return dEdx  The converted dE/dx
   */
  double ModBoxCorrection(const double &dQdx, const bool &tuned = false, const double &eField = calib::kEField);

} // calib
#endif
