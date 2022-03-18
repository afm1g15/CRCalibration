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
  
} // calib
#endif
