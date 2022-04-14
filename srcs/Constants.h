// Physical constants used throughout the analysis
// Many taken from 
//    $LARCOREOBJ_DIR/.../SimpleTypesAndConstants/PhysicalConstants.h

#ifndef CONSTANTS_H
#define CONSTANTS_H

namespace calib{

  /* @brief Conversion for energy deposited in GeV to number of ionization electrons produced */
  constexpr double kGeVToElectrons = 4.237e7; ///< 23.6eV per ion pair, 1e9 eV/GeV
  constexpr double kElectronsToADC = 5e-3; // Conversion from number of electrons to ADC for charge
  constexpr double kElectronsToADCTuned = 5.80906e-3; // Conversion from number of electrons to ADC for charge from my calibration work
  //constexpr double kElectronsToADC = 6.8906513e-3; // Conversion from number of electrons to ADC for charge from fcl

  /* @name Detector properties, some found in detectorproperties_dune.fcl
   *
   * @param Argon density in g/cm^3 
   * @param EField DUNE's electric field in kV/cm
   *
   */
  constexpr double kArDensity = 1.3954;
  constexpr double kEField    = 0.5;


  /**
   * @name Recombination factor coefficients (modified box, ArguNeuT JINST).
   * @see sim::ISCalculationSeparate::CalculateIonizationAndScintillation()
   *
   * Recombination factor coefficients come from Nucl.Instrum.Meth.A523:275-286,2004
   * * @f$ dE/dx @f$ is given by the voxel energy deposition, but have to convert it to MeV/cm from GeV/ voxel width
   * * electric field: @f$ E @f$ in kV/cm
   * * `kModBoxB` needs to be scaled with the electric field.
   *
   */
  constexpr double kModBoxA        = 0.930;    ///< Modified Box Alpha
  constexpr double kModBoxB        = 0.212;    ///< Modified Box Beta in g/(MeV cm&sup2;)*kV/cm
  constexpr double kWion           = 23.6e-6;  ///< ionization potenial in LAr, 23.6 eV = 1e, Wion in  MeV/e

}

#endif
