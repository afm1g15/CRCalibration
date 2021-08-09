#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <cmath>
#include <algorithm>
#include <string>
#include <iostream>
#include <vector>
#include "TVector3.h"
#include "Plane.h"
#include "Setup.h"

namespace calib{
  
  /**
   * @brief  Geometry class
   */
  class Geometry{

    public  :

      typedef std::vector<Plane> PlaneList;
  
      /**
       * @brief  Constructor
       */
      Geometry(std::vector<double> minX,
               std::vector<double> minY,
               std::vector<double> minZ,
               std::vector<double> maxX,
               std::vector<double> maxY,
               std::vector<double> maxZ,
               const bool fiducial);

      /**
       * @brief Get the list of planes in this geometry
       *
       * @return  vector of Plane objects
       */
      PlaneList GetPlaneList() const;

      /**
       * @brief Get the list of external planes in this geometry 
       *
       * @return  vector of Plane objects
       */
      PlaneList GetExternalPlaneList() const;

      /**
       * @brief Get the list of internal planes in this geometry 
       *
       * @param all  The full plane list to check
       * @param ext  The external plane list to check against 
       *
       * @return  vector of Plane objects
       */
      PlaneList GetInternalPlaneList(const PlaneList &all, const PlaneList &ext) const;
      
      /**
       * @brief  Get whether this object is a fiducial or active geometry
       *
       * @return True if fiducial
       */
      bool GetIsFiducial() const;

      /**
       * @brief  Get the number of TPCs based on the size of the vectors
       *
       * @return Number of TPCs
       */
      unsigned int GetNTPCs() const;
      
      /**
       * @brief  Get list of minimum x values
       *
       * @return vector of minimum values
       */
      const std::vector<double> &GetMinX() const;
      
      /**
       * @brief  Get list of minimum y values
       *
       * @return vector of minimum values
       */
      const std::vector<double> &GetMinY() const;
      
      /**
       * @brief  Get list of minimum z values
       *
       * @return vector of minimum values
       */
      const std::vector<double> &GetMinZ() const;
      
      /**
       * @brief  Get list of maximum x values
       *
       * @return vector of maximum values
       */
      const std::vector<double> &GetMaxX() const;
      
      /**
       * @brief  Get list of maximum y values
       *
       * @return vector of maximum values
       */
      const std::vector<double> &GetMaxY() const;
      
      /**
       * @brief  Get list of maximum z values
       *
       * @return vector of maximum values
       */
      const std::vector<double> &GetMaxZ() const;

    private :

      // Member variables
      PlaneList m_planes;         ///< List of planes associated with this detector
      PlaneList m_ext_planes;     ///< List of planes associated with this detector, external only (no central gaps)
      PlaneList m_int_planes;     ///< List of planes associated with this detector, internal only
      bool m_fiducial;            ///< Boolean for if the geometry is fiducial (true) or active (false)
      unsigned int m_n_tpcs;      ///< Number of TPCs in this geometry
      std::vector<double> m_min_x; ///< Vector of minimum values of x, vector size corresponds to the number of TPCs
      std::vector<double> m_min_y; ///< Vector of minimum values of y, vector size corresponds to the number of TPCs
      std::vector<double> m_min_z; ///< Vector of minimum values of z, vector size corresponds to the number of TPCs
      std::vector<double> m_max_x; ///< Vector of maximum values of x, vector size corresponds to the number of TPCs
      std::vector<double> m_max_y; ///< Vector of maximum values of y, vector size corresponds to the number of TPCs
      std::vector<double> m_max_z; ///< Vector of maximum values of z, vector size corresponds to the number of TPCs

  }; // Geometry
} // calib
#endif
