#include "Geometry.h"

namespace calib{

  // Geometry constructor
  Geometry::Geometry(std::vector<double> minX,
                     std::vector<double> minY,
                     std::vector<double> minZ,
                     std::vector<double> maxX,
                     std::vector<double> maxY,
                     std::vector<double> maxZ,
                     const bool fiducial): 
  m_min_x(minX),
  m_min_y(minY),
  m_min_z(minZ),
  m_max_x(maxX),
  m_max_y(maxY),
  m_max_z(maxZ),
  m_fiducial(fiducial)
  {
    const std::vector< std::vector<double> > limits = {m_min_x, m_min_y, m_min_z, m_max_x, m_max_y, m_max_z};
    const bool areSameSize = std::all_of(limits.begin(), limits.end(), [&](const std::vector<double> &x){return x.size() == limits.front().size();});
    if(!areSameSize){
      std::cerr << " Error: You are constructing a Geometry object with varying numbers of TPCs:" << std::endl;
      std::cerr << " Number of minimum x values: " << m_min_x.size() << std::endl;
      std::cerr << " Number of minimum y values: " << m_min_y.size() << std::endl;
      std::cerr << " Number of minimum z values: " << m_min_z.size() << std::endl;
      std::cerr << " Number of maximum x values: " << m_max_x.size() << std::endl;
      std::cerr << " Number of maximum y values: " << m_max_y.size() << std::endl;
      std::cerr << " Number of maximum z values: " << m_max_z.size() << std::endl;
      std::exit(1);
    }
   m_n_tpcs = m_min_x.size();
   if(m_n_tpcs == 0){
    std::cerr << " Error: Can't construct an emtpy Geometry object, exiting " << std::endl;
    std::exit(1);
   }

   if(m_fiducial)
     std::cout << " Constructing a fiducial volume Geometry object with " << m_n_tpcs << " TPCs " << std::endl;
   else
     std::cout << " Constructing an active volume Geometry object with " << m_n_tpcs << " TPCs " << std::endl;
   
  }

  //------------------------------------------------------------------------------------------ 
  
  Geometry::PlaneList Geometry::GetPlaneList() const {
    // Define the planes of the detector for each TPC
    PlaneList planes;

    //    Find the minimum element of the min_x vectors
    //    Find the maximum element of the max_x vectors
    double min_x = *std::min_element(m_min_x.begin(),m_min_x.end());
    double max_x = *std::max_element(m_max_x.begin(),m_max_x.end());
    
    //  The order: +/-x, +/-y, +/-z
    //
    // To take into account the fiducial border define 
    //    Half with fiducial   = Half width  - width fiducial border 
    //    Half height fiducial = Half height - height fiducial border
    //    Half length fiducial = Half length - length fiducial border
    for(unsigned int n = 0; n < m_n_tpcs; ++n){
      double h = (m_max_y.at(n)-m_min_y.at(n))*0.5; // half full detector height (y)
      double l = (m_max_z.at(n)-m_min_z.at(n))*0.5; // half full detector length (z)
      double w = (max_x-min_x)*0.5;

      // Define the planes of the detector, make sure there are no duplicates
      Plane tempPl1(TVector3(m_max_x.at(n), h, l), TVector3(0, 2*h, 0), TVector3(0, 0, -2*l), "h"+std::to_string(n+1));
      Plane tempPl2(TVector3(m_min_x.at(n), h, l), TVector3(0, 2*h, 0), TVector3(0, 0,  2*l), "h"+std::to_string(n));
      Plane tempPl3(TVector3(w, m_max_y.at(n), l), TVector3(2*w, 0, 0), TVector3(0, 0,  2*l), "t");
      Plane tempPl4(TVector3(w, m_min_y.at(n), l), TVector3(2*w, 0, 0), TVector3(0, 0, -2*l), "bo");
      Plane tempPl5(TVector3(w, h, m_min_z.at(n)), TVector3(2*w, 0, 0), TVector3(0, 2*h,  0), "f");
      Plane tempPl6(TVector3(w, h, m_max_z.at(n)), TVector3(2*w, 0, 0), TVector3(0,-2*h,  0), "ba");

      PlaneList newPlanes{tempPl1,tempPl2,tempPl3,tempPl4,tempPl5,tempPl6};

      for(const Plane &newPl : newPlanes){
        if(planes.size() > 0){
          bool found = false;
          for(const Plane &pl : planes){
            if(pl == newPl){
              found = true;
              break;
            }
          }
          if(!found)
            planes.push_back(newPl);
        }
        else{
          planes.push_back(newPl);
        } // IfPlanes
      } // NewPlanes
    } // TPCs
    return planes;
  }

  //------------------------------------------------------------------------------------------ 
  
  Geometry::PlaneList Geometry::GetExternalPlaneList() const {
    // Define the planes of the detector for each TPC
    PlaneList planes;

    // Only get the planes which correspond to the outer faces of the detector, this still requires looking at each TPC
    // First:
    //    Find the minimum element of the min_x vectors
    //    Find the maximum element of the max_x vectors
    double min_x = *std::min_element(m_min_x.begin(),m_min_x.end());
    double max_x = *std::max_element(m_max_x.begin(),m_max_x.end());
    double min_y = *std::min_element(m_min_y.begin(),m_min_y.end());
    double max_y = *std::max_element(m_max_y.begin(),m_max_y.end());
    double min_z = *std::min_element(m_min_z.begin(),m_min_z.end());
    double max_z = *std::max_element(m_max_z.begin(),m_max_z.end());
    
    double w = (max_x-min_x)*0.5; // half full detector width  (x)
    double h = (max_y-min_y)*0.5; // half full detector height (y)
    double l = (max_z-min_z)*0.5; // half full detector length (z)

    planes.emplace_back(TVector3(max_x, h, l), TVector3(0, 2*h, 0), TVector3(0, 0, -2*l), "h"+std::to_string(m_n_tpcs));
    planes.emplace_back(TVector3(min_x, h, l), TVector3(0, 2*h, 0), TVector3(0, 0,  2*l), "h0");
    planes.emplace_back(TVector3(w, max_y, l), TVector3(2*w, 0, 0), TVector3(0, 0,  2*l), "t");
    planes.emplace_back(TVector3(w, min_y, l), TVector3(2*w, 0, 0), TVector3(0, 0, -2*l), "bo");
    planes.emplace_back(TVector3(w, h, min_z), TVector3(2*w, 0, 0), TVector3(0, 2*h,  0), "f");
    planes.emplace_back(TVector3(w, h, max_z), TVector3(2*w, 0, 0), TVector3(0,-2*h,  0), "ba");

    return planes;
  }

  //------------------------------------------------------------------------------------------ 
  
  Geometry::PlaneList Geometry::GetInternalPlaneList(const PlaneList &all, const PlaneList &ext) const {
    // Using the planelist and external plane list, determine the internal planes
    PlaneList planes;

    // Loop over all planes and determine which are already defined as external
    for(const Plane &allPl : all){
      bool isExternal = false;
      for(const Plane &extPl : ext){
        if(allPl == extPl){
          isExternal = true;
          break;
        }
      }
      // If we did not find a matching external plane, this is an internal plane
      if(!isExternal)
        planes.push_back(allPl);
    }
    return planes;
  }

  //------------------------------------------------------------------------------------------ 
  
  bool Geometry::GetIsFiducial() const {return m_fiducial;}
  
  unsigned int Geometry::GetNTPCs() const {return m_n_tpcs;}

  const std::vector<double> &Geometry::GetMinX() const {return m_min_x;}
  const std::vector<double> &Geometry::GetMinY() const {return m_min_y;}
  const std::vector<double> &Geometry::GetMinZ() const {return m_min_z;}
  const std::vector<double> &Geometry::GetMaxX() const {return m_max_x;}
  const std::vector<double> &Geometry::GetMaxY() const {return m_max_y;}
  const std::vector<double> &Geometry::GetMaxZ() const {return m_max_z;}

  //------------------------------------------------------------------------------------------ 
  
} // Selection
