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

    //  The order: +/-x, +/-y, +/-z
    //
    // To take into account the fiducial border define 
    //    Half with fiducial   = Half width  - width fiducial border 
    //    Half height fiducial = Half height - height fiducial border
    //    Half length fiducial = Half length - length fiducial border
    for(unsigned int n = 0; n < m_n_tpcs; ++n){
      float w = (m_max_x.at(n)-m_min_x.at(n))*0.5; // half full detector width (x)
      float h = (m_max_y.at(n)-m_min_y.at(n))*0.5; // half full detector height (y)
      float l = (m_max_z.at(n)-m_min_z.at(n))*0.5; // half full detector length (z)

      // Adjust for relative locations
      w = m_min_x.at(n)+w;
      h = m_min_y.at(n)+h;
      l = m_min_z.at(n)+l;

      // Define the planes of the detector, make sure there are no duplicates
      Plane tempPl1(TVector3(m_max_x.at(n), h, l), TVector3(0, h, 0), TVector3(0, 0,  l));
      Plane tempPl2(TVector3(m_min_x.at(n), h, l), TVector3(0, h, 0), TVector3(0, 0, -l));
      Plane tempPl3(TVector3(w, m_max_y.at(n), l), TVector3(w, 0, 0), TVector3(0, 0, -l));
      Plane tempPl4(TVector3(w, m_min_y.at(n), l), TVector3(w, 0, 0), TVector3(0, 0,  l));
      Plane tempPl5(TVector3(w, h, m_min_z.at(n)), TVector3(w, 0, 0), TVector3(0,  h, 0));
      Plane tempPl6(TVector3(w, h, m_max_z.at(n)), TVector3(w, 0, 0), TVector3(0, -h, 0));

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

    // This time, only define the outermost planes, ignoring any internal gaps
    //  The order: +/-x, +/-y, +/-z
    //
    // To take into account the fiducial border define 
    //    Half with fiducial   = Half width  - width fiducial border 
    //    Half height fiducial = Half height - height fiducial border
    //    Half length fiducial = Half length - length fiducial border
    //
    // First:
    //    Find the minimum element of the min_[] vectors
    //    Find the maximum element of the max_[] vectors
    float min_x = *std::min_element(m_min_x.begin(),m_min_x.end());
    float min_y = *std::min_element(m_min_y.begin(),m_min_y.end());
    float min_z = *std::min_element(m_min_z.begin(),m_min_z.end());
    float max_x = *std::max_element(m_max_x.begin(),m_max_x.end());
    float max_y = *std::max_element(m_max_y.begin(),m_max_y.end());
    float max_z = *std::max_element(m_max_z.begin(),m_max_z.end());
    float w = (max_x-min_x)*0.5; // half full detector width (x)
    float h = (max_y-min_y)*0.5; // half full detector height (y)
    float l = (max_z-min_z)*0.5; // half full detector length (z)

    // Define the planes of the detector
    planes.emplace_back(TVector3( w, 0, l),         TVector3(0, h, 0), TVector3(0, 0,  l));
    planes.emplace_back(TVector3(-w, 0, l),         TVector3(0, h, 0), TVector3(0, 0, -l));
    planes.emplace_back(TVector3(0,  h, l),         TVector3(w, 0, 0), TVector3(0, 0, -l));
    planes.emplace_back(TVector3(0, -h, l),         TVector3(w, 0, 0), TVector3(0, 0,  l));
    planes.emplace_back(TVector3(0, 0, 0),          TVector3(w, 0, 0), TVector3(0,  h, 0));
    planes.emplace_back(TVector3(0, 0, (2*l)),      TVector3(w, 0, 0), TVector3(0, -h, 0));

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
