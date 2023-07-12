#include "GeometryHelpers.h"

namespace calib{
  
  //------------------------------------------------------------------------------------------ 
 
  float GetDistanceToPlane(const Plane &plane, const TVector3 &vtx, const TVector3 &end){
    // Get the value of the unit vector of the particle dotted with the normal to the plane
    // If this is zero, the particle is parallel so throw an exception and catch it in the main
    TVector3 track_direction   = (end-vtx).Unit();
    float direction_from_plane = track_direction.Dot(plane.GetUnitN());

    if(std::abs(direction_from_plane) <= std::numeric_limits<float>::epsilon()) return std::numeric_limits<float>::max();

    /*
    std::cout << " Plane: " << plane.GetLabel() << std::endl;
    std::cout << " A          : (" << plane.GetA().X() << ", " << plane.GetA().Y() << ", " << plane.GetA().Z() << ") " << std::endl;
    std::cout << " B          : (" << plane.GetB().X() << ", " << plane.GetB().Y() << ", " << plane.GetB().Z() << ") " << std::endl;
    std::cout << " V          : (" << plane.GetV().X() << ", " << plane.GetV().Y() << ", " << plane.GetV().Z() << ") " << std::endl;
    std::cout << " vtx        : (" << vtx.X() << ", " << vtx.Y() << ", " << vtx.Z() << ") " << std::endl;
    std::cout << " end        : (" << end.X() << ", " << end.Y() << ", " << end.Z() << ") " << std::endl;
    std::cout << " V - vtx    : (" << (plane.GetV() - vtx).X() << ", " << (plane.GetV() - vtx).Y() << ", " << (plane.GetV() - vtx).Z() << ") " << std::endl;
    std::cout << " 1/Direction:" << 1/direction_from_plane << std::endl;
    std::cout << " N          : (" << plane.GetUnitN().X() << ", " << plane.GetUnitN().Y() << ", " << plane.GetUnitN().Z() << ") " << std::endl;
    */

    // Multiply by 1 over the direction to the normal to account for the angle of the track relative to the plane
    return (1./direction_from_plane)*((plane.GetV() - vtx).Dot(plane.GetUnitN()));
  } // Get distance to plane
  
  //------------------------------------------------------------------------------------------ 

  Plane GetClosestPlane(const PlaneList &planes, const TVector3 &vtx, const TVector3 &end){
    float minDist = 999999.;
    unsigned int planeID = 0;
    unsigned int minID = 0;
    for(const Plane &pl : planes){
      float dist = GetDistanceToPlane(pl, vtx, end);
      if(abs(dist) < minDist){
        minDist = abs(dist);
        minID   = planeID;
      }
      planeID++;
    }
    return planes.at(minID);
  }
  
  //------------------------------------------------------------------------------------------ 

  bool CheckIfIntersectsPlane(const Plane &plane, const TVector3 &vtx, const TVector3 &end, const float &length){
    float d = GetDistanceToPlane(plane, vtx, end);

    // ignore if parallel
    if(abs(d - std::numeric_limits<float>::max()) <= std::numeric_limits<float>::epsilon()){
      return false;
    } 

    if(d < 0 || d > length) return false;
    TVector3 track_direction    = (end - vtx).Unit();
    TVector3 intersection_point = vtx + d * track_direction;

    bool intersects = IsProjectedPointInPlaneBounds(intersection_point, plane);
    return intersects;
  }

  //------------------------------------------------------------------------------------------ 

  unsigned int GetNExtPlanesCrossed(const double &length, const TVector3 &vtx, const TVector3 &end, const PlaneList &ext, const PlaneList &fidExt){
    // First, check if we come sufficiently close to and external plane
    Plane enteringPlane     = GetClosestPlane(ext, vtx, end);
    double distFromEntrance = GetDistanceToPlane(enteringPlane, vtx, end);

    // Find the closest plane to the end vertex and count it as a crossing plane
    Plane exitingPlane  = GetClosestPlane(ext, end, vtx);
    double distFromExit = GetDistanceToPlane(exitingPlane, end, vtx);

    unsigned int nExtCrossed = 0;
    for(const Plane &pl : ext){
      if(enteringPlane.GetLabel() == pl.GetLabel()){
        if(distFromEntrance < 30.){
          nExtCrossed++;
        }
      }
      else if(exitingPlane.GetLabel() == pl.GetLabel()){
        if(distFromExit < 30.){
          nExtCrossed++;
        }
      }
      // Counter for the number of external planes this track has crossed
      else if(CheckIfIntersectsPlane(pl,vtx,end,length)){
        nExtCrossed++;
      } // Intersects
    } // Planes
    return nExtCrossed;
  }

  //------------------------------------------------------------------------------------------ 

  void GetPlaneBoundaryCoordinates(const Plane &plane, std::vector<TVector3> &coords){

    // V should be the central point in each plane, alpha and beta are the full 2D dimenstions
    // V +/- A/2 and B/2 should give the boundaries
    coords.push_back(plane.GetV()+(plane.GetA())+(plane.GetB()));
    coords.push_back(plane.GetV()+(plane.GetA())-(plane.GetB()));
    coords.push_back(plane.GetV()-(plane.GetA())-(plane.GetB()));
    coords.push_back(plane.GetV()-(plane.GetA())+(plane.GetB()));

  }

  //------------------------------------------------------------------------------------------ 

  bool IsThroughGoing(const double &length, const TVector3 &vtx, const TVector3 &end, const PlaneList &ext, const PlaneList &fidExt){
    int nExtCrossed = GetNExtPlanesCrossed(length,vtx,end,ext,fidExt);
    if(nExtCrossed >= 2) return true;
    else return false;

  }
  
  //------------------------------------------------------------------------------------------ 

  bool IsStopping(const double &length, const TVector3 &vtx, const TVector3 &end, const PlaneList &ext, const PlaneList &fidExt){
    int nExtCrossed = GetNExtPlanesCrossed(length,vtx,end,ext,fidExt);
    if(nExtCrossed == 1) return true;
    else return false;

  }
  
  //------------------------------------------------------------------------------------------ 

  bool IsTrueThroughGoing(const TVector3 &vtx, const TVector3 &end, const TVector3 &vtxAV, const TVector3 &endAV){
      
    float dx = abs(endAV.X()-end.X());
    float dy = abs(endAV.Y()-end.Y());
    float dz = abs(endAV.Z()-end.Z());

    float dxS = abs(vtxAV.X()-vtx.X());
    float dyS = abs(vtxAV.Y()-vtx.Y());
    float dzS = abs(vtxAV.Z()-vtx.Z());

    float dE = dx+dy+dz;
    float dS = dxS+dyS+dzS;

    // If these match, the TPC end point and general end point are the same, therefore the particle stops
    if(dE < 1e-10 || dS < 1e-10) return false;
    else return true;
  
  }
  
  //------------------------------------------------------------------------------------------ 

  bool IsTrueStopping(const TVector3 &vtx, const TVector3 &end, const TVector3 &vtxAV, const TVector3 &endAV){
      
    // If these match, the TPC end point and general end point are the same, therefore the particle stops
    if(!IsTrueThroughGoing(vtx,end,vtxAV,endAV)) return true;
    else return false;
  
  }
  
  //------------------------------------------------------------------------------------------ 
  bool IsProjectedPointInPlaneBounds(const TVector3 &point, const Plane &plane){
    // Check if the point lies within the bound plane
    return (std::abs((point-plane.GetV()).Dot(plane.GetUnitA())) <= plane.GetAlpha() && std::abs((point-plane.GetV()).Dot(plane.GetUnitB())) <= plane.GetBeta());
  }
  
  //------------------------------------------------------------------------------------------ 
 
  bool CheckExternal(const Geometry &geom, const Plane &pl){
    PlaneList extPlanes = geom.GetExternalPlaneList();
    if(extPlanes.size() == 0){
      std::cerr << " Error: No external planes defined" << std::endl;
      return false;
    }
    for(const Plane &ext : extPlanes){
      if(pl.GetLabel() == ext.GetLabel())
        return true;
    }
    return false;
  }
  
  //------------------------------------------------------------------------------------------ 

  int GetBestPlane(const std::vector<int> &nHitsPerPlane){
    int bestPlane = -1;
    int currHits  = -999;
    for(int iPlane = 0; iPlane < 3; ++iPlane){
      if(nHitsPerPlane.at(iPlane) > currHits){
        currHits  = nHitsPerPlane.at(iPlane);
        bestPlane = iPlane; 
      } // CurrHits
    } // Planes
    if(bestPlane < 0){
      std::cerr << " Error: Best plane hasn't been set" << std::endl;
      std::exit(1);
    }
    return bestPlane;
  }
  
  // --------------------------------------------------------------------------------------------------------------------------------------------------
  
  void GetRecoBestPlane(const int &iTrk, const anatree *evt, int &bP, std::vector<int> &hits){
    // Get the best plane
    int currHits  = -999;
    for(int iPlane = 0; iPlane < 3; ++iPlane){
      hits.at(iPlane) = evt->ntrkhits_pandoraTrack[iTrk][iPlane];
      if(evt->ntrkhits_pandoraTrack[iTrk][iPlane] > currHits){
        currHits = evt->ntrkhits_pandoraTrack[iTrk][iPlane];
        bP       = iPlane; 
      } // CurrHits
    } // Planes
  }
  
  // --------------------------------------------------------------------------------------------------------------------------------------------------
  
  double GetCosTheta(const int &i, const TVector3 &vtx, const TVector3 &end){
  
    // First, get the direction of the wires on the plane
    TVector3 wireDir = wireDirections.at(i);

    // Now get the direction from start and end 
    TVector3 dir = (end-vtx);
    //dir *= 1/static_cast<double>(end.Mag()*vtx.Mag());

    // Now get the angle between them
    double costheta = wireDir.Dot(dir)/static_cast<double>(wireDir.Mag()*dir.Mag());

    return costheta;
  }

  // --------------------------------------------------------------------------------------------------------------------------------------------------
  
  double GetSinTheta(const int &i, const TVector3 &vtx, const TVector3 &end){
  
    // First, get the direction of the wires on the plane
    TVector3 wireDir = wireDirections.at(i);

    // Now get the direction from start and end 
    TVector3 dir = (end-vtx);
    dir *= 1/static_cast<double>(end.Mag()*vtx.Mag());

    // Now get the angle between them
    double sintheta = (wireDir.Cross(dir)).Mag()/static_cast<double>(wireDir.Mag()*dir.Mag());

    return sintheta;
  }

  // --------------------------------------------------------------------------------------------------------------------------------------------------
  
  double GetAngleToAPAs(const TVector3 &norm, const TVector3 &vtx, const TVector3 &end){

    // Get the unit direction of the track 
    TVector3 trackDir = (end-vtx);
    trackDir *= 1/static_cast<double>(end.Mag()*vtx.Mag());

    // Get the angle of the track to the wire planes
    double costoplane = norm.Dot(trackDir)/static_cast<double>(norm.Mag()*trackDir.Mag());

    return costoplane;

  }

  // --------------------------------------------------------------------------------------------------------------------------------------------------
  
  double GetHitPitch(const int &plane, const TVector3 &currXYZ, const TVector3 &nextXYZ){

    // pitch = 0.48/sintheta
    double sinTheta = GetSinTheta(plane,currXYZ,nextXYZ);

    // Get the pitch
    double hitPitch = 0.48/sinTheta;

    return hitPitch;
  }
  
  // --------------------------------------------------------------------------------------------------------------------------------------------------
  
  double GetCosDrift(const TVector3 &currXYZ, const TVector3 &nextXYZ){
    
    // Get the angle of the hit to the x axis
    TVector3 xDir(1,0,0);
    TVector3 hitVec = (nextXYZ-currXYZ);

    double cosDrift = hitVec.Dot(xDir)/static_cast<double>(hitVec.Mag()*xDir.Mag());
    return cosDrift;
  }
  
  // --------------------------------------------------------------------------------------------------------------------------------------------------
  
  double GetAngleToPlane(const TVector3 &startXYZ, const TVector3 &endXYZ, const TVector3 planeDir){
 
    // Get the track direction
    TVector3 trackDir = (endXYZ-startXYZ);
    double cosPlane = trackDir.Dot(planeDir)/static_cast<double>(trackDir.Mag()*planeDir.Mag());
    return cosPlane;
  }
} // calib
