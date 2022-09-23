/*
 * @brief Helper functions for using the DUNE FD HD SP detector geometry
 *
 */

#ifndef GEOMETRYHELPERS_H
#define GEOMETRYHELPERS_H

#include "../Setup/Setup.h"

namespace calib{

  std::vector<TVector3> wireDirections{
    TVector3(0, 0.812, 0.584),
    TVector3(0,-0.812, 0.584),
    TVector3(0,     1,     0)
  };

  /**
   * @brief  The the distance to a defined plane
   *
   * @param  plane 
   * @param  vtx
   * @param  end
   * 
   * @return the distance from the neutrino vertex to the plane
   */
  float GetDistanceToPlane(const Plane &plane, const TVector3 &vtx, const TVector3 &end);

  /**
   * @brief Find the closest plane to the vertex to determine where the track entered
   *
   * @param planes List of planes to check
   * @param vtx Start position of current track
   * @param end End position of current track
   *
   * @return plane Closest plane to vertex
   */
  Plane GetClosestPlane(const PlaneList &planes, const TVector3 &vtx, const TVector3 &end);

  /**
   * @brief  Check if it intersects a given plane
   *
   * @param  plane
   * @param  vtx
   * @param  end
   * @param  length
   *
   * @return true or false
   */
  bool CheckIfIntersectsPlane(const Plane &plane, const TVector3 &vtx, const TVector3 &end, const float &length);

  /**
   * @brief Count the number of external planes crossed by a track
   *
   * @param length  Length of the track
   * @param vtx     Vertex location of the start position
   * @param end     Vertex location of the end position
   * @param ext     List of external detector planes
   * @param fidExt  List of external fiducial detector planes
   *
   * @return true if track leaves the detector
   */
  unsigned int GetNExtPlanesCrossed(const double &length, const TVector3 &vtx, const TVector3 &end, const PlaneList &ext, const PlaneList &fidExt);
  
  /**
   * @brief Check if the track is through-going
   *
   * @param length  Length of the track
   * @param vtx     Vertex location of the start position
   * @param end     Vertex location of the end position
   * @param ext     List of external detector planes
   * @param fidExt  List of external fiducial detector planes
   *
   * @return true if track leaves the detector
   */
  bool IsThroughGoing(const double &length, const TVector3 &vtx, const TVector3 &end, const PlaneList &ext, const PlaneList &fidExt);

  /**
   * @brief Check if the track stops in the detector
   *
   * @param length  Length of the track
   * @param vtx     Vertex location of the start position
   * @param end     Vertex location of the end position
   * @param ext     List of external detector planes
   * @param fidExt  List of external fiducial detector planes
   *
   * @return true if track leaves the detector
   */
  bool IsStopping(const double &length, const TVector3 &vtx, const TVector3 &end, const PlaneList &ext, const PlaneList &fidExt);

  /**
   * @brief Check if the true G4 track is through-going
   *
   * @param vtx   Vertex location of the start position
   * @param end   Vertex location of the end position
   * @param vtxAV Vertex location of the start position in the AV
   * @param endAV Vertex location of the end position in the AV
   *
   * @return true if true track leaves the detector
   */
  bool IsTrueThroughGoing(const TVector3 &vtx, const TVector3 &end, const TVector3 &vtxAV, const TVector3 &endAV);

  /**
   * @brief Check if the true G4 track stops in the detector
   *
   * @param vtx   Vertex location of the start position
   * @param end   Vertex location of the end position
   * @param vtxAV Vertex location of the start position in the AV
   * @param endAV Vertex location of the end position in the AV
   *
   * @return true if true track stops in the detector
   */
  bool IsTrueStopping(const TVector3 &vtx, const TVector3 &end, const TVector3 &vtxAV, const TVector3 &endAV);
  
  /**
   * @brief  Check if the projected point is within the bounds of the given plane
   *
   * @param  point
   * @param  plane
   *
   * @return true or false
   */
  bool IsProjectedPointInPlaneBounds(const TVector3 &point, const Plane &plane);

  /**
   * @brief Check if the current plane is in the external list
   *
   * @param geom The geometry to check
   * @param pl Current plane to check
   *
   * @return if external
   */
  bool CheckExternal(const Geometry &geom, const Plane &pl);
  
  /**
   * @brief Function to get the plane with the most hits deposited for the current track
   *
   * @param nHitsPerPlane  Number of hits deposited by the track on each plane
   *
   * @return Iterator for the best plane
   */
  int GetBestPlane(const std::vector<int> &nHitsPerPlane);

  /**
   * @brief Get the identity of the best plane for the current reco track
   *
   * @param iTrk  Iterator of the current reco track
   * @param evt   Anatree event object
   * @param bP    BestPlane to allocate
   *
   */
  void GetRecoBestPlane(const int &iTrk, const anatree *evt, int &bP, std::vector<int> &hits);

  /**
   * @brief Get the angle between two points and the wires in the current plane
   *
   * @param i current plane index
   * @param vtx Vertex of the track
   * @param end End of the track
   *
   * @return costheta
   */
  double GetCosTheta(const int &i, const TVector3 &vtx, const TVector3 &end);
  
  /**
   * @brief Get the angle between two points and the wires in the current plane
   *
   * @param i current plane index
   * @param vtx Vertex of the track
   * @param end End of the track
   *
   * @return sintheta
   */
  double GetSinTheta(const int &i, const TVector3 &vtx, const TVector3 &end);
  
  /**
   * @brief Get the angle between the track and the wire plane
   *
   * @param norm Normal to the wire planes
   * @param vtx Vertex of the track
   * @param end End of the track
   *
   * @return costheta
   */
  double GetAngleToAPAs(const TVector3 &norm, const TVector3 &vtx, const TVector3 &end);

  /**
   * @brief Get the pitch of the current hit from the current and following hit location
   *
   * @param plane the wire plane 
   * @param currHitXYZ location of the current hit
   * @param nextHitXYZ location of the adjacent hit
   *
   * @return hit pitch
   */
  double GetHitPitch(const int &plane, const TVector3 &currXYZ, const TVector3 &nextXYZ);
  
  /**
   * @brief Get the angle of the current segment to the drift direction
   *
   * @param currHitXYZ location of the current hit
   * @param nextHitXYZ location of another hit hit
   *
   * @return cosDrift
   */
  double GetCosDrift(const TVector3 &currXYZ, const TVector3 &nextXYZ);

  /**
   * @brief Get the angle of the track to a given plane in radians
   *
   * @param startXYZ location of the current hit
   * @param nextXYZ location of next hit to assess
   * @param planeDir given plane co-ordinates
   *
   * @return angle to plane
   */
  double GetAngleToPlane(const TVector3 &startXYZ, const TVector3 &endXYZ, const TVector3 planeDir);

} // calib
#endif
