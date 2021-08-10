#include <string>
#include <vector>
#include <cmath>
#include "TVector3.h"
#include "Plane.h"

namespace calib{

  // Plane constructors
  Plane::Plane() :
    m_alpha(1.f),
    m_beta(1.f),
    m_label("")
  {}

  Plane::Plane(const TVector3 &V, const TVector3 &A, const TVector3 &B, const std::string &l) :
    m_V(V),
    m_A(A),
    m_B(B),
    m_label(l){
      m_a = ((1/m_A.Mag()) * m_A);
      m_b = ((1/m_B.Mag()) * m_B);
      m_n = m_a.Cross(m_b);
      m_alpha = m_A.Mag();
      m_beta  = m_B.Mag();
  }

  // Overloaded functions
  bool operator==(const Plane &lhs, const Plane &rhs){
    return (lhs.m_V     == rhs.m_V && 
            lhs.m_label == rhs.m_label);
  }
  
  bool operator!=(const Plane &lhs, const Plane &rhs){
    return !(lhs == rhs);  
  }

  // Plane object getters
  TVector3 Plane::GetV() const {return m_V;}
  TVector3 Plane::GetA() const {return m_A;}
  TVector3 Plane::GetB() const {return m_B;}
  TVector3 Plane::GetUnitA() const {return m_a;}
  TVector3 Plane::GetUnitB() const {return m_b;}
  TVector3 Plane::GetUnitN() const {return m_n;}
  float Plane::GetAlpha() const {return m_alpha;}
  float Plane::GetBeta()  const {return m_beta;}
  std::string Plane::GetLabel() const {return m_label;}

  // Set the label
  void Plane::SetLabel(const std::string &label){
    m_label = label;
  }

} // calib
