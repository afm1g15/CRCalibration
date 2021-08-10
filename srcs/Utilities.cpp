#include "Utilities.h"

namespace calib{
  
  //------------------------------------------------------------------------------------------ 
 
  float GetDistanceToPlane(const Plane &plane, const TVector3 &vtx, const TVector3 &end){
    // Get the value of the unit vector of the particle dotted with the normal to the plane
    // If this is zero, the particle is parallel so throw an exception and catch it in the main
    // then continue looping through the list of planes
    TVector3 track_direction   = (end - vtx).Unit();
    float direction_from_plane = track_direction.Dot(plane.GetUnitN());

    if(std::abs(direction_from_plane) <= std::numeric_limits<float>::epsilon()){
      throw std::domain_error("The track is parallel to the plane");
    }
   
    std::cout << " A          : (" << plane.GetA().X() << ", " << plane.GetA().Y() << ", " << plane.GetA().Z() << ") " << std::endl;
    std::cout << " B          : (" << plane.GetB().X() << ", " << plane.GetB().Y() << ", " << plane.GetB().Z() << ") " << std::endl;
    std::cout << " V          : (" << plane.GetV().X() << ", " << plane.GetV().Y() << ", " << plane.GetV().Z() << ") " << std::endl;
    std::cout << " vtx        : (" << vtx.X() << ", " << vtx.Y() << ", " << vtx.Z() << ") " << std::endl;
    std::cout << " end        : (" << end.X() << ", " << end.Y() << ", " << end.Z() << ") " << std::endl;
    std::cout << " V - vtx    : (" << (plane.GetV() - vtx).X() << ", " << (plane.GetV() - vtx).Y() << ", " << (plane.GetV() - vtx).Z() << ") " << std::endl;
    std::cout << " 1/Direction:" << 1/direction_from_plane << std::endl;
    std::cout << " N          : (" << plane.GetUnitN().X() << ", " << plane.GetUnitN().Y() << ", " << plane.GetUnitN().Z() << ") " << std::endl;

    return (1/direction_from_plane)*((plane.GetV() - vtx).Dot(plane.GetUnitN()));
  } // Get distance to plane
  
  //------------------------------------------------------------------------------------------ 

  Plane GetClosestPlane(const PlaneList &planes, const TVector3 &vtx, const TVector3 &end){
    float minDist = 999.;
    unsigned int planeID = 0;
    unsigned int minID = 0;
    for(const Plane &pl : planes){
      float dist = GetDistanceToPlane(pl, vtx, end);
      std::cout << " Plane: " << pl.GetLabel() << ", dist: " << dist << std::endl;
      if(abs(dist) < minDist){
        minDist = dist;
        minID   = planeID;
      }
      planeID++;
    }
    return planes.at(minID);
  }
  
  //------------------------------------------------------------------------------------------ 

  bool CheckIfIntersectsPlane(const Plane &plane, const TVector3 &vtx, const TVector3 &end, const float &length){
    float d = -std::numeric_limits<float>::max();
    try{
      // Will throw exception if the particle is parallel to the plane
      d = GetDistanceToPlane(plane, vtx, end);
    }
    // If caught, the particle is parallel to the plane and does not intersect
    catch(const std::domain_error&) {return false;}

//    if(d < 0 ) 
//      std::cout << "Backwards" << std::endl;

    if(d < 0 || d > length) {
      std::cout << " Doesn't intersect, d: " << d << ", length: " << length << std::endl;
      return false;
    }
 //   if(d > length) return false;
    TVector3 track_direction    = (end - vtx).Unit();
    TVector3 intersection_point = vtx + d * track_direction;

    bool intersects = IsProjectedPointInPlaneBounds(intersection_point, plane);
    if(!intersects){
      std::cout << "Intersection point: (" << intersection_point.X() << ", " << intersection_point.Y() << ", " << intersection_point.Z() << ") " << std::endl;
    }
    return intersects;
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
  
  void SetCanvasStyle(TCanvas *c, const double &l, const double &r, const double &t, const double &b, const bool logX, const bool logY, const bool logZ){
    c->SetLeftMargin(l);
    c->SetRightMargin(r);
    c->SetTopMargin(t);
    c->SetBottomMargin(b);

    if(logX)
      c->SetLogx();
    if(logY)
      c->SetLogy();
    if(logZ)
      c->SetLogz();
  } // Canvas Style
  
  // --------------------------------------------------------------------------------------------------------------------------------------------------
  
  void SetHistogramStyle1D(TH1D *h, const char *xLabel, const char *yLabel){
    h->GetXaxis()->SetTitle(xLabel);
    h->GetYaxis()->SetTitle(yLabel);
    h->GetXaxis()->SetTitleFont(132);
    h->GetXaxis()->SetLabelFont(132);
    h->GetYaxis()->SetTitleFont(132);
    h->GetYaxis()->SetLabelFont(132);
    h->GetXaxis()->SetTitleSize(0.055);
    h->GetXaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetTitleSize(0.055);
    h->GetYaxis()->SetLabelSize(0.045);
    h->GetXaxis()->SetMaxDigits(3);
    h->GetYaxis()->SetMaxDigits(3);
    h->SetStats(0);
  } // 1D Histogram Style

  // --------------------------------------------------------------------------------------------------------------------------------------------------
  
  void SetHistogramStyle2D(TH2D *h, const char *xLabel, const char *yLabel){
    h->GetXaxis()->SetTitle(xLabel);
    h->GetYaxis()->SetTitle(yLabel);
    h->GetXaxis()->SetTitleFont(132);
    h->GetXaxis()->SetLabelFont(132);
    h->GetYaxis()->SetTitleFont(132);
    h->GetYaxis()->SetLabelFont(132);
    h->GetXaxis()->SetTitleSize(0.055);
    h->GetXaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetTitleSize(0.055);
    h->GetYaxis()->SetLabelSize(0.045);
    h->GetXaxis()->SetMaxDigits(3);
    h->GetYaxis()->SetMaxDigits(3);
    h->GetYaxis()->SetTitleOffset(0.9);
    h->SetContour(99);
    h->SetStats(0);
  } // 2D Histogram Style
  
  // --------------------------------------------------------------------------------------------------------------------------------------------------
  
  void SetHistogramStyle3D(TH3D *h, const char *xLabel, const char *yLabel, const char *zLabel){
    h->GetXaxis()->SetTitle(xLabel);
    h->GetYaxis()->SetTitle(yLabel);
    h->GetZaxis()->SetTitle(zLabel);
    h->GetXaxis()->SetTitleFont(132);
    h->GetXaxis()->SetLabelFont(132);
    h->GetYaxis()->SetTitleFont(132);
    h->GetYaxis()->SetLabelFont(132);
    h->GetZaxis()->SetTitleFont(132);
    h->GetZaxis()->SetLabelFont(132);
    h->GetXaxis()->SetTitleSize(0.055);
    h->GetXaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetTitleSize(0.055);
    h->GetYaxis()->SetLabelSize(0.045);
    h->GetZaxis()->SetTitleSize(0.055);
    h->GetZaxis()->SetLabelSize(0.045);
    h->GetXaxis()->SetMaxDigits(3);
    h->GetYaxis()->SetMaxDigits(3);
    h->GetZaxis()->SetMaxDigits(3);
    h->SetContour(99);
    h->SetStats(0);
  } // 2D Histogram Style
  
  // --------------------------------------------------------------------------------------------------------------------------------------------------
  
  void FormatLatex(const double &x, const double &y, const char *s){
    // Setup the Latex object
    TLatex l;
    l.SetTextAlign(11); // Align at bottom
    l.SetTextSize(0.05);
    l.SetTextFont(132);
    l.DrawLatex(x,y,s);
  }

  // --------------------------------------------------------------------------------------------------------------------------------------------------
  void WriteStatsToTeX(ofstream   &file,
                       const int  &nFiles,
                       const std::vector<std::string> &contents,
                       const std::vector<unsigned int> &rates,
                       const double &denom){
 
    // Calculate the approximate number of days from the number of files
    // nFiles * 500 events per file / 14101 events per day
    double nDays = (nFiles*500.)/14118.;
    // Start by writing the first few lines of the tex file
    file << "\\begin{document} " << std::endl;
    file << "  \\thispagestyle{empty}" << std::endl;
    file << "  \\renewcommand{\\arraystretch}{1.2}" << std::endl;

    // Setup the table
    file << "  \\begin{table}[h!]" << std::endl;
    file << "    \\centering" << std::endl;
    file << "    \\begin{tabular}{ m{4cm} * {2}{ >{\\centering\\arraybackslash}m{4cm} } }" << std::endl;
    file << "      \\toprule" << std::endl;
    file << "      Statistic & Rate / " << std::setprecision(4) << nDays << " Days & \\% All Tracks \\\\" << std::endl;
    file << "      \\midrule" << std::endl;

    for(unsigned int nRate = 0; nRate < rates.size(); ++nRate){
      std::string str = contents.at(nRate);
      unsigned int val = rates.at(nRate);
      if(nRate == 0){
        file << "      " << str << " & " << "\\num{" << std::setprecision(4) << val << "} & - \\\\" << std::endl;
        file << "      \\midrule" << std::endl;
      }
      else
        file << "      " <<  str << " & " << "\\num{" << std::setprecision(4) << val << "} & " << std::setprecision(5) << (val/denom)*100 << "~\\%  \\\\" << std::endl;
    }
    file << "      \\bottomrule" << std::endl;
    file << "    \\end{tabular}" << std::endl;
    file << "  \\end{table}" << std::endl;
    file << "\\end{document}" << std::endl;
  }
  
} // calib
