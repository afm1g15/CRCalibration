#include "Setup.h"

namespace calib{

  // --------------------------------------------------------------------------------------------------------------------------------------------------
  void GetTime(time_t &rawtime){
    struct tm * timeinfo;
    time (&rawtime);
    timeinfo = localtime (&rawtime);
    std::cout << " Local time and date:  " << asctime(timeinfo)                << std::endl;
  } // GetTime
  // --------------------------------------------------------------------------------------------------------------------------------------------------
  void GetTotalTime(time_t &starttime, time_t &endtime){
    double timediff = difftime(endtime,starttime);

    // Calculate the time in a nice format
    int seconds_per_minute = 60;
    int seconds_per_hour   = 3600;
    int hours              = timediff / seconds_per_hour;
    int minutes            = (timediff - (hours * seconds_per_hour)) / (seconds_per_minute);
    int seconds            = timediff - (hours * seconds_per_hour) - (minutes * seconds_per_minute);

    // Now print
    std::cout << " Total time taken: "    << std::setw(4) << hours << " hours, ";
    std::cout                             << std::setw(4) << minutes << " minutes, ";
    std::cout                             << std::setw(4) << seconds << " seconds." << std::endl;
  } // GetTotalTime
  // --------------------------------------------------------------------------------------------------------------------------------------------------
  void SetCanvasStyle2D(TCanvas *c, const double &l, const double &r, const double &t, const double &b, const bool logX, const bool logY, const bool logZ){
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

} // calib (namespace)
