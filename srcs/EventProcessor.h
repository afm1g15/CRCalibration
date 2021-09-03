#ifndef EVENTPROCESSOR_H
#define EVENTPROCESSOR_H

#include "ReadFiles.h"

namespace calib{
  class EventProcessor{
    public:
      EventProcessor(const std::vector<TString> &allowedBranches, const std::string &inputList, int &nFiles);
      ~EventProcessor(){}

      void Initialize();
      void Finalize();

      int GetFiles() const;
      anatree * GetEvents() const;
      TChain  * GetTree() const;

      bool SelectEvent(anatree* evt) const;
      bool SelectTrack(anatree* evt, int iTrk) const;
      bool SelectHit(anatree* evt, int iTrk, int iPlane, int iHit) const;

      int WhichTPC(double x);

      const double APA_X_POSITIONS[3] = {-726.7681, 0., 726.7681};
      const double CPA_X_POSITIONS[2] = {-363.38405, 363.38405};
      const double TPC_Z_SIZE = 232.39;
      const double kXtoT = 1./160.563; // converts cm to ms, calculated for field 0.5kV/cm and temperature 87K

    private:
      TChain* fTree;
      anatree* fEvent;
      const char *fInputList;
      int fFiles;
      std::vector<TString> fAllowed;

  }; // EventProcessor
} // calib
#endif
