/************************************************************************
 * 
 * A macro to understand the truth-level hit and track distributions for 
 * dE/dx calibration studies
 *
 *
 * Example file list located here:
 *   /home/jones/work/cosmics/LArSoft-v08_50_00/work/files/anafiles.list
 *
 *
 *************************************************************************/

#include "EventProcessor.h"
#include "ConfigReader.h"

#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TStyle.h>
#include <iostream>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <fstream>
#include <vector>
#include <cmath>
#include <math.h>
#include <string>
#include <TImage.h>
#include <iomanip>
#include <TSpline.h>
#include <TText.h>
#include <TFrame.h>
#include "TMinuit.h"
#include <sys/stat.h>
#include <string>

using namespace calib;
using namespace cppsecrets;

// Allowed branches to read from the tree
std::vector<TString> allowed = {
   "run",
   "event",
   "geant_list_size",
   "inTPCActive",
   "TrackId",
   "pdg",
   "Mother",
   "ntracks_pandoraTrack",
   "trkId_pandoraTrack",
   "trkidtruth_pandoraTrack",
   "trkpdgtruth_pandoraTrack",
   "trkg4id_pandoraTrack",
   "ntrkhits_pandoraTrack",
   "trkdqdx_pandoraTrack",
   "trkdedx_pandoraTrack",
   "trkke_pandoraTrack",
   "trkxyz_pandoraTrack",
   "trkstartx_pandoraTrack",
   "trkstarty_pandoraTrack",
   "trkstartz_pandoraTrack",
   "trkendx_pandoraTrack",
   "trkendy_pandoraTrack",
   "trkendz_pandoraTrack",
   "trklen_pandoraTrack",
   "no_hits_stored",
   "hit_tpc",
   "hit_plane",
   "hit_charge",
   "hit_energy",
   "hit_nelec",
   "hit_trueX",
   "hit_trkid",
   "hit_peakT",
   "hit_startT",
   "hit_endT",
   "StartPointx",
   "StartPointy",
   "StartPointz",
   "StartE_tpcAV",
   "EndPointx",
   "EndPointy",
   "EndPointz",
   "StartPointx_tpcAV",
   "StartPointy_tpcAV",
   "StartPointz_tpcAV",
   "EndPointx_tpcAV",
   "EndPointy_tpcAV",
   "EndPointz_tpcAV",
   "NumberDaughters",
   "P",
   "Eng",
   "EndE",
   "pathlen",
   "trkresrg_pandoraTrack"
 };

// A translation list from plane labels to longer labels for plotting
std::map<std::string, std::string> planeLabels = {
  {"h0", "APA 1"},
  {"h1", "CPA 1"},
  {"h2", "APA 2"},
  {"h3", "CPA 2"},
  {"h4", "APA 3"},
  {"t",  "Top"},
  {"bo", "Bot."},
  {"f",  "Fro."},
  {"ba", "Back"},
};
  

typedef std::vector<Plane> PlaneList;

template <typename T>
ostream &operator<<(ostream &out, const vector<T> &v)
{
  out << "{";
  size_t last = v.size() - 1;
  for (size_t i = 0; i < v.size(); ++i)
  {
    out << v[i];
    if (i != last)
      out << ", ";
  }
  out << "}";
  return out;
}
     
int praveenUpdated(const char *config){

  // First, setup timing information so we can monitor the run
  time_t rawtime;
  std::cout << "-----------------------------------------------------------" << std::endl;
  GetTime(rawtime);
  std::cout << "-----------------------------------------------------------" << std::endl;

  //------------------------------------------------------------------------------------------
  //                                    Configure
  //------------------------------------------------------------------------------------------
  // Create object of the class ConfigReader
  // Parse the configuration file
  // Dump map on the console after parsing it
  ConfigReader* p = ConfigReader::getInstance();
  p->parseFile(config);
  std::cout << " Variables from configuration file: " << std::endl;
  p->dumpFileValues();
  std::cout << "-----------------------------------------------------------" << std::endl;

  // Get configuration variables and initiate the relevant ones
  int n = -1;   // How many files from the file list to run. Default: All (-1)
  std::string input_list = "";
  std::string location="";
  std::string tag="";

  // Access corresponding parameter in the configuration file
  p->getValue("InputList", input_list);
  p->getValue("Location",  location);
  p->getValue("Tag",       tag);
  p->getValue("NFiles",    n);


  // Get constrants
  int hitplane = 2;
  double Efield = 0.50;      // kV/cm protoDUNE electric filed
  double betap = 0.212;      //(kV/cm)(g/cm^2)/MeV // taken from ArgoNeuT experiment
  double Rho = 1.3936;       // g/cm^3 (liquid argon density at temp 87.596 K 18.0 psia)
  const double rho = 1.3936; // g/cm^3 at 87.596 K // density of liquid argon// to calculate the density with different temperature
  double alpha = 0.93;       // parameter from ArgoNeuT experiment at 0.481kV/cm
  double Wion = 23.6e-6;     // work function of argon // parameter from ArgoNeuT experiment at 0.481kV/cm
  const int Z = 18;          // Atomic number of Argon
  const double A = 39.948;   // g/mol Atomic mass of Argon
  const double I = 188.0e-6; // Mev , I comes with ln factor of logarithmic  so
  // to make factor unitless number value is taken in MeV  ( 188.0 ev )

  const double K = 0.307;     // Mev.cm^2 / mol  // standard valur taken from pdg paper
  const double Mmu = 105.658; // Mev for Mu
  const double Me = 0.51;     // Mev for electron


 
  // Sort out the file tag by adding an underscore
  if(tag != "")
    tag = "_"+tag;

  //--------------------------------------------------------------------------------- ---------
  //                                    Initialise
  //--------------------------------------------------------------------------------- ---------

  // Setup TTree from input file list
  std::cout << " Reading files and filling tree" << std::endl;

  EventProcessor evtProc(allowed, input_list, n);
  evtProc.Initialize();

  // Now setup the tree and event objects to work with
  TChain *tree = evtProc.GetTree();
  anatree *evt = evtProc.GetEvents();

  
  // Start of analysis (loop over chain and events
  std::cout << " Running analysis..." << std::endl;

  double calib_factor = 0.00531941;

  // Then setup the histograms, counters and any other variables to add to
  // Setup histograms
  TH2F *h_reco_dQdx_RR = new TH2F("h_reco_dQdx_RR", "; Residual range [cm];dQ/dx [ADC/cm]", 200, 0, 200, 100, 0, 1000);
  TH2F *h_reco_dQdx_RR_e3ms = new TH2F("h_reco_dQdx_RR_e3ms", "; Residual range [cm];dQ/dx [ADC/cm]", 200, 0, 200, 100, 0, 1000);
  //TH2F *h_reco_dEdx_RR_uncal = new TH2F("h_reco_dEdx_RR_uncal", "; Residual range [cm]; dE/dx [MeV/cm]", 200, 0, 200, 100, 0, 5);
  TH2F *h_reco_dEdx_RR_cal = new TH2F("h_reco_dEdx_RR_cal", "; Residual range [cm]; dE/dx [MeV/cm]", 200, 0, 200, 100, 0, 5);

  TH2F *fhist_dedx = new TH2F("cal_dedx_vs_rr", ";Residual range [cm]; dE/dx [MeV/cm]", 200, 0, 200, 100, 0, 5);
  TH1F *dqdx_rat = new TH1F("dQdx_ratio", "plane_2; dQ/dx ratio;Number of entries", 100, 0, 10);
  TH1F *my_hist = new TH1F("my_hist", " ;Kinetic energy [MeV]; MPV dE/dx [MeV/cm]", 10, 0, 500);

  TH2F *fhist_dqdxcal = new TH2F("cal_dqdx_rr", ";Residual range [cm]; dQ/dx [ADC/cm]", 200, 0, 200, 100, 0, 1000);
  //TH2F *fhist_dqdxuncal = new TH2F("uncal_dqdx_rr", ";Residual range [cm]; dQ/dx [ADC/cm]", 200, 0, 200, 100, 0, 1000);
  TH1F *dqdx_uncal = new TH1F("uncal_dqdx", ";dQ/dx [ADC/cm]; Number of entries", 100, 0, 1e3);
  TH1F *dqdx_cal = new TH1F("cal_dqdx", "Calibrated dQ/dx;dQ/dx [ADC/cm]; Number of entries", 100, 0, 1e3);

  //Set up a couple of binned histos
  int nbin = 20;
  int binsize = 10;

  TH1D *dedx[nbin];

  for (int i = 0; i < nbin; ++i)
  {
    if (i == 0)
      dedx[i] = new TH1D(Form("dedx_%d", i), "; dE/dx [MeV/cm]; Number of entries", 100, 0.0, 10);

    if (i != 0)
      dedx[i] = new TH1D(Form("dedx_%d", i), "; dE/dx [MeV/cm]; Number of entries", 100, 0.0, 10);

    dedx[i]->SetLineColor(kBlack);
    dedx[i]->Sumw2(); // also store errors
    std::cout << "For i = " << i << std::endl;
  }

  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(1111);
  
  // Now loop over the events
  unsigned int nentries = tree->GetEntries();
  //Long64_t nbytes = 0, nb = 0;

  std::cout << " |";
  std::cout << "nentries = " << nentries << std::endl;
  for (Long64_t jentry = 0; jentry < nentries; jentry++) {
    tree->GetEntry(jentry);
    if(!evtProc.SelectEvent(evt)) continue;

    //nbytes += nb;
    if (jentry % 100 == 0)
      cout << jentry << "/" << nentries << endl;

    std::cout << "jentry = " << jentry << std::endl;

    bool found_stopping_muons = false;

    int muonsevent = 0;
    for (int iG4 = 0; iG4 < evt->geant_list_size; ++iG4)
    {
      if (evt->inTPCActive[iG4] == 1 && abs(evt->pdg[iG4]) == 13)
        muonsevent = 1;
    }

    if (!evt->inTPCActive[0])
      continue;

    int n_muons = 0;
    int n_muons_associated_with_reco_tracks = 0;
    int n_muons_associated_with_only_one_reco_tracks = 0;

    double dx_avg_per_event = 0;

    for (int iG4 = 0; iG4 < evt->geant_list_size; ++iG4)
    {
      TVector3 vtx(evt->StartPointx[iG4], evt->StartPointy[iG4], evt->StartPointz[iG4]);
      TVector3 end(evt->EndPointx[iG4], evt->EndPointy[iG4], evt->EndPointz[iG4]);

      if (!evt->inTPCActive[iG4])
        continue;

      TVector3 vtxAV(evt->StartPointx_tpcAV[iG4], evt->StartPointy_tpcAV[iG4], evt->StartPointz_tpcAV[iG4]);
      TVector3 endAV(evt->EndPointx_tpcAV[iG4], evt->EndPointy_tpcAV[iG4], evt->EndPointz_tpcAV[iG4]);


      float ilengthAV = (endAV - vtxAV).Mag();
      float length = (end - vtx).Mag();
      float totalL = evt->pathlen[iG4];
      float iEndPointx_tpcAV = evt->EndPointx_tpcAV[iG4];
      float iEndPointy_tpcAV = evt->EndPointy_tpcAV[iG4];
      float iEndPointz_tpcAV = evt->EndPointz_tpcAV[iG4];
      float iEndPointx = evt->EndPointx[iG4];
      float iEndPointy = evt->EndPointy[iG4];
      float iEndPointz = evt->EndPointz[iG4];
      int inTPCActive = evt->inTPCActive[iG4];
      int iTrackId = evt->TrackId[iG4];

      bool stoppingmuons = false;
      float dx_end = abs(endAV.X() - evt->EndPointx[iG4]);
      float dx_start = abs(vtxAV.X() - evt->StartPointx[iG4]);
      float dy_end = abs(endAV.Y() - evt->EndPointy[iG4]);
      float dy_start = abs(vtxAV.Y() - evt->StartPointy[iG4]);
      float dz_end = abs(endAV.Z() - evt->EndPointz[iG4]);
      float dz_start = abs(vtxAV.Z() - evt->StartPointz[iG4]);

      if (dx_end + dy_end + dz_end == 0 && (dx_start != 0 || dy_start != 0 || dz_start != 0))
        stoppingmuons = true;

      int id = evt->TrackId[iG4];
      int ipdg = evt->pdg[iG4];
      int iMother = evt->Mother[iG4];

      if (abs(ipdg) != 13) // for muons
        continue;

      if (evt->Mother[iG4] != 0)
        continue;

      if (!stoppingmuons)
        continue;

      n_muons++;

      found_stopping_muons = true;

      int nHits = evt->no_hits_stored;

      double costheta = GetCosTheta(2, vtxAV, endAV);
      double cosCollection = GetCosTheta(2, vtxAV, endAV);

      vector<float> res, dq, first5dq, last5dq;

      int n_trks_per_muons_three_wire_plane = 0;
      int n_trks_per_muons_one_wire_plane = 0;
      int n_trks_per_muons_wire_plane0 = 0;
      int n_trks_per_muons_wire_plane1 = 0;
      int n_trks_per_muons_wire_plane2 = 0;

      bool found_muons_with_track = false;

      int n_tracks_per_muons = 0;
      int itrack_matched = -999;


      vector<double> matched_track_numbers;

      for (int iTrk = 0; iTrk < evt->ntracks_pandoraTrack; ++iTrk)
      {

        int trueID1 = evt->TrackId[iG4];
        for (int iWire = 0; iWire < 3; ++iWire) //?
        {
          if (evt->trkidtruth_pandoraTrack[iTrk][iWire] == trueID1) //?
          {
            n_tracks_per_muons++;
            matched_track_numbers.push_back(iTrk);
          }
        }
      } //iTrk loop


     if (n_tracks_per_muons == 0)
      {
        matched_track_numbers.push_back(-999);
        continue;
      }


      int greatest_track_number = *max_element(matched_track_numbers.begin(), matched_track_numbers.end());
      int lowest_track_number = *min_element(matched_track_numbers.begin(), matched_track_numbers.end());

      for (int iTrk = 0; iTrk < evt->ntracks_pandoraTrack; ++iTrk)
      {
        int trueID1 = evt->TrackId[iG4];


        TVector3 startVtx(evt->trkstartx_pandoraTrack[iTrk],
                          evt->trkstarty_pandoraTrack[iTrk],
                          evt->trkstartz_pandoraTrack[iTrk]);
        TVector3 endVtx(evt->trkendx_pandoraTrack[iTrk],
                        evt->trkendy_pandoraTrack[iTrk],
                        evt->trkendz_pandoraTrack[iTrk]);

        int currHits = -999;
        std::vector<int> hitsOnPlane(3, 0);
        for (int iPlane = 0; iPlane < 3; ++iPlane)
          float eng = -1.;

        int trueID = evt->TrackId[iG4];

        //wire

        for (int iWire = 0; iWire < 3; ++iWire)
        {
          if (abs(evt->trkpdgtruth_pandoraTrack[iTrk][iWire]) != 13)
            continue;

          if (evt->trkidtruth_pandoraTrack[iTrk][iWire] != trueID1)
            continue;

          int nHitsR = evt->ntrkhits_pandoraTrack[iTrk][iWire];
          if (nHitsR <= 0)
          {
            continue;
          }

          float energy = evt->trkke_pandoraTrack[iTrk][iWire] / 1000.; // GeV

          for (int j = 0; j < evt->ntrkhits_pandoraTrack[iTrk][iWire]; j++) // loop over pandora track hits
          {
            res.push_back(evt->trkresrg_pandoraTrack[iTrk][iWire][j]);
            dq.push_back(evt->trkdqdx_pandoraTrack[iTrk][iWire][j]);
          } // ntrkhits_pandoraTrack
          //end of buffer filling


          if (res.size() == 0)
            continue;           //removing empty tracks avoids seg fault


          int siz1 = evt->ntrkhits_pandoraTrack[iTrk][iWire];
          float max = *max_element(res.begin(), res.end());
          float min = *min_element(res.begin(), res.end());
          int siz = res.size();

          for (int iHit = 1; iHit < nHitsR - 1; ++iHit)
          {
            // Get the location of the current and following hits to determine the pitch
            TVector3 trkXYZ(evt->trkxyz_pandoraTrack[iTrk][iWire][iHit][0],
                            evt->trkxyz_pandoraTrack[iTrk][iWire][iHit][1],
                            evt->trkxyz_pandoraTrack[iTrk][iWire][iHit][2]);
            TVector3 nextXYZ(evt->trkxyz_pandoraTrack[iTrk][iWire][iHit + 1][0],
                             evt->trkxyz_pandoraTrack[iTrk][iWire][iHit + 1][1],
                             evt->trkxyz_pandoraTrack[iTrk][iWire][iHit + 1][2]);

            float x = trkXYZ.X();
            float t = x * evtProc.kXtoT;
            double dp = GetHitPitch(iWire, trkXYZ, nextXYZ);
            double cosDrift = GetCosDrift(trkXYZ, nextXYZ);

            int tpc = evtProc.WhichTPC(x) + 1;

            float hit_width = evt->hit_endT[iHit] - evt->hit_startT[iHit]; // In ticks; each tick =500 ns.
            float widthT = hit_width * 0.5e-6;                   // To s
            float widthX = widthT / static_cast<float>(evtProc.kXtoT);   // To cm

            float dx = (-1 + 2 * (tpc % 2)) * (x - evtProc.APA_X_POSITIONS[tpc / 2]); // calculate the distance(positive) between hitx and apa plane
            float dt = dx * evtProc.kXtoT;
            float corr = TMath::Exp(-dt / 2.88);
            float ecorr = TMath::Exp(-dt / 2.88);

            float hitWidth = evt->hit_endT[iHit] - evt->hit_startT[iHit];

            // std::cout << "hitWidth " << hitWidth << std::endl;
            if (iWire == 2)
            {

              float corrected_dq_dx = evt->trkdqdx_pandoraTrack[iTrk][iWire][iHit] / ecorr;

              float corrected_dE_dx = evt->trkdedx_pandoraTrack[iTrk][iWire][iHit] / ecorr;
	      std::cout << "corrected_dE_dx = " << corrected_dE_dx << std::endl;
              float scaled_corrected_dq_dx = float(corrected_dq_dx) / calib_factor;

             // float cal_de_dx = Dedx(scaled_corrected_dq_dx, Efield);
              float cal_de_dx = (exp(scaled_corrected_dq_dx * (betap / (Rho * Efield) * Wion)) - alpha) / (betap / (Rho * Efield));
              std::cout << "cal_de_dx = " << cal_de_dx << std::endl;
              int bin = int(evt->trkresrg_pandoraTrack[iTrk][iWire][iHit]) / binsize; //residual range
              fhist_dedx->Fill(evt->trkresrg_pandoraTrack[iTrk][iWire][iHit], cal_de_dx);
              fhist_dqdxcal->Fill(evt->trkresrg_pandoraTrack[iTrk][iWire][iHit], corrected_dq_dx);

              dqdx_cal->Fill(corrected_dq_dx / calib_factor); // scaled_corrected_dq_dx
              dqdx_uncal->Fill(evt->trkdqdx_pandoraTrack[iTrk][iWire][iHit] / calib_factor);

              h_reco_dQdx_RR->Fill(evt->trkresrg_pandoraTrack[iTrk][iWire][iHit], corrected_dq_dx);
              h_reco_dEdx_RR_cal->Fill(evt->trkresrg_pandoraTrack[iTrk][iWire][iHit], cal_de_dx);

              if (bin < nbin)
              {
                dedx[bin]->Fill(cal_de_dx);
              } // bin <40
            }
          } //iHit
        } //iWire
      } //iTrk
    } //G4 loop

  }// Event loop
  std::cout << " --- 100 % --- |" << std::endl;

  ofstream myfile1;

  double range_new[nbin];
  double dedxtheory[nbin];
  double keng;
  double krange;

  vector<double> range;
  vector<double> erange;
  vector<double> energy;
  vector<double> eenergy;
  float range_measured[nbin];
  float energy_measured[nbin];

  vector<double> Thke;
  vector<double> Theng;
  vector<double> Theke;
  vector<double> Theeng;

  vector<double> Thaveng;

  vector<double> chi_2_denominator;
  vector<double> chi_2_numerator;
  int dof = 0;


  myfile1.open("muon_mpv.txt");

  ifstream kerr;
  kerr.open("ke_rr.txt");

  ifstream dedxrr;
  dedxrr.open("dedx_rr.txt");

  for (int i = 0; i < nbin; i++)

  {

    std::cout << "Fitting bin ************** " << i << std::endl;
    double k1, d1, r1, r2;
    kerr >> k1 >> r1; //read in the values from Praveen's files
    dedxrr >> d1 >> r2;
    Thke.push_back(k1);
    Theng.push_back(r1);
    myfile1 << r1 << "  " << d1 << endl;
    range_new[i] = r1;
    dedxtheory[i] = d1;
    Theke.push_back(0);
    Theeng.push_back(0);
    std::cout << "Values from txt files" << std::endl;
    std::cout << d1 << "  " << r1 << endl;
    std::cout << k1 << "  " << r2 << endl;
    std::cout << "-------------------------------" << std::endl;

    TCanvas *d[i];
    d[i] = new TCanvas(Form("d_%d", i), Form("d_%d", i));

    //Double_t fr[2];
    Double_t fp[4], fpe[4];

    Double_t chisqr;
    Int_t ndf;
    //TF1 *fitsnr = this->langaufit(dedx[i], fr, sv, pllo, plhi, fp, fpe, &chisqr, &ndf);
    //langaufit
    Int_t i2;
    Char_t FunName[100];

    sprintf(FunName, "Fitfcn_%s", dedx[i]->GetName());

    TF1 *ffitold = (TF1 *)gROOT->GetListOfFunctions()->FindObject(FunName);
    if (ffitold)
      delete ffitold;

    double norm = dedx[i]->GetEntries() * dedx[i]->GetBinWidth(1);
    int maxbin          = dedx[i]->GetMaximumBin();
    double maxloc       = dedx[i]->GetBinCenter(maxbin);
    double sv[4] = {0.2*maxloc, maxloc, norm, 0.2*maxloc}; // starting values for parameters: Landau scale, Landau MPV, Norm, Gauss sigma
    double minR = dedx[i]->GetXaxis()->GetXmin();
    double maxR = dedx[i]->GetXaxis()->GetXmax();
    double nBinsFromPeak = 6;
    if(maxbin-nBinsFromPeak > 1)
      minR = dedx[i]->GetBinCenter(maxbin-nBinsFromPeak);
    if(maxbin+nBinsFromPeak < dedx[i]->GetNbinsX())
      maxR = dedx[i]->GetBinCenter(maxbin+nBinsFromPeak);


    TF1 *ffit;
    ffit = new TF1(FunName, langaufun, minR, maxR, 4);
    ffit->SetParameters(sv);
    ffit->SetParNames("Width", "MPV", "Area", "GSigma"); 

    //for (i2 = 0; i2 < 4; i2++)
   // {
      //std::cout << "low = " << pllo[i2] << "  and high = " << plhi[i2] << std::endl;
     // ffit->SetParLimits(i2, 1, 5);
      //ffit->SetParLimits(i2, 0, 4);
    //}
    
    if (dedx[i]->GetEntries() == 0)
	continue;

    auto result = dedx[i]->Fit(ffit, "QSMR", "");
    dedx[i]->Print();

    ffit->GetParameters(fp); // obtain fit parameters
    

    for (i2 = 0; i2 < 4; i2++)
    {
      fpe[i2] = ffit->GetParError(i2); // obtain fit parameter errors
    }

    chisqr = ffit->GetChisquare(); // obtain chi^2
    ndf = ffit->GetNDF();          // obt

    std::cout << "chisqr = " << chisqr << std::endl;
    std::cout << "ndf = " << ndf << std::endl;
    //----------------------
    //
    double mpv = result->Parameter(1);
    mpv = ffit->GetMaximumX(result->Parameter(1), result->Parameter(1) + result->Parameter(3));

    std::cout << "------------------------------------" << std::endl;
    for(unsigned int p = 0; p < result->NPar(); ++p){
      std::cout << " " << result->ParName(p) << " : " << result->Parameter(p) << std::endl;
    }
    std::cout << "------------------------------------" << std::endl;
    std::cout << " Peak pitch: " << maxloc << std::endl;
    std::cout << "------------------------------------" << std::endl;




    //-----------------
    
    ffit->SetLineColor(kRed);

    dedx[i]->SetMarkerStyle(8);
    dedx[i]->SetMarkerSize(0.5);
    dedx[i]->SetStats(1);
    dedx[i]->SetMaximum(dedx[i]->GetMaximum() * 1.25);

    gStyle->SetOptFit(1111);

    TPaveStats *ps = (TPaveStats *)dedx[i]->FindObject("stats");
    ps->SetX1NDC(0.6);
    ps->SetX2NDC(0.85);
    ps->SetY1NDC(0.65);
    ps->SetY2NDC(0.9);

    d[i]->Write();

    d[i]->SaveAs(Form("d_%d.pdf", i));

    d[i]->Close();

    std::cout << "ABOVE PASSES" << std::endl;
    if (gMinuit && dedx[i]->GetEntries() > 100)
    {
      std::cout << "PASS 100" << std::endl;
      TString test = gMinuit->fCstatu.Data();
      if (ffit->GetNDF() != 0)
      {
        std::cout << "PASS NDF" << std::endl;
        std::cout << "ffit->GetParError(1) " << ffit->GetParError(1) << std::endl;
        std::cout << "ffit->GetChisquare() / ffit->GetNDF() " << ffit->GetChisquare() / ffit->GetNDF() << std::endl;
        std::cout << "test " << test << std::endl;
        if ((test.EqualTo("CONVERGED ") or test.EqualTo("OK ") ) && (ffit->GetParError(1) < 1000) && ((ffit->GetChisquare() / ffit->GetNDF() < 10)))
        {
           std::cout << "PASS CONVERGE" << std::endl;
          range.push_back(r1);
          erange.push_back(0);
          range_measured[i] = r1;
          energy_measured[i] = ffit->GetParameter(1);
          energy.push_back(ffit->GetParameter(1));
          eenergy.push_back(ffit->GetParError(1));
	  std::cout << "energy measured in bin " <<ffit->GetParameter(1) <<std::endl;
          std::cout << "range measured in bin " << r1 <<std::endl;
        }
      }
    }
  }

  kerr.close();
  dedxrr.close();

  double sum = 0;

  TGraphErrors *th_mpv_graph = new TGraphErrors(Thke.size(), &Thke[0], &Theng[0], &Theke[0], &Theeng[0]);
  th_mpv_graph->SetMarkerStyle(7);
  th_mpv_graph->SetMarkerColor(kBlack);
  th_mpv_graph->SetLineColor(kBlack);

  TGraphErrors *th_ave_graph = new TGraphErrors(Thke.size(), &Thke[0], &Thaveng[0], &Theke[0], &Theeng[0]);
  th_ave_graph->SetMarkerStyle(21);
  th_ave_graph->SetMarkerColor(kGreen);

  TGraphErrors *e_graph = new TGraphErrors(range.size(), &range[0], &energy[0], &erange[0], &eenergy[0]);
  e_graph->GetXaxis()->SetTitle("Kinetic energy (MeV/cm)");
  e_graph->GetYaxis()->SetTitle("MPV Energy (MeV/cm)");
  e_graph->SetTitle("Kinetic energy Vs MPV energy");
  e_graph->SetMarkerStyle(7);
  e_graph->SetMarkerColor(kRed);
  e_graph->SetLineColor(kRed);

  TText *lable_1 = new TText();
  lable_1->SetNDC();
  lable_1->SetTextFont(1);
  lable_1->SetTextColor(kRed);
  lable_1->SetTextSize(0.03);
  lable_1->SetTextAlign(22);
  lable_1->SetTextAngle(0);

  TCanvas *c1 = new TCanvas();

  TGraph *theoretical = new TGraph(nbin, range_new, dedxtheory);
  std::cout << "range_new " << range_new[4] << "    dedxtheory " << dedxtheory[4] << std::endl;

  c1->cd();
  my_hist->Draw();

  e_graph->Draw("samePE");
  th_mpv_graph->Draw("samePE");

  auto legend0 = new TLegend(0.65, 0.7, 0.9, 0.85);

  legend0->SetTextSize(0.04);
  legend0->AddEntry(th_mpv_graph, "Theory", "lep");
  legend0->AddEntry(e_graph, "Measured", "lep");

  legend0->Draw("same");

  TLatex latex;
  latex.SetTextColor(kRed);
  latex.SetTextColor(kBlack);
  latex.SetTextSize(0.035);
  latex.SetTextAlign(13);

  c1->cd();
  c1->Draw();
  c1->Write();

  //c1->SaveAs(Form("MPV_ke.pdf"));
  c1->SaveAs(Form("MPV_ke.png"));

  c1->Clear();

  fhist_dedx->Write();
  fhist_dedx->Draw("COLZ");
  //c1->SaveAs(Form("fhist_dedx.pdf"));
  c1->SaveAs(Form("fhist_dedx.png"));
  c1->Clear();

  c1->Close();

  myfile1.close();

  TCanvas *cdedx = new TCanvas();
  cdedx->cd();

  fhist_dedx->Draw("COLZ");

  //cdedx->SaveAs(Form("fhist_dedx1.pdf"));
  cdedx->SaveAs(Form("fhist_dedx1.png"));

  theoretical->SetMarkerColor(kMagenta);;
  theoretical->SetMarkerColor(kMagenta);

  theoretical->SetMarkerStyle(1);
  theoretical->Draw("sameCP");

  TLegend *legdedx = new TLegend(0.2, 0.68, 0.83, 0.82);

  legdedx->SetTextSize(0.03);
  legdedx->AddEntry(theoretical, "Theoretical most probable value (Landau-Vavilov theory)", "LP");

  TGraph *gr3 = new TGraph(nbin, range_measured, energy_measured);

  gr3->SetLineColor(kRed);
  gr3->SetMarkerStyle(1);
  gr3->SetMarkerColor(kRed);
  gr3->Draw("sameCP");

  legdedx->AddEntry(gr3, "Measured using Gauss-Landau fit", "LP");
  legdedx->Draw("same");

  cdedx->Draw();
  cdedx->Write();

  //cdedx->SaveAs(Form("cdedx.pdf"));
  cdedx->SaveAs(Form("cdedx.png"));
  cdedx->Clear();
  cdedx->Close();

  TCanvas *c2 = new TCanvas();

  c2->cd();
  fhist_dqdxcal->Write();
  fhist_dqdxcal->Draw("COLZ");

  //c2->SaveAs(Form("fhist_dqdxcal.pdf"));
  c2->SaveAs(Form("fhist_dqdxcal.png"));
  c2->Clear();

  //fhist_dqdxuncal->Write();
  //fhist_dqdxuncal->Draw("COLZ");
  //c2->SaveAs(Form("fhist_dqdxuncal.pdf"));
  //c2->SaveAs(Form("fhist_dqdxuncal.png"));
  //c2->Clear();

  dqdx_cal->Write();
  dqdx_cal->Draw();
  //c2->SaveAs(Form("dqdx_cal.pdf"));
  c2->SaveAs(Form("dqdx_cal.png"));
  c2->Clear();

  dqdx_uncal->Write();
  dqdx_uncal->Draw();
  //c2->SaveAs(Form("dqdx_uncal.pdf"));
  c2->SaveAs(Form("dqdx_uncal.png"));
  c2->Clear();

  gr3->Write("measured_dedxrr");
  gr3->SetTitle(" ; Residual range [cm]; dE/dx [MeV/cm]");

  //gr3->Draw("ACP");
  gr3->Draw();

  //gPad->SaveAs(Form("measured_dEdx_rr.pdf"));
  gPad->SaveAs(Form("measured_dEdx_rr.png"));
  c2->Clear();

  theoretical->Write("theoretical");

  theoretical->SetTitle(" ; Residual range [cm]; dE/dx [MeV/cm]");

  theoretical->Draw("ACP");

  //gPad->SaveAs(Form("theory_dEdx_rr.pdf"));
  gPad->SaveAs(Form("theory_dEdx_rr.png"));
  c2->Clear();

  std::cout << "************************* Calibration.C has ended ***************************" << std::endl;

  gStyle->SetPalette(1, 0);
  gStyle->SetNumberContours(64);
  fhist_dqdxcal->Draw("colz");
  c2->SaveAs(Form("fhist_dqdxcal_last.png"));
  c2->Clear();
  c2->Close();

  TCanvas *c3 = new TCanvas();
  h_reco_dQdx_RR->Draw("COLZ");
  h_reco_dQdx_RR->Write(" h_reco_dQdx_RR");

  //c3->SaveAs("reco_dqdx_rr.pdf");
  c3->SaveAs("reco_dqdx_rr.png");
  c3->Write("reco_dQdx_RR");
  c3->Clear();

  //h_reco_dEdx_RR_uncal->Write();
  //h_reco_dEdx_RR_uncal->Draw("COLZ");

  //c3->SaveAs("reco_dEdx_rr_uncal.pdf");
  //c3->SaveAs("reco_dEdx_rr_uncal.png");
  //c3->Clear();

  h_reco_dEdx_RR_cal->Write();
  h_reco_dEdx_RR_cal->Draw("COLZ");

  //c3->SaveAs("reco_dEdx_rr_cal.pdf");
  c3->SaveAs("reco_dEdx_rr_cal.png");
  c3->Clear();

  fhist_dedx->Write();
  fhist_dedx->Draw("COLZ");
  //c3->SaveAs("fhist_dedx.pdf");
  c3->SaveAs("fhist_dedx.png");
  c3->Clear();

  c3->Close();

  gStyle->SetPalette(1, 0);
  gStyle->SetNumberContours(64);

  //tree->Print();

  // End of script
  std::cout << " ...finished analysis" << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;
  time_t rawtime_end;
  GetTime(rawtime_end);
  std::cout << "-----------------------------------------------------------" << std::endl;
  GetTotalTime(rawtime, rawtime_end);
  std::cout << "-----------------------------------------------------------" << std::endl;
  
  return 0;
}
