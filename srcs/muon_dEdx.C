////////////////////////// for stopping muons  //////////////
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
//#include <function.h>

using namespace std;

using std::ostream;

//A small function which seems to be used to print comments
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

//Make legends for plot
TLegend *MakeLegend(float left = 0.7, float bottom = 0.5, float right = 0.9, float top = 0.85)
{
  auto leg = new TLegend(left, bottom, right, top);
  leg->SetFillStyle(0); // unfortunately can't set this in TStyle :(
  return leg;
}

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
//
const double K = 0.307;     // Mev.cm^2 / mol  // standard valur taken from pdg paper
const double Mmu = 105.658; // Mev for Mu
const double Me = 0.51;     // Mev for electron

//function definition//

//Landau Gaussian function
Double_t langaufun(Double_t *x, Double_t *par)
{
  Double_t invsq2pi = 0.398942280401; // Control constants
  Double_t mpshift = -0.22278298;
  Double_t np = 500.0;
  Double_t sc = 5.0; // convolution extends to +-sc Gaussian sigmas
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow, xupp;
  Double_t step;
  Double_t i;

  // mpc = par[1]- mpshift * par[0];
  mpc = par[1];
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];
  step = (xupp - xlow) / np;

  for (i = 1.0; i <= np / 2; i++)
  {
    xx = xlow + (i - .5) * step;
    fland = TMath::Landau(xx, mpc, par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0], xx, par[3]);
    xx = xupp - (i - .5) * step;
    fland = TMath::Landau(xx, mpc, par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0], xx, par[3]);
  }

  return (par[2] * step * sum * invsq2pi / par[3]);
}

//Laudau Gaussian Fit, gets fit parameters, chi2
TF1 *langaufit(TH1D *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF)
{
  Int_t i;
  Char_t FunName[100];

  sprintf(FunName, "Fitfcn_%s", his->GetName()); // edited  complitation warning

  TF1 *ffitold = (TF1 *)gROOT->GetListOfFunctions()->FindObject(FunName);
  if (ffitold)
    delete ffitold;

  TF1 *ffit = new TF1(FunName, langaufun, fitrange[0], fitrange[1], 4);
  ffit->SetParameters(startvalues);
  ffit->SetParNames("Width", "MPV", "Area", "GSigma");

  for (i = 0; i < 4; i++)
  {
    ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
  }

  his->Fit(FunName, "RBO"); // fit within specified range, use ParLimits, do not plot ////////////////

  ffit->GetParameters(fitparams); // obtain fit parameters
  for (i = 0; i < 4; i++)
  {
    fiterrors[i] = ffit->GetParError(i); // obtain fit parameter errors
  }

  ChiSqr[0] = ffit->GetChisquare(); // obtain chi^2
  NDF[0] = ffit->GetNDF();          // obtain ndf

  return (ffit); // return fit function
}

Int_t langaupro(Double_t *params, Double_t &maxx, Double_t &FWHM)
{
  Double_t p, x, fy, fxr, fxl;
  Double_t step;
  Double_t l, lold;
  Int_t i = 0;
  Int_t MAXCALLS = 10000;

  // Search for maximum
  

  p = params[1] - 0.1 * params[0];
  step = 0.05 * params[0];
  lold = -2.0;
  l = -1.0;

  while ((l != lold) && (i < MAXCALLS))
  {
    i++;
    lold = l;
    x = p + step;
    l = langaufun(&x, params);
    if (l < lold)
      step = -step / 10;
    // step = -step / 4; //added
    p += step; 
  }

  if (i == MAXCALLS)
    return (-1);

  maxx = x;
  fy = l / 2;

  // Search for right x location of fy

  p = maxx + params[0];
  step = params[0];
  lold = -2.0;
  l = -1e300;
  i = 0;

  while ((l != lold) && (i < MAXCALLS))
  {
    i++;
    lold = l;
    x = p + step;
    l = TMath::Abs(langaufun(&x, params) - fy);
    if (l > lold)
      step = -step / 10;
     // step = -step / 4;
    p += step;
  }

  if (i == MAXCALLS)
    return (-2);

  fxr = x;

  // Search for left x location of fy


  step = -params[0];
  lold = -2.0;
  l = -1e300;
  i = 0;

  while ((l != lold) && (i < MAXCALLS))
  {
    i++;
    lold = l;
    x = p + step;
    l = TMath::Abs(langaufun(&x, params) - fy);
    if (l > lold)
      step = -step / 10;
    // step = -step / 4; //added
    p += step;
  }

  if (i == MAXCALLS)
    return (-3);

  fxl = x;
  FWHM = fxr - fxl;
  return (0);
}

//Function Def

//Appears to be a function to calculate dEdx from dQdx using modified box model
Float_t Dedx(float dqdx, float Ef)
{
  return (exp(dqdx * (betap / (Rho * Ef) * Wion)) - alpha) / (betap / (Rho * Ef));
}

void muon_dEdx::Loop()
{
  /////////////////////////////// Read only selected branch //////////////


  fChain->SetBranchStatus("*", 0);   // disable all branches
  fChain->SetBranchStatus("run", 1); // activate branchname
  fChain->SetBranchStatus("event", 1);
  fChain->SetBranchStatus("geant_list_size", 1);
  fChain->SetBranchStatus("inTPCActive", 1);
  fChain->SetBranchStatus("TrackId", 1);
  fChain->SetBranchStatus("pdg", 1);
  fChain->SetBranchStatus("Mother", 1);
  fChain->SetBranchStatus("ntracks_pandoraTrack", 1);
  fChain->SetBranchStatus("trkId_pandoraTrack", 1);
  fChain->SetBranchStatus("trkidtruth_pandoraTrack", 1);
  fChain->SetBranchStatus("trkpdgtruth_pandoraTrack", 1);
  fChain->SetBranchStatus("trkg4id_pandoraTrack", 1);
  fChain->SetBranchStatus("ntrkhits_pandoraTrack", 1);
  fChain->SetBranchStatus("trkresrg_pandoraTrack", 1);
  fChain->SetBranchStatus("trkdqdx_pandoraTrack", 1);
  fChain->SetBranchStatus("trkdedx_pandoraTrack", 1);
  fChain->SetBranchStatus("trkke_pandoraTrack", 1);
  fChain->SetBranchStatus("trkxyz_pandoraTrack", 1);
  fChain->SetBranchStatus("trkstartx_pandoraTrack", 1);
  fChain->SetBranchStatus("trkstarty_pandoraTrack", 1);
  fChain->SetBranchStatus("trkstartz_pandoraTrack", 1);
  fChain->SetBranchStatus("trkendx_pandoraTrack", 1);
  fChain->SetBranchStatus("trkendy_pandoraTrack", 1);
  fChain->SetBranchStatus("trkendz_pandoraTrack", 1);
  fChain->SetBranchStatus("trklen_pandoraTrack", 1);
  fChain->SetBranchStatus("no_hits_stored", 1);
  fChain->SetBranchStatus("hit_tpc", 1);
  fChain->SetBranchStatus("hit_plane", 1);
  fChain->SetBranchStatus("hit_charge", 1);
  fChain->SetBranchStatus("hit_energy", 1);
  fChain->SetBranchStatus("hit_nelec", 1);
  fChain->SetBranchStatus("hit_trueX", 1);
  fChain->SetBranchStatus("hit_trkid", 1);
  fChain->SetBranchStatus("hit_peakT", 1);
  fChain->SetBranchStatus("hit_startT", 1);
  fChain->SetBranchStatus("hit_endT", 1);
  fChain->SetBranchStatus("StartPointx", 1);
  fChain->SetBranchStatus("StartPointy", 1);
  fChain->SetBranchStatus("StartPointz", 1);
  fChain->SetBranchStatus("StartE_tpcAV", 1);
  fChain->SetBranchStatus("EndPointx", 1);
  fChain->SetBranchStatus("EndPointy", 1);
  fChain->SetBranchStatus("EndPointz", 1);
  fChain->SetBranchStatus("StartPointx_tpcAV", 1);
  fChain->SetBranchStatus("StartPointy_tpcAV", 1);
  fChain->SetBranchStatus("StartPointz_tpcAV", 1);
  fChain->SetBranchStatus("EndPointx_tpcAV", 1);
  fChain->SetBranchStatus("EndPointy_tpcAV", 1);
  fChain->SetBranchStatus("EndPointz_tpcAV", 1);
  fChain->SetBranchStatus("NumberDaughters", 1);
  fChain->SetBranchStatus("P", 1);
  fChain->SetBranchStatus("Eng", 1);
  fChain->SetBranchStatus("EndE", 1);
  fChain->SetBranchStatus("pathlen", 1);

  /////////////////////////////////////////////

  double calib_factor = 0.00531941;

  std::cout << "******************************* Calibration.C is running *******************************" << std::endl;

  // TH1F *h_energy_at_TPC = new TH1F("h_energy_at_TPC", "; Energy [GeV]; Number of events", 60, 0, 30);
  TH2F *h_reco_dQdx_RR = new TH2F("h_reco_dQdx_RR", "; Residual range [cm];dQ/dx [ADC/cm]", 200, 0, 200, 100, 0, 1000);
  TH2F *h_reco_dQdx_RR_e3ms = new TH2F("h_reco_dQdx_RR_e3ms", "; Residual range [cm];dQ/dx [ADC/cm]", 200, 0, 200, 100, 0, 1000);
  TH2F *h_reco_dEdx_RR_uncal = new TH2F("h_reco_dEdx_RR_uncal", "; Residual range [cm]; dE/dx [MeV/cm]", 200, 0, 200, 100, 0, 5);
  TH2F *h_reco_dEdx_RR_cal = new TH2F("h_reco_dEdx_RR_cal", "; Residual range [cm]; dE/dx [MeV/cm]", 200, 0, 200, 100, 0, 5);

  TH2F *fhist_dedx = new TH2F("cal_dedx_vs_rr", ";Residual range [cm]; dE/dx [MeV/cm]", 200, 0, 200, 100, 0, 5);
  TH1F *dqdx_rat = new TH1F("dQdx_ratio", "plane_2; dQ/dx ratio;Number of entries", 100, 0, 10);
  TH1F *my_hist = new TH1F("my_hist", " ;Kinetic energy [MeV]; MPV dE/dx [MeV/cm]", 10, 0, 500);

  TH2F *fhist_dqdxcal = new TH2F("cal_dqdx_rr", ";Residual range [cm]; dQ/dx [ADC/cm]", 200, 0, 200, 100, 0, 1000);
  TH2F *fhist_dqdxuncal = new TH2F("uncal_dqdx_rr", ";Residual range [cm]; dQ/dx [ADC/cm]", 200, 0, 200, 100, 0, 1000);
  TH1F *dqdx_uncal = new TH1F("uncal_dqdx", ";dQ/dx [ADC/cm]; Number of entries", 100, 0, 1e3);
  TH1F *dqdx_cal = new TH1F("cal_dqdx", "Calibrated dQ/dx;dQ/dx [ADC/cm]; Number of entries", 100, 0, 1e3);

  //////////////////////////////
  //
  //  //////////////////
  //
  //    // my_hist->SetFillColor(kCyan);
  //      // my_hist->SetLineColor(kCyan);
  //        // my_hist->SetBinContent(6, 10);
  //          // my_hist->SetBinContent(7, 10);
  //            // my_hist->SetBinContent(8, 10);
  //              // my_hist->SetBinContent(9, 10);
  //
  

  ofstream myfile;
  myfile.open("trkinfo.txt");
  int nbin = 20;
  int binsize = 10;

  TH1D *dedx[nbin];
  // TCanvas *d[nbin];

  for (int i = 0; i < nbin; ++i)
  {
    if (i == 0)
      dedx[i] = new TH1D(Form("dedx_%d", i), "; dE/dx [MeV/cm]; Number of entries", 300, 0.0, 15);

    if (i != 0)
      dedx[i] = new TH1D(Form("dedx_%d", i), "; dE/dx [MeV/cm]; Number of entries", 200, 0.0, 10);

    dedx[i]->SetLineColor(kBlack);
    dedx[i]->Sumw2(); // also store errors
  }

  gStyle->SetOptStat(1111);
  // gStyle->SetOptFit(111);
  gStyle->SetOptFit(1111);

  if (fChain == 0)
    return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;

  //EVENT LOOP
  for (Long64_t jentry = 0; jentry < nentries; jentry++)
  // for (Long64_t jentry = 0; jentry < 2; jentry++)
  {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0)
      break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;
    if (jentry % 100 == 0)
      cout << jentry << "/" << nentries << endl;

    bool found_stopping_muons = false;

    int muonsevent = 0;
    for (int iG4 = 0; iG4 < geant_list_size; ++iG4)
    {
      if (inTPCActive[iG4] == 1 && abs(pdg[iG4]) == 13)
        muonsevent = 1;
    }

    if (!inTPCActive[0])
      continue;

    int n_muons = 0;
    int n_muons_associated_with_reco_tracks = 0;
    int n_muons_associated_with_only_one_reco_tracks = 0;

    double dx_avg_per_event = 0;

    //G4 LOOP
    

    for (int iG4 = 0; iG4 < geant_list_size; ++iG4)
    {

      TVector3 vtx(StartPointx[iG4], StartPointy[iG4], StartPointz[iG4]);
      TVector3 end(EndPointx[iG4], EndPointy[iG4], EndPointz[iG4]);

      // Check the particle enters the TPC volume
      if (!inTPCActive[iG4])
        continue;

      TVector3 vtxAV(StartPointx_tpcAV[iG4], StartPointy_tpcAV[iG4], StartPointz_tpcAV[iG4]);
      TVector3 endAV(EndPointx_tpcAV[iG4], EndPointy_tpcAV[iG4], EndPointz_tpcAV[iG4]);

      // int ipdg = pdg[iG4];
     

      float ilengthAV = (endAV - vtxAV).Mag();
      float length = (end - vtx).Mag();
      float totalL = pathlen[iG4];
      float iEndPointx_tpcAV = EndPointx_tpcAV[iG4];
      float iEndPointy_tpcAV = EndPointy_tpcAV[iG4];
      float iEndPointz_tpcAV = EndPointz_tpcAV[iG4];
      float iEndPointx = EndPointx[iG4];
      float iEndPointy = EndPointy[iG4];
      float iEndPointz = EndPointz[iG4];
      int iinTPCActive = inTPCActive[iG4];
      int iTrackId = TrackId[iG4];

      bool stoppingmuons = false;
      float dx_end = abs(endAV.X() - EndPointx[iG4]);
      float dx_start = abs(vtxAV.X() - StartPointx[iG4]);
      float dy_end = abs(endAV.Y() - EndPointy[iG4]);
      float dy_start = abs(vtxAV.Y() - StartPointy[iG4]);
      float dz_end = abs(endAV.Z() - EndPointz[iG4]);
      float dz_start = abs(vtxAV.Z() - StartPointz[iG4]);

      if (dx_end + dy_end + dz_end == 0 && (dx_start != 0 || dy_start != 0 || dz_start != 0))
        stoppingmuons = true;

      int id = TrackId[iG4];
      int ipdg = pdg[iG4];
      int iMother = Mother[iG4];

      if (abs(ipdg) != 13) // for muons
        continue;

      if (Mother[iG4] != 0)
        continue;

      if (!stoppingmuons)
        continue;

      n_muons++;

      // stopping_muons++;
      found_stopping_muons = true;
      // int imuons =1;
      //


      int nHits = no_hits_stored;

      double costheta = GetCosTheta(2, vtxAV, endAV);
      double cosCollection = GetCosTheta(2, vtxAV, endAV);

      vector<float> res, dq, first5dq, last5dq;

      // cross check////////
      //


      int n_trks_per_muons_three_wire_plane = 0;
      int n_trks_per_muons_one_wire_plane = 0;
      int n_trks_per_muons_wire_plane0 = 0;
      int n_trks_per_muons_wire_plane1 = 0;
      int n_trks_per_muons_wire_plane2 = 0;

      bool found_muons_with_track = false;

      int n_tracks_per_muons = 0;
      int itrack_matched = -999;

      vector<double> matched_track_numbers;

      for (int iTrk = 0; iTrk < ntracks_pandoraTrack; ++iTrk)
      {

        // cout << "iTrk = "<<iTrk<<endl;
        int trueID1 = TrackId[iG4];
        for (int iWire = 0; iWire < 3; ++iWire) //?
        {

          if (trkidtruth_pandoraTrack[iTrk][iWire] == trueID1) //?
          {
            n_tracks_per_muons++;

            // cout << "iTrk = "<<a<<endl;
            matched_track_numbers.push_back(iTrk);
          }
        }
      }


      if (n_tracks_per_muons == 0)
      {
        matched_track_numbers.push_back(-999);
        continue;
      }

      int greatest_track_number = *max_element(matched_track_numbers.begin(), matched_track_numbers.end());
      int lowest_track_number = *min_element(matched_track_numbers.begin(), matched_track_numbers.end());

      for (int iTrk = 0; iTrk < ntracks_pandoraTrack; ++iTrk)
      {

        int trueID1 = TrackId[iG4];

        TVector3 startVtx(trkstartx_pandoraTrack[iTrk],
                          trkstarty_pandoraTrack[iTrk],
                          trkstartz_pandoraTrack[iTrk]);
        TVector3 endVtx(trkendx_pandoraTrack[iTrk],
                        trkendy_pandoraTrack[iTrk],
                        trkendz_pandoraTrack[iTrk]);

        int currHits = -999;
        std::vector<int> hitsOnPlane(3, 0);
        for (int iPlane = 0; iPlane < 3; ++iPlane)

          float eng = -1.;

        int trueID = TrackId[iG4];


        //Wire Loop

        for (int iWire = 0; iWire < 3; ++iWire)
        {

          // Only look at muons
          if (abs(trkpdgtruth_pandoraTrack[iTrk][iWire]) != 13)
            continue;

          if (trkidtruth_pandoraTrack[iTrk][iWire] != trueID1)
            continue;

          int nHitsR = ntrkhits_pandoraTrack[iTrk][iWire];
          if (nHitsR <= 0)
          {
            continue;
          }

          float energy = trkke_pandoraTrack[iTrk][iWire] / 1000.; // GeV
          for (int j = 0; j < ntrkhits_pandoraTrack[iTrk][iWire]; j++) // loop over pandora track hits
          {
            res.push_back(trkresrg_pandoraTrack[iTrk][iWire][j]);
            dq.push_back(trkdqdx_pandoraTrack[iTrk][iWire][j]);
          } // ntrkhits_pandoraTrack
          /**********************end of buffer filling*****************************/

          if (res.size() == 0)
            continue;
          /***********************removed empty tracks to avoid segmentation fault******************/


          int siz1 = ntrkhits_pandoraTrack[iTrk][iWire];
          float max = *max_element(res.begin(), res.end());
          float min = *min_element(res.begin(), res.end());
          int siz = res.size();

          for (int iHit = 1; iHit < nHitsR - 1; ++iHit)
          {

            // Get the location of the current and following hits to determine the pitch
            TVector3 trkXYZ(trkxyz_pandoraTrack[iTrk][iWire][iHit][0],
                            trkxyz_pandoraTrack[iTrk][iWire][iHit][1],
                            trkxyz_pandoraTrack[iTrk][iWire][iHit][2]);
            TVector3 nextXYZ(trkxyz_pandoraTrack[iTrk][iWire][iHit + 1][0],
                             trkxyz_pandoraTrack[iTrk][iWire][iHit + 1][1],
                             trkxyz_pandoraTrack[iTrk][iWire][iHit + 1][2]);

            float x = trkXYZ.X();
            float t = x * kXtoT;
            double dp = GetHitPitch(iWire, trkXYZ, nextXYZ);
            double cosDrift = GetCosDrift(trkXYZ, nextXYZ);

            int tpc = WhichTPC(x) + 1;


            float hit_width = hit_endT[iHit] - hit_startT[iHit]; // In ticks; each tick =500 ns.
            float widthT = hit_width * 0.5e-6;                   // To s
            float widthX = widthT / static_cast<float>(kXtoT);   // To cm

            float dx = (-1 + 2 * (tpc % 2)) * (x - APA_X_POSITIONS[tpc / 2]); // calculate the distance(positive) between hitx and apa plane
            float dt = dx * kXtoT;
            float corr = TMath::Exp(-dt / 2.88);
            float ecorr = TMath::Exp(-dt / 2.88);

            float hitWidth = hit_endT[iHit] - hit_startT[iHit];

            if (iWire == 2)
            {

              float corrected_dq_dx = trkdqdx_pandoraTrack[iTrk][iWire][iHit] / ecorr;

              float corrected_dE_dx = trkdedx_pandoraTrack[iTrk][iWire][iHit] / ecorr;
              float scaled_corrected_dq_dx = float(corrected_dq_dx) / calib_factor;

              float cal_de_dx = Dedx(scaled_corrected_dq_dx, Efield);

              int bin = int(trkresrg_pandoraTrack[iTrk][iWire][iHit]) / binsize;
              fhist_dedx->Fill(trkresrg_pandoraTrack[iTrk][iWire][iHit], cal_de_dx);
              fhist_dqdxcal->Fill(trkresrg_pandoraTrack[iTrk][iWire][iHit], corrected_dq_dx);

              dqdx_cal->Fill(corrected_dq_dx / calib_factor); // scaled_corrected_dq_dx
              dqdx_uncal->Fill(trkdqdx_pandoraTrack[iTrk][iWire][iHit] / calib_factor);

              h_reco_dQdx_RR->Fill(trkresrg_pandoraTrack[iTrk][iWire][iHit], corrected_dq_dx);
              h_reco_dEdx_RR_cal->Fill(trkresrg_pandoraTrack[iTrk][iWire][iHit], cal_de_dx);



              if (bin < nbin)
              {
                dedx[bin]->Fill(cal_de_dx);
              } // bin <40
            }

          } // loop over hits....

        } // iWire

      } // track loop

    } // G4 particle loop

  } // entries...........


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

    std::cout << "Fitting ************** " << i << std::endl;
    double k1, d1, r1, r2;
    kerr >> k1 >> r1;
    dedxrr >> d1 >> r2;
    Thke.push_back(k1);
    Theng.push_back(r1);
    myfile1 << r1 << "  " << d1 << endl;
    range_new[i] = r1;
    dedxtheory[i] = d1;
    Theke.push_back(0);
    Theeng.push_back(0);

    TCanvas *d[i];
    d[i] = new TCanvas(Form("d_%d", i), Form("d_%d", i));

    Double_t fr[2];
    Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];

    Double_t chisqr;
    Int_t ndf;
    TF1 *fitsnr = langaufit(dedx[i], fr, sv, pllo, plhi, fp, fpe, &chisqr, &ndf);
    fitsnr->SetLineColor(kRed);
    // fitsnr->SetMarkerStyle(2); //added
    dedx[i]->SetMarkerStyle(8);
    dedx[i]->SetMarkerSize(0.5);
    dedx[i]->SetStats(1);
    dedx[i]->SetMaximum(dedx[i]->GetMaximum() * 1.25);

    gStyle->SetOptFit(1111);

    // positon the stat box
        TPaveStats *ps = (TPaveStats *)dedx[i]->FindObject("stats");
    ps->SetX1NDC(0.6);
    ps->SetX2NDC(0.85);
    ps->SetY1NDC(0.65);
    ps->SetY2NDC(0.9);


    d[i]->Write();

    // d[i]->Clear();
    d[i]->SaveAs(Form("d_%d.pdf", i));
    // d[i]->Modified();
    // d[i]->Update();
    d[i]->Close();

    if (gMinuit && dedx[i]->GetEntries() > 100)
    {
      TString test = gMinuit->fCstatu.Data();
      if (fitsnr->GetNDF() != 0)
      {
        if (test.EqualTo("CONVERGED ") && (fitsnr->GetParError(1) < 1000) && ((fitsnr->GetChisquare() / fitsnr->GetNDF() < 10)) /* && (fitsnr->GetParError(0)<0.1)*/)
        {
          range.push_back(r1);
          erange.push_back(0);
          range_measured[i] = r1;
          energy_measured[i] = fitsnr->GetParameter(1);
          energy.push_back(fitsnr->GetParameter(1));
          eenergy.push_back(fitsnr->GetParError(1));
        }
      }
    }
  }


  kerr.close();
  dedxrr.close();

  // rootdedx->cd(); // change root directory to rootdedx



  double sum = 0;

  TGraphErrors *th_mpv_graph = new TGraphErrors(Thke.size(), &Thke[0], &Theng[0], &Theke[0], &Theeng[0]);
  th_mpv_graph->SetMarkerStyle(7);
  th_mpv_graph->SetMarkerColor(kBlack);
  th_mpv_graph->SetLineColor(kBlack);
  // th_mpv_graph->SetMarkerSize(.5);
  //
  //

  TGraphErrors *th_ave_graph = new TGraphErrors(Thke.size(), &Thke[0], &Thaveng[0], &Theke[0], &Theeng[0]);
  th_ave_graph->SetMarkerStyle(21);
  th_ave_graph->SetMarkerColor(kGreen);

  TGraphErrors *e_graph = new TGraphErrors(range.size(), &range[0], &energy[0], &erange[0], &eenergy[0]);
  // e_graph->SetMaximum(e_graph->GetMaximum() * 1.25);
  e_graph->GetXaxis()->SetTitle("Kinetic energy (MeV/cm)");
  e_graph->GetYaxis()->SetTitle("MPV Energy (MeV/cm)");
  e_graph->SetTitle("Kinetic energy Vs MPV energy");
  e_graph->SetMarkerStyle(7);
  e_graph->SetMarkerColor(kRed);
  e_graph->SetLineColor(kRed);
  // e_graph->SetMarkerSize(.5);
  

  TText *lable_1 = new TText();
  lable_1->SetNDC();
  lable_1->SetTextFont(1);
  lable_1->SetTextColor(kRed);
  lable_1->SetTextSize(0.03);
  lable_1->SetTextAlign(22);
  lable_1->SetTextAngle(0);

  TCanvas *c1 = new TCanvas();

  TGraph *theoretical = new TGraph(nbin, range_new, dedxtheory);

  c1->cd();
  my_hist->Draw();

  e_graph->Draw("samePE");
  th_mpv_graph->Draw("samePE");

  TLegend *legend0 = MakeLegend(0.65, 0.7, 0.9, 0.85);

  legend0->SetTextSize(0.04);
  legend0->AddEntry(th_mpv_graph, "Theory", "lep");
  legend0->AddEntry(e_graph, "Measured", "lep");

  legend0->Draw("same");

  TLatex latex;
  // latex.SetTextFont(1);
    latex.SetTextColor(kRed);
  latex.SetTextColor(kBlack);
  latex.SetTextSize(0.035);
  latex.SetTextAlign(13);

  c1->Draw();
  c1->Write();

  c1->SaveAs(Form("MPV_ke.pdf"));
  c1->SaveAs(Form("MPV_ke.png"));

  c1->Clear();

  /////

  fhist_dedx->Write();
  fhist_dedx->Draw("COLZ");
  c1->SaveAs(Form("fhist_dedx.pdf"));
  c1->SaveAs(Form("fhist_dedx.png"));
  c1->Clear();

  c1->Close();

  myfile1.close();

  TCanvas *cdedx = new TCanvas();
  cdedx->cd();

  fhist_dedx->Draw("COLZ");

  cdedx->SaveAs(Form("fhist_dedx1.pdf"));
  cdedx->SaveAs(Form("fhist_dedx1.png"));

  theoretical->SetMarkerColor(kBlack);
  theoretical->SetMarkerColor(kBlack);

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

  cdedx->SaveAs(Form("cdedx.pdf"));
  cdedx->SaveAs(Form("cdedx.png"));
  cdedx->Clear();
  cdedx->Close();

  TCanvas *c2 = new TCanvas();

  fhist_dqdxcal->Write();
  fhist_dqdxcal->Draw("COLZ");

  c2->SaveAs(Form("fhist_dqdxcal.pdf"));
  c2->SaveAs(Form("fhist_dqdxcal.png"));
  c2->Clear();

  fhist_dqdxuncal->Write();
  fhist_dqdxuncal->Draw("COLZ");
  c2->SaveAs(Form("fhist_dqdxuncal.pdf"));
  c2->SaveAs(Form("fhist_dqdxuncal.png"));
  c2->Clear();

  dqdx_cal->Write();
  dqdx_cal->Draw();
  c2->SaveAs(Form("dqdx_cal.pdf"));
  c2->SaveAs(Form("dqdx_cal.png"));
  c2->Clear();

  dqdx_uncal->Write();
  dqdx_uncal->Draw();
  c2->SaveAs(Form("dqdx_uncal.pdf"));
  c2->SaveAs(Form("dqdx_uncal.png"));
  c2->Clear();

  gr3->Write("measured_dedxrr");
  gr3->SetTitle(" ; Residual range [cm]; dE/dx [MeV/cm]");

  gr3->Draw("ACP");

  gPad->SaveAs(Form("measured_dEdx_rr.pdf"));
  gPad->SaveAs(Form("measured_dEdx_rr.png"));
  c2->Clear();

  theoretical->Write("theoretical");

  theoretical->SetTitle(" ; Residual range [cm]; dE/dx [MeV/cm]");

  theoretical->Draw("ACP");

  gPad->SaveAs(Form("theory_dEdx_rr.pdf"));
  gPad->SaveAs(Form("theory_dEdx_rr.png"));
  c2->Clear();

  std::cout << "************************* Calibration.C has ended ***************************" << std::endl;

  gStyle->SetPalette(1, 0);
  gStyle->SetNumberContours(64);
  fhist_dqdxcal->Draw("colz");
  c2->SaveAs(Form("fhist_dqdxcal_last.png"));
  c2->Clear();
  c2->Close();

  /////DRAW HISTOGRAM
  //
  TCanvas *c3 = new TCanvas();
  h_reco_dQdx_RR->Draw("COLZ");
  h_reco_dQdx_RR->Write(" h_reco_dQdx_RR");

  c3->SaveAs("plots_reco/plots_pdf/reco_dqdx_rr.pdf");
  c3->SaveAs("plots_reco/plots_pdf/reco_dqdx_rr.png");
  c3->Write("reco_dQdx_RR");
  c3->Clear();

  h_reco_dEdx_RR_uncal->Write();
  h_reco_dEdx_RR_uncal->Draw("COLZ");

  c3->SaveAs("plots_reco/plots_pdf/reco_dEdx_rr_uncal.pdf");
  c3->SaveAs("plots_reco/plots_pdf/reco_dEdx_rr_uncal.png");
  c3->Clear();

  h_reco_dEdx_RR_cal->Write();
  h_reco_dEdx_RR_cal->Draw("COLZ");

  c3->SaveAs("plots_reco/plots_pdf/reco_dEdx_rr_cal.pdf");
  c3->SaveAs("plots_reco/plots_pdf/reco_dEdx_rr_cal.png");
  c3->Clear();

  fhist_dedx->Write();
  fhist_dedx->Draw("COLZ");
  c3->SaveAs("plots_reco/plots_pdf/fhist_dedx.pdf");
  c3->SaveAs("plots_reco/plots_pdf/fhist_dedx.png");
  c3->Clear();

  c3->Close();

  gStyle->SetPalette(1, 0);
  gStyle->SetNumberContours(64);
}
