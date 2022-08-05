/************************************************************************
 *
 * A macro to run studies of the Pandora reconstructed objects in the 
 * sample.
 * 
 * Note (July 2022):
 *    This has not yet been updated to reflect the additional parameters 
 *    included in the production.
 *
 * Example file list located here:
 *   /home/jones/work/cosmics/LArSoft-v08_50_00/work/files/v09_41_00_02_files.list
 *
 *************************************************************************/

#include "Classes/EventProcessor.h"
#include "Classes/ConfigReader.h"

using namespace calib;
using namespace cppsecrets;

// Allowed branches to read from the tree
std::vector<TString> allowed = {
   "run",
   "event",
   "genie_primaries_pdg",
   "inTPCActive",
   "pdg",
   "TrackId",
   "Mother",
   "Eng",
   "pathlen",
   "StartPointx",
   "StartPointy",
   "StartPointz",
   "EndPointx",
   "EndPointy",
   "EndPointz",
   "StartPointx_tpcAV",
   "StartPointy_tpcAV",
   "StartPointz_tpcAV",
   "EndPointx_tpcAV",
   "EndPointy_tpcAV",
   "EndPointz_tpcAV",
   "StartE_tpcAV",
   "NumberDaughters",
   "hit_plane",
   "hit_trkid",
   "hit_nelec",
   "process_primary",
   "trkpdgtruth_pandoraTrack",
   "trkg4id_pandoraTrack",
   "trkidtruth_pandoraTrack",
   "ntracks_pandoraTrack",
   "trkId_pandoraTrack",
   "ntrkhits_pandoraTrack",
   "trkdqdx_pandoraTrack",
   "trkdedx_pandoraTrack",
   "trkxyz_pandoraTrack",
   "trkstartx_pandoraTrack",
   "trkstarty_pandoraTrack",
   "trkstartz_pandoraTrack",
   "trkendx_pandoraTrack",
   "trkendy_pandoraTrack",
   "trkendz_pandoraTrack",
   "trklen_pandoraTrack",
   "trkke_pandoraTrack",
   "nPFParticles",
   "pfp_selfID",
   "pfp_isPrimary",
   "pfp_isTrack",
   "pfp_trackID",
   "pfp_numClusters",
   "pfp_clusterIDs",
   "pfp_numDaughters",
   "pfp_daughterIDs"
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
     
int pandoraStudies(const char *config){

  // First, setup timing information so we can monitor the run
  time_t rawtime;
  std::cout << "-----------------------------------------------------------" << std::endl;
  GetTime(rawtime);
  std::cout << "-----------------------------------------------------------" << std::endl;

  //------------------------------------------------------------------------------------------
  //                                    Configure
  //------------------------------------------------------------------------------------------
  // Create object of the classConfigReader
  // Parse the configuration file
  // Dump map on the console after parsing it
  ConfigReader* p = ConfigReader::getInstance();
  p->parseFile(config);
  std::cout << " Variables from configuration file: " << std::endl;
  p->dumpFileValues();
  std::cout << "-----------------------------------------------------------" << std::endl;

  // Get configuration variables
  int n = -1;
  int yCut = 1;
  int thru = 0;
  std::string input_list = "";
  std::string location="";
  std::string tag="";
  std::vector<double> minx_fid, miny_fid, minz_fid;
  std::vector<double> maxx_fid, maxy_fid, maxz_fid;
  std::vector<double> minx_av, miny_av, minz_av;
  std::vector<double> maxx_av, maxy_av, maxz_av;

  p->getValue("InputList", input_list);
  p->getValue("Location",  location);
  p->getValue("Tag",       tag);
  p->getValue("NFiles",    n);
  p->getValue("YCut",      yCut);
  p->getValue("Thru",      thru);
  p->getValue("MinXFid",   minx_fid);
  p->getValue("MinYFid",   miny_fid);
  p->getValue("MinZFid",   minz_fid);
  p->getValue("MaxXFid",   maxx_fid);
  p->getValue("MaxYFid",   maxy_fid);
  p->getValue("MaxZFid",   maxz_fid);
  p->getValue("MinXAV",    minx_av);
  p->getValue("MinYAV",    miny_av);
  p->getValue("MinZAV",    minz_av);
  p->getValue("MaxXAV",    maxx_av);
  p->getValue("MaxYAV",    maxy_av);
  p->getValue("MaxZAV",    maxz_av);

  // Get the active and fiducial geometry objects
  Geometry fiducial(minx_fid,miny_fid,minz_fid,maxx_fid,maxy_fid,maxz_fid,true);
  Geometry active(minx_av,miny_av,minz_av,maxx_av,maxy_av,maxz_av,false);
  PlaneList extPlanes = active.GetExternalPlaneList();
  PlaneList allPlanes = active.GetPlaneList();
  PlaneList intPlanes = active.GetInternalPlaneList(allPlanes,extPlanes);
  PlaneList fidExtPlanes = fiducial.GetExternalPlaneList();
  PlaneList fidAllPlanes = fiducial.GetPlaneList();

  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Total number of planes in the active volume of the DUNE SP module: " << allPlanes.size() << std::endl;
  std::cout << " Consisting of " << extPlanes.size() << " external planes and " << intPlanes.size() << " internal planes" << std::endl; 
  std::cout << "-----------------------------------------------------------" << std::endl;
 
  // Sort out the file tag
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
  n = evtProc.GetFiles();
  std::cout << " Number of files: " << n << std::endl;
  
  // Start of analysis (loop over chain and events
  std::cout << " Running analysis..." << std::endl;

  // Then setup the histograms, counters and any other variables to add to
  // Setup histograms
  // 1D
  TH1D *h_true_length        = new TH1D("h_true_length","",100,0,30);
  TH1D *h_reco_length        = new TH1D("h_reco_length","",100,0,30);
  TH1D *h_reco_true_length   = new TH1D("h_reco_true_length","",100,0,30);
  TH1D *h_true_energy        = new TH1D("h_true_energy","",100,0.3,5000);
  TH1D *h_reco_energy        = new TH1D("h_reco_energy","",100,0.3,5000);
  TH1D *h_reco_true_energy   = new TH1D("h_reco_true_energy","",100,0.3,5000);
  TH1D *h_daughter_diff      = new TH1D("h_daughter_diff","",200,-100,100);

  // Log scale for energies
  SetLogX(h_true_energy);
  SetLogX(h_reco_energy);
  SetLogX(h_reco_true_energy);

  // Setup counters
  // Numerators
  unsigned int nPrimary_reco_assoc       = 0;
  unsigned int nLength_reco_assoc        = 0;
  unsigned int nPrimaryLength_reco_assoc = 0;
 
  // Denominators
  unsigned int nPrimary_true             = 0; // Efficiency
  unsigned int nPrimary_TG_true          = 0; // Efficiency, through-going
  unsigned int nPrimaryLength_reco       = 0; // Purity
  unsigned int nLength_reco              = 0; // Purity
  unsigned int nPrimary_reco             = 0; // Purity

  // Now loop over the events
  unsigned int nEvts = tree->GetEntries();
  unsigned int iIt = 1;

  std::cout << " |";
  for(unsigned int iEvt = 0; iEvt < nEvts; ++iEvt){
    tree->GetEntry(iEvt);
    if(!evtProc.SelectEvent(evt)) continue;
    int nGeant         = evt->geant_list_size;
    unsigned int nTrks = evt->ntracks_pandoraTrack;
    unsigned int nPfps = evt->nPFParticles;
    
    // Print the processing rate
    double evtFrac  = iEvt/static_cast<double>(nEvts);

    if((std::abs(0.1*iIt)-evtFrac) < std::numeric_limits<double>::epsilon()){
      std::cout << " --- " << evtFrac*100 << " %";
      std::cout.flush();
      iIt++;
    }
    //
    // Truth-level studies
    //
    // Loop over geant tracks before looping through reconstructed tracks and determining efficiencies
    for(int iG4 = 0; iG4 < nGeant; ++iG4){

      // Check the particle enters the TPC volume
      if(!evt->inTPCActive[iG4]) continue;

      TVector3 vtx(evt->StartPointx[iG4],evt->StartPointy[iG4],evt->StartPointz[iG4]);
      TVector3 end(evt->EndPointx[iG4],evt->EndPointy[iG4],evt->EndPointz[iG4]);
      
      TVector3 vtxAV(evt->StartPointx_tpcAV[iG4],evt->StartPointy_tpcAV[iG4],evt->StartPointz_tpcAV[iG4]);
      TVector3 endAV(evt->EndPointx_tpcAV[iG4],evt->EndPointy_tpcAV[iG4],evt->EndPointz_tpcAV[iG4]);

      // If these don't match, the TPC start and end point and general start and end point are not same, 
      bool throughGoing = IsTrueThroughGoing(vtx,end,vtxAV,endAV);
      
      int pdg        = evt->pdg[iG4];
      int id         = evt->TrackId[iG4];
      float length   = evt->pathlen[iG4]/100.; // [m]
      float energy   = evt->StartE_tpcAV[iG4];
      
      // Find our true primary muons and make sure they are through-going if desired
      if(abs(pdg) != 13) continue;
      if(evt->Mother[iG4] != 0) continue;
      nPrimary_true++;

      if(thru != throughGoing) continue; 
      nPrimary_TG_true++;

      // Fill true histograms
      h_true_length->Fill(length);
      h_true_energy->Fill(energy);

      // Now loop over reconstructed quantities and find associations with truth
      for(unsigned int iPfp = 0; iPfp < nPfps; ++iPfp){
        int PFPTrackID = evt->pfp_trackID[iPfp];
        int PFPID      = evt->pfp_selfID[iPfp];

        // Is primary 
        if(evt->pfp_isPrimary[iPfp])
          nPrimary_reco++;

        // Check that the current track ID is associated with the current G4 particle
        bool trueAssoc       = false;
        bool trueAssocLength = false;
        float reco_length = -99999.;
        float reco_energy = -99999.;
        for(unsigned int iTrk = 0; iTrk < nTrks; ++iTrk){

          // Get the best plane
          int bestPlane = GetBestPlane(std::vector<int>(std::begin(evt->ntrkhits_pandoraTrack[iTrk]),std::end(evt->ntrkhits_pandoraTrack[iTrk])));

          // Make sure we're looking at sensible reconstructed tracks
          if(evt->trkId_pandoraTrack[iTrk] != PFPTrackID) continue;

          // Is the track/pfp associated with the true G4 particle regardless of its length
          if(evt->trkg4id_pandoraTrack[iTrk] == id)
            trueAssoc = true;

          // Is long
          if(!evtProc.SelectTrack(evt,iTrk)) continue;
          nLength_reco++;

          // Is primary and long
          if(evt->pfp_isPrimary[iPfp])
            nPrimaryLength_reco++;

          // Now check if this track corresponds to the G4 track
          if(evt->trkg4id_pandoraTrack[iTrk] == id) {
            trueAssocLength = true;
            
            // Fill reco quantities
            reco_length = evt->trklen_pandoraTrack[iTrk]/100.; //[m]
            reco_energy = evt->trkke_pandoraTrack[iTrk][bestPlane]/1000.;

            break;
          }
        } // iTrk
        if(!trueAssoc) continue;

        // Counter for how many associations we have found
        unsigned int nAssoc = 0;
        if(evt->pfp_isPrimary[iPfp]){
          nPrimary_reco_assoc++;
          nAssoc++;
        }

        if(trueAssocLength){
          nLength_reco_assoc++;
          nAssoc++;
        }

        if(evt->pfp_isPrimary[iPfp] && trueAssocLength){
          nPrimaryLength_reco_assoc++;
          nAssoc++;

          // Fill histograms
          h_reco_length->Fill(reco_length);
          h_reco_energy->Fill(reco_energy);
          h_reco_true_energy->Fill(energy);
          h_reco_true_length->Fill(length);
          h_daughter_diff->Fill(evt->pfp_numDaughters[iPfp] - evt->NumberDaughters[iG4]);

        } // Reco primary and length

        if(nAssoc > 0)
          break;

      }// iPfp
    }// iG4
  }// Event loop
  std::cout << " --- 100 % --- |" << std::endl;

  // PLOTS
  // Get the efficiency histograms
  TH1D *h_eff_energy = static_cast<TH1D*>(h_reco_true_energy->Clone("h_eff_energy"));
  h_eff_energy->Divide(h_true_energy);

  TH1D *h_eff_length = static_cast<TH1D*>(h_reco_true_length->Clone("h_eff_length"));
  h_eff_length->Divide(h_true_length);

  // Overlay true, reco and efficiency with efficiency axis on RHS
  // Energy
  TCanvas *c0 = new TCanvas("c0","",900,800);
  SetCanvasStyle(c0, 0.12,0.12,0.08,0.12,0,0,0);
  c0->SetLogx();
  //c0->SetLogy();

  SetHistogramStyle1D(h_true_energy,"Muon energy [GeV]", "Rate");
  SetHistogramStyle1D(h_reco_true_energy,"Muon energy [GeV]", "Rate");

  double leftMax = 1.1*std::max(h_true_energy->GetMaximum(),h_reco_true_energy->GetMaximum());

  h_true_energy->GetYaxis()->SetRangeUser(0,leftMax);
  h_true_energy->Draw("hist");
  h_reco_true_energy->Draw("same hist");
  h_true_energy->SetMarkerColor(kWhite); // For the legend
  h_true_energy->SetLineWidth(3);
  h_reco_true_energy->SetLineWidth(3);
  h_true_energy->SetLineColor(kTeal-5);
  h_reco_true_energy->SetLineColor(kViolet-5);
  h_reco_true_energy->SetLineStyle(7);

  // Now setup the RHS
  double rightMax = 1.1*h_eff_energy->GetMaximum();
  double scale    = leftMax/rightMax;
  std::cout << " Efficiency max: " << h_eff_energy->GetMaximum() << std::endl;
  
  h_eff_energy->Scale(scale);
  h_eff_energy->Draw("hist same");
  h_eff_energy->SetLineWidth(3);
  h_eff_energy->SetLineStyle(2);
  h_eff_energy->SetLineColor(kOrange+5);
  
  TGaxis *rhsAxis = new TGaxis(h_eff_energy->GetXaxis()->GetXmax(),0,h_eff_energy->GetXaxis()->GetXmax(),leftMax,0,rightMax,510,"+L");
  rhsAxis->SetTitle("Efficiency");
  rhsAxis->SetLineColor(kOrange+5);
  rhsAxis->SetLabelColor(kOrange+5);
  rhsAxis->SetLabelFont(132);
  rhsAxis->SetLabelSize(0.045);
  rhsAxis->SetLabelOffset(0.01);
  rhsAxis->SetTitleOffset(1.1);
  rhsAxis->SetTitleFont(132);
  rhsAxis->SetTitleSize(0.055);
  rhsAxis->SetTitleColor(kOrange+5);
  rhsAxis->Draw();

  TLegend *l = new TLegend(0.18,0.92,0.94,0.995);
  l->SetNColumns(2);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->SetTextFont(132);

  l->AddEntry(h_true_energy,"True energy","L");
  l->AddEntry(h_reco_true_energy,"Reconstructed muon, true energy","L");
  l->AddEntry(h_true_energy," ","P");
  l->AddEntry(h_eff_energy,"Reconstruction efficiency","L");

  l->Draw();
  
  c0->SaveAs((location+"/energy_eff"+tag+".png").c_str());
  c0->SaveAs((location+"/energy_eff"+tag+".root").c_str());
  c0->Clear();
  l->Clear();
  
  // Length
  TCanvas *c1 = new TCanvas("c1","",900,800);
  SetCanvasStyle(c1, 0.12,0.12,0.08,0.12,0,0,0);

  SetHistogramStyle1D(h_true_length,"Muon length [m]", "Rate");
  SetHistogramStyle1D(h_reco_length,"Muon length [m]", "Rate");
  SetHistogramStyle1D(h_reco_true_length,"Muon length [m]", "Rate");

  leftMax = 1.1*std::max(h_true_length->GetMaximum(), std::max(h_reco_true_length->GetMaximum(),h_reco_length->GetMaximum()));

  h_true_length->GetYaxis()->SetRangeUser(0,leftMax);
  h_true_length->Draw("hist");
  h_reco_length->Draw("hist same");
  h_reco_true_length->Draw("same hist");
  h_true_length->SetLineWidth(3);
  h_reco_length->SetLineWidth(3);
  h_reco_true_length->SetLineWidth(3);
  h_true_length->SetLineColor(kTeal-5);
  h_reco_length->SetLineColor(kAzure-5);
  h_reco_true_length->SetLineColor(kViolet-5);
  h_reco_true_length->SetLineStyle(7);

  // Now setup the RHS
  rightMax = 1.1*h_eff_length->GetMaximum();
  scale    = leftMax/rightMax;
  std::cout << " Efficiency max: " << h_eff_length->GetMaximum() << std::endl;
  
  h_eff_length->Scale(scale);
  h_eff_length->Draw("hist same");
  h_eff_length->SetLineWidth(3);
  h_eff_length->SetLineStyle(2);
  h_eff_length->SetLineColor(kOrange+4);
  
  rhsAxis->DrawAxis(h_eff_length->GetXaxis()->GetXmax(),0,h_eff_length->GetXaxis()->GetXmax(),leftMax,0,rightMax,510,"+L");

  l->AddEntry(h_true_length,"True length","L");
  l->AddEntry(h_reco_true_length,"Reconstructed muon, true length","L");
  l->AddEntry(h_reco_length,"Reco length","L");
  l->AddEntry(h_eff_length,"Reconstruction efficiency","L");

  l->Draw();
  
  c1->SaveAs((location+"/length_eff"+tag+".png").c_str());
  c1->SaveAs((location+"/length_eff"+tag+".root").c_str());
  c1->Clear();
  l->Clear();
  
  // Single-axis plots
  // Number of daughters reco-truth
  SetHistogramStyle1D(h_daughter_diff,"N_{Daughters,Reco}-N_{Daughters,True}", "Rate");
  h_daughter_diff->Draw("hist");
  h_daughter_diff->SetLineWidth(3);
  h_daughter_diff->SetLineColor(kTeal-5);

  // Fit a Gaussian to the distribution
  // Find the min and max from the centre of the histogram, fit 1/8 of the total number of bins either side
  int maxbin    = h_daughter_diff->GetMaximumBin();
  int nbins     = h_daughter_diff->GetNbinsX();
  double maxloc = h_daughter_diff->GetBinCenter(maxbin);
  double fitMin = h_daughter_diff->GetBinCenter(maxbin-(nbins/8.)); 
  double fitMax = h_daughter_diff->GetBinCenter(maxbin+(nbins/8.)); 

  TF1 *fit = fit = new TF1("fit","gaus",fitMin,fitMax);
  auto result = h_daughter_diff->Fit(fit, "LEQSMR", "");
  fit->Draw("hist same");
  fit->SetLineWidth(3);
  fit->SetLineStyle(3);
  fit->SetLineColor(kViolet-5);

  l->AddEntry(h_daughter_diff,"#Delta N_{Daughters}","L");
  l->AddEntry(fit,"Gaussian fit","L");

  l->SetTextSize(0.038);
  l->Draw();

  c1->SaveAs((location+"/nDaughter_diff"+tag+".png").c_str());
  c1->SaveAs((location+"/nDaughter_diff"+tag+".root").c_str());
  c1->Clear();

  // COUNTERS
  ofstream txtFile;
  txtFile.open(location+"/statistics"+tag+".txt");
  txtFile << "----------------------------------------------------------------------------------------" << std::endl;
  txtFile << " Number of primary particles: " << std::endl;
  txtFile << " True:          " << nPrimary_true << std::endl;
  if(thru)
    txtFile << " True through-going: " << nPrimary_TG_true << std::endl;
  txtFile << " Reconstructed:      " << nPrimary_reco << std::endl;
  txtFile << "----------------------------------------------------------------------------------------" << std::endl;
  txtFile << " Efficiency of reconstructing all long muons:          " << (nLength_reco_assoc/static_cast<double>(nPrimary_true))*100. << " %" << std::endl;
  txtFile << " Efficiency of reconstructing all primary muons:       " << (nPrimary_reco_assoc/static_cast<double>(nPrimary_true))*100. << " %" << std::endl;
  txtFile << " Efficiency of reconstructing all primary, long muons: " << (nPrimaryLength_reco_assoc/static_cast<double>(nPrimary_true))*100. << " %" << std::endl;
  txtFile << std::endl;
  txtFile << " Efficiency of reconstructing through-going long muons:          " << (nLength_reco_assoc/static_cast<double>(nPrimary_TG_true))*100. << " %" << std::endl;
  txtFile << " Efficiency of reconstructing through-going primary muons:       " << (nPrimary_reco_assoc/static_cast<double>(nPrimary_TG_true))*100. << " %" << std::endl;
  txtFile << " Efficiency of reconstructing through-going primary, long muons: " << (nPrimaryLength_reco_assoc/static_cast<double>(nPrimary_TG_true))*100. << " %" << std::endl;
  txtFile << std::endl;
  txtFile << " Purity of reconstructing long muons:          " << (nLength_reco_assoc/static_cast<double>(nLength_reco))*100. << " %" << std::endl;
  txtFile << " Purity of reconstructing primary muons:       " << (nPrimary_reco_assoc/static_cast<double>(nPrimary_reco))*100. << " %" << std::endl;
  txtFile << " Purity of reconstructing primary, long muons: " << (nPrimaryLength_reco_assoc/static_cast<double>(nPrimaryLength_reco))*100. << " %" << std::endl;
  txtFile << "----------------------------------------------------------------------------------------" << std::endl;
  txtFile << " Daughter difference Gaussian fit results" << std::endl;
  txtFile << " Mean        = " << fit->GetParameter(1) << std::endl;
  txtFile << " StdDev      = " << fit->GetParameter(2) << std::endl;
  txtFile << " Chi^2/ndof  = " << fit->GetChisquare()/static_cast<double>(fit->GetNDF()) << std::endl;
  txtFile << "----------------------------------------------------------------------------------------" << std::endl;

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
