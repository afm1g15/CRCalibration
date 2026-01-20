/************************************************************************
 * 
 * A macro to make an event selection on true and reco stopping muons for 
 * dE/dx calibration studies
 *
 *
 * Input is a list of ana files.
 * Example file list located here:
 *   /exp/dune/app/users/amoor/duneCalibration/anafiles.list
 *
 *
 *************************************************************************/

#include "EventProcessor.h"
#include "ConfigReader.h"
#include "TTree.h"

using namespace calib;
using namespace cppsecrets;

// Allowed branches to read from the tree (all branches in afmanatree_core.h)
std::vector<TString> allowed = {
   "runID",
   "eventID",
   "nPFParticles",
   "nMCParticles",
   "mcParticlePdgCode",
   "mcParticleTrackID",
   "mcParticleParentTrackID",
   "mcParticleTrueEnergy",
   "mcParticleStartPositionX",
   "mcParticleStartPositionY",
   "mcParticleStartPositionZ",
   "mcParticleEndPositionX",
   "mcParticleEndPositionY",
   "mcParticleEndPositionZ",
   "pfpTrueParticleMatchedID",
   "pfpIsShower",
   "pfpID",
   "pfpShowerEnergy",
   "pfpShowerDirectionX",
   "pfpShowerDirectionY",
   "pfpShowerDirectionZ",
   "pfpShowerStartX",
   "pfpShowerStartY",
   "pfpShowerStartZ",
   "pfpShowerLength"
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
     
int stoppingPi0Study(const char *config){

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
  int thru = 0; // Do we want to select only through-going muons? Default: No (0)
  int stop = 0; // Do we want to select only stopping muons? Default: No (0)
  std::string input_list = "";
  std::string location="";
  std::string tag="";
  std::vector<double> minx_fid, miny_fid, minz_fid;
  std::vector<double> maxx_fid, maxy_fid, maxz_fid;
  std::vector<double> minx_av, miny_av, minz_av;
  std::vector<double> maxx_av, maxy_av, maxz_av;

  // Access corresponding parameter in the configuration file
  p->getValue("InputList", input_list);
  p->getValue("Location",  location);
  p->getValue("Tag",       tag);
  p->getValue("NFiles",    n);
  p->getValue("Thru",      thru);
  p->getValue("Stopping",  stop);
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

  // Sanity check the geometry definitions
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Total number of planes in the active volume of the DUNE SP module: " << allPlanes.size() << std::endl;
  std::cout << " Consisting of " << extPlanes.size() << " external planes and " << intPlanes.size() << " internal planes" << std::endl; 
  std::cout << "-----------------------------------------------------------" << std::endl;
 
  // Sort out the file tag by adding an underscore
  if(tag != "")
    tag = "_"+tag;

  //--------------------------------------------------------------------------------- ---------
  //                                    Initialise
  //--------------------------------------------------------------------------------- ---------

  // Setup TTree from input file list
  std::cout << " Reading files and filling tree..." << std::endl;

  EventProcessor evtProc(allowed, input_list, n);
  evtProc.Initialize();

  // Now setup the tree and event objects to work with
  TChain *tree = evtProc.GetTree();
  afmanatree *evt = evtProc.GetEvents();

  
  // Start of analysis (loop over chain and events
  std::cout << " Running analysis..." << std::endl;

  // Then setup the histograms, counters and any other variables to add to
  // Setup histograms if wanted
  TH1D *h_truth_inv_pi0_mass   = new TH1D("h_truth_inv_pi0_mass","",50,0,0.5);
  TH1D *h_reco_inv_pi0_mass   = new TH1D("h_reco_inv_pi0_mass","",50,0,0.5);
  TH1D *h_truth_shower_energy   = new TH1D("h_truth_shower_energy","",50,0,0.5);
  TH1D *h_reco_shower_energy   = new TH1D("h_reco_shower_energy","",50,0,0.5);
  TH1D *h_truth_cos_angle   = new TH1D("h_truth_cos_angle","",20,-1,1);
  TH1D *h_reco_cos_angle   = new TH1D("h_reco_cos_angle","",20,-1,1);
  TH1D *h_truth_dirX   = new TH1D("h_truth_dirX","",20,-1,1);
  TH1D *h_reco_dirX   = new TH1D("h_reco_dirX","",20,-1,1);
  TH1D *h_truth_dirY   = new TH1D("h_truth_dirY","",20,-1,1);
  TH1D *h_reco_dirY   = new TH1D("h_reco_dirY","",20,-1,1);
  TH1D *h_truth_dirZ   = new TH1D("h_truth_dirZ","",20,-1,1);
  TH1D *h_reco_dirZ   = new TH1D("h_reco_dirZ","",20,-1,1);
  
  // Setup counters
  unsigned int totalParticlesTrue = 0;
  unsigned int trueSignalPi0s = 0;
  unsigned int totalParticlesReco = 0;
  unsigned int recoSelectedPhotons = 0;
  unsigned int recoSelectedSignalPi0s = 0;
  
  // Now loop over the events
  unsigned int nEvts = tree->GetEntries();
  unsigned int iIt = 1;
  unsigned int trackRepeats = 0;
  unsigned int eventNum = 0;

  std::cout << " |";
  for(unsigned int iEvt = 0; iEvt < nEvts; ++iEvt){
    tree->GetEntry(iEvt);
    if(!evtProc.SelectEvent(evt)) continue;
    
    // Get the total number of true and reconstructed tracks to loop over
    int nPfps = evt->nPFParticles;            //reco
    int nMcs = evt->nMCParticles;                //true
    
    std::cout << "There are " << nPfps << " reconstructed particles and " << nMcs << " true particles in this event." << std::endl;


    // Print the processing rate
    double evtFrac  = iEvt/static_cast<double>(nEvts);

    // Prints out how much has been completed so far
    if(std::abs(0.1*iIt-evtFrac) < std::numeric_limits<double>::epsilon()){
      std::cout << " --- " << evtFrac*100 << " %";
      std::cout.flush();
      iIt++;
    }
  
    ///////////////////////////////////
    //          TRUTH                //
    ///////////////////////////////////          
    std::vector<int> trueMCPassId; //vector of track IDs that pass true signal cuts
    std::vector<int> signalpi0IDs;
    // Now loop over the true tracks
    for(int iMc = 0; iMc < nMcs; ++iMc){
      // Count tracks
      totalParticlesTrue++;

      //Look for true pdg
      int trupdg = evt->mcParticlePdgCode[iMc];
      if ( !(abs(trupdg) == 111) )
        continue;

      int signalpi0ID = evt->mcParticleTrackID[iMc];

      trueSignalPi0s++;
      signalpi0IDs.push_back(signalpi0ID);

    } // iMc, truth loop

    std::cout << "In this event, " << signalpi0IDs.size() << " pi0s were found in MC." << std::endl;

    std::vector<int> signalphotonIDs;
    std::vector<int> signalphotonNs;
    for(int iMc = 0; iMc < nMcs; ++iMc){

      int trupdg = evt->mcParticlePdgCode[iMc];
      int trunextpdg = evt->mcParticlePdgCode[iMc+1];
      if ( !(abs(trupdg) == 22 && abs(trunextpdg) == 22) )
        continue;

      
      int motherID = evt->mcParticleParentTrackID[iMc];
      int mothernextID = evt->mcParticleParentTrackID[iMc+1];
      if(!(CheckTrueIDAssoc(motherID,signalpi0IDs) && (motherID == mothernextID)))
        continue;


      int signalp1ID = evt->mcParticleTrackID[iMc];
      int signalp2ID = evt->mcParticleTrackID[iMc+1];

      signalphotonIDs.push_back(signalp1ID);
      signalphotonIDs.push_back(signalp2ID);
      signalphotonNs.push_back(iMc);
      signalphotonNs.push_back(iMc+1);
      std::cout << "A photon pair has found as a child of pi0 " << motherID << " with the track IDs " << signalp1ID << " and " << signalp2ID << std::endl;

      float trueE1 = evt->mcParticleTrueEnergy[iMc];
      float trueE2 = evt->mcParticleTrueEnergy[iMc+1];
      std::cout << "They have true energies of " << trueE1 << " and " << trueE2 << std::endl;
      h_truth_shower_energy->Fill(trueE1);
      h_truth_shower_energy->Fill(trueE2);

      TVector3 start1(evt->mcParticleStartPositionX[iMc],evt->mcParticleStartPositionY[iMc],evt->mcParticleStartPositionZ[iMc]);
      TVector3 end1(evt->mcParticleEndPositionX[iMc],evt->mcParticleEndPositionY[iMc],evt->mcParticleEndPositionZ[iMc]);

      TVector3 start2(evt->mcParticleStartPositionX[iMc+1],evt->mcParticleStartPositionY[iMc+1],evt->mcParticleStartPositionZ[iMc+1]);
      TVector3 end2(evt->mcParticleEndPositionX[iMc+1],evt->mcParticleEndPositionY[iMc+1],evt->mcParticleEndPositionZ[iMc+1]);

      float length1 = (end1 - start1).Mag();
      float length2 = (end2 - start2).Mag();
      TVector3 dir1 = (end1 - start1).Unit();
      TVector3 dir2 = (end2 - start2).Unit();
      float crosslength = (end1 - end2).Mag();
      std::cout << "They have true lengths of " << length1 << " and " << length2 << ", crosslength " << crosslength << std::endl;
      h_truth_dirX->Fill(dir1.X());
      h_truth_dirX->Fill(dir2.X());
      h_truth_dirY->Fill(dir1.Y());
      h_truth_dirY->Fill(dir2.Y());
      h_truth_dirZ->Fill(dir1.Z());
      h_truth_dirZ->Fill(dir2.Z());


      //float cos_angle = ((length1*length1)+(length2*length2)-(crosslength*crosslength))/(2.0*length1*length2);
      float opening_angle = dir1.Angle(dir2);
      float cos_angle = TMath::Cos(opening_angle); 

      std::cout << "cos(A) = " << cos_angle << std::endl;
      h_truth_cos_angle->Fill(cos_angle);     

      float pi0_mass = std::sqrt(2*trueE1*trueE2*(1 - cos_angle));
      std::cout << "truth invarient pi0 mass = " << pi0_mass << std::endl;
      h_truth_inv_pi0_mass->Fill(pi0_mass);
      std::cout << "-------------------------------" << std::endl;

    }

    std::cout << "-------------------------------" << std::endl;

    ///////////////////////////////////
    //            RECO               //
    ///////////////////////////////////
    std::vector<int> recoPfpPassId; //vector of true track IDs that pass reco cuts
    std::vector<int> recoSignalPhotonNs;
    std::vector<int> recoSignalPhotonMatchIDs;
    std::vector<std::pair<float, float>> matchedPairs;
    for(int iPfp = 0; iPfp < nPfps; ++iPfp){

      // Count tracks
      totalParticlesReco++;

      int pfpID = evt->pfpID[iPfp];
      int pfpTrueMatchID = evt->pfpTrueParticleMatchedID[iPfp];
      int pfpIsShower = evt->pfpIsShower[iPfp];

      if (!pfpIsShower)
        continue;


     if(!CheckTrueIDAssoc(pfpTrueMatchID,signalphotonIDs))
         continue;

      std::cout << "TrackID of true MC matched to this Pfp: " << pfpTrueMatchID << std::endl;

      auto photonMatches = std::make_pair(iPfp, pfpTrueMatchID);

      matchedPairs.push_back(photonMatches);

     recoSelectedPhotons++;

    } // iTrk, reco loop

   std::vector<int> pairNums;
   for (long unsigned int i = 0; i < matchedPairs.size(); i = i+1) {
     for (long unsigned int j = 0; j < matchedPairs.size(); j = j+1) {
       auto photon1 = matchedPairs[i].second;
       auto photon2 = matchedPairs[j].second;
         bool is_unique(true);
         if (pairNums.size() != 0) {
           for (long unsigned int k = 0; k < pairNums.size(); k = k+1)  {
             if ((photon1 == matchedPairs[k].second) || (photon2 == matchedPairs[k].second)) { 
                is_unique = false;
             }
           }
         }
       if ((photon1 != photon2) && (is_unique)) {
         pairNums.push_back(matchedPairs[i].first);
         pairNums.push_back(matchedPairs[j].first);
       }
    }
  }


    std::vector<float> reco_pi0_masses;
    for (long unsigned int i = 0; i < pairNums.size(); i = i+2) {

       int iPfp1 = pairNums[i];
       int iPfp2 = pairNums[i+1];
     

       float energy1 = evt->pfpShowerEnergy[iPfp1];
       float energy2 = evt->pfpShowerEnergy[iPfp2];

       h_reco_shower_energy->Fill(energy1);
       h_reco_shower_energy->Fill(energy2);

       TVector3 dir1(evt->pfpShowerDirectionX[iPfp1],evt->pfpShowerDirectionY[iPfp1],evt->pfpShowerDirectionZ[iPfp1]);
       TVector3 start1(evt->pfpShowerStartX[iPfp1],evt->pfpShowerStartY[iPfp1],evt->pfpShowerStartZ[iPfp1]);
       float length1 = evt->pfpShowerLength[iPfp1];

       TVector3 dir2(evt->pfpShowerDirectionX[iPfp2],evt->pfpShowerDirectionY[iPfp2],evt->pfpShowerDirectionZ[iPfp2]);
       TVector3 start2(evt->pfpShowerStartX[iPfp2],evt->pfpShowerStartY[iPfp2],evt->pfpShowerStartZ[iPfp2]);
       float length2 = evt->pfpShowerLength[iPfp2];
       h_reco_dirX->Fill(dir1.X());
       h_reco_dirX->Fill(dir2.X());
       h_reco_dirY->Fill(dir1.Y());
       h_reco_dirY->Fill(dir2.Y());
       h_reco_dirZ->Fill(dir1.Z());
       h_reco_dirZ->Fill(dir2.Z());

       TVector3 end1 = start1+(dir1*length1);
       TVector3 end2 = start2+(dir2*length2);

       std::cout << "-------------------------------------" << std::endl;
       std::cout << "Inputs: " << std::endl;
       std::cout << "- Energies = " << energy1 << " and " << energy2 << std::endl;
       std::cout << "- directions = (" << dir1.X()<< "," << dir1.Y() << "," << dir1.Z() << ") and " << dir2.X() << "," << dir2.Y() << "," << dir2.Z() << std::endl;

      float crosslength = (end1 - end2).Mag();
      std::cout << "They have true lengths of " << length1 << " and " << length2 << ", crosslength " << crosslength << std::endl;

      //float cos_angle = ((length1*length1)+(length2*length2)-(crosslength*crosslength))/(2.0*length1*length2);


       float opening_angle = dir1.Angle(dir2);
       float cos_angle = TMath::Cos(opening_angle);
       std::cout << "cos(A) = " << cos_angle << std::endl;
       h_reco_cos_angle->Fill(cos_angle);

       float pi0_mass = std::sqrt(2*energy1*energy2*std::abs(1 - cos_angle));
       std::cout << "reco invarient pi0 mass = " << pi0_mass << std::endl;
       h_reco_inv_pi0_mass->Fill(pi0_mass);
       reco_pi0_masses.push_back(pi0_mass);
    }



    if (recoPfpPassId.size() > 1) {
      std::sort(recoPfpPassId.begin(), recoPfpPassId.end());
      auto i1 = std::adjacent_find(recoPfpPassId.begin(), recoPfpPassId.end());
      bool isUnique = (i1 == recoPfpPassId.end());
      if (isUnique == 0) {
        trackRepeats++;
      }  //is unique
    } //recopassID

   std::cout << "------------------------------------" << std::endl;

  eventNum++;
  }// Event loop



  std::cout << " --- 100 % --- |" << std::endl;


  //Calculate the efficiency and purity of the selection
  float purity = 0;
  float efficiency = 100;
  if (recoSelectedPhotons != 0)
    purity = ((float)recoSelectedSignalPi0s/(float)recoSelectedPhotons)*100;
  if (trueSignalPi0s != 0)
    efficiency = ((float)recoSelectedSignalPi0s/(float)trueSignalPi0s)*100;

  // Print Stats
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Results..." << std::endl;
  std::cout << " True Particle #        = " << totalParticlesTrue << std::endl;
  std::cout << " True Signal #          = " << trueSignalPi0s << std::endl;
  std::cout << " Reco Particles #       = " << totalParticlesReco << std::endl;
  std::cout << " Reco Selected #        = " << recoSelectedPhotons << std::endl;
  std::cout << " Reco Selected Signal # = " << recoSelectedSignalPi0s << std::endl;
  std::cout << " Selection Purity       = " << purity << std::endl;
  std::cout << " Selection Efficiency   = " << efficiency << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;

  // Now write the histograms if wanted
  TCanvas *c1 = new TCanvas();
  h_truth_inv_pi0_mass->Draw("hist");
  h_truth_inv_pi0_mass->GetXaxis()->SetTitle("Pi0 Mass (Truth) [GeV]");
  h_truth_inv_pi0_mass->GetYaxis()->SetTitle("Counts");
  c1->SaveAs("truth_inv_pi0_mass.png");
  c1->Clear();

  h_reco_inv_pi0_mass->Draw("hist");
  h_reco_inv_pi0_mass->GetXaxis()->SetTitle("Pi0 Mass (Reco) [GeV]");
  h_reco_inv_pi0_mass->GetYaxis()->SetTitle("Counts");
  c1->SaveAs("reco_inv_pi0_mass.png");
  c1->Clear();

  h_truth_shower_energy->Draw("hist");
  h_truth_shower_energy->GetXaxis()->SetTitle("Shower Energy (Truth) [GeV]");
  h_truth_shower_energy->GetYaxis()->SetTitle("Counts");
  c1->SaveAs("truth_shower_energy.png");
  c1->Clear();

  h_reco_shower_energy->Draw("hist");
  h_reco_shower_energy->GetXaxis()->SetTitle("Shower Energy (Reco) [GeV]");
  h_reco_shower_energy->GetYaxis()->SetTitle("Counts");
  c1->SaveAs("reco_shower_energy.png");
  c1->Clear();

  h_truth_cos_angle->Draw("hist");
  h_truth_cos_angle->GetXaxis()->SetTitle("Cos Opening Angle (Truth)");
  h_truth_cos_angle->GetYaxis()->SetTitle("Counts");
  c1->SaveAs("truth_cos_angle.png");
  c1->Clear();

  h_reco_cos_angle->Draw("hist");
  h_reco_cos_angle->GetXaxis()->SetTitle("Cos Opening Angle (Reco)");
  h_reco_cos_angle->GetYaxis()->SetTitle("Counts");
  c1->SaveAs("reco_cos_angle.png");
  c1->Clear();

  h_truth_dirX->Draw("hist");
  h_truth_dirX->GetXaxis()->SetTitle("Shower Direction X (Truth)");
  h_truth_dirX->GetYaxis()->SetTitle("Counts");
  c1->SaveAs("truth_dirX.png");
  c1->Clear();

  h_reco_dirX->Draw("hist");
  h_reco_dirX->GetXaxis()->SetTitle("Shower Direction X (Reco)");
  h_reco_dirX->GetYaxis()->SetTitle("Counts");
  c1->SaveAs("reco_dirX.png");
  c1->Clear();

  h_truth_dirY->Draw("hist");
  h_truth_dirY->GetXaxis()->SetTitle("Shower Direction Y (Truth)");
  h_truth_dirY->GetYaxis()->SetTitle("Counts");
  c1->SaveAs("truth_dirY.png");
  c1->Clear();

  h_reco_dirY->Draw("hist");
  h_reco_dirY->GetXaxis()->SetTitle("Shower Direction Y (Reco)");
  h_reco_dirY->GetYaxis()->SetTitle("Counts");
  c1->SaveAs("reco_dirY.png");
  c1->Clear();

  h_truth_dirZ->Draw("hist");
  h_truth_dirZ->GetXaxis()->SetTitle("Shower Direction Z (Truth)");
  h_truth_dirZ->GetYaxis()->SetTitle("Counts");
  c1->SaveAs("truth_dirZ.png");
  c1->Clear();

  h_reco_dirZ->Draw("hist");
  h_reco_dirZ->GetXaxis()->SetTitle("Shower Direction Z (Reco)");
  h_reco_dirZ->GetYaxis()->SetTitle("Counts");
  c1->SaveAs("reco_dirZ.png");
  c1->Clear();


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
