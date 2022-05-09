//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Apr 25 15:27:04 2022 by ROOT version 6.22/08
// from TTree anatree/analysis tree
// found on file: /data/jones/dune/calibration/production_cr_chain/dunetpc_v09_41_00_02/20006000/20006000/MUSUN_dunefd_20006000_gen_g4_detsim_reco_ana.root
//////////////////////////////////////////////////////////

#ifndef anatree_h
#define anatree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

const int MAX_HITS = 100000;
const int MAX_VERTICES = 100;
const int MAX_CLUSTERS = 10000;
const int MAX_MCEVTS = 10;
const int MAX_GENIE_PRIMARIES = 10;
const int MAX_G4_PARTICLES = 250000;
const int MAX_SHOWERS = 100;
const int MAX_TRACKS = 5000;
const int MAX_TRACK_HITS = 2000;
const int MAX_MCTRACK_HITS = 2000;
const int MAX_MCTRACKS = 5000;
const int MAX_MCSHOWERS = 5000;
const int MAX_PFPARTICLES = 5000;

class anatree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Int_t           subrun;
   Int_t           event;
   Double_t        evttime;
   Double_t        beamtime;
   Double_t        pot;
   Char_t          isdata;
   Double_t        taulife;
   UInt_t          triggernumber;
   Double_t        triggertime;
   Double_t        beamgatetime;
   UInt_t          triggerbits;
   Double_t        potbnb;
   Double_t        potnumitgt;
   Double_t        potnumi101;
   Int_t           no_hits;
   Int_t           no_hits_stored;
   Short_t         hit_tpc[MAX_HITS];   //[no_hits_stored]
   Short_t         hit_plane[MAX_HITS];   //[no_hits_stored]
   Short_t         hit_wire[MAX_HITS];   //[no_hits_stored]
   Short_t         hit_channel[MAX_HITS];   //[no_hits_stored]
   Float_t         hit_peakT[MAX_HITS];   //[no_hits_stored]
   Float_t         hit_charge[MAX_HITS];   //[no_hits_stored]
   Float_t         hit_ph[MAX_HITS];   //[no_hits_stored]
   Float_t         hit_startT[MAX_HITS];   //[no_hits_stored]
   Float_t         hit_endT[MAX_HITS];   //[no_hits_stored]
   Float_t         hit_rms[MAX_HITS];   //[no_hits_stored]
   Float_t         hit_trueX[MAX_HITS];   //[no_hits_stored]
   Float_t         hit_goodnessOfFit[MAX_HITS];   //[no_hits_stored]
   Short_t         hit_multiplicity[MAX_HITS];   //[no_hits_stored]
   Short_t         hit_trkid[MAX_HITS];   //[no_hits_stored]
   Short_t         hit_trkKey[MAX_HITS];   //[no_hits_stored]
   Short_t         hit_clusterid[MAX_HITS];   //[no_hits_stored]
   Short_t         hit_clusterKey[MAX_HITS];   //[no_hits_stored]
   Short_t         hit_spacepointid[MAX_HITS];   //[no_hits_stored]
   Short_t         hit_spacepointKey[MAX_HITS];   //[no_hits_stored]
   Float_t         hit_nelec[MAX_HITS];   //[no_hits_stored]
   Float_t         hit_energy[MAX_HITS];   //[no_hits_stored]
   Short_t         nnuvtx;
   Float_t         nuvtxx[MAX_VERTICES];   //[nnuvtx]
   Float_t         nuvtxy[MAX_VERTICES];   //[nnuvtx]
   Float_t         nuvtxz[MAX_VERTICES];   //[nnuvtx]
   Short_t         nuvtxpdg[MAX_VERTICES];   //[nnuvtx]
   Short_t         nclusters;
   Short_t         clusterId[MAX_CLUSTERS];   //[nclusters]
   Short_t         clusterView[MAX_CLUSTERS];   //[nclusters]
   Float_t         cluster_StartCharge[MAX_CLUSTERS];   //[nclusters]
   Float_t         cluster_StartAngle[MAX_CLUSTERS];   //[nclusters]
   Float_t         cluster_EndCharge[MAX_CLUSTERS];   //[nclusters]
   Float_t         cluster_EndAngle[MAX_CLUSTERS];   //[nclusters]
   Float_t         cluster_Integral[MAX_CLUSTERS];   //[nclusters]
   Float_t         cluster_IntegralAverage[MAX_CLUSTERS];   //[nclusters]
   Float_t         cluster_SummedADC[MAX_CLUSTERS];   //[nclusters]
   Float_t         cluster_SummedADCaverage[MAX_CLUSTERS];   //[nclusters]
   Float_t         cluster_MultipleHitDensity[MAX_CLUSTERS];   //[nclusters]
   Float_t         cluster_Width[MAX_CLUSTERS];   //[nclusters]
   Short_t         cluster_NHits[MAX_CLUSTERS];   //[nclusters]
   Short_t         cluster_StartWire[MAX_CLUSTERS];   //[nclusters]
   Short_t         cluster_StartTick[MAX_CLUSTERS];   //[nclusters]
   Short_t         cluster_EndWire[MAX_CLUSTERS];   //[nclusters]
   Short_t         cluster_EndTick[MAX_CLUSTERS];   //[nclusters]
   Short_t         cluncosmictags_tagger[MAX_CLUSTERS];   //[nclusters]
   Float_t         clucosmicscore_tagger[MAX_CLUSTERS];   //[nclusters]
   Short_t         clucosmictype_tagger[MAX_CLUSTERS];   //[nclusters]
   Char_t          kNTracker;
   Char_t          kNVertexAlgos;
   Char_t          kNShowerAlgos;
   Int_t           mcevts_truth;
   Int_t           nuPDG_truth[MAX_MCEVTS];   //[mcevts_truth]
   Int_t           ccnc_truth[MAX_MCEVTS];   //[mcevts_truth]
   Int_t           mode_truth[MAX_MCEVTS];   //[mcevts_truth]
   Float_t         enu_truth[MAX_MCEVTS];   //[mcevts_truth]
   Float_t         Q2_truth[MAX_MCEVTS];   //[mcevts_truth]
   Float_t         W_truth[MAX_MCEVTS];   //[mcevts_truth]
   Float_t         X_truth[MAX_MCEVTS];   //[mcevts_truth]
   Float_t         Y_truth[MAX_MCEVTS];   //[mcevts_truth]
   Int_t           hitnuc_truth[MAX_MCEVTS];   //[mcevts_truth]
   Float_t         nuvtxx_truth[MAX_MCEVTS];   //[mcevts_truth]
   Float_t         nuvtxy_truth[MAX_MCEVTS];   //[mcevts_truth]
   Float_t         nuvtxz_truth[MAX_MCEVTS];   //[mcevts_truth]
   Float_t         nu_dcosx_truth[MAX_MCEVTS];   //[mcevts_truth]
   Float_t         nu_dcosy_truth[MAX_MCEVTS];   //[mcevts_truth]
   Float_t         nu_dcosz_truth[MAX_MCEVTS];   //[mcevts_truth]
   Float_t         lep_mom_truth[MAX_MCEVTS];   //[mcevts_truth]
   Float_t         lep_dcosx_truth[MAX_MCEVTS];   //[mcevts_truth]
   Float_t         lep_dcosy_truth[MAX_MCEVTS];   //[mcevts_truth]
   Float_t         lep_dcosz_truth[MAX_MCEVTS];   //[mcevts_truth]
   Float_t         vx_flux[MAX_MCEVTS];   //[mcevts_truth]
   Float_t         vy_flux[MAX_MCEVTS];   //[mcevts_truth]
   Float_t         vz_flux[MAX_MCEVTS];   //[mcevts_truth]
   Float_t         pdpx_flux[MAX_MCEVTS];   //[mcevts_truth]
   Float_t         pdpy_flux[MAX_MCEVTS];   //[mcevts_truth]
   Float_t         pdpz_flux[MAX_MCEVTS];   //[mcevts_truth]
   Float_t         ppdxdz_flux[MAX_MCEVTS];   //[mcevts_truth]
   Float_t         ppdydz_flux[MAX_MCEVTS];   //[mcevts_truth]
   Float_t         pppz_flux[MAX_MCEVTS];   //[mcevts_truth]
   Int_t           ptype_flux[MAX_MCEVTS];   //[mcevts_truth]
   Float_t         ppvx_flux[MAX_MCEVTS];   //[mcevts_truth]
   Float_t         ppvy_flux[MAX_MCEVTS];   //[mcevts_truth]
   Float_t         ppvz_flux[MAX_MCEVTS];   //[mcevts_truth]
   Float_t         muparpx_flux[MAX_MCEVTS];   //[mcevts_truth]
   Float_t         muparpy_flux[MAX_MCEVTS];   //[mcevts_truth]
   Float_t         muparpz_flux[MAX_MCEVTS];   //[mcevts_truth]
   Float_t         mupare_flux[MAX_MCEVTS];   //[mcevts_truth]
   Int_t           tgen_flux[MAX_MCEVTS];   //[mcevts_truth]
   Int_t           tgptype_flux[MAX_MCEVTS];   //[mcevts_truth]
   Float_t         tgppx_flux[MAX_MCEVTS];   //[mcevts_truth]
   Float_t         tgppy_flux[MAX_MCEVTS];   //[mcevts_truth]
   Float_t         tgppz_flux[MAX_MCEVTS];   //[mcevts_truth]
   Float_t         tprivx_flux[MAX_MCEVTS];   //[mcevts_truth]
   Float_t         tprivy_flux[MAX_MCEVTS];   //[mcevts_truth]
   Float_t         tprivz_flux[MAX_MCEVTS];   //[mcevts_truth]
   Float_t         dk2gen_flux[MAX_MCEVTS];   //[mcevts_truth]
   Float_t         gen2vtx_flux[MAX_MCEVTS];   //[mcevts_truth]
   Float_t         tpx_flux[MAX_MCEVTS];   //[mcevts_truth]
   Float_t         tpy_flux[MAX_MCEVTS];   //[mcevts_truth]
   Float_t         tpz_flux[MAX_MCEVTS];   //[mcevts_truth]
   Int_t           tptype_flux[MAX_MCEVTS];   //[mcevts_truth]
   Int_t           genie_no_primaries;
   Int_t           genie_primaries_pdg[MAX_GENIE_PRIMARIES];   //[genie_no_primaries]
   Float_t         genie_Eng[MAX_GENIE_PRIMARIES];   //[genie_no_primaries]
   Float_t         genie_Px[MAX_GENIE_PRIMARIES];   //[genie_no_primaries]
   Float_t         genie_Py[MAX_GENIE_PRIMARIES];   //[genie_no_primaries]
   Float_t         genie_Pz[MAX_GENIE_PRIMARIES];   //[genie_no_primaries]
   Float_t         genie_P[MAX_GENIE_PRIMARIES];   //[genie_no_primaries]
   Int_t           genie_status_code[MAX_GENIE_PRIMARIES];   //[genie_no_primaries]
   Float_t         genie_mass[MAX_GENIE_PRIMARIES];   //[genie_no_primaries]
   Int_t           genie_trackID[MAX_GENIE_PRIMARIES];   //[genie_no_primaries]
   Int_t           genie_ND[MAX_GENIE_PRIMARIES];   //[genie_no_primaries]
   Int_t           genie_mother[MAX_GENIE_PRIMARIES];   //[genie_no_primaries]
   Int_t           no_primaries;
   Int_t           geant_list_size;
   Int_t           geant_list_size_in_tpcAV;
   Int_t           pdg[MAX_G4_PARTICLES];   //[geant_list_size]
   Int_t           status[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         Mass[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         Eng[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         EndE[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         Px[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         Py[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         Pz[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         P[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         StartPointx[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         StartPointy[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         StartPointz[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         StartT[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         EndPointx[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         EndPointy[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         EndPointz[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         EndT[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         theta[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         phi[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         theta_xz[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         theta_yz[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         pathlen[MAX_G4_PARTICLES];   //[geant_list_size]
   Int_t           inTPCActive[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         StartPointx_tpcAV[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         StartPointy_tpcAV[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         StartPointz_tpcAV[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         StartT_tpcAV[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         StartE_tpcAV[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         StartP_tpcAV[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         StartPx_tpcAV[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         StartPy_tpcAV[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         StartPz_tpcAV[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         EndPointx_tpcAV[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         EndPointy_tpcAV[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         EndPointz_tpcAV[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         EndT_tpcAV[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         EndE_tpcAV[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         EndP_tpcAV[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         EndPx_tpcAV[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         EndPy_tpcAV[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         EndPz_tpcAV[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         pathlen_drifted[MAX_G4_PARTICLES];   //[geant_list_size]
   Int_t           inTPCDrifted[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         StartPointx_drifted[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         StartPointy_drifted[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         StartPointz_drifted[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         StartT_drifted[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         StartE_drifted[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         StartP_drifted[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         StartPx_drifted[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         StartPy_drifted[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         StartPz_drifted[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         EndPointx_drifted[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         EndPointy_drifted[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         EndPointz_drifted[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         EndT_drifted[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         EndE_drifted[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         EndP_drifted[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         EndPx_drifted[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         EndPy_drifted[MAX_G4_PARTICLES];   //[geant_list_size]
   Float_t         EndPz_drifted[MAX_G4_PARTICLES];   //[geant_list_size]
   Int_t           NumberDaughters[MAX_G4_PARTICLES];   //[geant_list_size]
   Int_t           Mother[MAX_G4_PARTICLES];   //[geant_list_size]
   Int_t           TrackId[MAX_G4_PARTICLES];   //[geant_list_size]
   Int_t           MergedId[MAX_G4_PARTICLES];   //[geant_list_size]
   Int_t           origin[MAX_G4_PARTICLES];   //[geant_list_size]
   Int_t           MCTruthIndex[MAX_G4_PARTICLES];   //[geant_list_size]
   Int_t           process_primary[MAX_G4_PARTICLES];   //[geant_list_size]
   vector<string>  *processname;
   Short_t         nPFParticles;
   Short_t         pfp_selfID[MAX_PFPARTICLES];   //[nPFParticles]
   Short_t         pfp_isPrimary[MAX_PFPARTICLES];   //[nPFParticles]
   Short_t         pfp_numDaughters[MAX_PFPARTICLES];   //[nPFParticles]
   Short_t         pfp_daughterIDs[MAX_PFPARTICLES][10];   //[nPFParticles]
   Short_t         pfp_parentID[MAX_PFPARTICLES];   //[nPFParticles]
   Short_t         pfp_vertexID[MAX_PFPARTICLES];   //[nPFParticles]
   Short_t         pfp_isShower[MAX_PFPARTICLES];   //[nPFParticles]
   Short_t         pfp_isTrack[MAX_PFPARTICLES];   //[nPFParticles]
   Short_t         pfp_trackID[MAX_PFPARTICLES];   //[nPFParticles]
   Short_t         pfp_showerID[MAX_PFPARTICLES];   //[nPFParticles]
   Int_t           pfp_pdgCode[MAX_PFPARTICLES];   //[nPFParticles]
   Short_t         pfp_numClusters[MAX_PFPARTICLES];   //[nPFParticles]
   Short_t         pfp_clusterIDs[MAX_PFPARTICLES][10];   //[nPFParticles]
   Short_t         pfp_isNeutrino[MAX_PFPARTICLES];   //[nPFParticles]
   Short_t         pfp_numNeutrinos;
   Short_t         pfp_neutrinoIDs[5];
   Short_t         nshowers_emshower;
   Short_t         showerID_emshower[MAX_SHOWERS];   //[nshowers_emshower]
   Short_t         shwr_bestplane_emshower[MAX_SHOWERS];   //[nshowers_emshower]
   Float_t         shwr_length_emshower[MAX_SHOWERS];   //[nshowers_emshower]
   Float_t         shwr_startdcosx_emshower[MAX_SHOWERS];   //[nshowers_emshower]
   Float_t         shwr_startdcosy_emshower[MAX_SHOWERS];   //[nshowers_emshower]
   Float_t         shwr_startdcosz_emshower[MAX_SHOWERS];   //[nshowers_emshower]
   Float_t         shwr_startx_emshower[MAX_SHOWERS];   //[nshowers_emshower]
   Float_t         shwr_starty_emshower[MAX_SHOWERS];   //[nshowers_emshower]
   Float_t         shwr_startz_emshower[MAX_SHOWERS];   //[nshowers_emshower]
   Float_t         shwr_totEng_emshower[MAX_SHOWERS][3];   //[nshowers_emshower]
   Float_t         shwr_dedx_emshower[MAX_SHOWERS][3];   //[nshowers_emshower]
   Float_t         shwr_mipEng_emshower[MAX_SHOWERS][3];   //[nshowers_emshower]
   Short_t         shwr_hasPFParticle_emshower[MAX_SHOWERS];   //[nshowers_emshower]
   Short_t         shwr_PFParticleID_emshower[MAX_SHOWERS];   //[nshowers_emshower]
   Float_t         shwr_pidmvamu_emshower[MAX_SHOWERS];   //[nshowers_emshower]
   Float_t         shwr_pidmvae_emshower[MAX_SHOWERS];   //[nshowers_emshower]
   Float_t         shwr_pidmvapich_emshower[MAX_SHOWERS];   //[nshowers_emshower]
   Float_t         shwr_pidmvaphoton_emshower[MAX_SHOWERS];   //[nshowers_emshower]
   Float_t         shwr_pidmvapr_emshower[MAX_SHOWERS];   //[nshowers_emshower]
   Short_t         ntracks_pandoraTrack;
   Short_t         trkId_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Short_t         trkncosmictags_tagger_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Float_t         trkcosmicscore_tagger_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Short_t         trkcosmictype_tagger_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Short_t         trkncosmictags_containmenttagger_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Float_t         trkcosmicscore_containmenttagger_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Short_t         trkcosmictype_containmenttagger_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Short_t         trkncosmictags_flashmatch_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Float_t         trkcosmicscore_flashmatch_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Short_t         trkcosmictype_flashmatch_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Float_t         trkke_pandoraTrack[MAX_TRACKS][3];   //[ntracks_pandoraTrack]
   Float_t         trkrange_pandoraTrack[MAX_TRACKS][3];   //[ntracks_pandoraTrack]
   Int_t           trkidtruth_pandoraTrack[MAX_TRACKS][3];   //[ntracks_pandoraTrack]
   Short_t         trkorigin_pandoraTrack[MAX_TRACKS][3];   //[ntracks_pandoraTrack]
   Int_t           trkpdgtruth_pandoraTrack[MAX_TRACKS][3];   //[ntracks_pandoraTrack]
   Float_t         trkefftruth_pandoraTrack[MAX_TRACKS][3];   //[ntracks_pandoraTrack]
   Float_t         trkpurtruth_pandoraTrack[MAX_TRACKS][3];   //[ntracks_pandoraTrack]
   Float_t         trkpitchc_pandoraTrack[MAX_TRACKS][3];   //[ntracks_pandoraTrack]
   Short_t         ntrkhits_pandoraTrack[MAX_TRACKS][3];   //[ntracks_pandoraTrack]
   Float_t         trkdedx_pandoraTrack[MAX_TRACKS][3][MAX_TRACK_HITS];   //[ntracks_pandoraTrack]
   Float_t         trkdqdx_pandoraTrack[MAX_TRACKS][3][MAX_TRACK_HITS];   //[ntracks_pandoraTrack]
   Float_t         trkresrg_pandoraTrack[MAX_TRACKS][3][MAX_TRACK_HITS];   //[ntracks_pandoraTrack]
   Int_t           trktpc_pandoraTrack[MAX_TRACKS][3][MAX_TRACK_HITS];   //[ntracks_pandoraTrack]
   Float_t         trkxyz_pandoraTrack[MAX_TRACKS][3][MAX_TRACK_HITS][3];   //[ntracks_pandoraTrack]
   Float_t         trkstartx_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Float_t         trkstarty_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Float_t         trkstartz_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Float_t         trkstartd_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Float_t         trkendx_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Float_t         trkendy_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Float_t         trkendz_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Float_t         trkendd_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Float_t         trkflashT0_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Float_t         trktrueT0_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Int_t           trkg4id_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Int_t           trkorig_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Float_t         trkpurity_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Float_t         trkcompleteness_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Float_t         trktheta_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Float_t         trkphi_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Float_t         trkstartdcosx_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Float_t         trkstartdcosy_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Float_t         trkstartdcosz_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Float_t         trkenddcosx_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Float_t         trkenddcosy_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Float_t         trkenddcosz_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Float_t         trkthetaxz_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Float_t         trkthetayz_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Float_t         trkmom_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Float_t         trkmomrange_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Float_t         trkmommschi2_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Float_t         trkmommsllhd_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Float_t         trklen_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Short_t         trksvtxid_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Short_t         trkevtxid_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Float_t         trkpidmvamu_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Float_t         trkpidmvae_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Float_t         trkpidmvapich_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Float_t         trkpidmvaphoton_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Float_t         trkpidmvapr_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Int_t           trkpidpdg_pandoraTrack[MAX_TRACKS][3];   //[ntracks_pandoraTrack]
   Float_t         trkpidchi_pandoraTrack[MAX_TRACKS][3];   //[ntracks_pandoraTrack]
   Float_t         trkpidchipr_pandoraTrack[MAX_TRACKS][3];   //[ntracks_pandoraTrack]
   Float_t         trkpidchika_pandoraTrack[MAX_TRACKS][3];   //[ntracks_pandoraTrack]
   Float_t         trkpidchipi_pandoraTrack[MAX_TRACKS][3];   //[ntracks_pandoraTrack]
   Float_t         trkpidchimu_pandoraTrack[MAX_TRACKS][3];   //[ntracks_pandoraTrack]
   Float_t         trkpidpida_pandoraTrack[MAX_TRACKS][3];   //[ntracks_pandoraTrack]
   Short_t         trkpidbestplane_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Short_t         trkhasPFParticle_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Short_t         trkPFParticleID_pandoraTrack[MAX_TRACKS];   //[ntracks_pandoraTrack]
   Short_t         nvtx_linecluster;
   Short_t         vtxId_linecluster[MAX_VERTICES];   //[nvtx_linecluster]
   Float_t         vtxx_linecluster[MAX_VERTICES];   //[nvtx_linecluster]
   Float_t         vtxy_linecluster[MAX_VERTICES];   //[nvtx_linecluster]
   Float_t         vtxz_linecluster[MAX_VERTICES];   //[nvtx_linecluster]
   Short_t         vtxhasPFParticle_linecluster[MAX_VERTICES];   //[nvtx_linecluster]
   Short_t         vtxPFParticleID_linecluster[MAX_VERTICES];   //[nvtx_linecluster]
   Short_t         nvtx_pandora;
   Short_t         vtxId_pandora[MAX_VERTICES];   //[nvtx_pandora]
   Float_t         vtxx_pandora[MAX_VERTICES];   //[nvtx_pandora]
   Float_t         vtxy_pandora[MAX_VERTICES];   //[nvtx_pandora]
   Float_t         vtxz_pandora[MAX_VERTICES];   //[nvtx_pandora]
   Short_t         vtxhasPFParticle_pandora[MAX_VERTICES];   //[nvtx_pandora]
   Short_t         vtxPFParticleID_pandora[MAX_VERTICES];   //[nvtx_pandora]
   Int_t           no_mcshowers;
   Int_t           mcshwr_origin[MAX_MCSHOWERS];   //[no_mcshowers]
   Int_t           mcshwr_pdg[MAX_MCSHOWERS];   //[no_mcshowers]
   Int_t           mcshwr_TrackId[MAX_MCSHOWERS];   //[no_mcshowers]
   vector<string>  *mcshwr_Process;
   Float_t         mcshwr_startX[MAX_MCSHOWERS];   //[no_mcshowers]
   Float_t         mcshwr_startY[MAX_MCSHOWERS];   //[no_mcshowers]
   Float_t         mcshwr_startZ[MAX_MCSHOWERS];   //[no_mcshowers]
   Float_t         mcshwr_endX[MAX_MCSHOWERS];   //[no_mcshowers]
   Float_t         mcshwr_endY[MAX_MCSHOWERS];   //[no_mcshowers]
   Float_t         mcshwr_endZ[MAX_MCSHOWERS];   //[no_mcshowers]
   Float_t         mcshwr_CombEngX[MAX_MCSHOWERS];   //[no_mcshowers]
   Float_t         mcshwr_CombEngY[MAX_MCSHOWERS];   //[no_mcshowers]
   Float_t         mcshwr_CombEngZ[MAX_MCSHOWERS];   //[no_mcshowers]
   Float_t         mcshwr_CombEngPx[MAX_MCSHOWERS];   //[no_mcshowers]
   Float_t         mcshwr_CombEngPy[MAX_MCSHOWERS];   //[no_mcshowers]
   Float_t         mcshwr_CombEngPz[MAX_MCSHOWERS];   //[no_mcshowers]
   Float_t         mcshwr_CombEngE[MAX_MCSHOWERS];   //[no_mcshowers]
   Float_t         mcshwr_dEdx[MAX_MCSHOWERS];   //[no_mcshowers]
   Float_t         mcshwr_dQdx[MAX_MCSHOWERS][3];   //[no_mcshowers]
   Float_t         mcshwr_StartDirX[MAX_MCSHOWERS];   //[no_mcshowers]
   Float_t         mcshwr_StartDirY[MAX_MCSHOWERS];   //[no_mcshowers]
   Float_t         mcshwr_StartDirZ[MAX_MCSHOWERS];   //[no_mcshowers]
   Int_t           mcshwr_isEngDeposited[MAX_MCSHOWERS];   //[no_mcshowers]
   Int_t           mcshwr_Motherpdg[MAX_MCSHOWERS];   //[no_mcshowers]
   Int_t           mcshwr_MotherTrkId[MAX_MCSHOWERS];   //[no_mcshowers]
   vector<string>  *mcshwr_MotherProcess;
   Float_t         mcshwr_MotherstartX[MAX_MCSHOWERS];   //[no_mcshowers]
   Float_t         mcshwr_MotherstartY[MAX_MCSHOWERS];   //[no_mcshowers]
   Float_t         mcshwr_MotherstartZ[MAX_MCSHOWERS];   //[no_mcshowers]
   Float_t         mcshwr_MotherendX[MAX_MCSHOWERS];   //[no_mcshowers]
   Float_t         mcshwr_MotherendY[MAX_MCSHOWERS];   //[no_mcshowers]
   Float_t         mcshwr_MotherendZ[MAX_MCSHOWERS];   //[no_mcshowers]
   Int_t           mcshwr_Ancestorpdg[MAX_MCSHOWERS];   //[no_mcshowers]
   Int_t           mcshwr_AncesotorTrkId[MAX_MCSHOWERS];   //[no_mcshowers]
   vector<string>  *mcshwr_AncesotorProcess;
   Float_t         mcshwr_AncestorstartX[MAX_MCSHOWERS];   //[no_mcshowers]
   Float_t         mcshwr_AncestorstartY[MAX_MCSHOWERS];   //[no_mcshowers]
   Float_t         mcshwr_AncestorstartZ[MAX_MCSHOWERS];   //[no_mcshowers]
   Float_t         mcshwr_AncestorendX[MAX_MCSHOWERS];   //[no_mcshowers]
   Float_t         mcshwr_AncestorendY[MAX_MCSHOWERS];   //[no_mcshowers]
   Float_t         mcshwr_AncestorendZ[MAX_MCSHOWERS];   //[no_mcshowers]
   Int_t           no_mctracks;
   Int_t           mctrk_origin[MAX_MCTRACKS];   //[no_mctracks]
   Int_t           mctrk_pdg[MAX_MCTRACKS];   //[no_mctracks]
   Int_t           mctrk_TrackId[MAX_MCTRACKS];   //[no_mctracks]
   vector<string>  *mctrk_Process;
   Float_t         mctrk_startX[MAX_MCTRACKS];   //[no_mctracks]
   Float_t         mctrk_startY[MAX_MCTRACKS];   //[no_mctracks]
   Float_t         mctrk_startZ[MAX_MCTRACKS];   //[no_mctracks]
   Float_t         mctrk_endX[MAX_MCTRACKS];   //[no_mctracks]
   Float_t         mctrk_endY[MAX_MCTRACKS];   //[no_mctracks]
   Float_t         mctrk_endZ[MAX_MCTRACKS];   //[no_mctracks]
   Float_t         mctrk_startX_drifted[MAX_MCTRACKS];   //[no_mctracks]
   Float_t         mctrk_startY_drifted[MAX_MCTRACKS];   //[no_mctracks]
   Float_t         mctrk_startZ_drifted[MAX_MCTRACKS];   //[no_mctracks]
   Float_t         mctrk_endX_drifted[MAX_MCTRACKS];   //[no_mctracks]
   Float_t         mctrk_endY_drifted[MAX_MCTRACKS];   //[no_mctracks]
   Float_t         mctrk_endZ_drifted[MAX_MCTRACKS];   //[no_mctracks]
   Float_t         mctrk_len_drifted[MAX_MCTRACKS];   //[no_mctracks]
   Float_t         mctrk_p_drifted[MAX_MCTRACKS];   //[no_mctracks]
   Float_t         mctrk_px_drifted[MAX_MCTRACKS];   //[no_mctracks]
   Float_t         mctrk_py_drifted[MAX_MCTRACKS];   //[no_mctracks]
   Float_t         mctrk_pz_drifted[MAX_MCTRACKS];   //[no_mctracks]
   Int_t           nmctrksteps[MAX_MCTRACKS];   //[no_mctracks]
   Float_t         mctrk_dEdx[MAX_MCTRACKS][MAX_MCTRACK_HITS];   //[no_mctracks]
   Float_t         mctrk_dQdx[MAX_MCTRACKS][3][MAX_MCTRACK_HITS];   //[no_mctracks]
   Int_t           mctrk_Motherpdg[MAX_MCTRACKS];   //[no_mctracks]
   Int_t           mctrk_MotherTrkId[MAX_MCTRACKS];   //[no_mctracks]
   vector<string>  *mctrk_MotherProcess;
   Float_t         mctrk_MotherstartX[MAX_MCTRACKS];   //[no_mctracks]
   Float_t         mctrk_MotherstartY[MAX_MCTRACKS];   //[no_mctracks]
   Float_t         mctrk_MotherstartZ[MAX_MCTRACKS];   //[no_mctracks]
   Float_t         mctrk_MotherendX[MAX_MCTRACKS];   //[no_mctracks]
   Float_t         mctrk_MotherendY[MAX_MCTRACKS];   //[no_mctracks]
   Float_t         mctrk_MotherendZ[MAX_MCTRACKS];   //[no_mctracks]
   Int_t           mctrk_Ancestorpdg[MAX_MCTRACKS];   //[no_mctracks]
   Int_t           mctrk_AncesotorTrkId[MAX_MCTRACKS];   //[no_mctracks]
   vector<string>  *mctrk_AncesotorProcess;
   Float_t         mctrk_AncestorstartX[MAX_MCTRACKS];   //[no_mctracks]
   Float_t         mctrk_AncestorstartY[MAX_MCTRACKS];   //[no_mctracks]
   Float_t         mctrk_AncestorstartZ[MAX_MCTRACKS];   //[no_mctracks]
   Float_t         mctrk_AncestorendX[MAX_MCTRACKS];   //[no_mctracks]
   Float_t         mctrk_AncestorendY[MAX_MCTRACKS];   //[no_mctracks]
   Float_t         mctrk_AncestorendZ[MAX_MCTRACKS];   //[no_mctracks]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_subrun;   //!
   TBranch        *b_event;   //!
   TBranch        *b_evttime;   //!
   TBranch        *b_beamtime;   //!
   TBranch        *b_pot;   //!
   TBranch        *b_isdata;   //!
   TBranch        *b_taulife;   //!
   TBranch        *b_triggernumber;   //!
   TBranch        *b_triggertime;   //!
   TBranch        *b_beamgatetime;   //!
   TBranch        *b_triggerbits;   //!
   TBranch        *b_potbnb;   //!
   TBranch        *b_potnumitgt;   //!
   TBranch        *b_potnumi101;   //!
   TBranch        *b_no_hits;   //!
   TBranch        *b_no_hits_stored;   //!
   TBranch        *b_hit_tpc;   //!
   TBranch        *b_hit_plane;   //!
   TBranch        *b_hit_wire;   //!
   TBranch        *b_hit_channel;   //!
   TBranch        *b_hit_peakT;   //!
   TBranch        *b_hit_charge;   //!
   TBranch        *b_hit_ph;   //!
   TBranch        *b_hit_startT;   //!
   TBranch        *b_hit_endT;   //!
   TBranch        *b_hit_rms;   //!
   TBranch        *b_hit_trueX;   //!
   TBranch        *b_hit_goodnessOfFit;   //!
   TBranch        *b_hit_multiplicity;   //!
   TBranch        *b_hit_trkid;   //!
   TBranch        *b_hit_trkKey;   //!
   TBranch        *b_hit_clusterid;   //!
   TBranch        *b_hit_clusterKey;   //!
   TBranch        *b_hit_spacepointid;   //!
   TBranch        *b_hit_spacepointKey;   //!
   TBranch        *b_hit_nelec;   //!
   TBranch        *b_hit_energy;   //!
   TBranch        *b_nnuvtx;   //!
   TBranch        *b_nuvtxx;   //!
   TBranch        *b_nuvtxy;   //!
   TBranch        *b_nuvtxz;   //!
   TBranch        *b_nuvtxpdg;   //!
   TBranch        *b_nclusters;   //!
   TBranch        *b_clusterId;   //!
   TBranch        *b_clusterView;   //!
   TBranch        *b_cluster_StartCharge;   //!
   TBranch        *b_cluster_StartAngle;   //!
   TBranch        *b_cluster_EndCharge;   //!
   TBranch        *b_cluster_EndAngle;   //!
   TBranch        *b_cluster_Integral;   //!
   TBranch        *b_cluster_IntegralAverage;   //!
   TBranch        *b_cluster_SummedADC;   //!
   TBranch        *b_cluster_SummedADCaverage;   //!
   TBranch        *b_cluster_MultipleHitDensity;   //!
   TBranch        *b_cluster_Width;   //!
   TBranch        *b_cluster_NHits;   //!
   TBranch        *b_cluster_StartWire;   //!
   TBranch        *b_cluster_StartTick;   //!
   TBranch        *b_cluster_EndWire;   //!
   TBranch        *b_cluster_EndTick;   //!
   TBranch        *b_cluncosmictags_tagger;   //!
   TBranch        *b_clucosmicscore_tagger;   //!
   TBranch        *b_clucosmictype_tagger;   //!
   TBranch        *b_kNTracker;   //!
   TBranch        *b_kNVertexAlgos;   //!
   TBranch        *b_kNShowerAlgos;   //!
   TBranch        *b_mcevts_truth;   //!
   TBranch        *b_nuPDG_truth;   //!
   TBranch        *b_ccnc_truth;   //!
   TBranch        *b_mode_truth;   //!
   TBranch        *b_enu_truth;   //!
   TBranch        *b_Q2_truth;   //!
   TBranch        *b_W_truth;   //!
   TBranch        *b_X_truth;   //!
   TBranch        *b_Y_truth;   //!
   TBranch        *b_hitnuc_truth;   //!
   TBranch        *b_nuvtxx_truth;   //!
   TBranch        *b_nuvtxy_truth;   //!
   TBranch        *b_nuvtxz_truth;   //!
   TBranch        *b_nu_dcosx_truth;   //!
   TBranch        *b_nu_dcosy_truth;   //!
   TBranch        *b_nu_dcosz_truth;   //!
   TBranch        *b_lep_mom_truth;   //!
   TBranch        *b_lep_dcosx_truth;   //!
   TBranch        *b_lep_dcosy_truth;   //!
   TBranch        *b_lep_dcosz_truth;   //!
   TBranch        *b_vx_flux;   //!
   TBranch        *b_vy_flux;   //!
   TBranch        *b_vz_flux;   //!
   TBranch        *b_pdpx_flux;   //!
   TBranch        *b_pdpy_flux;   //!
   TBranch        *b_pdpz_flux;   //!
   TBranch        *b_ppdxdz_flux;   //!
   TBranch        *b_ppdydz_flux;   //!
   TBranch        *b_pppz_flux;   //!
   TBranch        *b_ptype_flux;   //!
   TBranch        *b_ppvx_flux;   //!
   TBranch        *b_ppvy_flux;   //!
   TBranch        *b_ppvz_flux;   //!
   TBranch        *b_muparpx_flux;   //!
   TBranch        *b_muparpy_flux;   //!
   TBranch        *b_muparpz_flux;   //!
   TBranch        *b_mupare_flux;   //!
   TBranch        *b_tgen_flux;   //!
   TBranch        *b_tgptype_flux;   //!
   TBranch        *b_tgppx_flux;   //!
   TBranch        *b_tgppy_flux;   //!
   TBranch        *b_tgppz_flux;   //!
   TBranch        *b_tprivx_flux;   //!
   TBranch        *b_tprivy_flux;   //!
   TBranch        *b_tprivz_flux;   //!
   TBranch        *b_dk2gen_flux;   //!
   TBranch        *b_gen2vtx_flux;   //!
   TBranch        *b_tpx_flux;   //!
   TBranch        *b_tpy_flux;   //!
   TBranch        *b_tpz_flux;   //!
   TBranch        *b_tptype_flux;   //!
   TBranch        *b_genie_no_primaries;   //!
   TBranch        *b_genie_primaries_pdg;   //!
   TBranch        *b_genie_Eng;   //!
   TBranch        *b_genie_Px;   //!
   TBranch        *b_genie_Py;   //!
   TBranch        *b_genie_Pz;   //!
   TBranch        *b_genie_P;   //!
   TBranch        *b_genie_status_code;   //!
   TBranch        *b_genie_mass;   //!
   TBranch        *b_genie_trackID;   //!
   TBranch        *b_genie_ND;   //!
   TBranch        *b_genie_mother;   //!
   TBranch        *b_no_primaries;   //!
   TBranch        *b_geant_list_size;   //!
   TBranch        *b_geant_list_size_in_tpcAV;   //!
   TBranch        *b_pdg;   //!
   TBranch        *b_status;   //!
   TBranch        *b_Mass;   //!
   TBranch        *b_Eng;   //!
   TBranch        *b_EndE;   //!
   TBranch        *b_Px;   //!
   TBranch        *b_Py;   //!
   TBranch        *b_Pz;   //!
   TBranch        *b_P;   //!
   TBranch        *b_StartPointx;   //!
   TBranch        *b_StartPointy;   //!
   TBranch        *b_StartPointz;   //!
   TBranch        *b_StartT;   //!
   TBranch        *b_EndPointx;   //!
   TBranch        *b_EndPointy;   //!
   TBranch        *b_EndPointz;   //!
   TBranch        *b_EndT;   //!
   TBranch        *b_theta;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_theta_xz;   //!
   TBranch        *b_theta_yz;   //!
   TBranch        *b_pathlen;   //!
   TBranch        *b_inTPCActive;   //!
   TBranch        *b_StartPointx_tpcAV;   //!
   TBranch        *b_StartPointy_tpcAV;   //!
   TBranch        *b_StartPointz_tpcAV;   //!
   TBranch        *b_StartT_tpcAV;   //!
   TBranch        *b_StartE_tpcAV;   //!
   TBranch        *b_StartP_tpcAV;   //!
   TBranch        *b_StartPx_tpcAV;   //!
   TBranch        *b_StartPy_tpcAV;   //!
   TBranch        *b_StartPz_tpcAV;   //!
   TBranch        *b_EndPointx_tpcAV;   //!
   TBranch        *b_EndPointy_tpcAV;   //!
   TBranch        *b_EndPointz_tpcAV;   //!
   TBranch        *b_EndT_tpcAV;   //!
   TBranch        *b_EndE_tpcAV;   //!
   TBranch        *b_EndP_tpcAV;   //!
   TBranch        *b_EndPx_tpcAV;   //!
   TBranch        *b_EndPy_tpcAV;   //!
   TBranch        *b_EndPz_tpcAV;   //!
   TBranch        *b_pathlen_drifted;   //!
   TBranch        *b_inTPCDrifted;   //!
   TBranch        *b_StartPointx_drifted;   //!
   TBranch        *b_StartPointy_drifted;   //!
   TBranch        *b_StartPointz_drifted;   //!
   TBranch        *b_StartT_drifted;   //!
   TBranch        *b_StartE_drifted;   //!
   TBranch        *b_StartP_drifted;   //!
   TBranch        *b_StartPx_drifted;   //!
   TBranch        *b_StartPy_drifted;   //!
   TBranch        *b_StartPz_drifted;   //!
   TBranch        *b_EndPointx_drifted;   //!
   TBranch        *b_EndPointy_drifted;   //!
   TBranch        *b_EndPointz_drifted;   //!
   TBranch        *b_EndT_drifted;   //!
   TBranch        *b_EndE_drifted;   //!
   TBranch        *b_EndP_drifted;   //!
   TBranch        *b_EndPx_drifted;   //!
   TBranch        *b_EndPy_drifted;   //!
   TBranch        *b_EndPz_drifted;   //!
   TBranch        *b_NumberDaughters;   //!
   TBranch        *b_Mother;   //!
   TBranch        *b_TrackId;   //!
   TBranch        *b_MergedId;   //!
   TBranch        *b_origin;   //!
   TBranch        *b_MCTruthIndex;   //!
   TBranch        *b_process_primary;   //!
   TBranch        *b_processname;   //!
   TBranch        *b_nPFParticles;   //!
   TBranch        *b_pfp_selfID;   //!
   TBranch        *b_pfp_isPrimary;   //!
   TBranch        *b_pfp_numDaughters;   //!
   TBranch        *b_pfp_daughterIDs;   //!
   TBranch        *b_pfp_parentID;   //!
   TBranch        *b_pfp_vertexID;   //!
   TBranch        *b_pfp_isShower;   //!
   TBranch        *b_pfp_isTrack;   //!
   TBranch        *b_pfp_trackID;   //!
   TBranch        *b_pfp_showerID;   //!
   TBranch        *b_pfp_pdgCode;   //!
   TBranch        *b_pfp_numClusters;   //!
   TBranch        *b_pfp_clusterIDs;   //!
   TBranch        *b_pfp_isNeutrino;   //!
   TBranch        *b_pfp_numNeutrinos;   //!
   TBranch        *b_pfp_neutrinoIDs;   //!
   TBranch        *b_nshowers_emshower;   //!
   TBranch        *b_showerID_emshower;   //!
   TBranch        *b_shwr_bestplane_emshower;   //!
   TBranch        *b_shwr_length_emshower;   //!
   TBranch        *b_shwr_startdcosx_emshower;   //!
   TBranch        *b_shwr_startdcosy_emshower;   //!
   TBranch        *b_shwr_startdcosz_emshower;   //!
   TBranch        *b_shwr_startx_emshower;   //!
   TBranch        *b_shwr_starty_emshower;   //!
   TBranch        *b_shwr_startz_emshower;   //!
   TBranch        *b_shwr_totEng_emshower;   //!
   TBranch        *b_shwr_dedx_emshower;   //!
   TBranch        *b_shwr_mipEng_emshower;   //!
   TBranch        *b_shwr_hasPFParticle_emshower;   //!
   TBranch        *b_shwr_PFParticleID_emshower;   //!
   TBranch        *b_shwr_pidmvamu_emshower;   //!
   TBranch        *b_shwr_pidmvae_emshower;   //!
   TBranch        *b_shwr_pidmvapich_emshower;   //!
   TBranch        *b_shwr_pidmvaphoton_emshower;   //!
   TBranch        *b_shwr_pidmvapr_emshower;   //!
   TBranch        *b_ntracks_pandoraTrack;   //!
   TBranch        *b_trkId_pandoraTrack;   //!
   TBranch        *b_trkncosmictags_tagger_pandoraTrack;   //!
   TBranch        *b_trkcosmicscore_tagger_pandoraTrack;   //!
   TBranch        *b_trkcosmictype_tagger_pandoraTrack;   //!
   TBranch        *b_trkncosmictags_containmenttagger_pandoraTrack;   //!
   TBranch        *b_trkcosmicscore_containmenttagger_pandoraTrack;   //!
   TBranch        *b_trkcosmictype_containmenttagger_pandoraTrack;   //!
   TBranch        *b_trkncosmictags_flashmatch_pandoraTrack;   //!
   TBranch        *b_trkcosmicscore_flashmatch_pandoraTrack;   //!
   TBranch        *b_trkcosmictype_flashmatch_pandoraTrack;   //!
   TBranch        *b_trkke_pandoraTrack;   //!
   TBranch        *b_trkrange_pandoraTrack;   //!
   TBranch        *b_trkidtruth_pandoraTrack;   //!
   TBranch        *b_trkorigin_pandoraTrack;   //!
   TBranch        *b_trkpdgtruth_pandoraTrack;   //!
   TBranch        *b_trkefftruth_pandoraTrack;   //!
   TBranch        *b_trkpurtruth_pandoraTrack;   //!
   TBranch        *b_trkpitchc_pandoraTrack;   //!
   TBranch        *b_ntrkhits_pandoraTrack;   //!
   TBranch        *b_trkdedx_pandoraTrack;   //!
   TBranch        *b_trkdqdx_pandoraTrack;   //!
   TBranch        *b_trkresrg_pandoraTrack;   //!
   TBranch        *b_trktpc_pandoraTrack;   //!
   TBranch        *b_trkxyz_pandoraTrack;   //!
   TBranch        *b_trkstartx_pandoraTrack;   //!
   TBranch        *b_trkstarty_pandoraTrack;   //!
   TBranch        *b_trkstartz_pandoraTrack;   //!
   TBranch        *b_trkstartd_pandoraTrack;   //!
   TBranch        *b_trkendx_pandoraTrack;   //!
   TBranch        *b_trkendy_pandoraTrack;   //!
   TBranch        *b_trkendz_pandoraTrack;   //!
   TBranch        *b_trkendd_pandoraTrack;   //!
   TBranch        *b_trkflashT0_pandoraTrack;   //!
   TBranch        *b_trktrueT0_pandoraTrack;   //!
   TBranch        *b_trkg4id_pandoraTrack;   //!
   TBranch        *b_trkorig_pandoraTrack;   //!
   TBranch        *b_trkpurity_pandoraTrack;   //!
   TBranch        *b_trkcompleteness_pandoraTrack;   //!
   TBranch        *b_trktheta_pandoraTrack;   //!
   TBranch        *b_trkphi_pandoraTrack;   //!
   TBranch        *b_trkstartdcosx_pandoraTrack;   //!
   TBranch        *b_trkstartdcosy_pandoraTrack;   //!
   TBranch        *b_trkstartdcosz_pandoraTrack;   //!
   TBranch        *b_trkenddcosx_pandoraTrack;   //!
   TBranch        *b_trkenddcosy_pandoraTrack;   //!
   TBranch        *b_trkenddcosz_pandoraTrack;   //!
   TBranch        *b_trkthetaxz_pandoraTrack;   //!
   TBranch        *b_trkthetayz_pandoraTrack;   //!
   TBranch        *b_trkmom_pandoraTrack;   //!
   TBranch        *b_trkmomrange_pandoraTrack;   //!
   TBranch        *b_trkmommschi2_pandoraTrack;   //!
   TBranch        *b_trkmommsllhd_pandoraTrack;   //!
   TBranch        *b_trklen_pandoraTrack;   //!
   TBranch        *b_trksvtxid_pandoraTrack;   //!
   TBranch        *b_trkevtxid_pandoraTrack;   //!
   TBranch        *b_trkpidmvamu_pandoraTrack;   //!
   TBranch        *b_trkpidmvae_pandoraTrack;   //!
   TBranch        *b_trkpidmvapich_pandoraTrack;   //!
   TBranch        *b_trkpidmvaphoton_pandoraTrack;   //!
   TBranch        *b_trkpidmvapr_pandoraTrack;   //!
   TBranch        *b_trkpidpdg_pandoraTrack;   //!
   TBranch        *b_trkpidchi_pandoraTrack;   //!
   TBranch        *b_trkpidchipr_pandoraTrack;   //!
   TBranch        *b_trkpidchika_pandoraTrack;   //!
   TBranch        *b_trkpidchipi_pandoraTrack;   //!
   TBranch        *b_trkpidchimu_pandoraTrack;   //!
   TBranch        *b_trkpidpida_pandoraTrack;   //!
   TBranch        *b_trkpidbestplane_pandoraTrack;   //!
   TBranch        *b_trkhasPFParticle_pandoraTrack;   //!
   TBranch        *b_trkPFParticleID_pandoraTrack;   //!
   TBranch        *b_nvtx_linecluster;   //!
   TBranch        *b_vtxId_linecluster;   //!
   TBranch        *b_vtxx_linecluster;   //!
   TBranch        *b_vtxy_linecluster;   //!
   TBranch        *b_vtxz_linecluster;   //!
   TBranch        *b_vtxhasPFParticle_linecluster;   //!
   TBranch        *b_vtxPFParticleID_linecluster;   //!
   TBranch        *b_nvtx_pandora;   //!
   TBranch        *b_vtxId_pandora;   //!
   TBranch        *b_vtxx_pandora;   //!
   TBranch        *b_vtxy_pandora;   //!
   TBranch        *b_vtxz_pandora;   //!
   TBranch        *b_vtxhasPFParticle_pandora;   //!
   TBranch        *b_vtxPFParticleID_pandora;   //!
   TBranch        *b_no_mcshowers;   //!
   TBranch        *b_mcshwr_origin;   //!
   TBranch        *b_mcshwr_pdg;   //!
   TBranch        *b_mcshwr_TrackId;   //!
   TBranch        *b_mcshwr_Process;   //!
   TBranch        *b_mcshwr_startX;   //!
   TBranch        *b_mcshwr_startY;   //!
   TBranch        *b_mcshwr_startZ;   //!
   TBranch        *b_mcshwr_endX;   //!
   TBranch        *b_mcshwr_endY;   //!
   TBranch        *b_mcshwr_endZ;   //!
   TBranch        *b_mcshwr_CombEngX;   //!
   TBranch        *b_mcshwr_CombEngY;   //!
   TBranch        *b_mcshwr_CombEngZ;   //!
   TBranch        *b_mcshwr_CombEngPx;   //!
   TBranch        *b_mcshwr_CombEngPy;   //!
   TBranch        *b_mcshwr_CombEngPz;   //!
   TBranch        *b_mcshwr_CombEngE;   //!
   TBranch        *b_mcshwr_dEdx;   //!
   TBranch        *b_mcshwr_dQdx;   //!
   TBranch        *b_mcshwr_StartDirX;   //!
   TBranch        *b_mcshwr_StartDirY;   //!
   TBranch        *b_mcshwr_StartDirZ;   //!
   TBranch        *b_mcshwr_isEngDeposited;   //!
   TBranch        *b_mcshwr_Motherpdg;   //!
   TBranch        *b_mcshwr_MotherTrkId;   //!
   TBranch        *b_mcshwr_MotherProcess;   //!
   TBranch        *b_mcshwr_MotherstartX;   //!
   TBranch        *b_mcshwr_MotherstartY;   //!
   TBranch        *b_mcshwr_MotherstartZ;   //!
   TBranch        *b_mcshwr_MotherendX;   //!
   TBranch        *b_mcshwr_MotherendY;   //!
   TBranch        *b_mcshwr_MotherendZ;   //!
   TBranch        *b_mcshwr_Ancestorpdg;   //!
   TBranch        *b_mcshwr_AncesotorTrkId;   //!
   TBranch        *b_mcshwr_AncesotorProcess;   //!
   TBranch        *b_mcshwr_AncestorstartX;   //!
   TBranch        *b_mcshwr_AncestorstartY;   //!
   TBranch        *b_mcshwr_AncestorstartZ;   //!
   TBranch        *b_mcshwr_AncestorendX;   //!
   TBranch        *b_mcshwr_AncestorendY;   //!
   TBranch        *b_mcshwr_AncestorendZ;   //!
   TBranch        *b_no_mctracks;   //!
   TBranch        *b_mctrk_origin;   //!
   TBranch        *b_mctrk_pdg;   //!
   TBranch        *b_mctrk_TrackId;   //!
   TBranch        *b_mctrk_Process;   //!
   TBranch        *b_mctrk_startX;   //!
   TBranch        *b_mctrk_startY;   //!
   TBranch        *b_mctrk_startZ;   //!
   TBranch        *b_mctrk_endX;   //!
   TBranch        *b_mctrk_endY;   //!
   TBranch        *b_mctrk_endZ;   //!
   TBranch        *b_mctrk_startX_drifted;   //!
   TBranch        *b_mctrk_startY_drifted;   //!
   TBranch        *b_mctrk_startZ_drifted;   //!
   TBranch        *b_mctrk_endX_drifted;   //!
   TBranch        *b_mctrk_endY_drifted;   //!
   TBranch        *b_mctrk_endZ_drifted;   //!
   TBranch        *b_mctrk_len_drifted;   //!
   TBranch        *b_mctrk_p_drifted;   //!
   TBranch        *b_mctrk_px_drifted;   //!
   TBranch        *b_mctrk_py_drifted;   //!
   TBranch        *b_mctrk_pz_drifted;   //!
   TBranch        *b_nmctrksteps;   //!
   TBranch        *b_mctrk_dEdx;   //!
   TBranch        *b_mctrk_dQdx;   //!
   TBranch        *b_mctrk_Motherpdg;   //!
   TBranch        *b_mctrk_MotherTrkId;   //!
   TBranch        *b_mctrk_MotherProcess;   //!
   TBranch        *b_mctrk_MotherstartX;   //!
   TBranch        *b_mctrk_MotherstartY;   //!
   TBranch        *b_mctrk_MotherstartZ;   //!
   TBranch        *b_mctrk_MotherendX;   //!
   TBranch        *b_mctrk_MotherendY;   //!
   TBranch        *b_mctrk_MotherendZ;   //!
   TBranch        *b_mctrk_Ancestorpdg;   //!
   TBranch        *b_mctrk_AncesotorTrkId;   //!
   TBranch        *b_mctrk_AncesotorProcess;   //!
   TBranch        *b_mctrk_AncestorstartX;   //!
   TBranch        *b_mctrk_AncestorstartY;   //!
   TBranch        *b_mctrk_AncestorstartZ;   //!
   TBranch        *b_mctrk_AncestorendX;   //!
   TBranch        *b_mctrk_AncestorendY;   //!
   TBranch        *b_mctrk_AncestorendZ;   //!

   anatree(TTree *tree=0);
   virtual ~anatree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef anatree_cxx
anatree::anatree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/data/jones/dune/calibration/production_cr_chain/dunetpc_v09_41_00_02/20006000/20006000/MUSUN_dunefd_20006000_gen_g4_detsim_reco_ana.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/data/jones/dune/calibration/production_cr_chain/dunetpc_v09_41_00_02/20006000/20006000/MUSUN_dunefd_20006000_gen_g4_detsim_reco_ana.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/data/jones/dune/calibration/production_cr_chain/dunetpc_v09_41_00_02/20006000/20006000/MUSUN_dunefd_20006000_gen_g4_detsim_reco_ana.root:/analysistree");
      dir->GetObject("anatree",tree);

   }
   Init(tree);
}

anatree::~anatree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t anatree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t anatree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void anatree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   processname = 0;
   mcshwr_Process = 0;
   mcshwr_MotherProcess = 0;
   mcshwr_AncesotorProcess = 0;
   mctrk_Process = 0;
   mctrk_MotherProcess = 0;
   mctrk_AncesotorProcess = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("subrun", &subrun, &b_subrun);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("evttime", &evttime, &b_evttime);
   fChain->SetBranchAddress("beamtime", &beamtime, &b_beamtime);
   fChain->SetBranchAddress("pot", &pot, &b_pot);
   fChain->SetBranchAddress("isdata", &isdata, &b_isdata);
   fChain->SetBranchAddress("taulife", &taulife, &b_taulife);
   fChain->SetBranchAddress("triggernumber", &triggernumber, &b_triggernumber);
   fChain->SetBranchAddress("triggertime", &triggertime, &b_triggertime);
   fChain->SetBranchAddress("beamgatetime", &beamgatetime, &b_beamgatetime);
   fChain->SetBranchAddress("triggerbits", &triggerbits, &b_triggerbits);
   fChain->SetBranchAddress("potbnb", &potbnb, &b_potbnb);
   fChain->SetBranchAddress("potnumitgt", &potnumitgt, &b_potnumitgt);
   fChain->SetBranchAddress("potnumi101", &potnumi101, &b_potnumi101);
   fChain->SetBranchAddress("no_hits", &no_hits, &b_no_hits);
   fChain->SetBranchAddress("no_hits_stored", &no_hits_stored, &b_no_hits_stored);
   fChain->SetBranchAddress("hit_tpc", hit_tpc, &b_hit_tpc);
   fChain->SetBranchAddress("hit_plane", hit_plane, &b_hit_plane);
   fChain->SetBranchAddress("hit_wire", hit_wire, &b_hit_wire);
   fChain->SetBranchAddress("hit_channel", hit_channel, &b_hit_channel);
   fChain->SetBranchAddress("hit_peakT", hit_peakT, &b_hit_peakT);
   fChain->SetBranchAddress("hit_charge", hit_charge, &b_hit_charge);
   fChain->SetBranchAddress("hit_ph", hit_ph, &b_hit_ph);
   fChain->SetBranchAddress("hit_startT", hit_startT, &b_hit_startT);
   fChain->SetBranchAddress("hit_endT", hit_endT, &b_hit_endT);
   fChain->SetBranchAddress("hit_rms", hit_rms, &b_hit_rms);
   fChain->SetBranchAddress("hit_trueX", hit_trueX, &b_hit_trueX);
   fChain->SetBranchAddress("hit_goodnessOfFit", hit_goodnessOfFit, &b_hit_goodnessOfFit);
   fChain->SetBranchAddress("hit_multiplicity", hit_multiplicity, &b_hit_multiplicity);
   fChain->SetBranchAddress("hit_trkid", hit_trkid, &b_hit_trkid);
   fChain->SetBranchAddress("hit_trkKey", hit_trkKey, &b_hit_trkKey);
   fChain->SetBranchAddress("hit_clusterid", hit_clusterid, &b_hit_clusterid);
   fChain->SetBranchAddress("hit_clusterKey", hit_clusterKey, &b_hit_clusterKey);
   fChain->SetBranchAddress("hit_spacepointid", hit_spacepointid, &b_hit_spacepointid);
   fChain->SetBranchAddress("hit_spacepointKey", hit_spacepointKey, &b_hit_spacepointKey);
   fChain->SetBranchAddress("hit_nelec", hit_nelec, &b_hit_nelec);
   fChain->SetBranchAddress("hit_energy", hit_energy, &b_hit_energy);
   fChain->SetBranchAddress("nnuvtx", &nnuvtx, &b_nnuvtx);
   fChain->SetBranchAddress("nuvtxx", nuvtxx, &b_nuvtxx);
   fChain->SetBranchAddress("nuvtxy", nuvtxy, &b_nuvtxy);
   fChain->SetBranchAddress("nuvtxz", nuvtxz, &b_nuvtxz);
   fChain->SetBranchAddress("nuvtxpdg", nuvtxpdg, &b_nuvtxpdg);
   fChain->SetBranchAddress("nclusters", &nclusters, &b_nclusters);
   fChain->SetBranchAddress("clusterId", clusterId, &b_clusterId);
   fChain->SetBranchAddress("clusterView", clusterView, &b_clusterView);
   fChain->SetBranchAddress("cluster_StartCharge", cluster_StartCharge, &b_cluster_StartCharge);
   fChain->SetBranchAddress("cluster_StartAngle", cluster_StartAngle, &b_cluster_StartAngle);
   fChain->SetBranchAddress("cluster_EndCharge", cluster_EndCharge, &b_cluster_EndCharge);
   fChain->SetBranchAddress("cluster_EndAngle", cluster_EndAngle, &b_cluster_EndAngle);
   fChain->SetBranchAddress("cluster_Integral", cluster_Integral, &b_cluster_Integral);
   fChain->SetBranchAddress("cluster_IntegralAverage", cluster_IntegralAverage, &b_cluster_IntegralAverage);
   fChain->SetBranchAddress("cluster_SummedADC", cluster_SummedADC, &b_cluster_SummedADC);
   fChain->SetBranchAddress("cluster_SummedADCaverage", cluster_SummedADCaverage, &b_cluster_SummedADCaverage);
   fChain->SetBranchAddress("cluster_MultipleHitDensity", cluster_MultipleHitDensity, &b_cluster_MultipleHitDensity);
   fChain->SetBranchAddress("cluster_Width", cluster_Width, &b_cluster_Width);
   fChain->SetBranchAddress("cluster_NHits", cluster_NHits, &b_cluster_NHits);
   fChain->SetBranchAddress("cluster_StartWire", cluster_StartWire, &b_cluster_StartWire);
   fChain->SetBranchAddress("cluster_StartTick", cluster_StartTick, &b_cluster_StartTick);
   fChain->SetBranchAddress("cluster_EndWire", cluster_EndWire, &b_cluster_EndWire);
   fChain->SetBranchAddress("cluster_EndTick", cluster_EndTick, &b_cluster_EndTick);
   fChain->SetBranchAddress("cluncosmictags_tagger", cluncosmictags_tagger, &b_cluncosmictags_tagger);
   fChain->SetBranchAddress("clucosmicscore_tagger", clucosmicscore_tagger, &b_clucosmicscore_tagger);
   fChain->SetBranchAddress("clucosmictype_tagger", clucosmictype_tagger, &b_clucosmictype_tagger);
   fChain->SetBranchAddress("kNTracker", &kNTracker, &b_kNTracker);
   fChain->SetBranchAddress("kNVertexAlgos", &kNVertexAlgos, &b_kNVertexAlgos);
   fChain->SetBranchAddress("kNShowerAlgos", &kNShowerAlgos, &b_kNShowerAlgos);
   fChain->SetBranchAddress("mcevts_truth", &mcevts_truth, &b_mcevts_truth);
   fChain->SetBranchAddress("nuPDG_truth", nuPDG_truth, &b_nuPDG_truth);
   fChain->SetBranchAddress("ccnc_truth", ccnc_truth, &b_ccnc_truth);
   fChain->SetBranchAddress("mode_truth", mode_truth, &b_mode_truth);
   fChain->SetBranchAddress("enu_truth", enu_truth, &b_enu_truth);
   fChain->SetBranchAddress("Q2_truth", Q2_truth, &b_Q2_truth);
   fChain->SetBranchAddress("W_truth", W_truth, &b_W_truth);
   fChain->SetBranchAddress("X_truth", X_truth, &b_X_truth);
   fChain->SetBranchAddress("Y_truth", Y_truth, &b_Y_truth);
   fChain->SetBranchAddress("hitnuc_truth", hitnuc_truth, &b_hitnuc_truth);
   fChain->SetBranchAddress("nuvtxx_truth", nuvtxx_truth, &b_nuvtxx_truth);
   fChain->SetBranchAddress("nuvtxy_truth", nuvtxy_truth, &b_nuvtxy_truth);
   fChain->SetBranchAddress("nuvtxz_truth", nuvtxz_truth, &b_nuvtxz_truth);
   fChain->SetBranchAddress("nu_dcosx_truth", nu_dcosx_truth, &b_nu_dcosx_truth);
   fChain->SetBranchAddress("nu_dcosy_truth", nu_dcosy_truth, &b_nu_dcosy_truth);
   fChain->SetBranchAddress("nu_dcosz_truth", nu_dcosz_truth, &b_nu_dcosz_truth);
   fChain->SetBranchAddress("lep_mom_truth", lep_mom_truth, &b_lep_mom_truth);
   fChain->SetBranchAddress("lep_dcosx_truth", lep_dcosx_truth, &b_lep_dcosx_truth);
   fChain->SetBranchAddress("lep_dcosy_truth", lep_dcosy_truth, &b_lep_dcosy_truth);
   fChain->SetBranchAddress("lep_dcosz_truth", lep_dcosz_truth, &b_lep_dcosz_truth);
   fChain->SetBranchAddress("vx_flux", vx_flux, &b_vx_flux);
   fChain->SetBranchAddress("vy_flux", vy_flux, &b_vy_flux);
   fChain->SetBranchAddress("vz_flux", vz_flux, &b_vz_flux);
   fChain->SetBranchAddress("pdpx_flux", pdpx_flux, &b_pdpx_flux);
   fChain->SetBranchAddress("pdpy_flux", pdpy_flux, &b_pdpy_flux);
   fChain->SetBranchAddress("pdpz_flux", pdpz_flux, &b_pdpz_flux);
   fChain->SetBranchAddress("ppdxdz_flux", ppdxdz_flux, &b_ppdxdz_flux);
   fChain->SetBranchAddress("ppdydz_flux", ppdydz_flux, &b_ppdydz_flux);
   fChain->SetBranchAddress("pppz_flux", pppz_flux, &b_pppz_flux);
   fChain->SetBranchAddress("ptype_flux", ptype_flux, &b_ptype_flux);
   fChain->SetBranchAddress("ppvx_flux", ppvx_flux, &b_ppvx_flux);
   fChain->SetBranchAddress("ppvy_flux", ppvy_flux, &b_ppvy_flux);
   fChain->SetBranchAddress("ppvz_flux", ppvz_flux, &b_ppvz_flux);
   fChain->SetBranchAddress("muparpx_flux", muparpx_flux, &b_muparpx_flux);
   fChain->SetBranchAddress("muparpy_flux", muparpy_flux, &b_muparpy_flux);
   fChain->SetBranchAddress("muparpz_flux", muparpz_flux, &b_muparpz_flux);
   fChain->SetBranchAddress("mupare_flux", mupare_flux, &b_mupare_flux);
   fChain->SetBranchAddress("tgen_flux", tgen_flux, &b_tgen_flux);
   fChain->SetBranchAddress("tgptype_flux", tgptype_flux, &b_tgptype_flux);
   fChain->SetBranchAddress("tgppx_flux", tgppx_flux, &b_tgppx_flux);
   fChain->SetBranchAddress("tgppy_flux", tgppy_flux, &b_tgppy_flux);
   fChain->SetBranchAddress("tgppz_flux", tgppz_flux, &b_tgppz_flux);
   fChain->SetBranchAddress("tprivx_flux", tprivx_flux, &b_tprivx_flux);
   fChain->SetBranchAddress("tprivy_flux", tprivy_flux, &b_tprivy_flux);
   fChain->SetBranchAddress("tprivz_flux", tprivz_flux, &b_tprivz_flux);
   fChain->SetBranchAddress("dk2gen_flux", dk2gen_flux, &b_dk2gen_flux);
   fChain->SetBranchAddress("gen2vtx_flux", gen2vtx_flux, &b_gen2vtx_flux);
   fChain->SetBranchAddress("tpx_flux", tpx_flux, &b_tpx_flux);
   fChain->SetBranchAddress("tpy_flux", tpy_flux, &b_tpy_flux);
   fChain->SetBranchAddress("tpz_flux", tpz_flux, &b_tpz_flux);
   fChain->SetBranchAddress("tptype_flux", tptype_flux, &b_tptype_flux);
   fChain->SetBranchAddress("genie_no_primaries", &genie_no_primaries, &b_genie_no_primaries);
   fChain->SetBranchAddress("genie_primaries_pdg", &genie_primaries_pdg, &b_genie_primaries_pdg);
   fChain->SetBranchAddress("genie_Eng", &genie_Eng, &b_genie_Eng);
   fChain->SetBranchAddress("genie_Px", &genie_Px, &b_genie_Px);
   fChain->SetBranchAddress("genie_Py", &genie_Py, &b_genie_Py);
   fChain->SetBranchAddress("genie_Pz", &genie_Pz, &b_genie_Pz);
   fChain->SetBranchAddress("genie_P", &genie_P, &b_genie_P);
   fChain->SetBranchAddress("genie_status_code", &genie_status_code, &b_genie_status_code);
   fChain->SetBranchAddress("genie_mass", &genie_mass, &b_genie_mass);
   fChain->SetBranchAddress("genie_trackID", &genie_trackID, &b_genie_trackID);
   fChain->SetBranchAddress("genie_ND", &genie_ND, &b_genie_ND);
   fChain->SetBranchAddress("genie_mother", &genie_mother, &b_genie_mother);
   fChain->SetBranchAddress("no_primaries", &no_primaries, &b_no_primaries);
   fChain->SetBranchAddress("geant_list_size", &geant_list_size, &b_geant_list_size);
   fChain->SetBranchAddress("geant_list_size_in_tpcAV", &geant_list_size_in_tpcAV, &b_geant_list_size_in_tpcAV);
   fChain->SetBranchAddress("pdg", pdg, &b_pdg);
   fChain->SetBranchAddress("status", status, &b_status);
   fChain->SetBranchAddress("Mass", Mass, &b_Mass);
   fChain->SetBranchAddress("Eng", Eng, &b_Eng);
   fChain->SetBranchAddress("EndE", EndE, &b_EndE);
   fChain->SetBranchAddress("Px", Px, &b_Px);
   fChain->SetBranchAddress("Py", Py, &b_Py);
   fChain->SetBranchAddress("Pz", Pz, &b_Pz);
   fChain->SetBranchAddress("P", P, &b_P);
   fChain->SetBranchAddress("StartPointx", StartPointx, &b_StartPointx);
   fChain->SetBranchAddress("StartPointy", StartPointy, &b_StartPointy);
   fChain->SetBranchAddress("StartPointz", StartPointz, &b_StartPointz);
   fChain->SetBranchAddress("StartT", StartT, &b_StartT);
   fChain->SetBranchAddress("EndPointx", EndPointx, &b_EndPointx);
   fChain->SetBranchAddress("EndPointy", EndPointy, &b_EndPointy);
   fChain->SetBranchAddress("EndPointz", EndPointz, &b_EndPointz);
   fChain->SetBranchAddress("EndT", EndT, &b_EndT);
   fChain->SetBranchAddress("theta", theta, &b_theta);
   fChain->SetBranchAddress("phi", phi, &b_phi);
   fChain->SetBranchAddress("theta_xz", theta_xz, &b_theta_xz);
   fChain->SetBranchAddress("theta_yz", theta_yz, &b_theta_yz);
   fChain->SetBranchAddress("pathlen", pathlen, &b_pathlen);
   fChain->SetBranchAddress("inTPCActive", inTPCActive, &b_inTPCActive);
   fChain->SetBranchAddress("StartPointx_tpcAV", StartPointx_tpcAV, &b_StartPointx_tpcAV);
   fChain->SetBranchAddress("StartPointy_tpcAV", StartPointy_tpcAV, &b_StartPointy_tpcAV);
   fChain->SetBranchAddress("StartPointz_tpcAV", StartPointz_tpcAV, &b_StartPointz_tpcAV);
   fChain->SetBranchAddress("StartT_tpcAV", StartT_tpcAV, &b_StartT_tpcAV);
   fChain->SetBranchAddress("StartE_tpcAV", StartE_tpcAV, &b_StartE_tpcAV);
   fChain->SetBranchAddress("StartP_tpcAV", StartP_tpcAV, &b_StartP_tpcAV);
   fChain->SetBranchAddress("StartPx_tpcAV", StartPx_tpcAV, &b_StartPx_tpcAV);
   fChain->SetBranchAddress("StartPy_tpcAV", StartPy_tpcAV, &b_StartPy_tpcAV);
   fChain->SetBranchAddress("StartPz_tpcAV", StartPz_tpcAV, &b_StartPz_tpcAV);
   fChain->SetBranchAddress("EndPointx_tpcAV", EndPointx_tpcAV, &b_EndPointx_tpcAV);
   fChain->SetBranchAddress("EndPointy_tpcAV", EndPointy_tpcAV, &b_EndPointy_tpcAV);
   fChain->SetBranchAddress("EndPointz_tpcAV", EndPointz_tpcAV, &b_EndPointz_tpcAV);
   fChain->SetBranchAddress("EndT_tpcAV", EndT_tpcAV, &b_EndT_tpcAV);
   fChain->SetBranchAddress("EndE_tpcAV", EndE_tpcAV, &b_EndE_tpcAV);
   fChain->SetBranchAddress("EndP_tpcAV", EndP_tpcAV, &b_EndP_tpcAV);
   fChain->SetBranchAddress("EndPx_tpcAV", EndPx_tpcAV, &b_EndPx_tpcAV);
   fChain->SetBranchAddress("EndPy_tpcAV", EndPy_tpcAV, &b_EndPy_tpcAV);
   fChain->SetBranchAddress("EndPz_tpcAV", EndPz_tpcAV, &b_EndPz_tpcAV);
   fChain->SetBranchAddress("pathlen_drifted", pathlen_drifted, &b_pathlen_drifted);
   fChain->SetBranchAddress("inTPCDrifted", inTPCDrifted, &b_inTPCDrifted);
   fChain->SetBranchAddress("StartPointx_drifted", StartPointx_drifted, &b_StartPointx_drifted);
   fChain->SetBranchAddress("StartPointy_drifted", StartPointy_drifted, &b_StartPointy_drifted);
   fChain->SetBranchAddress("StartPointz_drifted", StartPointz_drifted, &b_StartPointz_drifted);
   fChain->SetBranchAddress("StartT_drifted", StartT_drifted, &b_StartT_drifted);
   fChain->SetBranchAddress("StartE_drifted", StartE_drifted, &b_StartE_drifted);
   fChain->SetBranchAddress("StartP_drifted", StartP_drifted, &b_StartP_drifted);
   fChain->SetBranchAddress("StartPx_drifted", StartPx_drifted, &b_StartPx_drifted);
   fChain->SetBranchAddress("StartPy_drifted", StartPy_drifted, &b_StartPy_drifted);
   fChain->SetBranchAddress("StartPz_drifted", StartPz_drifted, &b_StartPz_drifted);
   fChain->SetBranchAddress("EndPointx_drifted", EndPointx_drifted, &b_EndPointx_drifted);
   fChain->SetBranchAddress("EndPointy_drifted", EndPointy_drifted, &b_EndPointy_drifted);
   fChain->SetBranchAddress("EndPointz_drifted", EndPointz_drifted, &b_EndPointz_drifted);
   fChain->SetBranchAddress("EndT_drifted", EndT_drifted, &b_EndT_drifted);
   fChain->SetBranchAddress("EndE_drifted", EndE_drifted, &b_EndE_drifted);
   fChain->SetBranchAddress("EndP_drifted", EndP_drifted, &b_EndP_drifted);
   fChain->SetBranchAddress("EndPx_drifted", EndPx_drifted, &b_EndPx_drifted);
   fChain->SetBranchAddress("EndPy_drifted", EndPy_drifted, &b_EndPy_drifted);
   fChain->SetBranchAddress("EndPz_drifted", EndPz_drifted, &b_EndPz_drifted);
   fChain->SetBranchAddress("NumberDaughters", NumberDaughters, &b_NumberDaughters);
   fChain->SetBranchAddress("Mother", Mother, &b_Mother);
   fChain->SetBranchAddress("TrackId", TrackId, &b_TrackId);
   fChain->SetBranchAddress("MergedId", MergedId, &b_MergedId);
   fChain->SetBranchAddress("origin", origin, &b_origin);
   fChain->SetBranchAddress("MCTruthIndex", MCTruthIndex, &b_MCTruthIndex);
   fChain->SetBranchAddress("process_primary", process_primary, &b_process_primary);
   fChain->SetBranchAddress("processname", &processname, &b_processname);
   fChain->SetBranchAddress("nPFParticles", &nPFParticles, &b_nPFParticles);
   fChain->SetBranchAddress("pfp_selfID", pfp_selfID, &b_pfp_selfID);
   fChain->SetBranchAddress("pfp_isPrimary", pfp_isPrimary, &b_pfp_isPrimary);
   fChain->SetBranchAddress("pfp_numDaughters", pfp_numDaughters, &b_pfp_numDaughters);
   fChain->SetBranchAddress("pfp_daughterIDs", pfp_daughterIDs, &b_pfp_daughterIDs);
   fChain->SetBranchAddress("pfp_parentID", pfp_parentID, &b_pfp_parentID);
   fChain->SetBranchAddress("pfp_vertexID", pfp_vertexID, &b_pfp_vertexID);
   fChain->SetBranchAddress("pfp_isShower", pfp_isShower, &b_pfp_isShower);
   fChain->SetBranchAddress("pfp_isTrack", pfp_isTrack, &b_pfp_isTrack);
   fChain->SetBranchAddress("pfp_trackID", pfp_trackID, &b_pfp_trackID);
   fChain->SetBranchAddress("pfp_showerID", pfp_showerID, &b_pfp_showerID);
   fChain->SetBranchAddress("pfp_pdgCode", pfp_pdgCode, &b_pfp_pdgCode);
   fChain->SetBranchAddress("pfp_numClusters", pfp_numClusters, &b_pfp_numClusters);
   fChain->SetBranchAddress("pfp_clusterIDs", pfp_clusterIDs, &b_pfp_clusterIDs);
   fChain->SetBranchAddress("pfp_isNeutrino", pfp_isNeutrino, &b_pfp_isNeutrino);
   fChain->SetBranchAddress("pfp_numNeutrinos", &pfp_numNeutrinos, &b_pfp_numNeutrinos);
   fChain->SetBranchAddress("pfp_neutrinoIDs", pfp_neutrinoIDs, &b_pfp_neutrinoIDs);
   fChain->SetBranchAddress("nshowers_emshower", &nshowers_emshower, &b_nshowers_emshower);
   fChain->SetBranchAddress("showerID_emshower", showerID_emshower, &b_showerID_emshower);
   fChain->SetBranchAddress("shwr_bestplane_emshower", shwr_bestplane_emshower, &b_shwr_bestplane_emshower);
   fChain->SetBranchAddress("shwr_length_emshower", shwr_length_emshower, &b_shwr_length_emshower);
   fChain->SetBranchAddress("shwr_startdcosx_emshower", shwr_startdcosx_emshower, &b_shwr_startdcosx_emshower);
   fChain->SetBranchAddress("shwr_startdcosy_emshower", shwr_startdcosy_emshower, &b_shwr_startdcosy_emshower);
   fChain->SetBranchAddress("shwr_startdcosz_emshower", shwr_startdcosz_emshower, &b_shwr_startdcosz_emshower);
   fChain->SetBranchAddress("shwr_startx_emshower", shwr_startx_emshower, &b_shwr_startx_emshower);
   fChain->SetBranchAddress("shwr_starty_emshower", shwr_starty_emshower, &b_shwr_starty_emshower);
   fChain->SetBranchAddress("shwr_startz_emshower", shwr_startz_emshower, &b_shwr_startz_emshower);
   fChain->SetBranchAddress("shwr_totEng_emshower", shwr_totEng_emshower, &b_shwr_totEng_emshower);
   fChain->SetBranchAddress("shwr_dedx_emshower", shwr_dedx_emshower, &b_shwr_dedx_emshower);
   fChain->SetBranchAddress("shwr_mipEng_emshower", shwr_mipEng_emshower, &b_shwr_mipEng_emshower);
   fChain->SetBranchAddress("shwr_hasPFParticle_emshower", shwr_hasPFParticle_emshower, &b_shwr_hasPFParticle_emshower);
   fChain->SetBranchAddress("shwr_PFParticleID_emshower", shwr_PFParticleID_emshower, &b_shwr_PFParticleID_emshower);
   fChain->SetBranchAddress("shwr_pidmvamu_emshower", shwr_pidmvamu_emshower, &b_shwr_pidmvamu_emshower);
   fChain->SetBranchAddress("shwr_pidmvae_emshower", shwr_pidmvae_emshower, &b_shwr_pidmvae_emshower);
   fChain->SetBranchAddress("shwr_pidmvapich_emshower", shwr_pidmvapich_emshower, &b_shwr_pidmvapich_emshower);
   fChain->SetBranchAddress("shwr_pidmvaphoton_emshower", shwr_pidmvaphoton_emshower, &b_shwr_pidmvaphoton_emshower);
   fChain->SetBranchAddress("shwr_pidmvapr_emshower", shwr_pidmvapr_emshower, &b_shwr_pidmvapr_emshower);
   fChain->SetBranchAddress("ntracks_pandoraTrack", &ntracks_pandoraTrack, &b_ntracks_pandoraTrack);
   fChain->SetBranchAddress("trkId_pandoraTrack", trkId_pandoraTrack, &b_trkId_pandoraTrack);
   fChain->SetBranchAddress("trkncosmictags_tagger_pandoraTrack", trkncosmictags_tagger_pandoraTrack, &b_trkncosmictags_tagger_pandoraTrack);
   fChain->SetBranchAddress("trkcosmicscore_tagger_pandoraTrack", trkcosmicscore_tagger_pandoraTrack, &b_trkcosmicscore_tagger_pandoraTrack);
   fChain->SetBranchAddress("trkcosmictype_tagger_pandoraTrack", trkcosmictype_tagger_pandoraTrack, &b_trkcosmictype_tagger_pandoraTrack);
   fChain->SetBranchAddress("trkncosmictags_containmenttagger_pandoraTrack", trkncosmictags_containmenttagger_pandoraTrack, &b_trkncosmictags_containmenttagger_pandoraTrack);
   fChain->SetBranchAddress("trkcosmicscore_containmenttagger_pandoraTrack", trkcosmicscore_containmenttagger_pandoraTrack, &b_trkcosmicscore_containmenttagger_pandoraTrack);
   fChain->SetBranchAddress("trkcosmictype_containmenttagger_pandoraTrack", trkcosmictype_containmenttagger_pandoraTrack, &b_trkcosmictype_containmenttagger_pandoraTrack);
   fChain->SetBranchAddress("trkncosmictags_flashmatch_pandoraTrack", trkncosmictags_flashmatch_pandoraTrack, &b_trkncosmictags_flashmatch_pandoraTrack);
   fChain->SetBranchAddress("trkcosmicscore_flashmatch_pandoraTrack", trkcosmicscore_flashmatch_pandoraTrack, &b_trkcosmicscore_flashmatch_pandoraTrack);
   fChain->SetBranchAddress("trkcosmictype_flashmatch_pandoraTrack", trkcosmictype_flashmatch_pandoraTrack, &b_trkcosmictype_flashmatch_pandoraTrack);
   fChain->SetBranchAddress("trkke_pandoraTrack", trkke_pandoraTrack, &b_trkke_pandoraTrack);
   fChain->SetBranchAddress("trkrange_pandoraTrack", trkrange_pandoraTrack, &b_trkrange_pandoraTrack);
   fChain->SetBranchAddress("trkidtruth_pandoraTrack", trkidtruth_pandoraTrack, &b_trkidtruth_pandoraTrack);
   fChain->SetBranchAddress("trkorigin_pandoraTrack", trkorigin_pandoraTrack, &b_trkorigin_pandoraTrack);
   fChain->SetBranchAddress("trkpdgtruth_pandoraTrack", trkpdgtruth_pandoraTrack, &b_trkpdgtruth_pandoraTrack);
   fChain->SetBranchAddress("trkefftruth_pandoraTrack", trkefftruth_pandoraTrack, &b_trkefftruth_pandoraTrack);
   fChain->SetBranchAddress("trkpurtruth_pandoraTrack", trkpurtruth_pandoraTrack, &b_trkpurtruth_pandoraTrack);
   fChain->SetBranchAddress("trkpitchc_pandoraTrack", trkpitchc_pandoraTrack, &b_trkpitchc_pandoraTrack);
   fChain->SetBranchAddress("ntrkhits_pandoraTrack", ntrkhits_pandoraTrack, &b_ntrkhits_pandoraTrack);
   fChain->SetBranchAddress("trkdedx_pandoraTrack", trkdedx_pandoraTrack, &b_trkdedx_pandoraTrack);
   fChain->SetBranchAddress("trkdqdx_pandoraTrack", trkdqdx_pandoraTrack, &b_trkdqdx_pandoraTrack);
   fChain->SetBranchAddress("trkresrg_pandoraTrack", trkresrg_pandoraTrack, &b_trkresrg_pandoraTrack);
   fChain->SetBranchAddress("trktpc_pandoraTrack", trktpc_pandoraTrack, &b_trktpc_pandoraTrack);
   fChain->SetBranchAddress("trkxyz_pandoraTrack", trkxyz_pandoraTrack, &b_trkxyz_pandoraTrack);
   fChain->SetBranchAddress("trkstartx_pandoraTrack", trkstartx_pandoraTrack, &b_trkstartx_pandoraTrack);
   fChain->SetBranchAddress("trkstarty_pandoraTrack", trkstarty_pandoraTrack, &b_trkstarty_pandoraTrack);
   fChain->SetBranchAddress("trkstartz_pandoraTrack", trkstartz_pandoraTrack, &b_trkstartz_pandoraTrack);
   fChain->SetBranchAddress("trkstartd_pandoraTrack", trkstartd_pandoraTrack, &b_trkstartd_pandoraTrack);
   fChain->SetBranchAddress("trkendx_pandoraTrack", trkendx_pandoraTrack, &b_trkendx_pandoraTrack);
   fChain->SetBranchAddress("trkendy_pandoraTrack", trkendy_pandoraTrack, &b_trkendy_pandoraTrack);
   fChain->SetBranchAddress("trkendz_pandoraTrack", trkendz_pandoraTrack, &b_trkendz_pandoraTrack);
   fChain->SetBranchAddress("trkendd_pandoraTrack", trkendd_pandoraTrack, &b_trkendd_pandoraTrack);
   fChain->SetBranchAddress("trkflashT0_pandoraTrack", trkflashT0_pandoraTrack, &b_trkflashT0_pandoraTrack);
   fChain->SetBranchAddress("trktrueT0_pandoraTrack", trktrueT0_pandoraTrack, &b_trktrueT0_pandoraTrack);
   fChain->SetBranchAddress("trkg4id_pandoraTrack", trkg4id_pandoraTrack, &b_trkg4id_pandoraTrack);
   fChain->SetBranchAddress("trkorig_pandoraTrack", trkorig_pandoraTrack, &b_trkorig_pandoraTrack);
   fChain->SetBranchAddress("trkpurity_pandoraTrack", trkpurity_pandoraTrack, &b_trkpurity_pandoraTrack);
   fChain->SetBranchAddress("trkcompleteness_pandoraTrack", trkcompleteness_pandoraTrack, &b_trkcompleteness_pandoraTrack);
   fChain->SetBranchAddress("trktheta_pandoraTrack", trktheta_pandoraTrack, &b_trktheta_pandoraTrack);
   fChain->SetBranchAddress("trkphi_pandoraTrack", trkphi_pandoraTrack, &b_trkphi_pandoraTrack);
   fChain->SetBranchAddress("trkstartdcosx_pandoraTrack", trkstartdcosx_pandoraTrack, &b_trkstartdcosx_pandoraTrack);
   fChain->SetBranchAddress("trkstartdcosy_pandoraTrack", trkstartdcosy_pandoraTrack, &b_trkstartdcosy_pandoraTrack);
   fChain->SetBranchAddress("trkstartdcosz_pandoraTrack", trkstartdcosz_pandoraTrack, &b_trkstartdcosz_pandoraTrack);
   fChain->SetBranchAddress("trkenddcosx_pandoraTrack", trkenddcosx_pandoraTrack, &b_trkenddcosx_pandoraTrack);
   fChain->SetBranchAddress("trkenddcosy_pandoraTrack", trkenddcosy_pandoraTrack, &b_trkenddcosy_pandoraTrack);
   fChain->SetBranchAddress("trkenddcosz_pandoraTrack", trkenddcosz_pandoraTrack, &b_trkenddcosz_pandoraTrack);
   fChain->SetBranchAddress("trkthetaxz_pandoraTrack", trkthetaxz_pandoraTrack, &b_trkthetaxz_pandoraTrack);
   fChain->SetBranchAddress("trkthetayz_pandoraTrack", trkthetayz_pandoraTrack, &b_trkthetayz_pandoraTrack);
   fChain->SetBranchAddress("trkmom_pandoraTrack", trkmom_pandoraTrack, &b_trkmom_pandoraTrack);
   fChain->SetBranchAddress("trkmomrange_pandoraTrack", trkmomrange_pandoraTrack, &b_trkmomrange_pandoraTrack);
   fChain->SetBranchAddress("trkmommschi2_pandoraTrack", trkmommschi2_pandoraTrack, &b_trkmommschi2_pandoraTrack);
   fChain->SetBranchAddress("trkmommsllhd_pandoraTrack", trkmommsllhd_pandoraTrack, &b_trkmommsllhd_pandoraTrack);
   fChain->SetBranchAddress("trklen_pandoraTrack", trklen_pandoraTrack, &b_trklen_pandoraTrack);
   fChain->SetBranchAddress("trksvtxid_pandoraTrack", trksvtxid_pandoraTrack, &b_trksvtxid_pandoraTrack);
   fChain->SetBranchAddress("trkevtxid_pandoraTrack", trkevtxid_pandoraTrack, &b_trkevtxid_pandoraTrack);
   fChain->SetBranchAddress("trkpidmvamu_pandoraTrack", trkpidmvamu_pandoraTrack, &b_trkpidmvamu_pandoraTrack);
   fChain->SetBranchAddress("trkpidmvae_pandoraTrack", trkpidmvae_pandoraTrack, &b_trkpidmvae_pandoraTrack);
   fChain->SetBranchAddress("trkpidmvapich_pandoraTrack", trkpidmvapich_pandoraTrack, &b_trkpidmvapich_pandoraTrack);
   fChain->SetBranchAddress("trkpidmvaphoton_pandoraTrack", trkpidmvaphoton_pandoraTrack, &b_trkpidmvaphoton_pandoraTrack);
   fChain->SetBranchAddress("trkpidmvapr_pandoraTrack", trkpidmvapr_pandoraTrack, &b_trkpidmvapr_pandoraTrack);
   fChain->SetBranchAddress("trkpidpdg_pandoraTrack", trkpidpdg_pandoraTrack, &b_trkpidpdg_pandoraTrack);
   fChain->SetBranchAddress("trkpidchi_pandoraTrack", trkpidchi_pandoraTrack, &b_trkpidchi_pandoraTrack);
   fChain->SetBranchAddress("trkpidchipr_pandoraTrack", trkpidchipr_pandoraTrack, &b_trkpidchipr_pandoraTrack);
   fChain->SetBranchAddress("trkpidchika_pandoraTrack", trkpidchika_pandoraTrack, &b_trkpidchika_pandoraTrack);
   fChain->SetBranchAddress("trkpidchipi_pandoraTrack", trkpidchipi_pandoraTrack, &b_trkpidchipi_pandoraTrack);
   fChain->SetBranchAddress("trkpidchimu_pandoraTrack", trkpidchimu_pandoraTrack, &b_trkpidchimu_pandoraTrack);
   fChain->SetBranchAddress("trkpidpida_pandoraTrack", trkpidpida_pandoraTrack, &b_trkpidpida_pandoraTrack);
   fChain->SetBranchAddress("trkpidbestplane_pandoraTrack", trkpidbestplane_pandoraTrack, &b_trkpidbestplane_pandoraTrack);
   fChain->SetBranchAddress("trkhasPFParticle_pandoraTrack", trkhasPFParticle_pandoraTrack, &b_trkhasPFParticle_pandoraTrack);
   fChain->SetBranchAddress("trkPFParticleID_pandoraTrack", trkPFParticleID_pandoraTrack, &b_trkPFParticleID_pandoraTrack);
   fChain->SetBranchAddress("nvtx_linecluster", &nvtx_linecluster, &b_nvtx_linecluster);
   fChain->SetBranchAddress("vtxId_linecluster", vtxId_linecluster, &b_vtxId_linecluster);
   fChain->SetBranchAddress("vtxx_linecluster", vtxx_linecluster, &b_vtxx_linecluster);
   fChain->SetBranchAddress("vtxy_linecluster", vtxy_linecluster, &b_vtxy_linecluster);
   fChain->SetBranchAddress("vtxz_linecluster", vtxz_linecluster, &b_vtxz_linecluster);
   fChain->SetBranchAddress("vtxhasPFParticle_linecluster", vtxhasPFParticle_linecluster, &b_vtxhasPFParticle_linecluster);
   fChain->SetBranchAddress("vtxPFParticleID_linecluster", vtxPFParticleID_linecluster, &b_vtxPFParticleID_linecluster);
   fChain->SetBranchAddress("nvtx_pandora", &nvtx_pandora, &b_nvtx_pandora);
   fChain->SetBranchAddress("vtxId_pandora", vtxId_pandora, &b_vtxId_pandora);
   fChain->SetBranchAddress("vtxx_pandora", vtxx_pandora, &b_vtxx_pandora);
   fChain->SetBranchAddress("vtxy_pandora", vtxy_pandora, &b_vtxy_pandora);
   fChain->SetBranchAddress("vtxz_pandora", vtxz_pandora, &b_vtxz_pandora);
   fChain->SetBranchAddress("vtxhasPFParticle_pandora", vtxhasPFParticle_pandora, &b_vtxhasPFParticle_pandora);
   fChain->SetBranchAddress("vtxPFParticleID_pandora", vtxPFParticleID_pandora, &b_vtxPFParticleID_pandora);
   fChain->SetBranchAddress("no_mcshowers", &no_mcshowers, &b_no_mcshowers);
   fChain->SetBranchAddress("mcshwr_origin", mcshwr_origin, &b_mcshwr_origin);
   fChain->SetBranchAddress("mcshwr_pdg", mcshwr_pdg, &b_mcshwr_pdg);
   fChain->SetBranchAddress("mcshwr_TrackId", mcshwr_TrackId, &b_mcshwr_TrackId);
   fChain->SetBranchAddress("mcshwr_Process", &mcshwr_Process, &b_mcshwr_Process);
   fChain->SetBranchAddress("mcshwr_startX", mcshwr_startX, &b_mcshwr_startX);
   fChain->SetBranchAddress("mcshwr_startY", mcshwr_startY, &b_mcshwr_startY);
   fChain->SetBranchAddress("mcshwr_startZ", mcshwr_startZ, &b_mcshwr_startZ);
   fChain->SetBranchAddress("mcshwr_endX", mcshwr_endX, &b_mcshwr_endX);
   fChain->SetBranchAddress("mcshwr_endY", mcshwr_endY, &b_mcshwr_endY);
   fChain->SetBranchAddress("mcshwr_endZ", mcshwr_endZ, &b_mcshwr_endZ);
   fChain->SetBranchAddress("mcshwr_CombEngX", mcshwr_CombEngX, &b_mcshwr_CombEngX);
   fChain->SetBranchAddress("mcshwr_CombEngY", mcshwr_CombEngY, &b_mcshwr_CombEngY);
   fChain->SetBranchAddress("mcshwr_CombEngZ", mcshwr_CombEngZ, &b_mcshwr_CombEngZ);
   fChain->SetBranchAddress("mcshwr_CombEngPx", mcshwr_CombEngPx, &b_mcshwr_CombEngPx);
   fChain->SetBranchAddress("mcshwr_CombEngPy", mcshwr_CombEngPy, &b_mcshwr_CombEngPy);
   fChain->SetBranchAddress("mcshwr_CombEngPz", mcshwr_CombEngPz, &b_mcshwr_CombEngPz);
   fChain->SetBranchAddress("mcshwr_CombEngE", mcshwr_CombEngE, &b_mcshwr_CombEngE);
   fChain->SetBranchAddress("mcshwr_dEdx", mcshwr_dEdx, &b_mcshwr_dEdx);
   fChain->SetBranchAddress("mcshwr_dQdx", mcshwr_dQdx, &b_mcshwr_dQdx);
   fChain->SetBranchAddress("mcshwr_StartDirX", mcshwr_StartDirX, &b_mcshwr_StartDirX);
   fChain->SetBranchAddress("mcshwr_StartDirY", mcshwr_StartDirY, &b_mcshwr_StartDirY);
   fChain->SetBranchAddress("mcshwr_StartDirZ", mcshwr_StartDirZ, &b_mcshwr_StartDirZ);
   fChain->SetBranchAddress("mcshwr_isEngDeposited", mcshwr_isEngDeposited, &b_mcshwr_isEngDeposited);
   fChain->SetBranchAddress("mcshwr_Motherpdg", mcshwr_Motherpdg, &b_mcshwr_Motherpdg);
   fChain->SetBranchAddress("mcshwr_MotherTrkId", mcshwr_MotherTrkId, &b_mcshwr_MotherTrkId);
   fChain->SetBranchAddress("mcshwr_MotherProcess", &mcshwr_MotherProcess, &b_mcshwr_MotherProcess);
   fChain->SetBranchAddress("mcshwr_MotherstartX", mcshwr_MotherstartX, &b_mcshwr_MotherstartX);
   fChain->SetBranchAddress("mcshwr_MotherstartY", mcshwr_MotherstartY, &b_mcshwr_MotherstartY);
   fChain->SetBranchAddress("mcshwr_MotherstartZ", mcshwr_MotherstartZ, &b_mcshwr_MotherstartZ);
   fChain->SetBranchAddress("mcshwr_MotherendX", mcshwr_MotherendX, &b_mcshwr_MotherendX);
   fChain->SetBranchAddress("mcshwr_MotherendY", mcshwr_MotherendY, &b_mcshwr_MotherendY);
   fChain->SetBranchAddress("mcshwr_MotherendZ", mcshwr_MotherendZ, &b_mcshwr_MotherendZ);
   fChain->SetBranchAddress("mcshwr_Ancestorpdg", mcshwr_Ancestorpdg, &b_mcshwr_Ancestorpdg);
   fChain->SetBranchAddress("mcshwr_AncesotorTrkId", mcshwr_AncesotorTrkId, &b_mcshwr_AncesotorTrkId);
   fChain->SetBranchAddress("mcshwr_AncesotorProcess", &mcshwr_AncesotorProcess, &b_mcshwr_AncesotorProcess);
   fChain->SetBranchAddress("mcshwr_AncestorstartX", mcshwr_AncestorstartX, &b_mcshwr_AncestorstartX);
   fChain->SetBranchAddress("mcshwr_AncestorstartY", mcshwr_AncestorstartY, &b_mcshwr_AncestorstartY);
   fChain->SetBranchAddress("mcshwr_AncestorstartZ", mcshwr_AncestorstartZ, &b_mcshwr_AncestorstartZ);
   fChain->SetBranchAddress("mcshwr_AncestorendX", mcshwr_AncestorendX, &b_mcshwr_AncestorendX);
   fChain->SetBranchAddress("mcshwr_AncestorendY", mcshwr_AncestorendY, &b_mcshwr_AncestorendY);
   fChain->SetBranchAddress("mcshwr_AncestorendZ", mcshwr_AncestorendZ, &b_mcshwr_AncestorendZ);
   fChain->SetBranchAddress("no_mctracks", &no_mctracks, &b_no_mctracks);
   fChain->SetBranchAddress("mctrk_origin", mctrk_origin, &b_mctrk_origin);
   fChain->SetBranchAddress("mctrk_pdg", mctrk_pdg, &b_mctrk_pdg);
   fChain->SetBranchAddress("mctrk_TrackId", mctrk_TrackId, &b_mctrk_TrackId);
   fChain->SetBranchAddress("mctrk_Process", &mctrk_Process, &b_mctrk_Process);
   fChain->SetBranchAddress("mctrk_startX", mctrk_startX, &b_mctrk_startX);
   fChain->SetBranchAddress("mctrk_startY", mctrk_startY, &b_mctrk_startY);
   fChain->SetBranchAddress("mctrk_startZ", mctrk_startZ, &b_mctrk_startZ);
   fChain->SetBranchAddress("mctrk_endX", mctrk_endX, &b_mctrk_endX);
   fChain->SetBranchAddress("mctrk_endY", mctrk_endY, &b_mctrk_endY);
   fChain->SetBranchAddress("mctrk_endZ", mctrk_endZ, &b_mctrk_endZ);
   fChain->SetBranchAddress("mctrk_startX_drifted", mctrk_startX_drifted, &b_mctrk_startX_drifted);
   fChain->SetBranchAddress("mctrk_startY_drifted", mctrk_startY_drifted, &b_mctrk_startY_drifted);
   fChain->SetBranchAddress("mctrk_startZ_drifted", mctrk_startZ_drifted, &b_mctrk_startZ_drifted);
   fChain->SetBranchAddress("mctrk_endX_drifted", mctrk_endX_drifted, &b_mctrk_endX_drifted);
   fChain->SetBranchAddress("mctrk_endY_drifted", mctrk_endY_drifted, &b_mctrk_endY_drifted);
   fChain->SetBranchAddress("mctrk_endZ_drifted", mctrk_endZ_drifted, &b_mctrk_endZ_drifted);
   fChain->SetBranchAddress("mctrk_len_drifted", mctrk_len_drifted, &b_mctrk_len_drifted);
   fChain->SetBranchAddress("mctrk_p_drifted", mctrk_p_drifted, &b_mctrk_p_drifted);
   fChain->SetBranchAddress("mctrk_px_drifted", mctrk_px_drifted, &b_mctrk_px_drifted);
   fChain->SetBranchAddress("mctrk_py_drifted", mctrk_py_drifted, &b_mctrk_py_drifted);
   fChain->SetBranchAddress("mctrk_pz_drifted", mctrk_pz_drifted, &b_mctrk_pz_drifted);
   fChain->SetBranchAddress("nmctrksteps", nmctrksteps, &b_nmctrksteps);
   fChain->SetBranchAddress("mctrk_dEdx", mctrk_dEdx, &b_mctrk_dEdx);
   fChain->SetBranchAddress("mctrk_dQdx", mctrk_dQdx, &b_mctrk_dQdx);
   fChain->SetBranchAddress("mctrk_Motherpdg", mctrk_Motherpdg, &b_mctrk_Motherpdg);
   fChain->SetBranchAddress("mctrk_MotherTrkId", mctrk_MotherTrkId, &b_mctrk_MotherTrkId);
   fChain->SetBranchAddress("mctrk_MotherProcess", &mctrk_MotherProcess, &b_mctrk_MotherProcess);
   fChain->SetBranchAddress("mctrk_MotherstartX", mctrk_MotherstartX, &b_mctrk_MotherstartX);
   fChain->SetBranchAddress("mctrk_MotherstartY", mctrk_MotherstartY, &b_mctrk_MotherstartY);
   fChain->SetBranchAddress("mctrk_MotherstartZ", mctrk_MotherstartZ, &b_mctrk_MotherstartZ);
   fChain->SetBranchAddress("mctrk_MotherendX", mctrk_MotherendX, &b_mctrk_MotherendX);
   fChain->SetBranchAddress("mctrk_MotherendY", mctrk_MotherendY, &b_mctrk_MotherendY);
   fChain->SetBranchAddress("mctrk_MotherendZ", mctrk_MotherendZ, &b_mctrk_MotherendZ);
   fChain->SetBranchAddress("mctrk_Ancestorpdg", mctrk_Ancestorpdg, &b_mctrk_Ancestorpdg);
   fChain->SetBranchAddress("mctrk_AncesotorTrkId", mctrk_AncesotorTrkId, &b_mctrk_AncesotorTrkId);
   fChain->SetBranchAddress("mctrk_AncesotorProcess", &mctrk_AncesotorProcess, &b_mctrk_AncesotorProcess);
   fChain->SetBranchAddress("mctrk_AncestorstartX", mctrk_AncestorstartX, &b_mctrk_AncestorstartX);
   fChain->SetBranchAddress("mctrk_AncestorstartY", mctrk_AncestorstartY, &b_mctrk_AncestorstartY);
   fChain->SetBranchAddress("mctrk_AncestorstartZ", mctrk_AncestorstartZ, &b_mctrk_AncestorstartZ);
   fChain->SetBranchAddress("mctrk_AncestorendX", mctrk_AncestorendX, &b_mctrk_AncestorendX);
   fChain->SetBranchAddress("mctrk_AncestorendY", mctrk_AncestorendY, &b_mctrk_AncestorendY);
   fChain->SetBranchAddress("mctrk_AncestorendZ", mctrk_AncestorendZ, &b_mctrk_AncestorendZ);
   Notify();
}

Bool_t anatree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void anatree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t anatree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef anatree_cxx
