//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jan 14 09:39:22 2026 by ROOT version 6.28/12
// from TTree pandoraOutput/Pandora Output Tree
// found on file: afm_ana_hist.root
//////////////////////////////////////////////////////////

#ifndef afmanatree_h
#define afmanatree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class afmanatree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          eventID;
   UInt_t          runID;
   UInt_t          subrunID;
   UInt_t          nMCParticles;
   UInt_t          nPFParticles;
   Bool_t          mcIsMCPrimary[1390];   //[nMCParticles]
   Int_t           mcParticlePdgCode[1390];   //[nMCParticles]
   Double_t        mcParticleTrueEnergy[1390];   //[nMCParticles]
   Int_t           mcParticleTrackID[1390];   //[nMCParticles]
   Int_t           mcParticleParentTrackID[1390];   //[nMCParticles]
   Char_t          mcParticleStartProcess[1390];   //[nMCParticles]
   Char_t          mcParticleEndProcess[1390];   //[nMCParticles]
   Int_t           mcParticleNTrajectoryPoints[1390];   //[nMCParticles]
   Double_t        mcParticleStartPositionX[1390];   //[nMCParticles]
   Double_t        mcParticleStartPositionY[1390];   //[nMCParticles]
   Double_t        mcParticleStartPositionZ[1390];   //[nMCParticles]
   Double_t        mcParticleStartPositionT[1390];   //[nMCParticles]
   Double_t        mcParticleStartMomentumX[1390];   //[nMCParticles]
   Double_t        mcParticleStartMomentumY[1390];   //[nMCParticles]
   Double_t        mcParticleStartMomentumZ[1390];   //[nMCParticles]
   Double_t        mcParticleStartMomentumE[1390];   //[nMCParticles]
   Double_t        mcParticleEndPositionX[1390];   //[nMCParticles]
   Double_t        mcParticleEndPositionY[1390];   //[nMCParticles]
   Double_t        mcParticleEndPositionZ[1390];   //[nMCParticles]
   Double_t        mcParticleEndPositionT[1390];   //[nMCParticles]
   Double_t        mcParticleEndMomentumX[1390];   //[nMCParticles]
   Double_t        mcParticleEndMomentumY[1390];   //[nMCParticles]
   Double_t        mcParticleEndMomentumZ[1390];   //[nMCParticles]
   Double_t        mcParticleEndMomentumE[1390];   //[nMCParticles]
   Double_t        mcParticleVertexTime[1390];   //[nMCParticles]
   Double_t        mcParticleEndTime[1390];   //[nMCParticles]
   Int_t           mcParticleNHits[1390];   //[nMCParticles]
   Int_t           mcParticleNHitsView[1390][3];   //[nMCParticles]
   Int_t           pfpTrueParticleMatchedID[175];   //[nPFParticles]
   Int_t           pfpTrueParticleMatchedPosition[175];   //[nPFParticles]
   Bool_t          pfpIsPrimary[175];   //[nPFParticles]
   Int_t           pfpID[175];   //[nPFParticles]
   Int_t           pfpParentID[175];   //[nPFParticles]
   Int_t           pfpPdgCode[175];   //[nPFParticles]
   Int_t           pfpNChildren[175];   //[nPFParticles]
   Int_t           pfpNClusters[175];   //[nPFParticles]
   Int_t           pfpNHits[175];   //[nPFParticles]
   Int_t           pfpNHitsView[175][3];   //[nPFParticles]
   Int_t           pfpNSharedTrueParticleHits[175];   //[nPFParticles]
   Int_t           pfpNSharedTrueParticleHitsView[175][3];   //[nPFParticles]
   Int_t           pfpTrueParticleMatchedIDView[175][3];   //[nPFParticles]
   Int_t           pfpTrueParticleMatchedPositionView[175][3];   //[nPFParticles]
   Bool_t          pfpIsTrack[175];   //[nPFParticles]
   Bool_t          pfpIsShower[175];   //[nPFParticles]
   Int_t           pfpTrackID[175];   //[nPFParticles]
   Double_t        pfpTrackLength[175];   //[nPFParticles]
   Double_t        pfpTrackStartX[175];   //[nPFParticles]
   Double_t        pfpTrackStartY[175];   //[nPFParticles]
   Double_t        pfpTrackStartZ[175];   //[nPFParticles]
   Double_t        pfpTrackVertexX[175];   //[nPFParticles]
   Double_t        pfpTrackVertexY[175];   //[nPFParticles]
   Double_t        pfpTrackVertexZ[175];   //[nPFParticles]
   Double_t        pfpTrackEndX[175];   //[nPFParticles]
   Double_t        pfpTrackEndY[175];   //[nPFParticles]
   Double_t        pfpTrackEndZ[175];   //[nPFParticles]
   Double_t        pfpTrackTheta[175];   //[nPFParticles]
   Double_t        pfpTrackPhi[175];   //[nPFParticles]
   Double_t        pfpTrackZenithAngle[175];   //[nPFParticles]
   Double_t        pfpTrackAzimuthAngle[175];   //[nPFParticles]
   Double_t        pfpTrackStartDirectionX[175];   //[nPFParticles]
   Double_t        pfpTrackStartDirectionY[175];   //[nPFParticles]
   Double_t        pfpTrackStartDirectionZ[175];   //[nPFParticles]
   Double_t        pfpTrackVertexDirectionX[175];   //[nPFParticles]
   Double_t        pfpTrackVertexDirectionY[175];   //[nPFParticles]
   Double_t        pfpTrackVertexDirectionZ[175];   //[nPFParticles]
   Double_t        pfpTrackEndDirectionX[175];   //[nPFParticles]
   Double_t        pfpTrackEndDirectionY[175];   //[nPFParticles]
   Double_t        pfpTrackEndDirectionZ[175];   //[nPFParticles]
   Float_t         pfpTrackChi2[175];   //[nPFParticles]
   Int_t           pfpTrackStartNdof[175];   //[nPFParticles]
   Int_t           pfpCluPlane[175][100];   //[nPFParticles]
   Int_t           pfpCluView[175][100];   //[nPFParticles]
   Int_t           pfpCluNHits[175][100];   //[nPFParticles]
   Double_t        pfpCluIntegral[175][100];   //[nPFParticles]
   Int_t           pfpShowerID[175];   //[nPFParticles]
   Int_t           pfpShowerBestPlane[175];   //[nPFParticles]
   Double_t        pfpShowerEnergy[175];   //[nPFParticles]
   Double_t        pfpShowerEnergyToTrueEnergyRatio[175];   //[nPFParticles]
   Double_t        pfpShowerEnergyToTrueMomentumRatio[175];   //[nPFParticles]
   Double_t        pfpShowerDirectionX[175];   //[nPFParticles]
   Double_t        pfpShowerDirectionY[175];   //[nPFParticles]
   Double_t        pfpShowerDirectionZ[175];   //[nPFParticles]
   Double_t        pfpShowerDirectionErrX[175];   //[nPFParticles]
   Double_t        pfpShowerDirectionErrY[175];   //[nPFParticles]
   Double_t        pfpShowerDirectionErrZ[175];   //[nPFParticles]
   Double_t        pfpShowerStartX[175];   //[nPFParticles]
   Double_t        pfpShowerStartY[175];   //[nPFParticles]
   Double_t        pfpShowerStartZ[175];   //[nPFParticles]
   Double_t        pfpShowerStartErrX[175];   //[nPFParticles]
   Double_t        pfpShowerStartErrY[175];   //[nPFParticles]
   Double_t        pfpShowerStartErrZ[175];   //[nPFParticles]
   Double_t        pfpShowerLength[175];   //[nPFParticles]
   Double_t        pfpShowerOpeningAngle[175];   //[nPFParticles]
   Double_t        pfpCompleteness[175];   //[nPFParticles]
   Double_t        pfpCompletenessView[175][3];   //[nPFParticles]
   Double_t        pfpPurity[175];   //[nPFParticles]
   Double_t        pfpPurityView[175][3];   //[nPFParticles]

   // List of branches
   TBranch        *b_eventID;   //!
   TBranch        *b_runID;   //!
   TBranch        *b_subrunID;   //!
   TBranch        *b_nMCParticles;   //!
   TBranch        *b_nPFParticles;   //!
   TBranch        *b_mcIsMCPrimary;   //!
   TBranch        *b_mcParticlePdgCode;   //!
   TBranch        *b_mcParticleTrueEnergy;   //!
   TBranch        *b_mcParticleTrackID;   //!
   TBranch        *b_mcParticleParentTrackID;   //!
   TBranch        *b_mcParticleStartProcess;   //!
   TBranch        *b_mcParticleEndProcess;   //!
   TBranch        *b_mcParticleNTrajectoryPoints;   //!
   TBranch        *b_mcParticleStartPositionX;   //!
   TBranch        *b_mcParticleStartPositionY;   //!
   TBranch        *b_mcParticleStartPositionZ;   //!
   TBranch        *b_mcParticleStartPositionT;   //!
   TBranch        *b_mcParticleStartMomentumX;   //!
   TBranch        *b_mcParticleStartMomentumY;   //!
   TBranch        *b_mcParticleStartMomentumZ;   //!
   TBranch        *b_mcParticleStartMomentumE;   //!
   TBranch        *b_mcParticleEndPositionX;   //!
   TBranch        *b_mcParticleEndPositionY;   //!
   TBranch        *b_mcParticleEndPositionZ;   //!
   TBranch        *b_mcParticleEndPositionT;   //!
   TBranch        *b_mcParticleEndMomentumX;   //!
   TBranch        *b_mcParticleEndMomentumY;   //!
   TBranch        *b_mcParticleEndMomentumZ;   //!
   TBranch        *b_mcParticleEndMomentumE;   //!
   TBranch        *b_mcParticleVertexTime;   //!
   TBranch        *b_mcParticleEndTime;   //!
   TBranch        *b_mcParticleNHits;   //!
   TBranch        *b_mcParticleNHitsView;   //!
   TBranch        *b_pfpTrueParticleMatchedID;   //!
   TBranch        *b_pfpTrueParticleMatchedPosition;   //!
   TBranch        *b_pfpIsPrimary;   //!
   TBranch        *b_pfpID;   //!
   TBranch        *b_pfpParentID;   //!
   TBranch        *b_pfpPdgCode;   //!
   TBranch        *b_pfpNChildren;   //!
   TBranch        *b_pfpNClusters;   //!
   TBranch        *b_pfpNHits;   //!
   TBranch        *b_pfpNHitsView;   //!
   TBranch        *b_pfpNSharedTrueParticleHits;   //!
   TBranch        *b_pfpNSharedTrueParticleHitsView;   //!
   TBranch        *b_pfpTrueParticleMatchedIDView;   //!
   TBranch        *b_pfpTrueParticleMatchedPositionView;   //!
   TBranch        *b_pfpIsTrack;   //!
   TBranch        *b_pfpIsShower;   //!
   TBranch        *b_pfpTrackID;   //!
   TBranch        *b_pfpTrackLength;   //!
   TBranch        *b_pfpTrackStartX;   //!
   TBranch        *b_pfpTrackStartY;   //!
   TBranch        *b_pfpTrackStartZ;   //!
   TBranch        *b_pfpTrackVertexX;   //!
   TBranch        *b_pfpTrackVertexY;   //!
   TBranch        *b_pfpTrackVertexZ;   //!
   TBranch        *b_pfpTrackEndX;   //!
   TBranch        *b_pfpTrackEndY;   //!
   TBranch        *b_pfpTrackEndZ;   //!
   TBranch        *b_pfpTrackTheta;   //!
   TBranch        *b_pfpTrackPhi;   //!
   TBranch        *b_pfpTrackZenithAngle;   //!
   TBranch        *b_pfpTrackAzimuthAngle;   //!
   TBranch        *b_pfpTrackStartDirectionX;   //!
   TBranch        *b_pfpTrackStartDirectionY;   //!
   TBranch        *b_pfpTrackStartDirectionZ;   //!
   TBranch        *b_pfpTrackVertexDirectionX;   //!
   TBranch        *b_pfpTrackVertexDirectionY;   //!
   TBranch        *b_pfpTrackVertexDirectionZ;   //!
   TBranch        *b_pfpTrackEndDirectionX;   //!
   TBranch        *b_pfpTrackEndDirectionY;   //!
   TBranch        *b_pfpTrackEndDirectionZ;   //!
   TBranch        *b_pfpTrackChi2;   //!
   TBranch        *b_pfpTrackStartNdof;   //!
   TBranch        *b_pfpCluPlane;   //!
   TBranch        *b_pfpCluView;   //!
   TBranch        *b_pfpCluNHits;   //!
   TBranch        *b_pfpCluIntegral;   //!
   TBranch        *b_pfpShowerID;   //!
   TBranch        *b_pfpShowerBestPlane;   //!
   TBranch        *b_pfpShowerEnergy;   //!
   TBranch        *b_pfpShowerEnergyToTrueEnergyRatio;   //!
   TBranch        *b_pfpShowerEnergyToTrueMomentumRatio;   //!
   TBranch        *b_pfpShowerDirectionX;   //!
   TBranch        *b_pfpShowerDirectionY;   //!
   TBranch        *b_pfpShowerDirectionZ;   //!
   TBranch        *b_pfpShowerDirectionErrX;   //!
   TBranch        *b_pfpShowerDirectionErrY;   //!
   TBranch        *b_pfpShowerDirectionErrZ;   //!
   TBranch        *b_pfpShowerStartX;   //!
   TBranch        *b_pfpShowerStartY;   //!
   TBranch        *b_pfpShowerStartZ;   //!
   TBranch        *b_pfpShowerStartErrX;   //!
   TBranch        *b_pfpShowerStartErrY;   //!
   TBranch        *b_pfpShowerStartErrZ;   //!
   TBranch        *b_pfpShowerLength;   //!
   TBranch        *b_pfpShowerOpeningAngle;   //!
   TBranch        *b_pfpCompleteness;   //!
   TBranch        *b_pfpCompletenessView;   //!
   TBranch        *b_pfpPurity;   //!
   TBranch        *b_pfpPurityView;   //!

   afmanatree(TTree *tree=0);
   virtual ~afmanatree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef afmanatree_cxx
afmanatree::afmanatree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("afm_ana_hist.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("afm_ana_hist.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("afm_ana_hist.root:/ana");
      dir->GetObject("pandoraOutput",tree);

   }
   Init(tree);
}

afmanatree::~afmanatree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t afmanatree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t afmanatree::LoadTree(Long64_t entry)
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

void afmanatree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("eventID", &eventID, &b_eventID);
   fChain->SetBranchAddress("runID", &runID, &b_runID);
   fChain->SetBranchAddress("subrunID", &subrunID, &b_subrunID);
   fChain->SetBranchAddress("nMCParticles", &nMCParticles, &b_nMCParticles);
   fChain->SetBranchAddress("nPFParticles", &nPFParticles, &b_nPFParticles);
   fChain->SetBranchAddress("mcIsMCPrimary", mcIsMCPrimary, &b_mcIsMCPrimary);
   fChain->SetBranchAddress("mcParticlePdgCode", mcParticlePdgCode, &b_mcParticlePdgCode);
   fChain->SetBranchAddress("mcParticleTrueEnergy", mcParticleTrueEnergy, &b_mcParticleTrueEnergy);
   fChain->SetBranchAddress("mcParticleTrackID", mcParticleTrackID, &b_mcParticleTrackID);
   fChain->SetBranchAddress("mcParticleParentTrackID", mcParticleParentTrackID, &b_mcParticleParentTrackID);
   fChain->SetBranchAddress("mcParticleStartProcess", mcParticleStartProcess, &b_mcParticleStartProcess);
   fChain->SetBranchAddress("mcParticleEndProcess", mcParticleEndProcess, &b_mcParticleEndProcess);
   fChain->SetBranchAddress("mcParticleNTrajectoryPoints", mcParticleNTrajectoryPoints, &b_mcParticleNTrajectoryPoints);
   fChain->SetBranchAddress("mcParticleStartPositionX", mcParticleStartPositionX, &b_mcParticleStartPositionX);
   fChain->SetBranchAddress("mcParticleStartPositionY", mcParticleStartPositionY, &b_mcParticleStartPositionY);
   fChain->SetBranchAddress("mcParticleStartPositionZ", mcParticleStartPositionZ, &b_mcParticleStartPositionZ);
   fChain->SetBranchAddress("mcParticleStartPositionT", mcParticleStartPositionT, &b_mcParticleStartPositionT);
   fChain->SetBranchAddress("mcParticleStartMomentumX", mcParticleStartMomentumX, &b_mcParticleStartMomentumX);
   fChain->SetBranchAddress("mcParticleStartMomentumY", mcParticleStartMomentumY, &b_mcParticleStartMomentumY);
   fChain->SetBranchAddress("mcParticleStartMomentumZ", mcParticleStartMomentumZ, &b_mcParticleStartMomentumZ);
   fChain->SetBranchAddress("mcParticleStartMomentumE", mcParticleStartMomentumE, &b_mcParticleStartMomentumE);
   fChain->SetBranchAddress("mcParticleEndPositionX", mcParticleEndPositionX, &b_mcParticleEndPositionX);
   fChain->SetBranchAddress("mcParticleEndPositionY", mcParticleEndPositionY, &b_mcParticleEndPositionY);
   fChain->SetBranchAddress("mcParticleEndPositionZ", mcParticleEndPositionZ, &b_mcParticleEndPositionZ);
   fChain->SetBranchAddress("mcParticleEndPositionT", mcParticleEndPositionT, &b_mcParticleEndPositionT);
   fChain->SetBranchAddress("mcParticleEndMomentumX", mcParticleEndMomentumX, &b_mcParticleEndMomentumX);
   fChain->SetBranchAddress("mcParticleEndMomentumY", mcParticleEndMomentumY, &b_mcParticleEndMomentumY);
   fChain->SetBranchAddress("mcParticleEndMomentumZ", mcParticleEndMomentumZ, &b_mcParticleEndMomentumZ);
   fChain->SetBranchAddress("mcParticleEndMomentumE", mcParticleEndMomentumE, &b_mcParticleEndMomentumE);
   fChain->SetBranchAddress("mcParticleVertexTime", mcParticleVertexTime, &b_mcParticleVertexTime);
   fChain->SetBranchAddress("mcParticleEndTime", mcParticleEndTime, &b_mcParticleEndTime);
   fChain->SetBranchAddress("mcParticleNHits", mcParticleNHits, &b_mcParticleNHits);
   fChain->SetBranchAddress("mcParticleNHitsView", mcParticleNHitsView, &b_mcParticleNHitsView);
   fChain->SetBranchAddress("pfpTrueParticleMatchedID", pfpTrueParticleMatchedID, &b_pfpTrueParticleMatchedID);
   fChain->SetBranchAddress("pfpTrueParticleMatchedPosition", pfpTrueParticleMatchedPosition, &b_pfpTrueParticleMatchedPosition);
   fChain->SetBranchAddress("pfpIsPrimary", pfpIsPrimary, &b_pfpIsPrimary);
   fChain->SetBranchAddress("pfpID", pfpID, &b_pfpID);
   fChain->SetBranchAddress("pfpParentID", pfpParentID, &b_pfpParentID);
   fChain->SetBranchAddress("pfpPdgCode", pfpPdgCode, &b_pfpPdgCode);
   fChain->SetBranchAddress("pfpNChildren", pfpNChildren, &b_pfpNChildren);
   fChain->SetBranchAddress("pfpNClusters", pfpNClusters, &b_pfpNClusters);
   fChain->SetBranchAddress("pfpNHits", pfpNHits, &b_pfpNHits);
   fChain->SetBranchAddress("pfpNHitsView", pfpNHitsView, &b_pfpNHitsView);
   fChain->SetBranchAddress("pfpNSharedTrueParticleHits", pfpNSharedTrueParticleHits, &b_pfpNSharedTrueParticleHits);
   fChain->SetBranchAddress("pfpNSharedTrueParticleHitsView", pfpNSharedTrueParticleHitsView, &b_pfpNSharedTrueParticleHitsView);
   fChain->SetBranchAddress("pfpTrueParticleMatchedIDView", pfpTrueParticleMatchedIDView, &b_pfpTrueParticleMatchedIDView);
   fChain->SetBranchAddress("pfpTrueParticleMatchedPositionView", pfpTrueParticleMatchedPositionView, &b_pfpTrueParticleMatchedPositionView);
   fChain->SetBranchAddress("pfpIsTrack", pfpIsTrack, &b_pfpIsTrack);
   fChain->SetBranchAddress("pfpIsShower", pfpIsShower, &b_pfpIsShower);
   fChain->SetBranchAddress("pfpTrackID", pfpTrackID, &b_pfpTrackID);
   fChain->SetBranchAddress("pfpTrackLength", pfpTrackLength, &b_pfpTrackLength);
   fChain->SetBranchAddress("pfpTrackStartX", pfpTrackStartX, &b_pfpTrackStartX);
   fChain->SetBranchAddress("pfpTrackStartY", pfpTrackStartY, &b_pfpTrackStartY);
   fChain->SetBranchAddress("pfpTrackStartZ", pfpTrackStartZ, &b_pfpTrackStartZ);
   fChain->SetBranchAddress("pfpTrackVertexX", pfpTrackVertexX, &b_pfpTrackVertexX);
   fChain->SetBranchAddress("pfpTrackVertexY", pfpTrackVertexY, &b_pfpTrackVertexY);
   fChain->SetBranchAddress("pfpTrackVertexZ", pfpTrackVertexZ, &b_pfpTrackVertexZ);
   fChain->SetBranchAddress("pfpTrackEndX", pfpTrackEndX, &b_pfpTrackEndX);
   fChain->SetBranchAddress("pfpTrackEndY", pfpTrackEndY, &b_pfpTrackEndY);
   fChain->SetBranchAddress("pfpTrackEndZ", pfpTrackEndZ, &b_pfpTrackEndZ);
   fChain->SetBranchAddress("pfpTrackTheta", pfpTrackTheta, &b_pfpTrackTheta);
   fChain->SetBranchAddress("pfpTrackPhi", pfpTrackPhi, &b_pfpTrackPhi);
   fChain->SetBranchAddress("pfpTrackZenithAngle", pfpTrackZenithAngle, &b_pfpTrackZenithAngle);
   fChain->SetBranchAddress("pfpTrackAzimuthAngle", pfpTrackAzimuthAngle, &b_pfpTrackAzimuthAngle);
   fChain->SetBranchAddress("pfpTrackStartDirectionX", pfpTrackStartDirectionX, &b_pfpTrackStartDirectionX);
   fChain->SetBranchAddress("pfpTrackStartDirectionY", pfpTrackStartDirectionY, &b_pfpTrackStartDirectionY);
   fChain->SetBranchAddress("pfpTrackStartDirectionZ", pfpTrackStartDirectionZ, &b_pfpTrackStartDirectionZ);
   fChain->SetBranchAddress("pfpTrackVertexDirectionX", pfpTrackVertexDirectionX, &b_pfpTrackVertexDirectionX);
   fChain->SetBranchAddress("pfpTrackVertexDirectionY", pfpTrackVertexDirectionY, &b_pfpTrackVertexDirectionY);
   fChain->SetBranchAddress("pfpTrackVertexDirectionZ", pfpTrackVertexDirectionZ, &b_pfpTrackVertexDirectionZ);
   fChain->SetBranchAddress("pfpTrackEndDirectionX", pfpTrackEndDirectionX, &b_pfpTrackEndDirectionX);
   fChain->SetBranchAddress("pfpTrackEndDirectionY", pfpTrackEndDirectionY, &b_pfpTrackEndDirectionY);
   fChain->SetBranchAddress("pfpTrackEndDirectionZ", pfpTrackEndDirectionZ, &b_pfpTrackEndDirectionZ);
   fChain->SetBranchAddress("pfpTrackChi2", pfpTrackChi2, &b_pfpTrackChi2);
   fChain->SetBranchAddress("pfpTrackStartNdof", pfpTrackStartNdof, &b_pfpTrackStartNdof);
   fChain->SetBranchAddress("pfpCluPlane", pfpCluPlane, &b_pfpCluPlane);
   fChain->SetBranchAddress("pfpCluView", pfpCluView, &b_pfpCluView);
   fChain->SetBranchAddress("pfpCluNHits", pfpCluNHits, &b_pfpCluNHits);
   fChain->SetBranchAddress("pfpCluIntegral", pfpCluIntegral, &b_pfpCluIntegral);
   fChain->SetBranchAddress("pfpShowerID", pfpShowerID, &b_pfpShowerID);
   fChain->SetBranchAddress("pfpShowerBestPlane", pfpShowerBestPlane, &b_pfpShowerBestPlane);
   fChain->SetBranchAddress("pfpShowerEnergy", pfpShowerEnergy, &b_pfpShowerEnergy);
   fChain->SetBranchAddress("pfpShowerEnergyToTrueEnergyRatio", pfpShowerEnergyToTrueEnergyRatio, &b_pfpShowerEnergyToTrueEnergyRatio);
   fChain->SetBranchAddress("pfpShowerEnergyToTrueMomentumRatio", pfpShowerEnergyToTrueMomentumRatio, &b_pfpShowerEnergyToTrueMomentumRatio);
   fChain->SetBranchAddress("pfpShowerDirectionX", pfpShowerDirectionX, &b_pfpShowerDirectionX);
   fChain->SetBranchAddress("pfpShowerDirectionY", pfpShowerDirectionY, &b_pfpShowerDirectionY);
   fChain->SetBranchAddress("pfpShowerDirectionZ", pfpShowerDirectionZ, &b_pfpShowerDirectionZ);
   fChain->SetBranchAddress("pfpShowerDirectionErrX", pfpShowerDirectionErrX, &b_pfpShowerDirectionErrX);
   fChain->SetBranchAddress("pfpShowerDirectionErrY", pfpShowerDirectionErrY, &b_pfpShowerDirectionErrY);
   fChain->SetBranchAddress("pfpShowerDirectionErrZ", pfpShowerDirectionErrZ, &b_pfpShowerDirectionErrZ);
   fChain->SetBranchAddress("pfpShowerStartX", pfpShowerStartX, &b_pfpShowerStartX);
   fChain->SetBranchAddress("pfpShowerStartY", pfpShowerStartY, &b_pfpShowerStartY);
   fChain->SetBranchAddress("pfpShowerStartZ", pfpShowerStartZ, &b_pfpShowerStartZ);
   fChain->SetBranchAddress("pfpShowerStartErrX", pfpShowerStartErrX, &b_pfpShowerStartErrX);
   fChain->SetBranchAddress("pfpShowerStartErrY", pfpShowerStartErrY, &b_pfpShowerStartErrY);
   fChain->SetBranchAddress("pfpShowerStartErrZ", pfpShowerStartErrZ, &b_pfpShowerStartErrZ);
   fChain->SetBranchAddress("pfpShowerLength", pfpShowerLength, &b_pfpShowerLength);
   fChain->SetBranchAddress("pfpShowerOpeningAngle", pfpShowerOpeningAngle, &b_pfpShowerOpeningAngle);
   fChain->SetBranchAddress("pfpCompleteness", pfpCompleteness, &b_pfpCompleteness);
   fChain->SetBranchAddress("pfpCompletenessView", pfpCompletenessView, &b_pfpCompletenessView);
   fChain->SetBranchAddress("pfpPurity", pfpPurity, &b_pfpPurity);
   fChain->SetBranchAddress("pfpPurityView", pfpPurityView, &b_pfpPurityView);
   Notify();
}

Bool_t afmanatree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void afmanatree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t afmanatree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef afmanatree_cxx
