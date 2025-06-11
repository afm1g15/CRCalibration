#include <iostream>
#include <vector>
#include <limits>
 
#include "TChain.h"
#include "TCut.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TH1.h"
#include "TMath.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TROOT.h"
#include "TSystem.h"
 
#include "TMVA/GeneticAlgorithm.h"
#include "TMVA/GeneticFitter.h"
#include "TMVA/IFitterTarget.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Reader.h"

using std::vector, std::cout, std::endl;
 
using namespace TMVA;

//------------------------------------------------
//-----------Training-----------------------------
//------------------------------------------------

void Training(){
   std::string factoryOptions( "!V:!Silent:Transformations=I;D;P;G,D:AnalysisType=Classification" );
   //TString fname = "./tmva_example_multiple_background.root";

   std::unique_ptr<TFile> myFileSig( TFile::Open("signal_full.root") );
   std::unique_ptr<TFile> myFileBkg( TFile::Open("background_full.root") );
 
   TTree *signal      = myFileSig->Get<TTree>("sigtree");
   TTree *background0 = myFileBkg->Get<TTree>("bkgtree");

   //Global event weights per tree
   Double_t signalWeight      = 1.0;
   Double_t background0Weight = 1.0;

   // Create a new root output file.
   TString outfileName( "TMVASignalBackground0.root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   //background
   TMVA::Factory *factory = new TMVA::Factory( "TMVAMultiBkg0", outputFile, factoryOptions );
   TMVA::DataLoader *dataloader=new TMVA::DataLoader("datasetBkg0");

   dataloader->AddVariable( "nvtx", "N Vertex", "", 'F' );
   dataloader->AddVariable( "trkthetaxz", "Theta XZ", "", 'F' );
   dataloader->AddVariable( "trkthetayz", "Theta YZ", "", 'F' );
   dataloader->AddVariable( "length", "Length", "cm", 'F' );
   dataloader->AddVariable( "distEnter", "Dist Enter", "", 'F' );
   dataloader->AddVariable( "distExit", "Dist Exit", "", 'F' );
   dataloader->AddVariable( "trkstartd", "Start Dist", "", 'F' );
   dataloader->AddVariable( "starty", "Start Y", "", 'F' );
   dataloader->AddVariable( "endy", "End Y", "", 'F' );
   dataloader->AddVariable( "startx", "Start X", "", 'F' );
   dataloader->AddVariable( "endx", "End X", "", 'F' );
   dataloader->AddVariable( "startz", "Start Z", "", 'F' );
   dataloader->AddVariable( "endz", "End Z", "", 'F' );
   dataloader->AddVariable( "ke", "KE", "", 'F' );
   dataloader->AddVariable( "range", "range", "", 'F' );

   dataloader->AddSignalTree    ( signal,     signalWeight       );
   dataloader->AddBackgroundTree( background0, background0Weight );

   //cuts already applied before trees
   TCut mycuts = "";
   TCut mycutsb = "";

   // tell the factory to use all remaining events in the trees after training for testing:
   //dataloader->PrepareTrainingAndTestTree( mycuts, mycutsb,
     //                                   "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );
   dataloader->PrepareTrainingAndTestTree( mycuts,"nTrain_Signal=3585:nTrain_Background=1321:SplitMode=Random:!V" );  //274, 66
   // Boosted Decision Trees
      factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTG",  //G
         "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.30:UseBaggedBoost:BaggedSampleFraction=0.6:SeparationType=GiniIndex:nCuts=20:MaxDepth=2" );
   factory->TrainAllMethods();
   factory->TestAllMethods();
   factory->EvaluateAllMethods();
 
   outputFile->Close();
 
   delete factory;
   delete dataloader;
}


//---------------------------------------------
//-----------Application------------------------
//--------------------------------------------

//Create summary tree with classifer values

void ApplicationCreateCombinedTree(){

   // Create a new root output file.
   TString outfileName( "tmva_example_multiple_backgrounds__applied.root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
   TTree* outputTree = new TTree("multiBkg","multiple backgrounds tree");  


   float trkpurity, trkthetaxz, trkthetayz, completeness, length, ke, range;
   float nvtx;
   float distExit, trkstartd, starty, endy, startx, endx, startz, distEnter, endz;

   Int_t   classID = 0;
   Float_t weight = 1.f;

   Float_t classifier0;

   outputTree->Branch("classID", &classID, "classID/I");
   outputTree->Branch( "nvtx", &nvtx, "nvtx/F" );
   outputTree->Branch( "trkthetaxz", &trkthetaxz, "trkthetaxz/F");
   outputTree->Branch( "trkthetayz", &trkthetayz, "trkthetayz/F");
   outputTree->Branch( "length", &length, "length/F");
   outputTree->Branch( "distEnter", &distEnter, "distEnter/F");
   outputTree->Branch( "distExit", &distExit, "distExit/F");
   outputTree->Branch( "trkstartd", &trkstartd, "trkstartd/F");
   outputTree->Branch( "starty", &starty, "starty/F");
   outputTree->Branch( "endy", &endy, "endy/F");
   outputTree->Branch( "startx", &startx, "startx/F");
   outputTree->Branch( "endx", &endx, "endx/F");
   outputTree->Branch( "startz", &startz, "startz/F");
   outputTree->Branch( "endz", &endz, "endz/F");
   outputTree->Branch( "ke", &ke, "ke/F");
   outputTree->Branch( "range", &range, "range/F");
   outputTree->Branch("weight", &weight, "weight/F");
   outputTree->Branch("cls0", &classifier0, "cls0/F");

   //add readers
   TMVA::Reader *reader0 = new TMVA::Reader( "!Color:!Silent" );
   reader0->AddVariable( "nvtx", &nvtx);
   reader0->AddVariable( "trkthetaxz", &trkthetaxz);
   reader0->AddVariable( "trkthetayz", &trkthetayz);
   reader0->AddVariable( "length", &length);
   reader0->AddVariable( "distEnter", &distEnter);
   reader0->AddVariable( "distExit", &distExit);
   reader0->AddVariable( "trkstartd", &trkstartd);
   reader0->AddVariable( "starty", &starty);
   reader0->AddVariable( "endy", &endy);
   reader0->AddVariable( "startx", &startx);
   reader0->AddVariable( "endx", &endx);
   reader0->AddVariable( "startz", &startz);
   reader0->AddVariable( "endz", &endz);
   reader0->AddVariable( "ke", &ke);
   reader0->AddVariable( "range", &range);

   //load readers
   TString method =  "BDT method";
   reader0->BookMVA( "BDT method", "datasetBkg0/weights/TMVAMultiBkg0_BDTG.weights.xml" );

   //load file
   std::unique_ptr<TFile> myFileSig( TFile::Open("signal_full.root") );
   std::unique_ptr<TFile> myFileBkg( TFile::Open("background_full.root") );

   TTree* theTree = NULL;

   //loop through all trees
   for( int treeNumber = 0; treeNumber < 2; ++treeNumber ) {
      if( treeNumber == 0 ){
    theTree = (TTree*)myFileSig->Get("sigtree");
    std::cout << "--- Select signal sample" << std::endl;
    weight = 1;
    classID = 0;
      }else if( treeNumber == 1 ){
    theTree = (TTree*)myFileBkg->Get("bkgtree");
    std::cout << "--- Select background 0 sample" << std::endl;
    weight = 1;
    classID = 1;
   }


   theTree->SetBranchAddress( "nvtx", &nvtx);
   theTree->SetBranchAddress( "trkthetaxz", &trkthetaxz);
   theTree->SetBranchAddress( "trkthetayz", &trkthetayz);
   theTree->SetBranchAddress( "length", &length);
   theTree->SetBranchAddress( "distEnter", &distEnter);
   theTree->SetBranchAddress( "distExit", &distExit);
   theTree->SetBranchAddress( "trkstartd", &trkstartd);   //filtered to here
   theTree->SetBranchAddress( "starty", &starty);
   theTree->SetBranchAddress( "endy", &endy);
   theTree->SetBranchAddress( "startx", &startx);
   theTree->SetBranchAddress( "endx", &endx);
   theTree->SetBranchAddress( "startz", &startz);
   theTree->SetBranchAddress( "endz", &endz);
   theTree->SetBranchAddress( "ke", &ke);
   theTree->SetBranchAddress( "range", &range);
   

    std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
    TStopwatch sw;
    sw.Start();
    Int_t nEvent = theTree->GetEntries();
    for (Long64_t ievt=0; ievt<nEvent; ievt++) {
 
    if (ievt%1000 == 0){
       std::cout << "--- ... Processing event: " << ievt << std::endl;
    }
 
    theTree->GetEntry(ievt);

    classifier0 = reader0->EvaluateMVA( method );
    outputTree->Fill();
    }

      sw.Stop();
      std::cout << "--- End of event loop: "; sw.Print();
   }
   myFileSig->Close();
   myFileBkg->Close();

   outputFile->Write();
 
   outputFile->Close();
 
   std::cout << "--- Created root file: \"" << outfileName.Data() << "\" containing the MVA output histograms" << std::endl;
 
   delete reader0;

   std::cout << "==> Application of readers is done! combined tree created" << std::endl << std::endl;
 
}

//-------------------------------------------------
//------Genertic Alg Fitness Def-------------------
//-------------------------------------------------

class MyFitness : public IFitterTarget {
public:
   // constructor
      MyFitness( TChain* _chain ) : IFitterTarget() {
      chain = _chain;
 
      hSignal = new TH1F("hsignal","hsignal",100,-1,1);
      hFP = new TH1F("hfp","hfp",100,-1,1);
      hTP = new TH1F("htp","htp",100,-1,1);
 
      TString cutsAndWeightSignal  = "weight*(classID==0)";
      nSignal = chain->Draw("Entry$/Entries$>>hsignal",cutsAndWeightSignal,"goff");
      weightsSignal = hSignal->Integral();
 
   }

   // the output of this function will be minimized
      Double_t EstimatorFunction( std::vector<Double_t> & factors ){
 
      TString cutsAndWeightTruePositive  = Form("weight*((classID==0) && cls0>%f)",factors.at(0));
      TString cutsAndWeightFalsePositive = Form("weight*((classID >0) && cls0>%f)",factors.at(0));
 
      // Entry$/Entries$ just draws something reasonable. Could in principle anything
            Float_t nTP = chain->Draw("Entry$/Entries$>>htp",cutsAndWeightTruePositive,"goff");
      Float_t nFP = chain->Draw("Entry$/Entries$>>hfp",cutsAndWeightFalsePositive,"goff");
 
      weightsTruePositive = hTP->Integral();
      weightsFalsePositive = hFP->Integral();
 
      efficiency = 0;
      if( weightsSignal > 0 )
    efficiency = weightsTruePositive/weightsSignal;
 
      purity = 0;
      if( weightsTruePositive+weightsFalsePositive > 0 )
    purity = weightsTruePositive/(weightsTruePositive+weightsFalsePositive);
 
      Float_t effTimesPur = efficiency*purity;
 
      Float_t toMinimize = std::numeric_limits<float>::max(); // set to the highest existing number
      if( effTimesPur > 0 ) // if larger than 0, take 1/x. This is the value to minimize
    toMinimize = 1./(effTimesPur); // we want to minimize 1/efficiency*purity

      return toMinimize;
   }
 
 
   void Print(){
      std::cout << std::endl;
      std::cout << "======================" << std::endl
      << "Efficiency : " << efficiency << std::endl
      << "Purity     : " << purity << std::endl << std::endl
      << "True positive weights : " << weightsTruePositive << std::endl
      << "False positive weights: " << weightsFalsePositive << std::endl
      << "Signal weights        : " << weightsSignal << std::endl;
   }
 
   Float_t nSignal;
 
   Float_t efficiency;
   Float_t purity;
   Float_t weightsTruePositive;
   Float_t weightsFalsePositive;
   Float_t weightsSignal;
 
 
private:
   TChain* chain;
   TH1F* hSignal;
   TH1F* hFP;
   TH1F* hTP;
 
};

//------------------------------------------------
//-------------Call generic alg-------------------
//------------------------------------------------

void MaximizeSignificance(){
 
        // define all the parameters by their minimum and maximum value
        // in this example 3 parameters (=cuts on the classifiers) are defined.
        vector<Interval*> ranges;
        ranges.push_back( new Interval(-1,1) ); // for some classifiers (especially LD) the ranges have to be taken larger
 
   std::cout << "Classifier ranges (defined by the user)" << std::endl;
        for( std::vector<Interval*>::iterator it = ranges.begin(); it != ranges.end(); it++ ){
           std::cout << " range: " << (*it)->GetMin() << "   " << (*it)->GetMax() << std::endl;
        }
 
   TChain* chain = new TChain("multiBkg");
   chain->Add("tmva_example_multiple_backgrounds__applied.root");
 
        IFitterTarget* myFitness = new MyFitness( chain );

        const TString name( "multipleBackgroundGA" );
        const TString opts( "PopSize=100:Steps=30" );
 
        GeneticFitter mg( *myFitness, name, ranges, opts);

        std::vector<Double_t> result;
        Double_t estimator = mg.Run(result);
 
   dynamic_cast<MyFitness*>(myFitness)->Print();
   std::cout << std::endl;
 
   int n = 0;
   for( std::vector<Double_t>::iterator it = result.begin(); it<result.end(); it++ ){
      std::cout << "  cutValue[" << n << "] = " << (*it) << ";"<< std::endl;
      n++;
   }
 
 
}


//------------------------------------------------
//----------------Run-----------------------------
//-----------------------------------------------

void trainTMVABDT()
{   cout << "Start Test TMVAGAexample" << endl
        << "========================" << endl
        << endl;
 
   cout << endl;
   cout << "========================" << endl;
   cout << "--- Training" << endl;
   Training();
 
   cout << endl;
   cout << "========================" << endl;
   cout << "--- Application & create combined tree" << endl;
   ApplicationCreateCombinedTree();
 
   cout << endl;
   cout << "========================" << endl;
   cout << "--- maximize significance" << endl;
   MaximizeSignificance();
}
 
//int trainTVMABDT( int argc, char** argv ) {
//   TMVAMultipleBackgroundExample();
//}
