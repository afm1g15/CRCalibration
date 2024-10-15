#include <iostream>
using namespace TMVA::Experimental;

void tmvaReader()
{
   //Load in the model from the TMMA xml file
   RReader model("datasetBkg0/weights/TMVAMultiBkg0_BDTG.weights.xml");

   // In case you need a reminder of the names and order of the variables during
   // training, you can ask the model for it.
   auto variables = model.GetVariableNames();
   //std::cout << "vairables = " << variables << std::endl;

   //Apply model
   //Call as a lambda function to make the inference on a dataframe
   auto make_histo = [&](const std::string &treename, const std::string &filename) {
      ROOT::RDataFrame df(treename, filename);
      auto df2 = df.Define("y", Compute<7, float>(model), {"nvtx","trkthetaxz", "trkthetayz","length", "distEnter" ,"distExit", "trkstartd"});
      return df2.Histo1D({treename.c_str(), ";BDT score;N_{Events}", 60, -1, 1}, "y");
   };

   auto sig = make_histo("sigtree", "signal_full.root");
   auto bkg = make_histo("bkgtree", "background_full.root");

   //make a plot of it
   gStyle->SetOptStat(0);
   auto c = new TCanvas("", "", 800, 800);
 
   sig->SetLineColor(kBlue);
   bkg->SetLineColor(kRed);
   sig->SetLineWidth(2);
   bkg->SetLineWidth(2);
   bkg->Draw("HIST");
   sig->Draw("HIST SAME");
 
   TLegend legend(0.7, 0.7, 0.89, 0.89);
   legend.SetBorderSize(0);
   legend.AddEntry("TreeS", "Signal (Blue)", "l");
   legend.AddEntry("TreeB", "Background (Red)", "l");
   legend.Draw();
 
   c->DrawClone();
   c->SaveAs("TMVA.png");
}


