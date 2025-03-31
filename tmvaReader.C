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
      auto df2 = df.Define("y", Compute<15, float>(model), {"nvtx","trkthetaxz", "trkthetayz","length", "distEnter","distExit", "trkstartd", "starty", "endy", "startx", "endx", "startz", "endz", "ke", "range"}); //"length", "distEnter", "endz"
      return df2.Histo1D({treename.c_str(), ";BDT score;N_{Events}", 60, -1, 1}, "y");
   };

   auto sig = make_histo("sigtree", "signal_full.root");
   auto bkg = make_histo("bkgtree", "background_full.root");

   //make a plot of it
   gStyle->SetOptStat(0);
   auto c = new TCanvas("", "", 1000, 1000);
   c->SetLeftMargin(0.12);
   c->SetRightMargin(0.06);
   c->SetTopMargin(0.06);
   c->SetBottomMargin(0.12);
 
   sig->SetLineColor(kBlue);
   bkg->SetLineColor(kRed);
   sig->SetLineWidth(2);
   bkg->SetLineWidth(2);
   sig->GetXaxis()->SetLabelSize(0.03);
   bkg->GetXaxis()->SetLabelSize(0.03);
   sig->GetYaxis()->SetLabelSize(0.03);
   bkg->GetYaxis()->SetLabelSize(0.03);
   sig->GetXaxis()->SetTitleSize(0.04);
   bkg->GetXaxis()->SetTitleSize(0.04);
   sig->GetYaxis()->SetTitleSize(0.04);
   bkg->GetYaxis()->SetTitleSize(0.04);
   sig->Draw("HIST");
   bkg->Draw("HIST SAME");
 
   TLegend legend(0.6, 0.6, 0.79, 0.79);
   legend.SetBorderSize(0);
   legend.SetTextSize(0.04);
   legend.AddEntry("sigtree", "Signal", "l");
   legend.AddEntry("bkgtree", "Background", "l");
   legend.Draw();
   //c->BuildLegend();
 
   c->DrawClone();
   c->SaveAs("TMVA.png");
}


