#include "CaloCalibration/SourceCalib/inc/SourcePlotter.hh"
#include <chrono>
using namespace std::chrono;
using namespace CaloSourceCalib;

/* function to make global plots of the fit outputs*/
void SourcePlotter::ParamPlots(TTree* inputTree, TFile *outputFile) {
    float  crystalNo;
    float Peak;
    float ChiSq;
    float PeakErr;    
    inputTree-> SetBranchAddress("crystalNo", &crystalNo);//SetBranchAddress
    inputTree-> SetBranchAddress("Peak", &Peak); //SetBranchAddress
    inputTree-> SetBranchAddress("PeakErr", &PeakErr); //SetBranchAddress
    inputTree-> SetBranchAddress("ChiSq", &ChiSq); //SetBranchAddress  
    std::vector<Double_t> crystalNos, Peaks, PeakErrs, ChiSqs;
     
    // Fill the vectors with data from the TTree
    Long64_t nEntries = inputTree->GetEntries();
    for (Long64_t entry = 0; entry < nEntries; ++entry) {
        inputTree->GetEntry(entry);
        crystalNos.push_back(crystalNo);
        Peaks.push_back(Peak);
        PeakErrs.push_back(PeakErr);
        ChiSqs.push_back(ChiSq);        
    }

    // Create a TGraph for the scatter plot
    TGraphErrors graph(crystalNos.size(), &crystalNos[0], &Peaks[0], nullptr, &PeakErrs[0]);

    // Customize the plot
    graph.SetTitle("Cry Number vs Full Peak;Crystal Number;Full Peak");
    graph.SetMarkerStyle(20);
    graph.SetMarkerSize(0.8);
    graph.SetMarkerColor(kBlue);
    graph.SetLineColor(kBlue);
    TCanvas canvas("canvas", "Scatter Plot", 800, 600);
    graph.Draw("AP"); // A for axis, P for points
    outputFile->cd();
    graph.Write("Peaks");
    canvas.SaveAs("Peaks.root");

    //graph2SetTitle("Cry Number vs Chisq ;Crystal Number; Chi Square");
    //graph2.SetMarkerStyle(20);
    //graph2.SetMarkerSize(0.8);
    //graph2.SetMarkerColor(kBlue);
    // Create a canvas to draw the plot
    
    //TGraph graph2(crystalNos.size(), crystalNos.data(), ChiSqs.data());
		//TGraph *graph1,*graph2
		
    // Customize the plot
    //graph2.SetTitle("Cry Number vs Chi Square;Crystal Number;Chi Square");
    //graph2.SetMarkerStyle(20);
    //graph2.SetMarkerSize(0.8);
    //graph2.SetMarkerColor(kBlue);
    
    //graph2SetTitle("Cry Number vs Chisq ;Crystal Number; Chi Square");
    //graph2.SetMarkerStyle(20);
    //graph2.SetMarkerSize(0.8);
    //graph2.SetMarkerColor(kBlue);
    // Create a canvas to draw the plot
    //TCanvas canvas2("canvas2", "Scatter Plot2", 800, 600);
    //graph2.Draw("AP"); // A for axis, P for points

    // Save the canvas
    
    //canvas2.SaveAs("scatter_plot2.pdf");
 	  outputFile -> Write();
  	outputFile -> Close();
}
