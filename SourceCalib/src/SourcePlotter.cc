#include "CaloCalibration/SourceCalib/inc/SourcePlotter.hh"
#include <chrono>
using namespace std::chrono;
using namespace CaloSourceCalib;

void SourcePlotter::ParamPlots(TTree* t) {
    TString newouptName = "param_distribution_plots.root";
    TFile *newouptFile = new TFile(newouptName, "RECREATE");
    TH1F *newoutHist = new TH1F("paramTree","paramTree",100,0.05,0.08);
    // Vectors to hold data
   //std::vector<double> crystalNos, Peaks;
    float  crystalNo;
    float Peak;
    float ChiSq;    
    t -> SetBranchAddress("crystalNo", &crystalNo);//SetBranchAddress
    t -> SetBranchAddress("Peak", &Peak); //SetBranchAddress
    t -> SetBranchAddress("ChiSq", &ChiSq); //SetBranchAddress  
    std::vector<double> crystalNos, Peaks,ChiSqs;   
    // Fill the vectors with data from the TTree
    Long64_t nEntries = t->GetEntries();
    for (Long64_t entry = 0; entry < nEntries; ++entry) {

        t->GetEntry(entry);
        crystalNos.push_back(crystalNo);
        Peaks.push_back(Peak);
        ChiSqs.push_back(ChiSq);        
        newoutHist->Fill(Peak, crystalNo);//,ChiSq);
    }

    // Create a TGraph for the scatter plot
    TGraph graph(crystalNos.size(), crystalNos.data(), Peaks.data());
		//TGraph *graph1,*graph2
		
    // Customize the plot
    graph.SetTitle("Cry Number vs Full Peak;Crystal Number;Full Peak");
    graph.SetMarkerStyle(20);
    graph.SetMarkerSize(0.8);
    graph.SetMarkerColor(kBlue);
    
    //graph2SetTitle("Cry Number vs Chisq ;Crystal Number; Chi Square");
    //graph2.SetMarkerStyle(20);
    //graph2.SetMarkerSize(0.8);
    //graph2.SetMarkerColor(kBlue);
    // Create a canvas to draw the plot
    TCanvas canvas("canvas", "Scatter Plot", 800, 600);
    graph.Draw("AP"); // A for axis, P for points
    
    TGraph graph2(crystalNos.size(), crystalNos.data(), ChiSqs.data());
		//TGraph *graph1,*graph2
		
    // Customize the plot
    graph2.SetTitle("Cry Number vs Chi Square;Crystal Number;Chi Square");
    graph2.SetMarkerStyle(20);
    graph2.SetMarkerSize(0.8);
    graph2.SetMarkerColor(kBlue);
    
    //graph2SetTitle("Cry Number vs Chisq ;Crystal Number; Chi Square");
    //graph2.SetMarkerStyle(20);
    //graph2.SetMarkerSize(0.8);
    //graph2.SetMarkerColor(kBlue);
    // Create a canvas to draw the plot
    TCanvas canvas2("canvas2", "Scatter Plot2", 800, 600);
    graph2.Draw("AP"); // A for axis, P for points

    // Save the canvas
    canvas.SaveAs("scatter_plot2.pdf");
    canvas2.SaveAs("scatter_plot2.pdf");
 	  newouptFile -> Write();
  	newouptFile -> Close();

}
