#include "CaloCalibration/SourceCalib/inc/SourcePlotter.hh"
#include <chrono>
using namespace std::chrono;
using namespace CaloSourceCalib;

/* function to make global plots of the fit outputs*/
void SourcePlotter::ParamPlots(TTree* inputTree, TFile *inputFile, TFile *outputFile) {
    float crystalNo, Peak, ChiSq, PeakErr;
    int nEvents;    
    inputTree-> SetBranchAddress("crystalNo", &crystalNo);
    inputTree-> SetBranchAddress("Peak", &Peak);
    inputTree-> SetBranchAddress("PeakErr", &PeakErr);
    inputTree-> SetBranchAddress("ChiSq", &ChiSq);
    inputTree-> SetBranchAddress("nEvents", &nEvents); 
    std::vector<Double_t> crystalNos, Peaks, PeakErrs, ChiSqs, EventsErr, Events;
     
    // Fill the vectors with data from the TTree
    Long64_t nEntries = inputTree->GetEntries();
    for (Long64_t entry = 0; entry < nEntries; ++entry) {
        inputTree->GetEntry(entry);
        crystalNos.push_back(crystalNo);
        Peaks.push_back(1/(Peak/6.13));
        PeakErrs.push_back(PeakErr);
        ChiSqs.push_back(ChiSq);
        Events.push_back(nEvents);
        EventsErr.push_back(sqrt(nEvents));       
    }

    // Create a TGraph for the scatter plot
    TGraphErrors grpeaks(crystalNos.size(), &crystalNos[0], &Peaks[0], nullptr, &PeakErrs[0]);
    TGraph grchi2(crystalNos.size(), &crystalNos[0], &ChiSqs[0]);
    TGraphErrors grN(crystalNos.size(), &crystalNos[0], &Events[0], nullptr, &EventsErr[0]);

    // Customize the plot
    grpeaks.SetTitle("Cry Number vs MeV/ADC;Crystal Number;MeV/ADC");
    grpeaks.SetMarkerStyle(20);
    grpeaks.SetMarkerSize(0.8);
    grpeaks.SetMarkerColor(kBlue);
    grpeaks.SetLineColor(kBlue);
    TCanvas canvas("canvas", "Scatter Plot", 800, 600);
    grpeaks.Draw("AP"); // A for axis, P for points
    outputFile->cd();
    grpeaks.Write("Peaks");

    grchi2.SetTitle("Cry Number vs Chisq ;Crystal Number; Chi Square");
    grchi2.SetMarkerStyle(20);
    grchi2.SetMarkerSize(0.8);
    grchi2.SetMarkerColor(kRed);
    grchi2.SetLineColor(kRed);
    TCanvas canvas2("canvas", "Scatter Plot", 800, 600);
    grchi2.Draw("AP"); // A for axis, P for points
    outputFile->cd();
    grchi2.Write("Chi2");
    
    grN.SetTitle("Cry Number vs nEvents ;Crystal Number; nEvents");
    grN.SetMarkerStyle(20);
    grN.SetMarkerSize(0.8);
    grN.SetMarkerColor(kGreen);
    grN.SetLineColor(kGreen);
    TCanvas canvas3("canvas", "Scatter Plot", 800, 600);
    grN.Draw("AP"); // A for axis, P for points
    outputFile->cd();
    grN.Write("nEvents");
}
