#include "CaloCalibration/SourceCalib/inc/SourcePlotter.hh"
#include <chrono>
#include <numeric>
#include "TLine.h"
#include "TLegend.h"
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
		double Peaks_sum = std::accumulate(Peaks.begin(), Peaks.end(),0.0);
		double Peaks_avg = Peaks_sum / Peaks.size();
		double Peaks_stddev = TMath::StdDev(Peaks.begin(),Peaks.end());
		std::cout<<" nominal "<<0.0625<<std::endl;
		std::cout<<" average "<<Peaks_avg<<std::endl;
		std::cout<<" peak std dev "<<Peaks_stddev<<std::endl;
		double ChiSqs_sum = std::accumulate(ChiSqs.begin(), ChiSqs.end(),0.0);
		double ChiSqs_avg = ChiSqs_sum / ChiSqs.size();
		//double ChiSqs_stddev = TMath::StdDev(Peaks.begin(),Peaks.end());
	
		double Events_sum = std::accumulate(Events.begin(), Events.end(),0.0);
		double Events_avg = Events_sum / Events.size();
		double Events_stddev = TMath::StdDev(Events.begin(),Events.end());	
		std::cout<<" std dev of events "<<Events_stddev<<std::endl;			
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
    TLine *avgpeak = new TLine(674,Peaks_avg,1348,Peaks_avg);
    TLine *nominalpeak = new TLine(674,0.0625,1348,0.0625); 
    TLine *stddevpeak_abv = new TLine(674,Peaks_stddev,1348,Peaks_stddev);    
    TLine *stddevpeak_blw = new TLine(674,-Peaks_stddev,1348,-Peaks_stddev);       
    TCanvas canvas("canvas", "Scatter Plot", 800, 600);
    grpeaks.Draw("AP"); // A for axis, P for points
    avgpeak->SetLineColor(kCyan);
    nominalpeak->SetLineColor(kMagenta);
    stddevpeak_abv->SetLineColor(kOrange-3);
    stddevpeak_blw->SetLineColor(kOrange);        
    avgpeak->Draw("same");
    nominalpeak->Draw("same");
    auto legend = new TLegend(0.1,0.7,0.48,0.9);
    legend->AddEntry(avgpeak, "Average Value", "L");
    legend->AddEntry(nominalpeak,"Nominal Value","L");
    legend->AddEntry(stddevpeak_abv,"Std Dev","L");
    legend->AddEntry(stddevpeak_blw,"Std Dev","L");
    legend->Draw();
    outputFile->cd();
    grpeaks.Write("Peaks");
    canvas.SaveAs("Peaks.root");


    grchi2.SetTitle("Cry Number vs Chisq ;Crystal Number; Chi Square");
    grchi2.SetMarkerStyle(20);
    grchi2.SetMarkerSize(0.8);
    grchi2.SetMarkerColor(kRed);
    grchi2.SetLineColor(kRed);
    TLine *avgchi2 = new TLine(674,ChiSqs_avg,1348,ChiSqs_avg);
    TLine *optchi2 = new TLine(674,1,1348,1);
    //TLine *stddevchi2_abv = new TLine(674,ChiSqs_stddev,1348,ChiSqs_stddev);    
    //TLine *stddevchi2_blw = new TLine(674,-ChiSqs_stddev,1348,-ChiSqs_stddev);   
    TCanvas canvas2("canvas", "Scatter Plot", 800, 600);
    grchi2.Draw("AP"); // A for axis, P for points
    avgchi2->SetLineColor(kCyan);
    optchi2->SetLineColor(kMagenta);
    //stddevchi2_abv->SetLineColor(kOrange-3);
    //stddevchi2_blw->SetLineColor(kOrange-3);        
    avgchi2->Draw("same");
    optchi2->Draw("same");    
    //stddevchi2_abv->Draw("same");
    //stddevchi2_blw->Draw("same");
    auto legend_chi2 = new TLegend(0.1,0.7,0.48,0.9);
    legend_chi2->AddEntry(avgchi2, "Average Value", "L");
    legend_chi2->AddEntry(optchi2, "Optimal Value", "L");
    //legend_chi2->AddEntry(stddevchi2_abv,"Std Dev","L");
    legend_chi2->Draw();
    outputFile->cd();
    grchi2.Write("Chi2");
    canvas2.SaveAs("Chi2.root");    
    
    grN.SetTitle("Cry Number vs nEvents ;Crystal Number; nEvents");
    grN.SetMarkerStyle(20);
    grN.SetMarkerSize(0.8);
    grN.SetMarkerColor(kGreen);
    grN.SetLineColor(kGreen);
    TLine *avgevts = new TLine(674,Events_avg,1348,Events_avg);
    TLine *stddevevts_abv = new TLine(674,Events_stddev,1348,Events_stddev);    
    TLine *stddevevts_blw = new TLine(674,-Events_stddev,1348,-Events_stddev);   
    TCanvas canvas3("canvas", "Scatter Plot", 800, 600);
    grN.Draw("AP"); // A for axis, P for points
    avgevts->SetLineColor(kCyan);
    stddevevts_abv->SetLineColor(kOrange-3);
    stddevevts_blw->SetLineColor(kOrange-3);        
    avgevts->Draw("same");
    stddevevts_abv->Draw("same");
    stddevevts_blw->Draw("same");
    auto legend_evts = new TLegend(0.1,0.7,0.48,0.9);
    legend_evts->AddEntry(avgevts, "Average Value", "L");
    legend_evts->AddEntry(stddevevts_abv,"Std Dev","L");
    legend_evts->Draw();
    outputFile->cd();
    grN.Write("nEvents");
    canvas3.SaveAs("nEvents.root");    
}
