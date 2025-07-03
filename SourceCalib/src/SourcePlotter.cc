#include "CaloCalibration/SourceCalib/inc/SourcePlotter.hh"
#include <chrono>
#include <numeric>
#include "TF1.h"

using namespace std::chrono;
using namespace CaloSourceCalib;

/* function to make global plots of the fit outputs*/
void SourcePlotter::ParamPlots(TTree* inputTree, TFile *inputFile, TFile *outputFile,int cry_start, int cry_end) {//add param as plot range
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
        PeakErrs.push_back((1/(Peak/6.13))*(PeakErr/Peak));
        ChiSqs.push_back(ChiSq);
        Events.push_back(nEvents);
        EventsErr.push_back(sqrt(nEvents));       
    }
		double Peaks_sum = std::accumulate(Peaks.begin(), Peaks.end(),0.0);
 		double Peaks_avg = Peaks_sum / Peaks.size();
 		double Peaks_stddev = TMath::StdDev(Peaks.begin(),Peaks.end());
 		double stddev_prcnt = (Peaks_stddev/Peaks_avg)*100;
 		std::cout<<" peak 2 std dev "<<Peaks_stddev*2<<std::endl;
 		std::cout<<" peak std dev "<<Peaks_stddev<<std::endl;
 		std::cout<<" ------------ "<<std::endl; 		
 		int countHighPeakErrs = 0;
for (size_t i = 0; i < PeakErrs.size(); ++i) {
    if (PeakErrs[i] > 2 * Peaks_stddev) {
        ++countHighPeakErrs;
        std::cout << "Crystal " << crystalNos[i] 
                  << " has PeakErr = " << PeakErrs[i] 
                  << " > 2*Peaks_stddev (" << 2 * Peaks_stddev << ")\n";
    }
}
std::cout << "Number of PeakErrs > 2 stdev " << countHighPeakErrs << std::endl;

 		std::cout<<" ------------ "<<std::endl;
 		double ChiSqs_sum = std::accumulate(ChiSqs.begin(), ChiSqs.end(),0.0);
 		double ChiSqs_avg = ChiSqs_sum / ChiSqs.size();
 		double ChiSqs_stddev = TMath::StdDev(ChiSqs.begin(),ChiSqs.end());
 		std::cout<<" std dev of chi2 "<<ChiSqs_stddev<<std::endl;		
 		std::cout<<" 2 std dev of chi2 "<<2*ChiSqs_stddev<<std::endl;		
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
    TLine *avgpeak = new TLine(cry_start,Peaks_avg,cry_end,Peaks_avg);
    TLine *nominalpeak = new TLine(cry_start,0.0625,cry_end,0.0625); 
    TLine *stddevpeak_abv = new TLine(cry_start,Peaks_avg+2*Peaks_stddev,cry_end,Peaks_avg+2*Peaks_stddev);    
    TLine *stddevpeak_blw = new TLine(cry_start,Peaks_avg-2*Peaks_stddev,cry_end,Peaks_avg-2*Peaks_stddev);		
    TCanvas canvas("canvas", "Scatter Plot", 800, 600);
    grpeaks.Draw("AP"); // A for axis, P for points    
    TPaveText *avg = new TPaveText(0.15, 0.75, 0.35, 0.65, "brNDC");
   	avg -> SetFillStyle(0);
   	avg -> SetBorderSize(0);
   	avg -> SetTextSize(0.03);
   	avg -> SetTextColor(kBlack);
   	avg -> SetTextFont(72);
   	avg -> SetFillColor(kWhite);
    avg->AddText(Form("Peaks avg val = %4.8f#pm%.8f", Peaks_avg,Peaks_stddev));
		avg->AddText(Form("Std Dev = %.2f%%", stddev_prcnt));   
    avgpeak->SetLineColor(kCyan);
    avgpeak->SetLineWidth(2);
    nominalpeak->SetLineColor(kMagenta);
    nominalpeak->SetLineWidth(2);
    stddevpeak_abv->SetLineColor(kOrange-3);
    stddevpeak_blw->SetLineColor(kOrange);        
    stddevpeak_abv->SetLineWidth(2);    
    stddevpeak_blw->SetLineColor(kOrange);
    stddevpeak_blw->SetLineWidth(2);        
    avgpeak->Draw("same");
    nominalpeak->Draw("same");
    stddevpeak_abv->Draw("same");
    stddevpeak_blw->Draw("same");
    avg->Draw();
    auto legend = new TLegend();
    legend->AddEntry(avgpeak, "Average Value", "L");
    legend->AddEntry(nominalpeak,"Nominal Value","L");
    legend->AddEntry(stddevpeak_abv,"Std Dev","L");
    legend->AddEntry(stddevpeak_blw,"Std Dev","L");
     //legend->SetTextSize(0.4);
    legend->Draw();
    outputFile->cd();
    grpeaks.Write("Peaks");
    canvas.SaveAs("Peaks.root");
          
    //TCanvas canvas("canvas", "Scatter Plot", 800, 600);
    grpeaks.Draw("APsame"); // A for axis, P for points
    outputFile->cd();
    grpeaks.Write("Peaks");

    
    grchi2.SetTitle("SiPM Number vs Chisq ;SiPM Number; Chi Square");
    grchi2.SetMarkerStyle(20);
    grchi2.SetMarkerSize(0.8);
    grchi2.SetMarkerColor(kRed);
    grchi2.SetLineColor(kRed);
    TLine *avgchi2 = new TLine(cry_start,ChiSqs_avg,cry_end,ChiSqs_avg);
    TLine *optchi2 = new TLine(cry_start,1,cry_end,1);
    TLine *stddevchi2_abv = new TLine(cry_start,ChiSqs_avg+2*ChiSqs_stddev,cry_end,ChiSqs_avg+2*ChiSqs_stddev);    
    TLine *stddevchi2_blw = new TLine(cry_start,ChiSqs_avg-2*ChiSqs_stddev,cry_end,ChiSqs_avg-2*ChiSqs_stddev);
    TCanvas canvas2("canvas", "Scatter Plot", 800, 600);
    grchi2.Draw("AP"); // A for axis, P for points
    TPaveText *chi2text = new TPaveText(0.15, 0.75, 0.35, 0.65, "brNDC");
   	chi2text -> SetFillStyle(0);
   	chi2text -> SetBorderSize(0);
   	chi2text -> SetTextSize(0.03);
   	chi2text -> SetTextColor(kBlack);
   	chi2text -> SetTextFont(72);
   	chi2text -> SetFillColor(kWhite);
    chi2text->AddText(Form("Chi2 avg val = %4.5f#pm%.5f", ChiSqs_avg,ChiSqs_stddev));   
    avgchi2->SetLineColor(kCyan);
    avgchi2->SetLineWidth(2);
    optchi2->SetLineColor(kMagenta);
    optchi2->SetLineWidth(2);
    stddevchi2_abv->SetLineColor(kOrange-3);
    stddevchi2_blw->SetLineColor(kOrange-3);
    stddevchi2_blw->SetLineWidth(2);
    avgchi2->SetLineWidth(2);           
    avgchi2->Draw("same");
    optchi2->Draw("same");
    stddevchi2_abv->Draw("same");
    stddevchi2_blw->Draw("same");
    chi2text->Draw(); 
    auto legend_chi2 = new TLegend();
    legend_chi2->AddEntry(avgchi2, "Average Value", "L");
    legend_chi2->AddEntry(optchi2, "Optimal Value", "L");
    legend_chi2->AddEntry(stddevchi2_abv,"Std Dev","L");
    legend_chi2->Draw();
    outputFile->cd();
    grchi2.Write("Chi2");
    canvas2.SaveAs("Chi2.root");    

    grN.SetMarkerSize(0.8);
    grN.SetMarkerColor(kGreen);
    grN.SetLineColor(kGreen);
    TLine *avgevts = new TLine(cry_start,Events_avg,cry_end,Events_avg);
    TLine *stddevevts_abv = new TLine(cry_start,Events_avg+Events_stddev,cry_end,Events_avg+Events_stddev);        
    TLine *stddevevts_blw = new TLine(cry_start,Events_avg-Events_stddev,cry_end,Events_avg-Events_stddev);        
    TCanvas canvas3("canvas", "Scatter Plot", 800, 600);
    grN.Draw("AP"); // A for axis, P for points
    avgevts->SetLineColor(kCyan);
    avgevts->SetLineWidth(2);
    stddevevts_abv->SetLineColor(kOrange-3);
    stddevevts_blw->SetLineColor(kOrange-3);        
    stddevevts_abv->SetLineWidth(2);
    stddevevts_blw->SetLineColor(kOrange-3);
    stddevevts_blw->SetLineWidth(2);        
    avgevts->Draw("same");
    stddevevts_abv->Draw("same");
    stddevevts_blw->Draw("same");
    auto legend_evts = new TLegend();
    legend_evts->AddEntry(avgevts, "Average Value", "L");
    legend_evts->AddEntry(stddevevts_abv,"Std Dev","L");
    legend_evts->Draw();
    outputFile->cd();
    grN.Write("nEvents");
    canvas3.SaveAs("nEvents.root");    
    
// ---- Build SiPM pair maps by crystalNo using crystalNo parity (even/odd logic) ----
std::map<int, std::map<int, double>> sipm_peaks_by_crystal;

for (size_t i = 0; i < crystalNos.size(); ++i) {
    int sipm = static_cast<int>(crystalNos[i]);
    int crystal = sipm / 2;  // assumes sipm 2n and 2n+1 belong to same crystal
    int sipm_in_pair = sipm % 2; // 0 or 1

    sipm_peaks_by_crystal[crystal][sipm_in_pair] = Peaks[i];
}

// Prepare data for plots
std::vector<Double_t> sipm0_peaks, sipm1_peaks, pair_diffs, crystal_ids;

for (const auto& [crystal, sipm_map] : sipm_peaks_by_crystal) {
    if (sipm_map.size() != 2) {
        std::cerr << "Warning: crystalNo " << crystal << " does not have exactly 2 SiPMs (has " << sipm_map.size() << ")\n";
        continue;
    }

    if (sipm_map.count(0) == 0 || sipm_map.count(1) == 0) {
        std::cerr << "Warning: crystalNo " << crystal << " missing SiPM " << (sipm_map.count(0) ? 1 : 0) << "\n";
        continue;
    }

    double p0 = sipm_map.at(0);
    double p1 = sipm_map.at(1);

    sipm0_peaks.push_back(p0);
    sipm1_peaks.push_back(p1);
    pair_diffs.push_back(p0 - p1);
    crystal_ids.push_back(crystal);
}

if (sipm0_peaks.empty() || sipm1_peaks.empty()) {
    std::cerr << "Error: No valid crystal pairs found for plotting.\n";
    return;
}

// --- Scatter plot: SiPM 0 vs SiPM 1 peak ---
TGraph grScatter(sipm0_peaks.size(), &sipm0_peaks[0], &sipm1_peaks[0]);
grScatter.SetTitle("SiPM Pair Peaks;SiPM 0 Peak;SiPM 1 Peak");
grScatter.SetMarkerStyle(20);
grScatter.SetMarkerColor(kBlue);
grScatter.SetMarkerSize(0.9);

TCanvas canvasScatter("canvasScatter", "SiPM Scatter", 800, 600);
grScatter.Draw("AP");

double min_peak = std::min(*std::min_element(sipm0_peaks.begin(), sipm0_peaks.end()),
                           *std::min_element(sipm1_peaks.begin(), sipm1_peaks.end()));
double max_peak = std::max(*std::max_element(sipm0_peaks.begin(), sipm0_peaks.end()),
                           *std::max_element(sipm1_peaks.begin(), sipm1_peaks.end()));

TLine identity(min_peak, min_peak, max_peak, max_peak);
identity.SetLineColor(kRed);
identity.SetLineStyle(2);
identity.Draw("same");

outputFile->cd();
grScatter.Write("SiPM_Peak_Scatter");
canvasScatter.SaveAs("SiPM_Peak_Scatter.root");

// --- Difference plot: Peak0 - Peak1 vs crystal number ---
TGraph grDiff(crystal_ids.size(), &crystal_ids[0], &pair_diffs[0]);
grDiff.SetTitle("SiPM Peak Difference (0 - 1);Crystal Number;Peak Difference");
grDiff.SetMarkerStyle(21);
grDiff.SetMarkerColor(kGreen + 2);
grDiff.SetMarkerSize(0.9);

TCanvas canvasDiff("canvasDiff", "SiPM Peak Difference", 800, 600);
grDiff.Draw("AP");

TLine zeroLine(*std::min_element(crystal_ids.begin(), crystal_ids.end()), 0,
               *std::max_element(crystal_ids.begin(), crystal_ids.end()), 0);
zeroLine.SetLineColor(kRed);
zeroLine.SetLineStyle(2);
zeroLine.Draw("same");

outputFile->cd();
grDiff.Write("SiPM_Peak_Difference");
canvasDiff.SaveAs("SiPM_Peak_Difference.root");

// --- Projection: Histogram of Peak Differences and Gaussian Fit ---
int nbins = 50;
double diff_min = *std::min_element(pair_diffs.begin(), pair_diffs.end());
double diff_max = *std::max_element(pair_diffs.begin(), pair_diffs.end());

// Add some margin to bin edges
double margin = 0.05 * (diff_max - diff_min);
diff_min -= margin;
diff_max += margin;

// Create histogram
TH1D* hDiff = new TH1D("hPeakDiff", "Histogram of SiPM Peak Differences;Peak0 - Peak1;Counts", nbins, diff_min, diff_max);

// Fill histogram
for (double diff : pair_diffs) {
    hDiff->Fill(diff);
}

// Fit with Gaussian
TF1* gausFit = new TF1("gausFit", "gaus", diff_min, diff_max);
hDiff->Fit(gausFit, "Q");  // "Q" suppresses fit printing

// Draw
TCanvas canvasHist("canvasHist", "Peak Difference Histogram", 800, 600);
hDiff->SetLineColor(kBlue + 1);
hDiff->SetLineWidth(2);
hDiff->Draw();

gausFit->SetLineColor(kRed);
gausFit->SetLineWidth(2);
gausFit->Draw("same");

// Print fit parameters
double mean = gausFit->GetParameter(1);
double sigma = gausFit->GetParameter(2);
std::cout << "Gaussian Fit Results: Mean = " << mean << ", Sigma = " << sigma << std::endl;

// Save to file
outputFile->cd();
hDiff->Write("PeakDiff_Histogram");
gausFit->Write("PeakDiff_GaussianFit");
canvasHist.SaveAs("PeakDiff_Histogram.root");


}
