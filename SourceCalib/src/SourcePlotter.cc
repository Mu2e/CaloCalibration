#include "CaloCalibration/SourceCalib/inc/SourcePlotter.hh"
#include <chrono>
#include <numeric>
#include "TF1.h"

using namespace std::chrono;
using namespace CaloSourceCalib;

/* function to make global plots of the fit outputs*/
void SourcePlotter::ParamPlots(TTree* inputTree, TFile *inputFile, TFile *outputFile,int cry_start, int cry_end) {//add param as plot range
    float crystalNo, Peak, ChiSq, PeakErr,h_means,h_stds,frFull,frFrst,frScnd,convgstatus,firstPeak,secondPeak;
    int nEvents;    
    inputTree-> SetBranchAddress("crystalNo", &crystalNo);
    inputTree-> SetBranchAddress("Peak", &Peak);
    inputTree-> SetBranchAddress("PeakErr", &PeakErr);
    inputTree-> SetBranchAddress("ChiSq", &ChiSq);
    inputTree-> SetBranchAddress("nEvents", &nEvents); 
    inputTree-> SetBranchAddress("h_means", &h_means); 
    inputTree-> SetBranchAddress("h_stds", &h_stds); 
    inputTree-> SetBranchAddress("frFull", &frFull);
    inputTree-> SetBranchAddress("frFrst", &frFrst);
    inputTree-> SetBranchAddress("frScnd", &frScnd); 
    inputTree-> SetBranchAddress("convgstatus", &convgstatus); 
    inputTree-> SetBranchAddress("1stPeak", &firstPeak); 
    inputTree-> SetBranchAddress("2ndPeak", &secondPeak); 
    std::vector<Double_t> crystalNos, Peaks,ADCPeaks, PeakErrs,ADCPeakErrs, ChiSqs, EventsErr, Events,h_means_vec,h_stds_vec,frFull_vec,frFrst_vec,frScnd_vec,convgstatus_vec,firstPeak_vec,secondPeak_vec;
     
    // Fill the vectors with data from the TTree
    Long64_t nEntries = inputTree->GetEntries();
    for (Long64_t entry = 0; entry < nEntries; ++entry) {
        inputTree->GetEntry(entry);
        crystalNos.push_back(crystalNo);
        Peaks.push_back(1/(Peak/6.13));
        ADCPeaks.push_back(Peak);
        ADCPeakErrs.push_back(PeakErr);
        PeakErrs.push_back((1/(Peak/6.13))*(PeakErr/Peak));
        ChiSqs.push_back(ChiSq);
        Events.push_back(nEvents);
        EventsErr.push_back(sqrt(nEvents));       
        h_means_vec.push_back(h_means);       
        h_stds_vec.push_back(h_stds); 
        frFull_vec.push_back(frFull);
        frFrst_vec.push_back(frFrst);
        frScnd_vec.push_back(frScnd);
        convgstatus_vec.push_back(convgstatus);
        firstPeak_vec.push_back(firstPeak);
        secondPeak_vec.push_back(secondPeak);
    }
//----create value of avg and std dev for tlines for plots---
		double Peaks_sum = std::accumulate(Peaks.begin(), Peaks.end(),0.0);
 		double Peaks_avg = Peaks_sum / Peaks.size();
 		double Peaks_stddev = TMath::StdDev(Peaks.begin(),Peaks.end());
 		double stddev_prcnt = (Peaks_stddev/Peaks_avg)*100;
 		std::cout<<" peak 2 std dev "<<Peaks_stddev*2<<std::endl;
 		std::cout<<" peak std dev "<<Peaks_stddev<<std::endl;
 		std::cout<<" ------------ "<<std::endl; 
//----create loop to count and print how many crys/sipms have pear error larger than 2 stdev----		
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
//----resume creating values of avg and std dev for tlines for plots-----
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
    TGraph grpeakerrs(crystalNos.size(), &crystalNos[0], &PeakErrs[0]);
    TGraph grconvstatus(convgstatus_vec.size(), &convgstatus_vec[0], &PeakErrs[0]);

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
    canvas.Update();  // Must update to get access to Y axis
    grpeaks.GetYaxis()->SetRangeUser(0.05, 0.08);  // Y-axis limits  
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
//-----ADC Peaks graph that represent the full peak in ADC--------
		double ADCPeaks_sum = std::accumulate(ADCPeaks.begin(), ADCPeaks.end(),0.0);
 		double ADCPeaks_avg = ADCPeaks_sum / ADCPeaks.size();
 		double ADCPeaks_stddev = TMath::StdDev(ADCPeaks.begin(),ADCPeaks.end());
 		double ADCstddev_prcnt = (ADCPeaks_stddev/ADCPeaks_avg)*100;
 		std::cout<<" ADC peak 2 std dev "<<ADCPeaks_stddev*2<<std::endl;
 		std::cout<<" ADC peak std dev "<<ADCPeaks_stddev<<std::endl;
 		std::cout<<" ------------ "<<std::endl; 		
 		int ADCcountHighPeakErrs = 0;
for (size_t i = 0; i < ADCPeakErrs.size(); ++i) {
    if (ADCPeakErrs[i] > 2 * Peaks_stddev) {
        ++ADCcountHighPeakErrs;
        std::cout << "Crystal " << crystalNos[i] 
                  << " has ADCPeakErr = " << ADCPeakErrs[i] 
                  << " > 2*ADCPeaks_stddev (" << 2 * ADCPeaks_stddev << ")\n";
    }
}
std::cout << "Number of ADC PeakErrs > 2 stdev " << ADCcountHighPeakErrs << std::endl;
 		std::cout<<" ------------ "<<std::endl;
 			
    // Create a TGraph for the scatter plot
    TGraphErrors ADCgrpeaks(crystalNos.size(), &crystalNos[0], &ADCPeaks[0], nullptr, &ADCPeakErrs[0]);
    TGraph ADCgrpeakerrs(crystalNos.size(), &crystalNos[0], &ADCPeakErrs[0]);
    
    // Customize the plot
    ADCgrpeaks.SetTitle("Cry Number vs ADC;Crystal Number;ADC");
    ADCgrpeaks.SetMarkerStyle(20);
    ADCgrpeaks.SetMarkerSize(0.8);
    ADCgrpeaks.SetMarkerColor(kBlue);
    ADCgrpeaks.SetLineColor(kBlue);
    TLine *ADCavgpeak = new TLine(cry_start,ADCPeaks_avg,cry_end,ADCPeaks_avg);
    TLine *ADCnominalpeak = new TLine(cry_start,98,cry_end,98); 
    TLine *ADCstddevpeak_abv = new TLine(cry_start,ADCPeaks_avg+2*ADCPeaks_stddev,cry_end,ADCPeaks_avg+2*ADCPeaks_stddev);    
    TLine *ADCstddevpeak_blw = new TLine(cry_start,ADCPeaks_avg-2*ADCPeaks_stddev,cry_end,ADCPeaks_avg-2*ADCPeaks_stddev);		
    TCanvas ADCcanvas("ADCcanvas", "Scatter Plot", 800, 600);
    ADCgrpeaks.Draw("AP"); // A for axis, P for points  
    ADCcanvas.Update();  // Must update to get access to Y axis
    ADCgrpeaks.GetYaxis()->SetRangeUser(122, 77);  // Y-axis limits  
    TPaveText *ADCavg = new TPaveText(0.15, 0.75, 0.35, 0.65, "brNDC");
   	ADCavg -> SetFillStyle(0);
   	ADCavg -> SetBorderSize(0);
   	ADCavg -> SetTextSize(0.03);
   	ADCavg -> SetTextColor(kBlack);
   	ADCavg -> SetTextFont(72);
   	ADCavg -> SetFillColor(kWhite);
    ADCavg->AddText(Form("ADC Peaks avg val = %4.8f#pm%.8f", ADCPeaks_avg,ADCPeaks_stddev));
		ADCavg->AddText(Form("Std Dev = %.2f%%", ADCstddev_prcnt));   
    ADCavgpeak->SetLineColor(kCyan);
    ADCavgpeak->SetLineWidth(2);
    ADCnominalpeak->SetLineColor(kMagenta);
    ADCnominalpeak->SetLineWidth(2);
    ADCstddevpeak_abv->SetLineColor(kOrange-3);
    ADCstddevpeak_blw->SetLineColor(kOrange);        
    ADCstddevpeak_abv->SetLineWidth(2);    
    ADCstddevpeak_blw->SetLineColor(kOrange);
    ADCstddevpeak_blw->SetLineWidth(2);        
    ADCavgpeak->Draw("same");
    ADCnominalpeak->Draw("same");
    ADCstddevpeak_abv->Draw("same");
    ADCstddevpeak_blw->Draw("same");
    ADCavg->Draw();
    auto ADClegend = new TLegend();
    ADClegend->AddEntry(ADCavgpeak, "Average Value", "L");
    ADClegend->AddEntry(ADCnominalpeak,"Nominal Value","L");
    ADClegend->AddEntry(ADCstddevpeak_abv,"Std Dev","L");
    ADClegend->AddEntry(ADCstddevpeak_blw,"Std Dev","L");
     //legend->SetTextSize(0.4);
    ADClegend->Draw();
    outputFile->cd();
    ADCgrpeaks.Write("ADCPeaks");
    ADCcanvas.SaveAs("ADCPeaks.root");
          
    //TCanvas canvas("canvas", "Scatter Plot", 800, 600);
    ADCgrpeaks.Draw("APsame"); // A for axis, P for points
    outputFile->cd();
    ADCgrpeaks.Write("ADCPeaks");
  
//----Chi2 plot-----    
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
    
//----PeakErrors plot that represents the errors from the fit in MeV/ADC-----
    grpeakerrs.SetTitle("Cry Number vs Peak Errors;Crystal Number;Peak Errors");
	grpeakerrs.SetMarkerStyle(20);   // full circle (visible marker)
	grpeakerrs.SetMarkerSize(7.0);   // make points bigger
    grpeakerrs.SetMarkerColor(kGreen);
    grpeakerrs.SetLineColor(kGreen);        
    TCanvas canvas4("canvas", "Scatter Plot", 800, 600);
    grpeakerrs.Draw("AP"); // A for axis, P for points
    outputFile->cd();
    grpeakerrs.Write("peakerrs");
    canvas4.SaveAs("peakerrs.root");
    
    grpeakerrs.Draw("AP"); // A for axis, P for points
    outputFile->cd();
    grpeakerrs.Write("PeaksErrors");

//-----convergence status plots that come from the nll fit function-------- 
    grconvstatus.SetTitle("Peak Error vs Convg Status ;Convergence Status; PeakError");
    grconvstatus.SetMarkerStyle(20);
    grconvstatus.SetMarkerSize(0.8);
    grconvstatus.SetMarkerColor(kRed);
    grconvstatus.SetLineColor(kRed);
    TCanvas canvas_conv("canvas", "Scatter Plot", 800, 600);
    grconvstatus.Draw("AP"); // A for axis, P for points
    outputFile->cd();
    grconvstatus.Write("Convstatus");
    canvas_conv.SaveAs("convstatusvspeakerrs.root");  
   
// ---- Scatter plot of calib consts of sipm pairs laid agaisnt each other (i.e even on x axis odd on y axis)-----
//Build SiPM pair maps by crystalNo using crystalNo parity (even/odd logic) 
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

// --- Set axis limits ---
canvasScatter.Update();
grScatter.GetXaxis()->SetLimits(min_peak, max_peak); // set X range
grScatter.GetYaxis()->SetRangeUser(min_peak, max_peak); // set Y range


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


//----Plots for means and std devs of the histrogrammed data----
std::map<int, std::map<int, double>> sipm_hmeans_by_crystal;
std::map<int, std::map<int, double>> sipm_hstds_by_crystal;

for (size_t i = 0; i < crystalNos.size(); ++i) {
    int sipm = static_cast<int>(crystalNos[i]);
    int crystal = sipm / 2;
    int sipm_in_pair = sipm % 2;

    sipm_hmeans_by_crystal[crystal][sipm_in_pair] = h_means_vec[i];
    sipm_hstds_by_crystal[crystal][sipm_in_pair] = h_stds_vec[i];
}
std::vector<Double_t> sipm0_hmeans, sipm1_hmeans;
std::vector<Double_t> sipm0_hstds, sipm1_hstds;

for (const auto& [crystal, map_means] : sipm_hmeans_by_crystal) {
    if (map_means.size() != 2) continue;
    if (map_means.count(0) == 0 || map_means.count(1) == 0) continue;
    sipm0_hmeans.push_back(map_means.at(0));
    sipm1_hmeans.push_back(map_means.at(1));
    crystal_ids.push_back(crystal);
}

for (const auto& [crystal, map_stds] : sipm_hstds_by_crystal) {
    if (map_stds.size() != 2) continue;
    if (map_stds.count(0) == 0 || map_stds.count(1) == 0) continue;
    sipm0_hstds.push_back(map_stds.at(0));
    sipm1_hstds.push_back(map_stds.at(1));
}
// For h_means
TGraph grHMeans(sipm0_hmeans.size(), &sipm0_hmeans[0], &sipm1_hmeans[0]);
grHMeans.SetTitle("SiPM Pair h_means;SiPM 0 h_mean;SiPM 1 h_mean");
grHMeans.SetMarkerStyle(20);
grHMeans.SetMarkerColor(kBlue);
grHMeans.SetMarkerSize(0.9);

TCanvas canvasHMeans("canvasHMeans", "h_means Scatter", 800, 600);
grHMeans.Draw("AP");

double min_hmean = std::min(*std::min_element(sipm0_hmeans.begin(), sipm0_hmeans.end()),
                           *std::min_element(sipm1_hmeans.begin(), sipm1_hmeans.end()));
double max_hmean = std::max(*std::max_element(sipm0_hmeans.begin(), sipm0_hmeans.end()),
                           *std::max_element(sipm1_hmeans.begin(), sipm1_hmeans.end()));

TLine identityHMean(min_hmean, min_hmean, max_hmean, max_hmean);
identityHMean.SetLineColor(kRed);
identityHMean.SetLineStyle(2);
identityHMean.Draw("same");

outputFile->cd();
grHMeans.Write("SiPM_hmeans_Scatter");
canvasHMeans.SaveAs("SiPM_hmeans_Scatter.root");

// for h_stds
if (!sipm0_hstds.empty() && !sipm1_hstds.empty()) {
	TGraph grHStds(sipm0_hstds.size(), &sipm0_hstds[0], &sipm1_hstds[0]);
	grHStds.SetTitle("SiPM Pair h_stds;SiPM 0 h_std;SiPM 1 h_std");
	grHStds.SetMarkerStyle(20);
	grHStds.SetMarkerColor(kBlue);
	grHStds.SetMarkerSize(0.9);

	TCanvas canvasHStds("canvasHStds", "h_stds Scatter", 800, 600);
	grHStds.Draw("AP");

	double min_hstd = std::min(*std::min_element(sipm0_hstds.begin(), sipm0_hstds.end()),
                           *std::min_element(sipm1_hstds.begin(), sipm1_hstds.end()));
	double max_hstd = std::max(*std::max_element(sipm0_hstds.begin(), sipm0_hstds.end()),
                           *std::max_element(sipm1_hstds.begin(), sipm1_hstds.end()));

	TLine identityHStd(min_hstd, min_hstd, max_hstd, max_hstd);
	identityHStd.SetLineColor(kRed);
	identityHStd.SetLineStyle(2);
	identityHStd.Draw("same");

	outputFile->cd();
	grHStds.Write("SiPM_hstds_Scatter");
	canvasHStds.SaveAs("SiPM_hstds_Scatter.root");
}else {
    std::cout << "?? No h_stds data available for plotting!" << std::endl;
}

// ---- Scatter plot of full, first & scnd fractions of sipm pairs laid agaisnt each other (i.e even on x axis odd on y axis)-----
//Build SiPM pair maps for frFull, frFrst, frScnd 
std::map<int, std::map<int, double>> sipm_frFull_by_crystal;
std::map<int, std::map<int, double>> sipm_frFrst_by_crystal;
std::map<int, std::map<int, double>> sipm_frScnd_by_crystal;

for (size_t i = 0; i < crystalNos.size(); ++i) {
    int sipm = static_cast<int>(crystalNos[i]);
    int crystal = sipm / 2;       // assumes sipm 2n and 2n+1 belong to same crystal
    int sipm_in_pair = sipm % 2;  // 0 or 1

    sipm_frFull_by_crystal[crystal][sipm_in_pair] = frFull_vec[i];
    sipm_frFrst_by_crystal[crystal][sipm_in_pair] = frFrst_vec[i];
    sipm_frScnd_by_crystal[crystal][sipm_in_pair] = frScnd_vec[i];
}

// Prepare data vectors
std::vector<Double_t> sipm0_frFull, sipm1_frFull;
std::vector<Double_t> sipm0_frFrst, sipm1_frFrst;
std::vector<Double_t> sipm0_frScnd, sipm1_frScnd;

for (const auto& kv : sipm_frFull_by_crystal) {
    const auto& m = kv.second;
    if (m.count(0) && m.count(1)) {  // robust check
        sipm0_frFull.push_back(m.at(0));
        sipm1_frFull.push_back(m.at(1));
    }
}
for (const auto& kv : sipm_frFrst_by_crystal) {
    const auto& m = kv.second;
    if (m.count(0) && m.count(1)) {
        sipm0_frFrst.push_back(m.at(0));
        sipm1_frFrst.push_back(m.at(1));
    }
}
for (const auto& kv : sipm_frScnd_by_crystal) {
    const auto& m = kv.second;
    if (m.count(0) && m.count(1)) {
        sipm0_frScnd.push_back(m.at(0));
        sipm1_frScnd.push_back(m.at(1));
    }
}

    // Build graphs (graphs can be empty; that?s OK)
    TGraph grFrFull(sipm0_frFull.size(), sipm0_frFull.data(), sipm1_frFull.data());
    grFrFull.SetMarkerStyle(20);
    grFrFull.SetMarkerColor(kBlue+1);
    grFrFull.SetMarkerSize(1.0);

    TGraph grFrFrst(sipm0_frFrst.size(), sipm0_frFrst.data(), sipm1_frFrst.data());
    grFrFrst.SetMarkerStyle(24); // open circle for visibility
    grFrFrst.SetMarkerColor(kGreen+2);
    grFrFrst.SetMarkerSize(1.0);

    TGraph grFrScnd(sipm0_frScnd.size(), sipm0_frScnd.data(), sipm1_frScnd.data());
    grFrScnd.SetMarkerStyle(22); // triangle
    grFrScnd.SetMarkerColor(kMagenta+1);
    grFrScnd.SetMarkerSize(1.0);

    // Compute combined axis limits across ALL datasets (skip NaN/Inf)
    auto upd = [](double v, double& mn, double& mx) {
        if (std::isfinite(v)) { mn = std::min(mn, v); mx = std::max(mx, v); }
    };
    double xmin = +std::numeric_limits<double>::infinity();
    double xmax = -std::numeric_limits<double>::infinity();
    double ymin = +std::numeric_limits<double>::infinity();
    double ymax = -std::numeric_limits<double>::infinity();

    for (size_t i=0;i<sipm0_frFull.size(); ++i){ upd(sipm0_frFull[i], xmin, xmax); upd(sipm1_frFull[i], ymin, ymax); }
    for (size_t i=0;i<sipm0_frFrst.size(); ++i){ upd(sipm0_frFrst[i], xmin, xmax); upd(sipm1_frFrst[i], ymin, ymax); }
    for (size_t i=0;i<sipm0_frScnd.size(); ++i){ upd(sipm0_frScnd[i], xmin, xmax); upd(sipm1_frScnd[i], ymin, ymax); }

    // Fallback if one axis was totally empty
    if (!std::isfinite(xmin) || !std::isfinite(xmax)) { xmin = 0.; xmax = 1.; }
    if (!std::isfinite(ymin) || !std::isfinite(ymax)) { ymin = 0.; ymax = 1.; }

    // Add a small margin
    const double dx = (xmax - xmin) > 0 ? 0.05*(xmax - xmin) : 0.05;
    const double dy = (ymax - ymin) > 0 ? 0.05*(ymax - ymin) : 0.05;
    xmin -= dx; xmax += dx; ymin -= dy; ymax += dy;

    // Draw on a fixed frame so the later graphs cannot be clipped
    TCanvas canvasFrac("canvasFrac", "SiPM Fractions Scatter", 800, 600);
    TH1F* frame = gPad->DrawFrame(xmin, ymin, xmax, ymax);
    frame->SetTitle("SiPM Pair Fractions;SiPM 0;SiPM 1");

    // Identity line
    double lmin = std::min(xmin, ymin);
    double lmax = std::max(xmax, ymax);
    TLine identityFrac(lmin, lmin, lmax, lmax);
    identityFrac.SetLineColor(kRed);
    identityFrac.SetLineStyle(2);

    // Draw graphs
    grFrFull.Draw("P SAME");
    grFrFrst.Draw("P SAME");
    grFrScnd.Draw("P SAME");
    identityFrac.Draw("SAME");

    // Legend (only add entries that have points)
    TLegend frlegend(0.15, 0.75, 0.45, 0.9);
    if (grFrFull.GetN() > 0)  frlegend.AddEntry(&grFrFull, "frFull", "p");
    if (grFrFrst.GetN() > 0)  frlegend.AddEntry(&grFrFrst, "frFrst", "p");
    if (grFrScnd.GetN() > 0)  frlegend.AddEntry(&grFrScnd, "frScnd", "p");
    frlegend.Draw();

    outputFile->cd();
    grFrFull.Write("SiPM_FrFull_Scatter");
    grFrFrst.Write("SiPM_FrFrst_Scatter");
    grFrScnd.Write("SiPM_FrScnd_Scatter");
    canvasFrac.SaveAs("SiPM_Fractions_Scatter.root");

// --- Overlay plot: ADC Peak, 1stPeak, 2ndPeak vs Crystal Number ---
TGraph grADC(crystalNos.size(), &crystalNos[0], &ADCPeaks[0]);
TGraph gr1st(crystalNos.size(), &crystalNos[0], &firstPeak_vec[0]);
TGraph gr2nd(crystalNos.size(), &crystalNos[0], &secondPeak_vec[0]);

// Style
grADC.SetMarkerStyle(20); // filled circle
grADC.SetMarkerColor(kBlue);
grADC.SetLineColor(kBlue);
grADC.SetMarkerSize(0.9);

gr1st.SetMarkerStyle(21); // square
gr1st.SetMarkerColor(kRed);
gr1st.SetLineColor(kRed);
gr1st.SetMarkerSize(0.9);

gr2nd.SetMarkerStyle(22); // triangle
gr2nd.SetMarkerColor(kGreen+2);
gr2nd.SetLineColor(kGreen+2);
gr2nd.SetMarkerSize(0.9);

// Canvas
TCanvas cADC1st2nd("cADC1st2nd", "ADC, 1st, and 2nd Peaks", 800, 600);

// Compute combined y-range 
double ymin = std::numeric_limits<double>::infinity();
double ymax = -std::numeric_limits<double>::infinity();
auto updateMinMax = [&](const TGraph& g) {
    for (int i = 0; i < g.GetN(); ++i) {
        double x, y;
        g.GetPoint(i, x, y);
        if (std::isfinite(y)) {
            ymin = std::min(ymin, y);
            ymax = std::max(ymax, y);
        }
    }
};
updateMinMax(grADC);
updateMinMax(gr1st);
updateMinMax(gr2nd);

// Add margin
double dy = (ymax - ymin) * 0.05;
if (dy <= 0) dy = 1.0; // fallback
ymin -= dy;
ymax += dy;

// Draw empty frame first with full range
TH1F* frame = cADC1st2nd.DrawFrame(
    crystalNos.front(), ymin,
    crystalNos.back(), ymax
);
frame->SetTitle("Peak Locations;Crystal Number;ADC / Peak Value");

// Draw graphs
grADC.Draw("P SAME");
gr1st.Draw("P SAME");
gr2nd.Draw("P SAME");

// Legend
auto legADC = new TLegend(0.15, 0.75, 0.45, 0.9);
legADC->AddEntry(&grADC, "ADC Peak", "p");
legADC->AddEntry(&gr1st, "1st Peak", "p");
legADC->AddEntry(&gr2nd, "2nd Peak", "p");
legADC->Draw();

// Save
outputFile->cd();
grADC.Write("ADC_vs_Crystal");
gr1st.Write("FirstPeak_vs_Crystal");
gr2nd.Write("SecondPeak_vs_Crystal");
cADC1st2nd.SaveAs("ADC_1st_2nd_Peaks.root");

}
