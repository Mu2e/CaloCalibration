#include "CaloCalibration/SourceCalib/inc/SourcePlotter.hh"
#include <chrono>
#include <numeric>
#include "TF1.h"
#include "TH2.h"
#include "TProfile.h"
using namespace std::chrono;
using namespace CaloSourceCalib;

/* function to make global plots of the fit outputs*/
void SourcePlotter::ParamPlots(TTree* inputTree, TFile *inputFile, TFile *outputFile,int cry_start, int cry_end) {//add param as plot range
    float crystalNo, Peak, ChiSq, PeakErrHigh,PeakErrLo,h_means,h_stds,frFull,frFrst,frScnd,firstPeak,secondPeak,Width,WidthErrHigh,WidthErrLo,unreducedchi2;
    int nEvents,convgstatus,ndof;    
    inputTree-> SetBranchAddress("crystalNo", &crystalNo);
    inputTree-> SetBranchAddress("Peak", &Peak);
    inputTree-> SetBranchAddress("PeakErrLo", &PeakErrLo);
    inputTree-> SetBranchAddress("PeakErrHigh", &PeakErrHigh);
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
    inputTree-> SetBranchAddress("Width", &Width); 
    inputTree-> SetBranchAddress("WidthErrHigh", &WidthErrHigh); 
    inputTree-> SetBranchAddress("WidthErrLo", &WidthErrLo); 
    inputTree-> SetBranchAddress("unreducedchi2", &unreducedchi2); 
    inputTree-> SetBranchAddress("ndof", &ndof); 
    

    std::vector<Double_t> crystalNos, Peaks,ADCPeaks, PeakErrLos,PeakErrHighs,ADCPeakErrHighs,ADCPeakErrLos, ChiSqs, EventsErr, Events,h_means_vec,h_stds_vec,frFull_vec,frFrst_vec,frScnd_vec,convgstatus_vec,
firstPeak_vec,secondPeak_vec,Widths,WidthErrHighs,WidthErrLos,unreducedchi2s,ndofs;
     
    // Fill the vectors with data from the TTree
    Long64_t nEntries = inputTree->GetEntries();
    for (Long64_t entry = 0; entry < nEntries; ++entry) {
        inputTree->GetEntry(entry);
        crystalNos.push_back(crystalNo);
        Peaks.push_back(1/(Peak/6.13));
        ADCPeaks.push_back(Peak);
        ADCPeakErrHighs.push_back(fabs(PeakErrHigh));
        ADCPeakErrLos.push_back(fabs(PeakErrLo));
        PeakErrLos.push_back(fabs((6.13/Peak)*(PeakErrLo/Peak)));
        PeakErrHighs.push_back(fabs((6.13/Peak)*(PeakErrHigh/Peak)));
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
        Widths.push_back(Width);
        WidthErrHighs.push_back(WidthErrHigh);
        WidthErrLos.push_back(WidthErrLo);
        unreducedchi2s.push_back(unreducedchi2);
        ndofs.push_back(ndof);

    }
//----create value of avg and std dev for tlines for plots---
    double Peaks_sum = std::accumulate(Peaks.begin(), Peaks.end(), 0.0);
    double Peaks_avg = Peaks_sum / Peaks.size();
    double Peaks_stddev = TMath::StdDev(Peaks.begin(), Peaks.end());
    double stddev_prcnt = (Peaks_stddev / Peaks_avg) * 100;
    
    double peakerrhi_sum = std::accumulate(PeakErrHighs.begin(), PeakErrHighs.end(), 0.0);
    double peakerrlo_sum = std::accumulate(PeakErrLos.begin(), PeakErrLos.end(), 0.0);
    double avgpeakerr_hi = peakerrhi_sum / PeakErrHighs.size();
    double avgpeakerr_lo = peakerrlo_sum / PeakErrLos.size();
    
    std::cout << " peak 2 std dev " << Peaks_stddev * 2 << std::endl;
    std::cout << " peak std dev " << Peaks_stddev << std::endl;
    std::cout << " ------------ " << std::endl;

//----create loop to count and print how many crys/sipms have pear error larger than 2 stdev----		
    int countHighPeakErrHighs = 0;
    for (size_t i = 0; i < PeakErrHighs.size(); ++i) {
        if (PeakErrHighs[i] > 2 * Peaks_stddev) {
            ++countHighPeakErrHighs;
            std::cout << "Crystal " << crystalNos[i] 
                      << " has PeakErrHigh = " << PeakErrHighs[i] 
                      << " > 2*Peaks_stddev (" << 2 * Peaks_stddev << ")\n";
        }
    }
    int countHighPeakErrLos = 0;
    for (size_t i = 0; i < PeakErrLos.size(); ++i) {
        if (PeakErrLos[i] > 2 * Peaks_stddev) {
            ++countHighPeakErrLos;
            std::cout << "Crystal " << crystalNos[i] 
                      << " has PeakErrLos = " << PeakErrLos[i] 
                      << " > 2*Peaks_stddev (" << 2 * Peaks_stddev << ")\n";
        }
    }
    std::cout << "Number of PeakErrHighs > 2 stdev " << countHighPeakErrHighs << std::endl;
    std::cout << " ------------ " << std::endl;
    std::cout << "Number of PeakErrLos > 2 stdev " << countHighPeakErrLos << std::endl;
    std::cout << " ------------ " << std::endl;

//----resume creating values of avg and std dev for tlines for plots-----
    double ChiSqs_sum = std::accumulate(ChiSqs.begin(), ChiSqs.end(), 0.0);
    double ChiSqs_avg = ChiSqs_sum / ChiSqs.size();
    double ChiSqs_stddev = TMath::StdDev(ChiSqs.begin(), ChiSqs.end());
    std::cout << " std dev of chi2 " << ChiSqs_stddev << std::endl;
    std::cout << " 2 std dev of chi2 " << 2 * ChiSqs_stddev << std::endl;
    
    double Events_sum = std::accumulate(Events.begin(), Events.end(), 0.0);
    double Events_avg = Events_sum / Events.size();
    double Events_stddev = TMath::StdDev(Events.begin(), Events.end());
    std::cout << " std dev of events " << Events_stddev << std::endl;

    // =======================================================
    // SPLIT DATA INTO SUBSETS
    // Logic: The 8 specific IDs are LYSO. The rest are CsI.
    // =======================================================
    // Using long and lround to ensure ID matching works even if crystalNos are floats
    std::set<long> lyso_ids = {1164, 1165, 1220, 1221, 1218, 1219, 1274, 1275};
    
    // Vectors for LYSO (The 8 specific IDs)
    std::vector<double> crys_L, p_L, el_L, eh_L; 
    // Vectors for CsI (The rest)
    std::vector<double> crys_C, p_C, el_C, eh_C; 

    for(size_t i=0; i<crystalNos.size(); ++i) {
        // Use rounding to ensure 1120.00 matches 1120
        long id = std::lround(crystalNos[i]);
        
        if(lyso_ids.count(id)) {
            // It IS one of the 8 -> LYSO
            crys_L.push_back(crystalNos[i]);
            p_L.push_back(Peaks[i]);
            el_L.push_back(PeakErrLos[i]);
            eh_L.push_back(PeakErrHighs[i]);
        } else {
            // It is NOT one of the 8 -> CsI
            crys_C.push_back(crystalNos[i]);
            p_C.push_back(Peaks[i]);
            el_C.push_back(PeakErrLos[i]);
            eh_C.push_back(PeakErrHighs[i]);
        }
    }

    // =======================================================
    // GLOBAL PLOT DEFINITIONS
    // =======================================================
    TGraphAsymmErrors grpeaks(crystalNos.size(), &crystalNos[0], &Peaks[0], nullptr, nullptr, &PeakErrLos[0], &PeakErrHighs[0]);
    TGraph grchi2(crystalNos.size(), &crystalNos[0], &ChiSqs[0]);
    TGraphErrors grN(crystalNos.size(), &crystalNos[0], &Events[0], nullptr, &EventsErr[0]);
    TGraph grconvstatus(convgstatus_vec.size(), &convgstatus_vec[0], &PeakErrHighs[0]);
    TGraph grconvvschi2(convgstatus_vec.size(), &convgstatus_vec[0], &ChiSqs[0]);

    // =======================================================
    // CANVAS LAYOUT: TOP (All), BOTTOM-LEFT (CsI), BOTTOM-RIGHT (LYSO)
    // =======================================================
    TCanvas canvas("canvas", "Scatter Plot", 1000, 900); 
    
    TPad *padTop = new TPad("padTop", "All", 0.0, 0.5, 1.0, 1.0);
    TPad *padBL  = new TPad("padBL",  "CsI", 0.0, 0.0, 0.5, 0.5); // LEFT -> CsI
    TPad *padBR  = new TPad("padBR",  "LYSO",  0.5, 0.0, 1.0, 0.5); // RIGHT -> LYSO
    padTop->Draw();
    padBL->Draw();
    padBR->Draw();

    // --- 1. TOP PLOT (All SiPMs) ---
    padTop->cd();
    grpeaks.SetTitle("All SiPMs;SiPM Number;MeV/ADC");
    grpeaks.SetMarkerStyle(20);
    grpeaks.SetMarkerSize(0.8);
    grpeaks.SetMarkerColor(kBlue);
    grpeaks.SetLineColor(kBlue);
    grpeaks.Draw("AP");
    padTop->Update();
    grpeaks.GetYaxis()->SetRangeUser(0.05, 0.08);

    // Lines for Top
    TLine *avgpeak = new TLine(cry_start, Peaks_avg, cry_end, Peaks_avg);
    TLine *nominalpeak = new TLine(cry_start, 0.0625, cry_end, 0.0625); 
    TLine *stddevpeak_abv = new TLine(cry_start, Peaks_avg + 2*Peaks_stddev, cry_end, Peaks_avg + 2*Peaks_stddev);    
    TLine *stddevpeak_blw = new TLine(cry_start, Peaks_avg - 2*Peaks_stddev, cry_end, Peaks_avg - 2*Peaks_stddev);        
    
    avgpeak->SetLineColor(kCyan); avgpeak->SetLineWidth(2);
    nominalpeak->SetLineColor(kMagenta); nominalpeak->SetLineWidth(2);
    stddevpeak_abv->SetLineColor(kOrange-3); stddevpeak_abv->SetLineWidth(2);    
    stddevpeak_blw->SetLineColor(kOrange); stddevpeak_blw->SetLineWidth(2);        
    
    avgpeak->Draw("same");
    nominalpeak->Draw("same");
    stddevpeak_abv->Draw("same");
    stddevpeak_blw->Draw("same");

    // Text Box for Top
    TPaveText *avg = new TPaveText(0.15, 0.75, 0.35, 0.65, "brNDC");
    avg->SetFillStyle(0); avg->SetBorderSize(0); avg->SetTextSize(0.03);
    avg->SetTextColor(kBlack); avg->SetTextFont(72); avg->SetFillColor(kWhite);
    avg->AddText(Form("Avg = %4.8f#pm%.8f", Peaks_avg, Peaks_stddev));
    avg->AddText(Form("Std Dev = %.2f%%", stddev_prcnt));   
    avg->AddText(Form("Err Hi Avg = %4.8f", avgpeakerr_hi));
    avg->AddText(Form("Err Lo Avg = %4.8f", avgpeakerr_lo));
    avg->Draw();
    
    TLegend *legend = new TLegend();
    legend->AddEntry(avgpeak, "Average Value", "L");
    legend->AddEntry(nominalpeak, "Nominal Value", "L");
    legend->AddEntry(stddevpeak_abv, "Std Dev", "L");
    legend->Draw();

    // --- 2. BOTTOM LEFT PLOT (CsI Only) ---
    padBL->cd();
    // Stats for CsI
    double sum_C = std::accumulate(p_C.begin(), p_C.end(), 0.0);
    double avg_val_C = p_C.empty() ? 0 : sum_C / p_C.size();
    double std_val_C = TMath::StdDev(p_C.begin(), p_C.end());

    TGraphAsymmErrors *grCsI = new TGraphAsymmErrors(crys_C.size(), &crys_C[0], &p_C[0], nullptr, nullptr, &el_C[0], &eh_C[0]);
    grCsI->SetTitle("CsI Crystals;SiPM Number;MeV/ADC");
    grCsI->SetMarkerStyle(20); 
    grCsI->SetMarkerSize(0.8); // Set to 0.8
    grCsI->SetMarkerColor(kViolet+2); grCsI->SetLineColor(kViolet+2);
    grCsI->Draw("AP");
    padBL->Update();
    grCsI->GetYaxis()->SetRangeUser(0.05, 0.08);

    TLine *l_avg_C = new TLine(padBL->GetUxmin(), avg_val_C, padBL->GetUxmax(), avg_val_C);
    TLine *l_hi_C  = new TLine(padBL->GetUxmin(), avg_val_C + 2*std_val_C, padBL->GetUxmax(), avg_val_C + 2*std_val_C);
    TLine *l_lo_C  = new TLine(padBL->GetUxmin(), avg_val_C - 2*std_val_C, padBL->GetUxmax(), avg_val_C - 2*std_val_C);

    l_avg_C->SetLineColor(kCyan); l_avg_C->SetLineWidth(2);
    l_hi_C->SetLineColor(kOrange-3); l_hi_C->SetLineWidth(2);
    l_lo_C->SetLineColor(kOrange); l_lo_C->SetLineWidth(2);
    l_avg_C->Draw("same"); l_hi_C->Draw("same"); l_lo_C->Draw("same");

    TPaveText *ptC = new TPaveText(0.15, 0.75, 0.45, 0.60, "brNDC");
    ptC->SetFillStyle(0); ptC->SetBorderSize(0); ptC->SetTextSize(0.035);
    ptC->AddText(Form("Avg = %4.8f #pm %.8f", avg_val_C, std_val_C));
    ptC->AddText(Form("Std Dev = %.2f%%", (std_val_C/avg_val_C)*100));
    ptC->Draw();

    // --- 3. BOTTOM RIGHT PLOT (LYSO Only) ---
    padBR->cd();
    // Stats for LYSO
    double sum_L = std::accumulate(p_L.begin(), p_L.end(), 0.0);
    double avg_val_L = p_L.empty() ? 0 : sum_L / p_L.size();
    double std_val_L = TMath::StdDev(p_L.begin(), p_L.end());
    
    TGraphAsymmErrors *grLYSO = new TGraphAsymmErrors(crys_L.size(), &crys_L[0], &p_L[0], nullptr, nullptr, &el_L[0], &eh_L[0]);
    grLYSO->SetTitle("LYSO Crystals;SiPM Number;MeV/ADC");
    grLYSO->SetMarkerStyle(20); 
    grLYSO->SetMarkerSize(0.8); // Set to 0.8
    grLYSO->SetMarkerColor(kAzure-3); grLYSO->SetLineColor(kAzure-3);
    grLYSO->Draw("AP");
    padBR->Update();
    grLYSO->GetYaxis()->SetRangeUser(0.05, 0.08);

    TLine *l_avg_L = new TLine(padBR->GetUxmin(), avg_val_L, padBR->GetUxmax(), avg_val_L);
    TLine *l_hi_L  = new TLine(padBR->GetUxmin(), avg_val_L + 2*std_val_L, padBR->GetUxmax(), avg_val_L + 2*std_val_L);
    TLine *l_lo_L  = new TLine(padBR->GetUxmin(), avg_val_L - 2*std_val_L, padBR->GetUxmax(), avg_val_L - 2*std_val_L);
    
    l_avg_L->SetLineColor(kCyan); l_avg_L->SetLineWidth(2);
    l_hi_L->SetLineColor(kOrange-3); l_hi_L->SetLineWidth(2);
    l_lo_L->SetLineColor(kOrange); l_lo_L->SetLineWidth(2);
    l_avg_L->Draw("same"); l_hi_L->Draw("same"); l_lo_L->Draw("same");

    TPaveText *ptL = new TPaveText(0.15, 0.75, 0.45, 0.60, "brNDC");
    ptL->SetFillStyle(0); ptL->SetBorderSize(0); ptL->SetTextSize(0.035);
    ptL->AddText(Form("Avg = %4.8f #pm %.8f", avg_val_L, std_val_L));
    ptL->AddText(Form("Std Dev = %.2f%%", (std_val_L/avg_val_L)*100));
    ptL->Draw();

    outputFile->cd();
    grpeaks.Write("Peaks");
    canvas.SaveAs("Peaks.root");
// =======================================================
// SUBSET COMPARISON PLOTS (ODD / EVEN / RANGE1 / RANGE2)
// =======================================================

// -------- USER-EDITABLE RANGES --------
int cry_start1 = 0;   // Range 1 start
int cry_end1   = 400;     // Range 1 end

int cry_start2 = 900;         // <<< CHANGE RANGE 2 HERE
int cry_end2   = 1348;         // <<< CHANGE RANGE 2 HERE
// ------------------------------------

// --- Graphs ---
TGraphAsymmErrors *grOdd    = new TGraphAsymmErrors();
TGraphAsymmErrors *grEven   = new TGraphAsymmErrors();
TGraphAsymmErrors *grRange1 = new TGraphAsymmErrors();
TGraphAsymmErrors *grRange2 = new TGraphAsymmErrors();

// --- Storage for stats ---
std::vector<double> PeaksOdd, PeaksEven, PeaksRange1, PeaksRange2;

int iOdd = 0, iEven = 0, iR1 = 0, iR2 = 0;

// --- Fill graphs and stat vectors ---
for (size_t i = 0; i < crystalNos.size(); ++i) {

    double x   = crystalNos[i];
    double y   = Peaks[i];
    double eyL = PeakErrLos[i];
    double eyH = PeakErrHighs[i];

    if (((int)x) % 2 == 1) {
        grOdd->SetPoint(iOdd, x, y);
        grOdd->SetPointError(iOdd, 0, 0, eyL, eyH);
        PeaksOdd.push_back(y);
        ++iOdd;
    }

    if (((int)x) % 2 == 0) {
        grEven->SetPoint(iEven, x, y);
        grEven->SetPointError(iEven, 0, 0, eyL, eyH);
        PeaksEven.push_back(y);
        ++iEven;
    }

    if (x >= cry_start1 && x <= cry_end1) {
        grRange1->SetPoint(iR1, x, y);
        grRange1->SetPointError(iR1, 0, 0, eyL, eyH);
        PeaksRange1.push_back(y);
        ++iR1;
    }

    if (x >= cry_start2 && x <= cry_end2) {
        grRange2->SetPoint(iR2, x, y);
        grRange2->SetPointError(iR2, 0, 0, eyL, eyH);
        PeaksRange2.push_back(y);
        ++iR2;
    }
}

// --- Compute averages and std devs ---
auto Avg = [](const std::vector<double>& v){
    return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
};
auto Std = [](const std::vector<double>& v){
    return TMath::StdDev(v.begin(), v.end());
};

double avgOdd = Avg(PeaksOdd),     stdOdd = Std(PeaksOdd);
double avgEven = Avg(PeaksEven),   stdEven = Std(PeaksEven);
double avgR1 = Avg(PeaksRange1),   stdR1 = Std(PeaksRange1);
double avgR2 = Avg(PeaksRange2),   stdR2 = Std(PeaksRange2);

// --- Make lines ---
auto MakeLines = [](double xmin, double xmax, double avg, double std){
    TLine *avgL = new TLine(xmin, avg, xmax, avg);
    TLine *hiL  = new TLine(xmin, avg+2*std, xmax, avg+2*std);
    TLine *loL  = new TLine(xmin, avg-2*std, xmax, avg-2*std);

    avgL->SetLineColor(kCyan); avgL->SetLineWidth(2);
    hiL->SetLineColor(kOrange+1); hiL->SetLineWidth(2);
    loL->SetLineColor(kOrange+1); loL->SetLineWidth(2);

    return std::make_tuple(avgL, hiL, loL);
};

auto [avgOddL,  hiOddL,  loOddL ] = MakeLines(cry_start1, cry_end1, avgOdd,  stdOdd);
auto [avgEvenL, hiEvenL, loEvenL] = MakeLines(cry_start1, cry_end1, avgEven, stdEven);
auto [avgR1L,   hiR1L,   loR1L  ] = MakeLines(cry_start1, cry_end1, avgR1,   stdR1);
auto [avgR2L,   hiR2L,   loR2L  ] = MakeLines(cry_start2, cry_end2, avgR2,   stdR2);

// --- Text boxes ---
auto MakeText = [](double avg, double std){
    TPaveText *t = new TPaveText(0.15,0.75,0.35,0.65,"brNDC");
    t->SetFillStyle(0);
    t->SetBorderSize(0);
    t->SetTextSize(0.03);
    t->SetTextFont(72);
    t->AddText(Form("Avg = %4.8f", avg));
    t->AddText(Form("Std Dev = %.2f%%", (std/avg)*100));
    return t;
};

TPaveText *txtOdd  = MakeText(avgOdd, stdOdd);
TPaveText *txtEven = MakeText(avgEven, stdEven);
TPaveText *txtR1   = MakeText(avgR1, stdR1);
TPaveText *txtR2   = MakeText(avgR2, stdR2);

// --- Style graphs ---
grOdd->SetTitle("Odd SiPMs;SiPM Number;MeV/ADC");
grOdd->SetMarkerStyle(20); grOdd->SetMarkerColor(kRed+1);

grEven->SetTitle("Even SiPMs;SiPM Number;MeV/ADC");
grEven->SetMarkerStyle(20); grEven->SetMarkerColor(kBlue+1);

grRange1->SetTitle(Form("Range %d-%d;SiPM Number;MeV/ADC",cry_start1,cry_end1));
grRange1->SetMarkerStyle(20); grRange1->SetMarkerColor(kGreen+2);

grRange2->SetTitle(Form("Range %d-%d;SiPM Number;MeV/ADC",cry_start2,cry_end2));
grRange2->SetMarkerStyle(20); grRange2->SetMarkerColor(kMagenta+1);

// --- Canvas ---
TCanvas *cSplit = new TCanvas("cSplit","Subset Comparison",1600,400);
cSplit->Divide(4,1);

// --- Draw pads ---
auto DrawPad = [&](int pad, TGraphAsymmErrors* g,
                   TLine* a, TLine* h, TLine* l, TPaveText* t){
    cSplit->cd(pad);
    g->Draw("AP");
    g->GetYaxis()->SetRangeUser(0.05,0.08);
    a->Draw("same"); h->Draw("same"); l->Draw("same");
    t->Draw();
};

DrawPad(1, grOdd,    avgOddL,  hiOddL,  loOddL,  txtOdd);
DrawPad(2, grEven,   avgEvenL, hiEvenL, loEvenL, txtEven);
DrawPad(3, grRange1, avgR1L,   hiR1L,   loR1L,   txtR1);
DrawPad(4, grRange2, avgR2L,   hiR2L,   loR2L,   txtR2);

// --- Write output ---
outputFile->cd();
grOdd->Write("Peaks_Odd");
grEven->Write("Peaks_Even");
grRange1->Write("Peaks_Range1");
grRange2->Write("Peaks_Range2");
cSplit->SaveAs("Peaks_SubsetComparison.root");



//----------------------------------------//

    //TCanvas canvas("canvas", "Scatter Plot", 800, 600);
    grpeaks.Draw("APsame"); // A for axis, P for points
    outputFile->cd();
    grpeaks.Write("Peaks");
    // ---- save peaks and errors into a txt file ----
std::ofstream outFile("SiPM_FitParameters.txt");
if (outFile.is_open()) {
    outFile << "#SiPMNo\tPeak(MeV/ADC)\tPeakErrHigh\tPeakErrLo\tunredChi2\tndof\tnEvents\tWidth\tWidthErrHigh\tWidthErrLo\tChi2\n";
    for (size_t i = 0; i < crystalNos.size(); ++i) {
        outFile << std::fixed << std::setprecision(6)
                << crystalNos[i] << "\t"
                << Peaks[i]      << "\t"
                << PeakErrHighs[i]   << "\t"
                << PeakErrLos[i]   << "\t"
                << unreducedchi2s[i]<< "\t"
                << ndofs[i]<< "\t"
                << Events[i] << "\t"
                << Widths[i] << "\t"    
                << WidthErrHighs[i] << "\t" 
                << WidthErrLos[i] << "\t"    
                << ChiSqs[i] << "\t"
                << Events[i] << "\n";
    }
    outFile.close();
    std::cout << "All fit parameters written to SiPM_FitParameters.txt" << std::endl;
} else {
    std::cerr << "Error: could not open SiPM_FitParameters.txt for writing." << std::endl;
}

//create print out statment of convergance status percentages:
int total = convgstatus_vec.size();
int count0 = 0,count1=0,count2=0,count3=0,count4=0,count5=0,countx=0;
for (auto s:convgstatus_vec){
	if (s == 0) count0++;
	else if (s == 1) count1++;
	else if (s == 2) count2++;
	else if (s == 3) count3++;	
	else if (s == 4) count4++;			
	else if (s == 5) count5++;
	else if (s>5) countx++;
	}
if (total == 0) {
    std::cout << "No SiPMs in convgstatus_vec!" << std::endl;
} else {
double perc0 =100*(static_cast<double>(count0) /total);
double perc1 =100*(static_cast<double>(count1) /total);
double perc2 =100*(static_cast<double>(count2) /total);
double perc3 =100*(static_cast<double>(count3) /total);
double perc4 =100*(static_cast<double>(count4) /total);
double perc5 =100*(static_cast<double>(count5) /total);
double percx =100*(static_cast<double>(countx) /total);
std::cout << std::fixed << std::setprecision(2);
std::cout<< "% of sipms with convg status 0 ="<<perc0 <<std::endl;
std::cout<< "% of sipms with convg status 1 ="<< perc1<<std::endl;
std::cout<< "% of sipms with convg status 2 ="<< perc2<<std::endl;
std::cout<< "% of sipms with convg status 3 ="<< perc3<<std::endl;
std::cout<< "% of sipms with convg status 4 ="<< perc4<<std::endl;
std::cout<< "% of sipms with convg status 5 ="<< perc5<<std::endl;
std::cout<< "% of sipms with convg status unknown ="<< percx<<std::endl;}
  
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

//-----convergence status plots that come from the nll fit function vs peak error-------- 
    grconvstatus.SetTitle("Peak Error vs Convg Status ;Convergence Status; PeakErrorHigh");
    grconvstatus.SetMarkerStyle(20);
    grconvstatus.SetMarkerSize(0.8);
    grconvstatus.SetMarkerColor(kRed);
    grconvstatus.SetLineColor(kRed);
    TCanvas canvas_conv("canvas", "Scatter Plot", 800, 600);
    grconvstatus.Draw("AP"); // A for axis, P for points
    outputFile->cd();
    grconvstatus.Write("Convstatus");
    canvas_conv.SaveAs("convstatusvspeakerrs.root");  
    
//-----convergence status vs chi2 plots-------- 
    grconvvschi2.SetTitle("Peak Error vs Convg Status ;Convergence Status; Chi2");
    grconvvschi2.SetMarkerStyle(20);
    grconvvschi2.SetMarkerSize(0.8);
    grconvvschi2.SetMarkerColor(kRed);
    grconvvschi2.SetLineColor(kRed);
    TCanvas canvas_convvschi2("canvas", "Scatter Plot", 800, 600);
    grconvvschi2.Draw("AP"); // A for axis, P for points
    outputFile->cd();
    grconvvschi2.Write("Convstatusvschi2");
    canvas_convvschi2.SaveAs("convstatusvschi2.root");  
   
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
// =======================================================
// ODD vs EVEN + 45 DEGREE PROJECTION
// =======================================================

// --- Graphs ---
TGraph *grOddEven = new TGraph();
TH1D *hProj45 = new TH1D("hProj45",
    "Projection Along 45^{#circ};(Odd+Even)/#sqrt{2};Counts",
    60, 0.05*sqrt(2), 0.08*sqrt(2));

TH1D *hProjPerp = new TH1D("hProjPerp",
    "Projection Perpendicular to 45^{#circ};(Odd-Even)/#sqrt{2};Counts",
    60, -0.01, 0.01);

int p = 0;

// --- Build odd-even pairs ---
for (size_t i = 0; i + 1 < crystalNos.size(); ++i) {

    int sipm = (int)crystalNos[i];

    // enforce even-odd pairing
    if (sipm % 2 != 0) continue;

    double evenVal = Peaks[i];
    double oddVal  = Peaks[i+1];

    // Scatter plot
    grOddEven->SetPoint(p++, evenVal, oddVal);

    // Rotated coordinates
    double u = (evenVal + oddVal) / sqrt(2.0);   // along 45°
    double v = (oddVal  - evenVal) / sqrt(2.0);  // perpendicular

    hProj45->Fill(u);
    hProjPerp->Fill(v);
}
TCanvas *c45 = new TCanvas("c45","Odd-Even 45 Degree Analysis",1200,400);
c45->Divide(3,1);

// --- Odd vs Even scatter ---
c45->cd(1);
grOddEven->SetTitle("Odd vs Even SiPM Calibration;Even MeV/ADC;Odd MeV/ADC");
grOddEven->SetMarkerStyle(20);
grOddEven->SetMarkerSize(0.7);
grOddEven->Draw("AP");

// 45 degree reference line
TLine *diag = new TLine(0.05,0.05,0.08,0.08);
diag->SetLineColor(kRed);
diag->SetLineWidth(2);
diag->Draw("same");

// --- Projection along 45° ---
c45->cd(2);
hProj45->SetLineColor(kBlue+1);
hProj45->SetLineWidth(2);
hProj45->Draw();

// --- Projection perpendicular to 45° ---
c45->cd(3);
hProjPerp->SetLineColor(kMagenta+1);
hProjPerp->SetLineWidth(2);
hProjPerp->Draw();

// --- Save ---
outputFile->cd();
grOddEven->Write("OddVsEven");
hProj45->Write("Proj_45deg");
hProjPerp->Write("Proj_Perp45");
c45->SaveAs("OddEven_45deg.root");

//----------------------------------------------//

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
double ymin_adc = std::numeric_limits<double>::infinity();
double ymax_adc = -std::numeric_limits<double>::infinity();
auto updateMinMax = [&](const TGraph& g) {
    for (int i = 0; i < g.GetN(); ++i) {
        double x, y;
        g.GetPoint(i, x, y);
        if (std::isfinite(y)) {
            ymin_adc = std::min(ymin_adc, y);
            ymax_adc = std::max(ymax_adc, y);
        }
    }
};
updateMinMax(grADC);
updateMinMax(gr1st);
updateMinMax(gr2nd);

// Add margin
double dy_adc = (ymax_adc - ymin_adc) * 0.05;
if (dy_adc <= 0) dy_adc = 1.0; // fallback
ymin_adc -= dy_adc;
ymax_adc += dy_adc;

// Draw empty frame first with full range
TH1F* frame_adc = cADC1st2nd.DrawFrame(
    crystalNos.front(), ymin_adc,
    crystalNos.back(), ymax_adc
);
frame_adc->SetTitle("Peak Locations;Crystal Number;ADC / Peak Value");

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

/*// ---- Combined Fraction Plot: frFull, frFrst, frScnd ----

// Histogram settings
int nbins_frac = 50;
double frac_min = 0.0;
double frac_max = 1.0;

// Create histograms
TH1D* hFrFull = new TH1D("hFrFull", "SiPM Fraction Distributions;Fraction;Number of SiPMs", nbins_frac, frac_min, frac_max);
TH1D* hFrFrst = new TH1D("hFrFrst", "SiPM Fraction Distributions;Fraction;Number of SiPMs", nbins_frac, frac_min, frac_max);
TH1D* hFrScnd = new TH1D("hFrScnd", "SiPM Fraction Distributions;Fraction;Number of SiPMs", nbins_frac, frac_min, frac_max);
//;Normalized Density

// Fill histograms
for (double v : frFull_vec) hFrFull->Fill(v);
for (double v : frFrst_vec) hFrFrst->Fill(v);
for (double v : frScnd_vec) hFrScnd->Fill(v);


// ---- Overlay Plot ----
TCanvas* cOverlay = new TCanvas("cOverlay", "Overlay of SiPM Fraction Distributions", 800, 600);

hFrFull->SetLineColor(kBlue);
hFrFull->SetLineWidth(2);

hFrFrst->SetLineColor(kGreen + 2);
hFrFrst->SetLineWidth(2);

hFrScnd->SetLineColor(kMagenta + 1);
hFrScnd->SetLineWidth(2);

// Compute global max among all three histograms
double maxY = std::max({hFrFull->GetMaximum(), hFrFrst->GetMaximum(), hFrScnd->GetMaximum()});

// Apply it to the one drawn first (defines the axes)
hFrFull->SetMaximum(1.1 * maxY);  // add ~10% headroom
hFrFull->SetMinimum(0);

// Draw first and overlay others
hFrFull->Draw("HIST");
hFrFrst->Draw("HIST SAME");
hFrScnd->Draw("HIST SAME");

// Legend
TLegend* legOverlay = new TLegend(0.65, 0.7, 0.9, 0.9);
legOverlay->AddEntry(hFrFull, "frFull (Peak 1)", "l");
legOverlay->AddEntry(hFrFrst, "frFrst (Peak 2)", "l");
legOverlay->AddEntry(hFrScnd, "frScnd (Peak 3)", "l");
legOverlay->Draw();

// Save overlay plot
cOverlay->SaveAs("SiPM_Fractions_Overlay.root");

// ---- Stacked Plot ----
THStack* hStack = new THStack("hStack", "Stacked SiPM Fractions;Fraction; Number of SiPMs");

hFrFull->SetFillColor(kBlue - 9);
hFrFrst->SetFillColor(kGreen - 7);
hFrScnd->SetFillColor(kMagenta - 9);

hStack->Add(hFrFull);
hStack->Add(hFrFrst);
hStack->Add(hFrScnd);

TCanvas* cStack = new TCanvas("cStack", "Stacked SiPM Fractions", 800, 600);
hStack->Draw("HIST");

TLegend* legStack = new TLegend(0.65, 0.7, 0.9, 0.9);
legStack->AddEntry(hFrFull, "frFull (Peak 1)", "f");
legStack->AddEntry(hFrFrst, "frFrst (Peak 2)", "f");
legStack->AddEntry(hFrScnd, "frScnd (Peak 3)", "f");
legStack->Draw();

// ---- Save Histograms to Output File ----
outputFile->cd();
hFrFull->Write("SiPM_FrFull_Hist");
hFrFrst->Write("SiPM_FrFrst_Hist");
hFrScnd->Write("SiPM_FrScnd_Hist");
cStack->SaveAs("SiPM_Fractions_Hist.root");*/
// =======================================================
// Peak Location vs Fraction ? SCATTER (MeV space)
// =======================================================

const double E_min = 40.0;
const double E_max = 120.0;

// --- Graphs: one point per fit ---
TGraph* gMain   = new TGraph();
TGraph* gFirst  = new TGraph();
TGraph* gSecond = new TGraph();

gMain  ->SetName("gMainPeak_MeV_vs_Fraction");
gFirst ->SetName("gFirstPeak_MeV_vs_Fraction");
gSecond->SetName("gSecondPeak_MeV_vs_Fraction");

// ---- Fill graphs ----
int iMain = 0, iFirst = 0, iSecond = 0;

for (size_t i = 0; i < ADCPeaks.size(); ++i) {
    gMain  ->SetPoint(iMain++,   ADCPeaks[i],        frFull_vec[i]);
    gFirst ->SetPoint(iFirst++,  firstPeak_vec[i],    frFrst_vec[i]);
    gSecond->SetPoint(iSecond++, secondPeak_vec[i],   frScnd_vec[i]);
}

// ---- Styling ----
gMain->SetMarkerStyle(20);
gMain->SetMarkerSize(1.2);
gMain->SetMarkerColor(kOrange + 7);
gMain->SetLineColor(kOrange + 7);

gFirst->SetMarkerStyle(21);
gFirst->SetMarkerSize(1.2);
gFirst->SetMarkerColor(kMagenta + 1);
gFirst->SetLineColor(kMagenta + 1);

gSecond->SetMarkerStyle(22);
gSecond->SetMarkerSize(1.2);
gSecond->SetMarkerColor(kCyan + 1);
gSecond->SetLineColor(kCyan + 1);

// ---- Canvas ----
TCanvas* cPeakFrac = new TCanvas(
    "cPeakFrac",
    "Peak Location vs Fraction (per fit)",
    900, 650
);

// ---- Frame (DO NOT STORE POINTER) ----
cPeakFrac->DrawFrame(
    E_min, 0.0,
    E_max, 1.05,
    "Peak Location vs Fraction;Peak Energy [MeV];Fraction"
);

// ---- Draw graphs ----
gMain  ->Draw("P SAME");
gFirst ->Draw("P SAME");
gSecond->Draw("P SAME");

// ---- Legend ----
TLegend* leg = new TLegend(0.62, 0.72, 0.88, 0.88);
leg->AddEntry(gMain,   "Main Peak",   "p");
leg->AddEntry(gFirst,  "1st Escape",  "p");
leg->AddEntry(gSecond, "2nd Escape",  "p");
leg->Draw();

// ---- Save (ROOT-style, no extra files) ----
outputFile->cd();
gMain  ->Write();
gFirst ->Write();
gSecond->Write();
cPeakFrac->Write("PeakLocation_vs_Fraction");

cPeakFrac->SaveAs("PeakLocation_vs_Fraction.root");




}
