#include "CaloCalibration/SourceCalib/inc/mcinfo.hh"
#include "RooRealVar.h"
//#include "RooCBShape.h"
//#include "RooCrystalBall.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooFormulaVar.h"
#include "RooGenericPdf.h"
#include "RooMinimizer.h"
#include "TF1.h"

using namespace std::chrono;
using namespace CaloSourceCalib;
using namespace RooFit;

void mcinfo::RunMCTruth(TH1F* hist, int cryNum, int disk, TTree *trueinfo,
                        Int_t &cryNumparam, Int_t &tot_evts,
                        Int_t &mainpeak, Int_t &first_espeak, Int_t &second_espeak,
                        Int_t &background, Float_t &frmainpeak,
                        Float_t &frfirst_espeak, Float_t &frsecond_espeak, Float_t &frbackground,
                        Int_t &compton1, Int_t &compton2, Int_t &compton3,
                        Float_t &frcompton1, Float_t &frcompton2, Float_t &frcompton3)
{
    // 1. Define Energy Constants (MeV)
    // Using constexpr ensures these are known at compile time and efficient.
    static constexpr double E_main_MeV    = 6.13;
    static constexpr double E_1st_MeV     = 5.619;
    static constexpr double E_2nd_MeV     = 5.108;
    static constexpr double E_compton_MeV = 0.511;
    // Compton edge energies matching SourceFitter CE1/CE2/CE3
    static constexpr double CE1_MeV = 5.8842;
    static constexpr double CE2_MeV = 5.3741;
    static constexpr double CE3_MeV = 4.8640;

    // 2. Set Analysis Variables
    // Default to MeV scale (standard for MC Truth)
    int bin_lo = hist->FindBin(2.5);
    int bin_hi = hist->FindBin(6.2);

    double E_main    = E_main_MeV;
    double E_1st     = E_1st_MeV;
    double E_2nd     = E_2nd_MeV;
    double E_compton = E_compton_MeV;

    //uncomment if this file is being applied to the reconstructed hists
    /* int bin_lo = hist->FindBin(2.5)/0.0625;
    int bin_hi = hist->FindBin(6.2)/0.0625;
    double E_main  = 6.13/0.0625;
    double E_1st   = 5.619/0.0625;
    double E_2nd   = 5.108/0.0625;
    double E_compton = 0.511/0.0625;
    */

    // 3. Define Peak Windows (Bin Indices)
    // We calculate the bin index for the peak center and the compton edge.
    auto ax = hist->GetXaxis();
    int bin_main_center = ax->FindBin(E_main);
    int bin_1st_center  = ax->FindBin(E_1st);
    int bin_2nd_center  = ax->FindBin(E_2nd);
    // Logic: The window for a peak ends at its center, and starts right after the previous peak ends.
    
    // 2nd Escape Peak Window: [ (E_2nd - 0.511), E_2nd ]
    int bin_2nd_low  = ax->FindBin(E_2nd - E_compton);
    int bin_2nd_high = bin_2nd_center;

    // 1st Escape Peak Window: [ E_2nd + 1 bin, E_1st ]
    int bin_1st_low  = bin_2nd_high + 1; 
    int bin_1st_high = bin_1st_center;

    // Main Peak Window: [ E_1st + 1 bin, E_main ]
    int bin_main_low  = bin_1st_high + 1;
    int bin_main_high = bin_main_center;

    // 4. Compton edge bin positions (matching SourceFitter CE1/CE2/CE3)
    int bin_ce1 = ax->FindBin(CE1_MeV);
    int bin_ce2 = ax->FindBin(CE2_MeV);
    int bin_ce3 = ax->FindBin(CE3_MeV);

    // 5. Non-overlapping integrals — each peak window is split at its Compton edge:
    //   Compton window : [peak_window_low, CE]      (erfc plateau below edge)
    //   Peak window    : [CE+1bin, peak_window_high] (CB peak above edge)
    //
    //   Full partition of [4.597, 6.13]:
    //     compton3  [4.597,       CE3=4.864]
    //     2nd esc   [CE3+1,       5.108    ]
    //     compton2  [5.108+1,     CE2=5.374]
    //     1st esc   [CE2+1,       5.619    ]
    //     compton1  [5.619+1,     CE1=5.884]
    //     mainpeak  [CE1+1,       6.13     ]
    compton1     = hist->Integral(bin_main_low,  bin_ce1);          // [~5.619, 5.884]
    mainpeak     = hist->Integral(bin_ce1 + 1,   bin_main_high);    // [5.884,  6.13 ]
    compton2     = hist->Integral(bin_1st_low,   bin_ce2);          // [~5.108, 5.374]
    first_espeak = hist->Integral(bin_ce2 + 1,   bin_1st_high);     // [5.374,  5.619]
    compton3     = hist->Integral(bin_2nd_low,   bin_ce3);          // [4.597,  4.864]
    second_espeak= hist->Integral(bin_ce3 + 1,   bin_2nd_high);     // [4.864,  5.108]

    tot_evts = hist->Integral(bin_lo, bin_hi);
    double tot = (tot_evts > 0) ? static_cast<double>(tot_evts) : 1.0;

    frmainpeak      = mainpeak      / tot;
    frfirst_espeak  = first_espeak  / tot;
    frsecond_espeak = second_espeak / tot;
    frcompton1      = compton1      / tot;
    frcompton2      = compton2      / tot;
    frcompton3      = compton3      / tot;

    background   = tot_evts - (mainpeak + first_espeak + second_espeak + compton1 + compton2 + compton3);
    frbackground = background / tot;

    trueinfo->Fill();
}

void mcinfo::FinalizeMCSummary(TTree* trueinfo, TTree* truthfit,bool useCompton)
{
    std::cout << "[INFO] Building global MC truth histograms...\n";

    // --- 1. Initialize Histograms ---
    TH1F* h_tot      = new TH1F("h_tot_events",  "Total MC events per crystal", 200, 0, 40000);
    
    // Unweighted histograms
    TH1F* h_frMain   = new TH1F("h_frMain",   "Fraction Main Peak",  100, 0, 1);
    TH1F* h_fr1st    = new TH1F("h_fr1st",    "Fraction 1st Escape", 100, 0, 1);
    TH1F* h_fr2nd    = new TH1F("h_fr2nd",    "Fraction 2nd Escape", 100, 0, 1);
    TH1F* h_frBg     = new TH1F("h_frBg",     "Fraction Background", 100, 0, 1);
    TH1F* h_frC1     = new TH1F("h_frC1",     "Fraction Compton 1 (below CE1=5.884)", 100, 0, 1);
    TH1F* h_frC2     = new TH1F("h_frC2",     "Fraction Compton 2 (below CE2=5.374)", 100, 0, 1);
    TH1F* h_frC3     = new TH1F("h_frC3",     "Fraction Compton 3 (below CE3=4.864)", 100, 0, 1);

    // Weighted histograms (weighted by total events)
    TH1F* h_frMain_w = new TH1F("h_frMain_w", "Weighted Main Peak Fraction",    100, 0, 1);
    TH1F* h_fr1st_w  = new TH1F("h_fr1st_w",  "Weighted 1st Escape Fraction",   100, 0, 1);
    TH1F* h_fr2nd_w  = new TH1F("h_fr2nd_w",  "Weighted 2nd Escape Fraction",   100, 0, 1);
    TH1F* h_frC1_w   = new TH1F("h_frC1_w",   "Weighted Compton 1 Fraction",    100, 0, 1);
    TH1F* h_frC2_w   = new TH1F("h_frC2_w",   "Weighted Compton 2 Fraction",    100, 0, 1);
    TH1F* h_frC3_w   = new TH1F("h_frC3_w",   "Weighted Compton 3 Fraction",    100, 0, 1);

    // --- 2. Process TTree Data ---
    Int_t   tot_evts;
    Float_t frmainpeak, frfirst_espeak, frsecond_espeak, frbackground;
    Float_t frcompton1, frcompton2, frcompton3;

    trueinfo->SetBranchAddress("tot_evts",         &tot_evts);
    trueinfo->SetBranchAddress("frmainpeak",       &frmainpeak);
    trueinfo->SetBranchAddress("frfirst_espeak",   &frfirst_espeak);
    trueinfo->SetBranchAddress("frsecond_espeak",  &frsecond_espeak);
    trueinfo->SetBranchAddress("frbackground",     &frbackground);
    trueinfo->SetBranchAddress("frcompton1",       &frcompton1);
    trueinfo->SetBranchAddress("frcompton2",       &frcompton2);
    trueinfo->SetBranchAddress("frcompton3",       &frcompton3);

    Long64_t sumTot = 0, sumMain = 0, sum1st = 0, sum2nd = 0, sumBg = 0;
    Long64_t sumC1 = 0, sumC2 = 0, sumC3 = 0;
    Long64_t N = trueinfo->GetEntries();

    for (Long64_t i = 0; i < N; i++) {
        trueinfo->GetEntry(i);

        sumTot  += tot_evts;
        sumMain += frmainpeak * tot_evts;
        sum1st  += frfirst_espeak * tot_evts;
        sum2nd  += frsecond_espeak * tot_evts;
        sumBg   += frbackground * tot_evts;
        sumC1   += frcompton1 * tot_evts;
        sumC2   += frcompton2 * tot_evts;
        sumC3   += frcompton3 * tot_evts;

        h_tot->Fill(tot_evts);

        // Fill Unweighted
        h_frMain->Fill(frmainpeak);
        h_fr1st->Fill(frfirst_espeak);
        h_fr2nd->Fill(frsecond_espeak);
        h_frBg->Fill(frbackground);
        h_frC1->Fill(frcompton1);
        h_frC2->Fill(frcompton2);
        h_frC3->Fill(frcompton3);

        // Fill Weighted
        h_frMain_w->Fill(frmainpeak,      tot_evts);
        h_fr1st_w ->Fill(frfirst_espeak,  tot_evts);
        h_fr2nd_w ->Fill(frsecond_espeak, tot_evts);
        h_frC1_w  ->Fill(frcompton1,      tot_evts);
        h_frC2_w  ->Fill(frcompton2,      tot_evts);
        h_frC3_w  ->Fill(frcompton3,      tot_evts);
    }

    // --- 3. Calculate Global Averages ---
    double fMain = (sumTot > 0) ? double(sumMain) / sumTot : 0;
    double f1st  = (sumTot > 0) ? double(sum1st)  / sumTot : 0;
    double f2nd  = (sumTot > 0) ? double(sum2nd)  / sumTot : 0;
    double fBg   = (sumTot > 0) ? double(sumBg)   / sumTot : 0;
    double fC1   = (sumTot > 0) ? double(sumC1)   / sumTot : 0;
    double fC2   = (sumTot > 0) ? double(sumC2)   / sumTot : 0;
    double fC3   = (sumTot > 0) ? double(sumC3)   / sumTot : 0;

    // Normalize weighted histograms to integral 1 for comparison
    if(h_frMain_w->Integral() > 0) h_frMain_w->Scale(1.0 / h_frMain_w->Integral());
    if(h_fr1st_w->Integral()  > 0) h_fr1st_w ->Scale(1.0 / h_fr1st_w->Integral());
    if(h_fr2nd_w->Integral()  > 0) h_fr2nd_w ->Scale(1.0 / h_fr2nd_w->Integral());
    if(h_frC1_w->Integral()   > 0) h_frC1_w  ->Scale(1.0 / h_frC1_w->Integral());
    if(h_frC2_w->Integral()   > 0) h_frC2_w  ->Scale(1.0 / h_frC2_w->Integral());
    if(h_frC3_w->Integral()   > 0) h_frC3_w  ->Scale(1.0 / h_frC3_w->Integral());

    std::cout << std::fixed << std::setprecision(7);
    std::cout << "\n===== CUT & COUNT GLOBAL MC TRUTH SUMMARY =====\n";
    std::cout << "Total events:        " << sumTot << "\n";
    std::cout << "Main peak fraction:  " << fMain << "\n";
    std::cout << "1st escape fraction: " << f1st << "\n";
    std::cout << "2nd escape fraction: " << f2nd << "\n";
    std::cout << "Background fraction: " << fBg  << "\n";
    std::cout << "Compton1 fraction:   " << fC1  << "\n";
    std::cout << "Compton2 fraction:   " << fC2  << "\n";
    std::cout << "Compton3 fraction:   " << fC3  << "\n";
    std::cout << "Sum of fractions:    " << fBg+fMain+f1st+f2nd+fC1+fC2+fC3 << std::endl;
    std::cout << "==================================\n\n";

    TFile* outFile = new TFile("MC_Truth_Summary.root", "RECREATE");

    // =======================================================
    // PLOT 1 & 2: Overlay Histograms (Unweighted & Weighted)
    // =======================================================

    // Helper Lambda to draw the overlay plots (peaks + Compton components)
    // Peak colors match SourceFitter plotOn conventions; Compton use kGreen shades (dashed).
    auto draw_overlay = [&](const char* name, const char* title,
                            TH1F* hM,  TH1F* h1,  TH1F* h2,
                            TH1F* hC1, TH1F* hC2, TH1F* hC3,
                            const char* leg_suffix)
    {
        TCanvas* c = new TCanvas(name, title, 900, 600);

        // Peak styles (solid)
        hM ->SetLineColor(kBlue);      hM ->SetLineWidth(2); hM ->SetLineStyle(kSolid);
        h1 ->SetLineColor(kOrange+1);  h1 ->SetLineWidth(2); h1 ->SetLineStyle(kSolid);
        h2 ->SetLineColor(kMagenta+1); h2 ->SetLineWidth(2); h2 ->SetLineStyle(kSolid);
        // Compton styles (dashed, green shades matching SourceFitter)
        hC1->SetLineColor(kGreen);     hC1->SetLineWidth(2); hC1->SetLineStyle(kDashed);
        hC2->SetLineColor(kGreen+2);   hC2->SetLineWidth(2); hC2->SetLineStyle(kDashed);
        hC3->SetLineColor(kGreen+4);   hC3->SetLineWidth(2); hC3->SetLineStyle(kDashed);

        double max_val = std::max({ hM->GetMaximum(),  h1->GetMaximum(),  h2->GetMaximum(),
                                    hC1->GetMaximum(), hC2->GetMaximum(), hC3->GetMaximum() });
        hM->SetMaximum(1.1 * max_val);
        hM->SetMinimum(0);

        hM ->Draw("HIST");
        h1 ->Draw("HIST SAME");
        h2 ->Draw("HIST SAME");
        hC1->Draw("HIST SAME");
        hC2->Draw("HIST SAME");
        hC3->Draw("HIST SAME");

        TLegend* lg = new TLegend(0.60, 0.55, 0.92, 0.92);
        lg->SetBorderSize(0);
        lg->SetFillStyle(0);
        lg->AddEntry(hM,  Form("Main Peak %s",          leg_suffix), "l");
        lg->AddEntry(h1,  Form("1st Escape %s",         leg_suffix), "l");
        lg->AddEntry(h2,  Form("2nd Escape %s",         leg_suffix), "l");
        lg->AddEntry(hC1, Form("Compton 1 (CE=5.884) %s", leg_suffix), "l");
        lg->AddEntry(hC2, Form("Compton 2 (CE=5.374) %s", leg_suffix), "l");
        lg->AddEntry(hC3, Form("Compton 3 (CE=4.864) %s", leg_suffix), "l");
        lg->Draw();

        c->Write(name);
        return c;
    };

    draw_overlay("MC_Fractions_Overlay", "MC Truth Fraction Overlay",
                 h_frMain, h_fr1st, h_fr2nd, h_frC1, h_frC2, h_frC3, "");

    draw_overlay("MC_Fractions_Overlay_weighted", "MC Truth Fraction Overlay (Weighted)",
                 h_frMain_w, h_fr1st_w, h_fr2nd_w, h_frC1_w, h_frC2_w, h_frC3_w, "(wt)");

    // =======================================================
    // PLOT 3: Peak Location vs Fraction (SCATTER)
    // =======================================================

    // Helper Lambda to create and style a single point graph
    auto create_truth_graph = [&](double x, double y, int marker, int color) {
        TGraph* g = new TGraph(1);
        g->SetPoint(0, x, y);
        g->SetMarkerStyle(marker);
        g->SetMarkerSize(1.6);
        g->SetMarkerColor(color);
        return g;
    };

    // Peak positions and Compton edge positions (matching SourceFitter CE1/CE2/CE3)
    const double E_main = 6.13,  CE1 = 5.8842;
    const double E_1st  = 5.619, CE2 = 5.3741;
    const double E_2nd  = 5.108, CE3 = 4.8640;

    // Peak graphs (filled markers)
    TGraph* gTruthMain   = create_truth_graph(E_main, fMain, 20, kBlue);        // Filled circle
    TGraph* gTruthFirst  = create_truth_graph(E_1st,  f1st,  21, kOrange+1);    // Filled square
    TGraph* gTruthSecond = create_truth_graph(E_2nd,  f2nd,  22, kMagenta+1);   // Filled triangle
    // Compton edge graphs (open markers, green shades matching SourceFitter)
    TGraph* gTruthC1     = create_truth_graph(CE1,    fC1,   24, kGreen);        // Open circle
    TGraph* gTruthC2     = create_truth_graph(CE2,    fC2,   25, kGreen+2);      // Open square
    TGraph* gTruthC3     = create_truth_graph(CE3,    fC3,   26, kGreen+4);      // Open triangle

    TCanvas* cTruth = new TCanvas("cTruth", "MC Truth Peak Location vs Fraction", 900, 650);
    cTruth->DrawFrame(2.5, 0.0, 8.0, 1.05, "MC Truth: Peak Location vs Fraction;Peak Energy [MeV];Fraction of Events");

    gTruthMain->Draw("P SAME");
    gTruthFirst->Draw("P SAME");
    gTruthSecond->Draw("P SAME");
    gTruthC1->Draw("P SAME");
    gTruthC2->Draw("P SAME");
    gTruthC3->Draw("P SAME");

    TLegend* legTruth = new TLegend(0.60, 0.55, 0.92, 0.92);
    legTruth->SetBorderSize(0);
    legTruth->SetFillStyle(0);
    legTruth->AddEntry(gTruthMain,   "MC Main Peak",           "p");
    legTruth->AddEntry(gTruthFirst,  "MC 1st Escape",          "p");
    legTruth->AddEntry(gTruthSecond, "MC 2nd Escape",          "p");
    legTruth->AddEntry(gTruthC1,     "MC Compton 1 (CE=5.884)", "p");
    legTruth->AddEntry(gTruthC2,     "MC Compton 2 (CE=5.374)", "p");
    legTruth->AddEntry(gTruthC3,     "MC Compton 3 (CE=4.864)", "p");
    legTruth->Draw();

    cTruth->Write("McTruth_PeakLocation_vs_Fraction");

    // =======================================================
    // PLOT 4: Fitted Peak Location (MeV) vs Crystal Number
    // =======================================================
    if (truthfit && truthfit->GetEntries() > 0) {
        Int_t   tf_cryNum;
        Float_t tf_peak_mev;

        truthfit->SetBranchAddress("crystalNo", &tf_cryNum);
        truthfit->SetBranchAddress("peak_mev",  &tf_peak_mev);
        
        Float_t tf_fr_full, tf_fr_1st, tf_fr_2nd, tf_fr_c1, tf_fr_c2, tf_fr_c3,tf_fr_ebk;
        // THIS IS CORRECT - The strings match your MakeAnalysisTree.cc branch names
		truthfit->SetBranchAddress("fr_full", &tf_fr_full);
		truthfit->SetBranchAddress("fr_1st",  &tf_fr_1st);
		truthfit->SetBranchAddress("fr_2nd",  &tf_fr_2nd);
		truthfit->SetBranchAddress("fr_c1",   &tf_fr_c1);
		truthfit->SetBranchAddress("fr_c2",   &tf_fr_c2);
		truthfit->SetBranchAddress("fr_c3",   &tf_fr_c3);
		truthfit->SetBranchAddress("fr_ebk", &tf_fr_ebk);
        
		
		double sumFitFull = 0, sumFit1st = 0, sumFit2nd = 0;
        double sumFitC1 = 0, sumFitC2 = 0, sumFitC3 = 0, sumFitebk = 0;
        
        std::vector<double> cry_vec, peak_vec;
        Long64_t Nfit = truthfit->GetEntries();
        for (Long64_t i = 0; i < Nfit; i++) {
            truthfit->GetEntry(i);
            cry_vec.push_back(static_cast<double>(tf_cryNum));
            peak_vec.push_back(static_cast<double>(tf_peak_mev));
            
            // Sum up the fractions
            sumFitFull += tf_fr_full;
            sumFit1st  += tf_fr_1st;
            sumFit2nd  += tf_fr_2nd;
            sumFitC1   += tf_fr_c1;
            sumFitC2   += tf_fr_c2;
            sumFitC3   += tf_fr_c3;
            sumFitebk  += tf_fr_ebk;
        }
        // Calculate Averages for the printout
        double avgFitFull = sumFitFull / Nfit;
        double avgFit1st  = sumFit1st  / Nfit;
        double avgFit2nd  = sumFit2nd  / Nfit;
        double avgFitC1   = sumFitC1   / Nfit;
        double avgFitC2   = sumFitC2   / Nfit;
        double avgFitC3   = sumFitC3   / Nfit;
        double avgFitebk   = sumFitebk   / Nfit;
		if(useCompton){
        std::cout << std::fixed << std::setprecision(7);
        std::cout << "\n===== SHAPE-BASED (FIT) GLOBAL MC TRUTH SUMMARY =====\n";
        std::cout << "Total Crystals Fit:          " << Nfit << "\n";
        std::cout << "Average Main peak fraction:  " << avgFitFull << "\n";
        std::cout << "Average 1st escape fraction: " << avgFit1st << "\n";
        std::cout << "Average 2nd escape fraction: " << avgFit2nd << "\n";
        std::cout << "Average Compton1 fraction:   " << avgFitC1 << "\n";
        std::cout << "Average Compton2 fraction:   " << avgFitC2 << "\n";
        std::cout << "Average Compton3 fraction:   " << avgFitC3 << "\n";
        std::cout << "Sum of averages:             " << avgFitFull + avgFit1st + avgFit2nd + avgFitC1 + avgFitC2 + avgFitC3 << std::endl;
        std::cout << "=====================================================\n\n";
		}
		else{
		std::cout << std::fixed << std::setprecision(7);
        std::cout << "\n===== SHAPE-BASED (FIT) GLOBAL MC TRUTH SUMMARY =====\n";
        std::cout << "Total Crystals Fit:          " << Nfit << "\n";
        std::cout << "Average Main peak fraction:  " << avgFitFull << "\n";
        std::cout << "Average 1st escape fraction: " << avgFit1st << "\n";
        std::cout << "Average 2nd escape fraction: " << avgFit2nd << "\n";
        std::cout << "Average background:  		   "<< avgFitebk << "\n";
        std::cout << "Sum of averages:             " << avgFitFull + avgFit1st + avgFit2nd + avgFitebk << std::endl;
        std::cout << "=====================================================\n\n";
		}
        // Compute mean and std dev across all fitted peaks
        double sum_p = std::accumulate(peak_vec.begin(), peak_vec.end(), 0.0);
        double mean_p = sum_p / peak_vec.size();
        double sq_sum = 0.0;
        for (double v : peak_vec) sq_sum += (v - mean_p) * (v - mean_p);
        double stddev_p = (peak_vec.size() > 1) ? std::sqrt(sq_sum / peak_vec.size()) : 0.0;

        TGraph* grPeakMev = new TGraph(cry_vec.size(), cry_vec.data(), peak_vec.data());
        grPeakMev->SetTitle("MC Truth Fitted Peak Location;Crystal Number;Peak Energy [MeV]");
        grPeakMev->SetMarkerStyle(20);
        grPeakMev->SetMarkerSize(0.8);
        grPeakMev->SetMarkerColor(kBlue);
        grPeakMev->SetLineColor(kBlue);

        TCanvas* cPeakMev = new TCanvas("cPeakMev", "MC Truth Fitted Peak Location", 900, 600);
        grPeakMev->Draw("AP");
        grPeakMev->GetYaxis()->SetRangeUser(5.5, 6.5);

        double xlo = *std::min_element(cry_vec.begin(), cry_vec.end());
        double xhi = *std::max_element(cry_vec.begin(), cry_vec.end());

        TLine* refLine = new TLine(xlo, 6.13, xhi, 6.13);
        refLine->SetLineColor(kRed);
        refLine->SetLineStyle(kDashed);
        refLine->SetLineWidth(2);
        refLine->Draw();

        TLine* avgLine = new TLine(xlo, mean_p, xhi, mean_p);
        avgLine->SetLineColor(kCyan+1);
        avgLine->SetLineStyle(kSolid);
        avgLine->SetLineWidth(2);
        avgLine->Draw();

        TLine* sigUp = new TLine(xlo, mean_p + 2*stddev_p, xhi, mean_p + 2*stddev_p);
        sigUp->SetLineColor(kOrange+1);
        sigUp->SetLineStyle(kSolid);
        sigUp->SetLineWidth(2);
        sigUp->Draw();

        TLine* sigDn = new TLine(xlo, mean_p - 2*stddev_p, xhi, mean_p - 2*stddev_p);
        sigDn->SetLineColor(kOrange+1);
        sigDn->SetLineStyle(kSolid);
        sigDn->SetLineWidth(2);
        sigDn->Draw();

        TLegend* legPeak = new TLegend(0.50, 0.68, 0.92, 0.92);
        legPeak->SetBorderSize(0);
        legPeak->SetFillStyle(0);
        legPeak->AddEntry(grPeakMev, "Fitted peak (RooFit)",                           "p");
        legPeak->AddEntry(refLine,   "True peak 6.13 MeV",                             "l");
        legPeak->AddEntry(avgLine,   Form("Avg = %.3f MeV", mean_p),                   "l");
        legPeak->AddEntry(sigUp,     Form("#pm2#sigma (%.3f MeV)",2.0 * stddev_p),           "l");
        legPeak->Draw();
        TPaveText *ptL = new TPaveText(0.15, 0.75, 0.45, 0.60, "brNDC");
    	ptL->SetFillStyle(0); ptL->SetBorderSize(0); ptL->SetTextSize(0.035);
    	ptL->AddText(Form("Avg = %4.8f #pm %.8f", mean_p, stddev_p));
   	 	ptL->AddText(Form("Std Dev = %.2f%%", (stddev_p/mean_p)*100));
   		 ptL->Draw();

        cPeakMev->Write("McTruth_FittedPeak_vs_Crystal");
        grPeakMev->Write("grTruthPeakMev");
    }

    outFile->Close();
    std::cout << "[INFO] MC Truth summary histograms written.\n";
}

// =====================================================================
// FIT MC TRUTH SPECTRUM IN MEV SPACE
// =====================================================================
void mcinfo::FitMCTruthSpectrum(TH1F* hist_truth_mev, int crystalNo, TTree *truthfit_tree,
                                Float_t &truth_peak_mev, Float_t &truth_fr_full,
                                Float_t &truth_fr_1st, Float_t &truth_fr_2nd,
                                Float_t &truth_fr_c1, Float_t &truth_fr_c2, Float_t &truth_fr_c3,
                                Float_t &truth_fr_ebk,
                                Float_t &truth_width, Float_t &truth_alpha, Float_t &truth_chi2, Int_t &truth_ndof,
                                bool useCompton)
{
    std::cout << "======================================================================\n";
    std::cout << "Fitting MC Truth Spectrum (MeV) for crystal " << crystalNo << "\n";
    std::cout << "======================================================================\n";

    // --- Initialize parameters ---
    double initPeak_MeV   = 6.13;
    double peakLow_MeV    = 5.90;
    double peakHigh_MeV   = 6.25;
    //double initWidth_MeV  = 0.08;  // MC truth peaks are sharp
    //double initAlpha_MeV  = 1.2;   // smaller alpha = shorter CB tail, less background absorption

    // Initial fractions ordered to match physics: 1st escape largest, full peak smallest
    double initFrFull = 0.10;
    double initFr1st  = 0.21;
    double initFr2nd  = 0.17;
    double initFrC1   = 0.09;
    double initFrC2   = 0.13;
    double initFrC3   = 0.09;

    // --- Setup histogram with proper errors to avoid infinity errors ---
    hist_truth_mev->Sumw2();  // Enable errors for each bin (prevents chi2 INFINITY on empty bins)

    // --- Count events in fit range before defining yield parameters ---
    int fit_lo_init = hist_truth_mev->FindBin(2.5);
    int fit_hi_init = hist_truth_mev->FindBin(6.5);
    double nevents_init = hist_truth_mev->Integral(fit_lo_init, fit_hi_init);
    if (nevents_init <= 0) nevents_init = 1.0;

    // --- RooFit variables (MeV scale) ---
    RooRealVar energy("energy", "Reconstructed energy [MeV]", 2.5, 6.5);
    RooRealVar m_e("m_e", "electron mass [MeV]", 0.511);

    RooRealVar peak_full_mev("peak_full_mev", "Full peak [MeV]", initPeak_MeV, peakLow_MeV, peakHigh_MeV);
    RooRealVar width_mev("width_mev", "Width [MeV]", 0.01);
    width_mev.setConstant(kTRUE);
    // Left-side CB parameters (matching SourceFitter: alpha free, n free)
    //RooRealVar alpha_mev("alpha_mev", "Left CB alpha", initAlpha_MeV, 0.5,5.0);
    //RooRealVar n_deg("n_deg",         "Left CB n",     5.0, 5.0, 50.0);
    // Right-side CB parameters (matching SourceFitter: alpha fixed at 3, n free)
    //RooRealVar alpha_R("alpha_R", "Right CB alpha", 1.2,0.5,5.0);
    //RooRealVar n_R("n_R",         "Right CB n",     5.0, 1.0, 50.0);
    //unique alpha for mc info
    //RooRealVar alpha_full("alpha_full", "Left CB alpha for full", initAlpha_MeV, 0.5,5.0);
    //RooRealVar n_deg_full("n_deg_full",  "Left CB n for full",     10);

    // --- Define escape peak positions (derived) ---
    RooFormulaVar peak_1st_mev("peak_1st_mev", "First escape [MeV]",
                               "peak_full_mev - m_e", RooArgSet(peak_full_mev, m_e));
    RooFormulaVar peak_2nd_mev("peak_2nd_mev", "Second escape [MeV]",
                               "peak_full_mev - 2*m_e", RooArgSet(peak_full_mev, m_e));

    // --- Double-sided Crystal Ball PDFs (matching SourceFitter model) ---
    //RooCrystalBall full_cb("full_cb", "Full peak CB",   energy, peak_full_mev, width_mev, alpha_mev, n_deg, alpha_R, n_R);
    //RooCBShape full_cb("full_cb", "Full peak CB",   energy, peak_full_mev, width_mev, alpha_full, n_deg_full);
    //RooCrystalBall esc1_cb("esc1_cb", "1st escape CB",  energy, peak_1st_mev,  width_mev, alpha_mev, n_deg, alpha_R, n_R);
    //RooCrystalBall esc2_cb("esc2_cb", "2nd escape CB",  energy, peak_2nd_mev,  width_mev, alpha_mev, n_deg, alpha_R, n_R);
    RooGaussian full_cb("full_cb", "Full peak Gauss",    energy, peak_full_mev, 0.09);
    RooGaussian esc1_cb("esc1_cb", "1st escape Gauss",  energy, peak_1st_mev,  0.212);
    RooGaussian esc2_cb("esc2_cb", "2nd escape Gauss",  energy, peak_2nd_mev,  0.166);

    // --- CB peak yields (always needed) ---
    RooRealVar yield_full("yield_full", "Full peak yield",  initFrFull * nevents_init, 0.0, nevents_init);
    RooRealVar yield_1st ("yield_1st",  "1st escape yield", initFr1st  * nevents_init, 0.0, nevents_init);
    RooRealVar yield_2nd ("yield_2nd",  "2nd escape yield", initFr2nd  * nevents_init, 0.0, nevents_init);

    // --- Compton erfc background (useCompton=true) ---
    RooRealVar CE1_rv("CE1_rv", "Compton edge 1 [MeV]", 5.8842);
    RooRealVar CE2_rv("CE2_rv", "Compton edge 2 [MeV]", 5.3741);
    RooRealVar CE3_rv("CE3_rv", "Compton edge 3 [MeV]", 4.8640);
    CE1_rv.setConstant(kTRUE);
    CE2_rv.setConstant(kTRUE);
    CE3_rv.setConstant(kTRUE);
    TString erfc_form = "0.5 * TMath::Erfc((@0 - @1) / (TMath::Sqrt(2) * @2))";
    // RooArgList preserves insertion order so @0=energy, @1=CE_rv, @2=width_mev
    RooGenericPdf compton1_pdf("compton1_pdf", "Compton edge 1", erfc_form,
                               RooArgList(energy, CE1_rv, width_mev));
    RooGenericPdf compton2_pdf("compton2_pdf", "Compton edge 2", erfc_form,
                               RooArgList(energy, CE2_rv, width_mev));
    RooGenericPdf compton3_pdf("compton3_pdf", "Compton edge 3", erfc_form,
                               RooArgList(energy, CE3_rv, width_mev));
    // Compton yields tied to CB peak yields: yield_ci = ratio_ci * yield_cb_i
    // Physically: Compton background from each photon scales with that photon's flux (= CB yield).
    // This breaks the degeneracy between the 3 overlapping erfc steps.
    /*double initRatioC1 = (initFrFull > 0) ? initFrC1 / initFrFull : 0.5;
    double initRatioC2 = (initFr1st  > 0) ? initFrC2 / initFr1st  : 0.5;
    double initRatioC3 = (initFr2nd  > 0) ? initFrC3 / initFr2nd  : 0.5;
    RooRealVar ratio_c1("ratio_c1", "Compton1/full ratio", initRatioC1, 0.0, 10.0);
    RooRealVar ratio_c2("ratio_c2", "Compton2/1st ratio",  initRatioC2, 0.0, 10.0);
    RooRealVar ratio_c3("ratio_c3", "Compton3/2nd ratio",  initRatioC3, 0.0, 10.0);*/
    /*RooFormulaVar yield_c1_fv("yield_c1_fv", "Compton1 yield", "@0*@1", RooArgList(yield_full, ratio_c1));
    RooFormulaVar yield_c2_fv("yield_c2_fv", "Compton2 yield", "@0*@1", RooArgList(yield_1st,  ratio_c2));
    RooFormulaVar yield_c3_fv("yield_c3_fv", "Compton3 yield", "@0*@1", RooArgList(yield_2nd,  ratio_c3));*/
   /* // Define cumulative flux variables
RooFormulaVar flux_for_C2("flux_for_C2", "yield_full + yield_1st", RooArgList(yield_full, yield_1st));
RooFormulaVar flux_for_C3("flux_for_C3", "yield_full + yield_1st + yield_2nd", RooArgList(yield_full, yield_1st, yield_2nd));

// Tie ratios to the available flux
RooFormulaVar yield_c1_fv("yield_c1_fv", "@0*@1", RooArgList(yield_full, ratio_c1));
RooFormulaVar yield_c2_fv("yield_c2_fv", "@0*@1", RooArgList(flux_for_C2, ratio_c2));
RooFormulaVar yield_c3_fv("yield_c3_fv", "@0*@1", RooArgList(flux_for_C3, ratio_c3));*/
	RooRealVar yield_c1_fv("yield_c1_fv", "Compton 1 yield", initFrC1 * nevents_init, 0.02*nevents_init, 0.2*nevents_init);
    RooRealVar yield_c2_fv("yield_c2_fv", "Compton 2 yield", initFrC2 * nevents_init, 0.02*nevents_init, 0.2*nevents_init);
    RooRealVar yield_c3_fv("yield_c3_fv", "Compton 3 yield", initFrC3 * nevents_init, 0.02*nevents_init, 0.2*nevents_init);

    // --- Exponential background (useCompton=false) ---
    RooRealVar exp_slope("exp_slope", "Exp slope", -0.5, -5.0, -0.01);
    TString exp_form = "TMath::Exp(@1 * @0)";
    RooGenericPdf bkg_exp("bkg_exp", "Exponential background", exp_form,
                          RooArgList(energy, exp_slope));
    double initFrExp = initFrC1 + initFrC2 + initFrC3;
    RooRealVar yield_exp("yield_exp", "Exponential yield", initFrExp * nevents_init, 0.0, nevents_init);

    // --- Build combined model (toggle between Compton erfc and exponential background) ---
    std::cout << "[INFO] Background model: " << (useCompton ? "Compton erfc (yields tied to CB peaks)" : "Exponential") << "\n";
    RooAddPdf* model_ptr = nullptr;
    if (useCompton) {
        model_ptr = new RooAddPdf("model_truth", "MC Truth Model (Compton)",
            RooArgList(full_cb, esc1_cb, esc2_cb, compton1_pdf, compton2_pdf, compton3_pdf),
            RooArgList(yield_full, yield_1st, yield_2nd, yield_c1_fv, yield_c2_fv, yield_c3_fv));
    } else {
        model_ptr = new RooAddPdf("model_truth", "MC Truth Model (Exponential)",
            RooArgList(full_cb, esc1_cb, esc2_cb, bkg_exp),
            RooArgList(yield_full, yield_1st, yield_2nd, yield_exp));
    }
    RooAddPdf& model_truth = *model_ptr;

    // --- Fit the data ---
    int fit_lo = hist_truth_mev->FindBin(2.5);
    int fit_hi = hist_truth_mev->FindBin(6.5);
    int nbins = fit_hi - fit_lo + 1;
    int nevents = hist_truth_mev->Integral(fit_lo, fit_hi);

    std::cout << "[DEBUG] Fitting " << nevents << " events over " << nbins << " bins\n";

    RooDataHist datahist_truth("datahist_truth", "MC Truth Data", energy, hist_truth_mev);

    // Extended maximum likelihood fit — proper API for RooAddPdf with yield parameters
    RooFitResult* fitres_truth = model_truth.fitTo(datahist_truth,
        RooFit::Extended(kTRUE),
        RooFit::SumW2Error(kTRUE),
        RooFit::Save(kTRUE),
        RooFit::Minimizer("Minuit2", "Migrad"),
        RooFit::Strategy(1),
        RooFit::PrintLevel(-1),
        RooFit::Hesse(kTRUE)
    );

    // --- Extract results ---
    truth_peak_mev = peak_full_mev.getVal();
    truth_width    = width_mev.getVal();
    truth_alpha    = 0;//alpha_mev.getVal();

    // Fractions (normalize yields to 1)
    double yf  = yield_full.getVal();
    double y1  = yield_1st.getVal();
    double y2  = yield_2nd.getVal();
    double yc1, yc2, yc3, yebk;
    if (useCompton) {
        yc1 = yield_c1_fv.getVal();
        yc2 = yield_c2_fv.getVal();
        yc3 = yield_c3_fv.getVal();
        yebk = 0.0;
    } else {
        yc1 = 0.0;
        yc2 = 0.0;
        yc3 = 0.0;
        yebk = yield_exp.getVal(); // Exponential yield goes here
    }
    double ytot = yf + y1 + y2 + yc1 + yc2 + yc3 + yebk;

    truth_fr_full = (ytot > 0) ? yf  / ytot : 0;
    truth_fr_1st  = (ytot > 0) ? y1  / ytot : 0;
    truth_fr_2nd  = (ytot > 0) ? y2  / ytot : 0;
    truth_fr_c1   = (ytot > 0) ? yc1 / ytot : 0;
    truth_fr_c2   = (ytot > 0) ? yc2 / ytot : 0;
    truth_fr_c3   = (ytot > 0) ? yc3 / ytot : 0;
    truth_fr_ebk  = (ytot > 0) ? yebk / ytot : 0;

    // Chi-square at best-fit point (evaluate after fit so parameters are at converged values)
    //energy.setRange("eval_range", 2.5, 6.5);
    RooAbsReal* chi2_var = model_truth.createChi2(datahist_truth,
        RooFit::DataError(RooAbsData::Poisson),RooFit::Range(2.5, 6.5));
    truth_chi2 = static_cast<float>(chi2_var->getVal());
   /* ///
    // 3. Debugging prints (BEFORE deleting the pointer)
    std::cout << "chi2: " << truth_chi2 << "\n";
    std::cout << "Data entries: " << datahist_truth.sumEntries("", "eval_range") << std::endl;
    
    // Check if the model yield matches data yield
    double total_model_yield = yield_full.getVal() + yield_1st.getVal() + yield_2nd.getVal();
    if (useCompton) total_model_yield += (yield_c1.getVal() + yield_c2.getVal() + yield_c3.getVal());
    else total_model_yield += yield_exp.getVal();
    
    std::cout << "Model Yield: " << total_model_yield << " vs Data: " << nevents << std::endl;

    chi2_var->Print("v");
	///*/
    delete chi2_var;

    // peak+width+alpha + 3 CB yields + 3 Compton ratios = 9; peak+width+alpha + 3 CB yields + slope + exp yield = 8
    int nfree = useCompton ? 9 : 8;
    truth_ndof = nbins - nfree;
    std::cout << "chi2" << truth_chi2  << "\n";

    std::cout << "[RESULT] Peak (MeV): " << truth_peak_mev << " | Width: " << truth_width
              << " | Alpha: " << truth_alpha << "\n";
    std::cout << "[RESULT] Fractions: Full=" << truth_fr_full << " | 1st=" << truth_fr_1st
              << " | 2nd=" << truth_fr_2nd << "\n";
    if (useCompton) {
        std::cout << "[RESULT] Compton fractions: C1=" << truth_fr_c1
                  << " | C2=" << truth_fr_c2 << " | C3=" << truth_fr_c3 << "\n";
       // std::cout << "[RESULT] Compton/CB ratios: r1(C1/full)=" << ratio_c1.getVal()
        //          << " | r2(C2/1st)=" << ratio_c2.getVal()
         //         << " | r3(C3/2nd)=" << ratio_c3.getVal() << "\n";
    } else {
        std::cout << "[RESULT] Exp bkg fraction=" << truth_fr_c1
                  << " | slope=" << exp_slope.getVal() << "\n";
    }
    std::cout << "[RESULT] chi2/ndof = " << truth_chi2 << " / " << truth_ndof
              << " = " << (truth_ndof > 0 ? truth_chi2/truth_ndof : -1) << "\n";

    // --- Create plot and save to .root file ---
    gStyle->SetOptFit(1111);
    gStyle->SetOptStat(0);
    gStyle->SetPadBottomMargin(0.125);
    gStyle->SetPadTopMargin(0.075);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetTitleOffset(1.0, "x");
    gStyle->SetTitleOffset(1.75, "y");

    TCanvas* can_truth = new TCanvas("can_truth", "MC Truth Fit", 100, 100, 800, 600);
    can_truth->Draw();

    RooPlot* frame = energy.frame(RooFit::Title(Form("MC Truth Fit - Crystal %d", crystalNo)));
    datahist_truth.plotOn(frame, RooFit::MarkerColor(kBlack), RooFit::MarkerSize(0.7), RooFit::Name("data"));
    model_truth.plotOn(frame, RooFit::LineColor(kRed),      RooFit::LineWidth(2),                        RooFit::Name("model"));
    model_truth.plotOn(frame, RooFit::Components(full_cb),  RooFit::LineColor(kOrange),  RooFit::LineStyle(kDashed), RooFit::Name("full"));
    model_truth.plotOn(frame, RooFit::Components(esc1_cb),  RooFit::LineColor(kViolet),  RooFit::LineStyle(kDashed), RooFit::Name("1st_esc"));
    model_truth.plotOn(frame, RooFit::Components(esc2_cb),  RooFit::LineColor(kCyan),    RooFit::LineStyle(kDashed), RooFit::Name("2nd_esc"));
    if (useCompton) {
        model_truth.plotOn(frame, RooFit::Components(compton1_pdf), RooFit::LineColor(kGreen),   RooFit::LineStyle(kDashed), RooFit::Name("c1"));
        model_truth.plotOn(frame, RooFit::Components(compton2_pdf), RooFit::LineColor(kGreen+2), RooFit::LineStyle(kDashed), RooFit::Name("c2"));
        model_truth.plotOn(frame, RooFit::Components(compton3_pdf), RooFit::LineColor(kGreen+4), RooFit::LineStyle(kDashed), RooFit::Name("c3"));
    } else {
        model_truth.plotOn(frame, RooFit::Components(bkg_exp), RooFit::LineColor(kOrange+7), RooFit::LineStyle(kDashed), RooFit::Name("bkg"));
    }

    frame->Draw();

    // Add legend
    TLegend* leg = new TLegend(0.60, 0.50, 0.95, 0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry("data",    "MC Truth Data",          "ep");
    leg->AddEntry("model",   "Total Fit",              "l");
    leg->AddEntry("full",    "Full Peak (6.13 MeV)",   "l");
    leg->AddEntry("1st_esc", "1st Escape (5.619 MeV)", "l");
    leg->AddEntry("2nd_esc", "2nd Escape (5.108 MeV)", "l");
    if (useCompton) {
        leg->AddEntry("c1",  "Compton 1 (CE=5.884)",   "l");
        leg->AddEntry("c2",  "Compton 2 (CE=5.374)",   "l");
        leg->AddEntry("c3",  "Compton 3 (CE=4.864)",   "l");
    } else {
        leg->AddEntry("bkg", "Exponential background", "l");
    }
    leg->Draw();

    TString rootFileName = Form("MC_Truth_Fit_Crystal_%d.root", crystalNo);
    can_truth->SaveAs(rootFileName);

    std::cout << "[INFO] Saved MC Truth fit plot to " << rootFileName << "\n";

    truthfit_tree->Fill();

    delete can_truth;
    delete model_ptr;
    if (fitres_truth) delete fitres_truth;
}
