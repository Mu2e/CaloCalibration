#include "CaloCalibration/SourceCalib/inc/mcinfo.hh"

using namespace std::chrono;
using namespace CaloSourceCalib;

void mcinfo::RunMCTruth(TH1F* hist, int cryNum, int disk, TTree *trueinfo, 
                        Int_t &cryNumparam, Int_t &tot_evts, 
                        Int_t &mainpeak, Int_t &first_espeak, Int_t &second_espeak, 
                        Int_t &background, Float_t &frmainpeak, 
                        Float_t &frfirst_espeak, Float_t &frsecond_espeak, Float_t &frbackground)
{
    // 1. Define Energy Constants (MeV)
    // Using constexpr ensures these are known at compile time and efficient.
    static constexpr double E_main_MeV    = 6.13;
    static constexpr double E_1st_MeV     = 5.619;
    static constexpr double E_2nd_MeV     = 5.108;
    static constexpr double E_compton_MeV = 0.511;

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
    int bin_2nd_low  = ax->FindBin(E_2nd - E_compton);; 
    int bin_2nd_high = bin_2nd_center;

    // 1st Escape Peak Window: [ E_2nd + 1 bin, E_1st ]
    int bin_1st_low  = bin_2nd_high + 1; 
    int bin_1st_high = bin_1st_center;

    // Main Peak Window: [ E_1st + 1 bin, E_main ]
    int bin_main_low  = bin_1st_high + 1;
    int bin_main_high = bin_main_center;

    // 4. Calculate Integrals
    mainpeak      = hist->Integral(bin_main_low, bin_main_high);    // 5.619 -> 6.13
    first_espeak  = hist->Integral(bin_1st_low, bin_1st_high);      // 5.108 -> 5.619
    second_espeak = hist->Integral(bin_2nd_low, bin_2nd_high);      // Compton -> 5.108

    tot_evts      = hist->Integral(bin_lo, bin_hi);
    double tot = (tot_evts > 0) ? static_cast<double>(tot_evts) : 1.0; // Avoid div by zero

    frmainpeak      = mainpeak / tot;
    frfirst_espeak  = first_espeak / tot;
    frsecond_espeak = second_espeak / tot;

    background      = tot_evts - (mainpeak + first_espeak + second_espeak);
    frbackground    = background / tot;

    trueinfo->Fill();
}

void mcinfo::FinalizeMCSummary(TTree* trueinfo)
{
    std::cout << "[INFO] Building global MC truth histograms...\n";

    // --- 1. Initialize Histograms ---
    TH1F* h_tot      = new TH1F("h_tot_events",  "Total MC events per crystal", 200, 0, 40000);
    
    // Unweighted histograms
    TH1F* h_frMain   = new TH1F("h_frMain",      "Fraction Main Peak",          100, 0, 1);
    TH1F* h_fr1st    = new TH1F("h_fr1st",       "Fraction 1st Escape",         100, 0, 1);
    TH1F* h_fr2nd    = new TH1F("h_fr2nd",       "Fraction 2nd Escape",         100, 0, 1);
    TH1F* h_frBg     = new TH1F("h_frBg",        "Fraction Background",         100, 0, 1);

    // Weighted histograms (weighted by total events)
    TH1F* h_frMain_w = new TH1F("h_frMain_w", "Weighted Main Peak Fraction", 100, 0, 1);
    TH1F* h_fr1st_w  = new TH1F("h_fr1st_w",  "Weighted 1st Escape Fraction",100, 0, 1);
    TH1F* h_fr2nd_w  = new TH1F("h_fr2nd_w",  "Weighted 2nd Escape Fraction",100, 0, 1);

    // --- 2. Process TTree Data ---
    Int_t   tot_evts;
    Float_t frmainpeak, frfirst_espeak, frsecond_espeak, frbackground;

    trueinfo->SetBranchAddress("tot_evts",         &tot_evts);
    trueinfo->SetBranchAddress("frmainpeak",       &frmainpeak);
    trueinfo->SetBranchAddress("frfirst_espeak",   &frfirst_espeak);
    trueinfo->SetBranchAddress("frsecond_espeak",  &frsecond_espeak);
    trueinfo->SetBranchAddress("frbackground",     &frbackground);

    Long64_t sumTot = 0, sumMain = 0, sum1st = 0, sum2nd = 0, sumBg = 0;
    Long64_t N = trueinfo->GetEntries();

    for (Long64_t i = 0; i < N; i++) {
        trueinfo->GetEntry(i);
        
        sumTot  += tot_evts;
        sumMain += frmainpeak * tot_evts;
        sum1st  += frfirst_espeak * tot_evts;
        sum2nd  += frsecond_espeak * tot_evts;
        sumBg   += frbackground * tot_evts;

        h_tot->Fill(tot_evts);
        
        // Fill Unweighted
        h_frMain->Fill(frmainpeak);
        h_fr1st->Fill(frfirst_espeak);
        h_fr2nd->Fill(frsecond_espeak);
        h_frBg->Fill(frbackground);
        
        // Fill Weighted
        h_frMain_w->Fill(frmainpeak, tot_evts);
        h_fr1st_w ->Fill(frfirst_espeak, tot_evts);
        h_fr2nd_w ->Fill(frsecond_espeak, tot_evts);
    }

    // --- 3. Calculate Global Averages ---
    double fMain = (sumTot > 0) ? double(sumMain) / sumTot : 0;
    double f1st  = (sumTot > 0) ? double(sum1st)  / sumTot : 0;
    double f2nd  = (sumTot > 0) ? double(sum2nd)  / sumTot : 0;
    double fBg   = (sumTot > 0) ? double(sumBg)   / sumTot : 0;

    // Normalize weighted histograms to integral 1 for comparison
    if(h_frMain_w->Integral() > 0) h_frMain_w->Scale(1.0 / h_frMain_w->Integral());
    if(h_fr1st_w->Integral() > 0)  h_fr1st_w->Scale(1.0 / h_fr1st_w->Integral());
    if(h_fr2nd_w->Integral() > 0)  h_fr2nd_w->Scale(1.0 / h_fr2nd_w->Integral());

    std::cout << std::fixed << std::setprecision(7);
    std::cout << "\n===== GLOBAL MC TRUTH SUMMARY =====\n";
    std::cout << "Total events:        " << sumTot << "\n";
    std::cout << "Main peak fraction:  " << fMain << "\n";
    std::cout << "1st escape fraction: " << f1st << "\n";
    std::cout << "2nd escape fraction: " << f2nd << "\n";
    std::cout << "Background fraction: " << fBg  << "\n";
    std::cout << "Sum of fractions:    " << fBg+fMain+f1st+f2nd << std::endl;
    std::cout << "==================================\n\n";

    TFile* outFile = new TFile("MC_Truth_Summary.root", "RECREATE");

    // =======================================================
    // PLOT 1 & 2: Overlay Histograms (Unweighted & Weighted)
    // =======================================================

    // Helper Lambda to draw the overlay plots (Reduces redundancy)
    auto draw_overlay = [&](const char* name, const char* title, 
                            TH1F* hM, TH1F* h1, TH1F* h2, const char* leg_suffix) 
    {
        TCanvas* c = new TCanvas(name, title, 800, 600);
        
        // Style inputs
        hM->SetLineColor(kBlue);        hM->SetLineWidth(2);
        h1->SetLineColor(kGreen + 2);   h1->SetLineWidth(2);
        h2->SetLineColor(kMagenta + 1); h2->SetLineWidth(2);

        // Calculate Y Range to fit all
        double max_val = std::max({ hM->GetMaximum(), h1->GetMaximum(), h2->GetMaximum() });
        hM->SetMaximum(1.1 * max_val);
        hM->SetMinimum(0);

        hM->Draw("HIST");
        h1->Draw("HIST SAME");
        h2->Draw("HIST SAME");

        TLegend* lg = new TLegend(0.65, 0.7, 0.9, 0.9);
        lg->AddEntry(hM, Form("Main Peak %s", leg_suffix), "l");
        lg->AddEntry(h1, Form("1st Escape %s", leg_suffix), "l");
        lg->AddEntry(h2, Form("2nd Escape %s", leg_suffix), "l");
        lg->Draw();

        c->Write(name);
        return c; // Return if needed, though we Write() immediately
    };

    draw_overlay("MC_Fractions_Overlay", "MC Truth Fraction Overlay", 
                 h_frMain, h_fr1st, h_fr2nd, "");
                 
    draw_overlay("MC_Fractions_Overlay_weighted", "MC Truth Fraction Overlay (Weighted)", 
                 h_frMain_w, h_fr1st_w, h_fr2nd_w, "(wt)");

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

    // Constants for positions
    const double E_main = 6.13;
    const double E_1st  = 5.619;
    const double E_2nd  = 5.108;

    // Create graphs using helper
    TGraph* gTruthMain   = create_truth_graph(E_main, fMain, 24, kBlue);        // Open circle
    TGraph* gTruthFirst  = create_truth_graph(E_1st,  f1st,  25, kGreen + 2);   // Open square
    TGraph* gTruthSecond = create_truth_graph(E_2nd,  f2nd,  26, kMagenta + 1); // Open triangle

    TCanvas* cTruth = new TCanvas("cTruth", "MC Truth Peak Location vs Fraction", 900, 650);
    cTruth->DrawFrame(2.5, 0.0, 8.0, 1.05, "MC Truth: Peak Location vs Fraction;Peak Energy [MeV];Fraction of Events");

    gTruthMain->Draw("P SAME");
    gTruthFirst->Draw("P SAME");
    gTruthSecond->Draw("P SAME");

    TLegend* legTruth = new TLegend(0.65, 0.7, 0.9, 0.9);
    legTruth->AddEntry(gTruthMain,   "MC Main Peak",  "p");
    legTruth->AddEntry(gTruthFirst,  "MC 1st Escape", "p");
    legTruth->AddEntry(gTruthSecond, "MC 2nd Escape", "p");
    legTruth->Draw();

    cTruth->Write("McTruth_PeakLocation_vs_Fraction");

    outFile->Close();
    std::cout << "[INFO] MC Truth summary histograms written.\n";
}
