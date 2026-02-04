#include "CaloCalibration/SourceCalib/inc/MakeAnalysisTree.hh"
#include "CaloCalibration/SourceCalib/inc/SourceFitter.hh"
#include "CaloCalibration/SourceCalib/inc/SourcePlotter.hh"
#include "CaloCalibration/SourceCalib/inc/2dcontour.hh"
#include "CaloCalibration/SourceCalib/inc/mcinfo.hh"// extracting mc information -- need to be removed for real data
#include <chrono>
using namespace std::chrono;

using namespace CaloSourceCalib;

//method with one big root file
TString filepath_disk0 = "/exp/mu2e/app/home/mu2epro/sourcecalib/full_containment/newsinglephoton.root";
//"/exp/mu2e/app/home/mu2epro/sourcecalib/full_containment/newsinglephoton.root";
//"/exp/mu2e/app/home/mu2epro/sourcecalib/full_containment/reconew.root";
//"/exp/mu2e/app/home/mu2epro/sourcecalib/disk0/SourceAna1e9.root";
//"/exp/mu2e/app/home/mu2epro/sourcecalib/full_containment/singlephoton.root";
//"/exp/mu2e/app/home/mu2epro/sourcecalib/full_containment/multiphoton.root";

TString filepath_disk1 = "/exp/mu2e/app/home/mu2epro/sourcecalib/disk1/nts.mu2e.SourceCalibAna1e9.0.root";

//TString filepath_disk1 = "/exp/mu2e/app/users/hjafree/SourceFitDir/combined_disk1.root";

std::pair<TH1F*, TFile*> get_data_histogram(int cryNum, int disk) {
    TString filepath;
    TString histPath;
    // Currently the disk can be toggled between sipms (disk 0 or 1) vs crystals (disk 2 or 3-- representing 0 and 1 respectively) as input. This is a temporary switch being used for studies and will be removed before the final data run.
    if (disk == 0 || disk == 2) {
        filepath = filepath_disk0;
        //if "mc" flag is called
        histPath = (disk == 0) ? "SourceAna/sipm_ADC/sipm_" : "SourceAna/crystals_ADC/cry_";
        //uncomment line below and comment the one above if cut and count needs to be applied to the mc truth histograms
        //histPath = (disk == 0) ? "SourceAna/crystals_edep_truth/cry_" : "SourceAna/crystals_ADC/cry_";
    } else if (disk == 1 || disk == 3) {
        filepath = filepath_disk1;
        histPath = (disk == 1) ? "SourceAna/sipm_ADC/sipm_" : "SourceAna/crystals_ADC/cry_";
    }
 
    TFile *f = new TFile(filepath);
    TString crystalNumber = to_string(cryNum);
    TH1F* hist = (TH1F*)f->Get(histPath + crystalNumber); 
		hist->SetDirectory(0);  // Detach from file
    return std::make_pair(hist, f);
}
std::pair<double, double> ComputeHistogramStats(TH1F* hist) {
    double sum = 0;
    double weightedSum = 0;
    double weightedSumSq = 0;

    int nBins = hist->GetNbinsX();
    for (int i = 1; i <= nBins; ++i) {
        double content = hist->GetBinContent(i);
        double center = hist->GetBinCenter(i);

        sum += content;
        weightedSum += content * center;
        weightedSumSq += content * center * center;
    }

    double mean = (sum > 0) ? weightedSum / sum : 0;
    double variance = (sum > 0) ? (weightedSumSq / sum) - (mean * mean) : 0;
    double stddev = (variance > 0) ? std::sqrt(variance) : 0;

    return std::make_pair(mean, stddev);
}

/*main function allows a loop over all crystals or a choice of a single crystal*/
int main(int argc, char* argv[]) {
    std::cout << "========== Welcome to the Mu2e Source Calibration Analysis ==========" << std::endl;

    // 1. Safety Check: Ensure we have at least the 4 mandatory arguments
    // Usage: ./exe <start> <end> <alg> <disk> [optional flags...]
    if (argc < 5) {
        std::cerr << "[ERROR] Missing arguments.\n";
        std::cerr << "Usage: " << argv[0] << " <start_cry> <end_cry> <alg> <disk> [flags: overlay, contour, mc]\n";
        return 1;
    }

    // 2. Parse Fixed Positional Arguments
    int anacrys_start = std::atoi(argv[1]);
    int anacrys_end   = std::atoi(argv[2]);
    TString alg       = argv[3];
    int disk          = std::atoi(argv[4]);

    // 3. Initialize Flags
    bool doOverlay = false;
    bool contour   = false;
    bool isMC      = false;

    // 4. Parse Optional Flags
    // We start at i = 5 because indices 1, 2, 3, 4 are already used above.
    for (int i = 5; i < argc; i++) {
        TString arg = argv[i];
        arg.ToLower(); // Optional: makes it case-insensitive
        
        if (arg.Contains("overlay")) doOverlay = true;
        if (arg.Contains("contour")) contour   = true;
        if (arg.Contains("mc"))      isMC      = true;
    }
// --------------------------------------------------------
// MC TRUTH BLOCK
// --------------------------------------------------------
if (isMC) {
    std::cout << "[INFO] Running in MC mode, using MC file paths.\n";
    std::cout << "=======================================\n";
    std::cout << "        Running MC Truth Analysis        \n";
    std::cout << "=======================================\n";

    // 1. Create File and Tree INSIDE the block (so we don't create junk files for NLL fits)
    TFile *truthtable = new TFile("truthinfo.root", "RECREATE");
    TTree *trueinfo = new TTree("trueinfo", "MC Truth Information");

    // 2. Define Branch Variables
    Int_t cryNumparam, tot_evts, mainpeak, first_espeak, second_espeak, background;
    Float_t frmainpeak, frfirst_espeak, frsecond_espeak, frbackground;

    // 3. Link Branches
    trueinfo->Branch("cryNumparam",     &cryNumparam,     "cryNumparam/I");
    trueinfo->Branch("tot_evts",        &tot_evts,        "tot_evts/I");
    trueinfo->Branch("mainpeak",        &mainpeak,        "mainpeak/I");
    trueinfo->Branch("first_espeak",    &first_espeak,    "first_espeak/I");
    trueinfo->Branch("second_espeak",   &second_espeak,   "second_espeak/I");
    trueinfo->Branch("background",      &background,      "background/I");
    trueinfo->Branch("frmainpeak",      &frmainpeak,      "frmainpeak/F");
    trueinfo->Branch("frfirst_espeak",  &frfirst_espeak,  "frfirst_espeak/F");
    trueinfo->Branch("frsecond_espeak", &frsecond_espeak, "frsecond_espeak/F");
    trueinfo->Branch("frbackground",    &frbackground,    "frbackground/F");

    // 4. Loop Over Crystals
    for (int cryNum = anacrys_start; cryNum < anacrys_end; cryNum++) {
        
        // Helper to load histogram (make sure this returns the right MC histogram)
        auto [hist, file] = get_data_histogram(cryNum, disk);
        if (!hist) { 
            if (file) { file->Close(); delete file; }
            continue; 
        }
        hist->SetDirectory(0); // Detach from file so we can close the input file

        // Explicitly set the crystal number for the branch
        cryNumparam = cryNum; 

        // Run Logic
        mcinfo truth;
        truth.RunMCTruth(
            hist, cryNum, disk, trueinfo, 
            cryNumparam, tot_evts, mainpeak, first_espeak, second_espeak, background,
            frmainpeak, frfirst_espeak, frsecond_espeak, frbackground
        );

        // Cleanup Input
        file->Close();
        delete file;
        delete hist;
    }

    // 5. Finalize and Save
    mcinfo summary;
    summary.FinalizeMCSummary(trueinfo);

    truthtable->cd();
    trueinfo->Write(); // Write the tree explicitly
    truthtable->Write(); // Write any histograms created by FinalizeMCSummary
    truthtable->Close();
    
    std::cout << "[INFO] MC Truth analysis complete. Saved to truthinfo.root" << std::endl;
   
     return 0; 
}

// --------------------------------------------------------
// FITTING BLOCK
// --------------------------------------------------------
  TFile *table = new TFile("arXivTable.root", "RECREATE");
  Int_t nEvents;
  Int_t convergencestatus;
  Float_t fpeak, peakerrorhigh,peakerrorlo,fsigma, chiSq, fstpeak,
  scdpeak,fcbalphaparam,fcbndegparam,comCnstparam,
  combetaparam,frFullparam,frFrstparam,frScndparam,crystalNoparam,frBKGparam,
  pval,h_means,h_stddevs,unreducedchi2,fval,mparam,etaparam,widtherrorhigh,widtherrorlo,errbarhigh, errbarlo,
  evtfullerrorhigh,evtfullerrorlo,eventsFull,Esparam;//frBKGparam
  Int_t ndof;
  TTree *covar = new TTree("covar","Covariance Plot");
  covar->Branch("nEvents", &nEvents,"nEvents/I");
  covar->Branch("Peak", &fpeak,"fpeak/F");
  covar->Branch("PeakErrHigh", &peakerrorhigh,"peakerrorhigh/F");
  covar->Branch("PeakErrLo", &peakerrorlo,"peakerrorlo/F");
  covar->Branch("Width", &fsigma,"fsigma/F");
  covar->Branch("WidthErrHigh", &widtherrorhigh,"widtherrorhigh/F");
  covar->Branch("WidthErrLo", &widtherrorlo,"widtherrorlo/F");
  covar->Branch("ChiSq", &chiSq,"chiSq/F");
  covar->Branch("1stPeak", &fstpeak,"fstpeak/F");
  covar->Branch("2ndPeak", &scdpeak,"scdpeak/F");
  covar->Branch("Alpha", &fcbalphaparam,"fcbalphaparam/F");
  covar->Branch("Ndeg", &fcbndegparam,"fcbndegparam/F");
  covar->Branch("comCnst", &comCnstparam,"comCnstparam/F");
  covar->Branch("combeta", &combetaparam,"combetaparam/F");
  covar->Branch("frFull", &frFullparam,"frFullparam/F");
  covar->Branch("frFrst", &frFrstparam,"frFrstparam/F");
  covar->Branch("frScnd", &frScndparam,"frScndparam/F");
  covar->Branch("frBKG", &frBKGparam,"frBKGparam/F");
  covar->Branch("crystalNo", &crystalNoparam,"crystalNoparam/F");
  covar->Branch("convgstatus", &convergencestatus,"convergencestatus/I");
  covar->Branch("pval", &pval,"pval/F");
  covar->Branch("h_means", &h_means,"h_means/F");
  covar->Branch("h_stddevs", &h_stddevs,"h_stddevs/F");
  covar->Branch("unreducedchi2", &unreducedchi2,"unreducedchi2/F");
  covar->Branch("fval", &fval,"fval/F");
  covar->Branch("m", &mparam,"mparam/F");
  covar->Branch("eta", &etaparam,"etaparam/F");
  covar->Branch("ndof", &ndof,"ndof/I");
  covar->Branch("errbarhigh", &errbarhigh,"errbarhigh/I");
  covar->Branch("errbarlo", &errbarlo,"errbarlo/I");
  covar->Branch("evtfullerrorhigh", &evtfullerrorhigh,"evtfullerrorhigh/I");
  covar->Branch("evtfullerrorlo", &evtfullerrorlo,"evtfullerrorlo/I");
  covar->Branch("eventsFull", &eventsFull,"eventsFull/I");
  covar->Branch("Esparam", &Esparam,"Esparam/I");

  auto start_bin = high_resolution_clock::now();
  for(int cryNum=anacrys_start; cryNum<anacrys_end; cryNum++){
    auto [hSum, file] = get_data_histogram(cryNum, disk); // unpack pair    
    auto [mean, stddev] = ComputeHistogramStats(hSum);
		h_means   = mean;
		h_stddevs = stddev;
    SourceFitter *fit = new SourceFitter();
    fit->FitCrystal(hSum, alg, cryNum, covar, nEvents,convergencestatus, fpeak, peakerrorhigh,peakerrorlo, 
                    fsigma,widtherrorhigh,widtherrorlo, chiSq, fstpeak, scdpeak,
                    fcbalphaparam, fcbndegparam, comCnstparam, combetaparam,
                    frFullparam, frFrstparam, frScndparam, crystalNoparam, frBKGparam,
                    pval,h_means,h_stddevs,unreducedchi2,fval,mparam,etaparam,ndof,contour,errbarhigh, errbarlo,
                    evtfullerrorhigh,evtfullerrorlo,eventsFull,Esparam );
    file->Close();    
    delete file;     
    delete fit;      
    delete hSum;
}
// --------------------------------------------------------
// OVERLAY BLOCK (Optional)
// --------------------------------------------------------
if (doOverlay) {
	for (int cryNum = anacrys_start; cryNum + 1 < anacrys_end; cryNum += 2) {

    auto [hist_even, file_even] = get_data_histogram(cryNum, disk);
    auto [hist_odd, file_odd]  = get_data_histogram(cryNum + 1, disk);	 
     hist_even->SetDirectory(0);
     hist_odd->SetDirectory(0);
    // Style histograms
    hist_even->SetLineColor(kBlue);
    hist_even->SetLineWidth(2);
    hist_odd->SetLineColor(kRed);
    hist_odd->SetLineWidth(2);
    // Create residual histogram
    TH1F* residual = (TH1F*)hist_odd->Clone(Form("residual_%d_%d", cryNum, cryNum+1));
    residual->SetDirectory(0);
    residual->Reset();
    int nBins = hist_odd->GetNbinsX();
    	for (int i = 1; i <= nBins; ++i) {
        double odd   = hist_odd->GetBinContent(i);
        double even  = hist_even->GetBinContent(i);
        double denom = sqrt(odd + even);
        double value = (denom > 0) ? (odd - even) / denom : 0;
        residual->SetBinContent(i, value);

        double err_odd  = hist_odd->GetBinError(i);
        double err_even = hist_even->GetBinError(i);
        double err = denom > 0 ? sqrt(err_odd*err_odd + err_even*err_even) / denom : 0;
        residual->SetBinError(i, err);
    }
    // Create canvas for this pair
    TCanvas* cOverlay = new TCanvas(Form("cOverlay_%d_%d", cryNum, cryNum+1),
                                    Form("Even/Odd Overlay %d & %d", cryNum, cryNum+1),
                                    800, 800);
    cOverlay->Divide(1, 2, 0, 0);
    // --- Top pad: overlay ---
    cOverlay->cd(1);
    gPad->SetPad(0.0, 0.3, 1.0, 1.0);
    hist_even->Draw("hist");
    hist_even->GetYaxis()->SetRangeUser(0, 5000);
    hist_odd->Draw("hist same");
    hist_odd->GetYaxis()->SetRangeUser(0, 5000);

    TLegend* leg = new TLegend(0.6, 0.7, 0.88, 0.88);
    leg->AddEntry(hist_even, Form("Even SiPM (%d)", cryNum), "l");
    leg->AddEntry(hist_odd, Form("Odd SiPM (%d)", cryNum+1), "l");
    leg->Draw();

    TPaveText* statsText = new TPaveText(0.15, 0.7, 0.4, 0.88, "NDC");
    statsText->SetFillColor(0);
    statsText->SetTextSize(0.03);
    statsText->AddText(Form("Even Entries: %.0f", hist_even->GetEntries()));
    statsText->AddText(Form("Odd Entries: %.0f", hist_odd->GetEntries()));
    statsText->Draw();

    // --- Bottom pad: residual ---
    cOverlay->cd(2);
    gPad->SetPad(0.0, 0.0, 1.0, 0.3);
    residual->SetTitle("Normalised Residual (Residual / sqrt(Bin Count))");
    residual->GetXaxis()->SetTitle("ADC");
    residual->GetXaxis()->SetLabelSize(0.06);
    residual->GetXaxis()->SetTitleSize(0.07);
    residual->GetYaxis()->SetTitle("Counts");
    residual->GetYaxis()->SetLabelSize(0.06);
    residual->GetYaxis()->SetTitleSize(0.07);
    residual->Draw("hist");

    gPad->Update();
    if (auto stats = (TPaveStats*)residual->FindObject("stats")) {
        stats->SetX1NDC(0.7);
        stats->SetX2NDC(0.95);
        stats->SetY1NDC(0.15);
        stats->SetY2NDC(0.45);
        stats->SetTextSize(0.07);
    }
    // --- Save to a unique ROOT file ---
	TString oName = "mu2e_simu_hist_" + std::to_string(cryNum) + "_" + std::to_string(cryNum+1) + ".root";
    TFile* outputFile = new TFile(oName, "RECREATE");
    cOverlay->Write();
    outputFile->Close();
    delete outputFile;
    delete cOverlay;
    delete residual;
		delete leg;
		delete statsText;
		delete hist_even;
		delete hist_odd;
		}
}      
  auto end_bin = high_resolution_clock::now();
// --------------------------------------------------------
// SUMMARY & PLOTTING
// --------------------------------------------------------  
  std::cout<<" ******** Av. Time take to fit crystal: "<<duration_cast<seconds>((end_bin - start_bin)/(anacrys_end-anacrys_start))<<std::endl;
  TFile *globalPlots = new TFile("globalPlots.root", "RECREATE");
  SourcePlotter *plot = new SourcePlotter();
  plot->ParamPlots(covar, table, globalPlots, anacrys_start, anacrys_end);
  table->cd();
  table -> Write();
  table -> Close();
  globalPlots -> Write();
  globalPlots -> Close();

  std::cout << "================ Refit Summary =================" << std::endl;
  auto printConv = [](std::string label, int count, const std::vector<int>& list) {
    std::cout << label << ": " << count;
    if (!list.empty()) {
        std::cout << " (";
        for (size_t i = 0; i < list.size(); ++i) {
            std::cout << list[i];
            if (i + 1 < list.size()) std::cout << ", ";
        }
        std::cout << ")";
    }
    std::cout << std::endl;
};

  std::cout << "First fit converged: " << SourceFitter::nFirstFitConverged << std::endl;  
printConv("Second fit converged", SourceFitter::nSecondFitConverged, SourceFitter::crystalsSecondFitConverged);
printConv("Third fit converged",  SourceFitter::nThirdFitConverged,  SourceFitter::crystalsThirdFitConverged);
std::cout << "\nRetry attempts per crystal:\n";
for (const auto &p : SourceFitter::thirdFitRetryCount) {
    std::cout << "  Crystal " << p.first << ": " << p.second << " attempts\n";
}

if (!SourceFitter::convFailures.empty()) {
    std::cout << "Crystals with non-zero convergence status:\n";
    for (const auto& [cryNo, status] : SourceFitter::convFailures) {
        std::cout << "Crystal " << cryNo << " : status = " << status << std::endl;
    }
}

std::cout << "================================================" << std::endl;

  std::cout<<"Finished processing ..."<<std::endl;
  return 0;
}
