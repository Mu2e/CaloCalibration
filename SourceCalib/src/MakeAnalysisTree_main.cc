#include "CaloCalibration/SourceCalib/inc/MakeAnalysisTree.hh"
#include "CaloCalibration/SourceCalib/inc/SourceFitter.hh"
#include "CaloCalibration/SourceCalib/inc/SourcePlotter.hh"
#include <chrono>
using namespace std::chrono;

using namespace CaloSourceCalib;
double GetTreeMean(TTree* t, const char* branchName) {
    if (!t->GetBranch(branchName)) {
        std::cerr << "ERROR: Branch " << branchName << " not found in tree!\n";
        return -1;
    }
    t->Draw(branchName, "", "goff");
    Long64_t n = t->GetSelectedRows();
    if (n <= 0) return -1;
    return TMath::Mean(n, t->GetV1());
}
//method with one big root file
TString filepath_disk0 = "/exp/mu2e/app/home/mu2epro/sourcecalib/disk0/SourceAna1e9.root";



TString filepath_disk1 = "/exp/mu2e/app/home/mu2epro/sourcecalib/disk1/nts.mu2e.SourceCalibAna1e9.0.root";
//TString filepath_disk1 = "/exp/mu2e/app/users/hjafree/SourceFitDir/combined_disk1.root";

/*function to extract the TTree from the SourceCalibAna output*/
//original method to read sipms without cry option
/*std::pair<TH1F*, TFile*> get_data_histogram(int cryNum, int disk) {
    TString filepath;
    TString histPath;
    if (disk == 0) filepath = filepath_disk0;
    if (disk == 1) filepath = filepath_disk1;
    TFile *f = new TFile(filepath);
    TString crystalNumber = to_string(cryNum);
    TH1F* hist = (TH1F*)f->Get("SourceAna/sipm_ADC/sipm_" + crystalNumber); 
    return std::make_pair(hist, f);*/
    
//method to loop over several files
/*std::vector<std::string> ImportFileList(const char* filelist) {
    std::ifstream files(filelist);
    std::vector<std::string> fileNames;
    std::string line;
    while (std::getline(files, line)) {
        if (!line.empty()) fileNames.push_back(line);
    }
    return fileNames;
}
struct DiskConfig {
    std::vector<std::string> files;
    TString subdir; // "sipm_ADC" or "crystals_ADC"
    TString histPrefix; // e.g. "sipm_" or "cry_"
};
DiskConfig SetupFilesAndFolder(int disk) {
    DiskConfig cfg;
    // choose filelist file name and histogram folder/prefix depending on disk
    if (disk == 0 || disk == 2) {
        cfg.files = ImportFileList("disk0_files.txt");
        if (disk == 0) { cfg.subdir = "sipm_ADC"; cfg.histPrefix = "sipm_"; }
        else           { cfg.subdir = "crystals_ADC"; cfg.histPrefix = "cry_"; }
    } else if (disk == 1 || disk == 3) {
        cfg.files = ImportFileList("disk1_files.txt");
        if (disk == 1) { cfg.subdir = "sipm_ADC"; cfg.histPrefix = "sipm_"; }
        else           { cfg.subdir = "crystals_ADC"; cfg.histPrefix = "cry_"; }
    } else {
        throw std::runtime_error("Invalid disk number");
    }
    if (cfg.files.empty()) {
        std::cerr << "[SetupFilesAndFolder] WARNING: file list is empty for disk " << disk << "\n";
    }
    return cfg;
}   */
std::pair<TH1F*, TFile*> get_data_histogram(int cryNum, int disk) {
    TString filepath;
    TString histPath;
    if (disk == 0 || disk == 2) {
        filepath = filepath_disk0;
        histPath = (disk == 0) ? "SourceAna/sipm_ADC/sipm_" : "SourceAna/crystals_ADC/cry_";
    } else if (disk == 1 || disk == 3) {
        filepath = filepath_disk1;
        histPath = (disk == 1) ? "SourceAna/sipm_ADC/sipm_" : "SourceAna/crystals_ADC/cry_";
    }
 

// Example: collect every histogram in SourceAna/sipm_ADC
/*TH1F* SumHistogramAcrossFiles(const std::vector<std::string>& files, const TString& fullPath, const char* cloneName = nullptr) {
    TH1F* hSum = nullptr;
    for (const auto& fname : files) {
        TFile f(fname.c_str());
        if (f.IsZombie()) {
            std::cerr << "[SumHistogramAcrossFiles] cannot open: " << fname << "\n";
            continue;
        }
        TH1F* hTmp = (TH1F*)f.Get(fullPath);
        if (!hTmp) continue;
        if (!hSum) {
            // clone first found histogram (gives an independent object)
            if (cloneName)
                hSum = (TH1F*)hTmp->Clone(cloneName);
            else
                hSum = (TH1F*)hTmp->Clone(Form("sum_%s", fullPath.Data()));
            hSum->SetDirectory(0); // detach from file
        } else {
            hSum->Add(hTmp);
        }
    }
    return hSum;
}*/
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
int main(int argc, char* argv[]){
  std::cout<<"========== Welcome to the Mu2e Source Calibration Analysis =========="<<std::endl;
  int anacrys_start = std::atoi(argv[1]); //starting crystal//675
  int anacrys_end = std::atoi(argv[2]); //final crystal//680
  TString alg = argv[3]; // fitting alg (nll=NLL, chi2=chi2 fit)
  int disk = std::atoi(argv[4]); //disk number 0 or 1
  bool doOverlay = false;
    if (argc >= 6 && TString(argv[5]) == "overlay") {
        doOverlay = true;
    }
    // Read list of ROOT files (one per line)
 // Setup file list and which histogram folder/prefix to use
    //DiskConfig cfg = SetupFilesAndFolder(disk);
    //auto& files = cfg.files; // alias

  TFile *table = new TFile("arXivTable.root", "RECREATE");
  Int_t nEvents;
  Int_t convergencestatus;
  Float_t fpeak, dpeak, fsigma, chiSq, fstpeak, fstsigma, scdpeak,scdsigma,fcbalphaparam,fcbndegparam,Aparam,Bparam,Cparam,fullResparam,fstResparam,scdResparam,comCnstparam,
  combetaparam,frFullparam,frFrstparam,frScndparam,crystalNoparam,frBKGparam
  ,errbar,pval,h_means,h_stddevs,unreducedchi2,fval,mparam,etaparam,dsigma;//frBKGparam
  Int_t ndof;
  TTree *covar = new TTree("covar","Covariance Plot");
  covar->Branch("nEvents", &nEvents,"nEvents/I");
  covar->Branch("Peak", &fpeak,"fpeak/F");
  covar->Branch("PeakErr", &dpeak,"dpeak/F");
  covar->Branch("Width", &fsigma,"fsigma/F");
  covar->Branch("WidthErr", &dsigma,"dsigma/F");
  covar->Branch("ChiSq", &chiSq,"chiSq/F");
  covar->Branch("1stPeak", &fstpeak,"fstpeak/F");
  covar->Branch("1stWidth", &fstsigma,"fstsigma/F");
  covar->Branch("2ndPeak", &scdpeak,"scdpeak/F");
  covar->Branch("2ndWidth", &scdsigma,"scdsigma/F");
  covar->Branch("Alpha", &fcbalphaparam,"fcbalphaparam/F");
  covar->Branch("Ndeg", &fcbndegparam,"fcbndegparam/F");
  covar->Branch("A", &Aparam,"Aparam/F");
  covar->Branch("B", &Bparam,"Bparam/F");
  covar->Branch("C", &Cparam,"Cparam/F");
  covar->Branch("fullRes", &fullResparam,"fullResparam/F");
  covar->Branch("fstRes", &fstResparam,"fstResparam/F");
  covar->Branch("scdRes", &scdResparam,"scdResparam/F");
  covar->Branch("comCnst", &comCnstparam,"comCnstparam/F");
  covar->Branch("combeta", &combetaparam,"combetaparam/F");
  covar->Branch("frFull", &frFullparam,"frFullparam/F");
  covar->Branch("frFrst", &frFrstparam,"frFrstparam/F");
  covar->Branch("frScnd", &frScndparam,"frScndparam/F");
  covar->Branch("frBKG", &frBKGparam,"frBKGparam/F");
  covar->Branch("crystalNo", &crystalNoparam,"crystalNoparam/F");
  covar->Branch("convgstatus", &convergencestatus,"convergencestatus/I");
  covar->Branch("errorbar", &errbar,"errbar/F");
  covar->Branch("pval", &pval,"pval/F");
  covar->Branch("h_means", &h_means,"h_means/F");
  covar->Branch("h_stddevs", &h_stddevs,"h_stddevs/F");
  covar->Branch("unreducedchi2", &unreducedchi2,"unreducedchi2/F");
  covar->Branch("fval", &fval,"fval/F");
  covar->Branch("m", &mparam,"mparam/F");
  covar->Branch("eta", &etaparam,"etaparam/F");
  covar->Branch("ndof", &ndof,"ndof/I");

  auto start_bin = high_resolution_clock::now();
  /*for(int cryNum=anacrys_start; cryNum<anacrys_end; cryNum++){
    TH1F* h = get_data_histogram(cryNum, disk);
    SourceFitter *fit = new SourceFitter();
    fit->FitCrystal(h,alg, cryNum, covar, nEvents, fpeak, dpeak, fsigma,chiSq, fstpeak, fstsigma, scdpeak,scdsigma,fcbalphaparam,fcbndegparam,Aparam,Bparam,Cparam,
    fullResparam,fstResparam,scdResparam,comCnstparam,combetaparam,
    frFullparam,frFrstparam,frScndparam,crystalNoparam,frBKGparam);
  };*/


  for(int cryNum=anacrys_start; cryNum<anacrys_end; cryNum++){
    auto [hSum, file] = get_data_histogram(cryNum, disk); // unpack pair    
    auto [mean, stddev] = ComputeHistogramStats(hSum);
		h_means   = mean;
		h_stddevs = stddev;
    SourceFitter *fit = new SourceFitter();
    fit->FitCrystal(hSum, alg, cryNum, covar, nEvents,convergencestatus, fpeak, dpeak, 
                    fsigma,dsigma, chiSq, fstpeak, fstsigma, scdpeak, scdsigma,
                    fcbalphaparam, fcbndegparam, Aparam, Bparam, Cparam,
                    fullResparam, fstResparam, scdResparam, comCnstparam, combetaparam,
                    frFullparam, frFrstparam, frScndparam, crystalNoparam, frBKGparam
                    ,errbar,pval,h_means,h_stddevs,unreducedchi2,fval,mparam,etaparam,ndof);//migrad_status,hesse_status,minos_status
   bool need_refit = (convergencestatus > 0);
   double nEventsMean   = GetTreeMean(covar, "nEvents");
   double newPeakGuess   = GetTreeMean(covar, "Peak");
   double newSigmaGuess  = GetTreeMean(covar, "Width");
   double newAlphaGuess  = GetTreeMean(covar, "Alpha");
   double newCombetaGuess  = GetTreeMean(covar, "combeta");
   double newevtsFullGuess  = (GetTreeMean(covar, "frFull"))*nEventsMean;
   double newevtsFstGuess  = (GetTreeMean(covar, "frFrst"))*nEventsMean;
   double newevtsScndGuess  = (GetTreeMean(covar, "frScnd"))*nEventsMean;
   double newevtsBkgGuess  = (GetTreeMean(covar, "frBKG"))*nEventsMean;
   std::cout << "newPeakGuess = " << newPeakGuess << std::endl;
   std::cout << "newsigmaguess = " << newSigmaGuess << std::endl;
   std::cout << "newAlphaGuess = " << newAlphaGuess << std::endl;
   std::cout << "newCombetaGuess = " << newCombetaGuess << std::endl;
   std::cout << "newevtsFullGuess = " << newevtsFullGuess << std::endl;
   std::cout << "newevtsFstGuess = " << newevtsFstGuess << std::endl;
   std::cout << "newevtsScndGuess = " << newevtsScndGuess << std::endl;
   std::cout << "newevtsBkgGuess = " << newevtsBkgGuess << std::endl;
   std::cout << "nEventsMean = " << nEventsMean << std::endl;
   
   if (need_refit) {             
   fit->SetInitialGuesses(newPeakGuess, newSigmaGuess,newAlphaGuess,newCombetaGuess,newevtsFullGuess,newevtsFstGuess,newevtsScndGuess,newevtsBkgGuess);
   fit->FitCrystal(hSum, alg, cryNum, covar, nEvents, convergencestatus,
                    fpeak, dpeak, fsigma, dsigma, chiSq,
                    fstpeak, fstsigma, scdpeak, scdsigma,
                    fcbalphaparam, fcbndegparam, Aparam, Bparam, Cparam,
                    fullResparam, fstResparam, scdResparam,
                    comCnstparam, combetaparam,
                    frFullparam, frFrstparam, frScndparam,
                    crystalNoparam, frBKGparam,
                    errbar, pval, h_means, h_stddevs,
                    unreducedchi2, fval, mparam, etaparam,ndof);
}

    file->Close();    
    delete file;     
    delete fit;      
    delete hSum;
}
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
  std::cout<<" ******** Av. Time take to fit crystal: "<<duration_cast<seconds>((end_bin - start_bin)/(anacrys_end-anacrys_start))<<std::endl;
  
  TFile *globalPlots = new TFile("globalPlots.root", "RECREATE");
  SourcePlotter *plot = new SourcePlotter();
  plot->ParamPlots(covar, table, globalPlots, anacrys_start, anacrys_end);
  table->cd();
  table -> Write();
  table -> Close();
  globalPlots -> Write();
  globalPlots -> Close();
  std::cout<<"Finished processing ..."<<std::endl;
  return 0;
}
