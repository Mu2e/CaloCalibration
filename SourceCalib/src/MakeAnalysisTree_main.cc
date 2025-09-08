#include "CaloCalibration/SourceCalib/inc/MakeAnalysisTree.hh"
#include "CaloCalibration/SourceCalib/inc/SourceFitter.hh"
#include "CaloCalibration/SourceCalib/inc/SourcePlotter.hh"
#include <chrono>
using namespace std::chrono;

using namespace CaloSourceCalib;


TString filepath_disk0 = "/exp/mu2e/app/home/mu2epro/sourcecalib/disk0/SourceAna1e9.root";

TString filepath_disk1 = "/exp/mu2e/app/home/mu2epro/sourcecalib/disk1/nts.mu2e.SourceCalibAna1e9.0.root";

/*function to extract the TTree from the SourceCalibAna output*/
/*std::pair<TH1F*, TFile*> get_data_histogram(int cryNum, int disk) {
    TString filepath;
    TString histPath;
    if (disk == 0) filepath = filepath_disk0;
    if (disk == 1) filepath = filepath_disk1;
    TFile *f = new TFile(filepath);
    TString crystalNumber = to_string(cryNum);
    TH1F* hist = (TH1F*)f->Get("SourceAna/sipm_ADC/sipm_" + crystalNumber); 
    return std::make_pair(hist, f);*/
    
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
  //int nCry = anacrys_end - anacrys_start;  // number of crystals to analyze
  // Optional overlay flag (default = false)
  bool doOverlay = false;
    if (argc >= 6 && TString(argv[5]) == "overlay") {
        doOverlay = true;
    }

  TFile *table = new TFile("arXivTable.root", "RECREATE");
  Int_t nEvents;
  Float_t fpeak, dpeak, fsigma, chiSq, fstpeak, fstsigma, scdpeak,scdsigma,fcbalphaparam,fcbndegparam,Aparam,Bparam,Cparam,fullResparam,fstResparam,scdResparam,comCnstparam,
  combetaparam,frFullparam,frFrstparam,frScndparam,crystalNoparam,frBKGparam, convergencestatus,errbar,pval,h_means,h_stddevs;//frBKGparam
  TTree *covar = new TTree("covar","Covariance Plot");
	//TH1F* h_means = new TH1F("h_means", "Mean of Raw Histograms;Crystal;Mean ADC", nCry, anacrys_start, anacrys_end);
	//TH1F* h_stddevs = new TH1F("h_stddevs", "Width of Raw Histograms;Crystal;Std Dev (ADC)", nCry, anacrys_start, anacrys_end);
  covar->Branch("nEvents", &nEvents,"nEvents/I");
  covar->Branch("Peak", &fpeak,"fpeak/F");
  covar->Branch("PeakErr", &dpeak,"dpeak/F");
  covar->Branch("Width", &fsigma,"fsigma/F");
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
  covar->Branch("convgstatus", &convergencestatus,"convergencestatus/F");
  covar->Branch("errorbar", &errbar,"errbar/F");
  covar->Branch("pval", &pval,"pval/F");
  covar->Branch("h_means", &h_means,"h_means/F");
  covar->Branch("h_stddevs", &h_stddevs,"h_stddevs/F");
  auto start_bin = high_resolution_clock::now();
  /*for(int cryNum=anacrys_start; cryNum<anacrys_end; cryNum++){
    TH1F* h = get_data_histogram(cryNum, disk);
    SourceFitter *fit = new SourceFitter();
    fit->FitCrystal(h,alg, cryNum, covar, nEvents, fpeak, dpeak, fsigma,chiSq, fstpeak, fstsigma, scdpeak,scdsigma,fcbalphaparam,fcbndegparam,Aparam,Bparam,Cparam,
    fullResparam,fstResparam,scdResparam,comCnstparam,combetaparam,
    frFullparam,frFrstparam,frScndparam,crystalNoparam,frBKGparam);
  };*/
  for(int cryNum=anacrys_start; cryNum<anacrys_end; cryNum++){
    auto [h, file] = get_data_histogram(cryNum, disk); // unpack pair
    
    auto [mean, stddev] = ComputeHistogramStats(h);
		//h_means->SetBinContent(cryNum - anacrys_start + 1, mean);
		//h_stddevs->SetBinContent(cryNum - anacrys_start + 1, stddev);
		//std::cout << "cryNum: " << cryNum << ", bin = " << cryNum - anacrys_start + 1 << std::endl;
		h_means   = mean;
		h_stddevs = stddev;


    SourceFitter *fit = new SourceFitter();
    fit->FitCrystal(h, alg, cryNum, covar, nEvents, fpeak, dpeak, fsigma, chiSq, fstpeak, fstsigma, scdpeak, scdsigma,
                    fcbalphaparam, fcbndegparam, Aparam, Bparam, Cparam,
                    fullResparam, fstResparam, scdResparam, comCnstparam, combetaparam,
                    frFullparam, frFrstparam, frScndparam, crystalNoparam, frBKGparam, convergencestatus,errbar,pval,h_means,h_stddevs);

    file->Close();    // ? closes input file for current crystal
    delete file;      // ? avoids memory leak
    delete fit;       // (optional cleanup)
    delete h;
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

    // Cleanup
    file_even->Close();
    file_odd->Close();
    delete file_even;
    delete file_odd;
    delete cOverlay;
    delete residual;
		delete leg;
		delete statsText;
		delete hist_even;
		delete hist_odd;
		}
}      
    //file->Close();    // ? closes input file for current crystal
    //delete file;      // ? avoids memory leak
    //delete fit;       // (optional cleanup)
    //delete h;
//};


  auto end_bin = high_resolution_clock::now();
  std::cout<<" ******** Av. Time take to fit crystal: "<<duration_cast<seconds>((end_bin - start_bin)/(anacrys_end-anacrys_start))<<std::endl;
  
  TFile *globalPlots = new TFile("globalPlots.root", "RECREATE");
  SourcePlotter *plot = new SourcePlotter();
  plot->ParamPlots(covar, table, globalPlots, anacrys_start, anacrys_end);//add arguement of crynum
  table->cd();
  table -> Write();
  table -> Close();
  globalPlots -> Write();
  globalPlots -> Close();
  std::cout<<"Finished processing ..."<<std::endl;
  return 0;
}
