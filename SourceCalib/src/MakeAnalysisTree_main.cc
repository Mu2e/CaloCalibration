#include "CaloCalibration/SourceCalib/inc/MakeAnalysisTree.hh"
#include "CaloCalibration/SourceCalib/inc/SourceFitter.hh"
#include "CaloCalibration/SourceCalib/inc/SourcePlotter.hh"
#include <chrono>
using namespace std::chrono;

using namespace CaloSourceCalib;

//TString filepath_disk0 = "/pnfs/mu2e/tape/usr-nts/nts/sophie/SourceCalibSimAna/Disk0/root/9c/e5//nts.sophie.SourceCalibSimAna.Disk0.2.root";
TString filepath_disk0 = "/exp/mu2e/app/users/sophie/CaloCalib/nts.sophie.SourceCalibSimAna.v2.0.root";

TString filepath_disk1 = "/pnfs/mu2e/tape/usr-nts/nts/sophie/SourceCalibSimAna/Disk1/root/fc/49/nts.sophie.SourceCalibSimAna.Disk1.2.root";


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

    return std::make_pair(hist, f);
}

/*main function allows a loop over all crystals or a choice of a single crystal*/
int main(int argc, char* argv[]){
  std::cout<<"========== Welcome to the Mu2e Source Calibration Analysis =========="<<std::endl;
  int anacrys_start = std::atoi(argv[1]); //starting crystal//675
  int anacrys_end = std::atoi(argv[2]); //final crystal//680
  TString alg = argv[3]; // fitting alg (nll=NLL, chi2=chi2 fit)
  int disk = std::atoi(argv[4]); //disk number 0 or 1
  TFile *table = new TFile("arXivTable.root", "RECREATE");
  Int_t nEvents;
  Float_t fpeak, dpeak, fsigma, chiSq, fstpeak, fstsigma, scdpeak,scdsigma,fcbalphaparam,fcbndegparam,Aparam,Bparam,Cparam,fullResparam,fstResparam,scdResparam,comCnstparam,
  combetaparam,frFullparam,frFrstparam,frScndparam,crystalNoparam,frBKGparam, convergencestatus,errbar;//frBKGparam
  TTree *covar = new TTree("covar","Covariance Plot");
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
    SourceFitter *fit = new SourceFitter();
    fit->FitCrystal(h, alg, cryNum, covar, nEvents, fpeak, dpeak, fsigma, chiSq, fstpeak, fstsigma, scdpeak, scdsigma,
                    fcbalphaparam, fcbndegparam, Aparam, Bparam, Cparam,
                    fullResparam, fstResparam, scdResparam, comCnstparam, combetaparam,
                    frFullparam, frFrstparam, frScndparam, crystalNoparam, frBKGparam, 											convergencestatus,errbar);

// Example: Compare crystal 1042 (even) and 1043 (odd) from the same disk
auto [hist_even, file_even] = get_data_histogram(1042, disk);
auto [hist_odd, file_odd] = get_data_histogram(1043, disk);

// Prepare histograms for overlay
hist_even->SetLineColor(kBlue);
hist_even->SetLineWidth(2);
hist_odd->SetLineColor(kRed);
hist_odd->SetLineWidth(2);

// Create residual histogram: (odd - even)
TH1F* residual = (TH1F*)hist_odd->Clone("residual");
residual->SetDirectory(0);
residual->Sumw2();
residual->Add(hist_even, -1);  // residual = hist_odd - hist_even
// Don't set residual->SetEntries(...) ? let ROOT handle it

// Setup canvas
TCanvas* cOverlay = new TCanvas("cOverlay", "Even/Odd SiPM Overlay with Residual", 800, 800);
cOverlay->Divide(1, 2, 0, 0);

// ???????????????? Top pad: Overlay ????????????????
cOverlay->cd(1);
gPad->SetPad(0.0, 0.3, 1.0, 1.0);
hist_even->Draw("hist");
hist_odd->Draw("hist same");

TLegend* leg = new TLegend(0.6, 0.7, 0.88, 0.88);
leg->AddEntry(hist_even, "Even SiPM (1042)", "l");
leg->AddEntry(hist_odd, "Odd SiPM (1043)", "l");
leg->Draw();

// Add TPaveText to show entry counts
TPaveText* statsText = new TPaveText(0.15, 0.7, 0.4, 0.88, "NDC");
statsText->SetFillColor(0);
statsText->SetTextSize(0.03);
statsText->AddText(Form("Even Entries: %.0f", hist_even->GetEntries()));
statsText->AddText(Form("Odd Entries: %.0f", hist_odd->GetEntries()));
statsText->Draw();

// ???????????????? Bottom pad: Residual ????????????????
cOverlay->cd(2);
gPad->SetPad(0.0, 0.0, 1.0, 0.3);

// Show stat box and style it
gStyle->SetOptStat(1110);
residual->SetTitle("Residual (Odd - Even)");
residual->GetXaxis()->SetTitle("ADC");
residual->GetYaxis()->SetTitle("Counts");
residual->Draw("hist");

gPad->Update();  // Generate stats box
TPaveStats* stats = (TPaveStats*)residual->FindObject("stats");
if (stats) {
    stats->SetX1NDC(0.7);
    stats->SetX2NDC(0.95);
    stats->SetY1NDC(0.15);
    stats->SetY2NDC(0.45);
    stats->SetTextSize(0.04);
}

// Save to ROOT file
TFile* outputFile = new TFile("EvenOddOverlay.root", "RECREATE");
cOverlay->Write("cOverlay");
outputFile->Close();
delete outputFile;

// Cleanup
file_even->Close();
file_odd->Close();
delete file_even;
delete file_odd;
delete cOverlay;

                    
    file->Close();    // ? closes input file for current crystal
    delete file;      // ? avoids memory leak
    delete fit;       // (optional cleanup)
};


  auto end_bin = high_resolution_clock::now();
  std::cout<<" ******** Av. Time take to fit crystal: "<<duration_cast<seconds>((end_bin - start_bin)/(anacrys_end-anacrys_start))<<std::endl;
  
  TFile *globalPlots = new TFile("globalPlots.root", "RECREATE");
  SourcePlotter *plot = new SourcePlotter();
  plot->ParamPlots(covar, table, globalPlots, anacrys_start, anacrys_end);//add arguement of crynum
  table -> Write();
  table -> Close();
  globalPlots -> Write();
  globalPlots -> Close();
  std::cout<<"Finished processing ..."<<std::endl;
  return 0;
}
