#include "CaloCalibration/SourceCalib/inc/MakeAnalysisTree.hh"
#include "CaloCalibration/SourceCalib/inc/SourceFitter.hh"
#include "CaloCalibration/SourceCalib/inc/SourcePlotter.hh"
#include <chrono>
using namespace std::chrono;

using namespace CaloSourceCalib;

TString filepath = "/pnfs/mu2e/tape/usr-nts/nts/sophie/SourceCalibSimAna/ADCCuts/root/56/e6/nts.sophie.SourceCalibSimAna.ADCCuts.0.root";

/*function to extract the TTree from the SourceCalibAna output*/
TH1F* get_data_histogram(int cryNum){
    TFile *f =  new TFile(filepath);
    TString crystalNumber = to_string(cryNum);
    TH1F* hist = (TH1F*)f->Get("CaloSourceCalibDigiAna/crystals_ADC/hspec_" + crystalNumber); 
    return hist;
}

/*main function allows a loop over all crystals or a choice of a single crystal*/
int main(int argc, char* argv[]){
  std::cout<<"========== Welcome to the Mu2e Source Calibration Analysis =========="<<std::endl;
  int anacrys_start = 674; //starting crystal//675
  int anacrys_end = 1348; //final crystal//680
  TString alg = "nll"; // fitting alg (nll=NLL, chi2=chi2 fit)
  if(strcmp( argv[1], "chooseCrystal") == 0 ){
    cout<<"crystal to be analyzed (int) : "<<endl;
    cin>>anacrys_start;
    anacrys_end = anacrys_start+1;
  }
 
  TFile *table = new TFile("arXivTable.root", "RECREATE");
  Int_t nEvents;
  Float_t fpeak, dpeak, fsigma, chiSq, fstpeak, fstsigma, scdpeak,scdsigma,fcbalphaparam,fcbndegparam,Aparam,Bparam,Cparam,fullResparam,fstResparam,scdResparam,comCnstparam,
  combetaparam,frFullparam,frFrstparam,frScndparam,crystalNoparam,frBKGparam;//frBKGparam
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
  auto start_bin = high_resolution_clock::now();
  //for(int cryNum=anacrys_start; cryNum<anacrys_end; cryNum++){ //uncomment
	for(int cryNum=1330; cryNum<1348; cryNum++){//674 --> 1348 //delete before push
    TH1F* h = get_data_histogram(cryNum);
    SourceFitter *fit = new SourceFitter();
    fit->FitCrystal(h,alg, cryNum, covar, nEvents, fpeak, dpeak, fsigma,chiSq, fstpeak, fstsigma, scdpeak,scdsigma,fcbalphaparam,fcbndegparam,Aparam,Bparam,Cparam,
    fullResparam,fstResparam,scdResparam,comCnstparam,combetaparam,
    frFullparam,frFrstparam,frScndparam,crystalNoparam,frBKGparam);
  };

  auto end_bin = high_resolution_clock::now();
  std::cout<<" ******** Av. Time take to fit crystal: "<<duration_cast<seconds>((end_bin - start_bin)/(anacrys_end-anacrys_start))<<std::endl;

  
  TFile *globalPlots = new TFile("globalPlots.root", "RECREATE");
  SourcePlotter *plot = new SourcePlotter();
  plot->ParamPlots(covar, table, globalPlots);
  table -> Write();
  table -> Close();
  globalPlots -> Write();
  globalPlots -> Close();
  std::cout<<"Finished processing ..."<<std::endl;
  return 0;
}
