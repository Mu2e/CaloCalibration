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
  int anacrys_end = 680; //final crystal//680
  std::cout<<"choosing crystal "<<std::endl;
  if(strcmp( argv[1], "chooseCrystal") == 0 ){
    cout<<"crystal to be analyzed (int) : "<<endl;
    cin>>anacrys_start;
    anacrys_end = anacrys_start+1;
  }
 
  TFile *ouptFile = new TFile("paraFile.root", "RECREATE");
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

  for(int cryNum=anacrys_start; cryNum<anacrys_end; cryNum++){
    TH1F* h = get_data_histogram(cryNum);
    std::cout<<"Running fitter ....."<<std::endl;
    auto start_bin = high_resolution_clock::now();

    SourceFitter *fit = new SourceFitter();
    fit->FitCrystal(h,"nll", cryNum, covar, nEvents, fpeak, dpeak, fsigma,chiSq, fstpeak, fstsigma, scdpeak,scdsigma,fcbalphaparam,fcbndegparam,Aparam,Bparam,Cparam,
    fullResparam,fstResparam,scdResparam,comCnstparam,combetaparam,
    frFullparam,frFrstparam,frScndparam,crystalNoparam,frBKGparam);

    auto end_bin = high_resolution_clock::now();
    std::cout<<" ******** Time take to fit crystal: "<<cryNum<<" "<<duration_cast<seconds>(end_bin - start_bin)<<std::endl;
  };

  
  SourcePlotter *plot = new SourcePlotter();
  plot->ParamPlots(covar,ouptFile);
  ouptFile -> Write();
  ouptFile -> Close();
  std::cout<<"Finished processing ..."<<std::endl;
  return 0;
}
