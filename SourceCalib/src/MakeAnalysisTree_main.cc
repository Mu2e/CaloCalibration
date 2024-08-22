#include "CaloCalibration/SourceCalib/inc/MakeAnalysisTree.hh"
#include "CaloCalibration/SourceCalib/inc/SourceFitter.hh"
#include <chrono>
using namespace std::chrono;

using namespace CaloSourceCalib;

TString filepath = "/pnfs/mu2e/tape/usr-nts/nts/sophie/SourceCalibSimAna/ADCCuts/root/56/e6/nts.sophie.SourceCalibSimAna.ADCCuts.0.root";
//"/pnfs/mu2e/tape/usr-nts/nts/hjafree/SourceCalibAna/MDC2020ae/root/0a/1f/nts.hjafree.SourceCalibAna.MDC2020ae.0.root"
//"/exp/mu2e/app/users/hjafree/SourceFitDir/10M.root"

/*function to extract the TTree from the SourceCalibAna output*/
TH1F* get_data_histogram(int cryNum){
    TFile *f =  new TFile(filepath);
    TString crystalNumber = to_string(cryNum);
    TH1F* hist = (TH1F*)f->Get("CaloSourceCalibDigiAna/crystals_ADC/hspec_" + crystalNumber); 
    return hist;
}


/*function calls RooFit and fits a given crystal file*/
//void RunRooFit(TH1F* hist) {
//  SourceFitter *fit = new SourceFitter();
//  fit->FitCrystal(hist,"nll");
//}
//void RunRooFit(Tfile *f) {
 //SourceFitter *fit = new SourceFitter();
 // fit->FitCrystal(h_spec,"chi2");
//}

/*main function allows a loop over all crystals or a choice of a single crystal*/
//read in root file in main function and pass TFile
int main(int argc, char* argv[]){
  std::cout<<"========== Welcome to the Mu2e Source Calibration Analysis =========="<<std::endl;
  int anacrys_start = 675; //starting crystal
  int anacrys_end = 680; //final crystal
  std::cout<<"choosing crystal "<<std::endl;
  //RunRooFit(*t);
  if(strcmp( argv[1], "chooseCrystal") == 0 ){
    cout<<"crystal to be analyzed (int) : "<<endl;
    cin>>anacrys_start;
    anacrys_end = anacrys_start+1;

  }
  TFile *ouptFile = new TFile("paraFile.root", "RECREATE");
  Float_t fpeak;
  TTree *covar = new TTree("covar","Covariance Plot");
  covar->Branch("Peak", &fpeak,"fpeak/F");
  //add mean sigma and chisq
  //get python document and 
  
 for(int cryNum=674; cryNum<680; cryNum++){
  TH1F* h = get_data_histogram(cryNum);
  std::cout<<"start value"<<anacrys_start<<std::endl;
  std::cout<<"end value"<<anacrys_end<<std::endl;
  std::cout<<"Running crystal binning ....."<<std::endl;
  //MakeCrystalListOutputs(anacrys_start,anacrys_end);
  //MakeCrystalBinsOutputs(anacrys_start,anacrys_end);
  SourceFitter *fit = new SourceFitter();
  fit->FitCrystal(h,"nll", cryNum, covar, fpeak);//add another arg as ttree obj
 };
 	ouptFile -> Write();
  ouptFile -> Close();

  std::cout<<"Finished processing ..."<<std::endl;
  return 0;
}
