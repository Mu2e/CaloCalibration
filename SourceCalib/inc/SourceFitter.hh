#ifndef _SourceFitter_hh
#define _SourceFitter_hh

#include <fstream>
#include <iostream>
#include "TSystem.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1.h"

#include "TF1.h"
#include "TH2.h"
#include "TMath.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TMath.h"
#include "TMarker.h"
#include "TLine.h"
#include "TLegend.h"
#include <Riostream.h>

#include "TPaveStats.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPaveLabel.h"
#include "TAttFill.h"
#include "TObject.h"
// add roofit header files
#include "RooHist.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooKeysPdf.h"
#include "RooHistPdf.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooCBShape.h"
#include "RooGenericPdf.h"
#include "RooFFTConvPdf.h"
#include "RooFitResult.h"
#include "RooUniform.h"
#include "RooMinimizer.h"
#include "RooAbsReal.h"
#include "RooFormulaVar.h"
#include "RooAddition.h"
#include "RooChi2Var.h" 
#include "RooMsgService.h"
#include "Math/MinimizerOptions.h"
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include <random>
using namespace std;
using namespace TMath;
using namespace RooFit;

namespace CaloSourceCalib{
  class SourceFitter : public TObject {
      public:
        //explicit SourceFitter(){};
        SourceFitter() = default;
        //explicit SourceFitter(const SourceFitter &){};
        //SourceFitter& operator = (const SourceFitter &);
        virtual ~SourceFitter() = default;
        //#ifndef __CINT__
        void FitCrystal(TH1F* histogram, TString opt, int cryNum, TTree* covar, Int_t &nEvents,Int_t &convergencestatus, Float_t &fpeak, Float_t &peakerrorhigh,Float_t &peakerrorlo, Float_t &fsigma,Float_t &widtherrorhigh,Float_t &widtherrorlo, Float_t &chiSq, Float_t &fstpeak, Float_t &scdpeak,Float_t &fcbalphaparam,Float_t &fcbndegparam,Float_t &comCnstparam, Float_t &combetaparam, Float_t &frFullparam, Float_t &frFrstparam,Float_t &frScndparam,Float_t &crystalNoparam,Float_t &frBKGparam, Float_t &pval,Float_t &h_means,Float_t &h_stddevs,Float_t &unreducedchi2,Float_t &fval,Float_t &mparam,Float_t &etaparam,Int_t &ndof,bool contour,Float_t &errbarhigh, Float_t &errbarlo,Float_t &evtfullerrorhigh,Float_t &evtfullerrorlo,Float_t &eventsFull,Float_t &Esparam);
        void MCFitCrystal(int crystalNo, TString opt);
        void randomize_all_parameters();
        //counter for crystals that need a refit
        static int nSecondFits;
        static int nThirdFits;
        // vector of crystal numbers that need refit or fail
        static std::vector<int> crystalsSecondFit;
        static std::vector<int> crystalsThirdFit;
        static std::vector<std::pair<int,int>> convFailures; // {crystalNo, convStatus}
        //counters for how many crystals converged at which stage
        static int nFirstFitConverged;
        static int nSecondFitConverged;
        static int nThirdFitConverged;
        //vectors for crystal ID's that converged at different stages
        static std::vector<int> crystalsSecondFitConverged;
        static std::vector<int> crystalsThirdFitConverged;
        //map to show how many retries attempted, with crystal ID
        static std::map<int,int> thirdFitRetryCount; 

        //#endif

        ClassDef (SourceFitter,1);
    };
}
#endif  // _SourceFitter_hh
