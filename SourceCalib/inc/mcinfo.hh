#ifndef _mcinfo_hh
#define _mcinfo_hh

#include <fstream>
#include <iostream>
#include <chrono>
#include <numeric>
#include <iomanip>
#include <algorithm> 
#include "RooAbsReal.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TGraph.h"
#include "TLine.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TAxis.h"
#include <Riostream.h>
#include <vector>
#include "THStack.h"
#include "TStyle.h"


#include "TPaveText.h"
#include "TLegend.h"
#include "TPaveLabel.h"
#include "TAttFill.h"
using namespace std;
using namespace TMath;
using namespace RooFit;

namespace CaloSourceCalib{
  class mcinfo  {
      public:
        mcinfo() = default;
        mcinfo(const mcinfo &) = default;
        mcinfo& operator=(const mcinfo &) = default;
        virtual ~mcinfo() = default;
        #ifndef __CINT__
        void RunMCTruth(TH1F* hist, int cryNum, int disk,TTree *trueinfo,Int_t &cryNumparam,  Int_t &tot_evts,Int_t &mainpeak,Int_t &first_espeak,Int_t &second_espeak,Int_t &background,Float_t &frmainpeak,Float_t &frfirst_espeak,Float_t &frsecond_espeak,Float_t &frbackground,Int_t &compton1,Int_t &compton2,Int_t &compton3,Float_t &frcompton1,Float_t &frcompton2,Float_t &frcompton3);
        void FinalizeMCSummary(TTree* trueinfo, TTree* truthfit,bool useCompton);
        // MC Truth fitting (MeV space) for comparison with reconstructed fits (ADC space)
        void FitMCTruthSpectrum(TH1F* hist_truth_mev, int crystalNo, TTree *truthfit_tree,
                                Float_t &truth_peak_mev, Float_t &truth_fr_full,
                                Float_t &truth_fr_1st, Float_t &truth_fr_2nd,
                                Float_t &truth_fr_c1, Float_t &truth_fr_c2, Float_t &truth_fr_c3,
                                Float_t &truth_fr_ebk,
                                Float_t &truth_width, Float_t &truth_alpha, Float_t &truth_chi2, Int_t &truth_ndof,
                                bool useCompton = true);
        #endif
        ClassDef (mcinfo,1);
    };
}
#endif /* mcinfo.hh */
