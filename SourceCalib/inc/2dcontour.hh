#ifndef _2dcontour_hh
#define _2dcontour_hh

// ROOT and RooFit includes
#include "RooAbsPdf.h"
#include "RooDataHist.h"
#include "RooRealVar.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TMarker.h"
#include "TLegend.h"
#include "RooChi2Var.h"
#include "RooMinimizer.h"
#include "TStyle.h"
#include "TString.h"
#include <iostream>
namespace CaloSourceCalib {

    void MakeContourPlot(
        RooAbsPdf& fitFun, 
        RooDataHist& chSpec, 
        TString fitType, 
        int cryNum,
        // --- Selection Strings ---
        TString xName, 
        TString yName,
        // --- Special Variables (Manual Values passed from Fitter) ---
        RooRealVar& varPeak,  double valPeak,  double errLoPeak,  double errHiPeak,
        RooRealVar& varWidth, double valWidth, double errLoWidth, double errHiWidth,
        // --- Standard Variables (Values extracted automatically) ---
        RooRealVar& varAlpha,
        RooRealVar& varNFull,
        RooRealVar& varN1st,
        RooRealVar& varN2nd,
        RooRealVar& varNBkg,
        RooRealVar& varConst,
        RooRealVar& varBeta
    );

}

#endif
