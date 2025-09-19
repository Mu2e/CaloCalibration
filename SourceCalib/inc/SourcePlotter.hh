#ifndef _SourcePlotter_hh
#define _SourcePlotter_hh

#include <fstream>
#include <iostream>
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
#include "TAxis.h"
#include <Riostream.h>
#include <vector>

#include "TPaveText.h"
#include "TLegend.h"
#include "TPaveLabel.h"
#include "TAttFill.h"
using namespace std;
using namespace TMath;
using namespace RooFit;

namespace CaloSourceCalib{
  class SourcePlotter  {
      public:
        explicit SourcePlotter(){};
        explicit SourcePlotter(const SourcePlotter &){};
        SourcePlotter& operator = (const SourcePlotter &);
        virtual ~SourcePlotter() = default;
        #ifndef __CINT__
        void ParamPlots(TTree* t, TFile *inputFile, TFile *outputFile,int cry_start, int cry_end);        
        #endif
        ClassDef (SourcePlotter,1);
    };
}
#endif /* SourcePlotter.hh */
