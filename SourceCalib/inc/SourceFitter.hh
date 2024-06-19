#ifndef _SourceFitter_hh
#define _SourceFitter_hh

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
#include <Riostream.h>

using namespace std;
using namespace TMath;
using namespace RooFit;

namespace CaloSourceCalib{
  class SourceFitter  {
      public:
        explicit SourceFitter(){};
        explicit SourceFitter(const SourceFitter &){};
        SourceFitter& operator = (const SourceFitter &);
        virtual ~SourceFitter() = default;
        #ifndef __CINT__
        void SourceFitterOutputs();
        #endif
        ClassDef (SourceFitter,1);
    };
}
#endif /* SourceFitter.hh */
