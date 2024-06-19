#ifndef _MakeAnalysisTree_hh
#define _MakeAnalysisTree_hh

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
  class MakeAnalysisTree  {
      public:
        explicit MakeAnalysisTree(){};
        explicit MakeAnalysisTree(const MakeAnalysisTree &){};
        MakeAnalysisTree& operator = (const MakeAnalysisTree &);
        virtual ~MakeAnalysisTree() = default;
        #ifndef __CINT__
        void MakeAnalysisTreeOutputs();
        #endif
        ClassDef (MakeAnalysisTree,1);
    };
}
#endif /* MakeAnalysisTree.hh */
