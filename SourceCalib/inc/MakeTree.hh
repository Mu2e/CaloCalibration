#ifndef _MakeTree_hh
#define _MakeTree_hh

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
  class MakeTree  {
      public:
        explicit MakeTree(){};
        explicit MakeTree(const MakeTree &){};
        MakeTree& operator = (const MakeTree &);
        virtual ~MakeTree() = default;
        #ifndef __CINT__
        void MakeTreeOutputs();
        #endif
        ClassDef (MakeTree,1);
    };
}
#endif /* MakeTree.hh */
