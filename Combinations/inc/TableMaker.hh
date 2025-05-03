#ifndef _TableMaker_hh
#define _TableMaker_hh

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
  class TableMaker  {
      public:
        explicit TableMaker(){};
        explicit TableMaker(const TableMaker &){};
        TableMaker& operator = (const TableMaker &);
        virtual ~TableMaker() = default;
        #ifndef __CINT__
        void TableMakerOutputs();
        #endif
        ClassDef (TableMaker,1);
    };
}
#endif /* TableMaker.hh */
