/*#ifndef validationfitter_hh
#define validationfitter_hh

// Standard C++ libraries for file I/O and strings
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>

// ROOT and RooFit libraries
#include "TH1.h"
#include "TH1F.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooCBShape.h"
#include "RooGenericPdf.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "TCanvas.h"
#include "TAxis.h"

namespace CaloSourceCalib {

    class ValidationFitter {
    public:
        // Notice the 3rd argument is now std::string txtFilePath
        static void ValidateCrystal(TH1F* h_spec, int crystalNo, std::string txtFilePath);
    };

}

#endif*/
