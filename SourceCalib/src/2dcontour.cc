// =====================================================================
// 2dcontour.cc
// 
// Implementation of the 2D Contour Plotter for Source Calibration.
// This function performs a 2D scan of the likelihood (or Chi2) surface
// between any two selected variables (e.g., Peak vs Width, Alpha vs Peak).
//
// Key Features:
// 1. "Switch Board": Selects variables dynamically via string names.
// 2. Auto-Ranging: Automatically determines scan window based on fit errors.
// 3. Robust Plotting: Handles physical limits and draws 1-sigma contours.
// =====================================================================

#include "CaloCalibration/SourceCalib/inc/2dcontour.hh"
#include "RooChi2Var.h"
#include "RooNLLVar.h"
#include "RooMinimizer.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TMarker.h"
#include "TLine.h"
#include "TLegend.h"
#include "TStyle.h"
#include <iostream>

using namespace RooFit;

void CaloSourceCalib::MakeContourPlot(
    RooAbsPdf& fitFun, 
    RooDataHist& chSpec, 
    TString fitType, 
    int cryNum,
    // --- Selection Strings ---
    // These strings determine which variables are plotted on X and Y.
    TString xName, 
    TString yName,
    // --- Special Variables ---
    // Manual values passed for complex variables like Peak/Width
    RooRealVar& varPeak,  double valPeak,  double errLoPeak,  double errHiPeak,
    RooRealVar& varWidth, double valWidth, double errLoWidth, double errHiWidth,
    // --- Standard Variables ---
    // Direct references to RooFit variables for other parameters
    RooRealVar& varAlpha,
    RooRealVar& varNFull,
    RooRealVar& varN1st,
    RooRealVar& varN2nd,
    RooRealVar& varNBkg,
    RooRealVar& varConst,
    RooRealVar& varBeta
) {
    // ======================================================
    // 1. Internal Switch Board
    // ======================================================
    // This section maps the input strings (xName, yName) to the actual
    // variable pointers and values.
    
    RooRealVar* pVarX = nullptr;
    RooRealVar* pVarY = nullptr;
    double valX = 0, errXLo = 0, errXHi = 0;
    double valY = 0, errYLo = 0, errYHi = 0;

    // Helper lambda to select variables based on name
    auto selectParam = [&](TString name, RooRealVar*& ptr, double& val, double& eLo, double& eHi) {
        
        // --- CASE A: Special Manual Variables ---
        // These use the high-precision calculated values passed from the Fitter
        if (name == "Peak") {
            ptr = &varPeak; val = valPeak; eLo = errLoPeak; eHi = errHiPeak;
            return;
        }
        if (name == "Width") {
            ptr = &varWidth; val = valWidth; eLo = errLoWidth; eHi = errHiWidth;
            return;
        }

        // --- CASE B: Standard RooFit Variables ---
        // These map directly to the RooRealVar objects
        if      (name == "Alpha")  ptr = &varAlpha;
        else if (name == "N_Full") ptr = &varNFull;
        else if (name == "N_1st")  ptr = &varN1st;
        else if (name == "N_2nd")  ptr = &varN2nd;
        else if (name == "N_Bkg")  ptr = &varNBkg;
        else if (name == "Const")  ptr = &varConst;
        else if (name == "Beta")   ptr = &varBeta;

        // Extract values automatically for Case B
        if (ptr) {
            val = ptr->getVal();
            // Try asymmetric Minos errors first
            eLo = ptr->getErrorLo();
            eHi = ptr->getErrorHi();
            
            // If Minos didn't run (errors are 0), fallback to symmetric Hessian error
            if (eLo == 0 && eHi == 0) {
                double sym = ptr->getError();
                eLo = -sym; 
                eHi = sym;
            }
        } else {
            std::cout << "[Contour] ERROR: Unknown parameter '" << name << "'" << std::endl;
        }
    };

    // --- Apply Selection ---
    selectParam(xName, pVarX, valX, errXLo, errXHi);
    selectParam(yName, pVarY, valY, errYLo, errYHi);

    // Validate selection
    if (!pVarX || !pVarY) {
        std::cout << "[Contour] Skipping: Invalid parameters selected." << std::endl;
        return;
    }

    // ======================================================
    // 2. Setup Minimizer
    // ======================================================
    // Construct the objective function (NLL or Chi2) based on user choice.
    
    TString outName = Form("Contour_cry%d_%s_%s_vs_%s.root", cryNum, fitType.Data(), xName.Data(), yName.Data());
    std::cout << "[Contour] Scanning " << xName << " vs " << yName << " -> " << outName << std::endl;

    RooAbsReal* fitFunc = nullptr;
    double contourLevel = 1.0;
    TString zTitle = "";

    if (fitType == "nll") {
        fitFunc = fitFun.createNLL(chSpec, Range(40, 115.2), Extended(true), Verbose(false));
        contourLevel = 0.5;   // Delta NLL = 0.5 corresponds to 1 sigma
        zTitle = "#Delta NLL";
    } else {
        fitFunc = new RooChi2Var("chi2", "chi2", fitFun, chSpec, DataError(RooAbsData::SumW2), Range(40, 115.2));
        contourLevel = 1.0;   // Delta Chi2 = 1.0 corresponds to 1 sigma
        zTitle = "#Delta #chi^{2}";
    }

    RooMinimizer minimizer(*fitFunc);
    minimizer.setPrintLevel(-1);
    minimizer.setStrategy(2); // High precision strategy

    // ======================================================
    // 3. Scan Configuration (Auto-Ranging)
    // ======================================================
    // We determine the scan window by taking the best fit value +/- 4 sigma.
    
    int nX = 50, nY = 50;
    double rangeScale = 4.0; 

    // Safety fallback: If errors are zero (fixed param), default to 10% of value
    double safeErrX_lo = (fabs(errXLo) > 1e-6) ? fabs(errXLo) : 0.1 * fabs(valX);
    double safeErrX_hi = (fabs(errXHi) > 1e-6) ? fabs(errXHi) : 0.1 * fabs(valX);
    double safeErrY_lo = (fabs(errYLo) > 1e-6) ? fabs(errYLo) : 0.1 * fabs(valY);
    double safeErrY_hi = (fabs(errYHi) > 1e-6) ? fabs(errYHi) : 0.1 * fabs(valY);

    double xMin = valX - rangeScale * safeErrX_lo;
    double xMax = valX + rangeScale * safeErrX_hi;
    double yMin = valY - rangeScale * safeErrY_lo;
    double yMax = valY + rangeScale * safeErrY_hi;

    // Clamp to physical limits if they exist in the RooRealVar
    if (pVarX->hasMin() && xMin < pVarX->getMin()) xMin = pVarX->getMin();
    if (pVarX->hasMax() && xMax > pVarX->getMax()) xMax = pVarX->getMax();
    if (pVarY->hasMin() && yMin < pVarY->getMin()) yMin = pVarY->getMin();
    if (pVarY->hasMax() && yMax > pVarY->getMax()) yMax = pVarY->getMax();

    TH2F* h2Contour = new TH2F("h2Contour", Form("Crystal %d: %s vs %s", cryNum, xName.Data(), yName.Data()),
                               nX, xMin, xMax, nY, yMin, yMax);

    // ======================================================
    // 4. The Grid Loop
    // ======================================================
    // Perform the actual 2D scan. We fix X and Y, then minimize everything else.
    
    double minVal = 1e30;
    
    for (int ix = 0; ix < nX; ++ix) {
        double xv = xMin + ix * (xMax - xMin) / (nX - 1);
        pVarX->setVal(xv);
        pVarX->setConstant(true); // Fix X

        for (int iy = 0; iy < nY; ++iy) {
            double yv = yMin + iy * (yMax - yMin) / (nY - 1);
            pVarY->setVal(yv);
            pVarY->setConstant(true); // Fix Y

            minimizer.migrad(); // Minimize free parameters
            double val = fitFunc->getVal();
            
            // Handle convergence failures
            if (!std::isfinite(val)) val = 1e6;

            if (val < minVal) minVal = val;
            h2Contour->SetBinContent(ix + 1, iy + 1, val);

            pVarY->setConstant(false); // Release Y
        }
        pVarX->setConstant(false); // Release X
    }

    // ======================================================
    // 5. Draw
    // ======================================================
    // Normalize the plot (subtract minVal) so the Z-axis starts at 0.
    for (int i = 1; i <= nX*nY; ++i) {
        double v = h2Contour->GetBinContent(i);
        if (std::isfinite(v)) h2Contour->SetBinContent(i, v - minVal);
    }

    TCanvas* c2D = new TCanvas("c2D", "Contour", 800, 600);
    h2Contour->Draw("COLZ");

    // --- Draw Grid Minimum Marker ---
    int bx, by, bz;
    h2Contour->GetMinimumBin(bx, by, bz);
    TMarker* mBest = new TMarker(h2Contour->GetXaxis()->GetBinCenter(bx), 
                                 h2Contour->GetYaxis()->GetBinCenter(by), 34);
    mBest->SetMarkerColor(kGreen); mBest->SetMarkerSize(1.5); mBest->Draw("SAME");

    // --- Draw Fit Result Crosshairs ---
    // These show the best fit value +/- 1 sigma from the main fit
    TLine* lH = new TLine(valX - fabs(errXLo), valY, valX + fabs(errXHi), valY);
    TLine* lV = new TLine(valX, valY - fabs(errYLo), valX, valY + fabs(errYHi));
    lH->SetLineColor(kMagenta); lH->SetLineWidth(3); lH->Draw("SAME");
    lV->SetLineColor(kMagenta); lV->SetLineWidth(3); lV->Draw("SAME");

    // --- Draw Contour Line ---
    // This draws the specific Chi2=1 or NLL=0.5 iso-line
    TH2F* hC = (TH2F*)h2Contour->Clone("hC");
    double lvl[1] = { contourLevel };
    hC->SetContour(1, lvl);
    hC->SetLineColor(kCyan); hC->SetLineWidth(3); hC->Draw("CONT3 SAME");

    // --- Legend and Titles ---
    TLegend* leg = new TLegend(0.55, 0.75, 0.89, 0.89);
    leg->AddEntry(mBest, "Grid Min", "p");
    leg->AddEntry(lH, "Fit Result #pm 1#sigma", "l");
    leg->AddEntry(hC, Form("%s = %.1f", zTitle.Data(), contourLevel), "l");
    leg->Draw();

    h2Contour->GetZaxis()->SetRangeUser(0, 5); // Zoom Z-axis to relevant region
    h2Contour->SetXTitle(pVarX->GetTitle());
    h2Contour->SetYTitle(pVarY->GetTitle());
    h2Contour->GetZaxis()->SetTitle(zTitle);

    c2D->SaveAs(outName);

    // Cleanup
    delete c2D; delete h2Contour; delete hC; delete mBest; delete lH; delete lV; delete leg;
    if (fitFunc) delete fitFunc;
}
