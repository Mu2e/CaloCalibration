#include "CaloCalibration/SourceCalib/inc/SourceFitter.hh"
#include "RooConstVar.h"
using namespace TMath;
using namespace RooFit;
using namespace CaloSourceCalib;
void SourceFitter::SetInitialGuesses(double peak, double width, double alpha,double beta,double evtsfull,double evtsfst,double evtsscd,double evtsbkg) {
     overridePeak_  = peak;
     overrideWidth_ = width;
     overrideAlpha_ = alpha;
     overridebeta_ = beta;
     overrideevtsfull_ = evtsfull;
     overrideevtsfst_ = evtsfst;
     overrideevtsscd_ = evtsscd;
     overrideevtsbkg_ = evtsbkg;
}

void SourceFitter::FitCrystal(TH1F* h_spec, TString opt, int crystalNo,  TTree *covar, Int_t &nEvents,Int_t &convergencestatus, Float_t &fpeak, Float_t &dpeak, Float_t &fsigma,Float_t &dsigma, Float_t &chiSq, Float_t &fstpeak,Float_t &fstsigma, Float_t &scdpeak,Float_t &scdsigma,Float_t &fcbalphaparam,Float_t &fcbndegparam,Float_t &Aparam,Float_t &Bparam, Float_t &Cparam, Float_t &fullResparam, Float_t &fstResparam,Float_t &scdResparam,Float_t &comCnstparam, Float_t &combetaparam, Float_t &fr_fullparam, Float_t &fr_frstparam,Float_t &fr_scndparam,Float_t &crystalNoparam,Float_t &fr_bkgparam,Float_t &errbar, Float_t &pval,Float_t &h_means,Float_t &h_stddevs, Float_t &unreducedchi2,Float_t &fval,Float_t &mparam,Float_t &etaparam,Int_t &ndof ){
    
  // set stlye optionsr
  gStyle -> SetOptFit(1111);
  gStyle -> SetOptStat(0);
  gStyle -> SetPadBottomMargin(0.125);
  gStyle -> SetPadTopMargin(0.075);
  gStyle -> SetPadLeftMargin(0.15);
  gStyle -> SetTitleOffset(1.0, "x");
  gStyle -> SetTitleOffset(1.75, "y");

  TCanvas *can = new TCanvas("can", "", 100, 100, 600, 600);
  can -> Draw();
  TString cryNum = to_string(crystalNo);
  TString oName = "mu2e_simu_fitSpec_"+ opt+"_" + cryNum + ".root";
  TString title = "SiPM " + cryNum;
  // Use fallback defaults
  double initPeak  = 98.08;
  double initWidth = 0.5;
  double initAlpha = 0.5;
  double initBeta = 1.0;
  double initEvtFull = 100000;
  double initEvtFst = 100000;
  double initEvtScd = 90000;
  double initEvtBkg= 50000;
  // Override if user provided new guesses
  if (overridePeak_  > 0) initPeak  = overridePeak_;
  if (overrideWidth_ > 0) initWidth = overrideWidth_;
  if (overrideAlpha_ > 0) initAlpha = overrideAlpha_;
  if (overridebeta_ > 0) initBeta = overridebeta_;
  if (overrideevtsfull_ > 0) initEvtFull = overrideevtsfull_;
  if (overrideevtsfst_ > 0) initEvtFst = overrideevtsfst_;
  if (overrideevtsscd_ > 0) initEvtScd = overrideevtsscd_;
  if (overrideevtsbkg_ > 0) initEvtBkg = overrideevtsbkg_;
  
  
  
  RooRealVar m_e("m_e", "electron energy in MeV", 0.511);//0.511
  RooRealVar crysADC("crysADC", "ADC [counts]", 40, 120);
  RooRealVar E0("E0", "energy offset [MeV]", 0.0);

  // Crystal Ball shape params
  RooRealVar fcbalpha("fcbalpha", "alpha", initAlpha, 0.1, 5.0);
  RooRealVar fcbndeg("fcbndeg", "n", 5);//, 1.5, 30.0);

  RooRealVar A("A", "A parameter",0.05, 0.001, 3.0);
  RooRealVar B("B", "B parameter", 0.005, 0.0001, 1.0);
  RooRealVar C("C", "C parameter",0.002, 0.0001, 1.0);
  // Peak parameters in mev

  RooRealVar fullPeak("fullPeak", "Full peak [ADC]", initPeak, 85, 108);
  RooFormulaVar eta("eta", "ADC/MeV", "fullPeak/ (6.13-E0)", RooArgSet(fullPeak, E0));
  RooFormulaVar m("m","Mev/ADC", "(6.13-E0)/fullPeak", RooArgSet(fullPeak, E0));
  RooFormulaVar Es("Es", "calibrated energy [MeV]", "m*crysADC + E0", RooArgSet(m, crysADC,E0));//Linear energy calibration
  RooFormulaVar fstEsPeak("fstEsPeak", "First escape", "fullPeak - m_e*eta", RooArgSet(fullPeak, m_e,eta));
  RooFormulaVar scdEsPeak("scdEsPeak", "Second escape", "fullPeak - (2*m_e)*eta", RooArgSet(fullPeak, m_e,eta));

  //------------------------------------------------------------------
  // Resolution for each peak (in energy units)
  //------------------------------------------------------------------
  /*RooFormulaVar fullRes("fullRes", "Full peak resolution",
                          "sqrt(pow(A,2)/sqrt(6.13)+pow(B,2)+pow(C/6.13,2))", RooArgSet(A, B, C));
  RooFormulaVar fstRes("fstRes", "First escape resolution",
                         "sqrt(pow(A,2)/sqrt(5.62)+pow(B,2)+pow(C/5.62,2))", RooArgSet(A, B, C));
  RooFormulaVar scdRes("scdRes", "Second escape resolution",
                         "sqrt(pow(A,2)/sqrt(5.11)+pow(B,2)+pow(C/5.11,2))", RooArgSet(A, B, C));*/

  RooRealVar fullWidth("fullWidth", "Full width [MeV]",initWidth,0.2,1.5);//5,2,20);
  RooRealVar Egamma("Egamma", "Full peak [MeV]",6.13);
  RooRealVar fstesc("fstesc", "first peak [MeV]",5.619);
  RooRealVar scdesc("scdesc", "second peak [MeV]",5.108);
  //three peak crystal ball function
  RooCBShape fullErg("fullErg", "Full peak", Es, Egamma, fullWidth, fcbalpha, fcbndeg);
  RooCBShape firsErg("firsErg", "Single escape", Es,fstesc, fullWidth, fcbalpha, fcbndeg);
  RooCBShape secdErg("secdErg", "Double escape", Es, scdesc, fullWidth, fcbalpha, fcbndeg);

  // Yield
  nEvents = h_spec->GetEntries();
  RooRealVar evtsFull("evtsFull", "Full peak yield", initEvtFull, 0, nEvents);// try making this larger than nevents-- try 2*nevents
  RooRealVar evtsFrst("evtsFrst", "First escape yield", initEvtFst, 0, nEvents);
  RooRealVar evtsScnd("evtsScnd", "Second escape yield", initEvtScd, 0, nEvents);
  RooRealVar evtsbkg("evtsbkg", "Background yield", initEvtBkg, 0, nEvents);
  // Background (logistic)
  RooRealVar comCnst("comCnst", "Background const",4.0);//, 0.0, 10.0);
  RooRealVar combeta("combeta", "Background beta", initBeta, 0.001, 200.0);
  RooRealVar bkg_fixed("bkg_fixed", "Background offset", 0.0);//, -10.0, 10.0);
  //logistic background pdf
  RooGenericPdf comPdf("comPdf", "logistic background",
    "pow(1.0 + exp((Es - comCnst)/combeta), -1.0) + bkg_fixed",
    RooArgSet(Es, comCnst, combeta, comCnst,bkg_fixed));
  //combing four pdfs
  RooAddPdf fitFun("fitFun", "Total model",
                     RooArgList(fullErg, firsErg, secdErg, comPdf),
                     RooArgList(evtsFull, evtsFrst, evtsScnd,evtsbkg));
                  
  // Fix normalization set to crysADC
  fitFun.fixCoefNormalization(RooArgSet(crysADC));
  //preparing RooPlot
  RooPlot *chFrame = crysADC.frame(Title(title));
  h_spec->Sumw2();
  RooDataHist chSpec("crysADC","crysADC", crysADC, h_spec);
    // Width of the first bin (in keV if the axis is in keV)
	double binWidth = h_spec->GetXaxis()->GetBinWidth(500);  //remove before push
	std::cout << "Bin width = " << binWidth << " keV" << std::endl; //remove before push

  if(opt == "chi2"){ //binned chi2 fit
    RooFitResult *fitRes = fitFun.chi2FitTo(chSpec, Range(40,115.2),  Strategy(2),MaxCalls(10000),PrintLevel(1),Save(),DataError(RooAbsData::SumW2),Hesse(kTRUE));//,Minimizer("Fumili") Minos(kTRUE)
    fitRes->Print("v");
    fval = fitFun.chi2FitTo(chSpec, Range(40,115.2), Save())->minNll();
    convergencestatus =fitRes->status(); 
   }  
  if(opt == "nll"){ //binned nll fit
    RooAbsReal* nll = fitFun.createNLL(chSpec, Range(40,115.2));
    RooMinimizer m(*nll);
    m.migrad();  
    m.hesse();
    RooFitResult *nllRes = m.save();
    nll->Print("v");
	convergencestatus =nllRes->status(); 
	fval = nll->getVal();
  }
  //getting unreduced chisq
  RooChi2Var unred_chi2("chi2","unreduced chi2", fitFun, chSpec);
  unreducedchi2 = unred_chi2.getVal();  
  RooFitResult* fitResult = fitFun.fitTo(chSpec, Save());
  ndof = chSpec.numEntries() - fitResult->floatParsFinal().getSize();

  // plot components
  chSpec.plotOn(chFrame, MarkerColor(kBlack), LineColor(kBlack), MarkerSize(0.5), Name("chSpec"));
  fitFun.plotOn(chFrame, LineColor(kRed), LineStyle(1), Name("fit"));
  chiSq = chFrame->chiSquare(11);
  fitFun.plotOn(chFrame, Components(fullErg), LineColor(kOrange), LineStyle(5), Name("main"));
  fitFun.plotOn(chFrame, Components(firsErg), LineColor(kViolet), LineStyle(5), Name("fescape"));
  fitFun.plotOn(chFrame, Components(secdErg), LineColor(kCyan), LineStyle(5), Name("sescape"));
  fitFun.plotOn(chFrame, Components(comPdf), LineColor(kBlue), LineStyle(5), Name("background"));

  fpeak = fullPeak.getVal();
  dpeak = fullPeak.getError();
  fsigma = fullWidth.getVal();
  dsigma = fullWidth.getError();
  fstpeak = fstEsPeak.getVal();
  fstsigma = 0;//fstWidth.getVal();
  scdpeak = scdEsPeak.getVal();
  scdsigma = 0;//scdWidth.getVal();
  fcbalphaparam = fcbalpha.getVal();
  fcbndegparam = fcbndeg.getVal();
  Aparam = A.getVal();
  Bparam = B.getVal();
  Cparam = C.getVal();
  fullResparam = 0;//fullRes.getVal();
  fstResparam = 0;//fstRes.getVal();
  scdResparam = 0;//scdRes.getVal();
  comCnstparam =comCnst.getVal();
  combetaparam = combeta.getVal();
  fr_fullparam = evtsFull.getVal()/nEvents;
  fr_frstparam = evtsFull.getVal()/nEvents;
  fr_scndparam = evtsScnd.getVal()/nEvents;
  fr_bkgparam= evtsbkg.getVal()/nEvents;
  crystalNoparam = crystalNo;
  errbar = (1/(fpeak/6.13))*(dpeak/fpeak);                   
  pval = TMath::Prob(chiSq*11, 11);
  mparam = m.getVal();
  etaparam = eta.getVal();
  //make pretty plots
  TPaveLabel *ptitle = new TPaveLabel(0.80, 0.90, 0.85, 0.80, Form("Mu2e Simulation"), "brNDC");
  ptitle -> SetFillStyle(0);
  ptitle -> SetBorderSize(0);
  ptitle -> SetTextSize(0.4);
  ptitle -> SetTextColor(kBlack);
  ptitle -> SetTextFont(72);
  ptitle -> SetFillColor(kWhite);
  chFrame -> addObject(ptitle);
  TPaveLabel *pnentres = new TPaveLabel(0.15, 0.85, 0.25, 0.75, Form("Nentries = %4.2i", nEvents), "brNDC");
  pnentres -> SetFillStyle(0);
  pnentres -> SetBorderSize(0);
  pnentres -> SetTextSize(0.4);
  pnentres -> SetTextFont(42);
  pnentres -> SetTextColor(kBlack);
  pnentres -> SetFillColor(kWhite);
  chFrame -> addObject(pnentres);
  TPaveLabel *pchi2 = new TPaveLabel(0.15, 0.75, 0.25, 0.65, Form("#chi^{2}/ndf = %4.2f", chiSq), "brNDC");
  pchi2 -> SetFillStyle(0);
  pchi2 -> SetBorderSize(0);
  pchi2 -> SetTextSize(0.4);
  pchi2 -> SetTextFont(42);
  pchi2 -> SetTextColor(kBlack);
  pchi2 -> SetFillColor(kWhite);
  chFrame -> addObject(pchi2);
  TPaveLabel *fpk = new TPaveLabel(0.15, 0.65, 0.25, 0.55, Form("#mu_{main} = %.2f#pm%.2f", fpeak, dpeak), "brNDC");
  fpk -> SetFillStyle(0);
  fpk -> SetBorderSize(0);
  fpk -> SetTextSize(0.4);
  fpk -> SetTextFont(42);
  fpk -> SetTextColor(kBlack);
  fpk -> SetFillColor(kWhite);
  chFrame -> addObject(fpk);
  TPaveLabel *fsg = new TPaveLabel(0.15, 0.55, 0.25, 0.45, Form("#sigma_{main} =%.2f#pm%.4f", fsigma,dsigma), "brNDC");
  fsg -> SetFillStyle(0);
  fsg -> SetBorderSize(0);
  fsg -> SetTextSize(0.4);
  fsg -> SetTextFont(42);
  fsg -> SetTextColor(kBlack);
  fsg -> SetFillColor(kWhite);
  chFrame -> addObject(fsg);

  // Create top pad for fit
	TPad *pad1 = new TPad("pad1", "Top pad", 0, 0.25, 1, 1.0);
	pad1->SetBottomMargin(0.035); 
	pad1->Draw();
	pad1->cd();  // switch to top pad
  chFrame -> SetYTitle("Events per 0.4 keV");
  chFrame -> GetYaxis()->SetTitleOffset(1.0);
  chFrame -> GetYaxis()->SetRangeUser(0, 5000);
  chFrame -> Draw();
  TLegend* legend = new TLegend(0.5, 0.7);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry("main", "main (eqv. 6.13MeV)", "L");
  legend->AddEntry("fescape", "first escape", "L");
  legend->AddEntry("sescape", "second escape", "L");
  legend->AddEntry("background", "background", "L");
  legend->Draw();
  // Back to canvas
	can->cd();
	// Create bottom pad for residuals
	TPad *pad2 = new TPad("pad2", "Bottom pad", 0, 0.0, 1, 0.25);
	pad2->SetTopMargin(0.05);
	pad2->SetBottomMargin(0.3); // room for X-axis labels
	pad2->Draw();
	pad2->cd();
// Make residual histogram
double xMin = 40.0;
double xMax = h_spec->GetXaxis()->GetXmax();
int nBins   = h_spec->FindBin(xMax) - h_spec->FindBin(xMin) + 1;
int startBin = h_spec->FindBin(xMin);
int endBin = h_spec->FindBin(xMax);
double totalYield = h_spec->Integral(startBin,endBin); // total events in the histogram
TH1F* hresidual = new TH1F("hresidual","", nBins, xMin, xMax);

for (int i = startBin; i <= h_spec->GetNbinsX(); ++i) {
    double x     = h_spec->GetBinCenter(i);
    double yData = h_spec->GetBinContent(i);
    double yErr  = h_spec->GetBinError(i);
    double binW  = h_spec->GetBinWidth(i);
    RooArgSet vars(crysADC);
    crysADC.setVal(x);
    double muFit = fitFun.getVal(&vars) * (totalYield) * binW;
    double res = (yErr > 0.0) ? (yData - muFit)/ yErr : 0.0;
    // shift bin index so we start from 1 inside hresidual
    hresidual->SetBinContent(i - startBin + 1, res);
}
// Draw residual histogram
hresidual->SetTitle("");
hresidual->GetYaxis()->SetTitle("Normalized Residuals (data - fit)/#sigma");
hresidual->GetYaxis()->SetTitleSize(0.12);
hresidual->GetYaxis()->SetLabelSize(0.10);
hresidual->GetXaxis()->SetTitleSize(0.12);
hresidual->GetXaxis()->SetLabelSize(0.10);
hresidual->Draw("HIST");

  can -> SaveAs(oName); 
  can->Close();  
  delete can;    
  can = nullptr; 

  covar->Fill();
/*// ---- 2D contour scan ----
// Choose parameters to scan
// Build the NLL orc Chi2
//RooAbsReal* nll = fitFun.createNLL(chSpec);
RooChi2Var fitres("chi2", "chi2", fitFun, chSpec, DataError(RooAbsData::SumW2));

// Minuit interface
//RooMinimizer minimizer(*nll);
RooMinimizer minimizer(fitres);

int nX = 100;  // number of bins in X
int nY = 100;  // number of bins in Y
int minBinX, minBinY, minBinZ;

double xMincontour = 85;   // range of fullPeak
double xMaxcontour = 100;
double yMincontour = 0.6;  // range of fullWidth
double yMaxcontour = 0.8;
// histogram
TH2F* h2Contour = new TH2F("h2Contour", "2D scan of Chi2",//NLL",
                           nX, xMincontour, xMaxcontour, nY, yMincontour, yMaxcontour);

//double nllMin = 1e30;
double chi2Min = 1e30;

for (int ix = 0; ix < nX; ++ix) {
    double xVal = xMincontour + ix * (xMaxcontour - xMincontour) / (nX - 1);
    fullPeak.setVal(xVal);
    fullPeak.setConstant(true);   // FIX parameter during scan

    for (int iy = 0; iy < nY; ++iy) {
        double yVal = yMincontour + iy * (yMaxcontour - yMincontour) / (nY - 1);
        fullWidth.setVal(yVal);
        fullWidth.setConstant(true);  // FIX parameter during scan

        // Minimize over all other free parameters
        minimizer.migrad();
        minimizer.hesse();
        minimizer.setStrategy(2);
        //minimizer.minos();
        

        //double nllVal = nll->getVal();
        //if (!std::isfinite(nllVal)) nllVal = 1e6;

        //if (nllVal < nllMin) nllMin = nllVal;
        //h2Contour->SetBinContent(ix + 1, iy + 1, nllVal);
        
        double chi2Val = fitres.getVal();
        if (!std::isfinite(chi2Val)) chi2Val = 1e6;

        if (chi2Val < chi2Min) chi2Min = chi2Val;
        h2Contour->SetBinContent(ix + 1, iy + 1, chi2Val);//add end commend

        fullWidth.setConstant(false);  // unfix for next iteration
    }
    fullPeak.setConstant(false);  // unfix for next iteration
}

// normalize to ?NLL
for (int ix = 1; ix <= nX; ++ix) {
    for (int iy = 1; iy <= nY; ++iy) {
        double val = h2Contour->GetBinContent(ix, iy);
        if (std::isfinite(val)) {
            h2Contour->SetBinContent(ix, iy, val - chi2Min);//nllMin);
            //h2Contour->SetBinContent(ix, iy, val - nllMin);
        }
    }
}
//get minimum values
h2Contour->GetMinimumBin(minBinX, minBinY, minBinZ);
double bestX = h2Contour->GetXaxis()->GetBinCenter(minBinX);
double bestY = h2Contour->GetYaxis()->GetBinCenter(minBinY);

// Draw contour
TCanvas* c2D = new TCanvas("c2D", "2D Contour", 800, 600);
h2Contour->Draw("COLZ");

TMarker* gridBest = new TMarker(bestX, bestY, 3);
gridBest->SetMarkerColor(kGreen);
gridBest->SetMarkerStyle(34);
gridBest->SetMarkerSize(1.5);

// --- Uncertainty crosshair (±1?) ---
TLine* lHor = new TLine(fpeak - dpeak, fsigma, fpeak + dpeak, fsigma); // horizontal line
TLine* lVer = new TLine(fpeak, fsigma - dsigma, fpeak, fsigma + dsigma); // vertical line

for (auto L : {lHor, lVer}) {
    L->SetLineColor(kMagenta + 1);
    L->SetLineStyle(1);
    L->SetLineWidth(4);
    L->Draw("SAME");
}
// Clone the original 2D histogram
TH2F* h2Contour1 = (TH2F*)h2Contour->Clone("h2Contour1");

// Define the contour level (??² = 1)
double level1[1] = {1.0};  // 1-sigma for 1 parameter (??² = 1)
h2Contour1->SetContour(1, level1);

// Style settings
h2Contour1->SetLineColor(kCyan);
h2Contour1->SetLineWidth(3);
h2Contour1->SetLineStyle(1);

// Draw the contour on top of the existing color plot
h2Contour1->Draw("CONT3 SAME");
gridBest->Draw("SAME");
TLegend* leg = new TLegend(0.55, 0.75, 0.85, 0.9);
leg->SetBorderSize(0);
leg->SetFillStyle(0);
leg->AddEntry(gridBest, "Best Fit (grid)", "p");
leg->AddEntry(lVer, "Best Fit & Errors (standard fit)", "l");
//leg->AddEntry(h2Copy, "??^{2} = 1, 2.3, 4", "l");
leg->Draw();

//gPad->SetLogz();
h2Contour->GetZaxis()->SetRangeUser(0, 5);
h2Contour->SetXTitle("fullPeak");
h2Contour->SetYTitle("fullWidth");
h2Contour->GetZaxis()->SetTitle("#Delta#chi^{2}");
c2D->SaveAs("Contour_fullPeak_fullwidth.root");*/

/*//---------------------------------------------------------
// Refit with MINOS errors and save new plot
//---------------------------------------------------------

//TCanvas* hessecan = (TCanvas*)gROOT->FindObject("can");
//TCanvas *can_minos =(TCanvas*)hessecan->Clone("can_minos");
//can_minos->SetTitle("New comparative canvas");
//can_minos->Draw();
//TString minosName = "mu2e_minos_hesse"+ opt+"_" + cryNum + ".root";
TCanvas *can_minos = new TCanvas("can_minos", "Fit with MINOS errors", 100, 100, 600, 600);
TString minosName = "mu2e_minos_hesse" + opt + "_" + cryNum + ".root";
can_minos->cd();



// Redo the same model fit, but enable MINOS
if (opt == "chi2") {
    RooFitResult *fitResMinos = fitFun.chi2FitTo(chSpec,
                                   Range(40, 115.2),
                                   Strategy(0),        // higher strategy for robustness
                                   MaxCalls(10000),
                                   PrintLevel(1),
                                   Save(),
                                   DataError(RooAbsData::SumW2),
                                   Hesse(kTRUE),       // ?? Enable MINOS error evaluation
                                   Minos(kTRUE));
    fitResMinos->Print("v"); 
}
else if (opt == "nll") {
    RooAbsReal* nll = fitFun.createNLL(chSpec, Range(40,115.2));
    RooMinimizer m2(*nll);
    m2.migrad();
    m2.hesse();
    m2.minos();  // ?? Run MINOS explicitly
    //fitResMinos = m2.save();
}


// Plot the fit with MINOS error bands
RooPlot *frame_minos = crysADC.frame(Title("Fit with MINOS errors"));
chSpec.plotOn(frame_minos, MarkerColor(kBlack), LineColor(kBlack), MarkerSize(0.5));
//fitFun.plotOn(frame_minos, VisualizeError(fitResMinos, crysADC, 1, kFALSE), FillColor(kOrange), FillStyle(3001)); // 1? band
//fitFun.plotOn(frame_minos, LineColor(kRed));

 chSpec.plotOn(frame_minos, MarkerColor(kBlack), LineColor(kBlack), MarkerSize(0.5), Name("chSpec"));
 fitFun.plotOn(frame_minos, LineColor(kRed+2), LineStyle(1), Name("fit"));
  //chiSq = frame_minos->chiSquare(11);
  fitFun.plotOn(frame_minos, Components(fullErg), LineColor(kOrange+2), LineStyle(1), Name("main"));
  fitFun.plotOn(frame_minos, Components(firsErg), LineColor(kViolet+2), LineStyle(1), Name("fescape"));
  fitFun.plotOn(frame_minos, Components(secdErg), LineColor(kCyan+2), LineStyle(1), Name("sescape"));
  fitFun.plotOn(frame_minos, Components(comPdf), LineColor(kBlue+2), LineStyle(1), Name("background"));
// Draw and save
frame_minos->Draw();
can_minos->SaveAs(minosName);
can_minos->Close();
delete can_minos;
can->Close(); 
delete can;    // Delete the object
can = nullptr; // Safety (optional)*/

}
