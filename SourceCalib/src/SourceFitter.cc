#include "CaloCalibration/SourceCalib/inc/SourceFitter.hh"
#include "CaloCalibration/SourceCalib/inc/2dcontour.hh"
#include "RooFFTConvPdf.h"
#include "RooExponential.h"
#include "RooCrystalBall.h"
#include "TF1.h"
/*
Loops over all historgrams to fit data for each SiPM/crystal.
The code run the fit 4 iterative ways to find a fit with both successful convergence status (either with asymmetric or symmetric errors) AND reasonable chi2:
(1) with manually set initial param guesses,
(2) mean of successful fits as initial params,
(3) randomized guesses within param boundaries (40 times)
(4) if above has not procuded desrired assymetic errors AND chi2; fit will fall back on symmetric errors with reasonable chi2 and run the fit one last time 


Final params are saved in arxiv.root


Command line examples are in READ.ME
*/
using namespace TMath;
using namespace RooFit;
using namespace CaloSourceCalib;

int SourceFitter::nSecondFits = 0;
int SourceFitter::nThirdFits = 0;
std::vector<int> SourceFitter::crystalsSecondFit;
std::vector<int> SourceFitter::crystalsThirdFit;
int SourceFitter::nFirstFitConverged = 0;
int SourceFitter::nSecondFitConverged = 0;
int SourceFitter::nThirdFitConverged = 0;
std::vector<int> SourceFitter::crystalsSecondFitConverged;
std::vector<int> SourceFitter::crystalsThirdFitConverged;
std::map<int,int> SourceFitter::thirdFitRetryCount;
std::vector<std::pair<int,int>> SourceFitter::convFailures;
int SourceFitter::nAsymErrors = 0;
std::vector<int> SourceFitter::crystalsWithAsymErrors;
int SourceFitter::badchi2 = 0;
std::vector<int> SourceFitter::crystalswithbadchi2;
std::map<int,float> SourceFitter::badchi2Values;
int SourceFitter::nHesseFallbacks = 0;
std::vector<int> SourceFitter::crystalsHesseFallback;
struct CovarAccumulator {
   int count = 0;
   double sumPeak = 0.0;
   //double sumWidth = 0.0;
   double sumAlpha = 0.0;
   double sumn = 0.0;
   double sumAlphaR = 0.0;
   //double sumnR = 0.0;//comment for double sided
   double sumBeta = 0.0;
   double sumEvtFull = 0.0;
   double sumEvtFst = 0.0;
   double sumEvtScd = 0.0;
   double sumEvtBkg = 0.0;
   //double sumEvtCompton1 = 0;
   //double sumEvtCompton2 = 0;
   //double sumEvtCompton3 = 0;
} covarAcc;
// ---- SILENCE ROOFIT & MINUIT2 ----
static bool suppress_messages = [](){


   // 1. Silence RooFit INFO and DEBUG
   RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);


   // 2. Silence Minuit2 messages
   ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1); // GLOBAL SETTER
   ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(0); // optional


   // 3. Silence ROOT warnings (optional)
   gErrorIgnoreLevel = kError;   // hide Info + Warning


   return true;
}();
void SourceFitter::FitCrystal(TH1F* h_spec, TString opt, int crystalNo,  TTree *covar, Int_t &nEvents,Int_t &convergencestatus, Float_t &fpeak, Float_t &peakerrorhigh,Float_t &peakerrorlo,Float_t &redpeak, Float_t &fsigma,Float_t &widtherrorhigh,Float_t &widtherrorlo, Float_t &chiSq, Float_t &fstpeak, Float_t &scdpeak,Float_t &fcbalphaparam,Float_t &fcbndegparam,Float_t &comCnstparam, Float_t &combetaparam, Float_t &fr_fullparam, Float_t &fr_frstparam,Float_t &fr_scndparam,Float_t &crystalNoparam,Float_t &fr_bkgparam,Float_t &fr_comptonparam1,Float_t &fr_comptonparam2,Float_t &fr_comptonparam3, Float_t &pval,Float_t &h_means,Float_t &h_stddevs, Float_t &unreducedchi2,Float_t &fval,Float_t &mparam,Float_t &etaparam,Int_t &ndof,bool contour,TString xVar, TString yVar, Float_t &errbarhigh, Float_t &errbarlo,Float_t &evtfullerrorhigh,Float_t &evtfullerrorlo,Float_t &Esparam,Float_t &Aparam,Float_t &Aperr,Float_t &Bparam,Float_t &Berr,Float_t &Cparam,Float_t &Cerr,Float_t &fcbalphaRparam, Float_t &fcbndegRparam, bool islyso ){//Float_t &fr_comptonparam Float_t &fr_comptonparam1,Float_t &fr_comptonparam2,Float_t &fr_comptonparam3


std::cout << "**************************************************************\n";
std::cout << "Fitting crystal " << crystalNo << ".\n";
std::cout << "**************************************************************\n";
convergencestatus = -1; 
bool asymSuccess = false;    
double currentPeak  = 0.0;
//double currentWidth = 0.0;
double currentAlpha = 0.0;
double currentn= 0.0;
double currentAlphaR = 0.0;
//double currentnR= 0.0;
double currentBeta  = 0.0;
double currentEvtFull = 0.0;
double currentEvtFst  = 0.0;
double currentEvtScd  = 0.0;
double currentEvtBkg  = 0.0;
//double currentEvtCompton1  = 0.0;
//double currentEvtCompton2  = 0.0;
//double currentEvtCompton3  = 0.0;

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
 TString materialName = islyso ? "LYSO" : "CsI";
 TString oName = "mu2e_simu_fitSpec_"+ opt+"_"+ materialName+"_" + cryNum + ".root";
 TString title = "SiPM " + cryNum;


 double initPeak,initEvtFull,initEvtFst,initEvtScd,PeakLow,PeakHigh,EvtFullLow, EvtFullHigh,EvtFstLow, EvtFstHigh,EvtScdLow,EvtScdHigh,
 initEvtBkg,initBeta,initAlpha,initn,initAlphaR,initNdegR;
 //initEvtCompton1,initEvtCompton2,initEvtCompton3;;//,initWidth,
 int low = h_spec->FindBin(40);
 int high = h_spec->FindBin(115);
 //int low = h_spec->FindBin(140);
 //int high = h_spec->FindBin(220);
 int integral_evts = h_spec->Integral(low, high);
 nEvents = h_spec->GetEntries();
 float reducedchi2 = 0.0;
 int nPars = 8;
 if (islyso){
  // LYSO SPECIFIC PARAMETERS
  initPeak  = 95.08;
  PeakLow = 85;
  PeakHigh = 105;
  //initWidth = 0.5;
  initAlpha  = 0.8;
  initAlphaR = 3.0;//uncomment for double sided
  initNdegR  = 10;//uncomment for double sided
  initn = 10;
  //initn = 1;//for single sided
  initBeta = -0.09;//3.0;
  initEvtFull = 9000;
  EvtFullLow= 0.10*integral_evts;
  EvtFullHigh =  integral_evts;
  initEvtFst = 9000;
  EvtFstLow= 0.2*integral_evts;
  EvtFstHigh = 0.5*integral_evts;
  initEvtScd = 9000;
  EvtScdLow= 0.05*integral_evts;
  EvtScdHigh = 0.2*integral_evts;
  initEvtBkg= 5000;
  //initEvtCompton1= 5000;
  //initEvtCompton2= 5000;
  //initEvtCompton3= 5000;
 }
 else{
  // CsI SPECIFIC PARAMETERS
  initPeak  = 85;
  PeakLow = 85;
  PeakHigh = 108;
  //initWidth = 0.5;
  initAlpha  = 0.5;
  initAlphaR = 2.0;//uncomment for double sided
  initNdegR  = 10;//uncomment for double sided
  initn = 10;
  //initn = 1;//for single sided
  initBeta = -0.09;//3.0;
  initEvtFull = 9000;
  EvtFullLow= 0.05*integral_evts;//0.1
  EvtFullHigh =  0.6*integral_evts;
  //EvtFullHigh =  integral_evts;
  initEvtFst = 9000;
  EvtFstLow= 0.2*integral_evts;
  EvtFstHigh = integral_evts;
  //EvtFstHigh = 0.5*integral_evts;
  initEvtScd = 9000;
  EvtScdLow= 0.05*integral_evts;
  EvtScdHigh = integral_evts;
  initEvtBkg= 5000;
  //initEvtCompton1= 5000;
  //initEvtCompton2= 5000;
  //initEvtCompton3= 5000;
 }
  currentPeak  = initPeak;
  //currentWidth = initWidth;
  currentAlpha = initAlpha;
  currentn = initn;
  currentAlphaR = initAlphaR;
  //currentnR = initNdegR;
  currentBeta  = initBeta;
  currentEvtFull = initEvtFull;
  currentEvtFst  = initEvtFst;
  currentEvtScd  = initEvtScd;
  currentEvtBkg  = initEvtBkg;
  //currentEvtCompton1  = initEvtCompton1;
  //currentEvtCompton2  = initEvtCompton2;
  //currentEvtCompton3  = initEvtCompton3;
 RooRealVar m_e("m_e", "electron energy in MeV", 0.511);
 RooRealVar crysADC("crysADC", "ADC [counts]", 40, 115);
 //RooRealVar crysADC("crysADC", "ADC [counts]", 140, 220);
 RooRealVar E0("E0", "energy offset [MeV]", 0.0); //might need to remove if not 
 
 //RooRealVar fcbndeg("fcbndeg", "n", 5, 6, 7);//only for single sided
 //fcbndeg.setConstant(kTRUE);//only for single sided
 RooRealVar fullPeak("fullPeak", "Full peak [ADC]", initPeak, PeakLow, PeakHigh);
 RooFormulaVar eta("eta", "ADC/MeV", "fullPeak/ (6.13-E0)", RooArgSet(fullPeak, E0));
 RooFormulaVar m("m","Mev/ADC", "(6.13-E0)/fullPeak", RooArgSet(fullPeak, E0));
 RooFormulaVar Es("Es", "calibrated energy [MeV]", "m*crysADC + E0", RooArgSet(m, crysADC,E0));
 RooFormulaVar fstEsPeak("fstEsPeak", "First escape", "fullPeak - m_e*eta", RooArgSet(fullPeak, m_e,eta));
 RooFormulaVar scdEsPeak("scdEsPeak", "Second escape", "fullPeak - (2*m_e)*eta", RooArgSet(fullPeak, m_e,eta));


 //RooRealVar fullWidth("fullWidth", "Full width [MeV]",initWidth,0.2,1.5);
 //comment after test data test
 //RooRealVar fullWidth("fullWidth", "Full width [MeV]",0.5,0.2,1.5);
 RooRealVar Egamma("Egamma", "Full peak [MeV]",6.13);
 RooRealVar fstesc("fstesc", "first peak [MeV]",5.619);
 RooRealVar scdesc("scdesc", "second peak [MeV]",5.108);
 // Stochastic term: sigma = A * sqrt(E), where A [sqrt(MeV)] encodes photon collection efficiency
 RooRealVar A("A", "stochastic term [sqrt(MeV)]", 0.1, 0, 1.0);
 //RooRealVar B("B", "Constant term [(MeV)^2]", 0.1, 0, 1.0);
 //RooRealVar C("C", "Noise term", 0.1, 0, 1.0);
 
 //A and B term
 //RooFormulaVar fullWidth("fullWidth", "Full width [MeV]", "A*std::sqrt(Egamma) + B*Egamma", RooArgSet(A, Egamma, B));
 //RooFormulaVar fstwidth ("fstwidth",  "first width [MeV]",  "A*std::sqrt(fstesc) + B*fstesc",  RooArgSet(A, fstesc, B));
 //RooFormulaVar scdwidth ("scdwidth",  "second width [MeV]", "A*std::sqrt(scdesc) + B*scdesc", RooArgSet(A, scdesc, B));
 
 //B and C term
// RooFormulaVar fullWidth("fullWidth", "Full width [MeV]", "sqrt(B*Egamma + C)", RooArgSet(Egamma, B, C));
// RooFormulaVar fstwidth ("fstwidth",  "first width [MeV]",  "sqrt(B*fstesc + C)",  RooArgSet(fstesc, B, C));
 //RooFormulaVar scdwidth ("scdwidth",  "second width [MeV]", "sqrt(B*scdesc + C)", RooArgSet(scdesc, B, C));
 
 //only B term
 //RooFormulaVar fullWidth("fullWidth", "Full width [MeV]", "(B*Egamma)", RooArgSet(Egamma, B));
 //RooFormulaVar fstwidth ("fstwidth",  "first width [MeV]",  "(B*fstesc )",  RooArgSet(fstesc, B));
 //RooFormulaVar scdwidth ("scdwidth",  "second width [MeV]", "(B*scdesc )", RooArgSet(scdesc, B));
 
 // only A term
 RooFormulaVar fullWidth("fullWidth", "Full width [MeV]", "A*std::sqrt(Egamma)", RooArgSet(A, Egamma));
 RooFormulaVar fstwidth ("fstwidth",  "first width [MeV]",  "A*std::sqrt(fstesc)",  RooArgSet(A, fstesc));
 RooFormulaVar scdwidth ("scdwidth",  "second width [MeV]", "A*std::sqrt(scdesc)", RooArgSet(A, scdesc));
 // only C term-- not really needed bc sigma=c; so can just use the widths
// RooFormulaVar fullWidth("fullWidth", "Full width [MeV]", "C", RooArgSet(C, Egamma));
 //RooFormulaVar fstwidth ("fstwidth",  "first width [MeV]",  "C+fstesc",  RooArgSet(C, fstesc));
 //RooFormulaVar scdwidth ("scdwidth",  "second width [MeV]", "C+scdesc", RooArgSet(C, scdesc));
 
 //a b and c terma
 /*RooFormulaVar fullWidth("fullWidth", "Full width [MeV]", "A*std::sqrt(Egamma) + B*Egamma+C", RooArgSet(A, Egamma, B,C));
 RooFormulaVar fstwidth ("fstwidth",  "first width [MeV]",  "A*std::sqrt(fstesc) + B*fstesc+C",  RooArgSet(A, fstesc, B,C));
 RooFormulaVar scdwidth ("scdwidth",  "second width [MeV]", "A*std::sqrt(scdesc) + B*scdesc+C", RooArgSet(A, scdesc, B,C));*/


 //Single sided crystal ball w SAME width 
 //RooCBShape fullErg("fullErg", "Full peak", Es, Egamma, fullWidth, fcbalpha, fcbndeg);
// RooCBShape firsErg("firsErg", "Single escape", Es,fstesc, fullWidth, fcbalpha, fcbndeg);
 //RooCBShape secdErg("secdErg", "Double escape", Es, scdesc, fullWidth, fcbalpha, fcbndeg);
  //for single sided
 /*RooRealVar fcbalpha("fcbalpha", "alpha", initAlpha,0.1, 5.0);//0.01, 3.0);
 fcbalpha.setConstant(kFALSE);
 RooRealVar fcbndeg("fcbndeg", "n", initn, 5, 10);
 fcbndeg.setConstant(kFALSE);
 //Single sided crystal ball w different width using A/B/C
 RooCBShape fullErg("fullErg", "Full peak", Es, Egamma, fullWidth, fcbalpha, fcbndeg);
 RooCBShape firsErg("firsErg", "Single escape", Es,fstesc, fstwidth, fcbalpha, fcbndeg);
 RooCBShape secdErg("secdErg", "Double escape", Es, scdesc, scdwidth, fcbalpha, fcbndeg);*/

 
 //for Double sided
RooRealVar fcbalpha("fcbalpha", "alpha", initAlpha,0.1, 5.0);//0.01, 3.0);
 fcbalpha.setConstant(kTRUE);
 RooRealVar fcbndeg("fcbndeg", "n", initn, 5, 20);
 fcbndeg.setConstant(kFALSE);

 RooRealVar fcbalphaR("fcbalphaR", "alpha R", initAlphaR, 0.1, 3.5);
 fcbalphaR.setConstant(kFALSE);
 RooRealVar fcbndegR("fcbndegR", "n R", initNdegR, 1, 15);
 fcbndegR.setConstant(kTRUE);
 //Double sided crystal ball with Diferent widths and same alpha L/R an N l/r
 RooCrystalBall fullErg("fullErg", "Full peak", Es, Egamma, fullWidth, fcbalpha, fcbndeg, fcbalphaR, fcbndegR);
 //trying with one single sided cb
 //RooCBShape fullErg("fullErg", "Full peak", Es, Egamma, fullWidth, fcbalpha, fcbndeg);
 RooCrystalBall firsErg("firsErg", "Single escape", Es, fstesc, fstwidth, fcbalpha, fcbndeg, fcbalphaR, fcbndegR);
 RooCrystalBall secdErg("secdErg", "Double escape", Es, scdesc, scdwidth, fcbalpha, fcbndeg, fcbalphaR, fcbndegR);
/*
RooAbsPdf* fullErg = nullptr;
RooAbsPdf* firsErg = nullptr;
RooAbsPdf* secdErg = nullptr;
if (islyso) {
    // Instantiate Single-Sided Crystal Ball for LYSO
    fullErg = new RooCBShape("fullErg", "Full peak", Es, Egamma, fullWidth, fcbalpha, fcbndeg);
    firsErg = new RooCBShape("firsErg", "Single escape", Es, fstesc, fstwidth, fcbalpha, fcbndeg);
    secdErg = new RooCBShape("secdErg", "Double escape", Es, scdesc, scdwidth, fcbalpha, fcbndeg); 
}
else {
    // Instantiate Double-Sided Crystal Ball for CsI (Notice the fixed semicolon at the end!)
    fullErg = new RooCrystalBall("fullErg", "Full peak", Es, Egamma, fullWidth, fcbalpha, fcbndeg, fcbalphaR, fcbndegR);
    firsErg = new RooCrystalBall("firsErg", "Single escape", Es, fstesc, fstwidth, fcbalpha, fcbndeg, fcbalphaR, fcbndegR);
    secdErg = new RooCrystalBall("secdErg", "Double escape", Es, scdesc, scdwidth, fcbalpha, fcbndeg, fcbalphaR, fcbndegR);
}*/
 /*//unique alpha for full peak
 RooRealVar fcbalpha_full("fcbalpha_full", "alpha for full", 1.6,0.1,5.0);
 RooRealVar fcbndeg_full("fcbndeg_full", "n for full", 10);
 RooCBShape fullErg("fullErg", "Full peak", Es, Egamma, fullWidth, fcbalpha_full, fcbndeg_full);
 RooCBShape secdErg("secdErg", "Double escape", Es, scdesc, fullWidth, fcbalpha_full, fcbndeg_full);*/


 
  //Double sided crystal ball with same widths and alpha L/R an N l/r
 //RooCrystalBall fullErg("fullErg", "Full peak", Es, Egamma, fullWidth, fcbalpha, fcbndeg, fcbalphaR, fcbndegR);
 //RooCrystalBall firsErg("firsErg", "Single escape", Es, fstesc, fullWidth, fcbalpha, fcbndeg, fcbalphaR, fcbndegR);
 //RooCrystalBall secdErg("secdErg", "Double escape", Es, scdesc, fullWidth, fcbalpha, fcbndeg, fcbalphaR, fcbndegR);
 //Gaussian w same widths
 //RooGaussian fullErg("fullErg", "Full peak", Es, Egamma, fullWidth);
 //RooGaussian firsErg("firsErg", "Single escape", Es, fstesc, fullWidth);
 //RooGaussian secdErg("secdErg", "Double escape", Es, scdesc, fullWidth);
  /*RooRealVar fullWidth_ADC("fullWidth_ADC", "Full width [ADC]",initWidth,2,20);
  RooFormulaVar fstesc_ADC("fstesc_ADC", "fullPeak - (0.511/m)",RooArgSet(fullPeak,m));
 RooFormulaVar scdesc_ADC("scdesc_ADC","fullPeak - (2*0.511/m)",RooArgSet(fullPeak,m));
 RooCBShape fullErg("fullErg", "Full peak", crysADC, fullPeak, fullWidth_ADC, fcbalpha, fcbndeg);
 RooCBShape firsErg("firsErg", "Single escape", crysADC,fstesc_ADC, fullWidth_ADC, fcbalpha, fcbndeg);
 RooCBShape secdErg("secdErg", "Double escape", crysADC, scdesc_ADC, fullWidth_ADC, fcbalpha, fcbndeg);*/
 RooRealVar evtsFull("evtsFull", "Full peak yield", initEvtFull, EvtFullLow, EvtFullHigh);
 RooRealVar evtsFrst("evtsFrst", "First escape yield", initEvtFst, EvtFstLow, EvtFstHigh);
 RooRealVar evtsScnd("evtsScnd", "Second escape yield", initEvtScd,EvtScdLow,EvtScdHigh);
 RooRealVar evtsbkg("evtsbkg", "Background yield", initEvtBkg, 100, integral_evts);


 /*// Background (logistic)
 RooRealVar comCnst("comCnst", "Background const", 4.0);
 RooRealVar combeta("combeta", "Background beta", initBeta, 1,200);//0.001, 200);
 RooRealVar bkg_fixed("bkg_fixed", "Background offset", 0.0);//, -10.0, 10.0);
 //logistic background pdf
 //RooGenericPdf comPdf("comPdf", "logistic background",
 //  "pow(1.0 + exp((Es - comCnst)/combeta), -1.0) + bkg_fixed",
 //  RooArgSet(Es, comCnst, combeta,bkg_fixed));
 RooGenericPdf comPdf("comPdf", "logistic background",
   "pow(1.0 + exp((crysADC - comCnst)/combeta), -1.0) + bkg_fixed",
   RooArgSet(crysADC, comCnst, combeta,bkg_fixed));*/
   //exponential background
  // Lower bound -0.06: gives ~100:1 ratio from ADC=40 to ADC=115, already steep for Compton continuum.
  RooRealVar combeta("combeta", "Background beta", initBeta, -1.0, 1.0);
  RooExponential comPdf("comPdf", "exponential background", crysADC, combeta);
 RooAddPdf fitFun("fitFun", "Total Mu2e Calo Model",
       RooArgList(fullErg, firsErg, secdErg,comPdf ),
       RooArgList(evtsFull, evtsFrst, evtsScnd,evtsbkg ) );
       // Add them to the final list
 /*RooAddPdf fitFun("fitFun", "Total Mu2e Calo Model",
       RooArgList(fullErg, firsErg, secdErg ),
       RooArgList(evtsFull, evtsFrst, evtsScnd ) );*/




/*//adding an erfc function


//erfc centered at edge, with width sigmaADC
RooFormulaVar sigmaADC("sigmaADC","fullWidth/m",RooArgSet(fullWidth,m));
RooRealVar CE1_MeV("CE1_MeV", "Theory Edge 1",5.8842);
RooRealVar CE2_MeV("CE2_MeV","Theory Edge 2",5.3741);
RooRealVar CE3_MeV("CE3_MeV","Theory Edge 3", 4.8640);
RooFormulaVar edgeADC1("edgeADC1","CE1_MeV/m",RooArgSet(CE1_MeV,m));
RooFormulaVar edgeADC2("edgeADC2","CE2_MeV/m",RooArgSet(CE2_MeV,m));
RooFormulaVar edgeADC3("edgeADC3","CE3_MeV/m",RooArgSet(CE3_MeV,m));


RooRealVar evtsCompton1("evtsCompton1", "Yield Compton Full", initEvtCompton1, 0, integral_evts);
RooRealVar evtsCompton2("evtsCompton2", "Yield Compton 1st",  initEvtCompton2,  0, integral_evts);
RooRealVar evtsCompton3("evtsCompton3", "Yield Compton 2nd",  initEvtCompton3,  0, integral_evts);
//Mapping: @0 crystADC, @1 edgeADC,@2 sigmaADC
//Define the smearing (replaces the gaussian kernel)
//TString smearing = "0.5 * TMath::Erfc((@0 - @1) / (1.4142 * @2))";
//defining the smear
//TString rise = "1.0 / (pow(@1 - @0, 2) + 0.1)";
//combine two components
//String erfcFormula = "(1.0/(pow(@0-@1,2)+0.1)) + (0.5*TMath::Erfc((@0-@1)/(TMath::Sqrt(2)*@2))";
//TString fullformula = rise + " * " + smearing;


TString fullformula = "0.5 * TMath::Erfc((@0-@1)/(TMath::Sqrt(2)*@2))";
// RooArgList preserves insertion order so @0=crysADC, @1=edgeADC, @2=sigmaADC
RooGenericPdf compton1("compton1", "compton edge 1", fullformula, RooArgList(crysADC,edgeADC1,sigmaADC));
RooGenericPdf compton2("compton2", "compton edge 2", fullformula, RooArgList(crysADC,edgeADC2,sigmaADC));
RooGenericPdf compton3("compton3", "compton edge 3", fullformula, RooArgList(crysADC,edgeADC3,sigmaADC));


// Add them to the final list
 RooAddPdf fitFun("fitFun", "Total Mu2e Calo Model",
       RooArgList(fullErg, firsErg, secdErg,  compton1, compton2, compton3),
       RooArgList(evtsFull, evtsFrst, evtsScnd, evtsCompton1, evtsCompton2, evtsCompton3) );*/
        
 fitFun.fixCoefNormalization(RooArgSet(crysADC));
 RooPlot *chFrame = crysADC.frame(Title(title));
 h_spec->Sumw2();
 RooDataHist chSpec("crysADC","crysADC", crysADC, h_spec);
    asymSuccess = false;
    int migrad_status = -1;
auto run_one_fit = [&]() -> bool {
   fullPeak.setVal(currentPeak);
   //fullWidth.setVal(currentWidth);//fullWidth_ADC.setVal(currentWidth);
   //fcbalpha.setVal(initAlpha); 
   fcbndeg.setVal(currentn);
   fcbalphaR.setVal(initAlphaR); 
   //fcbndegR.setVal(currentnR);
   combeta.setVal(currentBeta);
   evtsFull.setVal(currentEvtFull);
   evtsFrst.setVal(currentEvtFst);
   evtsScnd.setVal(currentEvtScd);
   evtsbkg.setVal(currentEvtBkg);
   //evtsCompton1.setVal(currentEvtCompton1);
   //evtsCompton2.setVal(currentEvtCompton2);
   //evtsCompton3.setVal(currentEvtCompton3);
   RooFitResult *fitRes = nullptr;


   if (opt == "chi2") {
       // NOTE: DataError::Expected uses (O-E)^2/E; if E->0 but O>0 (e.g. near ADC=40),
       // chi2->inf. Switch to SumW2 to use (O-E)^2/O instead, which stays finite.
       //   RooFit::DataError(RooAbsData::SumW2),
      /*RooAbsReal* chi2Func = fitFun.createChi2(chSpec,
                              RooFit::DataError(RooAbsData::Expected),
                              RooFit::Range(40, 115.2),RooFit::Extended(true));*/
       RooAbsReal* chi2Func = fitFun.createChi2(chSpec,
                              RooFit::DataError(RooAbsData::Expected),
                              RooFit::Range(40, 115.2),RooFit::Extended(true));//RooFit::Range(40, 115.2),RooFit::Extended(true)


       RooMinimizer m(*chi2Func);
       m.setMinimizerType("Minuit2");
       m.setPrintLevel(-1);
       m.setStrategy(1);
       m.simplex();
       m.migrad();
      
       fitRes = m.save();
       if (fitRes->status() != 0) {
           m.setStrategy(2);
           m.migrad();
           delete fitRes;
           fitRes = m.save();
       }


       m.hesse();
       RooArgSet* fitParams = chi2Func->getVariables();
       RooRealVar* minosTarget = dynamic_cast<RooRealVar*>(fitParams->find("fullPeak"));
       if (minosTarget) {
           m.minos(RooArgSet(*minosTarget));
       }
       delete fitRes;
       fitRes = m.save();
       //=-------------------
       //delete once converging
       if (convergencestatus != 0) {
         std::cout << "[DEBUG] Fit Status History for Crystal " << crystalNo << ":" << std::endl;
           for (UInt_t i = 0; i < fitRes->numStatusHistory(); i++) {
               std::cout << "    Step " << i << ": "
                 << fitRes->statusLabelHistory(i) << " = "
                 << fitRes->statusCodeHistory(i) << std::endl;
   }
}
       //--------------------
      
       for (UInt_t i = 0; i < fitRes->numStatusHistory(); i++) {
           if (TString(fitRes->statusLabelHistory(i)).Contains("MIGRAD")) {
               migrad_status = fitRes->statusCodeHistory(i);
       }
   }
       unreducedchi2 = chi2Func->getVal();
       // chSpec.numEntries() returns ALL bins in [40,120], not just the fit range [40,115.2].
       // Count only bins whose centers fall inside the fit range.
       //int nBins = h_spec->FindBin(115.2) - h_spec->FindBin(40.0) + 1;
       int nBins = chSpec.numEntries();
       nPars = fitRes->floatParsFinal().getSize();
       ndof = nBins - nPars;
       reducedchi2= unreducedchi2/ndof;
       std::cout << "chi2/ndof for fit" << reducedchi2<< std::endl;
       if (chi2Func) delete chi2Func;
   }


   else if (opt == "nll") {
       RooAbsReal *nll = fitFun.createNLL(chSpec, Range(40,115.2));//Range(40,115.2)
       RooMinimizer m(*nll);
       m.migrad();
       m.hesse();
       fitRes = m.save();
       fval   = nll->getVal();
   }


   if (!fitRes) return false;


   convergencestatus = fitRes->status();
   RooRealVar* peak = dynamic_cast<RooRealVar*>(fitRes->floatParsFinal().find("fullPeak"));
   //coment after test data test is done
   //RooRealVar* widthpar = dynamic_cast<RooRealVar*>(fitRes->floatParsFinal().find("fullWidth"));


   if (peak) {
       fpeak = peak->getVal();
       //delete once converging-------
           std::cout << "Crystal Peak - hasAsym: " << peak->hasAsymError()
                     << " | ErrHi: " << peak->getErrorHi()
                     << " | ErrLo: " << peak->getErrorLo()
                     << " | SymErr: " << peak->getError() << std::endl;
       //-------------------


       if (peak->hasAsymError() && std::abs(peak->getErrorLo()) > 0.01 && std::abs(peak->getErrorHi())>0.01) {
           peakerrorhigh = peak->getErrorHi();
           peakerrorlo   = -peak->getErrorLo();
           asymSuccess   = true;
       } else {
           peakerrorhigh = peak->getError();
           peakerrorlo   = peak->getError();
           asymSuccess   = false;
       }


       fullPeak.setVal(fpeak);


   }
   // uncommend when only Apar is used
    /*RooRealVar* Apar = dynamic_cast<RooRealVar*>(fitRes->floatParsFinal().find("A"));
    if (Apar) { fsigma = fullWidth.getVal(); widtherrorhigh = Apar->getError()*std::sqrt(Egamma.getVal()); Aparam = Apar->getVal(); Aperr = Apar->getError();
    Bparam=0;
    Berr = 0;
    Cparam = 0;
    Cerr = 0;

    }*/
   //uncomment only when Bpar is used
  /* RooRealVar* Bpar = dynamic_cast<RooRealVar*>(fitRes->floatParsFinal().find("B"));
   if (Bpar) { Bparam = Bpar->getVal(); Berr = Bpar->getError();
        Aparam=0.0;
       Aperr=0.0;  }
   RooRealVar* Cpar = dynamic_cast<RooRealVar*>(fitRes->floatParsFinal().find("C"));
   if (Cpar) { Cparam = Cpar->getVal(); Cerr = Cpar->getError(); }
   else       { Cparam = 0.0;           Cerr  = 0.0; }*/
   //uncomment only when all Apar,Cpar and Bpar are used
   RooRealVar* Apar = dynamic_cast<RooRealVar*>(fitRes->floatParsFinal().find("A"));
   //RooRealVar* Bpar = dynamic_cast<RooRealVar*>(fitRes->floatParsFinal().find("B"));
   //RooRealVar* Cpar = dynamic_cast<RooRealVar*>(fitRes->floatParsFinal().find("C"));
   if (Apar) { fsigma = fullWidth.getVal(); widtherrorhigh = Apar->getError()*std::sqrt(Egamma.getVal()); Aparam = Apar->getVal(); Aperr = Apar->getError();}
   //if (Bpar) { fsigma = fullWidth.getVal(); widtherrorhigh = Bpar->getError()*Egamma.getVal(); Bparam = Bpar->getVal(); Berr = Bpar->getError(); }//fsigma = fullWidth.getVal(); widtherrorhigh = Bpar->getError()*Egamma.getVal();
   //if (Cpar) { Cparam = Cpar->getVal(); Cerr = Cpar->getError(); }
   /*if (widthpar) {
       fsigma         = widthpar->getVal();
       widtherrorhigh = widthpar->getError();
       widtherrorlo   = widthpar->getError();
       Aparam=0.0;
       Aperr=0.0;
   }*/
   


   RooRealVar* evtfull = dynamic_cast<RooRealVar*>(fitRes->floatParsFinal().find("evtsFull"));
   if (evtfull) {
       evtfullerrorhigh= evtfull->getError();
       evtfullerrorlo  = 0;
   }


   return (convergencestatus == 0||migrad_status == 0);
};
run_one_fit();
bool is_perfect = (convergencestatus == 0 && reducedchi2 <= 1.6 && asymSuccess);
if (convergencestatus == 0 && reducedchi2 <= 1.6) {
covarAcc.count++;
covarAcc.sumPeak     += fpeak;
//covarAcc.sumWidth    += fsigma;
covarAcc.sumAlpha    += fcbalphaparam;//comment for double sided
covarAcc.sumn    += fcbndegparam;
covarAcc.sumAlphaR    += fcbalphaRparam;
//covarAcc.sumnR    += fcbndegRparam;
covarAcc.sumBeta     += combetaparam;
covarAcc.sumEvtFull  += fr_fullparam*nEvents;
covarAcc.sumEvtFst   += fr_frstparam*nEvents;
covarAcc.sumEvtScd   += fr_scndparam*nEvents;
covarAcc.sumEvtBkg   += fr_bkgparam*nEvents;
//covarAcc.sumEvtCompton1   += fr_comptonparam1*nEvents;
//covarAcc.sumEvtCompton2   += fr_comptonparam2*nEvents;
//covarAcc.sumEvtCompton3   += fr_comptonparam3*nEvents;
}


//bool second_ok = true;
if (!is_perfect){
   std::cout << "======================================================================\n";
   std::cout << "needs refit (crystal " << crystalNo << ")\n";
   std::cout << "======================================================================\n";


   SourceFitter::nSecondFits++;
   SourceFitter::crystalsSecondFit.push_back(crystalNo);


   double newPeakGuess     = covarAcc.sumPeak     / covarAcc.count;
   //double newSigmaGuess    = covarAcc.sumWidth    / covarAcc.count;
   //double newAlphaGuess    = covarAcc.sumAlpha    / covarAcc.count;
   double newnGuess    = covarAcc.sumn    / covarAcc.count;
   double newAlphaRGuess    = covarAcc.sumAlphaR    / covarAcc.count;//uncomment for double sided
   //double newnRGuess    = covarAcc.sumnR    / covarAcc.count;//comment for double sided
   double newCombetaGuess  = covarAcc.sumBeta     / covarAcc.count;
   double newevtsFullGuess = covarAcc.sumEvtFull  / covarAcc.count;
   double newevtsFstGuess  = covarAcc.sumEvtFst   / covarAcc.count;
   double newevtsScdGuess  = covarAcc.sumEvtScd   / covarAcc.count;
   double newevtsBkgGuess  = covarAcc.sumEvtBkg   / covarAcc.count;
   //double newevtsComptonGuess1  = covarAcc.sumEvtCompton1   /covarAcc.count;
   //double newevtsComptonGuess2  = covarAcc.sumEvtCompton2   /covarAcc.count;
   //double newevtsComptonGuess3  = covarAcc.sumEvtCompton3   /covarAcc.count;


   currentPeak   = newPeakGuess;
   //currentWidth  = newSigmaGuess;
   //currentAlpha  = newAlphaGuess;
   currentn  = newnGuess;
   currentAlphaR = newAlphaRGuess;
   //currentnR = newnRGuess;
   currentBeta   = newCombetaGuess;
   currentEvtFull = newevtsFullGuess;
   currentEvtFst  = newevtsFstGuess;
   currentEvtScd  = newevtsScdGuess;
   currentEvtBkg  = newevtsBkgGuess;
   //currentEvtCompton1  = newevtsComptonGuess1;
   //currentEvtCompton2  = newevtsComptonGuess2;
   //currentEvtCompton3  = newevtsComptonGuess3;


   run_one_fit();


   if(convergencestatus == 0 && reducedchi2 <= 1.6 && asymSuccess){
       SourceFitter::nSecondFitConverged++;
       SourceFitter::crystalsSecondFitConverged.push_back(crystalNo);
       std::cout << "[SourceFitter] Second fit SUCCESS for crystal " << crystalNo << "\n";
   }
       else { 
       SourceFitter::nThirdFits++;
       SourceFitter::crystalsThirdFit.push_back(crystalNo);


       auto reset_to_defaults = [&]() {
           currentPeak    = initPeak;
           //currentWidth   = initWidth;
           //currentAlpha   = initAlpha;
           currentn   = initn;
           currentAlphaR   = initAlphaR;
           //currentnR   = initnR
           currentBeta    = initBeta;
           currentEvtFull = initEvtFull;
           currentEvtFst  = initEvtFst;
           currentEvtScd  = initEvtScd;      
           currentEvtBkg  = initEvtBkg;
           //currentEvtCompton1  = initEvtCompton1;
           //currentEvtCompton2  = initEvtCompton2;
           //currentEvtCompton3  = initEvtCompton3;
       };
       reset_to_defaults();


       auto randomize_all_parameters = [&]() {
           auto randomDouble = [](double min, double max) {
               return min + (max - min) * ((double)rand() / RAND_MAX);
           };
           currentPeak    = randomDouble(91.0, 108.0);
           //currentWidth   = randomDouble(0.2, 1.5);
           //currentAlpha   = randomDouble(0.6, 1.7);
           currentn   = randomDouble(5.0, 10);
           currentAlphaR = randomDouble(0.1, 3.5);//uncoment for double sided
           //currentnR = randomDouble(5.0, 200);
           currentBeta    = randomDouble(-1.0, -0.001);
           currentEvtFull = randomDouble(0.05 * integral_evts, integral_evts);
           currentEvtFst  = randomDouble(0.2 * integral_evts, integral_evts);
           //currentEvtFst  = randomDouble(0.2 * integral_evts, 0.5*integral_evts);
           currentEvtScd  = randomDouble(0.05 * integral_evts, integral_evts);
           currentEvtBkg  = randomDouble(0.05 * integral_evts, integral_evts);
           //currentEvtCompton1  = randomDouble(0.01 * integral_evts, 0.1*integral_evts);
           //currentEvtCompton2  = randomDouble(0.01 * integral_evts, 0.1*integral_evts);
           //currentEvtCompton3  = randomDouble(0.01 * integral_evts, 0.1*integral_evts);
       };


       struct FitCandidate {
           double chi2;
           std::map<std::string, double> p;
       };
       auto retry_third_fit = [&](int cId) -> bool {
           std::vector<FitCandidate> leaderboard;
           const int MAX_THIRD_TRIES = 40;
           int tries = 0;
            migrad_status = -1;


           for ( tries = 0; tries <= MAX_THIRD_TRIES; ++tries) {
               randomize_all_parameters();


               if (run_one_fit()) {
                   if (convergencestatus == 0 &&reducedchi2 <= 6 && asymSuccess) {
                       SourceFitter::thirdFitRetryCount[cId] = tries;
                       return true;
                   }


                   if (migrad_status == 0) {
                       FitCandidate cand;
                       cand.chi2 = reducedchi2;
                       cand.p["pk"] = currentPeak;  
                       //cand.p["wd"] = currentWidth;
                       //cand.p["al"] = currentAlpha;
                       cand.p["n"] = currentn;
                       cand.p["alR"] = currentAlphaR;//uncomment for double sided
                       //cand.p["nR"] = currentnR;
                       cand.p["bt"] = currentBeta;
                       cand.p["ef"] = currentEvtFull;
                       cand.p["e1"] = currentEvtFst;
                       cand.p["e2"] = currentEvtScd; 
                       cand.p["eb"] = currentEvtBkg;
                       //cand.p["ec1"] = currentEvtCompton1;
                       //cand.p["ec2"] = currentEvtCompton2;
                       //cand.p["ec3"] = currentEvtCompton3;
                       leaderboard.push_back(cand);
                   }
               }
           }
           SourceFitter::thirdFitRetryCount[cId] = MAX_THIRD_TRIES;


           std::cout << "[DEBUG] Crystal " << cId << " exhausted 40 tries. Leaderboard has "
                     << leaderboard.size() << " valid MIGRAD candidates to choose from.\n";


           if (leaderboard.empty()) return false;


           std::sort(leaderboard.begin(), leaderboard.end(), [](const FitCandidate& a, const FitCandidate& b) {
               return a.chi2 < b.chi2;
           });


           int to_try = std::min((int)leaderboard.size(), 10);
           std::cout << "[DEBUG] Crystal " << cId << " starting fallback: testing " << to_try << " candidates." << std::endl;


           double bestFallbackChi2 = std::numeric_limits<double>::max();
           int bestCandIdx = -1;


           for (int i = 0; i < to_try; i++) {
               auto& cand = leaderboard[i];
               currentPeak     = cand.p["pk"];
               //currentWidth    = cand.p["wd"];
               currentAlpha    = cand.p["al"];
               currentn        = cand.p["n"];
               //currentAlphaR   = cand.p["alR"];//uncomment for double sided
               //currentnR        = cand.p["nR"];
               currentBeta     = cand.p["bt"];
               currentEvtFull  = cand.p["ef"];
               currentEvtFst   = cand.p["e1"];
               currentEvtScd   = cand.p["e2"];
               currentEvtBkg   = cand.p["eb"];
               //cand.p["ec1"] = currentEvtCompton1;
               //cand.p["ec2"] = currentEvtCompton2;
               //cand.p["ec3"] = currentEvtCompton3;


               fullPeak.removeAsymError();


               bool fallback_ok = run_one_fit();
               if (fallback_ok) {
                   if (reducedchi2 <= 1.6) {
                       std::cout << "[DEBUG] Crystal " << cId << " converged on fallback candidate #"
                                 << i << " (chi2=" << reducedchi2 << ")" << std::endl;
                       convergencestatus = 0;
                       migrad_status = 0;
                       asymSuccess = false;
                       return true;
                   }
                   if (reducedchi2 < bestFallbackChi2) {
                       bestFallbackChi2 = reducedchi2;
                       bestCandIdx = i;
                   }
                   std::cout << "[DEBUG] Fallback #" << i << " succeeded but chi2 ("
                             << reducedchi2 << ") is still > 1.6. Trying next candidate...\n";
               }
           }


           // No candidate met <= 1.6 — restore and re-run the best found
           if (bestCandIdx >= 0) {
               auto& best = leaderboard[bestCandIdx];
               currentPeak     = best.p["pk"];
               //currentWidth    = best.p["wd"];
               //currentAlpha    = best.p["al"];
               currentn        = best.p["n"];
               currentAlphaR   = best.p["alR"];//uncomment for double sided
               //currentnR        = best.p["nR"];
               currentBeta     = best.p["bt"];
               currentEvtFull  = best.p["ef"];
               currentEvtFst   = best.p["e1"];
               currentEvtScd   = best.p["e2"];
               currentEvtBkg   = best.p["eb"];
               //currentEvtCompton1 = best.p["ec1"];
               //currentEvtCompton2 = best.p["ec2"];
               //currentEvtCompton3 = best.p["ec3"];
               fullPeak.removeAsymError();
               bool final_ok = run_one_fit();
               if (final_ok) {
                   std::cout << "[DEBUG] Crystal " << cId << " using best available fallback chi2="
                             << reducedchi2 << " (threshold not met, best was candidate #" << bestCandIdx << ")\n";
                   convergencestatus = 0;
                   migrad_status = 0;
                   asymSuccess = false;
                   return true;
               }
           }
           return false;
       }; // End of lambda
      
       bool third_ok = retry_third_fit(crystalNo);
       if (third_ok) {
           SourceFitter::nThirdFitConverged++;
           SourceFitter::crystalsThirdFitConverged.push_back(crystalNo);
           std::cout << "[SourceFitter] Third fit SUCCESS for crystal " << crystalNo << "\n";
       }
       else {
           std::cout << "======================================================================\n";
           std::cout << "[SourceFitter] Random-start refit FAILED for crystal "
                     << crystalNo << " after "
                     << SourceFitter::thirdFitRetryCount[crystalNo]
                     << " attempts.\n";
           std::cout << "======================================================================\n";
       }
   }


}


else {
   std::cout << "[SourceFitter] Fit ok for crystal " << crystalNo
             << " fpeak=" << fpeak << "\n";
   nFirstFitConverged++;
}


if (migrad_status>0){
   SourceFitter::convFailures.push_back({crystalNo, convergencestatus});
}
else if (asymSuccess && convergencestatus == 0) {
   SourceFitter::nAsymErrors++;
   SourceFitter::crystalsWithAsymErrors.push_back(crystalNo);
}
else if (!asymSuccess && convergencestatus == 0) {
   SourceFitter::nHesseFallbacks++;
   SourceFitter::crystalsHesseFallback.push_back(crystalNo);
}


TF1* totalModelFunc = fitFun.asTF(RooArgList(crysADC));
redpeak = totalModelFunc->GetMaximumX(40.0, 115.2);


 fstpeak = fstEsPeak.getVal();
 scdpeak = scdEsPeak.getVal();
 fcbalphaparam  = fcbalpha.getVal();
 fcbndegparam   = fcbndeg.getVal();
 fcbalphaRparam = fcbalphaR.getVal();
 fcbndegRparam  = fcbndegR.getVal();
 comCnstparam =0;// comCnst.getVal();
 combetaparam = combeta.getVal();
 fr_fullparam = evtsFull.getVal()/integral_evts;
 fr_frstparam = evtsFrst.getVal()/integral_evts;
 fr_scndparam = evtsScnd.getVal()/integral_evts;
 fr_bkgparam = evtsbkg.getVal()/integral_evts;
 fr_comptonparam1 =0;//evtsCompton1.getVal()/integral_evts;
 fr_comptonparam2 =0;//evtsCompton2.getVal()/integral_evts;
 fr_comptonparam3 =0;//evtsCompton3.getVal()/integral_evts;
 crystalNoparam = crystalNo;       
 etaparam = eta.getVal();
 errbarhigh = mparam*(peakerrorhigh/fpeak);
 errbarlo = mparam*(peakerrorlo/fpeak);
 Esparam = Es.getVal();
 //pval = TMath::Prob(chiSq*11, 11);
 if (ndof > 0) {
   pval = TMath::Prob(unreducedchi2, ndof); 
} else {
   pval = -1.0;  // flag for error
}


  chSpec.plotOn(chFrame, MarkerColor(kBlack), LineColor(kBlack), MarkerSize(0.5), Name("chSpec"));
 fitFun.plotOn(chFrame, LineColor(kRed), LineStyle(1), Name("fit"));
 fitFun.plotOn(chFrame, Components(fullErg), LineColor(kOrange), LineStyle(5), Name("main"));
 fitFun.plotOn(chFrame, Components(firsErg), LineColor(kViolet), LineStyle(5), Name("fescape"));
 fitFun.plotOn(chFrame, Components(secdErg), LineColor(kCyan), LineStyle(5), Name("sescape"));
 fitFun.plotOn(chFrame, Components(comPdf), LineColor(kBlue), LineStyle(5), Name("background"));
 //fitFun.plotOn(chFrame, Components(compton1), LineColor(kGreen), LineStyle(5), Name("compton1"));
 //fitFun.plotOn(chFrame, Components(compton2), LineColor(kGreen+2), LineStyle(5), Name("compton2"));
 //fitFun.plotOn(chFrame, Components(compton3), LineColor(kGreen+4), LineStyle(5), Name("compton3"));
 chiSq = chFrame->chiSquare("fit", "chSpec", nPars);
 std::cout << "[CHI2 COMPARE] SiPM " << crystalNo
           << "  reducedchi2 (minimizer, ndof=" << ndof << ") = " << reducedchi2
           << "  chiSq (chFrame, npars from 40-120) = " << chiSq << std::endl;
 if (reducedchi2 > 1.6){
   SourceFitter::badchi2++;
   SourceFitter::crystalswithbadchi2.push_back(crystalNo);
   SourceFitter::badchi2Values[crystalNo] = reducedchi2;
 }
 mparam = m.getVal();
 //draw vertical lines at the theorectical peak locations:
 double xLines[] = {98.08, 89.904, 81.728};
 int colors[] = {kOrange-2, kViolet-9, kCyan-9};
 for (int i = 0; i < 3; i++) {
   double yMax = chFrame->GetMaximum() > 0 ? chFrame->GetMaximum() : 10000;
   TLine *line = new TLine(xLines[i], 0, xLines[i], 0.5*yMax);
   line->SetLineColor(colors[i]);
   line->SetLineStyle(2);      
   line->SetLineWidth(2);
   chFrame->addObject(line);
}

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
 TPaveLabel *pchi2 = new TPaveLabel(0.15, 0.75, 0.25, 0.65, Form("#chi^{2}/ndf = %4.2f", reducedchi2), "brNDC");
 pchi2 -> SetFillStyle(0);
 pchi2 -> SetBorderSize(0);
 pchi2 -> SetTextSize(0.4);
 pchi2 -> SetTextFont(42);
 pchi2 -> SetTextColor(kBlack);
 pchi2 -> SetFillColor(kWhite);
 chFrame -> addObject(pchi2);
 TPaveLabel *fpk = new TPaveLabel(0.15, 0.65, 0.25, 0.55, Form("#mu_{main} = %.2f^{+%.2f}_{-%.2f}", fpeak, peakerrorhigh,peakerrorlo), "brNDC");
 fpk -> SetFillStyle(0);
 fpk -> SetBorderSize(0);
 fpk -> SetTextSize(0.4);
 fpk -> SetTextFont(42);
 fpk -> SetTextColor(kBlack);
 fpk -> SetFillColor(kWhite);
 chFrame -> addObject(fpk);
 TPaveLabel *fsg = new TPaveLabel(0.15, 0.55, 0.25, 0.45, Form("#sigma_{main} =%.2f #pm %.2f", fsigma,widtherrorhigh), "brNDC");
 fsg -> SetFillStyle(0);
 fsg -> SetBorderSize(0);
 fsg -> SetTextSize(0.4);
 fsg -> SetTextFont(42);
 fsg -> SetTextColor(kBlack);
 fsg -> SetFillColor(kWhite);
 chFrame -> addObject(fsg);


   TPad *pad1 = new TPad("pad1", "Top pad", 0, 0.25, 1, 1.0);
   pad1->SetBottomMargin(0.035);
   pad1->Draw();
   pad1->cd(); 
 chFrame -> SetYTitle("Event count");
 chFrame -> Draw();
 TLegend* legend = new TLegend(0.5, 0.7);
 legend->SetBorderSize(0);
 legend->SetFillStyle(0);
 legend->AddEntry("main", "main (eqv. 6.13MeV)", "L");
 legend->AddEntry("fescape", "first escape", "L");
 legend->AddEntry("sescape", "second escape", "L");
 legend->AddEntry("background", "background", "L");
 //legend->AddEntry("compton1", "compton1", "L");
 //legend->AddEntry("compton2", "compton2", "L");
 //legend->AddEntry("compton3", "compton3", "L");
 legend->Draw();
   can->cd();
   TPad *pad2 = new TPad("pad2", "Bottom pad", 0, 0.0, 1, 0.25);
   pad2->SetTopMargin(0.07);
   pad2->SetBottomMargin(0.3);
   pad2->Draw();
   pad2->cd();
  
// residual histogram
double xMin = 40.0;
double xMax = 115.0;//h_spec->GetXaxis()->GetXmax();
int nBins   = h_spec->FindBin(xMax) - h_spec->FindBin(xMin) + 1;
int startBin = h_spec->FindBin(xMin);
int endBin = h_spec->FindBin(xMax);
double totalYield = h_spec->Integral(startBin,endBin);
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
   hresidual->SetBinContent(i - startBin + 1, res);
}
hresidual->SetStats(0);
hresidual->SetTitle("");
hresidual->GetYaxis()->SetTitle("#splitline{Normalized Residuals}{(data - fit)/#sigma}");
hresidual->GetYaxis()->CenterTitle(true);
hresidual->GetYaxis()->SetTitleSize(0.08);
hresidual->GetYaxis()->SetLabelSize(0.10);
hresidual->GetYaxis()->SetTitleOffset(0.45);
hresidual->GetYaxis()->SetNdivisions(505);  
hresidual->GetXaxis()->SetTitleSize(0.12);
hresidual->GetXaxis()->SetLabelSize(0.10);
hresidual->GetXaxis()->SetTitle("ADC [counts]");
hresidual->Draw("HIST");


 can -> SaveAs(oName);
 can->Close(); 
 delete can;   
 can = nullptr;
 covar->Fill();


// =========================================================
// CONTOUR CONFIGURATION "SWITCH BOARD"
// =========================================================
if (contour) {
   std::cout << "[DEBUG] Initializing 2D Contour: " << xVar << " vs " << yVar << std::endl;


   CaloSourceCalib::MakeContourPlot(
       fitFun, chSpec, opt, crystalNo,
       xVar, yVar,               
       fullPeak,  fpeak,  peakerrorlo,  peakerrorhigh,
       //fullWidth, fsigma, widtherrorlo, widtherrorhigh,
       //fcbndeg,
       evtsFull,
       evtsFrst,
       evtsScnd
       //evtsbkg,
       //comCnst,
       //combeta
       //compton stuff added here
   );
}


}






