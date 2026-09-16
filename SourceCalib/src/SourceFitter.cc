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
bool SourceFitter::singlesided = false;
struct CovarAccumulator {
   int count = 0;
   double sumPeak = 0.0;
   double sumAlpha = 0.0;
   double sumn = 0.0;
   double sumAlphaR = 0.0; 
   double sumBeta = 0.0;
   double sumEvtFull = 0.0;
   double sumEvtFst = 0.0;
   double sumEvtScd = 0.0;
   double sumEvtBkg = 0.0;
} covarAcc;
// ---- SILENCE ROOFIT & MINUIT2 ----
static bool suppress_messages = [](){

   RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

   ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1); 
   ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(0); 
   gErrorIgnoreLevel = kError;  


   return true;
}();
void SourceFitter::FitCrystal(TH1F* h_spec, TString opt, int crystalNo,  TTree *covar, Int_t &nEvents,Int_t &convergencestatus, Float_t &fpeak, Float_t &peakerrorhigh,Float_t &peakerrorlo,Float_t &redpeak, Float_t &fsigma,Float_t &widtherrorhigh,Float_t &widtherrorlo, Float_t &chiSq, Float_t &fstpeak, Float_t &scdpeak,Float_t &fcbalphaparam,Float_t &fcbndegparam,Float_t &comCnstparam, Float_t &combetaparam, Float_t &fr_fullparam, Float_t &fr_frstparam,Float_t &fr_scndparam,Float_t &crystalNoparam,Float_t &fr_bkgparam,Float_t &fr_comptonparam1,Float_t &fr_comptonparam2,Float_t &fr_comptonparam3, Float_t &pval,Float_t &h_means,Float_t &h_stddevs, Float_t &unreducedchi2,Float_t &fval,Float_t &mparam,Float_t &etaparam,Int_t &ndof,bool contour,TString xVar, TString yVar, Float_t &errbarhigh, Float_t &errbarlo,Float_t &evtfullerrorhigh,Float_t &evtfullerrorlo,Float_t &Esparam,Float_t &Aparam,Float_t &Aperr,Float_t &Bparam,Float_t &Berr,Float_t &Cparam,Float_t &Cerr,Float_t &fcbalphaRparam, Float_t &fcbndegRparam, bool islyso ){

std::cout << "**************************************************************\n";
std::cout << "Fitting crystal " << crystalNo << ".\n";
std::cout << "**************************************************************\n";
convergencestatus = -1; 
bool asymSuccess = false;    
double currentPeak  = 0.0;
double currentAlpha = 0.0;
double currentn= 0.0;
double currentAlphaR = 0.0;
double currentBeta  = 0.0;
double currentEvtFull = 0.0;
double currentEvtFst  = 0.0;
double currentEvtScd  = 0.0;
double currentEvtBkg  = 0.0;

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
 initEvtBkg,initBeta,initAlpha,initn,initAlphaR = 0,initNdegR = 0 ;

 int low = h_spec->FindBin(40);
 int high = h_spec->FindBin(115);
 int integral_evts = h_spec->Integral(low, high);  
 nEvents = integral_evts;

 float reducedchi2 = 0.0;
 int nPars = 9;
 if (singlesided){
    if (islyso){
  // LYSO SPECIFIC PARAMETERS
  initPeak  = 95.08;
  PeakLow = 85;
  PeakHigh = 105;
  initAlpha  = 0.8;
  initn = 2;
  initBeta = -0.09;
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
 }
 else{
   // CsI SPECIFIC PARAMETERS
  initPeak  = 90;
  PeakLow = 85;
  PeakHigh = 108;
  initAlpha  = 0.5;
  initn = 1;
  initBeta = -0.09;
  initEvtFull = 9000;
  EvtFullLow= 0.05*integral_evts;
  EvtFullHigh =  0.6*integral_evts;
  initEvtFst = 9000;
  EvtFstLow= 0.2*integral_evts;
  EvtFstHigh = integral_evts;
  initEvtScd = 9000;
  EvtScdLow= 0.05*integral_evts;
  EvtScdHigh = integral_evts;
  initEvtBkg= 5000;
 }
 }
 else{
 if (islyso){
  // LYSO SPECIFIC PARAMETERS
  initPeak  = 95.08;
  PeakLow = 85;
  PeakHigh = 105;
  initAlpha  = 0.8;
  initAlphaR = 3.0;
  initNdegR  = 10;
  initn = 10; 
  initBeta = -0.09;
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
 }
 else{
  // CsI SPECIFIC PARAMETERS
  initPeak  = 90;
  PeakLow = 85;
  PeakHigh = 108;
  initAlpha  = 0.5;
  initAlphaR = 2.0;
  initNdegR  = 10;
  initn = 10;  
  initBeta = -0.09;
  initEvtFull = 9000;
  EvtFullLow= 0.05*integral_evts;
  EvtFullHigh =  0.6*integral_evts;
  initEvtFst = 9000;
  EvtFstLow= 0.2*integral_evts;
  EvtFstHigh = integral_evts;
  initEvtScd = 9000;
  EvtScdLow= 0.05*integral_evts;
  EvtScdHigh = integral_evts;
  initEvtBkg= 5000;

 }
}
if(singlesided){
  currentPeak  = initPeak;
  currentAlpha = initAlpha;
  currentn = initn;
  currentBeta  = initBeta;
  currentEvtFull = initEvtFull;
  currentEvtFst  = initEvtFst;
  currentEvtScd  = initEvtScd;
  currentEvtBkg  = initEvtBkg;
}
else{
   currentPeak  = initPeak;
  currentn = initn;
  currentAlphaR = initAlphaR;
  currentBeta  = initBeta;
  currentEvtFull = initEvtFull;
  currentEvtFst  = initEvtFst;
  currentEvtScd  = initEvtScd;
  currentEvtBkg  = initEvtBkg;
}
 RooRealVar m_e("m_e", "electron energy in MeV", 0.511);
 RooRealVar crysADC("crysADC", "ADC [counts]", 40, 115); 
 RooRealVar E0("E0", "energy offset [MeV]", 0.0); //test param

 RooRealVar fullPeak("fullPeak", "Full peak [ADC]", initPeak, PeakLow, PeakHigh);
 RooFormulaVar eta("eta", "ADC/MeV", "fullPeak/ (6.13-E0)", RooArgSet(fullPeak, E0));
 RooFormulaVar m("m","Mev/ADC", "(6.13-E0)/fullPeak", RooArgSet(fullPeak, E0));
 RooFormulaVar Es("Es", "calibrated energy [MeV]", "m*crysADC + E0", RooArgSet(m, crysADC,E0));
 RooFormulaVar fstEsPeak("fstEsPeak", "First escape", "fullPeak - m_e*eta", RooArgSet(fullPeak, m_e,eta));
 RooFormulaVar scdEsPeak("scdEsPeak", "Second escape", "fullPeak - (2*m_e)*eta", RooArgSet(fullPeak, m_e,eta));

 RooRealVar Egamma("Egamma", "Full peak [MeV]",6.13);
 RooRealVar fstesc("fstesc", "first peak [MeV]",5.619);
 RooRealVar scdesc("scdesc", "second peak [MeV]",5.108);
 // Stochastic term: sigma = A * sqrt(E), where A [sqrt(MeV)] encodes photon collection efficiency
 RooRealVar A("A", "stochastic term [sqrt(MeV)]", 0.1, 0, 1.0);
 RooFormulaVar fullWidth("fullWidth", "Full width [MeV]", "A*std::sqrt(Egamma)", RooArgSet(A, Egamma));
 RooFormulaVar fstwidth ("fstwidth",  "first width [MeV]",  "A*std::sqrt(fstesc)",  RooArgSet(A, fstesc));
 RooFormulaVar scdwidth ("scdwidth",  "second width [MeV]", "A*std::sqrt(scdesc)", RooArgSet(A, scdesc));

 RooRealVar *fcbalpha  = nullptr;
 RooRealVar *fcbndeg   = nullptr;
 RooRealVar *fcbalphaR = nullptr;
 RooRealVar *fcbndegR  = nullptr;
 RooAbsPdf  *fullErg   = nullptr;
 RooAbsPdf  *firsErg   = nullptr;
 RooAbsPdf  *secdErg   = nullptr;
if(singlesided){
 //Single sided crystal ball (RooCBShape) w different width using A/B/C
 fcbalpha = new RooRealVar("fcbalpha", "alpha", initAlpha, 0.1, 3.0);
 fcbalpha->setConstant(kFALSE);
 fcbndeg  = new RooRealVar("fcbndeg", "n", initn, 1, 4);
 fcbndeg->setConstant(kFALSE);
 fullErg = new RooCBShape("fullErg", "Full peak", Es, Egamma, fullWidth, *fcbalpha, *fcbndeg);
 firsErg = new RooCBShape("firsErg", "Single escape", Es, fstesc, fstwidth, *fcbalpha, *fcbndeg);
 secdErg = new RooCBShape("secdErg", "Double escape", Es, scdesc, scdwidth, *fcbalpha, *fcbndeg);
}
else{
 //Double sided crystal ball (RooCrystalBall) with different widths and same alpha L/R and N L/R
 fcbalpha = new RooRealVar("fcbalpha", "alpha", initAlpha, 0.1, 5.0);
 fcbalpha->setConstant(kTRUE);
 fcbndeg  = new RooRealVar("fcbndeg", "n", initn, 5, 20);
 fcbndeg->setConstant(kFALSE);
 fcbalphaR = new RooRealVar("fcbalphaR", "alpha R", initAlphaR, 0.1, 3.5);
 fcbalphaR->setConstant(kFALSE);
 fcbndegR = new RooRealVar("fcbndegR", "n R", initNdegR, 1, 15);
 fcbndegR->setConstant(kTRUE);
 fullErg = new RooCrystalBall("fullErg", "Full peak", Es, Egamma, fullWidth, *fcbalpha, *fcbndeg, *fcbalphaR, *fcbndegR);
 firsErg = new RooCrystalBall("firsErg", "Single escape", Es, fstesc, fstwidth, *fcbalpha, *fcbndeg, *fcbalphaR, *fcbndegR);
 secdErg = new RooCrystalBall("secdErg", "Double escape", Es, scdesc, scdwidth, *fcbalpha, *fcbndeg, *fcbalphaR, *fcbndegR);
}
 RooRealVar evtsFull("evtsFull", "Full peak yield", initEvtFull, EvtFullLow, EvtFullHigh);
 RooRealVar evtsFrst("evtsFrst", "First escape yield", initEvtFst, EvtFstLow, EvtFstHigh);
 RooRealVar evtsScnd("evtsScnd", "Second escape yield", initEvtScd,EvtScdLow,EvtScdHigh);
 RooRealVar evtsbkg("evtsbkg", "Background yield", initEvtBkg, 100, integral_evts);
   //exponential background
  RooRealVar combeta("combeta", "Background beta", initBeta, -1.0, 1.0);
  RooExponential comPdf("comPdf", "exponential background", crysADC, combeta);
 RooAddPdf fitFun("fitFun", "Total Mu2e Calo Model",
       RooArgList(*fullErg, *firsErg, *secdErg,comPdf ),
       RooArgList(evtsFull, evtsFrst, evtsScnd,evtsbkg ) );
        
 fitFun.fixCoefNormalization(RooArgSet(crysADC));
 RooPlot *chFrame = crysADC.frame(Title(title));
 h_spec->Sumw2();
 RooDataHist chSpec("crysADC","crysADC", crysADC, h_spec);
    asymSuccess = false;
    int migrad_status = -1;
auto run_one_fit = [&]() -> bool {
    if(singlesided){
        fcbalpha->setVal(currentAlpha);
        fcbndeg->setVal(currentn);
    }
    else{
        fcbndeg->setVal(currentn);
        fcbalphaR->setVal(currentAlphaR);
}
     fullPeak.setVal(currentPeak);
    combeta.setVal(currentBeta);
    evtsFull.setVal(currentEvtFull);
    evtsFrst.setVal(currentEvtFst);
    evtsScnd.setVal(currentEvtScd);
    evtsbkg.setVal(currentEvtBkg);
   RooFitResult *fitRes = nullptr;

   if (opt == "chi2") {
       RooAbsReal* chi2Func = fitFun.createChi2(chSpec,
                              RooFit::DataError(RooAbsData::Expected),
                              RooFit::Range(40, 115.2),RooFit::Extended(true));

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
      
       for (UInt_t i = 0; i < fitRes->numStatusHistory(); i++) {
           if (TString(fitRes->statusLabelHistory(i)).Contains("MIGRAD")) {
               migrad_status = fitRes->statusCodeHistory(i);
       }
   }
       unreducedchi2 = chi2Func->getVal();
       int nBins = chSpec.numEntries();
       nPars = fitRes->floatParsFinal().getSize();
       ndof = nBins - nPars;
       reducedchi2= unreducedchi2/ndof;
       if (chi2Func) delete chi2Func;
   }


   else if (opt == "nll") {
       RooAbsReal *nll = fitFun.createNLL(chSpec, Range(40,115.2));
       RooMinimizer m(*nll);
       m.migrad();
       m.hesse();
       fitRes = m.save();
       fval   = nll->getVal();
   }


   if (!fitRes) return false;


   convergencestatus = fitRes->status();
   RooRealVar* peak = dynamic_cast<RooRealVar*>(fitRes->floatParsFinal().find("fullPeak"));

   if (peak) {
       fpeak = peak->getVal();
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
   RooRealVar* Apar = dynamic_cast<RooRealVar*>(fitRes->floatParsFinal().find("A"));
   if (Apar) { fsigma = fullWidth.getVal(); widtherrorhigh = Apar->getError()*std::sqrt(Egamma.getVal()); Aparam = Apar->getVal(); Aperr = Apar->getError();}
   
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
    if (singlesided){
        covarAcc.sumAlpha    += fcbalphaparam;
    }
    else{
        covarAcc.sumAlphaR    += fcbalphaRparam;
    }
covarAcc.count++;
covarAcc.sumPeak     += fpeak;
covarAcc.sumn    += fcbndegparam;
covarAcc.sumBeta     += combetaparam;
covarAcc.sumEvtFull  += fr_fullparam*nEvents;
covarAcc.sumEvtFst   += fr_frstparam*nEvents;
covarAcc.sumEvtScd   += fr_scndparam*nEvents;
covarAcc.sumEvtBkg   += fr_bkgparam*nEvents;
}

if (!is_perfect){
   std::cout << "======================================================================\n";
   std::cout << "needs refit (crystal " << crystalNo << ")\n";
   std::cout << "======================================================================\n";


   SourceFitter::nSecondFits++;
   SourceFitter::crystalsSecondFit.push_back(crystalNo);
   double newPeakGuess     = covarAcc.sumPeak     / covarAcc.count;
   double newnGuess    = covarAcc.sumn    / covarAcc.count;
   double newCombetaGuess  = covarAcc.sumBeta     / covarAcc.count;
   double newevtsFullGuess = covarAcc.sumEvtFull  / covarAcc.count;
   double newevtsFstGuess  = covarAcc.sumEvtFst   / covarAcc.count;
   double newevtsScdGuess  = covarAcc.sumEvtScd   / covarAcc.count;
   double newevtsBkgGuess  = covarAcc.sumEvtBkg   / covarAcc.count;
   if (singlesided){
        currentAlpha  = covarAcc.sumAlpha  / covarAcc.count;
   }
   else{
        currentAlphaR = covarAcc.sumAlphaR / covarAcc.count;
   }


   currentPeak   = newPeakGuess;
   currentn  = newnGuess;
   currentBeta   = newCombetaGuess;
   currentEvtFull = newevtsFullGuess;
   currentEvtFst  = newevtsFstGuess;
   currentEvtScd  = newevtsScdGuess;
   currentEvtBkg  = newevtsBkgGuess;


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
            if (singlesided){
                currentAlpha   = initAlpha;
            }
            else{
                currentAlphaR   = initAlphaR;
            }
           currentPeak    = initPeak;
           currentn   = initn;
           currentBeta    = initBeta;
           currentEvtFull = initEvtFull;
           currentEvtFst  = initEvtFst;
           currentEvtScd  = initEvtScd;      
           currentEvtBkg  = initEvtBkg;
       };
       reset_to_defaults();


       auto randomize_all_parameters = [&]() {
           auto randomDouble = [](double min, double max) {
               return min + (max - min) * ((double)rand() / RAND_MAX);
           };
           if (singlesided){
            currentAlpha   = randomDouble(0.6, 1.7);
            currentn       = randomDouble(1.0, 4.0); 
           }
           else{
            currentAlphaR = randomDouble(0.1, 3.5);
            currentn      = randomDouble(5.0, 20.0); 
           }
           currentPeak    = randomDouble(91.0, 108.0);
           currentBeta    = randomDouble(-1.0, -0.001);
           currentEvtFull = randomDouble(0.05 * integral_evts, integral_evts);
           currentEvtFst  = randomDouble(0.2 * integral_evts, integral_evts);
           currentEvtScd  = randomDouble(0.05 * integral_evts, integral_evts);
           currentEvtBkg  = randomDouble(0.05 * integral_evts, integral_evts);
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
                       if (singlesided){
                        cand.p["al"] = currentAlpha;
                       }
                       else{
                        cand.p["alR"] = currentAlphaR; 
                       }
                       cand.chi2 = reducedchi2;
                       cand.p["pk"] = currentPeak;  
                       cand.p["n"] = currentn;
                       cand.p["bt"] = currentBeta;
                       cand.p["ef"] = currentEvtFull;
                       cand.p["e1"] = currentEvtFst;
                       cand.p["e2"] = currentEvtScd; 
                       cand.p["eb"] = currentEvtBkg;
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
               if (singlesided){
                currentAlpha    = cand.p["al"];
               }
               else{
                currentAlphaR   = cand.p["alR"]; 
               }
               currentPeak     = cand.p["pk"];
               currentn        = cand.p["n"];
               currentBeta     = cand.p["bt"];
               currentEvtFull  = cand.p["ef"];
               currentEvtFst   = cand.p["e1"];
               currentEvtScd   = cand.p["e2"];
               currentEvtBkg   = cand.p["eb"];
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
               if (singlesided){
                 currentAlpha    = best.p["al"];
               }
               else{
                currentAlphaR   = best.p["alR"];
               }
               currentPeak     = best.p["pk"];
               currentn        = best.p["n"];
               currentBeta     = best.p["bt"];
               currentEvtFull  = best.p["ef"];
               currentEvtFst   = best.p["e1"];
               currentEvtScd   = best.p["e2"];
               currentEvtBkg   = best.p["eb"];
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
if(singlesided){
    fcbalphaparam  = fcbalpha->getVal();
    fcbndegparam   = fcbndeg->getVal();
    fcbalphaRparam = 0;
    fcbndegRparam  = 0;
}
else{
    fcbalphaRparam = fcbalphaR->getVal();
    fcbndegRparam  = fcbndegR->getVal();
    fcbalphaparam  = fcbalpha->getVal();
    fcbndegparam   = fcbndeg->getVal();
}
 fstpeak = fstEsPeak.getVal();
 scdpeak = scdEsPeak.getVal();
 comCnstparam =0;
 combetaparam = combeta.getVal();
 fr_fullparam = evtsFull.getVal()/integral_evts;
 fr_frstparam = evtsFrst.getVal()/integral_evts;
 fr_scndparam = evtsScnd.getVal()/integral_evts;
 fr_bkgparam = evtsbkg.getVal()/integral_evts;
 crystalNoparam = crystalNo;       
 etaparam = eta.getVal();
 errbarhigh = mparam*(peakerrorhigh/fpeak);
 errbarlo = mparam*(peakerrorlo/fpeak);
 Esparam = Es.getVal();
 if (ndof > 0) {
   pval = TMath::Prob(unreducedchi2, ndof); 
} else {
   pval = -1.0; 
}


 chSpec.plotOn(chFrame, MarkerColor(kBlack), LineColor(kBlack), MarkerSize(0.5), Name("chSpec"));
 fitFun.plotOn(chFrame, LineColor(kRed), LineStyle(1), Name("fit"));
 fitFun.plotOn(chFrame, Components(*fullErg), LineColor(kOrange), LineStyle(5), Name("main"));
 fitFun.plotOn(chFrame, Components(*firsErg), LineColor(kViolet), LineStyle(5), Name("fescape"));
 fitFun.plotOn(chFrame, Components(*secdErg), LineColor(kCyan), LineStyle(5), Name("sescape"));
 fitFun.plotOn(chFrame, Components(comPdf), LineColor(kBlue), LineStyle(5), Name("background"));
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
 //delete for real data
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
 TPaveLabel *fsg = new TPaveLabel(0.15, 0.55, 0.25, 0.45, Form("#sigma_{main} =%.2f #pm %.3f", fsigma,widtherrorhigh), "brNDC");
 fsg -> SetFillStyle(0);
 fsg -> SetBorderSize(0);
 fsg -> SetTextSize(0.4);
 fsg -> SetTextFont(42);
 fsg -> SetTextColor(kBlack);
 fsg -> SetFillColor(kWhite);
 chFrame -> addObject(fsg);
 TString cbModelLabel = singlesided ? "Single-sided CB" : "Double-sided CB";
 TPaveLabel *pmodel = new TPaveLabel(0.15, 0.45, 0.28, 0.35, cbModelLabel, "brNDC");
 pmodel -> SetFillStyle(0);
 pmodel -> SetBorderSize(0);
 pmodel -> SetTextSize(0.4);
 pmodel -> SetTextFont(42);
 pmodel -> SetTextColor(kBlack);
 pmodel -> SetFillColor(kWhite);
 chFrame -> addObject(pmodel);


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
 legend->Draw();
   can->cd();
   TPad *pad2 = new TPad("pad2", "Bottom pad", 0, 0.0, 1, 0.25);
   pad2->SetTopMargin(0.07);
   pad2->SetBottomMargin(0.3);
   pad2->Draw();
   pad2->cd();
  
// residual histogram

double xMin = 40.0;
double xMax = 115.0;
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

       evtsFull,
       evtsFrst,
       evtsScnd

   );
}

delete fullErg;
delete firsErg;
delete secdErg;
delete fcbalpha;
delete fcbndeg;
delete fcbalphaR;   
delete fcbndegR;    

}

