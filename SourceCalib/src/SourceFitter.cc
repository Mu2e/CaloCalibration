#include "CaloCalibration/SourceCalib/inc/SourceFitter.hh"
#include "CaloCalibration/SourceCalib/inc/2dcontour.hh"
using namespace TMath;
using namespace RooFit;
using namespace CaloSourceCalib;

int SourceFitter::nSecondFits = 0; //counter for crystals that need a second fit
int SourceFitter::nThirdFits = 0; // counter for crystals that need a third fit
std::vector<int> SourceFitter::crystalsSecondFit;//vector of crystal number 
std::vector<int> SourceFitter::crystalsThirdFit;
int SourceFitter::nFirstFitConverged = 0;//counter for crystals that converged on first fit
int SourceFitter::nSecondFitConverged = 0; //counter for crystals that need a second fit
int SourceFitter::nThirdFitConverged = 0; // counter for crystals that need a third fit
std::vector<int> SourceFitter::crystalsSecondFitConverged;
std::vector<int> SourceFitter::crystalsThirdFitConverged;
std::map<int,int> SourceFitter::thirdFitRetryCount;
std::vector<std::pair<int,int>> SourceFitter::convFailures; //vector for crystal number of failing status crystals
struct CovarAccumulator {//calculated sum of previous crystals to take avg values to use for second fit
    int count = 0;
    double sumPeak = 0;
    double sumWidth = 0;
    double sumAlpha = 0;
    double sumBeta = 0;
    double sumEvtFull = 0;
    double sumEvtFst = 0;
    double sumEvtScd = 0;
    double sumEvtBkg = 0;
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
void SourceFitter::FitCrystal(TH1F* h_spec, TString opt, int crystalNo,  TTree *covar, Int_t &nEvents,Int_t &convergencestatus, Float_t &fpeak, Float_t &peakerrorhigh,Float_t &peakerrorlo, Float_t &fsigma,Float_t &widtherrorhigh,Float_t &widtherrorlo, Float_t &chiSq, Float_t &fstpeak, Float_t &scdpeak,Float_t &fcbalphaparam,Float_t &fcbndegparam,Float_t &comCnstparam, Float_t &combetaparam, Float_t &fr_fullparam, Float_t &fr_frstparam,Float_t &fr_scndparam,Float_t &crystalNoparam,Float_t &fr_bkgparam, Float_t &pval,Float_t &h_means,Float_t &h_stddevs, Float_t &unreducedchi2,Float_t &fval,Float_t &mparam,Float_t &etaparam,Int_t &ndof,bool contour, Float_t &errbarhigh, Float_t &errbarlo,Float_t &evtfullerrorhigh,Float_t &evtfullerrorlo,Float_t &eventsFull,Float_t &Esparam ){

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
  //Setting initial guess values for all params
  double initPeak  = 98.08;
  double initWidth = 0.5;
  double initAlpha = 0.5;
  double initBeta = 1.0;
  double initEvtFull = 100000;
  double initEvtFst = 90000;
  double initEvtScd = 90000;
  double initEvtBkg= 50000;

  double currentPeak  = initPeak;
  double currentWidth = initWidth;
  double currentAlpha = initAlpha;
  double currentBeta  = initBeta;
  double currentEvtFull = initEvtFull;
  double currentEvtFst  = initEvtFst;
  double currentEvtScd  = initEvtScd;
  double currentEvtBkg  = initEvtBkg;

  RooRealVar m_e("m_e", "electron energy in MeV", 0.511);
  RooRealVar crysADC("crysADC", "ADC [counts]", 40, 120);
  RooRealVar E0("E0", "energy offset [MeV]", 0.0); //might need to remove if not added

  // Crystal Ball shape params
  RooRealVar fcbalpha("fcbalpha", "alpha", initAlpha, 0.1, 5.0);
  RooRealVar fcbndeg("fcbndeg", "n", 5);
  
  // Peak parameters in mev
  RooRealVar fullPeak("fullPeak", "Full peak [ADC]", initPeak, 85, 108);
  RooFormulaVar eta("eta", "ADC/MeV", "fullPeak/ (6.13-E0)", RooArgSet(fullPeak, E0));
  RooFormulaVar m("m","Mev/ADC", "(6.13-E0)/fullPeak", RooArgSet(fullPeak, E0));
  RooFormulaVar Es("Es", "calibrated energy [MeV]", "m*crysADC + E0", RooArgSet(m, crysADC,E0));//Linear energy calibration
  RooFormulaVar fstEsPeak("fstEsPeak", "First escape", "fullPeak - m_e*eta", RooArgSet(fullPeak, m_e,eta));
  RooFormulaVar scdEsPeak("scdEsPeak", "Second escape", "fullPeak - (2*m_e)*eta", RooArgSet(fullPeak, m_e,eta));

  RooRealVar fullWidth("fullWidth", "Full width [MeV]",initWidth,0.2,1.5);
  RooRealVar Egamma("Egamma", "Full peak [MeV]",6.13);
  RooRealVar fstesc("fstesc", "first peak [MeV]",5.619);
  RooRealVar scdesc("scdesc", "second peak [MeV]",5.108);
 
  //three peak crystal ball function
  RooCBShape fullErg("fullErg", "Full peak", Es, Egamma, fullWidth, fcbalpha, fcbndeg);
  RooCBShape firsErg("firsErg", "Single escape", Es,fstesc, fullWidth, fcbalpha, fcbndeg);
  RooCBShape secdErg("secdErg", "Double escape", Es, scdesc, fullWidth, fcbalpha, fcbndeg);

  // Yield
  std::cout << "X axis range: " 
          << h_spec->GetXaxis()->GetXmin() 
          << " to " 
          << h_spec->GetXaxis()->GetXmax() 
          << std::endl;
  int low = h_spec->FindBin(40);
  int high = h_spec->FindBin(120);
  int inegral_evts = h_spec->Integral(low, high);
  nEvents = h_spec->GetEntries();
  RooRealVar evtsFull("evtsFull", "Full peak yield", initEvtFull, 0.05*inegral_evts, inegral_evts);
  RooRealVar evtsFrst("evtsFrst", "First escape yield", initEvtFst, 0, inegral_evts);
  RooRealVar evtsScnd("evtsScnd", "Second escape yield", initEvtScd, 0, inegral_evts);
  RooRealVar evtsbkg("evtsbkg", "Background yield", initEvtBkg, 0, inegral_evts);

  // Background (logistic)
  RooRealVar comCnst("comCnst", "Background const",4.0);//, 0.0, 10.0);
  RooRealVar combeta("combeta", "Background beta", initBeta, 0.001, 200.0);
  RooRealVar bkg_fixed("bkg_fixed", "Background offset", 0.0);//, -10.0, 10.0);
  //logistic background pdf
  RooGenericPdf comPdf("comPdf", "logistic background",
    "pow(1.0 + exp((Es - comCnst)/combeta), -1.0) + bkg_fixed",
    RooArgSet(Es, comCnst, combeta,bkg_fixed));

  // Combine all components
  RooAddPdf fitFun("fitFun", "Total Mu2e Calo Model",
                 RooArgList(fullErg, firsErg, secdErg, comPdf),
                 RooArgList(evtsFull, evtsFrst, evtsScnd, evtsbkg));

                  
  // Fix normalization set to crysADC
  fitFun.fixCoefNormalization(RooArgSet(crysADC));
  //preparing RooPlot
  RooPlot *chFrame = crysADC.frame(Title(title));
  h_spec->Sumw2();
  RooDataHist chSpec("crysADC","crysADC", crysADC, h_spec);
    // Width of the first bin (in keV if the axis is in keV)
	double binWidth = h_spec->GetXaxis()->GetBinWidth(500);  
	std::cout << "Bin width = " << binWidth << " keV" << std::endl; 
//main function that assignes the values to the params
auto run_one_fit = [&]() -> bool {
  // --- Assign starting values ---
    fullPeak.setVal(currentPeak);
    fullWidth.setVal(currentWidth);
    fcbalpha.setVal(currentAlpha);
    combeta.setVal(currentBeta);
    evtsFull.setVal(currentEvtFull);
    evtsFrst.setVal(currentEvtFst);
    evtsScnd.setVal(currentEvtScd);
    evtsbkg.setVal(currentEvtBkg);

    // --- RUN FIT ---
    RooFitResult *fitRes = nullptr;

    if (opt == "chi2") {
        RooChi2Var chi2Func("chi2", "chi2", fitFun, chSpec,
                    RooFit::DataError(RooAbsData::SumW2),
                    RooFit::Range(40, 115.2));
        RooMinimizer m(chi2Func);
        m.setStrategy(2);
        m.setPrintLevel(-1);
        m.migrad();
        m.hesse();
        m.minos();
        fitRes = m.save();
        fitFun.getParameters(chSpec)->assign(fitRes->floatParsFinal());
        unreducedchi2 = chi2Func.getVal(); 
        //calculating ndof
        int nBins = chSpec.numEntries();  // number of bins used in fit
        int nPars = fitRes->floatParsFinal().getSize();
        ndof = nBins - nPars;//defn for arxiv list
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

    // --- Extract peak & width with errors ---
    RooRealVar* peak = dynamic_cast<RooRealVar*>(fitRes->floatParsFinal().find("fullPeak"));
    RooRealVar* width= dynamic_cast<RooRealVar*>(fitRes->floatParsFinal().find("fullWidth"));

    if (peak) {
        fpeak         = peak->getVal(crysADC);
        peakerrorhigh = peak->getAsymErrorHi();
        peakerrorlo   = peak->getAsymErrorLo();
        fullPeak.setVal(peak->getVal(crysADC));        // Update value
    }


    if (width) {
        fsigma        = width->getVal();
        widtherrorhigh= width->getAsymErrorHi();
        widtherrorlo  = width->getAsymErrorLo();
    }

    // --- event yield extraction ---
    RooRealVar* evtfull = dynamic_cast<RooRealVar*>(fitRes->floatParsFinal().find("evtsFull"));
    if (evtfull) {
        eventsFull      = evtfull->getVal();
        evtfullerrorhigh= evtfull->getAsymErrorHi();
        evtfullerrorlo  = evtfull->getAsymErrorLo();
    }
    return (convergencestatus == 0);
};

  // -----------------------
  // 1) First attempt: start from init values 
  // -----------------------
bool first_ok = run_one_fit(); //run with init guesses
//add to counter to calculae avg of all values
covarAcc.count++;
covarAcc.sumPeak     += fpeak;
covarAcc.sumWidth    += fsigma;
covarAcc.sumAlpha    += fcbalphaparam;
covarAcc.sumBeta     += combetaparam;
covarAcc.sumEvtFull  += fr_fullparam*nEvents;
covarAcc.sumEvtFst   += fr_frstparam*nEvents;
covarAcc.sumEvtScd   += fr_scndparam*nEvents;
covarAcc.sumEvtBkg   += fr_bkgparam*nEvents;
bool second_ok = true;

/* -------------------------------------------------------------------
   FIRST FIT FAILED --> REFIT USING COVARIANCE MEAN
   ------------------------------------------------------------------- */
if (!first_ok)
{
    std::cout << "======================================================================\n";
    std::cout << "needs refit (crystal " << crystalNo << ")\n";
    std::cout << "======================================================================\n";

    SourceFitter::nSecondFits++;
    SourceFitter::crystalsSecondFit.push_back(crystalNo);

    double newPeakGuess     = covarAcc.sumPeak     / covarAcc.count;
    double newSigmaGuess    = covarAcc.sumWidth    / covarAcc.count;
    double newAlphaGuess    = covarAcc.sumAlpha    / covarAcc.count;
    double newCombetaGuess  = covarAcc.sumBeta     / covarAcc.count;
    double newevtsFullGuess = covarAcc.sumEvtFull  / covarAcc.count;
    double newevtsFstGuess  = covarAcc.sumEvtFst   / covarAcc.count;
    double newevtsScdGuess  = covarAcc.sumEvtScd   / covarAcc.count;
    double newevtsBkgGuess  = covarAcc.sumEvtBkg   / covarAcc.count;
    //set new values for the starting point of params
    currentPeak   = newPeakGuess;
    currentWidth  = newSigmaGuess;
    currentAlpha  = newAlphaGuess;
    currentBeta   = newCombetaGuess;
    currentEvtFull = newevtsFullGuess;
    currentEvtFst  = newevtsFstGuess;
    currentEvtScd  = newevtsScdGuess;
    currentEvtBkg  = newevtsBkgGuess;

    /* ---- RUN SECOND FIT ---- */
    second_ok = run_one_fit();//run with updated values

    if (second_ok) {
        SourceFitter::nSecondFitConverged++;
        SourceFitter::crystalsSecondFitConverged.push_back(crystalNo);
        std::cout << "[SourceFitter] Second fit SUCCESS for crystal " << crystalNo << "\n";
    } 
    else {  // second fit failed, try third fit
        SourceFitter::nThirdFits++;
        SourceFitter::crystalsThirdFit.push_back(crystalNo);

        // Lambda: randomize all parameters
        auto randomize_all_parameters = [&]() {
            auto randomDouble = [](double min, double max) {
                return min + (max - min) * ((double)rand() / RAND_MAX);
            };

            currentPeak     = randomDouble(85.0, 108.0);
            currentWidth    = randomDouble(0.2, 1.5);
            currentAlpha    = randomDouble(0.6, 1.7);
            currentBeta     = randomDouble(0.2, 1.6);
            currentEvtFull  = randomDouble(0.05 * inegral_evts, inegral_evts);
            currentEvtFst   = randomDouble(0.05 * inegral_evts, inegral_evts);
            currentEvtScd   = randomDouble(0.05 * inegral_evts, inegral_evts);
            currentEvtBkg   = randomDouble(0.05 * inegral_evts, inegral_evts);
        };

        // Lambda: retry third fit
        const int MAX_THIRD_TRIES = 40;
        auto retry_third_fit = [&](int crystal) -> bool {
            int tries = 0;
            bool converged = false;

            for (tries = 1; tries <= MAX_THIRD_TRIES; ++tries) {
                randomize_all_parameters();

                if (run_one_fit()) {
                    converged = true;
                    break;
                }
            }

            SourceFitter::thirdFitRetryCount[crystal] = tries;
            return converged;
        };

        // Run third fit
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

/* -------------------------------------------------------------------
   FIRST FIT OK (simple case)
   ------------------------------------------------------------------- */
else {
    std::cout << "[SourceFitter] Fit ok for crystal " << crystalNo
              << " fpeak=" << fpeak << "\n";
    nFirstFitConverged++;
}


/* -------------------------------------------------------------------
   Record convergence status
   ------------------------------------------------------------------- */
if (convergencestatus > 0) {
    SourceFitter::convFailures.push_back({crystalNo, convergencestatus});
}
//check to make sure the plot is cleared 
std::cout << "Before clear: " << chFrame->numItems() << std::endl;

  while (chFrame->numItems() > 0) {
    chFrame->remove();   // removes last-added item
}
  std::cout << "After clear: " << chFrame->numItems() << std::endl;
  //get values for parameters that will populate the arxiv table
  fstpeak = fstEsPeak.getVal();
  scdpeak = scdEsPeak.getVal();
  fcbalphaparam = fcbalpha.getVal();
  fcbndegparam = fcbndeg.getVal();
  comCnstparam =comCnst.getVal();
  combetaparam = combeta.getVal();
  fr_fullparam = evtsFull.getVal()/inegral_evts;
  fr_frstparam = evtsFrst.getVal()/inegral_evts;
  fr_scndparam = evtsScnd.getVal()/inegral_evts;
  fr_bkgparam= evtsbkg.getVal()/inegral_evts;
  crystalNoparam = crystalNo;                   
  pval = TMath::Prob(chiSq*11, 11);
  mparam = m.getVal();
  etaparam = eta.getVal();
  errbarhigh = mparam*(peakerrorhigh/fpeak);
  errbarlo = mparam*(peakerrorlo/fpeak);
  Esparam = Es.getVal();
  
  //make plot with components
  chSpec.plotOn(chFrame, MarkerColor(kBlack), LineColor(kBlack), MarkerSize(0.5), Name("chSpec"));
  fitFun.plotOn(chFrame, LineColor(kRed), LineStyle(1), Name("fit"));
  chiSq = chFrame->chiSquare("fit", "chSpec", 8);
  std::cout<<"reduced chi2 calculated by func"<<chiSq<<std::endl;
  fitFun.plotOn(chFrame, Components(fullErg), LineColor(kOrange), LineStyle(5), Name("main"));
  fitFun.plotOn(chFrame, Components(firsErg), LineColor(kViolet), LineStyle(5), Name("fescape"));
  fitFun.plotOn(chFrame, Components(secdErg), LineColor(kCyan), LineStyle(5), Name("sescape"));
  fitFun.plotOn(chFrame, Components(comPdf), LineColor(kBlue), LineStyle(5), Name("background"));

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
  TPaveLabel *fpk = new TPaveLabel(0.15, 0.65, 0.25, 0.55, Form("#mu_{main} = %.2f^{+%.2f}_{-%.2f}", fpeak, peakerrorhigh,peakerrorlo), "brNDC");
  std::cout<<"fpeak value at plotting" << fpeak<<std::endl;
  fpk -> SetFillStyle(0);
  fpk -> SetBorderSize(0);
  fpk -> SetTextSize(0.4);
  fpk -> SetTextFont(42);
  fpk -> SetTextColor(kBlack);
  fpk -> SetFillColor(kWhite);
  chFrame -> addObject(fpk);
  TPaveLabel *fsg = new TPaveLabel(0.15, 0.55, 0.25, 0.45, Form("#sigma_{main} =%.2f^{+%.4f}_{-%.4f}", fsigma,widtherrorhigh,widtherrorlo), "brNDC");
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
  chFrame -> GetYaxis()->SetRangeUser(0, 20000);
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
// =========================================================
// CONTOUR CONFIGURATION "SWITCH BOARD"
// =========================================================
// =========================================================
// CONTOUR CONFIGURATION "SWITCH BOARD"
// =========================================================
if (contour) {
    // Configure Plot Here -- Just change these two strings to pick your X and Y axes.
   // Options: "Peak", "Width", "Alpha", "N_Full", "N_1st", 
   //          "N_2nd", "N_Bkg", "Const", "Beta"
    TString xSelect = "Peak";    // <--- Change this to "Alpha", "N_Full", etc.
    TString ySelect = "Width";   // <--- Change this to "Peak", "Const", etc.

    // --------------------------------------------------------
    // Run Plotter
    // --------------------------------------------------------
    CaloSourceCalib::MakeContourPlot(
        fitFun, chSpec, opt, crystalNo,
        xSelect, ySelect,         // The selection
        // Special Manual Vars
        fullPeak,  fpeak,  peakerrorlo,  peakerrorhigh,
        fullWidth, fsigma, widtherrorlo, widtherrorhigh,
        // Standard Vars (Pass the RooRealVars directly)
        fcbalpha,
        evtsFull,
        evtsFrst,
        evtsScnd,
        evtsbkg,
        comCnst,
        combeta
    );
}

}

