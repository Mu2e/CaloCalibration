#include "CaloCalibration/SourceCalib/inc/SourceFitter.hh"

using namespace TMath;
using namespace RooFit;
using namespace CaloSourceCalib;


void SourceFitter::FitCrystal(int crystalNo, TString opt){

    // set stlye options
    gStyle -> SetOptFit(1111);
    gStyle -> SetOptStat(0);
    gStyle -> SetPadBottomMargin(0.125);
    gStyle -> SetPadTopMargin(0.075);
    gStyle -> SetPadRightMargin(0.05);
    gStyle -> SetPadLeftMargin(0.15);
    gStyle -> SetTitleOffset(1.0, "x");
    gStyle -> SetTitleOffset(1.75, "y");


    //for(unsigned int crystalNo = startCry; crystalNo < endCry; crystalNo++){

    TCanvas *can = new TCanvas("can", "", 100, 100, 600, 600);
    can -> Draw();
    TString cryNum = to_string(crystalNo);
    TString fName = "mu2e_caloSimu_crysEdep_"+cryNum+".root" ;
    TString oName = "mu2e_simu_fitSpec_"+ opt+"_" + cryNum + ".root";
    TString title = "Energy Spectrum of Crystal " + cryNum;
    TFile *inFile = TFile::Open(fName);

    TFile *ouptFile = new TFile("mu2e_caloSourceFit_" + cryNum + ".root", "RECREATE");

    
    // obtain the two initial values for the parameters
    double par1 = 100.;
    double par2 = 10.;

    //parameters
    RooRealVar crysEdep("crysEdep", "E_{reco} [MeV]", 3.0, 7.0);
    RooRealVar ergElec("ergElec", "electron energy", 0.511);
    RooRealVar fcbalpha("fcbalpha", "fcbalpha", 2.5, 0.05, 20.0);
    RooRealVar fcbndeg("fcbndeg", "fcbndeg", 10., 0.25, 80.);

    //Full peak:
    RooRealVar fullPeak("fullPeak", "full peak", 6.05, 5.0, 6.75);
    RooRealVar fullWidth("fullWidth", "width of the full peak", 0.35, 0.15, 0.70);

    // 1st escape peak:
    RooFormulaVar fstEsPeak("fstEsPeak", "first escape peak", "fullPeak - ergElec", RooArgSet(fullPeak, ergElec));
    RooRealVar fstWidth("fstWidth", "width of first escape peak", 0.35, 0.15, 0.70);

    // 2nd escape peak:
    RooFormulaVar scdEsPeak("scdEsPeak", "second escape peak", "fullPeak - 2*ergElec", RooArgSet(fullPeak, ergElec));
    RooRealVar scdWidth("scdWidth", "width of second escape peak", 0.35, 0.15, 0.70);

    // construct CBs:
    RooCBShape fullErg("fullErg", "full energy peak", crysEdep, fullPeak, fullWidth, fcbalpha, fcbndeg);
    RooCBShape firsErg("firsErg", "first escape peak", crysEdep, fstEsPeak, fstWidth, fcbalpha, fcbndeg);
    RooCBShape secdErg("secdErg", "second escape peak", crysEdep, scdEsPeak, scdWidth, fcbalpha, fcbndeg);

    // logistic background:
    RooRealVar comCnst("comCnst", "comCnst", par1, 10, 320.);
    RooRealVar combeta("combeta", "combeta", par2, 0.25, 30);
    RooGenericPdf comPdf("comPdf", "logistic", "1.0/(1.0+exp((crysEdep-comCnst)/combeta))",
               RooArgSet(crysEdep, comCnst, combeta));

    // Fraction of events in each peak
    RooRealVar frFull("frFull", "Fraction of full peak", 0.5, 0.25, 0.7);
    RooRealVar frFrst("frFrst", "Fraction of first escape peak", 0.3, 0.1, 0.5);
    RooRealVar frScnd("frScnd", "Fraction of second escape peak", 0.2, 0.1, 0.5);
    RooAddition constrSum("add","add", RooArgList(frFull,frFrst,frScnd));
    
    //RooRealVar frBKG("frBKG", "Fraction of BKG peak", 0.05, 0.1, 0.0);

    // the parameter is 11 after setting width of escape peaks opened
    Float_t fsigma, dsigma, fpeak, dpeak;
    Float_t chiSq = 0;
    //preparing RooPlot
    RooPlot *chFrame = crysEdep.frame(Title(title));

   if(opt == "chi2"){ //binned chi2 fit
      TH1F *h_spec = (TH1F*)inFile->Get("h_spec");
      RooDataHist chSpec("chSpec", "chSpec", crysEdep, h_spec);
      // combined fit function
      RooAddPdf fitFun("fitFun", "firsErg + (secdErg + (fullErg + comPdf))", RooArgList(firsErg, secdErg, fullErg, comPdf), RooArgList(frFrst, frScnd, frFull) );

      //chi2 fit
      RooFitResult *fitRes = fitFun.chi2FitTo(chSpec, Range(3.0, 7.2), Strategy(3), PrintLevel(1), Hesse(kTRUE), Extended(),Save());

      // plot components
      chSpec.plotOn(chFrame, MarkerColor(kBlack), LineColor(kBlack), MarkerSize(0.5), Name("chSpec"));
      fitFun.plotOn(chFrame, LineColor(kRed), LineStyle(1), Name("fullFit"));
      chiSq = chFrame->chiSquare(11);
      fitFun.plotOn(chFrame, Components(fullErg), LineColor(kOrange), LineStyle(5));
      fitFun.plotOn(chFrame, Components(firsErg), LineColor(kViolet), LineStyle(5));
      fitFun.plotOn(chFrame, Components(secdErg), LineColor(kCyan), LineStyle(5));
      fitFun.plotOn(chFrame, Components(comPdf), LineColor(kBlue), LineStyle(5));

      std::cout<<"Result "<< fitRes <<std::endl;
    } if(opt == "nll"){ //unbinned nll fit
      TTree *inTree = (TTree*)inFile->Get("crysTree");
      float crysEdepV;
      inTree -> SetBranchAddress("crysEdep", &crysEdepV);//ncalhitHit
    
      RooDataSet chSpec("chSpec", "chSpec", RooArgSet(crysEdep),  Import(*inTree));

      // combined fit function
      RooAddPdf fitFun("fitFun", "firsErg + (secdErg + (fullErg + comPdf))", RooArgList(firsErg, secdErg, fullErg, comPdf), RooArgList(frFrst, frScnd, frFull), kTRUE);

      //NLL minimizer
      /*RooAbsReal* nll = fitFun.createNLL(chSpec, Range(3.0, 7.2));
      RooMinimizer m(*nll);
      m.migrad();
      m.hesse();*/
      RooFitResult *fitRes = fitFun.fitTo(chSpec,Extended(kTRUE)) ;

      // plot components
      chSpec.plotOn(chFrame, MarkerColor(kBlack), LineColor(kBlack), MarkerSize(0.5), Name("chSpec"));
      fitFun.plotOn(chFrame, LineColor(kRed), LineStyle(1), Name("fullFit"));
      chiSq = chFrame->chiSquare(11);
      fitFun.plotOn(chFrame, Components(fullErg), LineColor(kOrange), LineStyle(5));
      fitFun.plotOn(chFrame, Components(firsErg), LineColor(kViolet), LineStyle(5));
      fitFun.plotOn(chFrame, Components(secdErg), LineColor(kCyan), LineStyle(5));
      fitFun.plotOn(chFrame, Components(comPdf), LineColor(kBlue), LineStyle(5));

      //RooFitResult *fitRes = m.save();
      std::cout<<"Result "<< fitRes <<std::endl;
    }
    
    //get fit parameters
    fpeak = fullPeak.getVal();
    dpeak = fullPeak.getError();
    fsigma = fullWidth.getVal();
    dsigma = fullWidth.getError();

    //make pretty plots
    TPaveLabel *pchi2 = new TPaveLabel(0.20, 0.70, 0.35, 0.80, Form("#chi^{2}/ndf = %4.2f", chiSq), "brNDC");
    pchi2 -> SetFillStyle(0);
    pchi2 -> SetBorderSize(0);
    pchi2 -> SetTextSize(0.25);
    pchi2 -> SetTextColor(kBlack);
    pchi2 -> SetFillColor(kWhite);
    chFrame -> addObject(pchi2);
    TPaveLabel *fpk = new TPaveLabel(0.20, 0.85, 0.35, 0.95, Form("peak = %.2f#pm%.2f", fpeak, dpeak), "brNDC");
    fpk -> SetFillStyle(0);
    fpk -> SetBorderSize(0);
    fpk -> SetTextSize(0.25);
    fpk -> SetTextColor(kBlack);
    fpk -> SetFillColor(kWhite);
    chFrame -> addObject(fpk);
    TPaveLabel *fsg = new TPaveLabel(0.20, 0.775, 0.35, 0.875, Form("sigma = %.2f#pm%.2f", fsigma, dsigma), "brNDC");
    fsg -> SetFillStyle(0);
    fsg -> SetBorderSize(0);
    fsg -> SetTextSize(0.25);
    fsg -> SetTextColor(kBlack);
    fsg -> SetFillColor(kWhite);
    chFrame -> addObject(fsg);
    std::cout << "chi2: " << chiSq << "; Probability: " << Prob(chiSq, 151) << std::endl;

    chFrame -> SetYTitle("Events per 25 keV");
    chFrame -> Draw();
    can -> SaveAs(oName);

    // the file must be closed
    inFile -> Close();

    ouptFile -> Write();
    ouptFile -> Close();
}
