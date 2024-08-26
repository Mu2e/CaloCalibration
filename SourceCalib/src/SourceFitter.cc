#include "CaloCalibration/SourceCalib/inc/SourceFitter.hh"

using namespace TMath;
using namespace RooFit;
using namespace CaloSourceCalib;


void SourceFitter::FitCrystal(TH1F* h_spec, TString opt, int crystalNo,  TTree *covar, Float_t &fpeak, Float_t &fsigma, Float_t &chiSq){
    
    // set stlye optionsr
    gStyle -> SetOptFit(1111);
    gStyle -> SetOptStat(0);
    gStyle -> SetPadBottomMargin(0.125);
    gStyle -> SetPadTopMargin(0.075);
    gStyle -> SetPadRightMargin(0.05);
    gStyle -> SetPadLeftMargin(0.15);
    gStyle -> SetTitleOffset(1.0, "x");
    gStyle -> SetTitleOffset(1.75, "y");

    //TCanvas *can = new TCanvas("can", "", 100, 100, 600, 600);--uncomment for pull
    //can -> Draw();--uncomment for pull
    TString cryNum = to_string(crystalNo);
    //TString fName = "mu2e_caloSimu_crysEdep_"+cryNum+".root" ;--uncomment for pull
   // TString oName = "mu2e_simu_fitSpec_"+ opt+"_" + cryNum + ".root";--uncomment for pull
    TString title = "Energy Spectrum of Crystal " + cryNum;
    //TFile *inFile = TFile::Open(fName);

   // TFile *ouptFile = new TFile("mu2e_caloSourceFit_" + cryNum + ".root", "RECREATE");--uncomment for pull

    
    // obtain the two initial values for the parameters
    double par1 = 1600;//100.;
    double par2 = 160;//10.;
    
    //parameters
    RooRealVar crysADC("crysADC", "ADC[counts]", 48, 112);
    RooRealVar ergElec("ergElec", "electron energy", 8.176);
    RooRealVar fcbalpha("fcbalpha", "fcbalpha", 40, 0.8, 320);
    RooRealVar fcbndeg("fcbndeg", "fcbndeg", 160., 4, 1280.);
    RooRealVar A("A", "coeff of E", 0.16, 0.08, 16);//0.01, 0.005, 1);
    RooRealVar B("B", "const B", 8, 0.08, 16);//0.5, 0.005, 1);
    RooRealVar C("C", "const C", 2.4);//Electronic noise in MeV
    std::cout<<"fcbndeg"<<fcbndeg<<std::endl;

    //Full peak:
    RooRealVar fullPeak("fullPeak", "full peak", 96.8, 80, 108);
    //RooFormulaVar fullRes("fullRes","full peak resolution","sqrt(pow(A/pow(fullPeak/1000,0.25),2))+B**2+pow(C/fullPeak,2)", RooArgSet(A, fullPeak, B, C));
    //RooFormulaVar fullWidth("fullWidth", "width of the full peak", "fullPeak*fullRes",RooArgSet(fullPeak,fullRes));
    RooRealVar fullWidth("fullWidth", "width of the full peak", 5.6, 2.4,11.2);//0.35, 0.15, 0.70);
    std::cout<<"fullWidth"<<fullWidth<<std::endl;
    // 1st escape peak:
    RooFormulaVar fstEsPeak("fstEsPeak", "first escape peak", "fullPeak - ergElec", RooArgSet(fullPeak, ergElec));
    //RooFormulaVar fstRes("fstRes","first peak resolution","sqrt(pow(A/pow(fstEsPeak/1000,0.25),2))+B**2+pow(C/fstEsPeak,2)", RooArgSet(A, fstEsPeak, B, C));
    //RooFormulaVar fstWidth("fstWidth", "width of first escape peak","fstEsPeak*fstRes",RooArgSet(fstEsPeak,fstRes));
    RooRealVar fstWidth("fstWidth", "width of first escape peak",5.6, 2.4,11.2);//0.35, 0.15, 0.70);
    std::cout<<"fstWidth"<<fstWidth<<std::endl;
    // 2nd escape peak:
    RooFormulaVar scdEsPeak("scdEsPeak", "second escape peak", "fullPeak - 2*ergElec", RooArgSet(fullPeak, ergElec));
    //RooFormulaVar scdRes("scdRes","second peak resolution","sqrt(pow(A/pow(scdEsPeak/1000,0.25),2))+B**2+pow(C/scdEsPeak,2)", RooArgSet(A, scdEsPeak, B, C));
    //RooFormulaVar scdWidth("scdWidth", "width of second escape peak","scdEsPeak*scdRes",RooArgSet(scdEsPeak,scdRes));
    RooRealVar scdWidth("scdWidth", "width of second escape peak", 5.6, 2.4,11.2);//0.35, 0.15, 0.70);
    std::cout<<"scdWidth"<<scdWidth<<std::endl;
    // construct CBs:
    RooCBShape fullErg("fullErg", "full energy peak", crysADC, fullPeak, fullWidth, fcbalpha, fcbndeg);
    RooCBShape firsErg("firsErg", "first escape peak", crysADC, fstEsPeak, fstWidth, fcbalpha, fcbndeg);
    RooCBShape secdErg("secdErg", "second escape peak", crysADC, scdEsPeak, scdWidth, fcbalpha, fcbndeg);
    std::cout<<"secdErg"<<secdErg<<std::endl;
    // logistic background:
    RooRealVar comCnst("comCnst", "comCnst", par1, 160, 5120);//10, 320.//160, 5120.
    RooRealVar combeta("combeta", "combeta", par2, 4, 480);//0.25, 30//4, 480.
    RooGenericPdf comPdf("comPdf", "logistic", "1.0/(1.0+exp((crysADC-comCnst)/combeta))",
               RooArgSet(crysADC, comCnst, combeta));
    std::cout<<"combeta"<<combeta<<std::endl;
    // Fraction of events in each peak
    RooRealVar frFull("frFull", "Fraction of full peak",0.5, 0.25, 0.7);// 8, 4, 11.2
    RooRealVar frFrst("frFrst", "Fraction of first escape peak", 0.3, 0.1, 0.5);//4.8, 1.6, 8
    RooRealVar frScnd("frScnd", "Fraction of second escape peak", 0.2, 0.1, 0.5);//3.2, 1.6, 8)
    std::cout<<"frScnd"<<frScnd<<std::endl;
    RooRealVar frBKG("frBKG", "Fraction of BKG peak", 0.05, 0.1, 0.0);
    RooAddition constrSum("add","add", RooArgList(frFull,frFrst,frScnd,frBKG));
    // the parameter is 11 after setting width of escape peaks opened
    Float_t dpeak,dsigma;//,dsigma,fpeak, fsigma
    //Float_t chiSq = 0;
    //preparing RooPlot
    RooPlot *chFrame = crysADC.frame(Title(title));
   if(opt == "nll"){ //binned chi2 fit
      //TH1F *h_spec = (TH1F*)inFile->Get("h_spec");
      RooDataHist chSpec("crysADC","crysADC", crysADC, h_spec);
      // combined fit function
      RooAddPdf fitFun("fitFun", "firsErg + (secdErg + (fullErg + comPdf))", RooArgList(firsErg, secdErg, fullErg, comPdf), RooArgList(frFrst, frScnd, frFull) );

      //chi2 fit
      //RooFitResult *fitRes = fitFun.chi2FitTo(chSpec, Range(48,115.2), Strategy(3), PrintLevel(1), Hesse(kTRUE), Extended(),Save());
      //nll fit
      RooFitResult *fitRes = fitFun.fitTo(chSpec,Extended(kTRUE)) ;
      // plot components
      chSpec.plotOn(chFrame, MarkerColor(kBlack), LineColor(kBlack), MarkerSize(0.5), Name("chSpec"));
      fitFun.plotOn(chFrame, LineColor(kRed), LineStyle(1), Name("fullFit"));
      chiSq = chFrame->chiSquare(11);
      fitFun.plotOn(chFrame, Components(fullErg), LineColor(kOrange), LineStyle(5));
      fitFun.plotOn(chFrame, Components(firsErg), LineColor(kViolet), LineStyle(5));
      fitFun.plotOn(chFrame, Components(secdErg), LineColor(kCyan), LineStyle(5));
      fitFun.plotOn(chFrame, Components(comPdf), LineColor(kBlue), LineStyle(5));
      std::cout<<"Result "<< fitRes <<std::endl;
    }
    
    //get fit parameters
    fpeak = fullPeak.getVal();
    dpeak = fullPeak.getError();
    fsigma = fullWidth.getVal();
    dsigma = fullWidth.getError();
    std::cout<<"dsigma"<<dsigma<<std::endl;
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
    TPaveLabel *fsg = new TPaveLabel(0.20, 0.775, 0.35, 0.875, Form("sigma = %.2f#pm%.2f", fsigma,dsigma), "brNDC");
    fsg -> SetFillStyle(0);
    fsg -> SetBorderSize(0);
    fsg -> SetTextSize(0.25);
    fsg -> SetTextColor(kBlack);
    fsg -> SetFillColor(kWhite);
    chFrame -> addObject(fsg);
    std::cout << "chi2: " << chiSq << "; Probability: " << Prob(chiSq, 151) << std::endl;

    chFrame -> SetYTitle("Events per 25 keV");
    chFrame -> Draw();
    //can -> SaveAs(oName);--uncomment for pull

    // the file must be closed
    //inFile -> Close();
		covar->Fill();
    //ouptFile -> Write();--uncomment for pull
    //ouptFile -> Close();--uncomment for pull
}
