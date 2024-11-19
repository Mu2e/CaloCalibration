#include "CaloCalibration/SourceCalib/inc/SourceFitter.hh"

using namespace TMath;
using namespace RooFit;
using namespace CaloSourceCalib;

void SourceFitter::FitCrystal(TH1F* h_spec, TString opt, int crystalNo,  TTree *covar, Float_t &fpeak, Float_t &fsigma, Float_t &chiSq, Float_t &fstpeak,Float_t &fstsigma, Float_t &scdpeak,Float_t &scdsigma,Float_t &fcbalphaparam,Float_t &fcbndegparam,Float_t &Aparam,Float_t &Bparam, Float_t &Cparam, Float_t &fullResparam, Float_t &fstResparam,Float_t &scdResparam,Float_t &comCnstparam, Float_t &combetaparam, Float_t &frFullparam, Float_t &frFrstparam,Float_t &frScndparam,Float_t &crystalNoparam,Float_t &frBKGparam){//Float_t &frBKGparam
    
    // set stlye optionsr
    gStyle -> SetOptFit(1111);
    gStyle -> SetOptStat(0);
    gStyle -> SetPadBottomMargin(0.125);
    gStyle -> SetPadTopMargin(0.075);
    gStyle -> SetPadRightMargin(0.05);
    gStyle -> SetPadLeftMargin(0.15);
    gStyle -> SetTitleOffset(1.0, "x");
    gStyle -> SetTitleOffset(1.75, "y");

    TCanvas *can = new TCanvas("can", "", 100, 100, 600, 600);
    can -> Draw();
    TString cryNum = to_string(crystalNo);
    TString oName = "mu2e_simu_fitSpec_"+ opt+"_" + cryNum + ".root";
    //how do i make anotehr branch in this file?
    // need to  make a branch that look at MC energy deposited and maybe another looking at 
    TString title = "Energy Spectrum of Crystal " + cryNum;
    //TFile *inFile = TFile::Open(fName);// can delete

   TFile *ouptFile = new TFile("mu2e_caloSourceFit_" + cryNum + ".root", "RECREATE");

    
    // obtain the two initial values for the parameters
    double par1 = 2500;//100.;1600
    double par2 = 50;//10.;
    //double ergElec_ADC = 8.176;
   // double b_adc = 2.4;
    //double c_adc = 3.2;
   // double b_mev = 0.2;
    //double c_mev = 0.15;
    //double ADC_conv = 0.0625;
    //add b and c as double, also define adc conversion as a constant
    //parameters
    RooRealVar crysADC("crysADC", "ADC[counts]", 48, 120);
    RooRealVar ergElec("ergElec", "electron energy in ADC", 8.176);
		RooRealVar fcbalpha("fcbalpha", "fcbalpha", 0.5, 0.1, 5);
    RooRealVar fcbndeg("fcbndeg", "fcbndeg", 10, 1, 40.);
    RooRealVar A("A", "coeff of E", 0.02,0.001,3);//0.16, 0.08, 16);0.006
    RooRealVar B("B", "const B",0.00405,0.0001,2);//8, 0.08, 16);0.00405
    RooRealVar C("C", "const C", 0.0027,0.0001,3);//Electronic noise in MeV 2.4  0.0027
    
//conv is GEV TO MEV then to ADC for A
    //Full peak:
    RooRealVar fullPeak("fullPeak", "full peak", 96.8, 80, 108);
    //RooFormulaVar fullRes("fullRes","full peak resolution","0.98650796*sqrt(((A)/pow(98.08/1000,0.5))+(B/0.0625)+(C/(98.08/1000)))", RooArgSet(A,B,C));
    RooFormulaVar fullRes("fullRes","full peak resolution","0.98650796*sqrt(pow(A/pow(fullPeak/1000,0.25),2)+pow(B/0.0625,2)+pow(C/(fullPeak/1000),2))", RooArgSet(A,B,C,fullPeak));    
    //RooFormulaVar fullRes("fullRes","full peak resolution","sqrt(pow((A*0.125)/pow(98.08/1000,0.25),2)+pow(B*0.0625,2)+pow(C/98.08,2))", RooArgSet(A,B,C));//old formula 
    RooFormulaVar fullWidth("fullWidth", "width of the full peak", "fullPeak*fullRes",RooArgSet(fullPeak, fullRes));
		//RooRealVar fullWidth("fullWidth", "width of the full peak", 5, 1, 10);
    // 1st escape peak:
    RooFormulaVar fstEsPeak("fstEsPeak", "first escape peak", "fullPeak - ergElec", RooArgSet(fullPeak,ergElec));
    //RooFormulaVar fstRes("fstRes","first peak resolution","0.98650796*sqrt(pow((A)/pow(89.904 /1000,0.5),2)+(B/0.0625)+(C/(89.904/1000) ))", RooArgSet(A,B,C));//,fstEsPeak));//89.904
    RooFormulaVar fstRes("fstRes","first peak resolution","0.98650796*sqrt(pow(A/pow(fstEsPeak /1000,0.25),2)+pow(B/0.0625,2)+pow(C/(fstEsPeak/1000),2))", RooArgSet(A,B,C,fstEsPeak));
    //RooFormulaVar fstRes("fstRes","first peak resolution","sqrt(pow((A*0.125)/pow(89.904 /1000,0.25),2)+pow(B*0.0625,2)+pow(C/89.904 ,2))", RooArgSet(A,B,C)); //old formula 
    RooFormulaVar fstWidth("fstWidth", "width of first escape peak","fstEsPeak*fstRes",RooArgSet(fstEsPeak,fstRes));

    // 2nd escape peak:
    RooFormulaVar scdEsPeak("scdEsPeak", "second escape peak", "fullPeak - 2*ergElec", RooArgSet(fullPeak,ergElec));
    //RooFormulaVar scdRes("scdRes","second peak resolution","0.98650796*sqrt(pow((A)/pow(81.72/1000,0.5),2)+(B/0.0625)+(C/(81.72/1000)))", RooArgSet(A,B,C));//,scdEsPeak));//81.72
    RooFormulaVar scdRes("scdRes","second peak resolution","0.98650796*sqrt(pow(A/pow(scdEsPeak/1000,0.25),2)+pow(B/0.0625,2)+pow(C/(scdEsPeak/1000),2))", RooArgSet(A,B,C,scdEsPeak));
    //RooFormulaVar scdRes("scdRes","second peak resolution","sqrt(pow((A*0.125)/pow(81.72/1000,0.25),2)+pow(B*0.0625,2)+pow(C/81.72,2))", RooArgSet(A,B,C));//old formula 
    RooFormulaVar scdWidth("scdWidth", "width of second escape peak","scdEsPeak*scdRes",RooArgSet(scdEsPeak,scdRes));

    // construct CBs:
    RooCBShape fullErg("fullErg", "full energy peak", crysADC, fullPeak, fullWidth, fcbalpha, fcbndeg);
    RooCBShape firsErg("firsErg", "first escape peak", crysADC, fstEsPeak, fstWidth, fcbalpha, fcbndeg);//fstWidth
    RooCBShape secdErg("secdErg", "second escape peak", crysADC, scdEsPeak, scdWidth, fcbalpha, fcbndeg);//scdWidth
    // logistic background:
    RooRealVar comCnst("comCnst", "comCnst", par1, 160, 5120);//10, 320.//160, 5120.
    RooRealVar combeta("combeta", "combeta", par2, 4, 480);//0.25, 30//4, 480.
    RooGenericPdf comPdf("comPdf", "logistic", "1.0/(1.0+exp((crysADC-comCnst)/combeta))",
               RooArgSet(crysADC, comCnst, combeta));

    // Fraction of events in each peak
    RooRealVar frFull("frFull", "Fraction of full peak",0.3,0.25, 1);
    RooRealVar frFrst("frFrst", "Fraction of first escape peak", 0.5,0.1, 1);
    RooRealVar frScnd("frScnd", "Fraction of second escape peak", 0.1, 0.1, 1);
   // RooRealVar frBKG("frBKG", "Fraction of BKG peak", 0.1, 0.01, 1);
    //RooAddition constrSum("add","add", RooArgList(frFull,frFull,frScnd));
    //RooFormulaVar frBKG("frBKG", "Fraction of BKG peak", "1 - (frFull+frFrst+frScnd)", RooArgSet(frFull,frFull,frScnd));
    //make the above line a double 1-frFull.getval() at the bottom with the other
    // the parameter is 11 after setting width of escape peaks opened
    //Float_t dpeak;//,dsigma,fpeak, fsigma
    //Float_t chiSq = 0;
    Float_t dpeak;//,valcombeta,dcombeta,valcomCnst,dcomCnst,fstpeak,fstsigma,scdpeak,scdsigma,
    //valfrFull,dfrFull,valfrFrst,dfrFrst,valfrScnd,dfrScnd,
    //valfcbndeg,dfcbndeg,valfcbalpha,dfcbalpha,valA,dA,valB,dB,valC,dC;//,valfrBKG,dfrBKG
    //preparing RooPlot
    RooPlot *chFrame = crysADC.frame(Title(title));
   if(opt == "chi2"){ //binned chi2 fit
      //TH1F *h_spec = (TH1F*)inFile->Get("h_spec");
      RooDataHist chSpec("crysADC","crysADC", crysADC, h_spec);
      // combined fit function
      RooAddPdf fitFun("fitFun", "firsErg + (secdErg + (fullErg + comPdf))", RooArgList(firsErg, secdErg, fullErg, comPdf), RooArgList(frFrst, frScnd, frFull) );
			//RooAddPdf fitFun("fitFun", "firsErg*frFrst + secdErg*frScnd + fullErg*frFull + comPdf(1-constrSum)", RooArgList(firsErg, secdErg, fullErg, comPdf), RooArgList(frFrst, frScnd, frFull) );
      //chi2 fit
      RooFitResult *fitRes = fitFun.chi2FitTo(chSpec, Range(48,115.2), Strategy(1),Hesse(kTRUE), PrintLevel(1),Save()) ;//Extended()
       fitRes->Print("v");
       //NLL fit
       //RooFitResult *fitRes = fitFun.fitTo(chSpec,Range(48,115.2), Strategy(2),Hesse(kTRUE), PrintLevel(1),Save());
       //fitRes->Print("v");
       
       //testing minimiser here
       //RooChi2Var chi2("chi2", "chi2", fitFun, chSpec, Range(48, 115.2), Extended());
       //RooMinimizer minimizer(chi2);
			// minimizer.minimize("Minuit2", "Migrad");
			 //minimizer.minos();
			 //minimizer.hesse();
			 //RooFitResult* fitRes = minimizer.save();
			 //fitRes->Print("v");
       //end test

      // plot components
      chSpec.plotOn(chFrame, MarkerColor(kBlack), LineColor(kBlack), MarkerSize(0.5), Name("chSpec"));
      fitFun.plotOn(chFrame, LineColor(kRed), LineStyle(1), Name("fullFit"));
      chiSq = chFrame->chiSquare(11);
      fitFun.plotOn(chFrame, Components(fullErg), LineColor(kOrange), LineStyle(5));
      fitFun.plotOn(chFrame, Components(firsErg), LineColor(kViolet), LineStyle(5));
      fitFun.plotOn(chFrame, Components(secdErg), LineColor(kCyan), LineStyle(5));
      fitFun.plotOn(chFrame, Components(comPdf), LineColor(kBlue), LineStyle(5));
    }
		
    fpeak = fullPeak.getVal();
    dpeak = fullPeak.getError();
    fsigma = fullWidth.getVal();
    fstpeak = fstEsPeak.getVal();
    fstsigma = fstWidth.getVal();
    scdpeak = scdEsPeak.getVal();
    scdsigma = scdWidth.getVal();
    fcbalphaparam = fcbalpha.getVal();
    fcbndegparam = fcbndeg.getVal();
    Aparam = A.getVal();
    Bparam = B.getVal();
    Cparam = C.getVal();
    fullResparam = fullRes.getVal();
    fstResparam = fstRes.getVal();
    scdResparam = scdRes.getVal();
    comCnstparam = comCnst.getVal();
    combetaparam = combeta.getVal();
    frFullparam = frFull.getVal();
    frFrstparam = frFrst.getVal();
    frScndparam = frScnd.getVal();
    frBKGparam= 1-(frFullparam+frFrstparam+frScndparam);
    std::cout<< "fraction of background" <<frBKGparam <<std::endl;
    crystalNoparam = crystalNo;

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
    TPaveLabel *fsg = new TPaveLabel(0.20, 0.775, 0.35, 0.875, Form("sigma = %4.2f", fsigma), "brNDC");//dsigma
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
    //inFile -> Close();
		covar->Fill();
    ouptFile -> Write();
    ouptFile -> Close();
}
