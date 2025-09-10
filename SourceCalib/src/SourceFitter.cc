#include "CaloCalibration/SourceCalib/inc/SourceFitter.hh"

using namespace TMath;
using namespace RooFit;
using namespace CaloSourceCalib;

void SourceFitter::FitCrystal(TH1F* h_spec, TString opt, int crystalNo,  TTree *covar, Int_t &nEvents, Float_t &fpeak, Float_t &dpeak, Float_t &fsigma, Float_t &chiSq, Float_t &fstpeak,Float_t &fstsigma, Float_t &scdpeak,Float_t &scdsigma,Float_t &fcbalphaparam,Float_t &fcbndegparam,Float_t &Aparam,Float_t &Bparam, Float_t &Cparam, Float_t &fullResparam, Float_t &fstResparam,Float_t &scdResparam,Float_t &comCnstparam, Float_t &combetaparam, Float_t &frFullparam, Float_t &frFrstparam,Float_t &frScndparam,Float_t &crystalNoparam,Float_t &frBKGparam, Float_t &convergencestatus,Float_t &errbar, Float_t &pval,Float_t &h_means,Float_t &h_stddevs){//Float_t &frBKGparam
    
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

  // obtain the two initial values for the parameters
  //double par1 = 2500;
  //double par2 = 50;
  //double ADC_conv = 0.0625;

  //parameters
  RooRealVar crysADC("crysADC", "ADC[counts]", 20, 120);//48
  RooRealVar ergElec("ergElec", "electron energy in ADC", 8.176);//0.511

	RooRealVar fcbalpha("fcbalpha", "fcbalpha",0.5, 0.1, 5);
  RooRealVar fcbndeg("fcbndeg", "fcbndeg",10, 1, 40);
  RooRealVar A("A", "coeff of E",0.02,0.001,3);
  RooRealVar B("B", "const B", 0.00405,0.0001,2);
  RooRealVar C("C", "const C",0.0027,0.0001,3);//Electronic noise in MeV 
  
  //Full peak:
  RooRealVar fullPeak("fullPeak", "full peak", 96.8, 80, 108);//6.13
  RooFormulaVar fullRes("fullRes","full peak resolution","0.98650796*sqrt(pow(A/pow(fullPeak/1000,0.25),2)+pow(B/0.0625,2)+pow(C/(fullPeak/1000),2))", RooArgSet(A,B,C,fullPeak));    
  RooFormulaVar fullWidth("fullWidth", "width of the full peak", "fullPeak*fullRes",RooArgSet(fullPeak, fullRes));
  
  // 1st escape peak:
  RooFormulaVar fstEsPeak("fstEsPeak", "first escape peak", "fullPeak - ergElec", RooArgSet(fullPeak,ergElec));
  //"fullPeak - 0.511*(fullpeak/6.13)"
  RooFormulaVar fstRes("fstRes","first peak resolution","0.98650796*sqrt(pow(A/pow(fstEsPeak /1000,0.25),2)+pow(B/0.0625,2)+pow(C/(fstEsPeak/1000),2))", RooArgSet(A,B,C,fstEsPeak));
  RooFormulaVar fstWidth("fstWidth", "width of first escape peak","fstEsPeak*fstRes",RooArgSet(fstEsPeak,fstRes));

  // 2nd escape peak:
  RooFormulaVar scdEsPeak("scdEsPeak", "second escape peak", "fullPeak - 2*ergElec", RooArgSet(fullPeak,ergElec));
  RooFormulaVar scdRes("scdRes","second peak resolution","0.98650796*sqrt(pow(A/pow(scdEsPeak/1000,0.25),2)+pow(B/0.0625,2)+pow(C/(scdEsPeak/1000),2))", RooArgSet(A,B,C,scdEsPeak));
  RooFormulaVar scdWidth("scdWidth", "width of second escape peak","scdEsPeak*scdRes",RooArgSet(scdEsPeak,scdRes));

  // construct CBs:
  RooCBShape fullErg("fullErg", "full energy peak", crysADC, fullPeak, fullWidth, fcbalpha, fcbndeg);
  RooCBShape firsErg("firsErg", "first escape peak", crysADC, fstEsPeak, fstWidth, fcbalpha, fcbndeg);
  RooCBShape secdErg("secdErg", "second escape peak", crysADC, scdEsPeak, scdWidth, fcbalpha, fcbndeg);
  
  // logistic background:
  RooRealVar comCnst("comCnst", "comCnst",  50, 10, 200);//par1, 160, 5120);//~100
  RooRealVar combeta("combeta", "combeta",50, 4, 480);//par2, 4, 480);
  RooGenericPdf comPdf("comPdf", "logistic", "1.0/(1.0+exp((crysADC-comCnst)/combeta))",
             RooArgSet(crysADC, comCnst, combeta));

  // Fraction of events in each peak-- fraction is used instead of number of events bc # of events gave issues in the fit
  RooRealVar frFull("frFull", "Fraction of full peak",0.3,0.1, 1);
  RooRealVar frFrst("frFrst", "Fraction of first escape peak",0.5,0.1, 1);
  RooRealVar frScnd("frScnd", "Fraction of second escape peak", 0.1, 0.1, 1);

  //preparing RooPlot
  RooPlot *chFrame = crysADC.frame(Title(title));
  RooDataHist chSpec("crysADC","crysADC", crysADC, h_spec);
  
  // combined fit function
	RooAddPdf fitFun("fitFun", "firsErg + (secdErg + (fullErg +comPdf))", RooArgList(firsErg, secdErg, fullErg,comPdf), RooArgList(frFrst, frScnd, frFull));//,frBKG) );
	
  if(opt == "chi2"){ //binned chi2 fit
    RooFitResult *fitRes = fitFun.chi2FitTo(chSpec, Range(40,115.2),Hesse(kTRUE),Minos(kTRUE), Strategy(1),MaxCalls(10000),PrintLevel(1),Save(),DataError(RooAbsData::SumW2));
    fitRes->Print("v");
   }  
  if(opt == "nll"){ //binned nll fit
    RooAbsReal* nll = fitFun.createNLL(chSpec, Range(40,115.2));
    RooMinimizer m(*nll);
    m.migrad();
    m.hesse();
    RooFitResult *fitRes = m.save();
    fitRes->Print("v");
	convergencestatus =fitRes->status();    
  }
  // plot components
  chSpec.plotOn(chFrame, MarkerColor(kBlack), LineColor(kBlack), MarkerSize(0.5), Name("chSpec"));
  fitFun.plotOn(chFrame, LineColor(kRed), LineStyle(1), Name("fit"));
  chiSq = chFrame->chiSquare(11);
  fitFun.plotOn(chFrame, Components(fullErg), LineColor(kOrange), LineStyle(5), Name("main"));
  fitFun.plotOn(chFrame, Components(firsErg), LineColor(kViolet), LineStyle(5), Name("fescape"));
  fitFun.plotOn(chFrame, Components(secdErg), LineColor(kCyan), LineStyle(5), Name("sescape"));
  fitFun.plotOn(chFrame, Components(comPdf), LineColor(kBlue), LineStyle(5), Name("background"));

  nEvents = h_spec->GetEntries();
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
  crystalNoparam = crystalNo;
  errbar = (1/(fpeak/6.13))*(dpeak/fpeak);                   
  pval = TMath::Prob(chiSq*11, 11);
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
  TPaveLabel *fsg = new TPaveLabel(0.15, 0.55, 0.25, 0.45, Form("#sigma_{main} = %4.2f", fsigma), "brNDC");
  fsg -> SetFillStyle(0);
  fsg -> SetBorderSize(0);
  fsg -> SetTextSize(0.4);
  fsg -> SetTextFont(42);
  fsg -> SetTextColor(kBlack);
  fsg -> SetFillColor(kWhite);
  chFrame -> addObject(fsg);

  // Create top pad for fit
	TPad *pad1 = new TPad("pad1", "Top pad", 0, 0.25, 1, 1.0);
	pad1->SetBottomMargin(0.035);  // no big gap between pads
	pad1->Draw();
	pad1->cd();  // switch to top pad
  chFrame -> SetYTitle("Events per 25 keV");
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
	RooHist *hpull = chFrame->pullHist();// // (data - fit)/sigma
	hpull->SetTitle("");
	hpull->GetYaxis()->SetTitle("Normalised Residuals (Residual / #sqrt{N})");
	hpull->GetYaxis()->SetTitleSize(0.12);
	hpull->GetYaxis()->SetLabelSize(0.10);
	hpull->GetXaxis()->SetTitleSize(0.12);
	hpull->GetXaxis()->SetLabelSize(0.10);

// Draw residuals in bottom pad
	hpull->Draw("AP");
  can -> SaveAs(oName); 
  can->Close();  // Close the associated file
  delete can;    // Delete the object
  can = nullptr; // Safety (optional)

  covar->Fill();
}
