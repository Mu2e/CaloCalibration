//to run:  root .L MyFitv2.C
#include <fstream>
#include <iostream>
#include "TSystem.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TMath.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TMath.h"
#include <Riostream.h>
#include "TLorentzVector.h"
#include "TPaveStats.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TPaveLabel.h"
#include "TAttFill.h"
// add roofit header files
#include "RooHist.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooKeysPdf.h"
#include "RooHistPdf.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooArgusBG.h"
#include "RooPolynomial.h"
#include "RooCBShape.h"
#include "RooGenericPdf.h"
#include "RooFFTConvPdf.h"
#include "RooFitResult.h"
using namespace std;
using namespace TMath;
using namespace RooFit;



int FitMacro()
{
  gROOT -> Reset();
  gSystem -> Load("libRooFit.so");
  gStyle -> SetOptFit(1111);
  gStyle -> SetOptStat(0);
  gStyle -> SetPadBottomMargin(0.125);
  gStyle -> SetPadTopMargin(0.075);
  gStyle -> SetPadRightMargin(0.05);
  gStyle -> SetPadLeftMargin(0.15);
  gStyle -> SetTitleOffset(1.0, "x");
  gStyle -> SetTitleOffset(1.75, "y");

  TFile *ouptFile = new TFile("paraFile.root", "RECREATE");
  Float_t fsigma, dsigma, fpeak, dpeak, chiSq, fullresolution, vfrFull, vfrFirst, vfrSecond;
  Float_t firstsigma, firstdsigma, firstpeak, firstdpeak, firstresolution;
  Float_t secondsigma, seconddsigma, secondpeak, seconddpeak, secondresolution;
  TTree *newTree = new TTree("newTree","newTree");
  newTree->Branch("FullPeak", &fpeak,"fpeak/F");
  newTree->Branch("FullWidth", &fsigma,"fsigma/F");
  newTree->Branch("FullDPeak", &dpeak,"dpeak/F");
  newTree->Branch("FullDSigma", &dsigma,"dsigma/F");
  newTree->Branch("Fullresolution", &fullresolution,"fullresolution/F");
  newTree->Branch("ChiSquared", &chiSq,"chiSq/F");
  newTree->Branch("vfrFull", &vfrFull,"vfrFull/F");
  newTree->Branch("vfrFirst", &vfrFirst,"vfrFirst/F");
  newTree->Branch("vfrSecond", &vfrSecond,"vfrSecond/F");

  newTree->Branch("FirstPeak", &firstpeak,"firstpeak/F");
  newTree->Branch("FirstWidth", &firstsigma,"firstsigma/F");
  newTree->Branch("FirstDPeak", &firstdpeak,"firstdpeak/F");
  newTree->Branch("FirstDSigma", &firstdsigma,"firstdsigma/F");
  newTree->Branch("Firstresolution", &firstresolution,"firstresolution/F");

  newTree->Branch("SecondPeak", &secondpeak,"secondpeak/F");
  newTree->Branch("SecondWidth", &secondsigma,"secondsigma/F");
  newTree->Branch("SecondDPeak", &seconddpeak,"seconddpeak/F");
  newTree->Branch("SecondDSigma", &seconddsigma,"seconddsigma/F");
  newTree->Branch("Secondresolution", &secondresolution,"secondresolution/F");

  for(unsigned int crystalNo = 674; crystalNo < 1348; crystalNo++){
    std::cout<<"================================="<<crystalNo<<"================================="<<std::endl;
    TCanvas *can = new TCanvas("can", "", 100, 100, 600, 600);
    can -> Draw();
    TString cryNum = to_string(crystalNo);
    TString fName = "caloSims/mu2e_caloSimu_crySpec" + cryNum + ".root";
    TString oName = "mu2e_simu_fitSpec" + cryNum + ".eps";
    TString tName = "mu2e_simu_ergTrth" + cryNum + ".eps";
    TString title = "Energy Spectrum of Crystal " + cryNum;
    TFile *fspec = TFile::Open(fName);
    TH1F *h_spec = (TH1F*)fspec->Get("h_spec");
    // obtain the two initial values for the parameters
    double par1 = 100.;
    double par2 = 10.;

    RooRealVar crySpec("crySpec", "E_{reco} [MeV]", 3.0, 7.0);
    RooDataHist chSpec("chSpec", "chSpec", crySpec, h_spec);

    RooRealVar ergElec("ergElec", "electron energy", 0.511);
    RooRealVar fcbalpha("fcbalpha", "fcbalpha", 2.5, 0.05, 20.0);
    RooRealVar fcbndeg("fcbndeg", "fcbndeg", 10., 0.25, 80.);

    //Full:
    RooRealVar fullPeak("fullPeak", "full peak", 6.05, 5.0, 6.75);
    RooRealVar fullWidth("fullWidth", "width of the full peak", 0.35, 0.15, 0.70);

    // 1st escape:
    RooFormulaVar fstEsPeak("fstEsPeak", "first escape peak", "fullPeak - ergElec", RooArgSet(fullPeak, ergElec));
    RooRealVar fstWidth("fstWidth", "width of first escape peak", 0.35, 0.15, 0.70);

    // 2nd escape:
    RooFormulaVar scdEsPeak("scdEsPeak", "second escape peak", "fullPeak - 2*ergElec", RooArgSet(fullPeak, ergElec));
    RooRealVar scdWidth("scdWidth", "width of second escape peak", 0.35, 0.15, 0.70);

    // CB construt:
    RooCBShape fullErg("fullErg", "full energy peak", crySpec, fullPeak, fullWidth, fcbalpha, fcbndeg);
    RooCBShape firsErg("firsErg", "first escape peak", crySpec, fstEsPeak, fstWidth, fcbalpha, fcbndeg);
    RooCBShape secdErg("secdErg", "second escape peak", crySpec, scdEsPeak, scdWidth, fcbalpha, fcbndeg);

    // logistic background
    RooRealVar comCnst("comCnst", "comCnst", par1, 10, 320.);
    RooRealVar combeta("combeta", "combeta", par2, 0.25, 30);
    RooGenericPdf comPdf("comPdf", "logistic", "1.0/(1.0+exp((crySpec-comCnst)/combeta))",
                         RooArgSet(crySpec, comCnst, combeta));

    // Fraction of events in eack peak
    RooRealVar frFull("frFull", "Fraction of full peak", 0.8, 0.25, 1.0);
    RooRealVar frFrst("frFrst", "Fraction of first escape peak", 0.3, 0.1, 0.7);
    RooRealVar frScnd("frScnd", "Fraction of second escape peak", 0.8, 0.25, 1.0);

    RooAddPdf fitFun("fitFun", "firsErg + (secdErg + (fullErg + comPdf))", RooArgList(firsErg, secdErg, fullErg, comPdf), RooArgList(frFrst, frScnd, frFull), kTRUE);

    RooFitResult *fitRes = fitFun.chi2FitTo(chSpec, Range(3.0, 7.2), Strategy(3), PrintLevel(1), Hesse(kTRUE), Save());

    RooPlot *chFrame = crySpec.frame(Title(title));
    chSpec.plotOn(chFrame, MarkerColor(kBlack), LineColor(kBlack), MarkerSize(0.5), Name("chSpec"));
    fitFun.plotOn(chFrame, LineColor(kRed), LineStyle(1), Name("fullFit"));

    // the parameter is 11 after setting width of escape peaks opened
    chiSq = chFrame->chiSquare(11);

    fitFun.plotOn(chFrame, Components(fullErg), LineColor(kOrange), LineStyle(5));
    fitFun.plotOn(chFrame, Components(firsErg), LineColor(kViolet), LineStyle(5));
    fitFun.plotOn(chFrame, Components(secdErg), LineColor(kCyan), LineStyle(5));
    fitFun.plotOn(chFrame, Components(comPdf), LineColor(kBlue), LineStyle(5));
    //double fpeak, fsigma;
    //double dpeak, dsigma;
    fpeak = fullPeak.getVal();
    dpeak = fullPeak.getError();
    fsigma = fullWidth.getVal();
    dsigma = fullWidth.getError();
    vfrFull = frFull.getVal();
    vfrFirst = frFrst.getVal();
    vfrSecond = frScnd.getVal();
    fullresolution = 100*(fsigma/fpeak);

    firstpeak = fstEsPeak.getVal();
    //firstdpeak = fstEsPeak.getError();
    firstsigma = fstWidth.getVal();
    firstdsigma = fstWidth.getError();
    firstresolution = 100*(firstsigma/firstpeak);

    secondpeak = scdEsPeak.getVal();
    //seconddpeak = scdEsPeak.getError();
    secondsigma = scdWidth.getVal();
    seconddsigma = scdWidth.getError();
    secondresolution = 100*(secondsigma/secondpeak);

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
    fspec -> Close();
    newTree->Fill();


  }
  ouptFile -> Write();
  ouptFile -> Close();
  return 0;
}

