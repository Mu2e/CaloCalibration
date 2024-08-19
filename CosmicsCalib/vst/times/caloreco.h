//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jul 26 15:18:29 2024 by ROOT version 6.28/06
// from TTree tree/tree
// found on file: cosmic_calibrated_run0.root
//////////////////////////////////////////////////////////

#ifndef caloreco_h
#define caloreco_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <iostream>

#include <string>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <sstream>

#define maxNcry     500
#define maxNsipm   1000
#define maxNsample  200

// Header file for the classes stored in the TTree if any.

class caloreco {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Int_t           subrun;
   Int_t           nevt;
   Int_t           dtcID;
   Long64_t        currentDTCEventWindow;
   Int_t           nhits;
   Int_t           boardID[1000];   //[nhits]
   Int_t           linkID[1000];   //[nhits]
   Int_t           chanID[1000];   //[nhits]
   Int_t           errflag[1000];   //[nhits]
   Int_t           fff[1000];   //[nhits]
   Int_t           timetot[1000];   //[nhits]
   Int_t           ewhit[1000];   //[nhits]
   Int_t           peakpos[1000];   //[nhits]
   Int_t           peakval[1000];   //[nhits]
   Int_t           nofsamples[1000];   //[nhits]
   Int_t           firstsample[1000];   //[nhits]
   Int_t           nsamples;
   Int_t           ADC[6300];   //[nsamples]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_subrun;   //!
   TBranch        *b_nevt;   //!
   TBranch        *b_dtcID;   //!
   TBranch        *b_currentDTCEventWindow;   //!
   TBranch        *b_nhits;   //!
   TBranch        *b_boardID;   //!
   TBranch        *b_linkID;   //!
   TBranch        *b_chanID;   //!
   TBranch        *b_errflag;   //!
   TBranch        *b_fff;   //!
   TBranch        *b_timetot;   //!
   TBranch        *b_ewhit;   //!
   TBranch        *b_peakpos;   //!
   TBranch        *b_peakval;   //!
   TBranch        *b_nofsamples;   //!
   TBranch        *b_firstsample;   //!
   TBranch        *b_nsamples;   //!
   TBranch        *b_ADC;   //!

   //caloreco(TTree *tree=0);
   caloreco(TString fName="",TTree *tree=0);

   virtual ~caloreco();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   //virtual void     Loop();
   virtual void     Loop(TString OutputFile="", int evflag=0, float xstart=-80., float xend=20.);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   // Additional user functions
   virtual void  BookOutput(TTree *sidet);
   virtual int   GetValues (int jChan, int evflag);
   virtual void  getTemplateFit(TGraphErrors* gt, float xstart, float xend);
   //virtual int   Get_iCry (int jRow, int jCol);
   //virtual int   Get_iSiPM (int jRow, int jCol, int jSiPM);

   //
   // Output nutple
   //

   TString         OutputFile;

   Int_t           nrun;
   Int_t           nsubrun;
   Int_t           evnum;
   Int_t           nHits;
   Int_t           _isLaser;
   Int_t           iDAQ[maxNsipm];
   Int_t           iBoa[maxNsipm];
   Int_t           iCha[maxNsipm];
   Int_t           iHit[maxNsipm];
   Int_t           iRow[maxNsipm];
   Int_t           iCol[maxNsipm];
   Int_t           SiPM[maxNsipm];
   Int_t           iCry[maxNsipm];
   Double_t        Xval[maxNsipm];
   Double_t        Yval[maxNsipm];
   Int_t           iMax[maxNsipm];
   Double_t        Qval[maxNsipm];
   Double_t        Tval[maxNsipm];
   Double_t        Vmax[maxNsipm];
   Int_t           nSamples[maxNsipm];
   Double_t        wave[maxNsipm][maxNsample];
   Double_t        tWave[maxNsipm][maxNsample];
   Double_t        bline[maxNsipm];
   Double_t        templTime[maxNsipm];
   Double_t        templChi2[maxNsipm];
   Double_t        templFit[maxNsipm][3];
   Double_t        templErr[maxNsipm][3];

   /*
   int         Itmp;
   int         ItmpLaser;
   float       Qtmp, Ttmp, Btmp, Vtmp, Xtmp, Ytmp;
   int         Stmp;
   Double_t    wavetmp[maxNsample];
   Double_t    twavetmp[maxNsample];
   int         BoaTmp, ChaTmp, DaqTmp, HitTmp, SipmTmp, CryTmp;
   int         DiskTmp, PhiTmp, CrateTmp, BoaIdx, MzbIdx, ConIdx, FEETmp, cryTmp, rouTmp;
   float       RowTmp, ColTmp, xTmp, yTmp;
   int         daqID[Nbrds][Nchan], rowID[Nbrds][Nchan], colID[Nbrds][Nchan];
   int         sipmID[Nbrds][Nchan], cryID[Nbrds][Nchan], boardIdx[Nbrds][Nchan];
   int         channIdx[Nbrds][Nchan], HitNum[Nbrds][Nchan];
   float       xVal[Nbrds][Nchan], yVal[Nbrds][Nchan];

   float       Time, Chi2;
   float       fitPar[3]={0.}, fitErr[3]={0.}, fitTmea=0., fitChi2=0.;
   */
};

#endif

#ifdef caloreco_cxx
//caloreco::caloreco(TTree *tree) : fChain(0) 
caloreco::caloreco(TString fileList, TTree *tree): fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
/*
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("cosmic_calibrated_run0.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("cosmic_calibrated_run0.root");
      }
      f->GetObject("tree",tree);

   }
*/
  if (tree == 0) {
     // Get list of files and chain ntuples
     TChain *chain = new TChain("tree","");
     ifstream file(fileList);
     string fileName;
     while( getline(file,fileName) ){
       fileName.append("/tree");
       chain->Add(&fileName[0]);  // string2char by assigning the 1st character 
       cout << fileName << endl;  // of string to a pointer to char
     }
     tree=chain;
  }
   Init(tree);
}

caloreco::~caloreco()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t caloreco::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t caloreco::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void caloreco::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("subrun", &subrun, &b_subrun);
   fChain->SetBranchAddress("nevt", &nevt, &b_nevt);
   fChain->SetBranchAddress("dtcID", &dtcID, &b_dtcID);
   fChain->SetBranchAddress("currentDTCEventWindow", &currentDTCEventWindow, &b_currentDTCEventWindow);
   fChain->SetBranchAddress("nhits", &nhits, &b_nhits);
   fChain->SetBranchAddress("boardID", boardID, &b_boardID);
   fChain->SetBranchAddress("linkID", linkID, &b_linkID);
   fChain->SetBranchAddress("chanID", chanID, &b_chanID);
   fChain->SetBranchAddress("errflag", errflag, &b_errflag);
   fChain->SetBranchAddress("fff", fff, &b_fff);
   fChain->SetBranchAddress("timetot", timetot, &b_timetot);
   fChain->SetBranchAddress("ewhit", ewhit, &b_ewhit);
   fChain->SetBranchAddress("peakpos", peakpos, &b_peakpos);
   fChain->SetBranchAddress("peakval", peakval, &b_peakval);
   fChain->SetBranchAddress("nofsamples", nofsamples, &b_nofsamples);
   fChain->SetBranchAddress("firstsample", firstsample, &b_firstsample);
   fChain->SetBranchAddress("nsamples", &nsamples, &b_nsamples);
   fChain->SetBranchAddress("ADC", ADC, &b_ADC);
   Notify();
}

Bool_t caloreco::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void caloreco::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t caloreco::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef caloreco_cxx
