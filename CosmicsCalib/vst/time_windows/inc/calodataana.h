//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jul 31 09:35:18 2024 by ROOT version 6.30/04
// from TTree sidet/sidet
// found on file: /exp/mu2e/data/users/pfedeli/CaloData/roottople/cosmic_calibrated_run0_reco.root
//////////////////////////////////////////////////////////

#ifndef calodataana_h
#define calodataana_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <fstream>
#include <iostream>  

// Header file for the classes stored in the TTree if any.

class calodataana {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           nrun;
   Int_t           nsubrun;
   Int_t           evnum;
   Int_t           nHits;
   Int_t           isLaser;
   Int_t           iDAQ[150];   //[nHits]
   Int_t           iBoa[150];   //[nHits]
   Int_t           iCha[150];   //[nHits]
   Int_t           iHit[150];   //[nHits]
   Int_t           iRow[150];   //[nHits]
   Int_t           iCol[150];   //[nHits]
   Int_t           SiPM[150];   //[nHits]
   Int_t           iCry[150];   //[nHits]
   Double_t        Xval[150];   //[nHits]
   Double_t        Yval[150];   //[nHits]
   Int_t           iMax[150];   //[nHits] //bin max dell'onda
   Double_t        Qval[150];   //[nHits]
   Double_t        Tval[150];   //[nHits]
   Double_t        Vmax[150];   //[nHits]
   Int_t           nSamples[150];   //[nHits]
   Double_t        wave[150][300];   //[nHits] //n entries in un bin
   Double_t        tWave[150][300];   //[nHits]
   Double_t        bline[150];   //[nHits] //baseline media primi 5 bin onda
   Double_t        templTime[150];   //[nHits]
   Double_t        templChi2[150];   //[nHits]
   Double_t        templFit[150][30];   //[nHits]
   Double_t        templErr[150][30];   //[nHits]

   // List of branches
   TBranch        *b_nrun;   //!
   TBranch        *b_nsubrun;   //!
   TBranch        *b_evnum;   //!
   TBranch        *b_nHits;   //!
   TBranch        *b_isLaser;   //!
   TBranch        *b_iDAQ;   //!
   TBranch        *b_iBoa;   //!
   TBranch        *b_iCha;   //!
   TBranch        *b_iHit;   //!
   TBranch        *b_iRow;   //!
   TBranch        *b_iCol;   //!
   TBranch        *b_SiPM;   //!
   TBranch        *b_iCry;   //!
   TBranch        *b_Xval;   //!
   TBranch        *b_Yval;   //!
   TBranch        *b_iMax;   //!
   TBranch        *b_Qval;   //!
   TBranch        *b_Tval;   //!
   TBranch        *b_Vmax;   //!
   TBranch        *b_nSamples;   //!
   TBranch        *b_wave;   //!
   TBranch        *b_tWave;   //!
   TBranch        *b_bline;   //!
   TBranch        *b_templTime;   //!
   TBranch        *b_templChi2;   //!
   TBranch        *b_templFit;   //!
   TBranch        *b_templErr;   //!

  calodataana(TString fName="", TTree *tree=0);
   virtual ~calodataana();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(const int nev);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef calodataana_cxx
calodataana::calodataana(TString fileList, TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  /*
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/exp/mu2e/data/users/pfedeli/CaloData/roottople/cosmic_calibrated_run0_reco.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/exp/mu2e/data/users/pfedeli/CaloData/roottople/cosmic_calibrated_run0_reco.root");
      }
      f->GetObject("sidet",tree);
   }
  */
  if (tree == 0) {
     // Get list of files and chain ntuples
     TChain *chain = new TChain("sidet","");
     ifstream file(fileList);
     string fileName;
     while( getline(file,fileName) ){
       fileName.append("/sidet");
       chain->Add(&fileName[0]);  // string2char by assigning the 1st character 
       std::cout << fileName << std::endl;  // of string to a pointer to char
     }
     tree=chain;
  }
   Init(tree);
}


  
calodataana::~calodataana()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t calodataana::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t calodataana::LoadTree(Long64_t entry)
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

void calodataana::Init(TTree *tree)
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

   fChain->SetBranchAddress("nrun", &nrun, &b_nrun);
   fChain->SetBranchAddress("nsubrun", &nsubrun, &b_nsubrun);
   fChain->SetBranchAddress("evnum", &evnum, &b_evnum);
   fChain->SetBranchAddress("nHits", &nHits, &b_nHits);
   fChain->SetBranchAddress("isLaser", &isLaser, &b_isLaser);
   fChain->SetBranchAddress("iDAQ", iDAQ, &b_iDAQ);
   fChain->SetBranchAddress("iBoa", iBoa, &b_iBoa);
   fChain->SetBranchAddress("iCha", iCha, &b_iCha);
   fChain->SetBranchAddress("iHit", iHit, &b_iHit);
   fChain->SetBranchAddress("iRow", iRow, &b_iRow);
   fChain->SetBranchAddress("iCol", iCol, &b_iCol);
   fChain->SetBranchAddress("SiPM", SiPM, &b_SiPM);
   fChain->SetBranchAddress("iCry", iCry, &b_iCry);
   fChain->SetBranchAddress("Xval", Xval, &b_Xval);
   fChain->SetBranchAddress("Yval", Yval, &b_Yval);
   fChain->SetBranchAddress("iMax", iMax, &b_iMax);
   fChain->SetBranchAddress("Qval", Qval, &b_Qval);
   fChain->SetBranchAddress("Tval", Tval, &b_Tval);
   fChain->SetBranchAddress("Vmax", Vmax, &b_Vmax);
   fChain->SetBranchAddress("nSamples", nSamples, &b_nSamples);
   fChain->SetBranchAddress("wave", wave, &b_wave);
   fChain->SetBranchAddress("tWave", tWave, &b_tWave);
   fChain->SetBranchAddress("bline", bline, &b_bline);
   fChain->SetBranchAddress("templTime", templTime, &b_templTime);
   fChain->SetBranchAddress("templChi2", templChi2, &b_templChi2);
   fChain->SetBranchAddress("templFit", templFit, &b_templFit);
   fChain->SetBranchAddress("templErr", templErr, &b_templErr);
   Notify();
}

Bool_t calodataana::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void calodataana::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t calodataana::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef calodataana_cxx
