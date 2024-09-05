#define caloreco_cxx
#include "caloreco.h"
#include <fstream>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TF1.h"
#include <TSpline.h>
#include <algorithm>
TSpline3 *sp3;

//=============================================================================
// Fit function for template
//=============================================================================

Double_t spline_fit(Double_t *x, Double_t *par)
{
  double f1;
  f1=par[0]*sp3->Eval(x[0]-par[1])+par[2];
  return f1;  
}

//=============================================================================
// User variables
//=============================================================================

int const Nbrds = 160;
int const Nchan =  20;

int         Itmp, Stmp;
float       Qtmp, Ttmp, Btmp, Vtmp, Xtmp, Ytmp;
Double_t    wavetmp[maxNsample];
Double_t    twavetmp[maxNsample];
int         BoaTmp, ChaTmp, DaqTmp, HitTmp, SipmTmp, CryTmp;
int         DiskTmp, PhiTmp, CrateTmp, BoaIdx, MzbIdx, ConIdx;
int         FEETmp, cryTmp, rouTmp, RowTmp, ColTmp;
float       xTmp, yTmp;
int         daqID[Nbrds][Nchan]={0}, rowID[Nbrds][Nchan]={0};
int         colID[Nbrds][Nchan]={0}, sipmID[Nbrds][Nchan]={0};
int         cryID[Nbrds][Nchan]={0}, boardIdx[Nbrds][Nchan]={0};
int         channIdx[Nbrds][Nchan]={0}, HitNum[Nbrds][Nchan]={0};
float       xOffl[Nbrds][Nchan]={0.}, yOffl[Nbrds][Nchan]={0.};

float       Time, Chi2;
float       fitPar[3]={0.}, fitErr[3]={0.}, fitTmea=0., fitChi2=0.;

//=============================================================================
// Loop function
//=============================================================================

//void caloreco::Loop()
void caloreco::Loop(TString OutputFile, int evflag, float xstart, float xend)
{
  // evflag: 0 for cosmics, 1 for laser runs
  // Template fit interval in ns: [peak+xstart,peak+xend] w.r.t. wave peak
  //
//   In a ROOT session, you can do:
//      root> .L caloreco.C
//      root> caloreco t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   //**************************************************************************
   // Settings for cosmics and laser runs
   //**************************************************************************

   int skipfit;
   float Qcut = 0;

   if( evflag == 0 ) {           // cosmics
     Qcut    = 2; // moved from 10 to 2;
     skipfit = 0;
   } 
   else if( evflag==1 ){         // laser runs
     Qcut = -100.;
     skipfit= 1;
   }
   else {
     cout << "Unknown evflag " << evflag << ": Fatal error" << endl;
     return;
   }

   //**************************************************************************
   // Input files: channel maps + template to fit timing
   //**************************************************************************

   ifstream dmap;
   dmap.open("dmap_bruno.dat");
   //if(dmap.eof()) break;

   string type;
   for ( int Iloop=0; Iloop<1349; Iloop++ ){
     dmap >> RowTmp >> ColTmp >> DiskTmp >> PhiTmp >> CrateTmp >> BoaTmp >> SipmTmp
          >> BoaIdx >> MzbIdx >> ConIdx >> ChaTmp >> FEETmp >> xTmp >> yTmp >> cryTmp >> rouTmp >> type;
     //cout << RowTmp << " " << ColTmp << " " << ChaTmp << " " << cryTmp << " " << xTmp << " " << rouTmp << " " << type << endl;
     boardIdx[BoaIdx][ChaTmp] = BoaIdx;
     channIdx[BoaIdx][ChaTmp] = ChaTmp;
     daqID[BoaIdx][ChaTmp]  = FEETmp;
     rowID[BoaIdx][ChaTmp]  = RowTmp;
     colID[BoaIdx][ChaTmp]  = ColTmp;
     sipmID[BoaIdx][ChaTmp] = SipmTmp;
     xOffl[BoaIdx][ChaTmp]  = xTmp;
     yOffl[BoaIdx][ChaTmp]  = yTmp;
     cryID[BoaIdx][ChaTmp]  = cryTmp;
   }
   dmap.close();

   TFile* f0 = new TFile("splines4_run39.root");
   TGraphErrors* gt = (TGraphErrors*)f0->Get("gtempl");

   //**************************************************************************
   // Output file
   //**************************************************************************

   cout << "Selected output file: " << OutputFile << endl;
   TFile *outFile = new TFile(OutputFile,"recreate");
   TTree *sidet = new TTree("sidet","sidet");
   BookOutput(sidet); // Book histos and ntuple

   //**************************************************************************
   // Loop on entries
   //**************************************************************************

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      //cout << "Event: " << nevt << endl;

     if( jentry%1000==0 ) cout << "Number of processed events: " << jentry << endl;
     nrun    = run;
     nsubrun = subrun; 
     evnum   = nevt;
     _isLaser=0;
     /*
     if( b1_charge7>2000 ){
       _isLaser=1;
       skipfit==1;
     }
     */
     int hitTot = 0; 

     // Zeroing number of hits per channel
     for( int jBoard=0; jBoard<Nbrds; jBoard++ ){
       for( int jChan=0; jChan<Nchan; jChan++ ){
	 HitNum[jBoard][jChan] = 0;  
       }
     }

     for( int jHit=0; jHit<nhits; jHit++ ){
              
       int Status = GetValues(jHit, evflag);
       if( Status!=0 ) cout << "GetValues error: " << Status << " " <<
			 ". Skipping event/DAQid: " << nevt << " " << DaqTmp << endl;

       if( Qtmp>Qcut ){                 ///// Check Qcut value!!!
	 iDAQ[hitTot] = DaqTmp;
	 iBoa[hitTot] = BoaTmp;
	 iCha[hitTot] = ChaTmp;
	 iHit[hitTot] = HitTmp;
	 iRow[hitTot] = RowTmp;
	 iCol[hitTot] = ColTmp;
	 SiPM[hitTot] = SipmTmp;
	 iCry[hitTot] = CryTmp;
	 Xval[hitTot] = Xtmp;
	 Yval[hitTot] = Ytmp;
	 Qval[hitTot] = Qtmp;
	 Tval[hitTot] = Ttmp;
	 Vmax[hitTot] = Vtmp;
	 iMax[hitTot] = Itmp;
	 nSamples[hitTot] = Stmp;
	 std::copy( wavetmp, wavetmp+maxNsample, wave[hitTot]);
	 std::copy(twavetmp,twavetmp+maxNsample,tWave[hitTot]);

	 /*
	 // Check wave 
	 cout << "Nsamples/iMax/Vmax " << Stmp << " " << Itmp << " " << Vtmp << endl;
	 for( int ii=0; ii<20; ii++ ){
	   cout << ii << " " << wave[hitTot][ii] << endl;
	 }
	 sleep(5);
	 */

	 bline[hitTot] = Btmp;

	 // Timing from template fit
	 if( skipfit==0){
	   //cout << "Fitting " << DaqTmp << endl;
	   getTemplateFit(gt,xstart,xend);
	 }
	 templTime[hitTot] = fitTmea;
	 templChi2[hitTot] = fitChi2;
	 std::copy(fitPar,fitPar+3,templFit[hitTot]);
	 std::copy(fitErr,fitErr+3,templErr[hitTot]);
	 
	 hitTot++;
	 if( hitTot>maxNsipm ){
	   cout << "ERROR: too many hits " << hitTot << ". You need to modify the size in the include file" << endl;
	 }
       }
     }

     nHits = hitTot;

     sidet->Fill();
      
   } // end loop on events

   cout << "Closing output file" << endl;

   outFile->cd();
   sidet->Write();
   outFile->Close();

}

//=============================================================================
// Ntuple booking
//=============================================================================
void caloreco::BookOutput(TTree *sidet)
{
  // Create ROOT tree output

  sidet->SetAutoSave(1000);

  sidet->Branch("nrun"      ,&nrun,      "nrun/I");
  sidet->Branch("nsubrun"   ,&nsubrun,   "nsubrun/I");
  sidet->Branch("evnum"     ,&evnum,     "evnum/I");
  sidet->Branch("nHits"     ,&nHits,     "nHits/I");
  sidet->Branch("isLaser"   ,&_isLaser,  "isLaser/I");       // number of SiPM seeing a laser pulse   
  sidet->Branch("iDAQ"      ,&iDAQ,      "iDAQ[nHits]/I");
  sidet->Branch("iBoa"      ,&iBoa,      "iBoa[nHits]/I");
  sidet->Branch("iCha"      ,&iCha,      "iCha[nHits]/I");
  sidet->Branch("iHit"      ,&iHit,      "iHit[nHits]/I");
  sidet->Branch("iRow"      ,&iRow,      "iRow[nHits]/I");
  sidet->Branch("iCol"      ,&iCol,      "iCol[nHits]/I");
  sidet->Branch("SiPM"      ,&SiPM,      "SiPM[nHits]/I");
  sidet->Branch("iCry"      ,&iCry,      "iCry[nHits]/I");
  sidet->Branch("Xval"      ,&Xval,      "Xval[nHits]/D");
  sidet->Branch("Yval"      ,&Yval,      "Yval[nHits]/D");
  sidet->Branch("iMax"      ,&iMax,      "iMax[nHits]/I");
  sidet->Branch("Qval"      ,&Qval,      "Qval[nHits]/D");
  sidet->Branch("Tval"      ,&Tval,      "Tval[nHits]/D");
  sidet->Branch("Vmax"      ,&Vmax,      "Vmax[nHits]/D");
  sidet->Branch("nSamples"  ,&nSamples,  "nSamples[nHits]/I");
  sidet->Branch("wave"      ,&wave,      "wave[nHits][200]/D");
  sidet->Branch("tWave"     ,&tWave,     "tWave[nHits][200]/D");
  sidet->Branch("bline"     ,&bline,     "bline[nHits]/D");
  sidet->Branch("templTime" ,&templTime, "templTime[nHits]/D");
  sidet->Branch("templChi2" ,&templChi2, "templChi2[nHits]/D");
  sidet->Branch("templFit"  ,&templFit,  "templFit[nHits][3]/D");
  sidet->Branch("templErr"  ,&templErr,  "templErr[nHits][3]/D");

}

//=============================================================================
// Get values for a given board/channel
//=============================================================================
int caloreco::GetValues(int jHit, int evflag)
{
  int ErrFlag = 0;

  int   nVal;
  float Time, bSum, Qsum;
  //int   QbinMin, QbinMax;
  float const camp = 1.;
  float dt = 1/camp;
  //float const Qpeak_min      =  60.; //  -60 ns from Tave
  //float const Qpeak_max      = 190.; // +190 ns from Tave

  BoaTmp  = boardID[jHit];
  ChaTmp  = chanID[jHit];
  // Check data consistency
  if( BoaTmp!=boardIdx[BoaTmp][ChaTmp] ){
    cout << "ERROR: BoardIdx inconsistency " << BoaTmp << " " << boardIdx[BoaTmp][ChaTmp] << endl;
    exit(0);
  }
  if( ChaTmp!=channIdx[BoaTmp][ChaTmp] ){
    cout << "ERROR: ChannIdx inconsistency " << ChaTmp << " " << channIdx[BoaTmp][ChaTmp] << endl;
    exit(0);
  }

  HitNum[BoaTmp][ChaTmp]++;
  HitTmp = HitNum[BoaTmp][ChaTmp];

  DaqTmp  = daqID[BoaTmp][ChaTmp];
  RowTmp  = rowID[BoaTmp][ChaTmp];
  ColTmp  = colID[BoaTmp][ChaTmp];
  SipmTmp = sipmID[BoaTmp][ChaTmp];
  CryTmp  = cryID[BoaTmp][ChaTmp];
  Xtmp  = xOffl[BoaTmp][ChaTmp];
  Ytmp  = yOffl[BoaTmp][ChaTmp];
  Qtmp = -999.;
  Btmp = 2048.;

  Ttmp = timetot[jHit]*5.;  // Time distance from event marker in ns
  //Itmp = peakpos[jHit]-1;   // WF index corresponding to the max. amplitude
  Stmp = nofsamples[jHit];
  if( Stmp>maxNsample ){
    cout << "ERROR: wave out of boundaries " << Stmp << 
      ".. reducing size to " << maxNsample << endl;
    Stmp = maxNsample;
  }

  int jStart = firstsample[jHit];
  int jEnd   = jStart + Stmp;
  
  // Get baseline
  
  nVal = 0;                 // Number of samples for baseline evaluation
  bSum = 0.;

  for(int iSample=jStart; iSample<jStart+5; iSample++ ){
    nVal++;
    bSum = bSum + ADC[iSample];
  }

  if( nVal!=0 ){
    Btmp = bSum / nVal;
  } 
  else{
    cout << "Run/Nev/DAQid " << run << " " << nevt << " " << DaqTmp << 
      ": no samples for baseline evaluation" << endl;
    ErrFlag = 2;
  }

  //
  //  Fill wave/twave & get charge, with baseline subtracted event-by-event
  //

  if( ErrFlag==0 ){

    Qsum  = 0.;
    float maxWave = 0.;
    int myIdx = 0;

    std::fill(wavetmp, wavetmp+maxNsample, 0.);
    std::fill(twavetmp, twavetmp+maxNsample, 0.);
    
    for(int iSample=jStart; iSample<jEnd; iSample++ ){
      twavetmp[myIdx] = myIdx*5.;
      wavetmp [myIdx] = ADC[iSample]-Btmp;
      if ( wavetmp[myIdx]>maxWave ){
	Itmp = myIdx;             // WF index corresponding to the max. amplitude
	maxWave = wavetmp[myIdx];
      }
      if ( iSample>jStart+5 ){    // Q integration starts from 60 ns before peak
 	Qsum = Qsum + wavetmp[myIdx];
      }
      myIdx++;
    }
    Qtmp = Qsum; //*dt/50;
    Vtmp = wavetmp[Itmp];
  }

  return ErrFlag;
}

//=============================================================================
// Get time from template fit to waveform
//=============================================================================
void caloreco::getTemplateFit(TGraphErrors* gt, float xstart, float xend)
{
  
  TGraphErrors *gwf;
  TF1          *fitf;

  double ex[maxNsample] = {0.05}; //X-error not useful .. just for display graphs
  double ey[maxNsample];
  std::fill_n(ey,maxNsample,3.6);
  int    t_peak = 0 ;
  double v_max  = 0.;
  double CF     = 0.2;

  fitTmea = 0.;
  fitChi2 = 0.;
  std::fill_n(fitPar,3,0.);
  std::fill_n(fitErr,3,0.);

  gwf = new TGraphErrors(160,twavetmp,wavetmp,ex,ey);
  v_max  = wavetmp[Itmp];
  t_peak = twavetmp[Itmp];
  fitf = new TF1("fitf",spline_fit,t_peak+xstart,t_peak+xend,3);
  sp3 = new TSpline3("sp3", gt,0,0.,800.);

  // Three parameter fit ( 0: norm, 1:t0, 2:baseline )
  //Set start parameters and limits
  fitf->SetParameter(0,v_max);
  fitf->SetParLimits(0,30.,1.2*v_max);
  float tstart = t_peak;
  fitf->SetParameter(1,tstart);
  fitf->SetParLimits(1,tstart-abs(xstart),tstart+abs(xend));
  fitf->SetParameter(2,  0.);
  fitf->SetParLimits(2,-10,10);
      
  gwf->Fit("fitf","RBOQ"); 
  //gwf->Fit("fitf","RBOQ");
  fitTmea   = fitf->GetX(v_max*CF);
  fitPar[0] = fitf->GetParameter(0);
  fitPar[1] = fitf->GetParameter(1);
  fitPar[2] = fitf->GetParameter(2);
  fitErr[0] = fitf->GetParError(0);
  fitErr[1] = fitf->GetParError(1);
  fitErr[2] = fitf->GetParError(2);
  //cout << "Fit: " << fitPar[1] << " ";
  if( fitf->GetNDF()!=0 ) fitChi2 = fitf->GetChisquare()/fitf->GetNDF();
  //cout << fitChi2<<endl;
  
  delete gwf;
  delete fitf;
  delete sp3;

}

