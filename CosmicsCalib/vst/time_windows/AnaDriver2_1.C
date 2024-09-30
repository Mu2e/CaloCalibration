#define AnaDriver2_cxx
#include "AnaDriver2_1.h"
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TF1.h"
#include <TSpline.h>
#include <algorithm>
#include <TGraphErrors.h>
TH1F* hDTdiff[3][20];
TH1F* hAsym[3][30];
TH1F* hCharge[3][20];
TH2F* hDTvsQ;
TH1F* hAsymall[3];
TH2F* hAsymvsE[3];
float MPV[3][20][2]; //estratti da calib_par.dat
TH1F* MPVL_r[3][20];
TH1F* MPVR_r[3][20];
TH1F *MPV_r = new TH1F("MPV_r", "Mpv RECALIBRATED distribution", 100, 0., 50.);
TH1F *MPV_r_sig = new TH1F("MPV_r_sig", "#sigma / Mpv RECALIBRATED distribution", 100, 0., 50.);

TF1 *Gauss_fit (TH1F * histo, float start = 9999.99, float end = 9999.99);

double AnaDriver2::Loop(TString OutputFile,
		      int timeflag, float MipCut, float Chi2Cut)
{
  //timeflag==0 correction for Tmarker (*5 in ns). <-
  //timeflag==1 no correction for Tmarker
   if (fChain == 0) return -1;
//------- initialization parts --------------
   int ndisplay = 0; 
   int ldebug = 1;

   //.. FIT Template
   TFile* f0 = new TFile("splines4_run39.root");
   TGraphErrors* gt = (TGraphErrors*)f0->Get("gtempl");

   // ------ output files ----------------------
   
   cout << "Selected output file: " << OutputFile << endl;
   TFile *outFile = new TFile(OutputFile,"recreate");
   BookHistos(); // Book histos for analysis

   // Read MIP file
   std::ifstream inputFile("calib_parameter.dat");
  if (!inputFile) {
    std::cerr << " Error opening Calib file" << std::endl;}

  std::string line;
  int board,chan, sipm;
  float mpv, sigma, chi2;
  
  while (std::getline(inputFile, line)) {
    std::istringstream iss(line);
    if (iss >> board  >> chan >> sipm >> mpv >> sigma >> chi2) {
      MPV[board][chan][sipm] = mpv;}
    cout << "MP B/Ch/Sipm " <<
      board << " " << chan << " " << sipm << " " << mpv << endl;
  }

   // Loop on Reconstructed Root File
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     if(jentry%10000==0){
       std::cout << "jentry: " << jentry<< std::endl;
     }
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      int Board = -1;
      // Just a list of matrixes: 3-Board , 20-Channel, L/R
      // This should become an ordered Structure 
      // with easy reporting vs IDAQ and/or vs X-Y
      float TVAL[3][20][2];
      float CHARGE[3][20][2];
      float TDIFF[3][20];
      //
      // Matrixes initialization. Set T and V different from 0
      //
      for (int jb=0;jb<3;jb++){
	for (int jch=0; jch<20; jch++){
	  TDIFF[jb][jch]= -999.;        
	  for (int js=0; js<2; js++){
	    TVAL[jb][jch][js]=-999;
	    CHARGE[jb][jch][js] = 0;
	  }
	}
      }
      //
      // Loop on hits and define int pointer to boards
      //
      for (int jhit=0; jhit<nHits; jhit++){
	if( iBoa[jhit]==120 || iBoa[jhit]==121 )Board=0;
	if( iBoa[jhit]==128 || iBoa[jhit]==129 )Board=1;
	if( iBoa[jhit]==136 || iBoa[jhit]==137 )Board=2;	
	if( Board!= -1 ){ // Consider only the 6 boards in use
	  int Chan = iCha[jhit];
	  int Sensor = SiPM[jhit];
	  float timing=-999;
	  // Define timing from Fit or Fit + Tmarker
	  if( templChi2[jhit]< Chi2Cut){
	    if(timeflag==0) timing= Tval[jhit]+templTime[jhit];
	    if(timeflag==1) timing= templTime[jhit];
	    TVAL[Board][Chan][Sensor] = timing;
	  }
	  // Fill charge block with Max Pulse Height
	  CHARGE[Board][Chan][Sensor] = Vmax[jhit]; 
	}
      }
      //
      // Now loop on matrices and calculate difference or Asym
      //
      for (int jb=0; jb<3; jb++){
	for (int jc=0; jc<20; jc++){
	  float MipL = MPV[jb][jc][0];
	  float MipR = MPV[jb][jc][1];
	  float Eleft=0, Eright=0;
	  if(MipL >0 && MipR >0){
	    Eleft = CHARGE[jb][jc][0]/MipL*20;
	    Eright= CHARGE[jb][jc][1]/MipR*20;
	    // std::cout << "eleft: " << Eleft << " ergiht: " << Eright << std::endl;
	  }
	  float Epeak = (Eleft+Eright)/2.;
	  float Asym = (Eleft-Eright)/(2.*Epeak);
	    
	  if( CHARGE[jb][jc][0]> MipCut && CHARGE[jb][jc][0]<1500 &&
	      CHARGE[jb][jc][1]> MipCut && CHARGE[jb][jc][0]<1500 ){
	   }
	  if( Epeak>0){
	    hCharge[jb][jc]->Fill(Epeak);
	    hAsym[jb][jc]->Fill(Asym);
	    hAsymall[jb]->Fill(Asym);
	    hAsymvsE[jb]->Fill(Epeak,Asym);

	  }
	  if(Epeak > MipCut){
	    //	hDTvsQ->Fill(Epeak, TDIFF[jb][jc]);
	    MPVL_r[jb][jc]->Fill(Eleft);
	    MPVR_r[jb][jc]->Fill(Eright);
	    if(TVAL[jb][jc][1]> -999 && TVAL[jb][jc][0] > -999){
	      TDIFF[jb][jc]=TVAL[jb][jc][1]-TVAL[jb][jc][0];
	      hDTdiff[jb][jc]->Fill(TDIFF[jb][jc]);
	    }
	  }

	}
      }
   }//END JENTRY

  // MEAND AND SIGMA OF TDIFF - Prepare histograms
  TH1F *h_mu = new TH1F("h_mu", "h_mu", 25, -10., 10.);
  TH1F *h_sigma = new TH1F("h_sig", "h_sig", 50, 0., 2.);
  // Cycle on boards
  for (int i_b = 0; i_b < 3; i_b ++ ){
    // Cycle on channels
    for (int i_c = 0; i_c < 20; i_c ++){
      TH1F *cur_hist = hDTdiff[i_b][i_c];
      TF1 *fit = Gauss_fit(cur_hist);
      double *params = fit -> GetParameters();
      // The parameters of gaus are constant, mean and sigma
      if(params[1]){
        h_mu -> Fill (params[1]);
      }
      if(params[2]){
        h_sigma -> Fill (params[2]);
      }
    }
  }
  // Fit mu and sigma DT histograms
  TF1 *fit_mu = Gauss_fit(h_mu);
  TF1 *fit_sigma = Gauss_fit(h_sigma);

// Save histograms
//
  outFile->cd();
  h_mu -> Write();
  h_sigma -> Write();
  
  for (int jb=0; jb<3; jb++){
    for (int jc=0; jc<20; jc++){
	    hDTdiff[jb][jc]->Write();
    }  
  }
   
  // Return the mean of the gaussian over the DT sigmas
  double mean_sigma = fit_sigma -> GetParameters()[1];
  return mean_sigma;
}

void AnaDriver2::BookHistos()
{ // ...  book histos ... 

  for(int jb=0; jb<3; jb++){
    hAsymall[jb] = new TH1F(Form("hAsymall_%d",jb),
			    Form("hAsymall_%d",jb),
			    200,-0.5,0.5);
    hAsymvsE[jb]= new TH2F(Form("hAsymvsE_%d",jb),
			   Form("hAsymvsE_%d",jb),
			   200,0.,2000,200,-0.5,0.5);
    for(int jch=0; jch<20; jch++){    
      hDTdiff[jb][jch] = new TH1F(Form("hDTdiff_%d_%d",jb,jch),
				  Form("hDTdiff_%d_%d",jb,jch),
				  200,-10.,10.);
      hCharge[jb][jch] = new TH1F(Form("hCharge_%d_%d",jb,jch),
				  Form("hCharge_%d_%d",jb,jch),
				  200,0.,2000.);
				  
      hAsym[jb][jch] = new TH1F(Form("hAsym_%d_%d",jb,jch),
				Form("hAsym_%d_%d",jb,jch),
				200,-0.5,0.5);
      MPVR_r[jb][jch] = new TH1F(Form("MPVR_r_%d_%d",jb,jch), Form("MPVR_r_%d_%d" ,jb,jch), 100, 0., 100);
      MPVL_r[jb][jch] = new TH1F(Form("MPVL_r_%d_%d",jb,jch), Form("MPVL_r_%d_%d" ,jb,jch), 100, 0., 100); 

    }
  }
}

TF1 *Gauss_fit (TH1F * histo){
  //Just a simple funciton that does a Gauss fit over all the histogram (or over the provided range)
  // The parameters of gaus are constant, mean and sigma
  float start = histo -> GetXaxis() -> GetXmax();
  float end = histo -> GetXaxis() -> GetXmin();
  TF1 *gauss = Gauss_fit (histo, start, end);
  return gauss;
}

TF1 *Gauss_fit (TH1F * histo, float start, float end){
  //Just a simple funciton that does a Gauss fit over all the histogram (or over the provided range)
  // The parameters of gaus are constant, mean and sigma
  TF1 *gauss = new TF1("Gauss", "gaus", start, end);
  histo -> Fit(gauss, "Q");
  return gauss;
}