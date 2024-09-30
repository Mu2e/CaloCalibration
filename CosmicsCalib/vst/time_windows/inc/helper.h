#ifndef HELPER_H
#define HELPER_H

#include "langaus.h"
#include <TH2.h>
#include <TH1F.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <fstream>
#include <iostream>
#include <TGraph.h>
#include <TMath.h>
#include <vector>
#include <TString.h>
#include <string>
#include <TMultiGraph.h>
#include <TFile.h>

struct event{
     int nhits;
     std::vector<std::vector<float>> posXY;
     std::vector<float> ene;
     std::vector<float> time;
     std::vector<std::vector<int>> boachasipm;
};

void printInfo(event evt){
    for(int i = 0; i< evt.nhits; i++){
      std::cout <<" iHit: "<< i+1 << "/" << evt.nhits <<std::endl;
      std::cout << " board: " << evt.boachasipm[i][0] << " chan: " << evt.boachasipm[i][1] << " SiPM: " << evt.boachasipm[i][2] << std::endl;
      std::cout << " energy: " << evt.ene[i] << " time " << evt.time[i] << std::endl;
	std::cout << " x: " << evt.posXY[i][0] << " y: " << evt.posXY[i][1] << std::endl << std::endl; 
  }
};


int trackselection(event evt, float enecut){
  TGraph err = new TGraphErrors();
  //flag for track type: 0 bad track, 1 vertical track 
  int tracktype=0;
  float err[evt.nhits] = {9.81}; //34/sqrt(12)
  std::vector<float> x;
  std::vector<float> y;
  //vertical track find
  for(int ihit =0; ihit < evt.nhits; ihit++){
    if(evt.ene[ihit] > enecut){
      x.push_back(evt.posXY[ihit][0]);
      y.push_back(evt.posXY[ihit][1]);
    }
  }
  for(int a =0; a < x.size(); a++){
    for(int b = a; b < x.size(); b++){
      if(1){}
    }
  }
  
  return tracktype;
};

event CutNClean(event evt, float enecut) {
  int ihit = 0;
  while (ihit < evt.nhits) {
    if (evt.ene[ihit] < enecut) {
      evt.ene.erase(evt.ene.begin() + ihit);
      evt.time.erase(evt.time.begin() + ihit);
      evt.posXY.erase(evt.posXY.begin() + ihit);
      evt.boachasipm.erase(evt.boachasipm.begin() + ihit);
      evt.nhits--;
    }
    else {
      ihit++;
    }
  }
  return evt;
};

std::vector<std::vector<std::vector<TH1F*>>> bookHistos(int nBoard, int nChan, int nbin, float xmin, float xmax, TString title){
  //initialize histots
  std::vector<std::vector<std::vector<TH1F*>>> histos(nBoard, std::vector<std::vector<TH1F*>>(nChan, std::vector<TH1F*>(2, nullptr)));
  std::string histtile;
  for(int iboard = 0;  iboard < nBoard; iboard++){
    for(int ichan = 0; ichan < nChan ; ichan++){
      for(int isipm =0; isipm<2; isipm++){
	TString histtitle = Form("%s in board %d channel %d SiPM %d", title.Data(), iboard, ichan, isipm);
	histos[iboard][ichan][isipm] = new TH1F(Form("Hist_%d_%d_%d", iboard, ichan, isipm), histtitle, nbin, xmin, xmax);
      }
    }
  }
  return histos;
};

std::vector<float> fitLangaus(TH1F* hist){
  std::vector<float> fitparameters;
  TF1* mylang = new TF1("mylang" ,langaufun,0.,50000., 4);
  TF1* myland = new TF1("myland", "landau");
  TF1* mygaus = new TF1("mygaus", "gaus");
  
  mylang->SetParName(0, "#sigma");
  mylang->SetParName(1, "MPV");
  mylang->SetParName(2, "Norm");
  mylang->SetParName(3, "g #sigma");
  
  int ibin = hist->GetMaximumBin();
  float Xmean = hist->GetBinCenter(ibin);
  float Xsigma= hist->GetRMS();
  float Xming = Xmean - 0.3*Xsigma;
  float Xmaxg = Xmean + 0.8*Xsigma;
  
  hist->Fit("myland", "Q0", "");

  float par0 = myland->GetParameter(2);
  float par1 = myland->GetParameter(1);
  float par2 = myland->GetParameter(0);
  
  hist->Fit("mygaus", "Q0", "", Xming, Xmaxg);
  float par3 = mygaus->GetParameter(2);

  mylang->SetParameters(par0, par1, par2, par3);

  float Xmin = par1 - 1*par0;
  float Xmax = par1 + 20*par0;
  
  hist->Fit("mylang","Q","",Xmax,Xmin);

  float chi2ndf= mylang->GetChisquare()/mylang->GetNDF();
  float sig= mylang->GetParameter(0);
  //if the fit goes wrong WE FIGHT!! FOR FRODO untill we get a reasonable fit
  int c= 2;
  if(sig < 10){ //sig < 200 for qval
    while(sig < 10 && c < 11){
      mylang->SetParameters(par0*c, par1, par2, par3);
      hist->Fit("mylang","Q","",Xmin,Xmax);
      sig = mylang->GetParameter(0);
      c++;
    }
  }

  fitparameters.push_back(mylang->GetParameter(0));
  fitparameters.push_back(mylang->GetParameter(1));
  fitparameters.push_back(mylang->GetParameter(2));
  fitparameters.push_back(mylang->GetParameter(3));
  fitparameters.push_back(mylang->GetChisquare()/mylang->GetNDF());
  fitparameters.push_back(mylang->GetParError(0));
  fitparameters.push_back(mylang->GetParError(1));
  fitparameters.push_back(mylang->GetParError(2));
  fitparameters.push_back(mylang->GetParError(3));

  return fitparameters;
};


void scatterEch(std::vector<std::vector<float>> fitspars){
  TH2F plotL("plotL", "MPV distribution per channels", 60, 0., 60., 100, 0., 15000.);
  TH2F plotR("plotR", "MPV distribution per channels", 60, 0., 60., 100, 0., 15000.);
  TGraph *spL = new TGraph(60);
  TGraph *spR = new TGraph(60);
  int counter=0;
  TMultiGraph *mg = new TMultiGraph();
  int size =  fitspars.size();
  for(int ch=0; ch <size; ch+=2){
    spL->AddPoint(counter, fitspars[ch+1][1]);
    spR->AddPoint(counter, fitspars[ch][1]);

    //since the channel fit parameters are saved alternatively sipm 0 and sipm 1 all even entries are sipm0 all odd are sipm one for thi we us a counter thta increase avery 2 step
    counter++;
  }
  spL->SetTitle("MPV Left SiPM");
  spL->GetXaxis()->SetTitle("Channel Number");
  spL->GetYaxis()->SetTitle("MPV");

  spR->SetTitle("MPV Right SiPM");
  spR->GetXaxis()->SetTitle("Channel Number");
  spR->GetYaxis()->SetTitle("MPV");
  spL->SetMarkerStyle(8);
  spR->SetMarkerStyle(22);
  spR->SetMarkerColor(kRed);
  mg->Add(spL, "p");
  mg->Add(spR, "p");

  auto cech = new TCanvas("cech", "cech");
  cech->cd();
  mg->Draw("A");
  mg->GetYaxis()->SetRangeUser(0., 15000);
  //cech->BuildLegend();

  cech->SaveAs("scatter_MPV_vmax.root");
};



void scatterSigmasumu(std::vector<std::vector<float>> fitspars){

  TGraph *spL = new TGraph(60);
  TGraph *spR = new TGraph(60);
  int counter=0;
  TMultiGraph *mg = new TMultiGraph();
  int size =  fitspars.size();
  for(int ch=0; ch <size; ch+=2){
    spL->AddPoint(counter, fitspars[ch+1][0]/fitspars[ch+1][1]);
    spR->AddPoint(counter, fitspars[ch][0]/fitspars[ch][1]);  

    //since the channel fit parameters are saved alternatively sipm 0 and sipm 1 all even entries are sipm0 all odd are sipm one for thi we us a counter thta increase avery 2 step
    
    counter++;
  }
  spL->SetTitle("#sigma / MPV Left SiPM");
  spL->GetXaxis()->SetTitle("Channel Number");
  spL->GetYaxis()->SetTitle("MPV");

  spR->SetTitle("#sigma / MPV Right SiPM");
  spR->GetXaxis()->SetTitle("Channel Number");
  spR->GetYaxis()->SetTitle("MPV");
  spL->SetMarkerStyle(8);
  spR->SetMarkerStyle(22);
  spR->SetMarkerColor(kRed);
  mg->Add(spL, "p");
  mg->Add(spR, "p");

  auto cech = new TCanvas("cechsm", "cechsm");
  cech->cd();
  mg->Draw("A");
  mg->GetYaxis()->SetRangeUser(0., 0.2);
  //cech->BuildLegend();

  cech->SaveAs("scatter_sigma_MPV_vmax.root");
}

void drawHisto(TH1F* histo, std::string title){
  auto chisto = new TCanvas("chisto", "chisto");
  chisto->cd();
  histo->Draw();
  chisto->Update();
  std::string tit = title+".root";
  chisto->SaveAs(tit.c_str());
}

void occupancy(std::vector<std::vector<std::vector<int>>> occu){
  
  TGraph *spL = new TGraph(60);
  TGraph *spR = new TGraph(60);
  int counter=0;
  TMultiGraph *mg = new TMultiGraph();
  int chansize = occu[0].size();
  int boasize = occu.size();
  for(int boa =0; boa < boasize; boa++){
    for(int ch=0; ch <chansize; ch++){
      spL->AddPoint(counter, occu[boa][ch][1]);
      spR->AddPoint(counter, occu[boa][ch][0]);  

      //since the channel fit parameters are saved alternatively sipm 0 and sipm 1 all even entries are sipm0 all odd are sipm one for thi we us a counter thta increase avery 2 step
    
      counter++;
    }
  }
   spL->SetTitle("Occupancy Left SiPM");
  spL->GetXaxis()->SetTitle("Channel Number");
  spL->GetYaxis()->SetTitle("MPV");

  spR->SetTitle("Occupancy Right SiPM");
  spR->GetXaxis()->SetTitle("Channel Number");
  spR->GetYaxis()->SetTitle("MPV");
  spL->SetMarkerStyle(8);
  spR->SetMarkerStyle(22);
  spR->SetMarkerColor(kRed);
  mg->Add(spL, "p");
  mg->Add(spR, "p");

  auto cech = new TCanvas("cocc", "cocc");
  cech->cd();
  mg->Draw("A");
  mg->GetYaxis()->SetRangeUser(0., 4000);
  //cech->BuildLegend();

  cech->SaveAs("occupancy_VMAX.root");
}

float Asym(float eL, float eR){
  float a = (eL - eR)/(eL + eR);
  return a;
};

std::vector<std::vector<float>> find2sipm(event evt, std::vector<float> MPV){
  std::vector<std::vector<float>> EeAsym;
  int diag=0;
  
  for(int ihit = 0; ihit< evt.nhits; ihit++){
    int board = evt.boachasipm[ihit][0];
    if(board == 120 || board ==121) board=0;
    else if(board == 128 || board ==129) board=1;
    else if(board == 136 || board ==137) board=2;

    int chann = evt.boachasipm[ihit][1];
    int sipm = evt.boachasipm[ihit][2];

    //looking for 2 simpm over trashold
    for(int ihit2 =ihit; ihit2< evt.nhits; ihit2++){
      int board2 = evt.boachasipm[ihit2][0];
	
      if(board2 == 120 || board2 ==121) board2=0;
      else if(board2 == 128 || board2 ==129) board2=1;
      else if(board2 == 136 || board2 ==137) board2=2;

      int chann2 = evt.boachasipm[ihit2][1];
      int sipm2 = evt.boachasipm[ihit2][2];

      if(board == board2 && chann == chann2 && sipm != sipm2){
	if(diag >0){
	std::cout << "found a compatibility: " << std::endl
		  << "board: " << board << " = " << board2 <<std::endl
		  << "chann: " << chann << " = " << chann2 <<std::endl
		  << "sipm: " << sipm << " = " << sipm2 <<std::endl;
	}
	float mpv_id = board2*40 + chann2*2 + sipm2; //sipm 1 == right energy odd board
	float E_r=0;
	float E_l=0;
	if(MPV[mpv_id-1] > 0){
	  E_r = evt.ene[ihit]/MPV[mpv_id]*20; //MeV mpv_id==1
	  E_l = evt.ene[ihit2]/MPV[mpv_id-1]*20; //MeV mpv_id-1 ==0;
	}
	float E = (E_r + E_l)/2;
	float asym = -999;
	if(E> 0) asym=Asym(E_l,E_r);
	std::vector<float> temp;
	temp.push_back(E);
	temp.push_back(asym);
	EeAsym.push_back(temp);
      }
    }
  }
  return EeAsym;
};

std::vector<float> fitGaus(TH1F* hist){
  int diag =0;
  std::vector<float> pars;
  TF1* mygaus = new TF1("mygaus", "gaus");
  hist->Fit("mygaus", "Q", "", -1.,1.);
  pars.push_back(mygaus->GetParameter(2)); //sigma
  pars.push_back(mygaus->GetParError(2)); //sigmaerr
  if(diag >0){
    std::cout << "mu : " << pars[0]
	      << " sigma: " <<  pars[1] << std::endl;
  }
  return pars;
};

std::vector<float> loadMPV(std::string namefile){

  struct DataPoint {
    int board;
    int chan;
    int sipm;
    double mpv;
    double sigma;
    double chi2_ndf;
    
};
  
  std::ifstream inputFile(namefile);
  if (!inputFile) {
    std::cerr << "Errore nell'apertura del file!" << std::endl;
  }
  std::vector<float> MPV;
  std::string line;
  
  while (std::getline(inputFile, line)) {
    std::istringstream iss(line);
    DataPoint dp;
    if (iss >> dp.board >> dp.chan >> dp.sipm >> dp.mpv >> dp.sigma >> dp.chi2_ndf) {
      MPV.push_back(dp.mpv);
    }
  }

  inputFile.close();
  return MPV;
};

void saveHistos(std::vector<std::vector<std::vector<TH1F*>>> histos, std::string title){
  TFile *file = TFile::Open(title.c_str(), "RECREATE");
   for(int iboard =0; iboard < histos.size(); iboard++){
     for(int ichan=0; ichan < histos[0].size(); ichan++){
       for(int isipm =0; isipm < histos[0][0].size(); isipm++){
	 histos[iboard][ichan][isipm]->Write(Form("board_%i_chan_%i_sipm_%i", iboard, ichan, isipm));
       }
     }
   }
};

void recalibratedE(std::vector<event> evts){
  int nboards =3, nchans =20, nsipms =2;
  std::vector<float> MPV = loadMPV("calib_parameter.dat");
  TH1F hqval("hqval", "MPV distribution recalibrated", 50, 0., 100);
  std::vector<std::vector<std::vector<TH1F*>>> histos = bookHistos(nboards, nchans, 50, 0., 250., "energy recalibrated deposited");
  for(int ievt =0; ievt< evts.size(); ievt++){

    for(int ihit = 0; ihit< evts[ievt].nhits; ihit++){
      int board = evts[ievt].boachasipm[ihit][0];
      if(board == 120 || board ==121) board=0;
      if(board == 128 || board ==129) board=1;
      if(board == 136 || board ==137) board=2;

      int chann = evts[ievt].boachasipm[ihit][1];

      int sipm = evts[ievt].boachasipm[ihit][2];

      int id_mpv = 40*board + 2*chann + sipm;
      // std::cout << "board: " << board << " chann: " << chann << " sipm: " << sipm << " id_mpv: " << id_mpv << std::endl;
      histos[board][chann][sipm]->Fill(evts[ievt].ene[ihit]/MPV[id_mpv]*20);

    }
  }
  for(int a = 0; a < nboards; a++){
    for(int b = 0; b < nchans; b++){
      for(int c =0; c<nsipms; c++){
	std::vector<float> fitpars = fitLangaus(histos[a][b][c]);
	  hqval.Fill(fitpars[1]);
	  if(fitpars[1] <= 0){
	    std::cout << "error mpv board: " << a << " chann: " << b << " sipm: " << c << std::endl;
	  }
      }
    }
  }
  saveHistos(histos, "langaus_fit_recalib.root");
  drawHisto(&hqval, "mpv_distr_recalib");


};

#endif
