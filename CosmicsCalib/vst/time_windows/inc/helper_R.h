#ifndef HELPER_R_H
#define HELPER_R_H

#include "crystalR.h"
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

void scatterEch_R(std::vector<std::vector<float>> fitspars){
  TH2F plotL("plotL", "MPV distribution per channels", 60, 0., 60., 100, 0., 15000.);
  TH2F plotR("plotR", "MPV distribution per channels", 60, 0., 60., 100, 0., 15000.);
  TGraph *spL = new TGraph(60);
  TGraph *spR = new TGraph(60);
  int counter=0;
  TMultiGraph *mg = new TMultiGraph();
  int size =  fitspars.size();
  for(int ch=0; ch <size; ch+=2){
    spL->AddPoint(chaR[counter], fitspars[ch+1][1]);
    spR->AddPoint(chaR[counter], fitspars[ch][1]);

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

  auto cech = new TCanvas("cMPV_R", "cMPV_R");
  cech->cd();
  mg->Draw("A");
  mg->GetYaxis()->SetRangeUser(0., 15000);
  //cech->BuildLegend();

  cech->SaveAs("scatter_MPV_qval_R.root");
};



void scatterSigmasumu_R(std::vector<std::vector<float>> fitspars){

  TGraph *spL = new TGraph(60);
  TGraph *spR = new TGraph(60);
  int counter=0;
  TMultiGraph *mg = new TMultiGraph();
  int size =  fitspars.size();
  for(int ch=0; ch <size; ch+=2){
    spL->AddPoint(chaR[counter], fitspars[ch+1][0]/fitspars[ch+1][1]);
    spR->AddPoint(chaR[counter], fitspars[ch][0]/fitspars[ch][1]);  

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

  auto cech = new TCanvas("cSMPV_R", "cSMPV_R");
  cech->cd();
  mg->Draw("A");
  mg->GetYaxis()->SetRangeUser(0., 0.2);
  //cech->BuildLegend();

  cech->SaveAs("scatter_sigma_MPV_qval_R.root");
}

void occupancy_R(std::vector<std::vector<std::vector<int>>> occu){
  
  TGraph *spL = new TGraph(60);
  TGraph *spR = new TGraph(60);
  int counter=0;
  TMultiGraph *mg = new TMultiGraph();
  int chansize = occu[0].size();
  int boasize = occu.size();
  for(int boa =0; boa < boasize; boa++){
    for(int ch=0; ch <chansize; ch++){
      spL->AddPoint(chaR[counter], occu[boa][ch][1]);
      spR->AddPoint(chaR[counter], occu[boa][ch][0]);  

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

  auto cech = new TCanvas("cOcc_R", "cOcc_R");
  cech->cd();
  mg->Draw("A");
  mg->GetYaxis()->SetRangeUser(0., 4000);
  //cech->BuildLegend();

  cech->SaveAs("occupancy_R.root");
}

#endif
