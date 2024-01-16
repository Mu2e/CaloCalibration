// =============================================================================
// extrapolate energy, position and time for offline calorimeter calibration
// 12-Dec-2023 - P.Fedeli & S.Giovannella
// 
// Starting code for retrieving calo information: 
// CaloEnergy & caloT0alig module from P.Fedeli & S.Giovannella
//
// =============================================================================

//includes

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Provenance.h"
#include "canvas/Utilities/InputTag.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "fhiclcpp/types/Sequence.h"

#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CalorimeterGeom/inc/CrystalMapper.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/DataProducts/inc/CaloConst.hh"

//C++ includes
#include <string>
#include <iostream>
#include <vector>
#include <fstream>

//ROOT includes
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TFitResult.h"

namespace mu2e{
  
  class CaloCosmicEnergy : public art::EDAnalyzer {

  public:

    explicit CaloCosmicEnergy(fhicl::ParameterSet const& pset);
    virtual ~CaloCosmicEnergy(){};
    virtual void beginJob();  
  
    void analyze(art::Event const& event) override;

    virtual void endJob();

    static double langaufun(Double_t *x, Double_t *par);

    float findpath(float m,float q, float x, float y);

  private:

    //calo parameters
    static constexpr int nSiPMs = CaloConst::_nSiPMPerCrystal;
    static constexpr int nCrystals = CaloConst::_nCrystalPerDisk;
    static constexpr int nDisks = CaloConst::_nDisk;
    static constexpr int nROchan = nDisks*nCrystals*nSiPMs;
    
    //fcl parameters
    int                   _diagLevel;
    int                   CutNCryHit; //>
    float                 CutEnergyDep;// >
    float                 CutChi2Norm;// <
  
    art::InputTag         _caloClusterTag;
    art::InputTag         _caloHitTag;
    TH1F *hSiPM[nDisks*nCrystals][nSiPMs];
    TH1F *hSiPMv[nDisks*nCrystals][nSiPMs];
    TH1F *hSiPMd[nDisks*nCrystals][nSiPMs];
    TH1F *hSiPMfp[nDisks*nCrystals][nSiPMs];
    TH1F *MPV,*MPVv,*MPVd,*MPVfp;


    std::ofstream              _outputfile;

    //helper variables
    int ncry =0; // need to count how many crystals pass the selections
    float Dx = 1400.;
    float Dy = 1400.; // need to calculate the spread alond x and y axis
    TF1* myPoly = new TF1("myPoly", "pol1", -600., 600.); //polinomial to fit the track
    TF1* mylang = new TF1("mylang" ,langaufun,15.,50., 4);//langaus function


    //we create arrays where we store informations of crystals that pass the selections
 
    float PosX[nCrystals];
    float PosY[nCrystals];
    float ErrPos[nCrystals];
    int IDs[nCrystals];
    float path[nCrystals];
    int whichHit[nCrystals];

  };//end class

  CaloCosmicEnergy::CaloCosmicEnergy(fhicl::ParameterSet const& pset) : 
    art::EDAnalyzer(pset), 
    _diagLevel(pset.get<int>("diagLevel", 0)),
    CutNCryHit(pset.get<int>("CutNCryHit", 3)),
    CutEnergyDep(pset.get<float>("CutEnergyDep", 10.)),
    CutChi2Norm(pset.get<float>("CutChi2Norm", 2.5)),
    _caloClusterTag(pset.get<art::InputTag>("caloClusterCollection")),
    _caloHitTag(pset.get<art::InputTag>("caloHitCollection")){}

  void CaloCosmicEnergy::beginJob(){
    if (_diagLevel > 0 ) std::cout<< "CaloCosmicEnergy: Entering beginJob" << std::endl;
   
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory tfdir = tfs->mkdir("All tracks");
    art::TFileDirectory tfdirv = tfs->mkdir("Vertical tracks");
    art::TFileDirectory tfdird = tfs->mkdir("Diagonal tracks");
    art::TFileDirectory tfdirfp = tfs->mkdir("All normalized tracks");
    art::TFileDirectory tfdir_res = tfs->mkdir("MPV distributions");

    for(int iCry=0; iCry < nDisks*nCrystals; iCry++){
      for(int iSiPM=0; iSiPM < nSiPMs; iSiPM++){
	hSiPM[iCry][iSiPM] = tfdir.make<TH1F>(Form("Crystal_%i_SiPM_%i",iCry, iSiPM),Form("Energy deposited in SiPM %i crystal %d",iSiPM, iCry),70, 15., 50.);
	hSiPMv[iCry][iSiPM] = tfdirv.make<TH1F>(Form("Crystal_%i_SiPM_%i_v",iCry, iSiPM),Form("Energy deposited in SiPM %i crystal %d only vertical tracks",iSiPM, iCry),70, 15., 50.);
	hSiPMd[iCry][iSiPM] = tfdird.make<TH1F>(Form("Crystal_%i_SiPM_%i_d",iCry, iSiPM),Form("Energy deposited in SiPM %i crystal %d only diagonal tracks",iSiPM, iCry),70, 15., 50.);
	hSiPMfp[iCry][iSiPM] = tfdirfp.make<TH1F>(Form("Crystal_%i_SiPM_%i_fp",iCry, iSiPM),Form("Energy deposited in SiPM %i crystal %d normalized track lenght",iSiPM, iCry),70, 15., 50.);
      }
    }
    MPV = tfdir_res.make<TH1F>("MPV_distr","MPV distribution",400, 10., 30.);
    MPVv = tfdir_res.make<TH1F>("MPV_v_distr","MPV distribution only vertical tracks",400, 10., 30.);
    MPVd = tfdir_res.make<TH1F>("MPV_d_distr","MPV distribution only diagonal tracks",400, 10., 30.);
    MPVfp = tfdir_res.make<TH1F>("MPV_fp_distr","MPV distribution all normalized tracks",400, 10., 30.);


    _outputfile.open("calibration_parameters.dat", std::ios::out);

  }//end begin job


  void CaloCosmicEnergy::analyze(art::Event const& event){
    art::ServiceHandle<GeometryService> geom;
    if ( !geom->hasElement<Calorimeter>() ){
      throw cet::exception("NO-CALO-GEOM") 
	<< "CaloEnergy: Calorimeter geometry not found";
    }

    const Calorimeter& cal = *(GeomHandle<Calorimeter>());

    art::Handle<CaloClusterCollection> caloClustersHandle;
    event.getByLabel(_caloClusterTag, caloClustersHandle);
    const CaloClusterCollection& caloClusters(*caloClustersHandle);

    int nCluster = caloClustersHandle -> size();

      for(int iClu = 0; iClu < nCluster; iClu++){
	ncry=0; //reset crystal counter each cluster
	Dx = 1400.;
	Dy = 1400.;
	if(_diagLevel > 0){
	  std::cout << "iClu: " << iClu << std::endl;
	}
	
	//resetting Array
	for(int f= 0; f< nCrystals; f++){ 
	  PosX[f]=0;
	  PosY[f]=0;
	  ErrPos[f]=9.81;// 34/sqrt(12)
	  IDs[f]=0;
	  path[f]=0;
	  whichHit[f]=0;
	}

	//loop over crystals
	for(long unsigned int iCry=0 ;iCry <caloClusters[iClu].caloHitsPtrVector().size(); iCry++){
	    
	  //cut on energy of crystal in cluster
	  if(caloClusters[iClu].caloHitsPtrVector()[iCry].get()->energyDep() > CutEnergyDep){
	    //counter n crystal above 10 MeV
	    PosX[ncry] = cal.geomUtil().mu2eToDiskFF(caloClusters[iClu].diskID(), cal.crystal(caloClusters[iClu].caloHitsPtrVector()[iCry].get()->crystalID()).position()).getX();
	    PosY[ncry] = cal.geomUtil().mu2eToDiskFF(caloClusters[iClu].diskID(), cal.crystal(caloClusters[iClu].caloHitsPtrVector()[iCry].get()->crystalID()).position()).getY();
	    IDs[ncry] = caloClusters[iClu].caloHitsPtrVector()[iCry].get()->crystalID();
	    whichHit[ncry] = iCry;
	    ncry++;
	  }

	}//end loop iCry in cluster
	//cut on ncrystalhit
	if(ncry > CutNCryHit){

	  float min_x =PosX[0];
	  float max_x =PosX[0]; //sbsitture xlictyposx->posX
	  
	  //bublbe sort to find Xmax and Xmin
	  for(int h =1; h< ncry; h++){
	    if(PosX[h] > max_x){
	      max_x = PosX[h];
	      //cout << "Max: " << max << endl;
	    }  
	    if(PosX[h] < min_x){
	      min_x = PosX[h];
	      //cout << "Min: " << min << endl;
	    }
	    Dx = abs(max_x-min_x);
	  }


	  float min_y =PosY[0];
	  float max_y =PosY[0];
	  
	  //bublbe sort to find Xmax and Xmin
	  for(int h =1; h< ncry; h++){
	    if(PosY[h] > max_y){
	      max_y = PosX[h];
	      //cout << "Max: " << max << endl;
	    }  
	    if(PosY[h] < min_y){
	      min_y = PosY[h];
	      //cout << "Min: " << min << endl;
	    }
	    Dy = abs(max_y-min_y);
	  }
	  
	  if(_diagLevel>0){
	    std::cout << "dx: " << Dx << " dy: "<< Dy << std::endl;
	  }
       
	  //geometry track
	  TGraph *gfirst = new TGraph(ncry, PosX, PosY);
	  TGraphErrors *gfit =  new TGraphErrors(ncry, PosX, PosY, ErrPos, ErrPos);	  
	  gfit->SetMarkerSize(0.9);
	  gfit->SetMarkerStyle(8);

	  float chi2norm =0.;


	  TFitResultPtr fitresult = gfirst->Fit("myPoly", "SQ");
	  TFitResultPtr fit = gfit->Fit("myPoly", "SQ");
	  //2 fits help convergence

	  chi2norm = fit->Chi2()/fit->Ndf();


	    //begin findpath

	  float m0 = -1;
	  m0 = myPoly->GetParameter(1);
	  if(Dx < 33){ //if vertical tracks fit could be wrong
	    m0 = 144.; //ca 89.6 deg
	  }


	  for(int h=0; h< ncry; h++){
	    if(_diagLevel>0){   
	      std::cout << "using m:" << m0 << " q: " << myPoly->GetParameter(0) << " x: " <<  PosX[h] << " y: "  << PosY[h]<< " ene: " << caloClusters[iClu].caloHitsPtrVector()[whichHit[h]].get()->energyDep() << std::endl;
	    }
	    path[h] = findpath(m0, myPoly->GetParameter(0), PosX[h], PosY[h]);
	
  }//end findpath

	  //only verical tracks
	  if(Dx < 34){
	    for(int iCry=0; iCry< ncry; iCry++){
	      for(int iSiPM =0; iSiPM< nSiPMs; iSiPM++){
		hSiPMv[IDs[iCry]][iSiPM]->Fill(caloClusters[iClu].caloHitsPtrVector()[whichHit[iCry]].get()->recoCaloDigis()[iSiPM].get()->energyDep());
	      }
	    }
	  }

	  //cut on diagonal tracks
	  if(chi2norm < 0.01){
	    for(int iCry=0; iCry< ncry; iCry++){
	      for(int iSiPM =0; iSiPM< nSiPMs; iSiPM++){
		hSiPMd[IDs[iCry]][iSiPM]->Fill(caloClusters[iClu].caloHitsPtrVector()[whichHit[iCry]].get()->recoCaloDigis()[iSiPM].get()->energyDep());
	      }
	    }
	  }

	  //cut on vertical tracks + good fit chi2
	  if(Dx < 35 || chi2norm < CutChi2Norm ){
	    for(int iCry=0; iCry< ncry; iCry++){
	      for(int iSiPM =0; iSiPM< nSiPMs; iSiPM++){
		hSiPM[IDs[iCry]][iSiPM]->Fill(caloClusters[iClu].caloHitsPtrVector()[whichHit[iCry]].get()->recoCaloDigis()[iSiPM].get()->energyDep());
		hSiPMfp[IDs[iCry]][iSiPM]->Fill(caloClusters[iClu].caloHitsPtrVector()[whichHit[iCry]].get()->recoCaloDigis()[iSiPM].get()->energyDep()*34.3/path[iCry]);
	      
	      }
	    }
	  }// chi2 cut

	  delete gfit;
	  delete gfirst; 
	  
	}//end ncrystal cut

      }//end loop iClu in ncluster
  }//end enalyze

  void CaloCosmicEnergy::endJob(){
    std::cout << "CaloCosmicEnergy: entering endjob" << std::endl;

    //fit params
    int redo[nDisks*nCrystals][nSiPMs];
    int ibin = 0;
    float Xmean = 0;
    float Xsigma= 0;
    float params[nDisks*nCrystals][nSiPMs][5];

    for(int k=0; k< nDisks*nCrystals; k++){
      for(int j =0; j < nSiPMs; j++){
	redo[k][j]=0;
	for(int i = 0;i< 5; i++){
	  params[k][j][i]=-1;
	}
      }
    }
    for(int iCry=0; iCry< nDisks*nCrystals; iCry++){
      for(int iSiPM=0; iSiPM < nSiPMs; iSiPM++){

	//set initial parameters

	ibin = hSiPM[iCry][iSiPM]->GetMaximumBin();
	Xmean = hSiPM[iCry][iSiPM]->GetBinCenter(ibin);
	Xsigma= hSiPM[iCry][iSiPM]->GetRMS();

	mylang->SetParameters(Xsigma,Xmean,hSiPM[iCry][iSiPM]->GetMaximum(), 1); 
	hSiPM[iCry][iSiPM]->Fit("mylang","Q","",10,50);


	if(_diagLevel>0){
	  std::cout<< "fit complete "<< std::endl;
	}
	if((mylang->GetChisquare()/mylang->GetNDF()) < 1.8 && (hSiPM[iCry][iSiPM]->GetMean() - mylang->GetParameter(1)) < 10.){
	  if(_diagLevel>0){
	    std::cout << "Crystal: " << iCry << " SiPM: "<< iSiPM << " MPV: " << mylang->GetParameter(1) << std::endl;
	  }

	 MPV->Fill(mylang->GetParameter(1));
	 
	 params[iCry][iSiPM][0]= mylang->GetParameter(1);
	 params[iCry][iSiPM][1]= mylang->GetParError(1);
	 params[iCry][iSiPM][2]= mylang->GetParameter(0);
	 params[iCry][iSiPM][3]= mylang->GetParError(0);
	 params[iCry][iSiPM][4]= mylang->GetChisquare()/mylang->GetNDF();

	}
	else { redo[iCry][iSiPM]=1; }

      }//end sipms loop
    }//end crystal loop
    if(_diagLevel >0){
      std::cout << " I am redoing fits" << std::endl;
}
  //redo wrong fit
  for(int iCry=0; iCry< nDisks*nCrystals; iCry++){
    for(int iSiPM=0; iSiPM< nSiPMs; iSiPM++){
      if(redo[iCry][iSiPM]==1){
	  ibin = hSiPM[iCry][iSiPM]->GetMaximumBin();
	  Xmean = hSiPM[iCry][iSiPM]->GetBinCenter(ibin);
	  Xsigma= hSiPM[iCry][iSiPM]->GetRMS();

	  //set wider parameter
	  mylang->SetParameters(Xsigma*1.2,Xmean,hSiPM[iCry][iSiPM]->GetMaximum(), 1.1); 

	 hSiPM[iCry][iSiPM]->Fit("mylang","Q","",10,50); 

	 MPV->Fill(mylang->GetParameter(1));

	 params[iCry][iSiPM][0]= mylang->GetParameter(1);
	 params[iCry][iSiPM][1]= mylang->GetParError(1);
	 params[iCry][iSiPM][2]= mylang->GetParameter(0);
	 params[iCry][iSiPM][3]= mylang->GetParError(0);
	 params[iCry][iSiPM][4]= mylang->GetChisquare()/mylang->GetNDF();

	 if(_diagLevel>0){
	   std::cout << "Crystal: " << iCry << " SiPM: "<< iSiPM << " MPV: " << mylang->GetParameter(1) << std::endl;
	 }	 

	}
    }
  }
  _outputfile  << "#Cry #SiPM #MPV #MPV err #sigm #sigma err #chi2/ndf \n\n"; 
  for(int iCry=0; iCry < nDisks*nCrystals; iCry++){
    for(int iSiPM=0; iSiPM < nSiPMs; iSiPM++){
      _outputfile << Form("%i %i %.3f %.3f %.3f %.3f %.3f\n", iCry, iSiPM, params[iCry][iSiPM][0], params[iCry][iSiPM][1], params[iCry][iSiPM][2], params[iCry][iSiPM][3], params[iCry][iSiPM][4]); 
      }
      	  _outputfile <<  "\n";
    }

    _outputfile.close();


    //fits for vertical diagonal and nromalized tracks
    for(int iCry=0; iCry< nDisks*nCrystals; iCry++){
      for(int iSiPM=0; iSiPM < nSiPMs; iSiPM++){

	//set initial parameters

	ibin = hSiPMv[iCry][iSiPM]->GetMaximumBin();
	Xmean = hSiPMv[iCry][iSiPM]->GetBinCenter(ibin);
	Xsigma= hSiPMv[iCry][iSiPM]->GetRMS();

	mylang->SetParameters(Xsigma,Xmean,hSiPMv[iCry][iSiPM]->GetMaximum(), 1); 
	hSiPMv[iCry][iSiPM]->Fit("mylang","Q","",10,50);
	MPVv->Fill(mylang->GetParameter(1));

	ibin = hSiPMd[iCry][iSiPM]->GetMaximumBin();
	Xmean = hSiPMd[iCry][iSiPM]->GetBinCenter(ibin);
	Xsigma= hSiPMd[iCry][iSiPM]->GetRMS();

	mylang->SetParameters(Xsigma,Xmean,hSiPMd[iCry][iSiPM]->GetMaximum(), 1); 
	hSiPMd[iCry][iSiPM]->Fit("mylang","Q","",10,50);
	MPVd->Fill(mylang->GetParameter(1));

	ibin = hSiPMfp[iCry][iSiPM]->GetMaximumBin();
	Xmean = hSiPMfp[iCry][iSiPM]->GetBinCenter(ibin);
	Xsigma= hSiPMfp[iCry][iSiPM]->GetRMS();

	mylang->SetParameters(Xsigma,Xmean,hSiPMfp[iCry][iSiPM]->GetMaximum(), 1); 
	hSiPMfp[iCry][iSiPM]->Fit("mylang","Q","",10,50);
	MPVfp->Fill(mylang->GetParameter(1));

      }
    }
    MPV->Fit("gaus","Q","",10,30.);
    MPVv->Fit("gaus","Q","",10,30.);
    MPVd->Fit("gaus","Q","",10,30.);
    MPVfp->Fit("gaus","Q","",10,30.);

  }//end endjob
 double CaloCosmicEnergy::langaufun(double *x, double *par) {

   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation), 
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.

      // Numeric constants
      Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      Double_t mpshift  = -0.22278298;       // Landau maximum location

      // Control constants
      Double_t np = 100.0;      // number of convolution steps
      Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

      // Variables
      Double_t xx;
      Double_t mpc;
      Double_t fland;
      Double_t sum = 0.0;
      Double_t xlow,xupp;
      Double_t step;
      Double_t i;


      // MP shift correction
      mpc = par[1] - mpshift * par[0]; 

      // Range of convolution integral
      xlow = x[0] - sc * par[3];
      xupp = x[0] + sc * par[3];

      step = (xupp-xlow) / np;

      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);

         xx = xupp - (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
      }

      return (par[2] * step * sum * invsq2pi / par[3]);
}//end langaus

  float CaloCosmicEnergy::findpath(float m,float q, float x, float y){
  float xup=0;
  float yright=0;
  float path=0;
  float xlow=0;
  float yleft=0;

  int diag =0;
  float halfcrysize = 17.15;
  
  if(m != 0){
    xup = ((y+halfcrysize)- q)/m;
    xlow = ((y-halfcrysize)- q)/m;
    yleft = m * (x - halfcrysize) + q;
    yright = m * (x + halfcrysize) + q;
  }
  if(xup <(x + halfcrysize) && xup > (x-halfcrysize)){ //hit the top of crystal
    if( xlow < (x+halfcrysize) && xlow > (x-halfcrysize)){
      path = halfcrysize*2 * TMath::Sqrt(1/(m*m) +1);
      if(diag==1){
	std::cout << "x0,y0: (" <<xup << ";" <<y+halfcrysize << ") x1,y1: (" <<xlow << ";" << y-halfcrysize<< ")" << 	std::endl;
      }
    }
    if(xlow < (x - halfcrysize)){
      path = TMath::Sqrt((xup-(x-halfcrysize))*(xup-(x-halfcrysize)) + ((y+halfcrysize) - yleft)*((y+halfcrysize) - yleft));
      if(diag==1){
		std::cout << "Track exit from left" << std::endl;
	std::cout << "x0,y0: (" << xup<< ";" <<y+halfcrysize << ") x1,y1: (" << x-halfcrysize<< ";" << yleft<< ")" << std::endl;
      }
    }
    if(xlow > (x + halfcrysize)){
      path = TMath::Sqrt((xup-(x+halfcrysize))*(xup-(x+halfcrysize)) + ((y+halfcrysize) - yright)*((y+halfcrysize) - yright));
      if(diag==1){
	std::cout << "x0,y0: (" <<xup << ";" << y+halfcrysize<< ") x1,y1: (" <<x+halfcrysize << ";" << yright<< ")" << std::endl;
	std::cout << "Track exit from right" << std::endl;
      }
    }
  }
  //enter from left
  else if(yleft < (y+halfcrysize) && yleft > (y-halfcrysize)){
    if(diag==1){
      std::cout << "track enter from left and ";
    }
    if(xlow < (x +halfcrysize) && xlow > (x -halfcrysize)){
      if(diag==1){
	std::cout << " it exit bot " << std::endl;
	std::cout << "x0,y0: (" <<x-halfcrysize << ";" << yleft<< ") x1,y1: (" <<xlow << ";" << y-halfcrysize<< ")" << std::endl;
      }
      path = TMath::Sqrt((x-halfcrysize - xlow)*(x-halfcrysize - xlow) +(yleft - (y-halfcrysize) )*(yleft - (y-halfcrysize)));
    }
    else{ //track exit from right
      if(diag==1){
	std::cout << " it exit from right"<< std::endl;
	std::cout << "x0,y0: (" <<x-halfcrysize << ";" << yleft<< ") x1,y1: (" <<x+halfcrysize << ";" << yright<< ")" << std::endl;
      }
      path = halfcrysize*2*TMath::Sqrt(1 + m*m);
    }
  }

  else if(yright < (y+halfcrysize) && yright > (y-halfcrysize)){
    if(diag==1){
      std::cout << "track enter from right and ";
    }
    if(xlow < (x +halfcrysize) && xlow > (x -halfcrysize)){
      if(diag==1){
	std::cout << " it exit bot " << std::endl;
	std::cout << "x0,y0: (" <<x+halfcrysize << ";" << yright<< ") x1,y1: (" <<xlow << ";" << y-halfcrysize<< ")" << std::endl;
      }
      path = TMath::Sqrt((x+halfcrysize - xlow)*(x+halfcrysize - xlow) +(yright - (y-halfcrysize) )*(yright - (y-halfcrysize)));
    }
  }

  else{
    if(diag==1){
      std::cout << "can't find a path" << std::endl;
    }
  }
  if(diag==1){
    std::cout << "path is: " << path << std::endl;
  }
  
  return path;
}//end findpath


}//end namespace

DEFINE_ART_MODULE(mu2e::CaloCosmicEnergy)
