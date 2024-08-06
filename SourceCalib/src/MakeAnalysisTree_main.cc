#include "CaloCalibration/SourceCalib/inc/MakeAnalysisTree.hh"
#include "CaloCalibration/SourceCalib/inc/SourceFitter.hh"
#include <chrono>
using namespace std::chrono;

using namespace CaloSourceCalib;

TString filepath = "/pnfs/mu2e/tape/phy-nts/nts/mu2e/SourceCalibAna/MDC2020t/root/26/3a/nts.mu2e.SourceCalibAna.MDC2020t.0.root";//"/pnfs/mu2e/tape/usr-nts/nts/hjafree/SourceCalibAna/MDC2020ae/root/0a/1f/nts.hjafree.SourceCalibAna.MDC2020ae.0.root"

/*function to extract the TTree from the SourceCalibAna output*/
TTree* get_data_tree(){
    TFile *f =  new TFile(filepath);
    TTree *t = (TTree*)f->Get("CaloExample/Calo");
    //TTree *t = (TTree*)f->Get("CaloRecoDigiMaker/SourceCalibDigiAna"); -->FOR ADC
    return t;
}

/* bin for a given crystal */
void AnalyzeCrystal(int crystalNo){
    
    TString cryNum;
    cryNum = to_string(crystalNo);

    // create the output file.
    TString ouptName = "mu2e_caloSimu_crysEdep_" + cryNum + ".root";
    TFile *ouptFile = new TFile(ouptName, "RECREATE");
    //create output tree
    Float_t crysEdep;//, ratio, ntrig, stim, time, tErg;
    TTree *outTree = new TTree("crysTree","crysTree");
    outTree->Branch("crysEdep", &crysEdep, "crysEdep/F"); //Energy Deposited
    
    // make histogram
    TH1F *h_spec = new TH1F("h_spec", "h_spec", 300, 0.0, 7.5);
    // input tree
    TTree *inTree = get_data_tree();

    // extract branches - these are arrays as there will be multiple crystal hits per event
    int  nCry;//number of crystal hits in event
    int cryId[10];
    float cryTime[10];
    float cryPosX[10], cryPosY[10], cryPosZ[10];
    float cryEdep[10];
    
    inTree -> SetBranchAddress("nCrystals", &nCry);//ncalhitHit
    inTree -> SetBranchAddress("cryId", &cryId);//same
    inTree -> SetBranchAddress("cryEdep", &cryEdep);//calhitRecoEdep
    inTree -> SetBranchAddress("cryTime", &cryTime);//calhitRecoTime
    inTree -> SetBranchAddress("cryPosX", &cryPosX);//calhitRecoPosX
    inTree -> SetBranchAddress("cryPosY", &cryPosY);//calhitRecoPosY
    inTree -> SetBranchAddress("cryPosZ", &cryPosZ);//calhitRecoPosZ

    unsigned int nEvt = (int)inTree -> GetEntries();
    unsigned int nEvtCrys = 0; //events in this crystal
    for(unsigned int iEvt=0; iEvt<nEvt; iEvt++)
    {
      // extract single entry in the Tree
      inTree -> GetEntry(iEvt);
      
      int idExist = 0;
      std::vector<float> sameCryTime;
      std::vector<float> sameCryEdep;
      
      for(int icry=0; icry<nCry; icry++)
      {
        
        if(cryId[icry] == crystalNo)
        {
          idExist += 1;
          sameCryTime.push_back(cryTime[icry]);
          sameCryEdep.push_back(cryEdep[icry]);
        }
      }

      if(idExist == 0) continue;
      sort(sameCryTime.begin(), sameCryTime.end());
      float edepTarget = 0.0;
      float edepOthers   = 0.0;

      std::vector<float> deltaTime;

      for(int icry=0; icry<nCry; icry++)
      {
        if(cryId[icry] == crystalNo)
        {
          edepTarget = cryEdep[icry];
          deltaTime.push_back(cryTime[icry]);

        }
        else
        {
          edepOthers += cryEdep[icry];
          deltaTime.push_back(cryTime[icry]);
        }
      }
      sort(deltaTime.begin(), deltaTime.end());
      float difTime = deltaTime.back() - deltaTime.front();
      
      if((edepTarget / (edepTarget + edepOthers)) >= 0.8 && difTime < 4)
      {
        crysEdep =edepTarget;
        h_spec->Fill(edepTarget);
        nEvtCrys+=1;
      }
      outTree->Fill();
      
    }
    std::cout<<" Events analyzed in this crystal : "<<nEvtCrys<<std::endl;
  ouptFile -> Write();
  ouptFile -> Close();

}

/*function calls RooFit and fits a given crystal file*/
void RunRooFit(int crystalNo) {
  SourceFitter *fit = new SourceFitter();
  fit->FitCrystal(crystalNo,"nll");
}


/* function to loop over ntuple from SourceCalibAna and bin by crystal*/
void MakeCrystalBinsOutputs( int start,  int end){
  for(int crystalNo = start; crystalNo < end; crystalNo++){
    std::cout<<" Finding Hits in Crystal # "<<crystalNo<<std::endl;
    auto start_bin = high_resolution_clock::now();
    AnalyzeCrystal(crystalNo);
    auto end_bin = high_resolution_clock::now();
    std::cout<<" Time take to filter crystal "<<duration_cast<seconds>(end_bin - start_bin)<<std::endl;
    std::cout<<" Fitting Crystal # "<<crystalNo<<std::endl;
    auto start_fit = high_resolution_clock::now();
    RunRooFit(crystalNo);
    auto end_fit = high_resolution_clock::now();
    std::cout<<" Time take to fit crystal "<<duration_cast<seconds>(end_fit - start_fit)<<std::endl;
 }
}


/*main function allows a loop over all crystals or a choice of a single crystal*/
int main(int argc, char* argv[]){
  std::cout<<"========== Welcome to the Mu2e Source Calibration Analysis =========="<<std::endl;
  int anacrys_start = 675; //starting crystal
  int anacrys_end = 680; //final crystal
  std::cout<<"choosing crystal "<<std::endl;
  if(strcmp( argv[1], "chooseCrystal") == 0 ){
    cout<<"crystal to be analyzed (int) : "<<endl;
    cin>>anacrys_start;
    anacrys_end = anacrys_start+1;
  }
  std::cout<<"Running crystal binning ....."<<std::endl;
  //MakeCrystalListOutputs(anacrys_start,anacrys_end);
  MakeCrystalBinsOutputs(anacrys_start,anacrys_end);
  std::cout<<"Finished processing ..."<<std::endl;
  return 0;
}
