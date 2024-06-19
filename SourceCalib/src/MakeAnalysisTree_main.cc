#include "CaloCalibration/SourceCalib/inc/MakeAnalysisTree.hh"
using namespace CaloSourceCalib;

TString filepath = "/pnfs/mu2e/tape/usr-nts/nts/hjafree/SourceCalibAna/MDC2020ae/root/0a/1f/nts.hjafree.SourceCalibAna.MDC2020ae.0.root";

/*function to extract the TTree from the SourceCalibAna output*/
TTree* get_data_tree(){
    TFile *f =  new TFile(filepath);
    TTree *t = (TTree*)f->Get("CaloExample/Calo");
    return t;
}

/* function to loop over ntuple from SourceCalibAna and bin by crystal*/
void MakeCrystalBinsOutputs(){

  for(unsigned int numVal = 674; numVal < 676; numVal++){
    std::cout<<" Analyzing Crystal # "<<numVal<<std::endl;

    TString cryNum;
    cryNum = to_string(numVal);

    // create the output file.
    TString ouptName = "mu2e_caloSimu_crySpec_" + cryNum + ".root";
    TFile *ouptFile = new TFile(ouptName, "RECREATE");
    //create output tree
    Float_t spec, ratio, ntrig, stim, time, tErg;
    TTree *outTree = new TTree("outTree","outTree");
    outTree->Branch("spec", &spec, "spec/F"); //Energy Deposited
    outTree->Branch("ratio", &ratio, "ratio/F"); //Ratio of Energy Deposited
    outTree->Branch("ntrig", &ntrig,"ntrig/F");//Number of Hits on Target Crystal
    outTree->Branch("stim", &stim, "stim/F");//Time Difference from Same Crystal
    outTree->Branch("time", &time, "time/F");//Maxium Time Difference for All Hits
    outTree->Branch("tErg", &tErg, "tErg/F");//Energy Spectrum with Special Time Difference

    // input tree
    TTree *inTree = get_data_tree();

    unsigned int nEvt = (int)inTree -> GetEntries();
    for(unsigned int iEvt=0; iEvt<nEvt; iEvt++)
    {
      inTree -> GetEntry(iEvt);
      int idExist = 0;
      std::vector<float> sameCryTime;
      std::vector<float> sameCryEdep;
      
      for(unsigned int icry=0; icry<nCry; icry++)
      {
        if(cryId[icry] == numVal)
        {
          idExist += 1;
          sameCryTime.push_back(cryTime[icry]);
          sameCryEdep.push_back(cryEdep[icry]);
        }
      }
      
      float cryEdep[100];
      int nCry;
      int cryId[100];
      float cryTime[100];
      float cryPosX[100], cryPosY[100], cryPosZ[100];

      inTree -> SetBranchAddress("nCry", &nCry);
      inTree -> SetBranchAddress("cryId", &cryId);
      inTree -> SetBranchAddress("cryEdep", &cryEdep);
      inTree -> SetBranchAddress("cryTime", &cryTime);
      inTree -> SetBranchAddress("cryPosX", &cryPosX);
      inTree -> SetBranchAddress("cryPosY", &cryPosY);
      inTree -> SetBranchAddress("cryPosZ", &cryPosZ);


      if(idExist == 0) continue;
      sort(sameCryTime.begin(), sameCryTime.end());
      float sameCryDeltaT = sameCryTime.back() - sameCryTime.front();
      stim = sameCryDeltaT;

      float ergTarget = 0.0;
      float ergRest   = 0.0;
      int nTarget = 0;
      std::vector<float> deltaTime;

      for(int icry=0; icry<nCry; icry++)
      {
        if(cryId[icry] == numVal)
        {
          ergTarget = cryEdep[icry];
          deltaTime.push_back(cryTime[icry]);
          nTarget += 1;
        }
        else
        {
          ergRest += cryEdep[icry];
          deltaTime.push_back(cryTime[icry]);
        }
      }
      sort(deltaTime.begin(), deltaTime.end());
      float difTime = deltaTime.back() - deltaTime.front();
      int numHit = sameCryEdep.size();
      if(difTime > 20)
      {
        for(unsigned int iHit = 0; iHit < numHit; iHit++)
        {
          tErg = sameCryEdep[iHit];
        }
      }
      if((ergTarget / (ergTarget + ergRest)) >= 0.8 && difTime < 4)
      {
        spec =ergTarget;
      }
      ratio= (ergTarget / (ergTarget + ergRest));
      ntrig = (nTarget);
      time  = (difTime);
      outTree->Fill();
    }

  ouptFile -> Write();
  ouptFile -> Close();

  }
 }
 
int main(){//int argc, char* argv[]){

  MakeCrystalBinsOutputs();

}
