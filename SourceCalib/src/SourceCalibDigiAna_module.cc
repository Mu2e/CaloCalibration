#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

#include "Offline/ConditionsService/inc/ConditionsHandle.hh"
#include "Offline/ConditionsService/inc/CalorimeterCalibrations.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/RecoDataProducts/inc/CaloDigi.hh"
#include "Offline/RecoDataProducts/inc/CaloRecoDigi.hh"
#include "Offline/RecoDataProducts/inc/ProtonBunchTime.hh"
#include "Offline/CaloReco/inc/CaloWaveformProcessor.hh"
#include "Offline/CaloReco/inc/CaloTemplateWFProcessor.hh"
#include "Offline/CaloReco/inc/CaloRawWFProcessor.hh"
#include "Offline/DAQConditions/inc/EventTiming.hh"
#include "TH1F.h"
#include <iostream>
#include <string>
#include <sstream>
#include "art_root_io/TFileService.h"
#include "TTree.h"

constexpr int ncrystals = 1348;
constexpr int nsipms = 2696;

int Contains(std::vector<int> v, int x)
{
  return std::count(v.begin(), v.end(), x);
}

//std::vector<int> badcrys = {774, 775,776,792, 793,794,795,796,798,814,815,816,818,819,820,822,851,865,881,941, 965,969,1006,1029,1039,1071,1097,1144,1183,1209,1243,1275,1281,1316,1336};
namespace mu2e {

  class SourceCalibDigiAna : public art::EDAnalyzer
  {
     public:
        enum processorStrategy {NoChoice, RawExtract, Template};

        struct Config
        {
           using Name    = fhicl::Name;
           using Comment = fhicl::Comment;
           fhicl::Table<mu2e::CaloRawWFProcessor::Config>      proc_raw_conf       { Name("RawProcessor"),        Comment("Raw processor config") };
           fhicl::Table<mu2e::CaloTemplateWFProcessor::Config> proc_templ_conf     { Name("TemplateProcessor"),   Comment("Log normal fit processor config") };
           fhicl::Atom<art::InputTag>                          caloDigiCollection  { Name("caloDigiCollection"),  Comment("Calo Digi module label") };
           fhicl::Atom<art::InputTag>                          pbttoken            { Name("ProtonBunchTimeTag"),  Comment("ProtonBunchTime producer")};
           fhicl::Atom<std::string>                            processorStrategy   { Name("processorStrategy"),   Comment("Digi reco processor name") };
           fhicl::Atom<double>                                 digiSampling        { Name("digiSampling"),        Comment("Calo ADC sampling time (ns)") };
           fhicl::Atom<double>                                 maxChi2Cut          { Name("maxChi2Cut"),          Comment("Chi2 cut for keeping reco digi") };
           fhicl::Atom<double>                                 ratioCut          { Name("ratioCut"),          Comment("ratio of energy in crystal to others"), 0.8 };
           fhicl::Atom<double>                                timeCut          { Name("timeCut"),          Comment("time for digis in same event"), 20 };
           fhicl::Atom<int>                                    diagLevel           { Name("diagLevel"),           Comment("Diagnosis level") };
           fhicl::Atom<double>                                    endBin           { Name("endBin"),           Comment("finalBin in ADC"),120 };
           fhicl::Atom<bool>                                    docuts           { Name("docuts"),           Comment("apply cuts"),true };
        };

        explicit SourceCalibDigiAna(const art::EDAnalyzer::Table<Config>& config) :
           EDAnalyzer{config},
           caloDigisToken_    {consumes<CaloDigiCollection>(config().caloDigiCollection())},
           pbttoken_          {consumes<ProtonBunchTime>(config().pbttoken())},
           processorStrategy_ (config().processorStrategy()),
           digiSampling_      (config().digiSampling()),
           maxChi2Cut_        (config().maxChi2Cut()),
           ratioCut_          (config().ratioCut()),
           timeCut_          (config().timeCut()),
           diagLevel_         (config().diagLevel()),
           endBin_         (config().endBin()),
           docuts_          (config().docuts())
        {
            std::map<std::string, processorStrategy> spmap;
            spmap["RawExtract"]  = RawExtract;
            spmap["TemplateFit"] = Template;
            switch (spmap[processorStrategy_])
            {
                case RawExtract:
                {
                    waveformProcessor_ = std::make_unique<CaloRawWFProcessor>(config().proc_raw_conf());
                    break;
                }
                case Template:
                {
                    waveformProcessor_ = std::make_unique<CaloTemplateWFProcessor>(config().proc_templ_conf());
                    break;
                }
                default:
                {
                    throw cet::exception("CATEGORY")<< "Unrecognized processor in SourceCalibDigiAna module";
                }
            }
        }
        //explicit SourceCalibDigiAna(const art::EDAnalyzer::Table<Config>& config);
        virtual ~SourceCalibDigiAna() {}
        virtual void beginJob();
        virtual void beginRun(const art::Run& aRunn);
        virtual void analyze(const art::Event& e);

     private:
        //std::ofstream badfile;
        //std::ofstream goodfile;
        void extractRecoDigi(const art::ValidHandle<CaloDigiCollection>&, double );

        const  art::ProductToken<CaloDigiCollection> caloDigisToken_;
        const  art::ProductToken<ProtonBunchTime>    pbttoken_;
        const  std::string                           processorStrategy_;
        double                                       digiSampling_;
        double                                       maxChi2Cut_;
        double                                       ratioCut_;
        double                                       timeCut_;
        int                                          diagLevel_;
        double                                       endBin_;
        bool                                          docuts_;
        std::unique_ptr<CaloWaveformProcessor>       waveformProcessor_;
        TH1F* list_of_crys_hists[ncrystals];
        TH1F* list_of_sipm_hists[nsipms];
  };

  
  void SourceCalibDigiAna::analyze(const art::Event& event)
  {
      if (diagLevel_ > 0) std::cout<<"[SourceCalibDigiAna::analyze] begin"<<std::endl;
      const auto& caloDigisH = event.getValidHandle(caloDigisToken_);
      auto pbtH = event.getValidHandle(pbttoken_);
      const ProtonBunchTime& pbt(*pbtH);
      double pbtOffset = pbt.pbtime_;
      extractRecoDigi(caloDigisH,  pbtOffset);
      if (diagLevel_ > 0) std::cout<<"[SourceCalibDigiAna::analyze] end"<<std::endl;
  }

  void SourceCalibDigiAna::beginJob(){
      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory crydir = tfs->mkdir( "crystals_ADC" );
      art::TFileDirectory sipmdir = tfs->mkdir( "sipm_ADC" );
      for(int i = 0; i < ncrystals ; i++){
        TString histname = "cry_"+std::to_string(i);
        list_of_crys_hists[i] = crydir.make<TH1F>( histname , histname, 300, 0.0, endBin_);
      }
      for(int i = 0; i < nsipms ; i++){
        TString histname = "sipm_"+std::to_string(i);
        list_of_sipm_hists[i] = sipmdir.make<TH1F>( histname , histname, 300, 0.0, endBin_);
      }
      //badfile.open("badcrys.csv");
      //goodfile.open("goodcrys.csv");
  }

  //--------------------------------------------------
  void SourceCalibDigiAna::beginRun(const art::Run& aRun)
  {
      waveformProcessor_->initialize();
  }

  //------------------------------------------------------------------------------------------------------------
  void SourceCalibDigiAna::extractRecoDigi(const art::ValidHandle<CaloDigiCollection>& caloDigisHandle, double pbtOffset)
  {
      
      const auto& caloDigis = *caloDigisHandle;
      ConditionsHandle<CalorimeterCalibrations> calorimeterCalibrations("ignored");

      double totEnergyReco(0);
      std::vector<double> x{},y{};
      

      std::vector<int> crystals_in_event;
      std::vector<int> sipms_in_event;
      std::vector<double> time;
      std::vector<double> total_energy_in_crystal(ncrystals, 0);
      std::vector<double> total_energy_in_sipm(nsipms, 0);
      for (const auto& caloDigi : caloDigis)
      {
          int    SiPMID   = caloDigi.SiPMID();
          double t0       = caloDigi.t0();
          time.push_back(t0);
          const std::vector<int>& waveform = caloDigi.waveform();
          //size_t index = &caloDigi - &caloDigis.front();
          //art::Ptr<CaloDigi> caloDigiPtr(caloDigisHandle, index);

          x.clear();y.clear();
          for (unsigned int i=0;i<waveform.size();++i)
          {
              x.push_back(t0 + (i+0.5)*digiSampling_); // add 0.5 to be in middle of bin
              y.push_back(waveform.at(i));
          }

          waveformProcessor_->reset();
          waveformProcessor_->extract(x,y);

          for (int i=0;i<waveformProcessor_->nPeaks();++i)
          {
              double chi2      = waveformProcessor_->chi2();
              int    ndf       = waveformProcessor_->ndf();
              double eDep      = waveformProcessor_->amplitude(i);
              if (SiPMID%2==0) totEnergyReco += eDep;
              if (chi2/float(ndf) > maxChi2Cut_) continue;
              
              if(crystals_in_event.size() !=0 ) {
                if (Contains(crystals_in_event, SiPMID/2) == 0) {
                  crystals_in_event.push_back(SiPMID/2);
                  total_energy_in_crystal[SiPMID/2] = eDep;
                } else{ 
                  total_energy_in_crystal[SiPMID/2]+= eDep;
                }
              } else crystals_in_event.push_back(SiPMID/2);
                
                if (Contains(sipms_in_event, caloDigi.SiPMID()) == 0) {
                  sipms_in_event.push_back(caloDigi.SiPMID());
                  total_energy_in_sipm[caloDigi.SiPMID()] = eDep;
                } else{
                  total_energy_in_sipm[caloDigi.SiPMID()]+= eDep;
                }
                sipms_in_event.push_back(caloDigi.SiPMID());
          }
      }

      // For crystals:
      bool passes_time_cry = true; 
      bool passes_ratio_cry = true;
      
      double difTime_cry = 0;
      if(time.size() !=0) {
        sort(time.begin(), time.end());
        difTime_cry = time.back() - time.front();
      }
      for(const auto& id : crystals_in_event){
        float edepTarget = 0.0;
        float edepOthers   = 0.0;
        int id_sipm1 = id*2;
        int id_sipm2 = id*2+1;
        for(const auto& id2 : crystals_in_event){
          if(id == id2) edepTarget += total_energy_in_crystal[id];
          else edepOthers += total_energy_in_crystal[id];
        }
        if(docuts_ and difTime_cry > timeCut_){
          passes_time_cry = false;
        }

        if(docuts_ and edepTarget / (edepTarget + edepOthers) < ratioCut_){
          passes_ratio_cry = false;
        }
        if( (passes_time_cry and passes_ratio_cry)){
          list_of_crys_hists[id]->Fill(total_energy_in_crystal[id]);
          if(total_energy_in_sipm[id_sipm1]!=0) list_of_sipm_hists[id_sipm1]->Fill(total_energy_in_sipm[id_sipm1]);
          if(total_energy_in_sipm[id_sipm2]!=0) list_of_sipm_hists[id_sipm2]->Fill(total_energy_in_sipm[id_sipm2]);
        }
        /*if (passes_time and passes_ratio and Contains(badcrys, id) == 1) { 
          badfile<<id<<","<<passes_time<<","<<passes_ratio<<","<<total_energy_in_crystal[id]<<std::endl;
        }*/
      }
      
      // For SiPMs
      /*bool passes_time_sipm = true; 
      bool passes_ratio_sipm = true;
      double difTime_sipm = 0;
      if(time.size() !=0) {
        sort(time.begin(), time.end());
        difTime_sipm = time.back() - time.front();
      }
      for(const auto& id : sipms_in_event){
        float edepTarget = 0.0;
        float edepOthers   = 0.0;
        for(const auto& id2 : sipms_in_event){
          if(id == id2) edepTarget += total_energy_in_sipm[id];
          else edepOthers += total_energy_in_sipm[id];
        }
        if(difTime_sipm > timeCut_){
          passes_time_sipm = false;
        }

        if(passes_time_sipm and passes_ratio_sipm and total_energy_in_sipm[id]!=0){
          //list_of_sipm_hists[id]->Fill(total_energy_in_sipm[id]);
          std::cout<<"ID "<<id<<" energy "<<total_energy_in_sipm[id]<<std::endl;
        }
      }*/
      if (diagLevel_ > 1) std::cout<<"[SourceCalibDigiAna] Total energy reco "<<totEnergyReco <<std::endl;
  }
}
DEFINE_ART_MODULE(mu2e::SourceCalibDigiAna)
