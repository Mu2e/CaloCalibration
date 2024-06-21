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
int Contains(std::vector<int> v, int x)
{
  return std::count(v.begin(), v.end(), x);
}
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
           diagLevel_         (config().diagLevel())
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
        void extractRecoDigi(const art::ValidHandle<CaloDigiCollection>&, double );

        const  art::ProductToken<CaloDigiCollection> caloDigisToken_;
        const  art::ProductToken<ProtonBunchTime>    pbttoken_;
        const  std::string                           processorStrategy_;
        double                                       digiSampling_;
        double                                       maxChi2Cut_;
        double                                       ratioCut_;
        double                                       timeCut_;
        int                                          diagLevel_;
        std::unique_ptr<CaloWaveformProcessor>       waveformProcessor_;
        TH1F* list_of_hists[ncrystals];
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
      art::TFileDirectory tfdir = tfs->mkdir( "crystals_ADC" );
      for(int i = 0; i < 1348 ; i++){
        TString histname = "hspec_"+std::to_string(i);
        list_of_hists[i] = tfdir.make<TH1F>( histname , histname, 300, 0.0, 150);
      }
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
      
      //int idExist = 0;
      std::vector<int> crystals_in_event;
      std::vector<double> time;
      double total_energy_in_crystal[ncrystals];

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
                } else{ total_energy_in_crystal[SiPMID/2]+= eDep; }
              } else crystals_in_event.push_back(SiPMID/2);
              //list_of_hists[SiPMID/2]->Fill(waveformProcessor_->amplitude(i));
          }
      }
      
      bool passes_time = true; 
      bool passes_ratio = true;
      double difTime = 0;
      if(time.size() !=0) {
        sort(time.begin(), time.end());
        difTime = time.back() - time.front();
      }
      for(const auto& id : crystals_in_event){
        float edepTarget = 0.0;
        float edepOthers   = 0.0;
        for(const auto& id2 : crystals_in_event){
          if(id == id2) edepTarget += total_energy_in_crystal[id];
          else edepOthers += total_energy_in_crystal[id];
        }
        std::cout<<difTime<<std::endl;
        if(difTime > timeCut_){
          passes_time = false;
        }

        if(edepTarget / (edepTarget + edepOthers) < ratioCut_){
          passes_ratio = false;
        }
        if(passes_time and passes_ratio){
          list_of_hists[id]->Fill(total_energy_in_crystal[id]);
        }
      }
      if (diagLevel_ > 1) std::cout<<"[SourceCalibDigiAna] Total energy reco "<<totEnergyReco <<std::endl;
  }
}
DEFINE_ART_MODULE(mu2e::SourceCalibDigiAna)
