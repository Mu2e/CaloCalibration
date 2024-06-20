/*// Simplied analyzer for source calibration analysis - will work with digi ADC waveforms
// Author: Sophie Middleton 2022, 2024
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Provenance.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Utilities/InputTag.h"

#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/ConditionsService/inc/ConditionsHandle.hh"

#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Offline/DataProducts/inc/CaloSiPMId.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/VirtualDetector.hh"
#include "Offline/CaloReco/inc/CaloWaveformProcessor.hh"
#include "Offline/CaloReco/inc/CaloTemplateWFProcessor.hh"
#include "Offline/CaloReco/inc/CaloRawWFProcessor.hh"

#include "Offline/RecoDataProducts/inc/CaloDigi.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"

#include "TTree.h"
namespace
{
   constexpr int ntupLen = 1638400;
}

int Contains(std::vector<int> v, int x)
{
  return std::count(v.begin(), v.end(), x);
}
namespace mu2e {

    class SourceCalibDigiAna : public art::EDAnalyzer {

      public:
      struct Config
      {
        using Name    = fhicl::Name;
        using Comment = fhicl::Comment;

        fhicl::Atom<art::InputTag> caloDigiCollection{Name("CaloDigiCollection"),Comment("Calo Digi collection name")};
        fhicl::Atom<double> digiSampling{ Name("digiSampling"),Comment("Calo ADC sampling time (ns)")};
        fhicl::Table<mu2e::CaloTemplateWFProcessor::Config> proc_templ_conf     { Name("TemplateProcessor"),   Comment("Log normal fit processor config") };
        fhicl::Atom<int> diagLevel{Name("diagLevel"),Comment("Diag Level"),0};
      };
      explicit SourceCalibDigiAna(const art::EDAnalyzer::Table<Config>& config);
      virtual ~SourceCalibDigiAna() {}

      virtual void beginJob();
      virtual void endJob() {};
      virtual void analyze(const art::Event& e);


      private:

        art::InputTag         caloDigiTag_;
        double digiSampling_;
        int                   diagLevel_;
        int                   nProcess_;

        TTree* Ntup_;
        unsigned int _evt,_run;
        unsigned int ncaldigi_, nCrystals_;
        unsigned int SiPMId_[ntupLen], cryId_[ntupLen],crySectionId_[ntupLen];
        float cryEtot_,calDigiRecoTime_[ntupLen], calDigiRecoADC_[ntupLen], calDigiRecoTimeErr_[ntupLen],calDigiRecoADCErr_[ntupLen];
        std::unique_ptr<CaloWaveformProcessor> waveformProcessor_;
        };


    SourceCalibDigiAna::SourceCalibDigiAna(const art::EDAnalyzer::Table<Config>& config):
      EDAnalyzer{config},
      caloDigiTag_    (config().caloDigiCollection()),
      digiSampling_    (config().digiSampling()),
      diagLevel_          (config().diagLevel()),
      nProcess_(0),
      Ntup_(0)
      {
        waveformProcessor_ = std::make_unique<CaloTemplateWFProcessor>(config().proc_templ_conf());
      
      }

    void SourceCalibDigiAna::beginJob(){

      art::ServiceHandle<art::TFileService> tfs;

      Ntup_  = tfs->make<TTree>("SourceCalibDigiAna","SourceCalibDigiAna");

      Ntup_->Branch("evt",          &_evt ,         "evt/I");
      Ntup_->Branch("run",          &_run ,         "run/I");

      // Reconstructed carystal hit info (from CaloDigiCollection):
      Ntup_->Branch("calDigiRecoEtot",      &cryEtot_ ,     "calDigiRecoEtot/F");
      Ntup_->Branch("ncaldigi",         &ncaldigi_ ,       "ncaldigi/I");
      Ntup_->Branch("nCrystals",         &nCrystals_ ,       "nCrystals/I");
      Ntup_->Branch("cryId",        &cryId_ ,       "cryId[ncaldigi]/I");
      Ntup_->Branch("SiPMId",        &SiPMId_ ,       "SiPMId[ncaldigi]/I");
      Ntup_->Branch("crySectionId", &crySectionId_, "crySectionId[ncaldigi]/I");
      Ntup_->Branch("calDigiRecoADC",      &calDigiRecoADC_ ,     "calDigiRecoADC[ncaldigi]/F");
      Ntup_->Branch("calDigiRecoADCErr",   &calDigiRecoADCErr_ ,  "calDigiRecoADCErr[ncaldigi]/F");
      Ntup_->Branch("calDigiRecoTime",      &calDigiRecoTime_ ,     "calDigiRecoTime[ncaldigi]/F");
      Ntup_->Branch("calDigiRecoTimeErr",   &calDigiRecoTimeErr_ ,  "calDigiRecoTimeErr[ncaldigi]/F");
    }

    void SourceCalibDigiAna::analyze(const art::Event& event){
        ++nProcess_;
        if (nProcess_%10==0 && diagLevel_ > 0) std::cout<<"Processing event from SourceCalibDigiAna =  "<<nProcess_ <<std::endl;

        //Calorimeter crystal hits (average from readouts)
        
        const auto& caloDigisH = event.getValidHandle(caloDigiTag_);
        const auto& caloDigis = *caloDigisH;
        auto pbtH = event.getValidHandle(pbttoken_);
        const ProtonBunchTime& pbt(*pbtH);
        double pbtOffset = pbt.pbtime_;

        _evt = event.id().event();
        _run = event.run();

        //--------------------------  Do calorimeter hits --------------------------------
        ncaldigi_ = 0;
        cryEtot_ = 0.0;

        double totEnergyReco(0);
        std::vector<double> x{},y{};
        std::vector<int> crystalsHit;
        int ic = 0;
        for (const auto& caloDigi : caloDigis)
        {
          int    SiPMID   = caloDigi.SiPMID();
          SiPMId_[ncaldigi_]  = SiPMID; //TODO should we be per siPM
          
          cryId_[ncaldigi_]        =  SiPMID/2;

          if(ic == 0) crystalsHit.push_back(SiPMID/2);
          else if(Contains(crystalsHit, SiPMID/2) == 0) crystalsHit.push_back(SiPMID/2);
          
          crySectionId_[ncaldigi_] = 0;//TODO diskId;
 
          double t0       = caloDigi.t0();
          //calDigiRecoTime_[ncaldigi_]  = t0;
          const std::vector<int>& waveform = caloDigi.waveform();

          size_t index = &caloDigi - &caloDigis.front();
          art::Ptr<CaloDigi> caloDigiPtr(caloDigisHandle, index);

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
              calDigiRecoADC_[ncaldigi_]      = waveformProcessor_->amplitude(i);
              calDigiRecoADCErr_[ncaldigi_]      = waveformProcessor_->amplitudeErr(i);
              calDigiRecoTime_[ncaldigi_] = waveformProcessor_->time(i);
              calDigiRecoTime_[ncaldigi_] = waveformProcessor_->timeErr(i);
              if (SiPMID%2==0) totEnergyReco += waveformProcessor_->amplitude(i);//combine two sipms
          }
          ++ncaldigi_;
          ic += 1;
        }
        
  }
}
DEFINE_ART_MODULE(mu2e::SourceCalibDigiAna)*/
