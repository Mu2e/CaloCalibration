#ifndef CaloCalibTableMaker_hh
#define CaloCalibTableMaker_hh
#include<vector>
#include <iostream>
#include <utility>
#include <tuple>
#include <fstream>
using namespace std;

struct ArchiveTable{
  std::vector<int> roid_;
  std::vector<double> peak_;
  std::vector<double> errpeak_;
  std::vector<double> width_;
  std::vector<double> errwidth_;
  std::vector<double> sigma_;
  std::vector<double> errsigma_;
  std::vector<double> chisq_;
  std::vector<int> Nhits_;
  std::string tag_;
  //const char* algtag_;
  ArchiveTable(){};
  ArchiveTable(
  std::vector<int> roid, 
  std::vector<double> peak, 
  std::vector<double> errpeak, 
  std::vector<double> width, 
  std::vector<double> errwidth, 
  std::vector<double> sigma, 
  std::vector<double> errsigma, 
  std::vector<double> chisq,
  std::vector<int> nhits,
  std::string tag) 
  : roid_(roid), peak_(peak), errpeak_(errpeak), width_(width), errwidth_(errwidth), sigma_(sigma), errsigma_(errsigma), chisq_(chisq), Nhits_(nhits), tag_(tag) {};
};

struct RecoTable{
  std::vector<int> roid_;
  std::vector<double> ADC2MeV_;
  RecoTable(){};
  RecoTable(
  std::vector<int> roid, 
  std::vector<double> ADC2MeV ) 
  : roid_(roid), ADC2MeV_(ADC2MeV) {};
};

class CaloCalibTableMaker
	{
      public:
        ArchiveTable inputTable(const char *ArchiveFile, const char *tag);
        RecoTable combineAlg(std::vector<ArchiveTable> tables);
        void WriteOutput(RecoTable table);
  };
#endif

