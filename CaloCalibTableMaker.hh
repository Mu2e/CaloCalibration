#ifndef CaloCalibTableMaker_hh
#define CaloCalibTableMaker_hh
#include<vector>
#include <iostream>
#include <utility>
#include <tuple>
using namespace std;

struct ArchiveTable{
  std::vector<int> roid_;
  std::vector<double> peak_;
  std::vector<double> errpeak_;
  std::vector<double> width_;
  std::vector<double> errwidth_;
  std::vector<double> chisq_;
  const char* algtag_;
  ArchiveTable(){};
  ArchiveTable(
  std::vector<int> roid, 
  std::vector<double> peak, 
  std::vector<double> errpeak, 
  std::vector<double> width, 
  std::vector<double> errwidth, 
  std::vector<double> chisq, 
  const char* tag) 
  : roid_(roid), peak_(peak), errpeak_(errpeak), width_(width), errwidth_(errwidth), chisq_(chisq), algtag_(tag) {};
};

struct RecoTable{
  std::vector<int> roid_;
  std::vector<double> peak_;
  std::vector<double> errpeak_;
  int algtag_;
  RecoTable(){};
  RecoTable(
  std::vector<int> roid, 
  std::vector<double> peak, 
  std::vector<double> errpeak, 
  int tag) 
  : roid_(roid), peak_(peak), errpeak_(errpeak), algtag_(tag) {};
};

class CaloCalibTableMaker
	{
      public:
        ArchiveTable inputTable(const char *ArchiveFile, const char *tag);
        RecoTable combineAlg(std::vector<ArchiveTable> tables);
        void WriteOutput(RecoTable table);
  };
#endif

