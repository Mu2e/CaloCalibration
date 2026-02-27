#ifndef CALOCALIBTABLEMAKER_HH
#define CALOCALIBTABLEMAKER_HH

#include <string>
#include <unordered_map>
#include <vector>

#include "CaloCalibration/Combinations/inc/CalEnergyCalibCombiner.hh"

struct ArchiveFitRow {
  int roid = -1;
  double peak = 0.0;
  double errPeak = 0.0;
  double chi2 = 0.0;
  int ndf = 0;
};

struct CalEnergyCalibRow {
  int roid = -1;
  double adc2MeV = 1.0;
};

struct CombinedCalibRow {
  int roid = -1;
  double adc2MeV = 1.0;
  double adc2MeVErr = 0.0;
  int statusCode = static_cast<int>(CalibStatus::kFallbackAllInvalid);
  std::string statusMessage = "FALLBACK: bad fits";
};

class CaloCalibTableMaker {
 public:
  static constexpr int kTotalChannels = 2720;

  std::unordered_map<int, ArchiveFitRow> readCosmicTable(const std::string& inputPath) const;
  std::unordered_map<int, ArchiveFitRow> readSourceTable(const std::string& inputPath) const;
  std::unordered_map<int, CalEnergyCalibRow> readNominalTable(const std::string& inputPath) const;

  std::vector<CombinedCalibRow> combineTables(
      const std::unordered_map<int, CalEnergyCalibRow>& nominal,
      const std::unordered_map<int, ArchiveFitRow>& cosmic,
      const std::unordered_map<int, ArchiveFitRow>& source,
      bool useCosmic,
      bool useSource,
      const CalEnergyCalibCombiner& combiner) const;

  bool writeRecoTable(const std::string& outputPath, const std::vector<CombinedCalibRow>& rows) const;
  bool writeInfoTable(const std::string& outputPath, const std::vector<CombinedCalibRow>& rows) const;
};

#endif
