#include "CaloCalibration/Combinations/inc/CaloCalibTableMaker.hh"

#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace {
constexpr double kCosmicReferenceMeV = 20.0;
constexpr double kSourceReferenceMeV = 6.13;

std::string trim(const std::string& value) {
  std::size_t begin = 0;
  while (begin < value.size() && std::isspace(static_cast<unsigned char>(value[begin]))) {
    ++begin;
  }
  std::size_t end = value.size();
  while (end > begin && std::isspace(static_cast<unsigned char>(value[end - 1]))) {
    --end;
  }
  return value.substr(begin, end - begin);
}

std::vector<std::string> splitCsv(const std::string& line) {
  std::vector<std::string> tokens;
  std::stringstream ss(line);
  std::string token;
  while (std::getline(ss, token, ',')) {
    tokens.push_back(trim(token));
  }
  return tokens;
}

bool isDataLine(const std::string& line) {
  const std::string stripped = trim(line);
  if (stripped.empty()) return false;
  if (stripped.rfind("TABLE", 0) == 0) return false;
  if (stripped[0] == '#') return false;
  return true;
}
}  // namespace

std::unordered_map<int, ArchiveFitRow> CaloCalibTableMaker::readCosmicTable(
    const std::string& inputPath) const {
  std::unordered_map<int, ArchiveFitRow> rows;
  std::ifstream input(inputPath);
  if (!input) {
    std::cerr << "Error: cannot open cosmic table " << inputPath << "\n";
    return rows;
  }

  std::string line;
  while (std::getline(input, line)) {
    if (!isDataLine(line)) continue;
    const auto columns = splitCsv(line);
    if (columns.size() < 9) continue;

    try {
      const int roid = std::stoi(columns[0]);
      rows[roid] = ArchiveFitRow{
          roid,
          std::stod(columns[1]),
          std::stod(columns[2]),
          std::stod(columns[7]),
          std::stoi(columns[8])};
    } catch (const std::exception&) {
    }
  }

  return rows;
}

std::unordered_map<int, ArchiveFitRow> CaloCalibTableMaker::readSourceTable(
    const std::string& inputPath) const {
  std::unordered_map<int, ArchiveFitRow> rows;
  std::ifstream input(inputPath);
  if (!input) {
    std::cerr << "Error: cannot open source table " << inputPath << "\n";
    return rows;
  }

  std::string line;
  while (std::getline(input, line)) {
    if (!isDataLine(line)) continue;
    const auto columns = splitCsv(line);
    // Expected order from CalSourceEnergyCalib:
    // roid,fullepeak,fullerrepeak,...,frsecond,chisq,ndf
    if (columns.size() < 18) continue;

    try {
      const int roid = std::stoi(columns[0]);
      rows[roid] = ArchiveFitRow{
          roid,
          std::stod(columns[1]),
          std::stod(columns[2]),
          std::stod(columns[16]),
          std::stoi(columns[17])};
    } catch (const std::exception&) {
    }
  }

  return rows;
}

std::unordered_map<int, CalEnergyCalibRow> CaloCalibTableMaker::readNominalTable(
    const std::string& inputPath) const {
  std::unordered_map<int, CalEnergyCalibRow> rows;
  std::ifstream input(inputPath);
  if (!input) {
    std::cerr << "Warning: cannot open nominal table " << inputPath
              << ", fallback R0 will default to 1.0\n";
    return rows;
  }

  std::string line;
  while (std::getline(input, line)) {
    if (!isDataLine(line)) continue;
    const auto columns = splitCsv(line);
    if (columns.size() < 2) continue;

    try {
      const int roid = std::stoi(columns[0]);
      rows[roid] = CalEnergyCalibRow{roid, std::stod(columns[1])};
    } catch (const std::exception&) {
    }
  }

  return rows;
}

std::vector<CombinedCalibRow> CaloCalibTableMaker::combineTables(
    const std::unordered_map<int, CalEnergyCalibRow>& nominal,
    const std::unordered_map<int, ArchiveFitRow>& cosmic,
    const std::unordered_map<int, ArchiveFitRow>& source,
    bool useCosmic,
    bool useSource,
    const CalEnergyCalibCombiner& combiner) const {
  std::vector<CombinedCalibRow> out;
  out.reserve(kTotalChannels);

  for (int roid = 0; roid < kTotalChannels; ++roid) {
    const auto nominalIt = nominal.find(roid);
    const double oldADC2MeV = (nominalIt != nominal.end()) ? nominalIt->second.adc2MeV : 1.0;

    double cosmicPeak = 0.0;
    double cosmicErrPeak = -1.0;
    double cosmicChi2 = 0.0;
    int cosmicNdf = 0;
    if (useCosmic) {
      const auto cosmicIt = cosmic.find(roid);
      if (cosmicIt != cosmic.end()) {
        // Cosmic table stores peak in ADC; convert to method calibration constant.
        cosmicPeak = cosmicIt->second.peak / kCosmicReferenceMeV;
        cosmicErrPeak = cosmicIt->second.errPeak / kCosmicReferenceMeV;
        cosmicChi2 = cosmicIt->second.chi2;
        cosmicNdf = cosmicIt->second.ndf;
      }
    }

    double sourcePeak = 0.0;
    double sourceErrPeak = -1.0;
    double sourceChi2 = 0.0;
    int sourceNdf = 0;
    if (useSource) {
      const auto sourceIt = source.find(roid);
      if (sourceIt != source.end()) {
        // Source table stores full-energy peak in ADC; convert similarly.
        sourcePeak = sourceIt->second.peak / kSourceReferenceMeV;
        sourceErrPeak = sourceIt->second.errPeak / kSourceReferenceMeV;
        sourceChi2 = sourceIt->second.chi2;
        sourceNdf = sourceIt->second.ndf;
      }
    }

    const CalibCombineResult result = combiner.combineChannel(
        oldADC2MeV,
        cosmicPeak,
        cosmicErrPeak,
        cosmicChi2,
        cosmicNdf,
        sourcePeak,
        sourceErrPeak,
        sourceChi2,
        sourceNdf);

    out.push_back(CombinedCalibRow{
        roid,
        result.adc2MeV,
        result.adc2MeVErr,
        static_cast<int>(result.status),
        result.statusMessage});
  }

  return out;
}

bool CaloCalibTableMaker::writeRecoTable(
    const std::string& outputPath, const std::vector<CombinedCalibRow>& rows) const {
  std::ofstream output(outputPath);
  if (!output) {
    std::cerr << "Error: cannot write " << outputPath << "\n";
    return false;
  }

  output << "TABLE CalEnergyCalib\n";
  output << std::fixed << std::setprecision(5);
  for (const auto& row : rows) {
    output << row.roid << "," << row.adc2MeV << "\n";
  }
  return output.good();
}

bool CaloCalibTableMaker::writeInfoTable(
    const std::string& outputPath, const std::vector<CombinedCalibRow>& rows) const {
  std::ofstream output(outputPath);
  if (!output) {
    std::cerr << "Error: cannot write " << outputPath << "\n";
    return false;
  }

  output << "TABLE CalEnergyCalibInfo\n";
  output << std::fixed << std::setprecision(5);
  for (const auto& row : rows) {
    output << row.roid << "," << row.adc2MeV << "," << row.adc2MeVErr << ","
           << row.statusCode << "," << row.statusMessage << "\n";
  }
  return output.good();
}

int main(int argc, char** argv) {
  bool useCosmic = true;
  bool useSource = true;
  std::string cosmicPath = "Cosmics.txt";
  std::string sourcePath = "Source.txt";
  std::string nominalPath = "nominal.txt";
  std::string recoOutputPath = "RecoTable.txt";
  std::string infoOutputPath = "CalEnergyCalibInfo.txt";

  if (argc >= 3) {
    useCosmic = (std::stoi(argv[1]) != 0);
    useSource = (std::stoi(argv[2]) != 0);
  }
  if (argc >= 4) cosmicPath = argv[3];
  if (argc >= 5) sourcePath = argv[4];
  if (argc >= 6) nominalPath = argv[5];
  if (argc >= 7) recoOutputPath = argv[6];
  if (argc >= 8) infoOutputPath = argv[7];

  if (!useCosmic && !useSource) {
    std::cerr << "Error: at least one input method must be enabled\n";
    return 1;
  }

  CaloCalibTableMaker maker;
  const auto nominal = maker.readNominalTable(nominalPath);
  const auto cosmic = useCosmic ? maker.readCosmicTable(cosmicPath)
                                : std::unordered_map<int, ArchiveFitRow>{};
  const auto source = useSource ? maker.readSourceTable(sourcePath)
                                : std::unordered_map<int, ArchiveFitRow>{};

  CalEnergyCalibCombiner::Config config;
  config.fitPvalueThreshold = 0.01;
  config.compatPvalueThreshold = 0.05;
  const CalEnergyCalibCombiner combiner(config);

  const auto combined = maker.combineTables(
      nominal, cosmic, source, useCosmic, useSource, combiner);

  if (!maker.writeRecoTable(recoOutputPath, combined)) return 1;
  if (!maker.writeInfoTable(infoOutputPath, combined)) return 1;

  std::cout << "Wrote " << combined.size() << " rows to " << recoOutputPath
            << " and " << infoOutputPath << "\n";
  return 0;
}
