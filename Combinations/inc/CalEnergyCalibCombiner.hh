#ifndef CALENERGYCALIBCOMBINER_HH
#define CALENERGYCALIBCOMBINER_HH

#include <string>
#include <utility>

enum class CalibStatus : int {
  kCombinedConsistent = 0,
  kFallbackInconsistent = 1,
  kFallbackAllInvalid = 2,
  kCosmicOnly = 101,
  kSourceOnly = 102
};

struct CalibCombineResult {
  double adc2MeV = 1.0;
  double adc2MeVErr = 0.0;
  CalibStatus status = CalibStatus::kFallbackAllInvalid;
  std::string statusMessage = "FALLBACK: bad fits";
};

class CalEnergyCalibCombiner {
 public:
  struct Config {
    double fitPvalueThreshold;
    double compatPvalueThreshold;
    Config() : fitPvalueThreshold(0.01), compatPvalueThreshold(0.05) {}
  };

  CalEnergyCalibCombiner() : config_() {}
  explicit CalEnergyCalibCombiner(const Config& config) : config_(config) {}

  CalibCombineResult combineChannel(
      double oldADC2MeV,
      double cosmicPeak,
      double cosmicErrPeak,
      double cosmicChi2,
      int cosmicNdf,
      double sourcePeak,
      double sourceErrPeak,
      double sourceChi2,
      int sourceNdf) const;

 private:
  Config config_;

  double computePvalue(double chi2, int ndf) const;
  double computeCompatPvalue(double mu1, double sigma1, double mu2, double sigma2) const;
  std::pair<double, double> weightedAverage(double mu1, double sigma1, double mu2, double sigma2) const;
  static std::string statusToMessage(CalibStatus status);
};

#endif
