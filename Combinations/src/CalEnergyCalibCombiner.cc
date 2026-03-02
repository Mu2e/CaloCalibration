#include "CaloCalibration/Combinations/inc/CalEnergyCalibCombiner.hh"

#include <cmath>

namespace {
double incompleteGamma(double a, double x) {
  if (x < 0.0 || a <= 0.0) return 0.0;
  if (x == 0.0) return 0.0;

  if (x < a + 1.0) {
    double ap = a;
    double sum = 1.0 / a;
    double del = sum;
    for (int n = 1; n <= 100; ++n) {
      ap += 1.0;
      del *= x / ap;
      sum += del;
      if (std::abs(del) < std::abs(sum) * 1e-10) {
        return sum * std::exp(-x + a * std::log(x) - std::lgamma(a));
      }
    }
  } else {
    double b = x + 1.0 - a;
    double c = 1.0 / 1e-30;
    double d = 1.0 / b;
    double h = d;
    for (int n = 1; n <= 100; ++n) {
      double an = -n * (n - a);
      b += 2.0;
      d = an * d + b;
      if (std::abs(d) < 1e-30) d = 1e-30;
      c = b + an / c;
      if (std::abs(c) < 1e-30) c = 1e-30;
      d = 1.0 / d;
      double del = d * c;
      h *= del;
      if (std::abs(del - 1.0) < 1e-10) {
        return 1.0 - std::exp(-x + a * std::log(x) - std::lgamma(a)) * h;
      }
    }
  }

  return (x < a + 1.0) ? 0.0 : 1.0;
}

double chi2CDF(double x, int ndf) {
  if (x <= 0.0 || ndf <= 0) return 0.0;
  return incompleteGamma(ndf / 2.0, x / 2.0);
}
}  // namespace

double CalEnergyCalibCombiner::computePvalue(double chi2, int ndf) const {
  if (ndf <= 0 || chi2 < 0.0) return 0.0;
  return 1.0 - chi2CDF(chi2, ndf);
}

double CalEnergyCalibCombiner::computeCompatPvalue(
    double mu1, double sigma1, double mu2, double sigma2) const {
  const double varianceSum = sigma1 * sigma1 + sigma2 * sigma2;
  if (varianceSum <= 0.0) return 0.0;
  const double delta = mu1 - mu2;
  const double chi2Compat = delta * delta / varianceSum;
  return 1.0 - chi2CDF(chi2Compat, 1);
}

std::pair<double, double> CalEnergyCalibCombiner::weightedAverage(
    double mu1, double sigma1, double mu2, double sigma2) const {
  const double w1 = 1.0 / (sigma1 * sigma1);
  const double w2 = 1.0 / (sigma2 * sigma2);
  const double weightSum = w1 + w2;
  return {(w1 * mu1 + w2 * mu2) / weightSum, std::sqrt(1.0 / weightSum)};
}

std::string CalEnergyCalibCombiner::statusToMessage(CalibStatus status) {
  switch (status) {
    case CalibStatus::kCombinedConsistent:
      return "OK: C+S";
    case CalibStatus::kFallbackInconsistent:
      return "FALLBACK: inconsistent";
    case CalibStatus::kFallbackAllInvalid:
      return "FALLBACK: bad fits";
    case CalibStatus::kCosmicOnly:
      return "OK: cosmic only";
    case CalibStatus::kSourceOnly:
      return "OK: source only";
    default:
      return "UNKNOWN";
  }
}

CalibCombineResult CalEnergyCalibCombiner::combineChannel(
    double oldADC2MeV,
    double cosmicPeak,
    double cosmicErrPeak,
    double cosmicChi2,
    int cosmicNdf,
    double sourcePeak,
    double sourceErrPeak,
    double sourceChi2,
    int sourceNdf) const {
  CalibCombineResult result;

  const bool cosmicValid = (cosmicErrPeak > 0.0);
  const bool sourceValid = (sourceErrPeak > 0.0);

  const double pCosmic = cosmicValid ? computePvalue(cosmicChi2, cosmicNdf) : 0.0;
  const double pSource = sourceValid ? computePvalue(sourceChi2, sourceNdf) : 0.0;

  const bool cosmicPass = cosmicValid && (pCosmic > config_.fitPvalueThreshold);
  const bool sourcePass = sourceValid && (pSource > config_.fitPvalueThreshold);

  if (!cosmicPass && !sourcePass) {
    result.adc2MeV = oldADC2MeV;
    result.adc2MeVErr = 0.0;
    result.status = CalibStatus::kFallbackAllInvalid;
    result.statusMessage = statusToMessage(result.status);
    return result;
  }

  if (cosmicPass && !sourcePass) {
    result.adc2MeV = cosmicPeak;
    result.adc2MeVErr = cosmicErrPeak;
    result.status = CalibStatus::kCosmicOnly;
    result.statusMessage = statusToMessage(result.status);
    return result;
  }

  if (!cosmicPass && sourcePass) {
    result.adc2MeV = sourcePeak;
    result.adc2MeVErr = sourceErrPeak;
    result.status = CalibStatus::kSourceOnly;
    result.statusMessage = statusToMessage(result.status);
    return result;
  }

  const auto [muCombined, sigmaCombined] =
      weightedAverage(cosmicPeak, cosmicErrPeak, sourcePeak, sourceErrPeak);
  const double pCompat =
      computeCompatPvalue(cosmicPeak, cosmicErrPeak, sourcePeak, sourceErrPeak);

  if (pCompat < config_.compatPvalueThreshold) {
    result.adc2MeV = oldADC2MeV;
    result.adc2MeVErr = 0.0;
    result.status = CalibStatus::kFallbackInconsistent;
    result.statusMessage = statusToMessage(result.status);
    return result;
  }

  result.adc2MeV = muCombined;
  result.adc2MeVErr = sigmaCombined;
  result.status = CalibStatus::kCombinedConsistent;
  result.statusMessage = statusToMessage(result.status);
  return result;
}
