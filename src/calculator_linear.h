#ifndef ECHMET_LEMNG_CALCULATOR_LINEAR_H
#define ECHMET_LEMNG_CALCULATOR_LINEAR_H

#include "lemng_p.h"
#include "calculator_types.h"
#include <vector>

namespace ECHMET {
namespace LEMNG {
namespace Calculator {

class Eigenzone {
public:
	Eigenzone(std::vector<double> &&constituentConcentrations) noexcept;
	Eigenzone(const double zoneMobility, std::vector<double> &&constituentConcentrations, SolutionProperties &&solutionProperties, const bool tainted, const bool isAnalyzeZone) noexcept;
	Eigenzone(const Eigenzone &other);
	Eigenzone(Eigenzone &&other) noexcept;

	const std::vector<double> constituentConcentrations;
	const SolutionProperties solutionProperties;
	const double zoneMobility;
	const bool tainted;
	const bool isAnalyzeZone;
};

class LinearResults {
public:
	LinearResults(std::vector<Eigenzone> &&eigenzones, const QLQRPack &QLQR, const EMMatrix &M1, const EMMatrix &M2, const bool allZonesValid);
	LinearResults(std::vector<Eigenzone> &&eigenzones, QLQRPack &&QLQR, EMMatrix &&M1, EMMatrix &&M2, const bool allZonesValid) noexcept;
	LinearResults(const LinearResults &other);
	LinearResults(LinearResults &&other) noexcept;

	const std::vector<Eigenzone> eigenzones;
	const QLQRPack QLQR;
	const EMMatrix M1;
	const EMMatrix M2;
	const bool allZonesValid;
};

LinearResults calculateLinear(const CalculatorSystemPack &systemPack, const DeltaPackVec &deltaPacks, const bool correctForIonicStrength);

} // namespace Calculator
} // namespace LEMNG
} // namespace ECHMET


#endif // ECHMET_LEMNG_CALCULATOR_LINEAR_H
