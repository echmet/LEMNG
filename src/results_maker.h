#ifndef ECHMET_LEMNG_RESULTS_MAKER_H
#define ECHMET_LEMNG_RESULTS_MAKER_H

#include "lemng.h"
#include "calculator_linear.h"
#include "calculator_nonlinear.h"

namespace ECHMET {
namespace LEMNG {

	void fillResults(const ChemicalSystemPtr &chemSystemBGE, const ChemicalSystemPtr &chemSystemFull, const Calculator::SolutionProperties &BGEProperties, const Calculator::SolutionProperties &BGELikeProperties, const Calculator::LinearResults &linResults, const Calculator::EigenzoneDispersionVec &ezDisps, const NonidealityCorrections corrections, Results &r);
	void fillResultsBGE(const ChemicalSystemPtr &chemSystemBGE, const Calculator::SolutionProperties &BGEProperties, const NonidealityCorrections corrections, Results &r);
	void fillResultsAnalytesDissociation(const ChemicalSystemPtr &chemSystemFull, const Calculator::SolutionProperties &BGELikeProperties, Results &r);
	Results prepareResults(const ChemicalSystemPtr &chemSystemBGE, const ChemicalSystemPtr &chemSystemFull, IsAnalyteFunc &isAnalyte);

} // namespace LEMNG
} // namespace ECHMET

#endif // ECHMET_LEMNG_RESULTS_MAKER_H
