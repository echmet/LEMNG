#ifndef ECHMET_LEMNG_RESULTS_MAKER_H
#define ECHMET_LEMNG_RESULTS_MAKER_H

#include "lemng.h"
#include "calculator_linear.h"
#include "calculator_nonlinear.h"

namespace ECHMET {
namespace LEMNG {

	void fillResults(const ChemicalSystemPtr &chemSystemBGE, const ChemicalSystemPtr &chemSystemFull, const Calculator::SolutionProperties &BGEProperties, const Calculator::LinearResults &linResults, const Calculator::EigenzoneDispersionVec &ezDisps, const bool correctForIonicStrength, Results &r);
	void fillResultsPartial(const ChemicalSystemPtr &chemSystemBGE, const Calculator::SolutionProperties &BGEProperties, const bool correctForIonicStrength, Results &r);
	Results prepareResults(const ChemicalSystemPtr &chemSystemBGE, const ChemicalSystemPtr &chemSystemFull);

} // namespace LEMNG
} // namespace ECHMET

#endif // ECHMET_LEMNG_RESULTS_MAKER_H
