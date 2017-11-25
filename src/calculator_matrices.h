#ifndef ECHMET_LEMNG_CALCULATOR_MATRICES_H
#define ECHMET_LEMNG_CALCULATOR_MATRICES_H

#include "calculator_types.h"
#include "lemng_p.h"

namespace ECHMET {

/* Forward-declarations of CAES types */
namespace CAES {
	class Solver;
} // namespace CAES

namespace LEMNG {
namespace Calculator {

typedef std::vector<EMMatrix> EMMatrixVec;

CalculatorSystemPack makeIonicFormPack(const ChemicalSystemPtr &BGESystem, const CalculatedPropertiesPtr &BGEcalcProps,
				       const ChemicalSystemPtr &fullSystem, CalculatedPropertiesPtr &fullCalcProps,
				       const RealVecPtr &analConcsBGE,
				       const std::function<bool (const std::string &s)> &isAnalyte,
				       const bool correctForIonicStrength);
EMMatrix makeMatrixM1(const CalculatorSystemPack &systemPack);
EMMatrix makeMatrixM2(const CalculatorSystemPack &systemPack, const DeltaPackVec &deltaPacks);
EMMatrix makeM1Derivative(const CalculatorSystemPack &systemPack, const DeltaPack &deltaPack) noexcept;
EMMatrix makeM2Derivative(const CalculatorSystemPack &systemPack, const RealVecPtr &analyticalConcentrations, const SysComp::Constituent *pivotalConstituent, CAES::Solver *solver, RealVec *derivatives);

} // namespace Calculator
} // namespace LEMNG
} // namespace ECHMET

#endif // ECHMET_LEMNG_CALCULATOR_MATRICES_H
