#ifndef ECHMET_LEMNG_CALCULATOR_NONLINEAR_H
#define ECHMET_LEMNG_CALCULATOR_NONLINEAR_H

#include "lemng_p.h"
#include "calculator_types.h"
#include <vector>

namespace ECHMET {
namespace LEMNG {
namespace Calculator {

class EigenzoneDispersion {
public:
	EigenzoneDispersion(const double a2t, const double uEMD);

	const double a2t;	/*!< Time-independent coefficient used to derive the HVL a2 value. */
	const double uEMD;	/*!< Electromigration dispersion velocity slope. */
};
typedef std::vector<EigenzoneDispersion> EigenzoneDispersionVec;

EigenzoneDispersionVec calculateNonlinear(const CalculatorSystemPack &systemPack, const RealVecPtr &analyticalConcentrations,
					  const DeltaPackVec &deltaPacks, const ConcentrationDeltasVec &concentrationDeltasVec,
					  const EMMatrix &M1, const EMMatrix &M2, const QLQRPack &QLQR,
					  const NonidealityCorrections corrections);

} // namespace Calculator
} // namespace LEMNG
} // namespace ECHMET

#endif // ECHMET_LEMNG_CALCULATOR_NONLINEAR_H
