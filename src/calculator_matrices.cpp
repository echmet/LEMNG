#include "calculator_matrices.h"
#include "calculator_common.h"
#include "helpers.h"

#ifndef ECHMET_IMPORT_INTERNAL
#define ECHMET_IMPORT_INTERNAL
#endif // ECHMET_IMPORT_INTERNAL
#include <echmetphchconsts.h>
#include <echmetsyscomp.h>
#include <echmetcaes_extended.h>

#include "tracing/lemng_tracer_impl.h"
#include <sstream>

namespace ECHMET {
namespace LEMNG {
namespace Calculator {

static
int multiplicity(const size_t i, const MultiplicityVec &mList)
{
	for (const auto &m : mList)
		if (i == m.first) return m.second;

	return 0;
}

EMMatrix makeM1Derivative(const CalculatorSystemPack &systemPack, const DeltaPack &deltaPack) noexcept
{
	const size_t ROWS = systemPack.constituents.size();
	const size_t COLS = systemPack.ionicForms.size();
	const size_t H3O_idx = COLS - 2;
	const size_t OH_idx = COLS - 1;
	EMMatrix MOneDer{ROWS, COLS};

	const CalculatorConstituentVec &ccVec = systemPack.constituents;
	const CalculatorIonicFormVec &ifVec = systemPack.ionicForms;
	const double baseConductivity = systemPack.conductivity * 1.0e9;
	const double dKdcJ = deltaPack.conductivityDelta * 1.0e9;

	ECHMET_TRACE(LEMNGTracing, CALC_MATRIX_DM1_INPUT, baseConductivity, std::cref(deltaPack.concentrationDeltas));

	for (size_t row = 0; row < ROWS; row++) {
		const CalculatorConstituent &c = ccVec.at(row);

		const double uIcIdKdC = [&c, baseConductivity, dKdcJ]() {
			double s = 0.0;

			for (const CalculatorIonicForm *iF : c.ionicForms) {
				const int d = (iF->internalIonicForm->ligand != nullptr) ? iF->internalIonicForm->ligandCount : 1;
				s += iF->concentration * iF->mobility * cxsgn(iF->charge) * d;
				ECHMET_TRACE(LEMNGTracing, CALC_MATRIX_DM1_UICIDKDC, s, iF->mobility, cxsgn(iF->charge), d);
			}

			return s * 2.0 / std::pow(baseConductivity, 2) * dKdcJ;
		}();

		const double uIdcIdcJ = [&c, &deltaPack, baseConductivity]() {
			double s = 0.0;

			for (const CalculatorIonicForm *iF : c.ionicForms) {
				const int d = (iF->internalIonicForm->ligand != nullptr) ? iF->internalIonicForm->ligandCount : 1;
				s += iF->mobility * cxsgn(iF->charge) * deltaPack.concentrationDeltas(iF->globalIonicFormConcentrationIdx) * d;
				ECHMET_TRACE(LEMNGTracing, CALC_MATRIX_DM1_UIDCIDCJ, std::cref(iF->name), iF->globalIonicFormConcentrationIdx,
					     deltaPack.concentrationDeltas(iF->globalIonicFormConcentrationIdx), s, d);
			}

			return -s / baseConductivity;
		}();

		const double termTwo = (uIcIdKdC + uIdcIdcJ) * PhChConsts::F;

		for (size_t col = 0; col < COLS - 2; col++) {
			const CalculatorIonicForm *iF = ifVec.at(col);
			const auto &multiplicities = iF->multiplicities;
			const int32_t charge = iF->charge;
			const double mobility = iF->mobility;
			const int d = multiplicity(row, multiplicities);

			double result = -d * cxsgn(charge) * dKdcJ / baseConductivity;
			result += termTwo * std::abs(charge);
			result *= mobility;

			ECHMET_TRACE(LEMNGTracing, CALC_MATRIX_DM1_KRD, d);

			MOneDer(row, col) = result;
		}

		const double H3OMob = ifVec.at(H3O_idx)->mobility;
		const double OHMob = ifVec.at(OH_idx)->mobility;

		ECHMET_TRACE(LEMNGTracing, CALC_MATRIX_DM1_ROW_BLOCK, uIcIdKdC, uIdcIdcJ, termTwo, H3OMob, OHMob);

		/* H3O+ and OH- are in the last two columns */
		MOneDer(row, H3O_idx) = termTwo * H3OMob;
		MOneDer(row, OH_idx) = termTwo * OHMob;
	}

	ECHMET_TRACE(LEMNGTracing, CALC_MATRIX_DM1_OUTPUT, std::cref(MOneDer));

	return MOneDer;
}

EMMatrix makeM2Derivative(const CalculatorSystemPack &systemPack, const RealVecPtr &analyticalConcentrations, const SysComp::Constituent *pivotalConstituent, CAES::Solver *solver, RealVec *derivatives)
{
	const ECHMETReal H = DELTA_H;

	const size_t ROWS = systemPack.ionicForms.size();
	const size_t COLS = systemPack.constituents.size();
	const size_t H3O_idx = ROWS - 2;
	const size_t OH_idx = ROWS - 1;
#ifdef ECHMET_LEMNG_SENSITIVE_NUMDERS
	(void)analyticalConcentrations;
	const RealVecPtr analyticalConcentrationsForDiffs = makeAnalyticalConcentrationsForDerivator(systemPack);
#else
	const RealVecPtr &analyticalConcentrationsForDiffs = analyticalConcentrations;
#endif
	EMMatrix MTwoDer{ROWS, COLS};

	SysComp::ChemicalSystem chemSystemRaw = *systemPack.chemSystemRaw;
	const SysComp::CalculatedProperties *calcPropsRaw = systemPack.calcPropsRaw;

	for (size_t col = 0; col < COLS; col++) {
		const SysComp::Constituent *cK = systemPack.constituents.at(col).internalConstituent;

		::ECHMET::RetCode tRet = CAES::calculateCrossConcentrationDerivatives_prepared(derivatives, solver, H, chemSystemRaw, analyticalConcentrationsForDiffs.get(), pivotalConstituent, cK, calcPropsRaw->ionicStrength);
		if (tRet != ::ECHMET::RetCode::OK)
			throw CalculationException{"Cannot calculate concentration derivatives for M2 derivative", coreLibsErrorToNativeError(tRet)};

		for (size_t row = 0; row < ROWS - 2; row++) {
			const CalculatorIonicForm *iF = systemPack.ionicForms.at(row);

			MTwoDer(row, col) = ECHMETRealToDouble(derivatives->at(iF->internalIonicFormConcentrationIdx));
		}

		MTwoDer(H3O_idx, col) = ECHMETRealToDouble(derivatives->at(0));
		MTwoDer(OH_idx, col) = ECHMETRealToDouble(derivatives->at(1));
	}

	ECHMET_TRACE(LEMNGTracing, CALC_MATRIX_DM2_OUTPUT, std::cref(MTwoDer));

	return MTwoDer;
}

EMMatrix makeMatrixD1(const CalculatorSystemPack &systemPack, const ERVector &diffusionCoefficients)
{
	const size_t ROWS = systemPack.constituents.size();
	const size_t COLS = systemPack.ionicForms.size();
	const size_t H3O_idx = COLS - 2;
	const size_t OH_idx = COLS - 1;
	EMMatrix DOne{ROWS, COLS};

	const CalculatorConstituentVec &ccVec = systemPack.constituents;
	const CalculatorIonicFormVec &ifVec = systemPack.ionicForms;

	ECHMET_TRACE(LEMNGTracing, CALC_MATRIX_D1_DIMS, ROWS, COLS);

	for (size_t row = 0; row < ROWS; row++) {
		const CalculatorConstituent &c = ccVec.at(row);

		const double uIcIFSum = [](const CalculatorConstituent &c, const double conductivity) {
			double s = 0.0;

			ECHMET_TRACE(LEMNGTracing, CALC_MATRIX_D1_UICIF_BLOCK, c.name.c_str());

			for (size_t idx = 0; idx < c.ionicForms.size(); idx++) {
				const CalculatorIonicForm *iF = c.ionicForms.at(idx);

				const int d = (iF->internalIonicForm->ligand != nullptr) ? iF->internalIonicForm->ligandCount : 1;
				s += iF->concentration * iF->mobility * cxsgn(iF->charge) * d;
				ECHMET_TRACE(LEMNGTracing, CALC_MATRIX_D1_UICIF_INTERMEDIATE, s, iF, d);
			}

			return s * PhChConsts::F / (conductivity * 1.0e9);
		}(c, systemPack.conductivity);

		for (size_t col = 0; col < COLS - 2; col++) {
			const CalculatorIonicForm *iF = ifVec.at(col);
			const auto &multiplicities = iF->multiplicities;
			const double diffCoeff = diffusionCoefficients.at(col);
			const int d = multiplicity(row, multiplicities);

			ECHMET_TRACE(LEMNGTracing, CALC_MATRIX_D1_ROW_BLOCK, std::cref(iF->name), d, col, diffCoeff);

			DOne(row, col) = (d - uIcIFSum * iF->charge) * diffCoeff;
		}

		/* H3O+ and OH- are in the last two columns */
		DOne(row, H3O_idx) = uIcIFSum * diffusionCoefficients.at(H3O_idx);
		DOne(row, OH_idx) = uIcIFSum * diffusionCoefficients.at(OH_idx);
	}

	ECHMET_TRACE(LEMNGTracing, CALC_MATRIX_D1_OUTPUT, std::cref(DOne), std::cref(systemPack.constituents), std::cref(systemPack.ionicForms));

	return DOne;
}

EMMatrix makeMatrixD2(const CalculatorSystemPack &systemPack, const DeltaPackVec &deltaPacks)
{
	const size_t ROWS = systemPack.ionicForms.size();
	const size_t COLS = systemPack.constituents.size();
	EMMatrix DTwo{ROWS, COLS};

	ECHMET_TRACE(LEMNGTracing, CALC_MATRIX_D2_DIMS, ROWS, COLS);

	for (size_t col = 0; col < COLS; col++)
		DTwo.col(col) = deltaPacks.at(col).concentrationDeltas;


	ECHMET_TRACE(LEMNGTracing, CALC_MATRIX_D2_OUTPUT, std::cref(DTwo), std::cref(systemPack.constituents), std::cref(systemPack.ionicForms));

	return DTwo;
}

EMMatrix makeMatrixM1(const CalculatorSystemPack &systemPack)
{
	const size_t ROWS = systemPack.constituents.size();
	const size_t COLS = systemPack.ionicForms.size();
	const size_t H3O_idx = COLS - 2;
	const size_t OH_idx = COLS - 1;
	EMMatrix MOne{ROWS, COLS};

	const CalculatorConstituentVec &ccVec = systemPack.constituents;
	const CalculatorIonicFormVec &ifVec = systemPack.ionicForms;

	ECHMET_TRACE(LEMNGTracing, CALC_MATRIX_M1_DIMS, ROWS, COLS);

	for (size_t row = 0; row < ROWS; row++) {
		const CalculatorConstituent &c = ccVec.at(row);

		const double uIcIFSum = [](const CalculatorConstituent &c, const double conductivity) {
			double s = 0.0;

			ECHMET_TRACE(LEMNGTracing, CALC_MATRIX_M1_UICIF_BLOCK, c.name.c_str());

			for (size_t idx = 0; idx < c.ionicForms.size(); idx++) {
				const CalculatorIonicForm *iF = c.ionicForms.at(idx);

				const int d = (iF->internalIonicForm->ligand != nullptr) ? iF->internalIonicForm->ligandCount : 1;
				s += iF->concentration * iF->mobility * cxsgn(iF->charge) * d;
				ECHMET_TRACE(LEMNGTracing, CALC_MATRIX_M1_UICIF_INTERMEDIATE, s, iF, d);
			}

			return s * PhChConsts::F / (conductivity * 1.0e9);
		}(c, systemPack.conductivity);

		for (size_t col = 0; col < COLS - 2; col++) {
			const CalculatorIonicForm *iF = ifVec.at(col);
			const auto &multiplicities = iF->multiplicities;
			const int32_t charge = iF->charge;
			const double mobility = iF->mobility;
			const int d = multiplicity(row, multiplicities);

			ECHMET_TRACE(LEMNGTracing, CALC_MATRIX_M1_ROW_BLOCK, std::cref(c.name), std::cref(iF->name), d, col, mobility);

			MOne(row, col) = (d * cxsgn(charge) - uIcIFSum * std::abs(charge)) * mobility;
		}

		/* H3O+ and OH- are in the last two columns */
		MOne(row, H3O_idx) = -uIcIFSum * ifVec.at(H3O_idx)->mobility;
		MOne(row, OH_idx) = -uIcIFSum * ifVec.at(OH_idx)->mobility;
	}

	ECHMET_TRACE(LEMNGTracing, CALC_MATRIX_M1_OUTPUT, std::cref(MOne), std::cref(systemPack.constituents), std::cref(systemPack.ionicForms));

	return MOne;
}

EMMatrix makeMatrixM2(const CalculatorSystemPack &systemPack, const DeltaPackVec &deltaPacks)
{
	const size_t ROWS = systemPack.ionicForms.size();
	const size_t COLS = systemPack.constituents.size();
	EMMatrix MTwo{ROWS, COLS};

	ECHMET_TRACE(LEMNGTracing, CALC_MATRIX_M2_DIMS, ROWS, COLS);

	for (size_t col = 0; col < COLS; col++)
		MTwo.col(col) = deltaPacks.at(col).concentrationDeltas;


	ECHMET_TRACE(LEMNGTracing, CALC_MATRIX_M2_OUTPUT, std::cref(MTwo), std::cref(systemPack.constituents), std::cref(systemPack.ionicForms));

	return MTwo;
}

} // namespace Calculator
} // namespace LEMNG

#ifndef ECHMET_TRACER_DISABLE_TRACING

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_MATRIX_M1_DIMS, "Matrix M1 dimensions")
ECHMET_BEGIN_MAKE_LOGGER(LEMNGTracing, CALC_MATRIX_M1_DIMS, const size_t rows, const size_t cols)
{
	std::ostringstream ss{};

	ss << "Calculating Matrix M1(" << rows << ", " << cols << ")";
	return ss.str();
}
ECHMET_END_MAKE_LOGGER

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_MATRIX_M2_DIMS, "Matrix M2 dimensions")
ECHMET_BEGIN_MAKE_LOGGER(LEMNGTracing, CALC_MATRIX_M2_DIMS, const size_t rows, const size_t cols)
{
	std::ostringstream ss{};

	ss << "Calculating Matrix M2(" << rows << ", " << cols << ")";
	return ss.str();
}
ECHMET_END_MAKE_LOGGER

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_MATRIX_M1_UICIF_BLOCK, "Matrix M1 uIcIF intermediate block beginning")
ECHMET_BEGIN_MAKE_LOGGER(LEMNGTracing, CALC_MATRIX_M1_UICIF_BLOCK, const char *block)
{
	std::ostringstream ss{};

	ss << "Entering Matrix M1 uIcIF block " << block;
	return ss.str();
}
ECHMET_END_MAKE_LOGGER

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_MATRIX_M1_UICIF_INTERMEDIATE, "Matrix M1 uIcIF intermediate block output")
ECHMET_BEGIN_MAKE_LOGGER(LEMNGTracing, CALC_MATRIX_M1_UICIF_INTERMEDIATE, const double s, const ECHMET::LEMNG::Calculator::CalculatorIonicForm *ccIf, const int d)
{
	std::ostringstream ss{};
	const std::string &name = ccIf->name;

	ss << "Entering Matrix M1 uIcIF intermediate(" << name << "):\n"
	   << "s(current) = " << s << "\n"
	   << "concentration = " << ccIf->concentration << "\n"
	   << "totalCharge = " << ccIf->charge << "\n"
	   << "mobility = " << ccIf->mobility << "\n"
	   << "ligand " << [](const auto iF) {
		   if (iF->ligand == nullptr)
			   return "(NO LIGAND)";
		   else
			   return iF->ligand->name->c_str();
	   }(ccIf->internalIonicForm)
	   << " multiplicity = " << d << "\n";

	return ss.str();
}
ECHMET_END_MAKE_LOGGER

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_MATRIX_M1_ROW_BLOCK, "Matrix M1 row block")
ECHMET_BEGIN_MAKE_LOGGER(LEMNGTracing, CALC_MATRIX_M1_ROW_BLOCK, const std::string &cName, const std::string &iFName, const int mul, const size_t col, const double mobility)
{
	std::ostringstream ss{};

	ss << "M1 row block[" << cName << "], multiplicity = " << mul << ", col = " << col << "[" << iFName << "], mobility = " << mobility;

	return ss.str();
}
ECHMET_END_MAKE_LOGGER

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_MATRIX_M1_OUTPUT, "Matrix M1 output")
ECHMET_BEGIN_MAKE_LOGGER(LEMNGTracing, CALC_MATRIX_M1_OUTPUT, const ECHMET::LEMNG::Calculator::EMMatrix &MOne, const ECHMET::LEMNG::Calculator::CalculatorConstituentVec &ccVec, const ECHMET::LEMNG::Calculator::CalculatorIonicFormVec &cIfVec)
{
	std::ostringstream ss{};

	ss << "-- Matrix M1 --\n";
	ss << "Columns -> ";
	for (auto *iF : cIfVec)
		ss << iF->name << "; ";
	ss << "\nRows -> ";
	for (auto &cc : ccVec)
		ss << cc.name << "; ";
	ss << "\n";

	ss << "---\n\n" << MOne << "\n\n---";

	return ss.str();
}
ECHMET_END_MAKE_LOGGER

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_MATRIX_M2_OUTPUT, "Matrix M2 output")
ECHMET_BEGIN_MAKE_LOGGER(LEMNGTracing, CALC_MATRIX_M2_OUTPUT, const ECHMET::LEMNG::Calculator::EMMatrix &MTwo, const ECHMET::LEMNG::Calculator::CalculatorConstituentVec &ccVec, const ECHMET::LEMNG::Calculator::CalculatorIonicFormVec &cIfVec)
{
	std::ostringstream ss{};

	ss << "-- Matrix M2 --\n";
	ss << "Rows -> ";
	for (auto &cc : ccVec)
		ss << cc.name << "; ";
	ss << "\nColumns ->";
	for (auto *iF : cIfVec)
		ss << iF->name << "; ";

	ss << "\n";

	ss << "---\n\n" << MTwo << "\n\n---";

	return ss.str();
}
ECHMET_END_MAKE_LOGGER

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_MATRIX_DM1_INPUT, "dM1/dC input")
ECHMET_BEGIN_MAKE_LOGGER(LEMNGTracing, CALC_MATRIX_DM1_INPUT, const double baseConductivity, const ECHMET::LEMNG::Calculator::EMVector &cDeltas)
{
	std::ostringstream ss{};

	ss << "dM1/dC imput\n"
	   << "base conductivity = " << baseConductivity << " (S/m)\n"
	   << "concentrations deltas (mmol/dm3)\n"
	   << cDeltas;

	return ss.str();
}
ECHMET_END_MAKE_LOGGER

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_MATRIX_DM1_UICIDKDC, "dM1/dC dConductivity/dC")
ECHMET_BEGIN_MAKE_LOGGER(LEMNGTracing, CALC_MATRIX_DM1_UICIDKDC, const double s, const double mobility, const int32_t absCharge, const int d)
{
	std::ostringstream ss{};

	ss << "dM1/dC dConductivity/dC\n"
	   << "s = " << s << ", mobility = " << mobility << ", charge = " << absCharge << ", multiplicity = " << d;

	return ss.str();
}
ECHMET_END_MAKE_LOGGER

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_MATRIX_DM1_UIDCIDCJ, "dM1/dC dCi/dCj")
ECHMET_BEGIN_MAKE_LOGGER(LEMNGTracing, CALC_MATRIX_DM1_UIDCIDCJ, const std::string &name, const size_t gIdx, const double cDelta, const double s, const double d)
{
	std::ostringstream ss{};

	ss << "dM1/dC dCi/dCj\n"
	   << name << ", Global IFIdx = " << gIdx << ", concDelta = " << cDelta << ", s = " << s << ", multiplicity = " << d;

	return ss.str();
}
ECHMET_END_MAKE_LOGGER

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_MATRIX_DM1_KRD, "dM1/dC Kroenecker delta")
ECHMET_BEGIN_MAKE_LOGGER(LEMNGTracing, CALC_MATRIX_DM1_KRD, const int KrD)
{
	return std::string{"Kroenecker delta = "} + std::to_string(KrD);
}
ECHMET_END_MAKE_LOGGER

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_MATRIX_DM1_ROW_BLOCK, "dM1/dC row block")
ECHMET_BEGIN_MAKE_LOGGER(LEMNGTracing, CALC_MATRIX_DM1_ROW_BLOCK, const double uIcIcKdC, const double uIdcIdcJ, const double termTwo, const double H3OMob, const double OHMob)
{
	std::ostringstream ss{};

	ss << "dM1/dC row block\n"
	   << "uIcIcKcD = " << uIcIcKdC << ", uIdcIdcJ = " << uIdcIdcJ << "\n"
	   << "termTwo = " << termTwo
	   << "water mobs = " << H3OMob << ", " << OHMob;

	return ss.str();
}
ECHMET_END_MAKE_LOGGER

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_MATRIX_DM1_OUTPUT, "dM1/dC output")
ECHMET_BEGIN_MAKE_LOGGER(LEMNGTracing, CALC_MATRIX_DM1_OUTPUT, const ECHMET::LEMNG::Calculator::EMMatrix &MOneDer)
{
	std::ostringstream ss{};

	ss << "-- Matrix dM1/dC --\n"
	   << "---\n\n" << MOneDer << "\n\n---";

	return ss.str();
}
ECHMET_END_MAKE_LOGGER

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_MATRIX_DM2_OUTPUT, "dM2/dC output")
ECHMET_BEGIN_MAKE_LOGGER(LEMNGTracing, CALC_MATRIX_DM2_OUTPUT, const ECHMET::LEMNG::Calculator::EMMatrix &MTwoDer)
{
	std::ostringstream ss{};

	ss << "-- Matrix dM2/dC --\n"
	   << "---\n\n" << MTwoDer << "\n\n---";

	return ss.str();
}
ECHMET_END_MAKE_LOGGER

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_MATRIX_D1_DIMS, "Matrix D1 dimensions")
ECHMET_BEGIN_MAKE_LOGGER(LEMNGTracing, CALC_MATRIX_D1_DIMS, const size_t rows, const size_t cols)
{
	std::ostringstream ss{};

	ss << "Calculating Matrix D1(" << rows << ", " << cols << ")";
	return ss.str();
}
ECHMET_END_MAKE_LOGGER

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_MATRIX_D2_DIMS, "Matrix D2 dimensions")
ECHMET_BEGIN_MAKE_LOGGER(LEMNGTracing, CALC_MATRIX_D2_DIMS, const size_t rows, const size_t cols)
{
	std::ostringstream ss{};

	ss << "Calculating Matrix D2(" << rows << ", " << cols << ")";
	return ss.str();
}
ECHMET_END_MAKE_LOGGER

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_MATRIX_D1_ROW_BLOCK, "Matrix D1 row block")
ECHMET_BEGIN_MAKE_LOGGER(LEMNGTracing, CALC_MATRIX_D1_ROW_BLOCK, const std::string &name, const int mul, const size_t col, const double diffCoeff)
{
	std::ostringstream ss{};

	ss << "D1 row block[" << name << "], multiplicity = " << mul << ", col = " << col << ", diffCoeff = " << diffCoeff;

	return ss.str();
}
ECHMET_END_MAKE_LOGGER

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_MATRIX_D1_OUTPUT, "Matrix D1 output")
ECHMET_BEGIN_MAKE_LOGGER(LEMNGTracing, CALC_MATRIX_D1_OUTPUT, const ECHMET::LEMNG::Calculator::EMMatrix &DOne, const ECHMET::LEMNG::Calculator::CalculatorConstituentVec &ccVec, const ECHMET::LEMNG::Calculator::CalculatorIonicFormVec &cIfVec)
{
	std::ostringstream ss{};

	ss << "-- Matrix D1 --\n";
	ss << "Columns -> ";
	for (auto *iF : cIfVec)
		ss << iF->name << "; ";
	ss << "\nRows -> ";
	for (auto &cc : ccVec)
		ss << cc.name << "; ";
	ss << "\n";

	ss << "---\n\n" << DOne << "\n\n---";

	return ss.str();
}
ECHMET_END_MAKE_LOGGER

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_MATRIX_D2_OUTPUT, "Matrix D2 output")
ECHMET_BEGIN_MAKE_LOGGER(LEMNGTracing, CALC_MATRIX_D2_OUTPUT, const ECHMET::LEMNG::Calculator::EMMatrix &DOne, const ECHMET::LEMNG::Calculator::CalculatorConstituentVec &ccVec, const ECHMET::LEMNG::Calculator::CalculatorIonicFormVec &cIfVec)
{
	std::ostringstream ss{};

	ss << "-- Matrix D2 --\n";
	ss << "Columns -> ";
	for (auto *iF : cIfVec)
		ss << iF->name << "; ";
	ss << "\nRows -> ";
	for (auto &cc : ccVec)
		ss << cc.name << "; ";
	ss << "\n";

	ss << "---\n\n" << DOne << "\n\n---";

	return ss.str();
}
ECHMET_END_MAKE_LOGGER

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_MATRIX_D1_UICIF_BLOCK, "Matrix D1 uIcIF intermediate block beginning")
ECHMET_BEGIN_MAKE_LOGGER(LEMNGTracing, CALC_MATRIX_D1_UICIF_BLOCK, const char *block)
{
	std::ostringstream ss{};

	ss << "Entering Matrix D1 uIcIF block " << block;
	return ss.str();
}
ECHMET_END_MAKE_LOGGER

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_MATRIX_D1_UICIF_INTERMEDIATE, "Matrix D1 uIcIF intermediate block output")
ECHMET_BEGIN_MAKE_LOGGER(LEMNGTracing, CALC_MATRIX_D1_UICIF_INTERMEDIATE, const double s, const ECHMET::LEMNG::Calculator::CalculatorIonicForm *ccIf, const int d)
{
	std::ostringstream ss{};
	const std::string &name = ccIf->name;

	ss << "Entering Matrix D1 uIcIF intermediate(" << name << "):\n"
	   << "s(current) = " << s << "\n"
	   << "concentration = " << ccIf->concentration << "\n"
	   << "totalCharge = " << ccIf->charge << "\n"
	   << "mobility = " << ccIf->mobility << "\n"
	   << "ligand " << [](const auto iF) {
		   if (iF->ligand == nullptr)
			   return "(NO LIGAND)";
		   else
			   return iF->ligand->name->c_str();
	   }(ccIf->internalIonicForm)
	   << " multiplicity = " << d << "\n";

	return ss.str();
}
ECHMET_END_MAKE_LOGGER

#endif // ECHMET_TRACER_DISABLE_TRACING

} // namespace ECHMET
