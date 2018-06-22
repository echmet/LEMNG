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
int M1KroeneckerDelta(const size_t i, const std::vector<size_t> &mList)
{
	for (const size_t m : mList)
		if (i == m) return 1;

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

	_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::CALC_MATRIX_DM1_INPUT, const double&, const EMVector&>(baseConductivity, deltaPack.concentrationDeltas);

	for (size_t row = 0; row < ROWS; row++) {
		const CalculatorConstituent &c = ccVec.at(row);

		const double uIcIdKdC = [&c, baseConductivity, dKdcJ]() {
			double s = 0.0;

			for (const CalculatorIonicForm *iF : c.ionicForms) {
				s += iF->concentration * iF->mobility * cxsgn(iF->charge);
				_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::CALC_MATRIX_DM1_UICIDKDC, const double&, const double&, const int32_t&>(s, iF->mobility, cxsgn(iF->charge));
			}

			return s * 2.0 / std::pow(baseConductivity, 2) * dKdcJ;
		}();

		const double uIdcIdcJ = [&c, &deltaPack, baseConductivity]() {
			double s = 0.0;

			for (const CalculatorIonicForm *iF : c.ionicForms) {
				s += iF->mobility * cxsgn(iF->charge) * deltaPack.concentrationDeltas(iF->globalIonicFormConcentrationIdx);
				_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::CALC_MATRIX_DM1_UIDCIDCJ, const std::string&, const size_t&, const double&, const double&>(iF->name, iF->globalIonicFormConcentrationIdx, deltaPack.concentrationDeltas(iF->globalIonicFormConcentrationIdx), s);
			}

			return -s / baseConductivity;
		}();

		const double termTwo = (uIcIdKdC + uIdcIdcJ) * PhChConsts::F;

		for (size_t col = 0; col < COLS - 2; col++) {
			const CalculatorIonicForm *iF = ifVec.at(col);
			const auto &ctuentList = iF->containedConstituents;
			const int32_t charge = iF->charge;
			const double mobility = iF->mobility;
			const int d = M1KroeneckerDelta(row, ctuentList);

			double result = -d * cxsgn(charge) * dKdcJ / baseConductivity;
			result += termTwo * std::abs(charge);
			result *= mobility;

			_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::CALC_MATRIX_DM1_KRD, const int&>(d);

			MOneDer(row, col) = result;
		}

		const double H3OMob = ifVec.at(H3O_idx)->mobility;
		const double OHMob = ifVec.at(OH_idx)->mobility;

		_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::CALC_MATRIX_DM1_ROW_BLOCK, const double&, const double&, const double&, const double&, const double&>(uIcIdKdC, uIdcIdcJ, termTwo, H3OMob, OHMob);

		/* H3O+ and OH- are in the last two columns */
		MOneDer(row, H3O_idx) = termTwo * H3OMob;
		MOneDer(row, OH_idx) = termTwo * OHMob;
	}

	_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::CALC_MATRIX_DM1_OUTPUT, const EMMatrix&>(MOneDer);

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

	_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::CALC_MATRIX_DM2_OUTPUT, const EMMatrix&>(MTwoDer);

	return MTwoDer;
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

	_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::CALC_MATRIX_M1_DIMS, const size_t&, const size_t&>(ROWS, COLS);

	for (size_t row = 0; row < ROWS; row++) {
		const CalculatorConstituent &c = ccVec.at(row);

		const double uIcIFSum = [](const CalculatorConstituent &c, const double conductivity) {
			double s = 0.0;

			_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::CALC_MATRIX_M1_UICIF_BLOCK, const char*>(c.name.c_str());

			for (size_t idx = 0; idx < c.ionicForms.size(); idx++) {
				const CalculatorIonicForm *iF = c.ionicForms.at(idx);

				s += iF->concentration * iF->mobility * cxsgn(iF->charge);
				_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::CALC_MATRIX_M1_UICIF_INTERMEDIATE, const double&, const CalculatorIonicForm*&>(s, iF);
			}

			return s * PhChConsts::F / (conductivity * 1.0e9);
		}(c, systemPack.conductivity);

		for (size_t col = 0; col < COLS - 2; col++) {
			const CalculatorIonicForm *iF = ifVec.at(col);
			const auto &ctuentList = iF->containedConstituents;
			const int32_t charge = iF->charge;
			const double mobility = iF->mobility;
			const int d = M1KroeneckerDelta(row, ctuentList);

			_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::CALC_MATRIX_M1_ROW_BLOCK, const std::string&, const int &, const size_t&, const double &>(iF->name, d, col, mobility);

			MOne(row, col) = (d * cxsgn(charge) - uIcIFSum * std::abs(charge)) * mobility;
		}

		/* H3O+ and OH- are in the last two columns */
		MOne(row, H3O_idx) = -uIcIFSum * ifVec.at(H3O_idx)->mobility;
		MOne(row, OH_idx) = -uIcIFSum * ifVec.at(OH_idx)->mobility;
	}

	_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::CALC_MATRIX_M1_OUTPUT, const EMMatrix&, const CalculatorConstituentVec&, const CalculatorIonicFormVec&>(MOne, systemPack.constituents, systemPack.ionicForms);

	return MOne;
}

EMMatrix makeMatrixM2(const CalculatorSystemPack &systemPack, const DeltaPackVec &deltaPacks)
{
	const size_t ROWS = systemPack.ionicForms.size();
	const size_t COLS = systemPack.constituents.size();
	EMMatrix MTwo{ROWS, COLS};

	_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::CALC_MATRIX_M2_DIMS, const size_t&, const size_t&>(ROWS, COLS);

	for (size_t col = 0; col < COLS; col++)
		MTwo.col(col) = deltaPacks.at(col).concentrationDeltas;


	_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::CALC_MATRIX_M2_OUTPUT, const EMMatrix&, const CalculatorConstituentVec&, const CalculatorIonicFormVec&>(MTwo, systemPack.constituents, systemPack.ionicForms);

	return MTwo;
}

} // namespace Calculator
} // namespace LEMNG

#ifndef ECHMET_TRACER_DISABLE_TRACING

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_MATRIX_M1_DIMS, "Matrix M1 dimensions")
ECHMET_MAKE_LOGGER(LEMNGTracing, CALC_MATRIX_M1_DIMS, const size_t &rows, const size_t &cols)
{
	std::ostringstream ss{};

	ss << "Calculating Matrix M1(" << rows << ", " << cols << ")";
	return ss.str();
}

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_MATRIX_M2_DIMS, "Matrix M2 dimensions")
ECHMET_MAKE_LOGGER(LEMNGTracing, CALC_MATRIX_M2_DIMS, const size_t &rows, const size_t &cols)
{
	std::ostringstream ss{};

	ss << "Calculating Matrix M2(" << rows << ", " << cols << ")";
	return ss.str();
}

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_MATRIX_M1_UICIF_BLOCK, "Matrix M1 uIcIF intermediate block beginning")
ECHMET_MAKE_LOGGER(LEMNGTracing, CALC_MATRIX_M1_UICIF_BLOCK, const char *block)
{
	std::ostringstream ss{};

	ss << "Entering Matrix M1 uIcIF block " << block;
	return ss.str();
}

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_MATRIX_M1_UICIF_INTERMEDIATE, "Matrix M1 uIcIF intermediate block output")
ECHMET_MAKE_LOGGER(LEMNGTracing, CALC_MATRIX_M1_UICIF_INTERMEDIATE, const double &s, const ECHMET::LEMNG::Calculator::CalculatorIonicForm *&ccIf)
{
	std::ostringstream ss{};
	const std::string &name = ccIf->name;

	ss << "Entering Matrix M1 uIcIF intermediate(" << name << "):\n"
	   << "s(current) = " << s << "\n"
	   << "concentration = " << ccIf->concentration << "\n"
	   << "totalCharge = " << ccIf->charge << "\n"
	   << "mobility = " << ccIf->mobility << "\n";

	return ss.str();
}

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_MATRIX_M1_ROW_BLOCK, "Matrix M1 row block")
ECHMET_MAKE_LOGGER(LEMNGTracing, CALC_MATRIX_M1_ROW_BLOCK, const std::string &name, const int &KrD, const size_t &col, const double &mobility)
{
	std::ostringstream ss{};

	ss << "M1 row block[" << name << "], KroeneckerDelta = " << KrD << ", col = " << col << ", mobility = " << mobility;

	return ss.str();
}

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_MATRIX_M1_OUTPUT, "Matrix M1 output")
ECHMET_MAKE_LOGGER(LEMNGTracing, CALC_MATRIX_M1_OUTPUT, const ECHMET::LEMNG::Calculator::EMMatrix &MOne, const ECHMET::LEMNG::Calculator::CalculatorConstituentVec &ccVec, const ECHMET::LEMNG::Calculator::CalculatorIonicFormVec &cIfVec)
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

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_MATRIX_M2_OUTPUT, "Matrix M2 output")
ECHMET_MAKE_LOGGER(LEMNGTracing, CALC_MATRIX_M2_OUTPUT, const ECHMET::LEMNG::Calculator::EMMatrix &MTwo, const ECHMET::LEMNG::Calculator::CalculatorConstituentVec &ccVec, const ECHMET::LEMNG::Calculator::CalculatorIonicFormVec &cIfVec)
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

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_MATRIX_DM1_INPUT, "dM1/dC input")
ECHMET_MAKE_LOGGER(LEMNGTracing, CALC_MATRIX_DM1_INPUT, const double &baseConductivity, const ECHMET::LEMNG::Calculator::EMVector &cDeltas)
{
	std::ostringstream ss{};

	ss << "dM1/dC imput\n"
	   << "base conductivity = " << baseConductivity << " (S/m)\n"
	   << "concentrations deltas (mmol/dm3)\n"
	   << cDeltas;

	return ss.str();
}

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_MATRIX_DM1_UICIDKDC, "dM1/dC dConductivity/dC")
ECHMET_MAKE_LOGGER(LEMNGTracing, CALC_MATRIX_DM1_UICIDKDC, const double &s, const double &mobility, const int32_t &absCharge)
{
	std::ostringstream ss{};

	ss << "dM1/dC dConductivity/dC\n"
	   << "s = " << s << ", mobility = " << mobility << ", charge = " << absCharge;

	return ss.str();
}

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_MATRIX_DM1_UIDCIDCJ, "dM1/dC dCi/dCj")
ECHMET_MAKE_LOGGER(LEMNGTracing, CALC_MATRIX_DM1_UIDCIDCJ, const std::string &name, const size_t &gIdx, const double &cDelta, const double &s)
{
	std::ostringstream ss{};

	ss << "dM1/dC dCi/dCj\n"
	   << name << ", Global IFIdx = " << gIdx << ", concDelta = " << cDelta << ", s = " << s;

	return ss.str();
}

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_MATRIX_DM1_KRD, "dM1/dC Kroenecker delta")
ECHMET_MAKE_LOGGER(LEMNGTracing, CALC_MATRIX_DM1_KRD, const int &KrD)
{
	return std::string{"Kroenecker delta = "} + std::to_string(KrD);
}

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_MATRIX_DM1_ROW_BLOCK, "dM1/dC row block")
ECHMET_MAKE_LOGGER(LEMNGTracing, CALC_MATRIX_DM1_ROW_BLOCK, const double &uIcIcKdC, const double &uIdcIdcJ, const double &termTwo, const double &H3OMob, const double &OHMob)
{
	std::ostringstream ss{};

	ss << "dM1/dC row block\n"
	   << "uIcIcKcD = " << uIcIcKdC << ", uIdcIdcJ = " << uIdcIdcJ << "\n"
	   << "termTwo = " << termTwo
	   << "water mobs = " << H3OMob << ", " << OHMob;

	return ss.str();
}

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_MATRIX_DM1_OUTPUT, "dM1/dC output")
ECHMET_MAKE_LOGGER(LEMNGTracing, CALC_MATRIX_DM1_OUTPUT, const ECHMET::LEMNG::Calculator::EMMatrix &MOneDer)
{
	std::ostringstream ss{};

	ss << "-- Matrix dM1/dC --\n"
	   << "---\n\n" << MOneDer << "\n\n---";

	return ss.str();
}

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_MATRIX_DM2_OUTPUT, "dM2/dC output")
ECHMET_MAKE_LOGGER(LEMNGTracing, CALC_MATRIX_DM2_OUTPUT, const ECHMET::LEMNG::Calculator::EMMatrix &MTwoDer)
{
	std::ostringstream ss{};

	ss << "-- Matrix dM2/dC --\n"
	   << "---\n\n" << MTwoDer << "\n\n---";

	return ss.str();
}

#endif // ECHMET_TRACER_DISABLE_TRACING

} // namespace ECHMET
