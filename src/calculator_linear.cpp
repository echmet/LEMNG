#include "calculator_linear.h"
#include "calculator_common.h"
#include "calculator_matrices.h"
#include "helpers.h"
#include <cmath>

#ifndef ECHMET_IMPORT_INTERNAL
#define ECHMET_IMPORT_INTERNAL
#endif // ECHMET_IMPORT_INTERNAL
#include <echmetionprops.h>

#include "tracing/lemng_tracer_impl.h"
#include <sstream>

namespace ECHMET {
namespace LEMNG {
namespace Calculator {

typedef Eigen::ComplexEigenSolver<EMMatrixC> EMSolverC;

Eigenzone::Eigenzone(std::vector<double> &&zeroCC) noexcept :
	constituentConcentrations(zeroCC),
	zoneMobility{0},
	tainted{true},
	isAnalyzeZone{false}
{
}

Eigenzone::Eigenzone(const double zoneMobility, std::vector<double> &&constituentConcentrations, SolutionProperties &&solutionProperties, const bool tainted, const bool isAnalyzeZone) noexcept :
	constituentConcentrations(constituentConcentrations),
	solutionProperties{solutionProperties},
	zoneMobility{zoneMobility},
	tainted{tainted},
	isAnalyzeZone{isAnalyzeZone}
{
}

Eigenzone::Eigenzone(const Eigenzone &other) :
	constituentConcentrations(other.constituentConcentrations),
	solutionProperties{other.solutionProperties},
	zoneMobility{other.zoneMobility},
	tainted{other.tainted},
	isAnalyzeZone{other.isAnalyzeZone}
{
}

Eigenzone::Eigenzone(Eigenzone &&other) noexcept :
	constituentConcentrations(std::move(other.constituentConcentrations)),
	solutionProperties{std::move(other.solutionProperties)},
	zoneMobility{other.zoneMobility},
	tainted{other.tainted},
	isAnalyzeZone{other.isAnalyzeZone}
{
}

LinearResults::LinearResults(std::vector<Eigenzone> &&eigenzones, const QLQRPack &QLQR, const EMMatrix &M1, const EMMatrix &M2, const bool allZonesValid) :
	eigenzones(eigenzones),
	QLQR{QLQR},
	M1{M1},
	M2{M2},
	allZonesValid{allZonesValid}
{
}

LinearResults::LinearResults(std::vector<Eigenzone> &&eigenzones, QLQRPack &&QLQR, EMMatrix &&M1, EMMatrix &&M2, const bool allZonesValid) noexcept :
	eigenzones(eigenzones),
	QLQR{QLQR},
	M1{M1},
	M2{M2},
	allZonesValid{allZonesValid}
{
}

LinearResults::LinearResults(const LinearResults &other) :
	eigenzones(other.eigenzones),
	QLQR{other.QLQR},
	M1{other.M1},
	M2{other.M2},
	allZonesValid{other.allZonesValid}
{
}

LinearResults::LinearResults(LinearResults &&other) noexcept :
	eigenzones(std::move(other.eigenzones)),
	QLQR{std::move(other.QLQR)},
	M1{std::move(other.M1)},
	M2{std::move(other.M2)},
	allZonesValid{other.allZonesValid}
{
}

std::vector<EMMatrixC> calculatePMatrices(const EMMatrixC &QL, const EMMatrixC &QR)
{
	std::vector<EMMatrixC> PMatrices{};

	PMatrices.reserve(QR.rows());

	for (int idx = 0; idx < QR.rows(); idx++) {
		const EMMatrixC qRi = QR.col(idx);
		const EMMatrixC qLi = QL.row(idx);

		PMatrices.emplace_back(qRi * qLi);
	}

	return PMatrices;
}

QLQRPack calculateQLQR(EMSolverC &ces)
{
	EMMatrixC QR = ces.eigenvectors();
	EMMatrixC QL = QR.inverse();

	return QLQRPack{QL, QR};
}

std::vector<std::tuple<std::vector<double>, bool, bool>> calculateEigenzoneCompositions(const std::vector<EMMatrixC> &PMatrices, const CalculatorSystemPack &systemPack)
{
	auto isAnalytePresent = [](const double cInZone, const double cInSample) {
		return cInZone >= cInSample * 0.9;
	};

	const size_t NCO = systemPack.constituents.size();
	bool tainted = false;

	EMMatrixC deltaCVec{NCO, 1};	/* Note that this is a row vector */
	std::vector<std::tuple<std::vector<double>, bool, bool>> ezPackVec;

	for (size_t idx = 0; idx < NCO; idx++) {
		const CalculatorConstituent &cc = systemPack.constituents.at(idx);

		deltaCVec(idx) = cc.concentrationSample - cc.concentrationBGE;
	}

	int zoneCtr = 0;
	for (const auto &p : PMatrices) {
		const EMMatrixC ezConcDeltas = p * deltaCVec;
		bool isAnalyzeZone = false;

		std::vector<double> ezConcs{};
		ezConcs.reserve(NCO);
		for (size_t idx = 0; idx < NCO; idx++) {

			auto effectiveZero = [&tainted, zoneCtr](const double c, const std::string &name) {
				if (std::abs(c) <= 1.0e-13)
					return 1.0e-13;
				else if (c < 0.0) {
					_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::CALC_LIN_ZONE_TAINTED, const int &, const std::string &, const double &>(zoneCtr, name, c);

					tainted = true;
					return 1.0e-13;
				}

				return c;
			};

			const CalculatorConstituent &cc = systemPack.constituents.at(idx);
			const double c = cc.concentrationBGE + ezConcDeltas(idx).real(); /* Convert the delta to actual concentration */

			if (cc.isAnalyte && isAnalytePresent(c, cc.concentrationSample))
				isAnalyzeZone = true;

			ezConcs.emplace_back(effectiveZero(c, cc.name));
		}
		zoneCtr++;

		ezPackVec.emplace_back(std::make_tuple(std::move(ezConcs), tainted, isAnalyzeZone));
	}

	return ezPackVec;
}

LinearResults calculateLinear(const CalculatorSystemPack &systemPack, const DeltaPackVec &deltaPacks, const NonidealityCorrections corrections)
{
	/* Calculate the mobility matrix. */
	EMMatrix M1{};
	EMMatrix M2{};
	EMMatrix MFin{};

	_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::CALC_LIN_PROGRESS, const char *>("Starting");

	try {
		M1 = makeMatrixM1(systemPack);
		M2 = makeMatrixM2(systemPack, deltaPacks);
		MFin = M1 * M2;

		_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::CALC_LIN_MFIN, const EMMatrix &>(MFin);
	} catch (std::bad_alloc &) {
		throw CalculationException{"Insufficient memory to calculate mobility matrix", RetCode::E_NO_MEMORY};
	}

	_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::CALC_LIN_PROGRESS, const char *>("Solving eigenzones' compositions");

	if (MFin.rows() < 1)
		return LinearResults{{}, QLQRPack{EMMatrixC{0,0}, EMMatrixC{0,0}}, std::move(M1), std::move(M2), true};


	/* Calculate eigenmobilites and zone compositions.
	 * Eigenmobilities are the eigenvalues of the MFin matrix.
	 * Zone composition are derived from the QL and QR eigenvectors.
	 */
	try {
		EMSolverC ces{MFin};
		EMVectorC eigenmobs = ces.eigenvalues();
		_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::CALC_EIGENMOBS, const EMVectorC &>(eigenmobs);
		if (isComplex(eigenmobs))
			throw CalculationException{"Detected complex eigenmobilities", RetCode::E_COMPLEX_EIGENMOBILITIES};

		const QLQRPack QLQR = calculateQLQR(ces);
		const std::vector<EMMatrixC> PMatrices = calculatePMatrices(QLQR.QL(), QLQR.QR());
		auto eigenzoneCompositions = calculateEigenzoneCompositions(PMatrices, systemPack);

		std::vector<Eigenzone> eigenzones{};
		eigenzones.reserve(eigenzoneCompositions.size());
		bool allZonesValid = true;
		for (size_t idx = 0; idx < eigenzoneCompositions.size(); idx++) {
			auto &&ez = std::get<0>(eigenzoneCompositions.at(idx));
			const bool tainted = std::get<1>(eigenzoneCompositions.at(idx));
			const bool isAnalyzeZone = std::get<2>(eigenzoneCompositions.at(idx));
			RealVecPtr zoneConcsVec = makeAnalyticalConcentrationsVec(systemPack.chemSystemRaw);
			CalculatedPropertiesPtr zoneCalcProps = makeCalculatedProperties(systemPack.chemSystemRaw);
			const double zoneMobility = eigenmobs(idx).real();

			/* Analytical concentrations in eigenzones are ordered by the CalculatorSystemPack
			 * ordering which may not correspond to the SysComp ordering.
			 * We need to make sure that we remap them back to SysComp ordering.
			 */
			for (size_t jdx = 0; jdx < systemPack.constituents.size(); jdx++) {
				const CalculatorConstituent &cc = systemPack.constituents.at(jdx);
				const size_t scIdx = cc.internalConstituent->analyticalConcentrationIndex;

				(*zoneConcsVec)[scIdx] = ez.at(jdx);
			}

			try {
				SolutionProperties zoneProps = calculateSolutionProperties(systemPack.chemSystemRaw, zoneConcsVec, zoneCalcProps.get(), corrections);
				eigenzones.emplace_back(zoneMobility, std::move(ez), std::move(zoneProps), tainted, isAnalyzeZone);
			} catch (CalculationException &) {
				std::vector<double> dummyCC;
				dummyCC.resize(eigenzoneCompositions.size());

				eigenzones.emplace_back(Eigenzone{std::move(dummyCC)});
				allZonesValid = false;
			}
		}

		_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::CALC_LIN_PROGRESS, const char *>("Done");

		return LinearResults{std::move(eigenzones), std::move(QLQR), std::move(M1), std::move(M2), allZonesValid};
	} catch (std::bad_alloc &) {
		throw CalculationException{"Insufficient memory to calculate eigenzone compositions", RetCode::E_NO_MEMORY};
	} catch (SysCompException &ex) {
		throw CalculationException{"SysComp library exception: " + std::string{ex.what()}, coreLibsErrorToNativeError(ex.errorCode())};
	}
}

} // namespace Calculator
} // namespace LEMNG

#ifndef ECHMET_TRACER_DISABLE_TRACING

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_LIN_PROGRESS, "Linear calculations progress reports")
ECHMET_MAKE_LOGGER(LEMNGTracing, CALC_LIN_PROGRESS, const char *stage)
{
	return std::string{"Linear calculations stage: "} + std::string{stage};
}

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_LIN_MFIN, "Linear mobility matrix")
ECHMET_MAKE_LOGGER(LEMNGTracing, CALC_LIN_MFIN, const ECHMET::LEMNG::Calculator::EMMatrix &MFin)
{
	std::ostringstream ss{};

	ss << "-- Matrix MFin --\n";
	ss << "---\n\n" << MFin << "\n\n---";

	return ss.str();
}

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_LIN_ZONE_TAINTED, "Eigenzone tainted")
ECHMET_MAKE_LOGGER(LEMNGTracing, CALC_LIN_ZONE_TAINTED, const int &zoneNum, const std::string &offendingConstituent, const double &calculatedConc)
{
	std::ostringstream ss{};

	ss << "Zone " << zoneNum << " is tainted, concentration of " << offendingConstituent << " was computed as " << calculatedConc << " (mmol/dm3), clapming to zero";

	return ss.str();
}

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_EIGENMOBS, "Complex eigenmobilities")
ECHMET_MAKE_LOGGER(LEMNGTracing, CALC_EIGENMOBS, const LEMNG::Calculator::EMVectorC &mobilities)
{
	std::ostringstream ss{};

	for (int idx = 0; idx < mobilities.cols(); idx++) {
		const auto &cu = mobilities(idx);
		ss << "Real: " << cu.real() << "; Imag: " << cu.imag() << "\n";
	}

	return ss.str();
}

#endif // ECHMET_TRACER_DISABLE_TRACING

} // namespace ECHMET
