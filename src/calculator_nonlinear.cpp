#include "calculator_nonlinear.h"
#include "calculator_common.h"
#include "calculator_matrices.h"
#include "calculator_linear.h"
#include "helpers.h"
#include <future>

#ifndef ECHMET_IMPORT_INTERNAL
#define ECHMET_IMPORT_INTERNAL
#endif // ECHMET_IMPORT_INTERNAL
#include <echmetcaes_extended.h>
#include <echmetphchconsts.h>

#include "tracing/lemng_tracer_impl.h"
#include <sstream>

namespace ECHMET {
namespace LEMNG {
namespace Calculator {

EigenzoneDispersion::EigenzoneDispersion(const double a2t, const double uEMD) :
	a2t{a2t},
	uEMD{uEMD}
{
}

static
EigenzoneDispersionVec calculateEigenzoneDispersion(const QLQRPack &QLQR, const EMMatrixVec &MDerivatives, const EMMatrix &concentrationDeltas, const EMMatrix &diffMatrix, const size_t NCO)
{
	const EMMatrixC &QL = QLQR.QL();
	const EMMatrixC &QR = QLQR.QR();
	EigenzoneDispersionVec ezDisps{};

	ezDisps.reserve(NCO);

	/* Precalculate QL * dM/dcK * QR products */
	std::vector<EMMatrixC> LMRs{};
	LMRs.reserve(NCO);
	for (size_t idx = 0; idx < NCO; idx++) {
		const EMMatrix MD = MDerivatives.at(idx);
		const EMMatrixC LMR = QL * MD * QR;

		LMRs.emplace_back(LMR);
	}

	/* Transformation to w domain.
	 * Hruška V, Riesová M, Gaš B, ELECTROPHORESIS 2012, Volume: 33, Pages: 923-930 (DOI: 10.1002/elps.201100554)
	 * states equation 18 in reverse order c = QR * w, we use QL to get w from concentration deltas.
	 */
	const EMMatrixC wVec = QL * concentrationDeltas;

	/* Diffusive parameters */
	const EMMatrixC LDiffR = QL * diffMatrix * QR;

	ECHMET_TRACE(LEMNGTracing, CALC_NONLIN_DIFF_PARAMS_MATRIX, std::cref(LDiffR));

	for (size_t idx = 0; idx < NCO; idx++) {
		/* Get diffusive parameter of the eigenzone */
		const double a2t = [&NCO, &LDiffR](int idx) {
			const double v = LDiffR(idx, idx).real();

			if (v <= 0)
				return 0.5;
			return v;
		}(idx);

		/* Calculate uEMD parameter of the eigenzone
		 * Step 1: Calculate dLambda/dW
		 * Uses equation 22 from Hruška V, Riesová M, Gaš B, ELECTROPHORESIS 2012, Volume: 33, Pages: 923-930 (DOI: 10.1002/elps.201100554)
		 */
		double dLdW = 0.0;

		/* Hruška V, Riesová M, Gaš B, ELECTROPHORESIS 2012, Volume: 33, Pages: 923-930 (DOI: 10.1002/elps.201100554)
		 * Equation 22 states that values from QR matrix shall be taken from positions i,k whereas they shall
		 * be taken from positions k,i
		 */
		for (size_t k = 0; k < NCO; k++)
			dLdW += QR(k, idx).real() * LMRs.at(k)(idx, idx).real();

		/* Step 2: Calculate uEMD */
		const double uEMD = dLdW * wVec(idx).real();

		ezDisps.emplace_back(a2t, uEMD);
	}

	return ezDisps;
}

static
EMMatrixVec calculateMDerivatives(const EMMatrix &MOne, const EMMatrix &MTwo, const EMMatrixVec &MOneDerivatives, const EMMatrixVec &MTwoDerivatives)
{
	EMMatrixVec MDerivatives{};
	MDerivatives.reserve(MOneDerivatives.size());

#ifdef ECHMET_LEMNG_PARALLEL_NUM_OPS
	const auto worker = [&](const size_t idx) -> EMMatrix {
		const EMMatrix MLeft = MOneDerivatives.at(idx) * MTwo;
		const EMMatrix MRight = MOne * MTwoDerivatives.at(idx);

		return MLeft + MRight;
	};

	const size_t N = MOneDerivatives.size();

	std::vector<std::future<EMMatrix>> results{};
	results.reserve(N);

	for (size_t idx = 0; idx < N; idx++)
		results.emplace_back(std::async(std::launch::async, worker, idx));

	for (auto &f: results)
		MDerivatives.emplace_back(f.get());
#else // ECHMET_LEMNG_PARALLEL_NUM_OPS
	for (size_t idx = 0; idx < MOneDerivatives.size(); idx++) {
		const EMMatrix MLeft = MOneDerivatives.at(idx) * MTwo;
		const EMMatrix MRight = MOne * MTwoDerivatives.at(idx);

		MDerivatives.emplace_back(MLeft + MRight);
	}
#endif // ECHMET_LEMNG_PARALLEL_NUM_OPS

	return MDerivatives;
}

static
EMMatrixVec calculateM1Derivatives(const CalculatorSystemPack &systemPack, const DeltaPackVec &deltaPacks)
{
	EMMatrixVec M1Derivatives{};

	M1Derivatives.reserve(deltaPacks.size());

#ifdef ECHMET_LEMNG_PARALLEL_NUM_OPS
	const auto worker = [&](const size_t idx) -> EMMatrix {
		return makeM1Derivative(systemPack, deltaPacks.at(idx));
	};

	const size_t N = deltaPacks.size();

	std::vector<std::future<EMMatrix>> results{};
	results.reserve(N);

	for (size_t idx = 0; idx < N; idx++)
		results.emplace_back(std::async(std::launch::async, worker, idx));

	for (auto &f : results)
		M1Derivatives.emplace_back(f.get());
#else // ECHMET_LEMNG_PARALLEL_NUM_OPS
	for (const auto &dPack : deltaPacks)
		M1Derivatives.emplace_back(makeM1Derivative(systemPack, dPack));
#endif // ECHMET_LEMNG_PARALLEL_NUM_OPS

	return M1Derivatives;
}

static
EMMatrixVec calculateM2Derivatives(const CalculatorSystemPack &systemPack, const RealVecPtr &analyticalConcentrations, const NonidealityCorrections corrections)
{
	const size_t NCO = systemPack.constituents.size();
	EMMatrixVec M2Derivatives{};

	M2Derivatives.reserve(NCO);

	const SysComp::ChemicalSystem &chemSystemRaw = *systemPack.chemSystemRaw;

	CAES::Solver *solver = nullptr;
	::ECHMET::RealVec *derivatives = nullptr;

	::ECHMET::RetCode tRet = CAES::prepareDerivatorContext(derivatives, solver, chemSystemRaw, corrections);
	if (tRet != ::ECHMET::RetCode::OK)
		throw CalculationException{std::string{"Cannot make derivator context: "} + std::string{errorToString(tRet)}, coreLibsErrorToNativeError(tRet)};

#ifdef ECHMET_LEMNG_PARALLEL_NUM_OPS
	const auto worker = [&](const size_t idx) -> EMMatrix {
		const SysComp::Constituent *pivotalConstituent = systemPack.constituents.at(idx).internalConstituent;
		RealVec *_derivatives = ::ECHMET::createRealVec(derivatives->size());
		if (_derivatives == nullptr)
			throw CalculationException{"Cannot allocate thread-local derivatives vector", RetCode::E_NO_MEMORY};
		if (_derivatives->resize(derivatives->size()) != ::ECHMET::RetCode::OK) {
			_derivatives->destroy();
			throw CalculationException{"Cannot resize thread-local derivatives vector", RetCode::E_NO_MEMORY};
		}

		const EMMatrix res = makeM2Derivative(systemPack, analyticalConcentrations, pivotalConstituent, solver, _derivatives);
		_derivatives->destroy();
		return res;
	};

	std::vector<std::future<EMMatrix>> results{};
	try {
		results.reserve(NCO);
	} catch (std::bad_alloc &) {
		solver->context()->destroy();
		solver->destroy();
		derivatives->destroy();
		throw;
	}

	try {
		for (size_t idx = 0; idx < NCO; idx++)
			results.emplace_back(std::async(std::launch::async, worker, idx));

		for (auto &f : results)
			M2Derivatives.emplace_back(f.get());
	} catch (CalculationException &) {
			for (auto &f : results) {
				if (f.valid())
					f.wait();
			}

			solver->context()->destroy();
			solver->destroy();
			derivatives->destroy();
			throw;
	}
#else // ECHMET_LEMNG_PARALLEL_NUM_OPS
	for (size_t idx = 0; idx < NCO; idx++) {
		const SysComp::Constituent *pivotalConstituent = systemPack.constituents.at(idx).internalConstituent;

		try {
			M2Derivatives.emplace_back(makeM2Derivative(systemPack, analyticalConcentrations, pivotalConstituent, solver, derivatives));
		} catch (CalculationException &ex) {
			solver->context()->destroy();
			solver->destroy();
			derivatives->destroy();
			throw CalculationException{"Cannot calculate concentration derivatives for M2 derivative", coreLibsErrorToNativeError(tRet)};
		}
	}
#endif // ECHMET_LEMNG_PARALLEL_NUM_OPS

	solver->context()->destroy();
	solver->destroy();
	derivatives->destroy();

	return M2Derivatives;
}

static
EMMatrix makeDiffusionMatrix(const CalculatorSystemPack &systemPackUncharged, const DeltaPackVec &deltaPackUncharged)
{
	static const auto calcDiffCoeff = [](const double mobility, const int32_t charge) {
		const int32_t _charge = charge == 0 ? 1 : charge;

		return mobility * ECHMET::PhChConsts::Tlab * ECHMET::PhChConsts::bk / (std::abs(_charge) * ECHMET::PhChConsts::e);
	};

	const CalculatorIonicFormVec &ionicForms = systemPackUncharged.ionicForms;
	const size_t NIF = ionicForms.size();
	const SysComp::CalculatedProperties *calcPropsRaw = systemPackUncharged.calcPropsRaw;

	ERVector diffusionCoefficients{};

	diffusionCoefficients.reserve(NIF);

	for (size_t ifIdx = 0; ifIdx < NIF; ifIdx++) {
		const SysComp::IonicForm *iF = ionicForms.at(ifIdx)->internalIonicForm;
		const ECHMETReal mobility = [calcPropsRaw](const SysComp::IonicForm *iF) -> double {
			if (iF->totalCharge != 0)
				return calcPropsRaw->ionicMobilities->at(iF->ionicMobilityIndex);

			const SysComp::Constituent *c = iF->nucleus;

			/* Our ionic form has zero charge. This does not work with our diffusion coefficient formula
			 * as it only works for charged particles. Approximate the diffusion coefficient from the
			 * charged ionic forms that are the most similar to our given form.
			 */
			const SysComp::IonicForm *chMinusOne = nullptr;
			const SysComp::IonicForm *chPlusOne = nullptr;
			auto findNextToZeroIFs = [](const SysComp::IonicForm *&chMinusOne, const SysComp::IonicForm *&chPlusOne, const SysComp::IonicForm *iF,  const SysComp::Constituent *c) {
				std::function<bool (const SysComp::IonicForm *, const SysComp::IonicForm *)> haveSameBuildingBlocks = [&haveSameBuildingBlocks](const SysComp::IonicForm *target, const SysComp::IonicForm *candidate) {
					/* Check whether both forms either are or are not a complex */
					const bool targetHasLigand = target->ligand != nullptr;
					const bool candidateHasLigand = candidate->ligand != nullptr;
					if (targetHasLigand != candidateHasLigand)
						return false;

					/* If we have complexes check if we have a same ligand */
					if (targetHasLigand && candidateHasLigand) {
						if (*(target->ligand->name) != *(candidate->ligand->name))
							return false;

						/* Check that both forms have the same ancestor structure. This is a fun one */
						const bool targetHasAncestor = target->ancestor != nullptr;
						const bool candidateHasAncestor = candidate->ancestor != nullptr;

						if (targetHasAncestor != candidateHasAncestor)
							return false;

						/* We have reached the bottom of the ancestor tree */
						if ((targetHasAncestor == false) && (candidateHasAncestor == false))
							return true;

						if (*(target->ancestor->nucleus->name) != *(candidate->ancestor->nucleus->name))
							return false;

						/* Dive down the ancestor tree with some recursive madness */
						return haveSameBuildingBlocks(target->ancestor, candidate->ancestor);
					}

					return *(target->nucleus->name) == *(candidate->nucleus->name);
				};

				for (size_t idx = 0; idx < c->ionicForms->size(); idx++) {
					const SysComp::IonicForm *ionicForm = c->ionicForms->at(idx);

					ECHMET_TRACE(LEMNGTracing, CALC_NONLIN_NEIGHBOUR_FORMS_LOOKUP, ionicForm->name->c_str(), iF->name->c_str());

					/* TODO: return early once we find all the ionic forms we need */

					if (ionicForm->totalCharge == -1 && chMinusOne == nullptr) {
						if (haveSameBuildingBlocks(iF, ionicForm))
							chMinusOne = ionicForm;
					} else if (ionicForm->totalCharge == 1 && chPlusOne == nullptr) {
						if (haveSameBuildingBlocks(iF, ionicForm))
							chPlusOne = ionicForm;
					}
				}
			};

			findNextToZeroIFs(chMinusOne, chPlusOne, iF, c);

			/* Acid */
			if (chMinusOne && !chPlusOne)
				return calcPropsRaw->ionicMobilities->at(chMinusOne->ionicMobilityIndex);
			/* Base */
			else if (!chMinusOne && chPlusOne)
				return calcPropsRaw->ionicMobilities->at(chPlusOne->ionicMobilityIndex);
			/* Ampholyte */
			else if (chMinusOne && chPlusOne) {
				const double mobMinusOne = calcPropsRaw->ionicMobilities->at(chMinusOne->ionicMobilityIndex);
				const double mobPlusOne = calcPropsRaw->ionicMobilities->at(chPlusOne->ionicMobilityIndex);

				return (mobMinusOne + mobPlusOne) / 2.0;
			} else
				return 20.0; /* Arbitrarily chosen mobility for constituents that really have no charge */
		}(iF);

		diffusionCoefficients.emplace_back(calcDiffCoeff(mobility, iF->totalCharge));
	}

	ECHMET_TRACE(LEMNGTracing, CALC_NONLIN_DIFFUSION_COEFFS, std::cref(systemPackUncharged), std::cref(diffusionCoefficients));

	const EMMatrix DOne = makeMatrixD1(systemPackUncharged, diffusionCoefficients);
	const EMMatrix DTwo = makeMatrixD2(systemPackUncharged, deltaPackUncharged);

	const EMMatrix diffMatrix = DOne * DTwo;

	ECHMET_TRACE(LEMNGTracing, CALC_NONLIN_DIFF_MATRIX, std::cref(diffMatrix));

	return diffMatrix;
}

static
EMMatrix makeConcentrationDeltas(const CalculatorSystemPack &systemPack)
{
	const size_t NCO = systemPack.constituents.size();
	EMMatrix deltaCVec{NCO, 1};

	for (size_t idx = 0; idx < NCO; idx++) {
		const CalculatorConstituent &cc = systemPack.constituents.at(idx);

		deltaCVec(idx) = cc.concentrationSample - cc.concentrationBGE;
	}

	return deltaCVec;
}

EigenzoneDispersionVec calculateNonlinear(const CalculatorSystemPack &systemPack, const CalculatorSystemPack &systemPackUncharged,
					  const RealVecPtr &analyticalConcentrations,
					  const DeltaPackVec &deltaPacks, const DeltaPackVec &deltaPacksUncharged,
					  const EMMatrix &M1, const EMMatrix &M2, const QLQRPack &QLQR,
					  const NonidealityCorrections corrections)
{
	ECHMET_TRACE(LEMNGTracing, CALC_NONLIN_PROGRESS, "Starting");

	const EMMatrixVec M1Derivatives = calculateM1Derivatives(systemPack, deltaPacks);
	const EMMatrixVec M2Derivatives = calculateM2Derivatives(systemPack, analyticalConcentrations, corrections);

	ECHMET_TRACE(LEMNGTracing, CALC_NONLIN_PROGRESS, "Individual matrix derivatives solved");

	const EMMatrix diffMatrix = makeDiffusionMatrix(systemPackUncharged, deltaPacksUncharged);
	const EMMatrixVec MDerivatives = calculateMDerivatives(M1, M2, M1Derivatives, M2Derivatives);
	const EMMatrix deltaCVec = makeConcentrationDeltas(systemPack);

	return calculateEigenzoneDispersion(QLQR, MDerivatives, deltaCVec, diffMatrix, systemPack.constituents.size());
}

} // namespace Calculator
} // namespace LEMNG

#ifndef ECHMET_TRACER_DISABLE_TRACING

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_NONLIN_PROGRESS, "Nonlinear calculations progress reports")
ECHMET_BEGIN_MAKE_LOGGER(LEMNGTracing, CALC_NONLIN_PROGRESS, const char *stage)
{
	return std::string{"Nonlinear calculations stage: "} + std::string{stage};
}
ECHMET_END_MAKE_LOGGER

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_NONLIN_NEIGHBOUR_FORMS_LOOKUP, "Neighbour ionic forms lookup for Nernst-Einstein")
ECHMET_BEGIN_MAKE_LOGGER(LEMNGTracing, CALC_NONLIN_NEIGHBOUR_FORMS_LOOKUP, const char *iF, const char *target)
{
	std::ostringstream ss{};

	ss << "Testing IF " << iF << " for target " << target;

	return ss.str();
}
ECHMET_END_MAKE_LOGGER

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_NONLIN_NERNST_EINST_INPUT, "Mobility and total charge of ionic form used in Nermst-Einstein equation")
ECHMET_BEGIN_MAKE_LOGGER(LEMNGTracing, CALC_NONLIN_NERNST_EINST_INPUT, const double mobility, const int32_t totalCharge)
{
	std::ostringstream ss{};

	ss << "u " << mobility << ", TC " << totalCharge;

	return ss.str();
}
ECHMET_END_MAKE_LOGGER

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_NONLIN_DIFFUSION_COEFFS, "Diffusion coefficients")
ECHMET_BEGIN_MAKE_LOGGER(LEMNGTracing, CALC_NONLIN_DIFFUSION_COEFFS, const ECHMET::LEMNG::Calculator::CalculatorSystemPack &systemPackUncharged, const ECHMET::LEMNG::Calculator::ERVector &diffCoeffs)
{
	std::ostringstream ss{};

	ss << "-- Diffusion coefficients --\n";
	for (size_t ifIdx = 0; ifIdx < systemPackUncharged.ionicForms.size(); ifIdx++)
		ss << systemPackUncharged.ionicForms.at(ifIdx)->name << "; " << diffCoeffs[ifIdx] << "\n";

	return ss.str();
}
ECHMET_END_MAKE_LOGGER

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_NONLIN_DIFF_MATRIX, "Diffusion matrix")
ECHMET_BEGIN_MAKE_LOGGER(LEMNGTracing, CALC_NONLIN_DIFF_MATRIX, const ECHMET::LEMNG::Calculator::EMMatrix &diffMatrix)
{
	std::ostringstream ss{};

	ss << "-- Diffusion matrix --\n"
	   << "---\n\n" << diffMatrix << "\n\n---";

	return ss.str();
}
ECHMET_END_MAKE_LOGGER

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_NONLIN_DIFF_PARAMS_MATRIX, "Diffusive parameters matrix")
ECHMET_BEGIN_MAKE_LOGGER(LEMNGTracing, CALC_NONLIN_DIFF_PARAMS_MATRIX, const ECHMET::LEMNG::Calculator::EMMatrixC &dpMatrix)
{
	std::ostringstream ss{};

	ss << "-- Diffusive parameters matrix --\n"
	   << "---\n\n" << dpMatrix << "\n\n---";

	return ss.str();
}
ECHMET_END_MAKE_LOGGER

#endif // ECHMET_TRACER_DISABLE_TRACING

} // namespace ECHMET
