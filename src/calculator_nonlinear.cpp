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
	 * incorrectly states that the concentration deltas vector shall be multiplied by QR instead. */
	const EMMatrixC wVec = QL * concentrationDeltas;

	/* Diffusive parameters */
	const EMMatrixC LDiffR = QL * diffMatrix * QR;

	_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::CALC_NONLIN_DIFF_PARAMS_MATRIX, const EMMatrixC&>(LDiffR);

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
EMMatrix makeDiffusionMatrix(const CalculatorSystemPack &systemPack, const RealVecPtr &analyticalConcentrations, const ConcentrationDeltasVec &concentrationDeltasVec)
{
	/* WARNING - This piece of code is a terrible read. */
	typedef std::vector<ERVector> ERVecVec;

	const SysComp::CalculatedProperties &calcPropsRaw = *systemPack.calcPropsRaw;
	const size_t NCO = systemPack.constituents.size();
	const double kappa = 1.0e9 * systemPack.conductivity / PhChConsts::F;
	EMMatrix diffMatrix = EMMatrix::Constant(NCO, NCO, 0);

	static const auto calcDiffCoeff = [](const double mobility, const int32_t charge) {
		const int32_t _charge = charge == 0 ? 1 : charge;

		return mobility * ECHMET::PhChConsts::Tlab * ECHMET::PhChConsts::bk / (std::abs(_charge) * ECHMET::PhChConsts::e);
	};

	/* Precalculate degree of dissociation and its derivative.
	 * Outer vector is for each constituent, inner vector for each of its ionic forms. */
	ERVecVec diffCoeffs;
	ERVecVec dissocDegrees;
	ERVecVec dDissocDegreesDX;

	diffCoeffs.reserve(NCO);
	dissocDegrees.reserve(NCO);
	dDissocDegreesDX.reserve(NCO);

	/* A note to your future self - DO NOT try to figure this out. */
	for (size_t idx = 0; idx < NCO; idx++) {
		const SysComp::Constituent *c = systemPack.constituents.at(idx).internalConstituent;
		const size_t IFCount = c->ionicForms->size();
		const ERVector &concentrationDeltas = concentrationDeltasVec.at(idx);
		ERVector diffCoeffsForCtuent;
		ERVector dissocDegreesForCtuent;
		ERVector dDissocDegreesDXForCtuent;

		diffCoeffsForCtuent.reserve(IFCount);
		dissocDegreesForCtuent.reserve(IFCount);
		dDissocDegreesDXForCtuent.reserve(IFCount);

		for (size_t ifIdx = 0; ifIdx < IFCount; ifIdx++) {
			const SysComp::IonicForm *iF = c->ionicForms->at(ifIdx);
			const ECHMETReal mobility = [&calcPropsRaw](const SysComp::IonicForm *iF, const SysComp::Constituent *c) -> double {
				if (iF->totalCharge != 0)
					return calcPropsRaw.ionicMobilities->at(iF->ionicMobilityIndex);

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

						_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::CALC_NONLIN_NEIGHBOUR_FORMS_LOOKUP, const char *, const char *>(ionicForm->name->c_str(), iF->name->c_str());

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
					return calcPropsRaw.ionicMobilities->at(chMinusOne->ionicMobilityIndex);
				/* Base */
				else if (!chMinusOne && chPlusOne)
					return calcPropsRaw.ionicMobilities->at(chPlusOne->ionicMobilityIndex);
				/* Ampholyte */
				else if (chMinusOne && chPlusOne) {
					const double mobMinusOne = calcPropsRaw.ionicMobilities->at(chMinusOne->ionicMobilityIndex);
					const double mobPlusOne = calcPropsRaw.ionicMobilities->at(chPlusOne->ionicMobilityIndex);

					return (mobMinusOne + mobPlusOne) / 2.0;
				} else
					return 20.0; /* Arbitrarily chosen mobility for constituents that really have no charge */
			}(iF, c);

			_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::CALC_NONLIN_NERNST_EINST_INPUT, const double&, const int32_t&>( mobility, iF->totalCharge);

			const double anC = analyticalConcentrations->at(c->analyticalConcentrationIndex);
			const double ionicC = calcPropsRaw.ionicConcentrations->at(iF->ionicConcentrationIndex);
			const double dissocDegree = ionicC / anC;

			diffCoeffsForCtuent.emplace_back(calcDiffCoeff(mobility, iF->totalCharge));
			dissocDegreesForCtuent.emplace_back(dissocDegree);

			/* The diffusion model expects the "equilibrium communication" expressed
			 * as derivatives of degrees of dissociation whereas CAES provides
			 * derivatives of concentrations.
			 * Here we convert concentration derivatives to derivatives of
			 * degrees of dissociation. */
                        const double cDer = concentrationDeltas.at(iF->ionicConcentrationIndex);
			const double dCHdCX = concentrationDeltas.at(0);
			const double dDissocDegreeDX = (cDer - dissocDegree) / (dCHdCX * anC);

			dDissocDegreesDXForCtuent.emplace_back(dDissocDegreeDX);
		}

		diffCoeffs.emplace_back(diffCoeffsForCtuent);
		dissocDegrees.emplace_back(dissocDegreesForCtuent);
		dDissocDegreesDX.emplace_back(dDissocDegreesDXForCtuent);
	}

	/* Precalculate terms that are used commonly in the final calculation */
	const auto calcTermA = [&diffCoeffs, &dissocDegrees]() {
		assert(diffCoeffs.size() == dissocDegrees.size());

		ERVector termA;
		termA.reserve(diffCoeffs.size());

		for (size_t idx = 0; idx < diffCoeffs.size(); idx++) {
			const ERVector &diffCoeffsForCtuent = diffCoeffs.at(idx);
			const ERVector &dissocDegreesForCtuent = dissocDegrees.at(idx);
			double t = 0.0;

			assert(diffCoeffsForCtuent.size() == dissocDegreesForCtuent.size());

			for (size_t ifIdx = 0; ifIdx < diffCoeffsForCtuent.size(); ifIdx++) {
				const ECHMET::ECHMETReal diff = diffCoeffsForCtuent.at(ifIdx);
				const ECHMET::ECHMETReal dissoc = dissocDegreesForCtuent.at(ifIdx);

				t += diff * dissoc;
			}

			termA.emplace_back(t);
		}

		return termA;

	};
	const auto calcTermAz = [&systemPack, NCO, &diffCoeffs, &dissocDegrees]() {
		assert(diffCoeffs.size() == dissocDegrees.size());
		assert(diffCoeffs.size() == NCO);

		ERVector termAz;
		termAz.reserve(NCO);

		for (size_t idx = 0; idx < NCO; idx++) {
			const ERVector &diffCoeffsForCtuent = diffCoeffs.at(idx);
			const ERVector &dissocDegreesForCtuent = dissocDegrees.at(idx);
			const SysComp::Constituent *c = systemPack.constituents.at(idx).internalConstituent;
			double t = 0.0;

			assert(diffCoeffsForCtuent.size() == dissocDegreesForCtuent.size());
			assert(diffCoeffsForCtuent.size() == c->ionicForms->size());

			for (size_t ifIdx = 0; ifIdx < c->ionicForms->size(); ifIdx++) {
				const SysComp::IonicForm *iF = c->ionicForms->at(ifIdx);
				const double diff = diffCoeffsForCtuent.at(ifIdx);
				const double dissoc = dissocDegreesForCtuent.at(ifIdx);

				t += diff * dissoc * iF->totalCharge;
			}

			termAz.emplace_back(t);
		}

		return termAz;

	};
	const auto calcTermB = [&systemPack, &calcPropsRaw, NCO, &dissocDegrees]() {
		ERVector termB;
		termB.reserve(NCO);

		for (size_t idx = 0; idx < NCO; idx++) {
			const ERVector &dissocDegreesForCtuent = dissocDegrees.at(idx);
			const SysComp::Constituent *c = systemPack.constituents.at(idx).internalConstituent;
			assert(dissocDegreesForCtuent.size() == c->ionicForms->size());

			double t = 0.0;

			for (size_t ifIdx = 0; ifIdx < c->ionicForms->size(); ifIdx++) {
				const SysComp::IonicForm *iF = c->ionicForms->at(ifIdx);
				const double mobility = calcPropsRaw.ionicMobilities->at(iF->ionicMobilityIndex);

				t += cxsgn(iF->totalCharge) * mobility * dissocDegreesForCtuent.at(ifIdx);
			}

			termB.emplace_back(t);
		}

		return termB;
	};
	const auto calcTermC = [&systemPack, &calcPropsRaw, NCO, &diffCoeffs, &dDissocDegreesDX, &kappa, &analyticalConcentrations](const ERVector &termB) {
		const double diffCoeffH3O = calcDiffCoeff(PhChConsts::mobilityH3O, 1);
		const double diffCoeffOH = calcDiffCoeff(PhChConsts::mobilityOH, -1);
		double sum = 0.0;

		/* Calculate the sum */
		for (size_t idx = 0; idx < NCO; idx++) {
			const ERVector &diffCoeffsForCtuent = diffCoeffs.at(idx);
			const ERVector &dDissocDegreesDXForCtuent = dDissocDegreesDX.at(idx);
			const SysComp::Constituent *c = systemPack.constituents.at(idx).internalConstituent;

			assert(dDissocDegreesDXForCtuent.size() == diffCoeffsForCtuent.size());
			assert(dDissocDegreesDXForCtuent.size()	== c->ionicForms->size());

			double t = 0.0;

			for (size_t ifIdx = 0; ifIdx < c->ionicForms->size(); ifIdx++) {
				const SysComp::IonicForm *iF = c->ionicForms->at(ifIdx);

				t += iF->totalCharge * diffCoeffsForCtuent.at(ifIdx) * dDissocDegreesDXForCtuent.at(ifIdx);
			}

			const double anC = analyticalConcentrations->at(c->analyticalConcentrationIndex);

			sum += t * anC;
		}

		/* PeakMaster 5.3 ignores the water terms here
		sum += calcPropsRaw.ionicConcentrations->at(0) * diffCoeffH3O;
		sum += calcPropsRaw.ionicConcentrations->at(1) * diffCoeffOH; */

		_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::CALC_NONLIN_DIFFCOEFF_TERM_C_SUM, const double&>(sum);

		const double v = sum + diffCoeffH3O + diffCoeffOH * (PhChConsts::KW_298 * 1.0e6) / (std::pow(calcPropsRaw.ionicConcentrations->at(0), 2));

		_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::CALC_NONLIN_DIFFCOEFF_TERM_C_V, const double&>(v);

		/* Calculate internal term */
		ERVector termInternal;		/* Zi3 in PeakMaster 5.3 */
		termInternal.reserve(NCO);

		for (size_t idx = 0; idx < NCO; idx++) {
			const SysComp::Constituent *c = systemPack.constituents.at(idx).internalConstituent;
			const ERVector &diffCoeffsForCtuent = diffCoeffs.at(idx);
			const ERVector &dDissocDegreesDXForCtuent = dDissocDegreesDX.at(idx);
			double t = 0;

			for (size_t ifIdx = 0; ifIdx < c->ionicForms->size(); ifIdx++) {
				_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::CALC_NONLIN_DIFFCOEFF_TERM_C_INTER2, const double&, const double&>(diffCoeffsForCtuent.at(ifIdx), dDissocDegreesDXForCtuent.at(ifIdx));
				t += diffCoeffsForCtuent.at(ifIdx) * dDissocDegreesDXForCtuent.at(ifIdx);
			}

			termInternal.emplace_back(t);
		}

		/* Now finally calculate the termC */
		ERVector termC;
		termC.reserve(NCO);

		for (size_t idx = 0; idx < NCO; idx++) {
			const SysComp::Constituent *c = systemPack.constituents.at(idx).internalConstituent;
			const double anC = analyticalConcentrations->at(c->analyticalConcentrationIndex);

			const double t = anC * (termInternal.at(idx) - (termB.at(idx) / kappa) * v);
			termC.emplace_back(t);
		}

		return termC;
	};
	const auto calcTermD = [&systemPack, &calcPropsRaw, NCO, &dissocDegrees, &dDissocDegreesDX, &analyticalConcentrations]() {
		double sum = 0.0;

		/* Calculate the sum */
		for (size_t idx = 0; idx < NCO; idx++) {
			const SysComp::Constituent *c = systemPack.constituents.at(idx).internalConstituent;
			const ERVector &dDissocDegreesDXForCtuent = dDissocDegreesDX.at(idx);

			for (size_t ifIdx = 0; ifIdx < c->ionicForms->size(); ifIdx++) {
				const SysComp::IonicForm *iF = c->ionicForms->at(ifIdx);
				const double anC = analyticalConcentrations->at(c->analyticalConcentrationIndex);

				sum += anC * iF->totalCharge * dDissocDegreesDXForCtuent.at(ifIdx);
			}
		}

		_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::CALC_NONLIN_DIFFCOEFF_TERM_D_SUM, const double&>(sum);

		/* Calculate termD */
		ERVector termD;
		termD.reserve(NCO);

		const double denominator = 1 + ((PhChConsts::KW_298 * 1.0e6) / std::pow(calcPropsRaw.ionicConcentrations->at(0), 2)) + sum;

		for (size_t idx = 0; idx < NCO; idx++) {
			const SysComp::Constituent *c = systemPack.constituents.at(idx).internalConstituent;
			const ERVector &dissocDegreesForCtuent = dissocDegrees.at(idx);
			double t = 0.0;

			for (size_t ifIdx = 0; ifIdx < c->ionicForms->size(); ifIdx++) {
				const SysComp::IonicForm *iF = c->ionicForms->at(ifIdx);

				t += -iF->totalCharge * dissocDegreesForCtuent.at(ifIdx) / denominator;
			}

			termD.emplace_back(t);
		}

		return termD;
	};

	/* You should want to kill yourself already right about now but we are not done yet... */

	const ERVector termA = calcTermA();		/* Zi1 in PeakMaster 5.3 */
	const ERVector termAz = calcTermAz();		/* Zk1 in PeakMaster 5.3 */
	const ERVector termB = calcTermB();		/* Zi2 in PeakMaster 5.3 */
	const ERVector termC = calcTermC(termB);	/* Zi4 in PeakMaster 5.3 */
	const ERVector termD = calcTermD();		/* Zk3 in PeakMaster 5.3 */

	/* Calculate the diffusion matrix */
	for (size_t row = 0; row < NCO; row++) {
		const SysComp::Constituent *c = systemPack.constituents.at(row).internalConstituent;

		diffMatrix(row, row) = termA.at(row);

		const double anC = analyticalConcentrations->at(c->analyticalConcentrationIndex);
		for (size_t col = 0; col < NCO; col++)
			diffMatrix(row, col) += -(anC / kappa) * termB.at(row) * termAz.at(col) + termC.at(row) * termD.at(col);
	}

	_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::CALC_NONLIN_DIFF_MATRIX, const EMMatrix&>(diffMatrix);

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

EigenzoneDispersionVec calculateNonlinear(const CalculatorSystemPack &systemPack, const RealVecPtr &analyticalConcentrations,
					  const DeltaPackVec &deltaPacks, const ConcentrationDeltasVec &concentrationDeltasVec,
					  const EMMatrix &M1, const EMMatrix &M2, const QLQRPack &QLQR,
					  const NonidealityCorrections corrections)
{
	_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::CALC_NONLIN_PROGRESS, const char*>("Starting");

	const EMMatrixVec M1Derivatives = calculateM1Derivatives(systemPack, deltaPacks);
	const EMMatrixVec M2Derivatives = calculateM2Derivatives(systemPack, analyticalConcentrations, corrections);

	_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::CALC_NONLIN_PROGRESS, const char*>("Individual matrix derivatives solved");

	const EMMatrix diffMatrix = makeDiffusionMatrix(systemPack, analyticalConcentrations, concentrationDeltasVec);
	const EMMatrixVec MDerivatives = calculateMDerivatives(M1, M2, M1Derivatives, M2Derivatives);
	const EMMatrix deltaCVec = makeConcentrationDeltas(systemPack);

	return calculateEigenzoneDispersion(QLQR, MDerivatives, deltaCVec, diffMatrix, systemPack.constituents.size());
}

} // namespace Calculator
} // namespace LEMNG

#ifndef ECHMET_TRACER_DISABLE_TRACING

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_NONLIN_PROGRESS, "Nonlinear calculations progress reports")
ECHMET_MAKE_LOGGER(LEMNGTracing, CALC_NONLIN_PROGRESS, const char *stage)
{
	return std::string{"Nonlinear calculations stage: "} + std::string{stage};
}

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_NONLIN_NEIGHBOUR_FORMS_LOOKUP, "Neighbour ionic forms lookup for Nernst-Einstein")
ECHMET_MAKE_LOGGER(LEMNGTracing, CALC_NONLIN_NEIGHBOUR_FORMS_LOOKUP, const char *iF, const char *target)
{
	std::ostringstream ss{};

	ss << "Testing IF " << iF << " for target " << target;

	return ss.str();
}

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_NONLIN_NERNST_EINST_INPUT, "Mobility and total charge of ionic form used in Nermst-Einstein equation")
ECHMET_MAKE_LOGGER(LEMNGTracing, CALC_NONLIN_NERNST_EINST_INPUT, const double &mobility, const int32_t &totalCharge)
{
	std::ostringstream ss{};

	ss << "u " << mobility << ", TC " << totalCharge;

	return ss.str();
}

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_NONLIN_DIFFCOEFF_TERM_C_SUM, "Sum used to calculate \"Term C\" in diffusion coefficient equation")
ECHMET_MAKE_LOGGER(LEMNGTracing, CALC_NONLIN_DIFFCOEFF_TERM_C_SUM, const double &sum)
{
	std::ostringstream ss{};

	ss << "Term C Sum " << sum;

	return ss.str();
}

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_NONLIN_DIFFCOEFF_TERM_C_V, "\"v\" value used to calculate \"Term C\" in diffusion coefficient equation")
ECHMET_MAKE_LOGGER(LEMNGTracing, CALC_NONLIN_DIFFCOEFF_TERM_C_V, const double &v)
{
	std::ostringstream ss{};

	ss << "Term C v " << v;

	return ss.str();
}

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_NONLIN_DIFFCOEFF_TERM_C_INTER2, "Intermediate values used to calculate \"Term C\" in diffusion coefficient equation")
ECHMET_MAKE_LOGGER(LEMNGTracing, CALC_NONLIN_DIFFCOEFF_TERM_C_INTER2, const double &diffCoeff, const double &dDiffCoeffDX)
{
	std::ostringstream ss{};

	ss << "diffCoeff " << diffCoeff << ", dDiffCoeffDX " << dDiffCoeffDX;

	return ss.str();
}

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_NONLIN_DIFFCOEFF_TERM_D_SUM, "Sum used to calculate \"Term D\" in diffusion coefficient equation")
ECHMET_MAKE_LOGGER(LEMNGTracing, CALC_NONLIN_DIFFCOEFF_TERM_D_SUM, const double &sum)
{
	std::ostringstream ss{};

	ss << "Term D Sum " << sum;

	return ss.str();
}

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_NONLIN_DIFF_MATRIX, "Diffusion matrix")
ECHMET_MAKE_LOGGER(LEMNGTracing, CALC_NONLIN_DIFF_MATRIX, const ECHMET::LEMNG::Calculator::EMMatrix &diffMatrix)
{
	std::ostringstream ss{};

	ss << "-- Diffusion matrix --\n"
	   << "---\n\n" << diffMatrix << "\n\n---";

	return ss.str();
}

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_NONLIN_DIFF_PARAMS_MATRIX, "Diffusive parameters matrix")
ECHMET_MAKE_LOGGER(LEMNGTracing, CALC_NONLIN_DIFF_PARAMS_MATRIX, const ECHMET::LEMNG::Calculator::EMMatrixC &dpMatrix)
{
	std::ostringstream ss{};

	ss << "-- Diffusive parameters matrix --\n"
	   << "---\n\n" << dpMatrix << "\n\n---";

	return ss.str();
}

#endif // ECHMET_TRACER_DISABLE_TRACING

} // namespace ECHMET
