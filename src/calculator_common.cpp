#include "calculator_common.h"
#include "helpers.h"
#include <vector>
#include <cassert>
#include <iterator>
#include <future>

#include "tracing/lemng_tracer_impl.h"
#include <sstream>

#ifndef ECHMET_IMPORT_INTERNAL
#define ECHMET_IMPORT_INTERNAL
#endif // ECHMET_IMPORT_INTERNAL
#include <echmetcaes_extended.h>
#include <echmetionprops.h>

namespace ECHMET {
namespace LEMNG {
namespace Calculator {

const size_t SOLVER_MAX_ITERATIONS = 5000;

CalculationException::CalculationException(const std::string &message, const RetCode errorCode) : std::exception{},
    m_errorCode{errorCode},
    m_message{message}
{
}

CalculationException::CalculationException(std::string &&message, const RetCode errorCode) noexcept : std::exception{},
    m_errorCode{errorCode},
    m_message{message}
{
}

RetCode CalculationException::errorCode() const noexcept
{
    return m_errorCode;
}

const char * CalculationException::what() const noexcept
{
	return m_message.c_str();
}

void buildSystemPackVectors(CalculatorConstituentVec &ccVec, CalculatorIonicFormVec &ifVec, const std::vector<const SysComp::Constituent *> &allConstituents,
			    const SysComp::IonicForm *internalIFH3O, const SysComp::IonicForm *internalIFOH,
			    const std::function<bool (const std::string &)> &isAnalyte)
{

	auto findInIfVec = [&ifVec](const char *name) -> CalculatorIonicForm * {
		const std::string _name{name};

		for (const auto iF : ifVec)
			if (iF->name == _name) return iF;

		return nullptr;
	};

	assert(internalIFH3O->ifType == SysComp::IonicFormType::H);
	assert(internalIFOH->ifType == SysComp::IonicFormType::OH);

	/* WARNING: Raw pointers ahead - use your head here! */
	for (size_t idx = 0; idx < allConstituents.size(); idx++) {
		const SysComp::Constituent *ctuent = allConstituents.at(idx);
		CalculatorIonicFormVec locIfVec{};
		const bool ctuentIsAnalyte = isAnalyte(ctuent->name->c_str());

		locIfVec.reserve(ctuent->ionicForms->size());

		for (size_t ifIdx = 0; ifIdx < ctuent->ionicForms->size(); ifIdx++) {
			const SysComp::IonicForm *iF = ctuent->ionicForms->at(ifIdx);
			bool iFisAnalyte = ctuentIsAnalyte;

			/* We do not need to represent uncharged ionic forms in the mobility matrices.
			 * Leaving them out will spare us some unnecessary multiplications.
			 */
			if (iF->totalCharge == 0)
				continue;

			/* We need to build a list of indices of all ligands that are present
			 * in a given ionic form. This is necessary to have a reasonably efficient
			 * function to calculate Kroenecker delta in makeMatrixM1().
			 */
			std::vector<size_t> containedConstituentsList = [&allConstituents, &isAnalyte, &iFisAnalyte](const SysComp::IonicForm *iF) {
				auto findLigandIdx = [&allConstituents](const SysComp::Constituent *c) {
					for (size_t idx = 0; idx < allConstituents.size(); idx++) {
						const SysComp::Constituent *otherC = allConstituents.at(idx);

						if (*(otherC->name) == *c->name)
							return idx;
					}

					throw CalculationException{"Ligand index not found", RetCode::E_INTERNAL_ERROR};
				};

				const SysComp::IonicForm *ancestor = iF;
				std::vector<size_t> idxList{};

				while (ancestor->ligand != nullptr) {
					const size_t idx = findLigandIdx(ancestor->ligand);
					if (isAnalyte(ancestor->ligand->name->c_str()))
						iFisAnalyte = true;

					idxList.emplace_back(idx);
					ancestor = ancestor->ancestor;
				};

				return idxList;
			}(iF);

			containedConstituentsList.emplace_back(idx); /* Ionic form always contains at least the nucleus. */

			/* This is a slightly less obvious part, focus now!
			 * The ifVec contains all *ionic forms* that are present
			 * in the system. Each ionic form has to be contained only once
			 * in this vector. This seems logical, right?
			 * However, each CalculatorConstituent contains its own vector
			 * of ionic forms that contain the said constituent.
			 * If a
			 */

			CalculatorIonicForm *locIF = findInIfVec(iF->name->c_str());
			if (locIF == nullptr) {
				try {
					locIF = new CalculatorIonicForm{std::string{iF->name->c_str()}, iF->totalCharge,
									iF,
									iF->ionicConcentrationIndex,
									ifVec.size(), /* This ionic form will be placed at this index in the global ifVec. */
									std::move(containedConstituentsList), iFisAnalyte};
					ifVec.emplace_back(locIF);
				} catch (std::bad_alloc &) {
					for (auto &&item : locIfVec)
						delete item;

					throw;
				}
			}

			try {
				locIfVec.emplace_back(locIF);
			} catch (std::bad_alloc &) {
				delete locIF;
				for (auto &&item : locIfVec)
					delete item;

				throw;
			}
		}

		locIfVec.shrink_to_fit();

		try {
			ccVec.emplace_back(std::string{ctuent->name->c_str()}, std::move(locIfVec), ctuent, ctuentIsAnalyte);
		} catch (std::bad_alloc &) {
			for (auto &&item : locIfVec)
				delete item;

			throw;
		}
	}

	CalculatorIonicForm *ifH3O = nullptr;
	CalculatorIonicForm *ifOH = nullptr;
	/* H3O+ and OH- forms go at the very end.
	 * The do not correspond to any constituent.
	 * Note that we cannot initialize their mobilities
	 * at this point as they are also effected by
	 * ionic strength.
	 */
	try {
		ifH3O = new CalculatorIonicForm{"H3O+", 1, internalIFH3O, 0, ifVec.size(), {}, false};

		ifVec.emplace_back(ifH3O);
	} catch (std::bad_alloc &) {
		delete ifH3O;
		throw;
	}
	try {
		ifOH = new CalculatorIonicForm{"OH-", -1, internalIFOH, 1, ifVec.size(), {}, false};

		ifVec.emplace_back(ifOH);
	} catch (std::bad_alloc &) {
		delete ifOH;
		throw;
	}
}

void calcIonicProperties(const SysComp::ChemicalSystem *chemSystem, const RealVecPtr &concentrations, SysComp::CalculatedProperties *calcProps, const bool correctForIonicStrength)
{
	::ECHMET::RetCode tRet;
	IonProps::ComputationContext *ctx = IonProps::makeComputationContext(*chemSystem, concentrations.get(), *calcProps);

	if (ctx == nullptr)
		throw CalculationException{"Cannot create IonProps computation context", RetCode::E_NO_MEMORY};

	if (correctForIonicStrength) {
		tRet = IonProps::correctMobilities(ctx);
		if (tRet != ::ECHMET::RetCode::OK) {
			ctx->destroy();
			throw CalculationException{"Cannot correct ionic mobilities for ionic strength: " + std::string(errorToString(tRet)), coreLibsErrorToNativeError(tRet)};
		}
	}
	_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::CALC_COMMON_CALC_SOLPROPS_ION_MOBS, const SysComp::IonicFormVec*, const SysComp::CalculatedProperties*>(chemSystem->ionicForms, calcProps);

	tRet = IonProps::calculateEffectiveMobilities(ctx);
	if (tRet != ::ECHMET::RetCode::OK) {
		ctx->destroy();
		throw CalculationException{"Cannot calculate effective mobilities: " + std::string(errorToString(tRet)), coreLibsErrorToNativeError(tRet)};
	}
	_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::CALC_COMMON_CALC_SOLPROPS_EFF_MOBS, const SysComp::ConstituentVec*, const SysComp::CalculatedProperties*>(chemSystem->constituents, calcProps);

	IonProps::calculateConductivity(ctx);
	_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::CALC_COMMON_CALC_SOLPROPS_CONDUCTIVITY, const double&>(calcProps->conductivity);

	ctx->destroy();
}

double calculateSolutionBufferCapacity(const SysComp::ChemicalSystem *chemSystem, const RealVecPtr &concentrations, const SysComp::CalculatedProperties *calcProps, const bool correctForIonicStrength)
{
	::ECHMET::RetCode tRet;
	double bufferCapacity;

	tRet = ::ECHMET::CAES::calculateBufferCapacity(bufferCapacity, *chemSystem, *calcProps, concentrations.get(), correctForIonicStrength);
	if (tRet != ::ECHMET::RetCode::OK)
		return -1.0;
	return bufferCapacity;
}

SolutionProperties calculateSolutionProperties(const SysComp::ChemicalSystem *chemSystem, const RealVecPtr &concentrations, SysComp::CalculatedProperties *calcProps, const bool correctForIonicStrength,
					       const bool calcBufferCapacity)
{
	ECHMET_TRACE(LEMNGTracing, CALC_COMMON_CALC_SOLPROPS_PROGRESS, ECHMET_S("Solving equilibrium"));
	solveChemicalSystem(chemSystem, concentrations, calcProps, correctForIonicStrength);

	auto analyticalConcentrations = [&concentrations]() {
		std::vector<double> anC{};
		anC.reserve(concentrations->size());
		for (size_t idx = 0; idx < concentrations->size(); idx++)
			anC.emplace_back(concentrations->at(idx));

		return anC;
	}();

	auto ionicConcentrations = [calcProps]() {
		std::vector<double> iC{};
		iC.reserve(calcProps->ionicConcentrations->size());
		for (size_t idx = 0; idx < calcProps->ionicConcentrations->size(); idx++)
			iC.emplace_back(calcProps->ionicConcentrations->at(idx));

		return iC;
	}();

	const double bufferCapacity = [&]() {
		if (calcBufferCapacity)
			return calculateSolutionBufferCapacity(chemSystem, concentrations, calcProps, correctForIonicStrength);
		return -1.0;
	}();

	return SolutionProperties{bufferCapacity,
				  calcProps->conductivity,
				  calcProps->ionicStrength,
				  std::move(analyticalConcentrations),
				  std::move(ionicConcentrations)};
}

SolutionProperties calculateSolutionProperties(const ChemicalSystemPtr &chemSystem, const RealVecPtr &concentrations, CalculatedPropertiesPtr &calcProps, const bool correctForIonicStrength, const bool calcBufferCapacity)
{
	return calculateSolutionProperties(chemSystem.get(), concentrations, calcProps.get(), correctForIonicStrength, calcBufferCapacity);
}

template <>
bool isComplex<EMMatrixC>(const EMMatrixC &I)
{
	for (int row = 0; row < I.rows(); row++) {
		for (int col = 0; col < I.cols(); col++) {
			const std::complex<double> &n = I(row, col);
			if (n.imag() != 0.0)
				return true;
		}
	}

	return false;
}

template <>
bool isComplex<EMVectorC>(const EMVectorC &I)
{
	for (int col = 0; col < I.cols(); col++) {
		const std::complex<double> &n = I(col);
		if (n.imag() != 0.0)
			return true;
	}

	return false;
}

CalculatorSystemPack makeSystemPack(const ChemicalSystemPtr &chemSystem, const CalculatedPropertiesPtr &calcProps,
				    const std::function<bool (const std::string &)> &isAnalyte)
{
	CalculatorConstituentVec ccVec{};
	CalculatorIonicFormVec ifVec{};
	const auto orderedConstituents = sysCompToLEMNGOrdering(chemSystem, isAnalyte);

	ccVec.reserve(orderedConstituents.size());
	ifVec.reserve(chemSystem->ionicForms->size());


	/* WARNING: Raw pointers inside! */
	try {
		buildSystemPackVectors(ccVec, ifVec, orderedConstituents,
				       chemSystem->ionicForms->at(0), chemSystem->ionicForms->at(1),
				       isAnalyte);
	} catch (std::bad_alloc &) {
		for (auto &&item : ifVec)
			delete item;

		throw;
	}

	try {
		return CalculatorSystemPack{std::move(ccVec), std::move(ifVec), chemSystem.get(), calcProps.get()};
	} catch (std::bad_alloc &) {
		for (auto &&item : ifVec)
			delete item;

		throw;
	}
}

void bindSystemPack(CalculatorSystemPack &systemPack, const RealVecPtr &analConcsBGELike, const RealVecPtr &analConcsSample)
{
	const SysComp::ChemicalSystem *chemSystem = systemPack.chemSystemRaw;
	const SysComp::CalculatedProperties *calcProps = systemPack.calcPropsRaw;

	/* Bind ionic forms to state of the fully resolved system */
	for (size_t idx = 0; idx < systemPack.ionicForms.size() - 2; idx++) {
		CalculatorIonicForm *iF = systemPack.ionicForms.at(idx);
		const std::string &name = iF->name;
		size_t mobilityIdx;

		if (chemSystem->ionicMobilitiesByName->at(mobilityIdx, name.c_str()) != ::ECHMET::RetCode::OK)
			throw CalculationException{"Cannot find ionic form mobility index", RetCode::E_INTERNAL_ERROR};
		iF->mobility = ECHMETRealToDouble(calcProps->ionicMobilities->at(mobilityIdx));

		if (iF->isAnalyte)
			iF->concentration = 0.0;
		else {
			size_t concentrationIdx;
			if (chemSystem->ionicConcentrationsByName->at(concentrationIdx, name.c_str()) != ::ECHMET::RetCode::OK)
				throw CalculationException{"Cannot find ionic form concentration index", RetCode::E_INTERNAL_ERROR};

			iF->concentration = ECHMETRealToDouble(calcProps->ionicConcentrations->at(concentrationIdx));
		}
	}

	/* H3O+ and OH- are at the end of the vector in LEMNG */
	const size_t H3O_idx = systemPack.ionicForms.size() - 2;
	const size_t OH_idx = systemPack.ionicForms.size() - 1;
	systemPack.ionicForms[H3O_idx]->mobility = ECHMETRealToDouble(calcProps->ionicMobilities->at(0));		/* H3O+ is always on top of the SysComp vector */
	systemPack.ionicForms[OH_idx]->mobility = ECHMETRealToDouble(calcProps->ionicMobilities->at(1));		/* OH- is always second in the SysComp vector */
	systemPack.ionicForms[H3O_idx]->concentration = ECHMETRealToDouble(calcProps->ionicConcentrations->at(0));	/* H3O+ is always on top of the SysComp vector */
	systemPack.ionicForms[OH_idx]->concentration = ECHMETRealToDouble(calcProps->ionicConcentrations->at(1));	/* OH- is always second in the SysComp vector */

	/* Bind constituents to analytical concentrations */
	for (CalculatorConstituent &cc : systemPack.constituents) {
		const size_t idx = cc.internalConstituent->analyticalConcentrationIndex;

		cc.concentrationSample = analConcsSample->at(idx);
		if (cc.isAnalyte)
			cc.concentrationBGE = 0.0;
		else
			cc.concentrationBGE = analConcsBGELike->at(idx);
	}
}

#ifdef ECHMET_LEMNG_SENSITIVE_NUMDERS
RealVecPtr makeAnalyticalConcentrationsForDerivator(const CalculatorSystemPack &systemPack)
{
	RealVec *v;

	::ECHMET::RetCode tRet = SysComp::makeAnalyticalConcentrationsVec(v, *systemPack.chemSystemRaw);
	if (tRet != ::ECHMET::RetCode::OK)
		throw CalculationException{"Cannot initialize analytical concentrations vector to calculate dervatives", coreLibsErrorToNativeError(tRet)};

	for (const CalculatorConstituent &cc : systemPack.constituents) {
		const size_t idx = cc.internalConstituent->analyticalConcentrationIndex;

		if (cc.isAnalyte)
			(*v)[idx] = ANALYTE_CONCENTRATION_NUMDERS; /* REALLY damn low concentration. Derivator can deal with them whereas the equilibirum solver cannot. */
		else
			(*v)[idx] = cc.concentrationBGE;
	}

	return std::unique_ptr<RealVec, decltype(&echmetRealVecDeleter)>{v, echmetRealVecDeleter};
}
#endif // ECHMET_LEMNG_SENSITIVE_NUMDERS

void precalculateConcentrationDeltas(CalculatorSystemPack &systemPack, DeltaPackVec &deltaPacks, const RealVecPtr &analyticalConcentrations, const bool correctForIonicStrength)
{
	static const ECHMETReal H = DELTA_H;

	const size_t NCO = systemPack.constituents.size();
	const size_t NIF = systemPack.ionicForms.size();
	const size_t H3O_idx = NIF - 2;
	const size_t OH_idx = NIF - 1;
	const SysComp::ChemicalSystem &chemSystemRaw = *systemPack.chemSystemRaw;
#ifdef ECHMET_LEMNG_SENSITIVE_NUMDERS
	(void)analyticalConcentrations;
	const RealVecPtr analyticalConcentrationsForDiffs = makeAnalyticalConcentrationsForDerivator(systemPack);
#else
	const RealVecPtr &analyticalConcentrationsForDiffs = analyticalConcentrations;
#endif // ECHMET_LEMNG_SENSITIVE_NUMDERS

	CAES::Solver *solver = nullptr;
	RealVec *derivatives = nullptr;

	deltaPacks.reserve(NCO);

	::ECHMET::RetCode tRet = CAES::prepareDerivatorContext(derivatives, solver, chemSystemRaw, correctForIonicStrength);
	if (tRet != ::ECHMET::RetCode::OK)
		throw CalculationException{std::string{"Cannot make derivator context: "} + std::string{errorToString(tRet)}, coreLibsErrorToNativeError(tRet)};

#if ECHMET_LEMNG_PARALLEL_NUM_OPS
	auto worker = [&](const SysComp::Constituent *perturbedConstituent) -> DeltaPack {
		const size_t ND = derivatives->size();
		EMVector deltas{NIF};
		ECHMETReal conductivityDerivative;
		RealVec *_derivatives = ::ECHMET::createRealVec(ND);
		if (_derivatives == nullptr)
			throw CalculationException{"Cannot allocate thread-local derivatives vector", RetCode::E_NO_MEMORY};
		if (_derivatives->resize(ND) != ::ECHMET::RetCode::OK) {
			_derivatives->destroy();
			throw CalculationException{"Cannot allocate thread-local derivatives vector", RetCode::E_NO_MEMORY};
		}

		::ECHMET::RetCode tRet = CAES::calculateFirstConcentrationDerivatives_prepared(_derivatives, conductivityDerivative,
											       solver, H,
											       chemSystemRaw, analyticalConcentrationsForDiffs.get(),
											       perturbedConstituent);
		if (tRet != ::ECHMET::RetCode::OK) {
			_derivatives->destroy();
			throw CalculationException{"Cannot calculate concentration derivatives for M2", coreLibsErrorToNativeError(tRet)};
		}

		/* Use LEMNG ordering for concentration deltas */
		for (size_t iFIdx = 0; iFIdx < NIF - 2; iFIdx++) {
			const CalculatorIonicForm *iF = systemPack.ionicForms.at(iFIdx);

			deltas(iFIdx) = ECHMETRealToDouble(_derivatives->at(iF->internalIonicFormConcentrationIdx));
		}

		deltas(H3O_idx) = ECHMETRealToDouble(_derivatives->at(0));
		deltas(OH_idx) = ECHMETRealToDouble(_derivatives->at(1));

		_derivatives->destroy();

		return DeltaPack{std::move(deltas), conductivityDerivative, perturbedConstituent};
	};

	std::vector<std::future<DeltaPack>> results{};

	try {
		results.reserve(NCO);
	} catch (std::bad_alloc &) {
		solver->context()->destroy();
		solver->destroy();
		derivatives->destroy();
		throw CalculationException{"Cannot allocate thread futures", RetCode::E_NO_MEMORY};
	}

	try {
		for (size_t cIdx = 0; cIdx < NCO; cIdx++) {
			const SysComp::Constituent *perturbedConstituent = systemPack.constituents.at(cIdx).internalConstituent;

			results.emplace_back(std::async(std::launch::async, worker, perturbedConstituent));
		}

		for (auto &f : results)
			deltaPacks.emplace_back(f.get());
	} catch (CalculationException &ex) {
		for (auto &f : results) {
			if (f.valid())
				f.wait();
		}

		solver->context()->destroy();
		solver->destroy();
		derivatives->destroy();
		throw ex;
	}
#else // ECHMET_LEMNG_PARALLEL_NUM_OPS
	for (size_t cIdx = 0; cIdx < NCO; cIdx++) {
		const SysComp::Constituent *perturbedConstituent = systemPack.constituents.at(cIdx).internalConstituent;
		ECHMETReal conductivityDerivative;
		EMVector deltas{NIF};

		tRet = CAES::calculateFirstConcentrationDerivatives_prepared(derivatives, conductivityDerivative,
									     solver, H,
									     chemSystemRaw, analyticalConcentrationsForDiffs.get(),
									     perturbedConstituent);
		if (tRet != ::ECHMET::RetCode::OK) {
			solver->context()->destroy();
			solver->destroy();
			derivatives->destroy();
			throw CalculationException{"Cannot calculate concentration derivatives for M2", coreLibsErrorToNativeError(tRet)};
		}

		/* Use LEMNG ordering for concentration deltas */
		for (size_t iFIdx = 0; iFIdx < NIF - 2; iFIdx++) {
			const CalculatorIonicForm *iF = systemPack.ionicForms.at(iFIdx);

			deltas(iFIdx) = ECHMETRealToDouble(derivatives->at(iF->internalIonicFormConcentrationIdx));
		}

		deltas(H3O_idx) = ECHMETRealToDouble(derivatives->at(0));
		deltas(OH_idx) = ECHMETRealToDouble(derivatives->at(1));

		try {
			deltaPacks.emplace_back(std::move(deltas), ECHMETRealToDouble(conductivityDerivative), perturbedConstituent);
		} catch (std::bad_alloc &) {
			solver->context()->destroy();
			solver->destroy();
			derivatives->destroy();
			throw;
		}
	}
#endif // ECHMET_LEMNG_PARALLEL_NUM_OPS
	solver->context()->destroy();
	solver->destroy();
	derivatives->destroy();
}

void prepareModelData(CalculatorSystemPack &systemPack, DeltaPackVec &deltaPacks, const RealVecPtr &analConcsBGELike, const RealVecPtr &analConcsSample, const bool correctForIonicStrength)
{
	/* Step 1 - Identify the target and its flaws, there are always flaws... oops, not this "step one"...
	 *
	 *

	█████████████████████████████████████████
	█████████████████████████████████████████
	████ ▄▄▄▄▄ █▀▄██▄ ▀▀██▀ ▀▀ █▀█ ▄▄▄▄▄ ████
	████ █   █ █  ▀ ▄▀▄▀█ ▄▄█ ▄▀██ █   █ ████
	████ █▄▄▄█ █ ▄██▄ █▄▀█▀ ▀ █ ██ █▄▄▄█ ████
	████▄▄▄▄▄▄▄█▄▀ ▀▄▀ ▀ █ ▀▄█▄▀ █▄▄▄▄▄▄▄████
	████▀█▄▀█▄▄ ▄█▀█▄▄▀▄▄ ▄  █ ▀ ▄▀▄▄  ▄▀████
	█████▀▄   ▄ ▀▀▀▄█▀▀█  █▀ ▄  ▀▄▄ ▀ ▄▄▄████
	████▄  ██▄▄  █▄▄▄██ ▄ ▀▄█▄ ▀▄▀  █▄██▄████
	████▀▀▄▀██▄   █ ▄█▀▄▀██ ██▄ ▄ █▄ ▄█▀ ████
	█████ ▀▀▀ ▄▄▀█▀▄██  ▀ ▄  ██  █  █ ▀█▀████
	█████ ▀▀█▄▄▀▀   ████ ▀▄█ █  ▄  ▀█▀▄▀▄████
	████████ ▄▄██ █▄█▄▄██▄ ▀▄▄▀▄▄▄▄▀████▀████
	████▄▄▀█▀█▄ ▄▀ ▀▀█▄██ ▀▄▀▄▄ ▄▄▄▀▄ ▄▀ ████
	████▄▄▄█▄▄▄▄▀▀▀▀█▄ ▀██▄▀▄▀█▀ ▄▄▄  ▀▀▀████
	████ ▄▄▄▄▄ █ ▄▄▀█ ▀▀▄▀▀█ ▀█  █▄█ ▀█▀ ████
	████ █   █ ████▄▀█ ▀▄▀▄▄█ █   ▄▄  ▀▀█████
	████ █▄▄▄█ █▄ ██▄▄█  ▄█▄█    █▄▄  █▄▄████
	████▄▄▄▄▄▄▄███▄██▄▄▄▄█▄▄██▄▄▄██▄▄▄██▄████
	█████████████████████████████████████████
	█████████████████████████████████████████

	 *
	 * Solve the almost-like-BGE system to get ionic concentrations and corrected ionic mobilities.
	 */
	solveChemicalSystem(systemPack.chemSystemRaw, analConcsBGELike, systemPack.calcPropsRaw, correctForIonicStrength);

	/* Step 2 - Bind the now known properties of the present ionic forms to the SystemPack.
	 */
	bindSystemPack(systemPack, analConcsBGELike, analConcsSample);

	/* Step 3 - Precalculate concentration derivatives */
	precalculateConcentrationDeltas(systemPack, deltaPacks, analConcsBGELike, correctForIonicStrength);
}

void solveChemicalSystem(const SysComp::ChemicalSystem *chemSystem, const RealVecPtr &concentrations, SysComp::CalculatedProperties *calcProps, const bool correctForIonicStrength)
{
	::ECHMET::RetCode tRet;
	CAES::SolverContext *solverCtx;
	CAES::Solver *solver;
	CAES::Solver::Options opts = CAES::Solver::defaultOptions() | CAES::Solver::Options::DISABLE_THREAD_SAFETY;

	if (correctForIonicStrength)
		opts |= CAES::Solver::IONIC_STRENGTH_CORRECTION;

	tRet = CAES::createSolverContext(solverCtx, *chemSystem);
	if (tRet != ::ECHMET::RetCode::OK) {
		throw CalculationException{"Failed to create solver context: " + std::string{errorToString(tRet)}, coreLibsErrorToNativeError(tRet)};
	}

	solver = CAES::createSolver(solverCtx, opts);
	if (solver == nullptr) {
		solverCtx->destroy();
		throw CalculationException{"Failed to create solver", RetCode::E_NO_MEMORY};
	}

	tRet = CAES::estimateDistribution(solver, concentrations.get(), *calcProps);
	if (tRet != ::ECHMET::RetCode::OK) {
		solver->destroy();
		solverCtx->destroy();
		throw CalculationException{"Failed to estimate distribution: " + std::string(errorToString(tRet)), coreLibsErrorToNativeError(tRet)};
	}

	CAES::SolverIterations solvIters{};
	tRet = solver->solve(concentrations.get(), *calcProps, SOLVER_MAX_ITERATIONS, &solvIters);
	if (tRet != ::ECHMET::RetCode::OK) {
		solver->destroy();
		solverCtx->destroy();
		throw CalculationException{"Solver was unable to calculate equilibrium composition: " + std::string(errorToString(tRet)), coreLibsErrorToNativeError(tRet)};
	}

	ECHMET_TRACE(LEMNGTracing, CALC_COMMON_CALC_SOLPROPS_PROGRESS, ECHMET_S("Equilibrium successfuly solved"));
	_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::CALC_COMMON_CALC_SOLPROPS_ITERS, const uint32_t&, const uint32_t&>(solvIters.outer, solvIters.total);
	_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::CALC_COMMON_CALC_SOLPROPS_EQ_COMP, const SysComp::IonicFormVec*, const SysComp::CalculatedProperties*>(chemSystem->ionicForms, calcProps);

	solver->destroy();
	solverCtx->destroy();

	/* Calculate ionic properties */
	calcIonicProperties(chemSystem, concentrations, calcProps, correctForIonicStrength);
}

void solveChemicalSystem(const ChemicalSystemPtr chemSystem, const RealVecPtr &concentrations, CalculatedPropertiesPtr &calcProps, const bool correctForIonicStrength)
{
	return solveChemicalSystem(chemSystem.get(), concentrations, calcProps.get(), correctForIonicStrength);
}

std::vector<const SysComp::Constituent *> sysCompToLEMNGOrdering(const ChemicalSystemPtr &chemSystem, const std::function<bool (const std::string &)> &isAnalyte)
{
	/*
	 * This reorders the constituents from the order they came in into
	 * "background first, analytes next, water last" order.
	 * Such ordering is introduced in a series of articles describing
	 * the original linear theory of electromigration:
	 *
	 * Štědrý M, Jaroš M, Gaš, B, JOURNAL OF CHROMATOGRAPHY A 2002, Volume: 960, Pages: 187-198
	 * Štědrý M, Jaroš M, Včeláková K, Gaš B, ELECTROPHORESIS 2003, Volume: 24, Pages: 536-547 (DOI: 10.1002/elps.200390061)
	 * Štědrý M, Jaroš M, Hruška V, Gaš G, ELECTROPHORESIS 2004, Volume: 25, Pages: 3071–3079 (DOI: 10.1002/elps.200405981)
	 * Jaroš M, Hruška V, Štědrý M, Zusková I, Gaš B, ELECTROPHORESIS 2004, Volume: 25, Pages: 3080-3085 (DOI: 10.1002/elps.200405982)
	 */
	std::vector<const SysComp::Constituent *> all{};
	std::vector<const SysComp::Constituent *> analytes{};

	all.reserve(chemSystem->constituents->size());
	analytes.reserve(chemSystem->constituents->size());

	for (size_t idx = 0; idx < chemSystem->constituents->size(); idx++) {
		const SysComp::Constituent *c = chemSystem->constituents->at(idx);

		if (isAnalyte(c->name->c_str()))
			analytes.emplace_back(c);
		else
			all.emplace_back(c);
	}

	std::copy(analytes.cbegin(), analytes.cend(), std::back_inserter(all));

	return all;
}

} // namespace Calculator
} // namespace LEMNG

#ifndef ECHMET_TRACER_DISABLE_TRACING

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_COMMON_CALC_SOLPROPS_PROGRESS, "Solution properties calculation progress")
ECHMET_MAKE_LOGGER(LEMNGTracing, CALC_COMMON_CALC_SOLPROPS_PROGRESS, const char *msg)
{
	std::ostringstream ss{};

	ss << "Calculating solution properties, stage: " << msg;
	return ss.str();
}

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_COMMON_CALC_SOLPROPS_ITERS, "Iterations needed to calculate concentration equilibrium")
ECHMET_MAKE_LOGGER(LEMNGTracing, CALC_COMMON_CALC_SOLPROPS_ITERS, const uint32_t &outer, const uint32_t &total)
{
	std::ostringstream ss{};

	ss << "Iterations needed to calculate equilibrium: outer (for IS correction): " << outer << ", total (NRS): " << total;
	return ss.str();
}

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_COMMON_CALC_SOLPROPS_EQ_COMP, "Equilibrium composition")
ECHMET_MAKE_LOGGER(LEMNGTracing, CALC_COMMON_CALC_SOLPROPS_EQ_COMP, const ECHMET::SysComp::IonicFormVec *ionicForms, const ECHMET::SysComp::CalculatedProperties *calcProps)
{
	std::ostringstream ss{};

	ss << "-- Equilibium composition --\n";

	for (size_t idx = 0; idx < ionicForms->size(); idx++) {
		const ECHMET::SysComp::IonicForm *iF = ionicForms->at(idx);
		const double c = ECHMET::ECHMETRealToDouble(calcProps->ionicConcentrations->at(iF->ionicConcentrationIndex));

		switch (iF->ifType) {
		case ECHMET::SysComp::IonicFormType::H:
			ss << "[H+]: ";
			break;
		case ECHMET::SysComp::IonicFormType::OH:
			ss << "[OH-]: ";
			break;
		default:
			ss << "[" << iF->name->c_str() << "]: ";
			break;
		}

		ss << c << " (mmol/dm3)\n";
	}

	return ss.str();
}

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_COMMON_CALC_SOLPROPS_ION_MOBS, "Ionic moblities corrected to ionic strength")
ECHMET_MAKE_LOGGER(LEMNGTracing, CALC_COMMON_CALC_SOLPROPS_ION_MOBS, const ECHMET::SysComp::IonicFormVec *ionicForms, const ECHMET::SysComp::CalculatedProperties *calcProps)
{
	std::ostringstream ss{};

	ss << "-- Ionic mobilities --\n";

	for (size_t idx = 0; idx < ionicForms->size(); idx++) {
		const ECHMET::SysComp::IonicForm *iF = ionicForms->at(idx);
		const double ionMob = ECHMET::ECHMETRealToDouble(calcProps->ionicConcentrations->at(iF->ionicMobilityIndex));

		switch (iF->ifType) {
		case ECHMET::SysComp::IonicFormType::H:
			ss << "[H+]: ";
			break;
		case ECHMET::SysComp::IonicFormType::OH:
			ss << "[OH-]: ";
			break;
		default:
			ss << "[" << iF->name->c_str() << "]: ";
			break;
		}

		ss << ionMob << " (m.m/V/s . e-9)\n";
	}

	return ss.str();
}

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_COMMON_CALC_SOLPROPS_EFF_MOBS, "Ionic moblities corrected to ionic strength")
ECHMET_MAKE_LOGGER(LEMNGTracing, CALC_COMMON_CALC_SOLPROPS_EFF_MOBS, const ECHMET::SysComp::ConstituentVec *constituents, const ECHMET::SysComp::CalculatedProperties *calcProps)
{
	std::ostringstream ss{};

	ss << "-- Effective mobilities --\n";

	for (size_t idx = 0; idx < constituents->size(); idx++) {
		const ECHMET::SysComp::Constituent *ctuent =  constituents->at(idx);
		const double effMob = ECHMET::ECHMETRealToDouble(calcProps->effectiveMobilities->at(ctuent->effectiveMobilityIndex));

		ss << "[" << ctuent->name->c_str() << "]: ";
		ss << effMob << " (m.m/V/s . e-9)\n";
	}

	return ss.str();
}

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_COMMON_CALC_SOLPROPS_CONDUCTIVITY, "Solution conductivity")
ECHMET_MAKE_LOGGER(LEMNGTracing, CALC_COMMON_CALC_SOLPROPS_CONDUCTIVITY, const double &conductivity)
{
	std::ostringstream ss{};

	ss << "Conductivity: " << conductivity << " (S/m)";
	return ss.str();
}

#endif // ECHMET_TRACER_DISABLE_TRACING

} // namespace ECHMET
