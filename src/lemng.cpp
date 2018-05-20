#include "lemng_p.h"
#include "calculator_common.h"
#include "calculator_linear.h"
#include "calculator_nonlinear.h"
#include "helpers.h"
#include "results_maker.h"
#include <new>

#define USE_ECHMET_CONTAINERS
#include <containers/echmetskmap_p.h>
#include <containers/echmetvec_p.h>

#include "tracing/lemng_tracer_impl.h"
#include <sstream>

#define _STRINGIFY(input) #input
#define ERROR_CODE_CASE(erCase) case RetCode::erCase: return _STRINGIFY(erCase)

#ifndef ECHMET_TRACER_DISABLE_TRACING

ECHMET_MAKE_TRACER(LEMNGTracing)

#else // ECHMET_TRACER_DISABLE_TRACING

ECHMET_MAKE_TRACER(ECHMET::__DUMMY_TRACER_CLASS)

#endif // ECHMET_TRACER_DISABLE_TRACING

namespace ECHMET {
namespace LEMNG {

class ConcentrationTooLowException : public std::runtime_error {
public:
	using std::runtime_error::runtime_error;
};

class CannotApplyConcentrationException : public std::runtime_error {
public:
	CannotApplyConcentrationException() : std::runtime_error{""}
	{}
};

CZESystemImpl::CZESystemImpl() :
	m_chemicalSystemBGE{std::unique_ptr<SysComp::ChemicalSystem, decltype(&chemicalSystemDeleter)>{new SysComp::ChemicalSystem, &chemicalSystemDeleter}},
	m_chemicalSystemFull{std::unique_ptr<SysComp::ChemicalSystem, decltype(&chemicalSystemDeleter)>{new SysComp::ChemicalSystem, &chemicalSystemDeleter}},
	m_calcPropsBGE{std::unique_ptr<SysComp::CalculatedProperties, decltype(&calculatedPropertiesDeleter)>{new SysComp::CalculatedProperties, &calculatedPropertiesDeleter}},
	m_calcPropsFull{std::unique_ptr<SysComp::CalculatedProperties, decltype(&calculatedPropertiesDeleter)>{new SysComp::CalculatedProperties, &calculatedPropertiesDeleter}}
{
}

CZESystemImpl::CZESystemImpl(CZESystemImpl &&other) noexcept :
	m_chemicalSystemBGE{std::move(other.m_chemicalSystemBGE)},
	m_chemicalSystemFull{std::move(other.m_chemicalSystemFull)},
	m_calcPropsBGE{std::move(other.m_calcPropsBGE)},
	m_calcPropsFull{std::move(other.m_calcPropsFull)},
	m_systemPack{std::move(other.m_systemPack)},
	m_isAnalyteMap{std::move(other.m_isAnalyteMap)}
{
}

CZESystemImpl::CZESystemImpl(const SysComp::ChemicalSystem &chemicalSystemBGE, const SysComp::CalculatedProperties &calcPropsBGE, const SysComp::ChemicalSystem &chemicalSystemFull, const SysComp::CalculatedProperties &calcPropsFull, const IsAnalyteMap &iaMap) :
	m_chemicalSystemBGE{std::unique_ptr<SysComp::ChemicalSystem, decltype(&chemicalSystemDeleter)>{new SysComp::ChemicalSystem, &chemicalSystemDeleter}},
	m_chemicalSystemFull{std::unique_ptr<SysComp::ChemicalSystem, decltype(&chemicalSystemDeleter)>{new SysComp::ChemicalSystem, &chemicalSystemDeleter}},
	m_calcPropsBGE{std::unique_ptr<SysComp::CalculatedProperties, decltype(&calculatedPropertiesDeleter)>{new SysComp::CalculatedProperties, &calculatedPropertiesDeleter}},
	m_calcPropsFull{std::unique_ptr<SysComp::CalculatedProperties, decltype(&calculatedPropertiesDeleter)>{new SysComp::CalculatedProperties, &calculatedPropertiesDeleter}},
	m_isAnalyteMap{iaMap}
{
	setupInternal(chemicalSystemBGE, calcPropsBGE, chemicalSystemFull, calcPropsFull);
}

CZESystemImpl::CZESystemImpl(const SysComp::ChemicalSystem &chemicalSystemBGE, const SysComp::CalculatedProperties &calcPropsBGE, const SysComp::ChemicalSystem &chemicalSystemFull, const SysComp::CalculatedProperties &calcPropsFull, IsAnalyteMap &&iaMap) :
	m_chemicalSystemBGE{std::unique_ptr<SysComp::ChemicalSystem, decltype(&chemicalSystemDeleter)>{new SysComp::ChemicalSystem, &chemicalSystemDeleter}},
	m_chemicalSystemFull{std::unique_ptr<SysComp::ChemicalSystem, decltype(&chemicalSystemDeleter)>{new SysComp::ChemicalSystem, &chemicalSystemDeleter}},
	m_calcPropsBGE{std::unique_ptr<SysComp::CalculatedProperties, decltype(&calculatedPropertiesDeleter)>{new SysComp::CalculatedProperties, &calculatedPropertiesDeleter}},
	m_calcPropsFull{std::unique_ptr<SysComp::CalculatedProperties, decltype(&calculatedPropertiesDeleter)>{new SysComp::CalculatedProperties, &calculatedPropertiesDeleter}},
	m_isAnalyteMap{iaMap}
{
	setupInternal(chemicalSystemBGE, calcPropsBGE, chemicalSystemFull, calcPropsFull);
}

CZESystemImpl::~CZESystemImpl() noexcept
{
}

RetCode ECHMET_CC CZESystemImpl::evaluate(const InAnalyticalConcentrationsMap *acBGE, const InAnalyticalConcentrationsMap *acSample,
					  const NonidealityCorrections corrections, Results &results) noexcept
{
	auto applyConcentrationMapping = [](RealVecPtr &acVec, const InAnalyticalConcentrationsMap *acMap, const ChemicalSystemPtr &chemSystem) {
		InAnalyticalConcentrationsMap::Iterator *it = acMap->begin();
		if (it == nullptr)
			throw CannotApplyConcentrationException{};

		 while (it->hasNext()) {
			const char *name = it->key();
			const double cAc = it->value();
			size_t idx;

			if (cAc < minimumSafeConcentration()) {
				it->destroy();
				throw ConcentrationTooLowException{name};
			}


			if (chemSystem->analyticalConcentrationsByName->at(idx, name) != ::ECHMET::RetCode::OK) {
				it->destroy();
				throw CannotApplyConcentrationException{};
			}

			(*acVec.get())[idx] = cAc;

			it->next();
		}
		it->destroy();
	};

	auto applyConcentrationMappingBGELike = [&acSample, &acBGE, this](RealVecPtr &acVec) {
		/* Map analytical concentrations from the BGE to full system.
		 * Concentrations of BGE components are the same as in the plain BGE,
		 * concentrations of analytes are "very small".
		 *
		 * This idea here is to solve an almost-like-BGE system so that
		 * we can account for ionic strength effects on the mobility of
		 * the analytes without having the analytes affect the overall
		 * properties of the system.
		 */
		InAnalyticalConcentrationsMap::Iterator *it = acSample->begin();
		if (it == nullptr)
			throw CannotApplyConcentrationException{};

		while (it->hasNext()) {
			const char *name = it->key();
			size_t idx ;

			if (m_chemicalSystemFull->analyticalConcentrationsByName->at(idx, name) != ::ECHMET::RetCode::OK) {
				it->destroy();
				throw CannotApplyConcentrationException{};
			}

			if (isAnalyte(name)) {
				(*acVec.get())[idx] = Calculator::ANALYTE_CONCENTRATION; /* This is our "very small" concentration */
			} else {
				double cAc;
				if (acBGE->at(cAc, name) != ::ECHMET::RetCode::OK) {
					it->destroy();
					throw CannotApplyConcentrationException{};
				}

				if (cAc < minimumSafeConcentration()) {
					it->destroy();
					throw ConcentrationTooLowException{name};
				}

				(*acVec.get())[idx] = cAc;
			}

			it->next();
		}
		it->destroy();
	};

	Calculator::DeltaPackVec deltaPacks{};				/* LEMNG ordering */
	Calculator::ConcentrationDeltasVec concentrationDeltasVec{};	/* SysComp ordering */
	/* Initialize vectors of concentrations */
	RealVecPtr analConcsBGE{nullptr, echmetRealVecDeleter};
	RealVecPtr analConcsBGELike{nullptr, echmetRealVecDeleter};
	RealVecPtr analConcsFull{nullptr, echmetRealVecDeleter};
	try {
		analConcsBGE = makeAnalyticalConcentrationsVec(m_chemicalSystemBGE);
		analConcsBGELike = makeAnalyticalConcentrationsVec(m_chemicalSystemFull);
		analConcsFull = makeAnalyticalConcentrationsVec(m_chemicalSystemFull);
	} catch (std::bad_alloc &) {
		m_lastErrorString = "Cannot make vectors of analytical concentrations";
		_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::EVAL_INIT_ERR, const char*, const char*>("Cannot make vectors of analytical concentrations", "Insufficient memory");

		return RetCode::E_NO_MEMORY;
	}

	try {
		applyConcentrationMapping(analConcsBGE, acBGE, m_chemicalSystemBGE);
		applyConcentrationMapping(analConcsFull, acSample, m_chemicalSystemFull);
		applyConcentrationMappingBGELike(analConcsBGELike);
	} catch (const CannotApplyConcentrationException &ex) {
		m_lastErrorString = "Cannot process input analytical concentrations";
		_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::EVAL_INIT_ERR, const char*, const char*>("Cannot process input analytical concentrations", "Malformed input data");

		return RetCode::E_INTERNAL_ERROR;
	} catch (const ConcentrationTooLowException &ex) {
		m_lastErrorString = "Concentration of " + std::string{ex.what()} + " is too low for the numerical solver";
		_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::EVAL_INIT_ERR, const char*, const char*>("Cannot process input analytical concentrations", "Concentration too low");

		return RetCode::E_CONCENTRATION_TOO_LOW;
	}

	auto isAnalyteFunc = [this](const std::string &s){ return this->isAnalyte(s); };

	/* Prepare output results */
	try {
		results = prepareResults(m_chemicalSystemBGE, m_chemicalSystemFull, isAnalyteFunc);
	} catch (std::bad_alloc &) {
		m_lastErrorString = "Insufficient memory to prepare results";
		_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::EVAL_INIT_ERR, const char*, const char*>("Cannot prepare Results data structures", "Insufficient memory");

		return RetCode::E_NO_MEMORY;
	}

	Calculator::SolutionProperties BGEProps;
	try {
		BGEProps = Calculator::calculateSolutionProperties(m_chemicalSystemBGE, analConcsBGE, m_calcPropsBGE, corrections, true);
	} catch (const Calculator::CalculationException &ex) {
		releaseResults(results);
		m_lastErrorString = std::string{"Unable to calculate BGE properties: "} + ex.what();
		_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::EVAL_PROGRESS_ERR, const char*, const char*>("Unable to calculate BGE properties", ex.what());

		return RetCode::E_CANNOT_SOLVE_BGE;
	}

	m_systemPack.conductivity = BGEProps.conductivity;

	Calculator::SolutionProperties BGELikeProps;
	/* Precalculate what is used in many places of the linear model */
	try {
		Calculator::prepareModelData(m_systemPack, deltaPacks, concentrationDeltasVec, analConcsBGELike, analConcsFull, BGELikeProps, corrections);
	} catch (std::bad_alloc &) {
		fillResultsBGE(m_chemicalSystemBGE, BGEProps, corrections, results);
		return RetCode::E_NO_MEMORY;
	} catch (Calculator::CalculationException &ex) {
		fillResultsBGE(m_chemicalSystemBGE, BGEProps, corrections, results);
		m_lastErrorString = ex.what();
		_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::EVAL_PROGRESS_ERR, const char*, const char*>("Cannot prepare model data", ex.what());

		return ex.errorCode();
	}

	/* Solve the linear model and first nonlinearity term */
	try {
		Calculator::LinearResults linResults = Calculator::calculateLinear(m_systemPack, deltaPacks, corrections);
		Calculator::EigenzoneDispersionVec ezDisps = Calculator::calculateNonlinear(m_systemPack, analConcsBGELike, deltaPacks, concentrationDeltasVec,
											    linResults.M1, linResults.M2, linResults.QLQR, corrections);

		fillResults(m_chemicalSystemBGE, m_chemicalSystemFull, BGEProps, BGELikeProps, linResults, ezDisps, corrections, results);
		/* There was a TODO for a case where not all zones were valid. Be prepared to
		 * revisit this in case we need to handle this somehow in the future. */
	} catch (std::bad_alloc &) {
		fillResultsBGE(m_chemicalSystemBGE, BGEProps, corrections, results);
		fillResultsAnalytesDissociation(m_chemicalSystemFull, BGELikeProps, results);
		_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::EVAL_PROGRESS_ERR, const char*, const char*>("Cannot evaluate linear model", "Insufficient memory");

		return RetCode::E_NO_MEMORY;
	} catch (Calculator::CalculationException &ex) {
		fillResultsBGE(m_chemicalSystemBGE, BGEProps, corrections, results);
		fillResultsAnalytesDissociation(m_chemicalSystemFull, BGELikeProps, results);
		m_lastErrorString = ex.what();
		_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::EVAL_PROGRESS_ERR, const char*, const char*>("Cannot evaluate linear model", ex.what());

		return ex.errorCode();
	}

	return RetCode::OK;
}

bool CZESystemImpl::isAnalyte(const std::string &name)
{
	return m_isAnalyteMap.at(name);
}

const char * ECHMET_CC CZESystemImpl::lastErrorString() const noexcept
{
	return m_lastErrorString.c_str();
}

CZESystemImpl * CZESystemImpl::make(const SysComp::InConstituentVec *inCtuentVecBGE, const SysComp::InConstituentVec *inCtuentVecSample)
{
	SysComp::ChemicalSystem chemSystemBGE{};
	SysComp::ChemicalSystem chemSystemFull{};
	SysComp::CalculatedProperties calcPropsBGE{};
	SysComp::CalculatedProperties calcPropsFull{};

	::ECHMET::RetCode tRet = SysComp::makeComposition(chemSystemBGE, calcPropsBGE, inCtuentVecBGE);
	if (tRet != ::ECHMET::RetCode::OK)
		throw SysCompException{"Cannot make BGE system composition", tRet};

	tRet = SysComp::makeComposition(chemSystemFull, calcPropsFull, inCtuentVecSample);
	if (tRet != ::ECHMET::RetCode::OK) {
		SysComp::releaseChemicalSystem(chemSystemBGE);
		throw SysCompException{"Cannot make full system composition", tRet};
	}

	try {
		IsAnalyteMap iaMap = makeIsAnalyteMap(inCtuentVecBGE, inCtuentVecSample);

		return new CZESystemImpl(chemSystemBGE, calcPropsBGE, chemSystemFull, calcPropsFull, std::move(iaMap));
	} catch (std::bad_alloc &up) {
		SysComp::releaseChemicalSystem(chemSystemBGE);
		SysComp::releaseChemicalSystem(chemSystemFull);
		throw up;
	}
}

void ECHMET_CC CZESystemImpl::toggleAllTracepoints(const bool state) noexcept
{
	if (state)
		TRACER_INSTANCE<LEMNGTracing>().enableAllTracepoints();
	else
		TRACER_INSTANCE<LEMNGTracing>().disableAllTracepoints();
}

void ECHMET_CC CZESystemImpl::toggleTracepoint(const int32_t TPID, const bool state) noexcept
{
	if (state)
		TRACER_INSTANCE<LEMNGTracing>().enableTracepoint(TPID);
	else
		TRACER_INSTANCE<LEMNGTracing>().disableTracepoint(TPID);
}

FixedString * ECHMET_CC CZESystemImpl::trace(const bool dontClear) noexcept
{
#ifdef ECHMET_TRACER_DISABLE_TRACING
	(void)dontClear;
	return nullptr;
#else
	return createFixedString(TRACER_INSTANCE<LEMNGTracing>().logged(dontClear).c_str());
#endif // ECHMET_TRACER_DISABLE_TRACING
}

TracepointInfoVec * ECHMET_CC CZESystemImpl::tracepointInfo() const noexcept
{
#ifdef ECHMET_TRACER_DISABLE_TRACING
	return nullptr;
#else // ECHMET_TRACER_DISABLE_TRACING
	VecImpl<TracepointInfo, false> *tpiVec = createECHMETVec<TracepointInfo, false>(0);
	if (tpiVec == nullptr)
		return nullptr;

	auto &tpiVecSTL = tpiVec->STL();
	auto &tracerInstance = TRACER_INSTANCE<LEMNGTracing>();
	auto tracepoints = tracerInstance.tracepoints();

	try {
		for (auto &tp : tracepoints) {
			FixedString *desc = createFixedString(std::get<1>(tp).c_str());
			if (desc == nullptr)
				throw std::bad_alloc{};
			tpiVecSTL.emplace_back(TracepointInfo{std::get<0>(tp), desc});
		}
	} catch (std::bad_alloc &) {
		tpiVec->destroy();
		return nullptr;
	}

	return tpiVec;
#endif // ECHMET_TRACER_DISABLE_TRACING
}

bool ECHMET_CC CZESystemImpl::tracepointState(const int32_t TPID) const noexcept
{
	return TRACER_INSTANCE<LEMNGTracing>().isTracepointEnabled(TPID);
}

RetCode ECHMET_CC CZESystemImpl::makeAnalyticalConcentrationsMaps(InAnalyticalConcentrationsMap *&acMapBGE, InAnalyticalConcentrationsMap *&acMapFull) const noexcept
{
	typedef std::unique_ptr<MutSKMapImpl<double>> InAnalyticalConcentrationsMapPtr;

	auto fillMap = [](InAnalyticalConcentrationsMapPtr &acMap, const SysComp::ConstituentVec *cVec) {
		SKMapImpl<double>::STLMap &stlMap = acMap->STL();

		for (size_t idx = 0; idx < cVec->size(); idx++) {
			const SysComp::Constituent *c = cVec->at(idx);
			const std::string name{c->name->c_str()};

			stlMap.emplace(name, 0.0);
		}
	};

	try {
		InAnalyticalConcentrationsMapPtr acMapBGEPtr{new MutSKMapImpl<double>{}};
		InAnalyticalConcentrationsMapPtr acMapFullPtr{new MutSKMapImpl<double>{}};

		fillMap(acMapBGEPtr, m_chemicalSystemBGE->constituents);
		fillMap(acMapFullPtr, m_chemicalSystemFull->constituents);

		acMapBGE = acMapBGEPtr.release();
		acMapFull = acMapFullPtr.release();

		return RetCode::OK;
	} catch (std::bad_alloc &) {
		acMapBGE = nullptr;
		acMapFull = nullptr;

		return RetCode::E_NO_MEMORY;
	}
}

void CZESystemImpl::setupInternal(const SysComp::ChemicalSystem &chemicalSystemBGE, const SysComp::CalculatedProperties &calcPropsBGE, const SysComp::ChemicalSystem &chemicalSystemFull, const SysComp::CalculatedProperties &calcPropsFull)
{
	m_chemicalSystemBGE->constituents = chemicalSystemBGE.constituents;
	m_chemicalSystemBGE->ionicForms = chemicalSystemBGE.ionicForms;
	m_chemicalSystemBGE->analyticalConcentrationsByName = chemicalSystemBGE.analyticalConcentrationsByName;
	m_chemicalSystemBGE->ionicConcentrationsByName = chemicalSystemBGE.ionicConcentrationsByName;
	m_chemicalSystemBGE->effectiveMobilitiesByName = chemicalSystemBGE.effectiveMobilitiesByName;
	m_chemicalSystemBGE->ionicMobilitiesByName = chemicalSystemBGE.ionicMobilitiesByName;

	m_calcPropsBGE->ionicConcentrations = calcPropsBGE.ionicConcentrations;
	m_calcPropsBGE->ionicMobilities = calcPropsBGE.ionicMobilities;
	m_calcPropsBGE->effectiveMobilities = calcPropsBGE.effectiveMobilities;
	m_calcPropsBGE->ionicStrength = 0;
	m_calcPropsBGE->conductivity = 0;

	m_chemicalSystemFull->constituents = chemicalSystemFull.constituents;
	m_chemicalSystemFull->ionicForms = chemicalSystemFull.ionicForms;
	m_chemicalSystemFull->analyticalConcentrationsByName = chemicalSystemFull.analyticalConcentrationsByName;
	m_chemicalSystemFull->ionicConcentrationsByName = chemicalSystemFull.ionicConcentrationsByName;
	m_chemicalSystemFull->effectiveMobilitiesByName = chemicalSystemFull.effectiveMobilitiesByName;
	m_chemicalSystemFull->ionicMobilitiesByName = chemicalSystemFull.ionicMobilitiesByName;

	m_calcPropsFull->ionicConcentrations = calcPropsFull.ionicConcentrations;
	m_calcPropsFull->ionicMobilities = calcPropsFull.ionicMobilities;
	m_calcPropsFull->effectiveMobilities = calcPropsFull.effectiveMobilities;
	m_calcPropsFull->ionicStrength = 0;
	m_calcPropsFull->conductivity = 0;

	m_systemPack = Calculator::makeSystemPack(m_chemicalSystemFull, m_calcPropsFull, [this](const std::string &s) { return this->isAnalyte(s); });
}

const char * ECHMET_CC LEMNGerrorToString(const RetCode tRet) noexcept
{
	switch (tRet) {
		ERROR_CODE_CASE(OK);
		ERROR_CODE_CASE(E_NO_MEMORY);
		ERROR_CODE_CASE(E_INVALID_ARGUMENT);
		ERROR_CODE_CASE(E_INVALID_CAPILLARY);
		ERROR_CODE_CASE(E_INVALID_DETECTOR_POSITION);
		ERROR_CODE_CASE(E_NOT_IMPLEMENTED);
		ERROR_CODE_CASE(E_CANNOT_SOLVE_BGE);
		ERROR_CODE_CASE(E_DATA_TOO_LARGE);
		ERROR_CODE_CASE(E_INVALID_CONSTITUENT);
		ERROR_CODE_CASE(E_INVALID_COMPLEXATION);
		ERROR_CODE_CASE(E_DUPLICIT_CONSTITUENTS);
		ERROR_CODE_CASE(E_UNKW_CORELIBS_ERROR);
		ERROR_CODE_CASE(E_CHEM_SYSTEM_UNSOLVABLE);
		ERROR_CODE_CASE(E_INTERNAL_ERROR);
		ERROR_CODE_CASE(E_COMPLEX_EIGENMOBILITIES);
	default:
		return "Unknown error code";
	}
}

RetCode ECHMET_CC makeCZESystem(SysComp::InConstituentVec *BGE, SysComp::InConstituentVec *sample, CZESystem *&czeSystem) noexcept
{
	try {
		czeSystem = CZESystemImpl::make(BGE, sample);
	} catch (std::bad_alloc &) {
		return RetCode::E_NO_MEMORY;
	} catch (SysCompException &ex) {
		_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::MAKE_CZE_SYSTEM_ERR, const char*>(ex.what());

		return coreLibsErrorToNativeError(ex.errorCode());
	}

	return RetCode::OK;
}

double ECHMET_CC minimumSafeConcentration() noexcept
{
	return Calculator::ANALYTE_CONCENTRATION * 10.0;
}

void ECHMET_CC releaseCZESystem(const CZESystem *czeSystem) noexcept
{
	const CZESystemImpl *czeSystemImpl{};

	if (czeSystem == nullptr)
		return;

	czeSystemImpl = dynamic_cast<const CZESystemImpl *>(czeSystem);
	if (czeSystemImpl != nullptr)
		delete czeSystemImpl;
}

CZESystem::~CZESystem() noexcept {}

} // namespace LEMNG

#ifndef ECHMET_TRACER_DISABLE_TRACING

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, MAKE_CZE_SYSTEM_ERR, "makeCZEsystemError")
ECHMET_MAKE_LOGGER(LEMNGTracing, MAKE_CZE_SYSTEM_ERR, const char *err)
{
	return std::string{err};
}

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, EVAL_INIT_ERR, "Evaluation initialization error")
ECHMET_MAKE_LOGGER(LEMNGTracing, EVAL_INIT_ERR, const char *where, const char *why)
{
	std::ostringstream ss{};
	ss << "Error during evaluation initialization: " << where << " (" << why << ")";

	return ss.str();
}

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, EVAL_PROGRESS_ERR, "Error during evaluation")
ECHMET_MAKE_LOGGER(LEMNGTracing, EVAL_PROGRESS_ERR, const char *where, const char *why)
{
	std::ostringstream ss{};
	ss << "Error during evaluation: " << where << " (" << why << ")";

	return ss.str();
}

#endif // ECHMET_TRACER_DISABLE_TRACING

} // namespace ECHMET
