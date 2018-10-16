#include "lemng_p.h"
#include "calculator_common.h"
#include "calculator_linear.h"
#include "calculator_nonlinear.h"
#include "helpers.h"
#include "results_maker.h"
#include "lemng_config.h"
#include <new>

#define USE_ECHMET_CONTAINERS
#include <containers/echmetskmap_p.h>
#include <containers/echmetvec_p.h>

#include "tracing/lemng_tracer_impl.h"
#include <sstream>

#define _STRINGIFY(input) #input
#define ERROR_CODE_CASE(erCase) case RetCode::erCase: return _STRINGIFY(erCase)
#define MK_VERSION_STRING(vmaj, vmin, vpatch) _STRINGIFY(vmaj) "." _STRINGIFY(vmin) "." _STRINGIFY(vpatch)

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

class InvalidComposition : public std::runtime_error {
public:
	enum Type {
		MISMATCHING_PARAMETERS,
		MISSING_IN_SAMPLE
	};

	InvalidComposition(const Type t) :
		std::runtime_error{"Input composition is invalid"},
		type(t)
	{}

	const Type type;
};

static
void complexesOnlyWithAnalytes(const SysComp::InCFVec *cfVec, const IsAnalyteMap &iaMap)
{
	for (size_t idx = 0; idx < cfVec->size(); idx++) {
		const auto &cf = cfVec->at(idx);

		for (size_t jdx = 0; jdx < cf.ligandGroups->size(); jdx++) {
			const auto &lgg = cf.ligandGroups->at(jdx);

			for (size_t kdx = 0; kdx < lgg.ligands->size(); kdx++) {
				if (!iaMap.at(lgg.ligands->at(kdx).ligandName->c_str()))
					throw InvalidComposition{InvalidComposition::Type::MISMATCHING_PARAMETERS};
			}
		}
	}
}

static
bool findMatchingLigandForm(const SysComp::InLigandForm &ligBGE, const SysComp::InComplexForm &cfSam)
{
	/* Go through all ligand groups and look for the same ligand
	 * with the same charge.
	 * It does not matter if the ligand groups in BGE and sample
	 * constituent are not in the same order. One specific ligand form
	 * must have the same complexation parameters regardless of the
	 * groups it belongs to. Allowing this would go against the laws
	 * of thermodynamics.
	 */

	for (size_t idx = 0; idx < cfSam.ligandGroups->size(); idx++) {
		const auto &lgg = cfSam.ligandGroups->at(idx);

		for (size_t jdx = 0; jdx < lgg.ligands->size(); jdx++) {
			if (SysComp::compareInLigandForms(ligBGE, lgg.ligands->at(jdx)))
				return true;
		}
	}

	return false;
}

static
void validateComplexForms(const SysComp::InCFVec *samVec, const SysComp::InCFVec *BGEVec)
{
	for (size_t idx = 0; idx < BGEVec->size(); idx++) {
		const auto &cfBGE = BGEVec->at(idx);

		/* Look for a complexForm with the same nucleus charge in
		 * sample. Fail if there is not one. */
		size_t jdx;
		for (jdx = 0; jdx < samVec->size(); jdx++) {
			const auto &cfSam = BGEVec->at(jdx);

			if (cfBGE.nucleusCharge != cfSam.nucleusCharge)
				continue;

			 /* Go through all ligand groups in BGE constituent and
			  * all ligands in them. */
			for (size_t kdx = 0; kdx < cfBGE.ligandGroups->size(); kdx++) {
				const auto &lggBGE = cfBGE.ligandGroups->at(kdx);

				for (size_t ldx = 0; ldx < lggBGE.ligands->size(); ldx++) {
					const auto &lig = lggBGE.ligands->at(ldx);

					if (!findMatchingLigandForm(lig, cfSam))
						throw InvalidComposition{InvalidComposition::Type::MISMATCHING_PARAMETERS};
				}
			}

			break;
		}

		if (jdx == samVec->size())
			throw InvalidComposition{InvalidComposition::Type::MISMATCHING_PARAMETERS};
	}
}

static
void validateCompositions(const SysComp::InConstituentVec *BGEVec, const SysComp::InConstituentVec *sampleVec, const IsAnalyteMap &iaMap)
{
	const size_t sizeSam = sampleVec->size();

	for (size_t idx = 0; idx < BGEVec->size(); idx++) {
		const auto &cBGE = BGEVec->at(idx);

		size_t jdx;
		for (jdx = 0; jdx < sizeSam; jdx++) {
			const auto &cSam = sampleVec->at(jdx);

			if (*cSam.name == *cBGE.name) {
				if (cBGE.ctype == SysComp::ConstituentType::NUCLEUS && cBGE.complexForms->size() > 0)
					validateComplexForms(cSam.complexForms, cBGE.complexForms);
				else {
					/* Sample component may have some complexations defined.
					 * This is okay iff these complexations are with analytes only. */
					if (cBGE.ctype == SysComp::ConstituentType::NUCLEUS) {
						assert(cSam.complexForms != nullptr);
						complexesOnlyWithAnalytes(cSam.complexForms, iaMap);
					}
					if (!SysComp::compareInConstituents(cBGE, cSam, false))
						throw InvalidComposition{InvalidComposition::Type::MISMATCHING_PARAMETERS};
				}
				break;
			}
		}

		if (jdx == sizeSam)
			throw InvalidComposition{InvalidComposition::Type::MISSING_IN_SAMPLE};
	}
}

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
	m_systemPackUncharged{std::move(other.m_systemPackUncharged)},
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

	Calculator::DeltaPackVec deltaPacks{};
	Calculator::DeltaPackVec deltaPacksUncharged{};
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
		ECHMET_TRACE(LEMNGTracing, EVAL_INIT_ERR, "Cannot make vectors of analytical concentrations", "Insufficient memory");

		return RetCode::E_NO_MEMORY;
	}

	try {
		applyConcentrationMapping(analConcsBGE, acBGE, m_chemicalSystemBGE);
		applyConcentrationMapping(analConcsFull, acSample, m_chemicalSystemFull);
		applyConcentrationMappingBGELike(analConcsBGELike);
	} catch (const CannotApplyConcentrationException &ex) {
		m_lastErrorString = "Cannot process input analytical concentrations";
		ECHMET_TRACE(LEMNGTracing, EVAL_INIT_ERR, "Cannot process input analytical concentrations", "Malformed input data");

		return RetCode::E_INTERNAL_ERROR;
	} catch (const ConcentrationTooLowException &ex) {
		m_lastErrorString = "Concentration of " + std::string{ex.what()} + " is too low for the numerical solver";
		ECHMET_TRACE(LEMNGTracing, EVAL_INIT_ERR, "Cannot process input analytical concentrations", "Concentration too low");

		return RetCode::E_CONCENTRATION_TOO_LOW;
	}

	auto isAnalyteFunc = [this](const std::string &s){ return this->isAnalyte(s); };

	/* Prepare output results */
	try {
		results = prepareResults(m_chemicalSystemBGE, m_chemicalSystemFull, isAnalyteFunc);
	} catch (std::bad_alloc &) {
		m_lastErrorString = "Insufficient memory to prepare results";
		ECHMET_TRACE(LEMNGTracing, EVAL_INIT_ERR, "Cannot prepare Results data structures", "Insufficient memory");

		return RetCode::E_NO_MEMORY;
	}

	Calculator::SolutionProperties BGEProps;
	try {
		BGEProps = Calculator::calculateSolutionProperties(m_chemicalSystemBGE, analConcsBGE, m_calcPropsBGE, corrections, true);
	} catch (const Calculator::CalculationException &ex) {
		releaseResults(results);
		m_lastErrorString = std::string{"Unable to calculate BGE properties: "} + ex.what();
		ECHMET_TRACE(LEMNGTracing, EVAL_PROGRESS_ERR, "Unable to calculate BGE properties", ex.what());

		return RetCode::E_CANNOT_SOLVE_BGE;
	}

	m_systemPack.conductivity = BGEProps.conductivity;

	Calculator::SolutionProperties BGELikeProps;
	/* Precalculate what is used in many places of the linear model */
	try {
		Calculator::prepareModelData(m_systemPack, m_systemPackUncharged, deltaPacks, deltaPacksUncharged, analConcsBGELike, analConcsFull, BGELikeProps, corrections);
	} catch (std::bad_alloc &) {
		fillResultsBGE(m_chemicalSystemBGE, BGEProps, corrections, results);
		return RetCode::E_NO_MEMORY;
	} catch (Calculator::CalculationException &ex) {
		fillResultsBGE(m_chemicalSystemBGE, BGEProps, corrections, results);
		m_lastErrorString = ex.what();
		ECHMET_TRACE(LEMNGTracing, EVAL_PROGRESS_ERR, "Cannot prepare model data", ex.what());

		return ex.errorCode();
	}

	/* Solve the linear model and first nonlinearity term */
	bool allZonesValid;
	try {
		Calculator::LinearResults linResults = Calculator::calculateLinear(m_systemPack, deltaPacks, corrections);
		Calculator::EigenzoneDispersionVec ezDisps = Calculator::calculateNonlinear(m_systemPack, m_systemPackUncharged, analConcsBGELike, deltaPacks, deltaPacksUncharged,
											    linResults.M1, linResults.M2, linResults.QLQR, corrections);

		fillResults(m_chemicalSystemBGE, m_chemicalSystemFull, BGEProps, BGELikeProps, linResults, ezDisps, corrections, results);
		allZonesValid = linResults.allZonesValid;
	} catch (std::bad_alloc &) {
		fillResultsBGE(m_chemicalSystemBGE, BGEProps, corrections, results);
		fillResultsAnalytesDissociation(m_chemicalSystemFull, BGELikeProps, results);
		ECHMET_TRACE(LEMNGTracing, EVAL_PROGRESS_ERR, "Cannot evaluate linear model", "Insufficient memory");

		return RetCode::E_NO_MEMORY;
	} catch (Calculator::CalculationException &ex) {
		fillResultsBGE(m_chemicalSystemBGE, BGEProps, corrections, results);
		fillResultsAnalytesDissociation(m_chemicalSystemFull, BGELikeProps, results);
		m_lastErrorString = ex.what();
		ECHMET_TRACE(LEMNGTracing, EVAL_PROGRESS_ERR, "Cannot evaluate linear model", ex.what());

		return ex.errorCode();
	}

	if (allZonesValid)
		return RetCode::OK;
	return RetCode::E_PARTIAL_EIGENZONES;
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

		validateCompositions(inCtuentVecBGE, inCtuentVecSample, iaMap);

		return new CZESystemImpl(chemSystemBGE, calcPropsBGE, chemSystemFull, calcPropsFull, std::move(iaMap));
	} catch (std::bad_alloc &up) {
		SysComp::releaseChemicalSystem(chemSystemBGE);
		SysComp::releaseChemicalSystem(chemSystemFull);
		throw up;
	}
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

	m_systemPack = Calculator::makeSystemPack(m_chemicalSystemFull, m_calcPropsFull, [this](const std::string &s) { return this->isAnalyte(s); }, false);
	m_systemPackUncharged = Calculator::makeSystemPack(m_chemicalSystemFull, m_calcPropsFull, [this](const std::string &s) { return this->isAnalyte(s); }, true);
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
		ERROR_CODE_CASE(E_PARTIAL_EIGENZONES);
		ERROR_CODE_CASE(E_INVALID_COMPOSITION_PARAMS);
		ERROR_CODE_CASE(E_INVALID_COMPOSITION_MISSING);
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
	} catch (InvalidComposition &ex) {
		switch (ex.type) {
		case InvalidComposition::Type::MISMATCHING_PARAMETERS:
			return RetCode::E_INVALID_COMPOSITION_PARAMS;
		case InvalidComposition::Type::MISSING_IN_SAMPLE:
			return RetCode::E_INVALID_COMPOSITION_MISSING;
		}
	} catch (SysCompException &ex) {
		ECHMET_TRACE(LEMNGTracing, MAKE_CZE_SYSTEM_ERR, ex.what());

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

void ECHMET_CC toggleAllTracepoints(const bool state) noexcept
{
	if (state)
		TRACER_INSTANCE<LEMNGTracing>().enableAllTracepoints();
	else
		TRACER_INSTANCE<LEMNGTracing>().disableAllTracepoints();
}

void ECHMET_CC toggleTracepoint(const int32_t TPID, const bool state) noexcept
{
	if (state)
		TRACER_INSTANCE<LEMNGTracing>().enableTracepoint(TPID);
	else
		TRACER_INSTANCE<LEMNGTracing>().disableTracepoint(TPID);
}

const char * ECHMET_CC versionString() noexcept
{
	return MK_VERSION_STRING(LEMNG_VERSION_MAJOR, LEMNG_VERSION_MINOR, LEMNG_VERSION_PATCH);
}

FixedString * ECHMET_CC trace(const bool dontClear) noexcept
{
#ifdef ECHMET_TRACER_DISABLE_TRACING
	(void)dontClear;
	return nullptr;
#else
	return createFixedString(TRACER_INSTANCE<LEMNGTracing>().logged(dontClear).c_str());
#endif // ECHMET_TRACER_DISABLE_TRACING
}

TracepointInfoVec * ECHMET_CC tracepointInfo() noexcept
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

bool ECHMET_CC tracepointState(const int32_t TPID) noexcept
{
	return TRACER_INSTANCE<LEMNGTracing>().isTracepointEnabled(TPID);
}

} // namespace LEMNG

#ifndef ECHMET_TRACER_DISABLE_TRACING

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, MAKE_CZE_SYSTEM_ERR, "makeCZEsystemError")
ECHMET_BEGIN_MAKE_LOGGER(LEMNGTracing, MAKE_CZE_SYSTEM_ERR, const char *err)
{
	return std::string{err};
}
ECHMET_END_MAKE_LOGGER

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, EVAL_INIT_ERR, "Evaluation initialization error")
ECHMET_BEGIN_MAKE_LOGGER(LEMNGTracing, EVAL_INIT_ERR, const char *where, const char *why)
{
	std::ostringstream ss{};
	ss << "Error during evaluation initialization: " << where << " (" << why << ")";

	return ss.str();
}
ECHMET_END_MAKE_LOGGER

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, EVAL_PROGRESS_ERR, "Error during evaluation")
ECHMET_BEGIN_MAKE_LOGGER(LEMNGTracing, EVAL_PROGRESS_ERR, const char *where, const char *why)
{
	std::ostringstream ss{};
	ss << "Error during evaluation: " << where << " (" << why << ")";

	return ss.str();
}
ECHMET_END_MAKE_LOGGER

#endif // ECHMET_TRACER_DISABLE_TRACING

} // namespace ECHMET
