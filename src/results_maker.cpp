#include "results_maker.h"

#ifndef ECHMET_IMPORT_INTERNAL
#define ECHMET_IMPORT_INTERNAL
#endif // ECHMET_IMPORT_INTERNAL
#include <echmetionprops.h>

#define USE_ECHMET_CONTAINERS
#include <containers/echmetvec_p.h>
#include <containers/echmetskmap_p.h>

namespace ECHMET {
namespace LEMNG {

void releaseRIon(RIon &ion) noexcept
{
	ion.name->destroy();
}

void releaseRForm(RForm &f) noexcept
{
	for (size_t idx = 0; idx < f.ions->size(); idx++)
		releaseRIon((*f.ions)[idx]);

	f.ions->destroy();
}

void releaseRConstituent(RConstituent &c)
{
	SKMapImpl<RForm> *formsImpl = dynamic_cast<SKMapImpl<RForm> *>(c.forms);
	if (formsImpl == nullptr)
		throw std::runtime_error{"Failed cast to SKMapImpl<RForm> *"};

	for (auto &&f : formsImpl->STL())
		releaseRForm(f.second);

	c.name->destroy();
	c.forms->destroy();
}

void teardownRConstituentMap(SKMapImpl<RConstituent> *map)
{
	for (auto &&c : map->STL())
		releaseRConstituent(c.second);

	map->destroy();
}

void releaseRSolutionProperties(RSolutionProperties &props)
{
	SKMapImpl<RConstituent> *compositionImpl = dynamic_cast<SKMapImpl<RConstituent> *>(props.composition);
	if (compositionImpl == nullptr)
		throw std::runtime_error{"Failed cast to SKMapImpl<RConstituent> *"};

	teardownRConstituentMap(compositionImpl);
}

void teardownRFormMap(SKMapImpl<RForm> *map)
{
	for (auto &&f : map->STL())
		releaseRForm(f.second);

	map->destroy();
}

void teardownRIonVec(VecImpl<RIon, false> *vec)
{
	for (auto &&i : vec->STL())
		releaseRIon(i);

	vec->destroy();
}

void teardownREigenzoneVec(VecImpl<REigenzone, false> *vec)
{
	for (auto &&ez: vec->STL())
		releaseRSolutionProperties(ez.solutionProperties);

	vec->destroy();
}

typedef std::unique_ptr<SKMapImpl<RConstituent>, decltype(&teardownRConstituentMap)> RConstituentMapWrapper;
typedef std::unique_ptr<VecImpl<REigenzone, false>, decltype(&teardownREigenzoneVec)> REigenzoneVecWrapper;
typedef std::unique_ptr<SKMapImpl<RForm>, decltype(&teardownRFormMap)> RFormMapWrapper;
typedef std::unique_ptr<VecImpl<RIon, false>, decltype(&teardownRIonVec)> RIonVecWrapper;

template <typename T>
void zeroize(typename std::enable_if<!std::is_pointer<T>::value, T>::type *t)
{
	const size_t SZ = sizeof(T);
	memset(t, 0, SZ);
}

RConstituentMapWrapper prepareComposition(const ChemicalSystemPtr &chemSystem)
{
	RConstituentMapWrapper compMap{new SKMapImpl<RConstituent>{}, teardownRConstituentMap};

	for (size_t idx = 0; idx < chemSystem->constituents->size(); idx++) {
		const SysComp::Constituent *c = chemSystem->constituents->at(idx);
		RFormMapWrapper fMap{new SKMapImpl<RForm>, teardownRFormMap};

		for (size_t jdx = 0; jdx < c->ionicForms->size(); jdx++) {
			RIonVecWrapper ions{new VecImpl<RIon, false>{}, teardownRIonVec};

			const SysComp::IonicForm *iF = c->ionicForms->at(jdx);

			/* Add the nucleus */
			RIon ion{createFixedString(iF->nucleus->name->c_str()),
				 iF->totalCharge,
				 1};

			if (ion.name == nullptr)
				throw std::bad_alloc{};

			if (ions->push_back(ion) != ::ECHMET::RetCode::OK)
				throw std::bad_alloc{};

			if (iF->ligand != nullptr) {
				const SysComp::IonicForm *currentIF = iF;
				const FixedString *lastName = iF->ligand->name;
				int32_t lastCharge = iF->ligandCharge;

				RIon ion{createFixedString(iF->ligand->name->c_str()),
					 iF->ligandCharge,
					 iF->ligandCount};

				if (ion.name == nullptr)
					throw std::bad_alloc{};

				if (ions->push_back(ion) != ::ECHMET::RetCode::OK)
					throw std::bad_alloc{};

				while (currentIF->ancestor->ligand != nullptr) {
					currentIF = currentIF->ancestor;

					if (*(currentIF->ligand->name) != *lastName || currentIF->ligandCharge != lastCharge) {
						RIon _ion{createFixedString(currentIF->ligand->name->c_str()),
							  currentIF->ligandCharge,
							  currentIF->ligandCount};

						if (ion.name == nullptr)
							throw std::bad_alloc{};

						if (ions->push_back(_ion) != ::ECHMET::RetCode::OK)
							throw std::bad_alloc{};
					}

					lastName = currentIF->ligand->name;
					lastCharge = currentIF->ligandCharge;
				}
			}

			RForm form{iF->totalCharge, 0, ions.release()};

			try {
				fMap->STL().emplace(iF->name->c_str(), form);
			} catch (std::bad_alloc &) {
				releaseRForm(form);
				throw;
			}
		}

		SKMapImpl<RForm> *fMapRaw = fMap.release();
		RConstituent ctuent{createFixedString(c->name->c_str()), 0, fMapRaw};
		if (ctuent.name == nullptr) {
			teardownRFormMap(fMapRaw);
			throw;
		}

		try {
			compMap->STL().emplace(c->name->c_str(), ctuent);
		} catch (std::bad_alloc &) {
			releaseRConstituent(ctuent);
			throw;
		}
	}

	return compMap;
}

REigenzoneVecWrapper prepareEigenzones(const ChemicalSystemPtr &chemSystem)
{
	REigenzoneVecWrapper eigenzones{new VecImpl<REigenzone, false>{}, teardownREigenzoneVec};

	for (size_t idx = 0; idx < chemSystem->constituents->size(); idx++) {
		REigenzone ez;

		zeroize<REigenzone>(&ez);

		auto composition = prepareComposition(chemSystem);
		SKMapImpl<RConstituent> *compositionRaw = composition.release();

		ez.solutionProperties.composition = compositionRaw;
		if (eigenzones->push_back(ez) != ::ECHMET::RetCode::OK) {
			teardownRConstituentMap(compositionRaw);
			throw;
		}
	}

	return eigenzones;
}

void fillSolutionProperties(const ChemicalSystemPtr &chemSystem, const Calculator::SolutionProperties &props, const NonidealityCorrections corrections, RSolutionProperties &rProps)
{
	auto H3OConcentration = [](const std::vector<double> &icVec) {
		return icVec.at(0);
	};

	auto mapComposition =  [](const ChemicalSystemPtr &chemSystem, const std::vector<double> &anConcs, const std::vector<double> &icConcs, RConstituentMap *composition) {
		auto &compositionSTL = static_cast<SKMapImpl<RConstituent> *>(composition)->STL();

		for (size_t idx = 0; idx < chemSystem->constituents->size(); idx++) {
			const SysComp::Constituent *c = chemSystem->constituents->at(idx);
			RConstituent &rCtuent = compositionSTL[c->name->c_str()];
			auto &rFormsSTL = static_cast<SKMapImpl<RForm> *>(rCtuent.forms)->STL();

			rCtuent.concentration = anConcs.at(c->analyticalConcentrationIndex);

			for (size_t jdx = 0; jdx < c->ionicForms->size(); jdx++) {
				const SysComp::IonicForm *iF = c->ionicForms->at(jdx);

				RForm &rForm = rFormsSTL[iF->name->c_str()];
				rForm.concentration = icConcs.at(iF->ionicConcentrationIndex);
			}
		}
	};

	const bool correctForIS = corrections & NonidealityCorrections::CORR_DEBYE_HUCKEL;
	const ECHMETReal cH = H3OConcentration(props.ionicConcentrations);
	rProps.bufferCapacity = props.bufferCapacity;
	rProps.conductivity = props.conductivity;
	rProps.ionicStrength = props.ionicStrength;
	rProps.pH = IonProps::calculatepH_direct(cH, (correctForIS ? props.ionicStrength : 0.0));
	mapComposition(chemSystem, props.analyticalConcentrations, props.ionicConcentrations, rProps.composition);
}

void fillResults(const ChemicalSystemPtr &chemSystemBGE, const ChemicalSystemPtr &chemSystemFull, const Calculator::SolutionProperties &BGEProperties, const Calculator::LinearResults &linResults, const Calculator::EigenzoneDispersionVec &ezDisps, const NonidealityCorrections corrections, Results &r)
{
	auto fillEigenzone = [corrections](const ChemicalSystemPtr &chemSystem, const Calculator::Eigenzone &ez, const Calculator::EigenzoneDispersion &disp, REigenzone &rEz) {
		fillSolutionProperties(chemSystem, ez.solutionProperties, corrections, rEz.solutionProperties);

		rEz.mobility = ez.zoneMobility;
		rEz.a2t = disp.a2t;
		rEz.uEMD = disp.uEMD;
		rEz.tainted = ez.tainted;
		rEz.ztype = [ez]() {
			if (ez.isAnalyzeZone)
				return EigenzoneType::ANALYTE;
			return EigenzoneType::SYSTEM;
		}();
	};

	/* Fill out BGE properties */
	fillSolutionProperties(chemSystemBGE, BGEProperties, corrections, r.BGEProperties);
	r.isBGEValid = true;

	/* Fill out all eigenzones */
	for (size_t idx = 0; idx < linResults.eigenzones.size(); idx++) {
		const Calculator::Eigenzone &ez = linResults.eigenzones.at(idx);
		const Calculator::EigenzoneDispersion &disp = ezDisps.at(idx);
		REigenzone &rEz = (*r.eigenzones)[idx];

		fillEigenzone(chemSystemFull, ez, disp, rEz);
	}
}

void fillResultsPartial(const ChemicalSystemPtr &chemSystemBGE, const Calculator::SolutionProperties &BGEProperties, const NonidealityCorrections corrections, Results &r)
{
	fillSolutionProperties(chemSystemBGE, BGEProperties, corrections, r.BGEProperties);
	r.isBGEValid = true;
}

Results prepareResults(const ChemicalSystemPtr &chemSystemBGE, const ChemicalSystemPtr &chemSystemFull)
{
	Results r;
	zeroize<Results>(&r);

	RConstituentMapWrapper BGEPropertiesWrapped = prepareComposition(chemSystemBGE);
	REigenzoneVecWrapper eigenzones = prepareEigenzones(chemSystemFull);

	r.BGEProperties.composition = BGEPropertiesWrapped.release();
	r.eigenzones = eigenzones.release();

	return r;
}

void ECHMET_CC releaseResults(Results &r) noexcept
{
	releaseRSolutionProperties(r.BGEProperties);

	VecImpl<REigenzone, false> *eigenzonesImpl = dynamic_cast<VecImpl<REigenzone, false> *>(r.eigenzones);
	if (eigenzonesImpl == nullptr)
		return;

	teardownREigenzoneVec(eigenzonesImpl);
}

} // namespace LEMNG
} // namespace ECHMET
