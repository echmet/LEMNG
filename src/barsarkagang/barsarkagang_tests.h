#ifndef BARSARKAGANG_TESTS_H
#define BARSARKAGANG_TESTS_H

#include <cstdlib>
#include <initializer_list>
#include <string>
#include <vector>

#include <lemng.h>

namespace ECHMET {
namespace Barsarkagang {

using CMapping = std::vector<std::pair<std::string, double>>;
using InConstituentList = std::initializer_list<SysComp::InConstituent>;

static const double TOLERANCE{1.0e-8};

static inline
bool numberMatches(const double got, const double expected, const double tol = TOLERANCE)
{
	const double TMAX = 1.0 + 0.5 * tol;
	const double TMIN = 1.0 - 0.5 * tol;
	const double norm_got = got / expected;

	fprintf(stderr, "%0.9g;%0.9g\n", got, norm_got);

	return norm_got <= TMAX && norm_got >= TMIN;
}

static inline
void failIfFalse(const bool b)
{
	if (!b)
		std::exit(EXIT_FAILURE);
}

static inline
void failIfMismatch(const double got, const double expected, const double tol = TOLERANCE)
{

	if (!numberMatches(got, expected, tol)) {
		fprintf(stderr, "Value mismatch: got %0.9g; expected %0.9g\n", got, expected);
		std::exit(EXIT_FAILURE);
	}
}

static inline
void failIfError(const ::ECHMET::RetCode tRet)
{
	if (tRet != ::ECHMET::RetCode::OK) {
		std::cerr << errorToString(tRet);
		std::exit(EXIT_FAILURE);
	}
}

static inline
void failIfError(const LEMNG::RetCode tRet)
{
	if (tRet != LEMNG::RetCode::OK) {
		std::cerr << LEMNG::LEMNGerrorToString(tRet) << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

static inline
void checkEigenzone(const LEMNG::REigenzoneVec *ezs, const double u, const double uEMD, const double pH)
{
	for (size_t idx = 0; idx < ezs->size(); idx++) {
		auto &ez = ezs->at(idx);

		if (numberMatches(ez.mobility, u, 0.1)) {
			failIfMismatch(ez.uEMD, uEMD);
			failIfMismatch(ez.solutionProperties.pH, pH);

			return;
		}
	}

	std::cerr << "Eigenzone not found" << std::endl;

	std::exit(EXIT_FAILURE);
}

static inline
void checkSolProps(const LEMNG::RSolutionProperties &props, const double pH, const double cond, const double ionicStrength)
{
	failIfMismatch(props.pH, pH);
	failIfMismatch(props.conductivity, cond);
	failIfMismatch(props.ionicStrength, ionicStrength);
}

static inline
SysComp::InConstituentVec * mkInConstVec(std::initializer_list<SysComp::InConstituent> l)
{
	auto v = SysComp::createInConstituentVec(0);

	for (auto &&i : l)
		v->push_back(i);

	return v;
}

static inline
LEMNG::Results calculate(const InConstituentList &bge, const InConstituentList &sample, const CMapping &bgeMaps, const CMapping &sampleMaps)
{
	LEMNG::CZESystem *czeSys;

	auto ret = LEMNG::makeCZESystem(mkInConstVec(bge), mkInConstVec(sample), czeSys);
	failIfError(ret);

	LEMNG::InAnalyticalConcentrationsMap *acBGEMap = nullptr;
	LEMNG::InAnalyticalConcentrationsMap *acSampleMap = nullptr;

	failIfError(czeSys->makeAnalyticalConcentrationsMaps(acBGEMap, acSampleMap));

	for (auto &&item : bgeMaps)
		acBGEMap->item(item.first.c_str()) = item.second;

	for (auto &&item : sampleMaps)
		acSampleMap->item(item.first.c_str()) = item.second;

	auto corrections = defaultNonidealityCorrections();

	LEMNG::Results results;
	failIfError(czeSys->evaluate(acBGEMap, acSampleMap, corrections, results));

	return results;
}

static inline
RealVec * mkRealVec(std::initializer_list<double> l)
{
	auto v = createRealVec(0);

	for (auto &&i : l)
		v->push_back(i);

	return v;
}

static inline
SysComp::InCFVec * noComplexes()
{
	return SysComp::createInCFVec(0);
}

} // namespace Barsarkagang
} // namespace ECHMET

#endif // BARSARKAGANG_TESTS_H
