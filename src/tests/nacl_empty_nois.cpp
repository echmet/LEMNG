#include <cstdlib>
#include "barsarkagang_tests.h"


using namespace ECHMET;
using namespace ECHMET::Barsarkagang;


int main(int , char ** )
{
	SysComp::InConstituent chloride{
		SysComp::ConstituentType::NUCLEUS,
		createFixedString("Chloride"),
		-1,
		0,
		mkRealVec( { -2.0 } ),
		mkRealVec( { 79.1, 0.0 } ),
		noComplexes(),
		0.0
	};

	SysComp::InConstituent sodium{
		SysComp::ConstituentType::NUCLEUS,
		createFixedString("Sodium"),
		0,
		1,
		mkRealVec( { 13.7 } ),
		mkRealVec( { 0.0, 51.9 } ),
		noComplexes(),
		0.0
	};

	CMapping cBGE = {
		{ "Chloride", 9.0 },
		{ "Sodium", 10.0 }
	};

	CMapping cSample = {
		{ "Chloride", 7.0 },
		{ "Sodium", 8.0 }
	};

	auto r = calculate(
		{
			chloride,
			sodium
		},
		{
			chloride,
			sodium
		},
		cBGE, cSample,
		false, false, false, false);

	checkBGE(r, 10.991436593, 0.13810663565, 0.0099804751555, 2.3024973608);

	checkEigenzone(1, r.eigenzones, 1.808142924e-07, 1.1114445903e-07, 1.4684061008, 10.891527377, 0.10981822032);

	checkEigenzone(2, r.eigenzones, -179.83608821, 8.0629168727, 1.4865248413, 11.073709041, 0.14119055221);

	LEMNG::releaseResults(r);
	SysComp::releaseInConstituent(chloride);
	SysComp::releaseInConstituent(sodium);

	return EXIT_SUCCESS;
}

