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

	const auto r = calculate(
		{
			chloride,
			sodium
		},
		{
			chloride,
			sodium
		},
		cBGE, cSample,
		true, true, false);

	checkBGE(r, 10.949715048, 0.1300633734, 0.0099839393407, 2.302525756);

	checkEigenzone(r.eigenzones, 2.048216106e-07, 1.258946201e-07, 10.85379174, 0.10403577448);

	checkEigenzone(r.eigenzones, -172.04923579, 7.8420114551, 11.031664059, 0.13302316692);

	return EXIT_SUCCESS;
}

