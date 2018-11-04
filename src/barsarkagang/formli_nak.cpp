#include <cstdlib>
#include <iostream>

#include "barsarkagang_tests.h"

using namespace ECHMET;
using namespace ECHMET::Barsarkagang;

int main(int, char **)
{
	/* BGE */
	SysComp::InConstituent formic_acid{
		SysComp::ConstituentType::NUCLEUS,
		createFixedString("Formic acid"),
		-1,
		0,
		mkRealVec( { 3.752 } ),
		mkRealVec( { 56.6, 0.0 } ),
		noComplexes(),
		0.0
	};

	SysComp::InConstituent lithium{
		SysComp::ConstituentType::NUCLEUS,
		createFixedString("Lithium"),
		0,
		1,
		mkRealVec( { 13.8 } ),
		mkRealVec( { 0.0, 40.1 } ),
		noComplexes(),
		0.0
	};

	/* Analytes */
	SysComp::InConstituent potassium{
		SysComp::ConstituentType::NUCLEUS,
		createFixedString("Potassium"),
		0,
		1,
		mkRealVec( { 13 } ),
		mkRealVec( { 0.0, 76.2 } ),
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

	CMapping cBGE = { { "Formic acid", 30.0 },
		          { "Lithium", 15.0 } };

	CMapping cSample = { { "Formic acid", 25.0 },
			     { "Lithium", 13.0 },
		             { "Potassium", 0.1 },
		             { "Sodium", 0.2 } };

	auto r = calculate({ formic_acid, lithium },
			   { formic_acid, lithium, potassium, sodium },
			   cBGE, cSample);

	checkBGE(r, 3.7620167177, 0.14694657546, 0.015172974976, 17.665381469);

	checkEigenzone(r.eigenzones, -1.449570266e-07, -5.7548838091e-08, 3.8215787468, 0.13011347082);
	checkEigenzone(r.eigenzones, 7.6900661281, 0.080672507344, 3.7574545028, 0.14653407349);
	checkEigenzone(r.eigenzones, 51.899999994, -0.34631776329, 3.7630511705, 0.14743687887);
	checkEigenzone(r.eigenzones, 76.199999956, -0.63761663008, 3.7630161769, 0.14756140682);

	return EXIT_SUCCESS;
}
