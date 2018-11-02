#include <cstdlib>
#include <iostream>

#include "barsarkagang_tests.h"

using namespace ECHMET;
using namespace ECHMET::Barsarkagang;

int main(int, char **)
{

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

	CMapping cBGE = { { "Sodium", 10.0 },
			  { "Chloride", 9.0 } };
	CMapping cSample = { { "Sodium", 8.0 },
			     { "Chloride", 7.0 } };

	auto r = calculate({ sodium, chloride }, { sodium, chloride },
			   cBGE, cSample);

	failIfFalse(r.isBGEValid);

	checkSolProps(r.BGEProperties, 10.9914366, 0.138106636, 0.00998047516);
	failIfMismatch(r.BGEProperties.bufferCapacity, 2.30249736);

	checkEigenzone(r.eigenzones, 1.8079625298577537e-7, 1.111370655142905e-7, 10.891527376537205);



	return EXIT_SUCCESS;
}
