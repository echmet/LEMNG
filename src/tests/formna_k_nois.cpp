#include <cstdlib>
#include "barsarkagang_tests.h"


using namespace ECHMET;
using namespace ECHMET::Barsarkagang;


SysComp::InCFVec * gen_complexforms_formic_acid()
{
	const ComplexDef cDef = {
		{ /* InComplexForm c-tor begin */
			-1,
			/* InLGVec */
			{
			}
		}, /* InComplexForm c-tor end */
		{ /* InComplexForm c-tor begin */
			0,
			/* InLGVec */
			{
			}
		} /* InComplexForm c-tor end */
	};

	return buildComplexes(cDef);
}

SysComp::InCFVec * gen_complexforms_na()
{
	const ComplexDef cDef = {
		{ /* InComplexForm c-tor begin */
			0,
			/* InLGVec */
			{
			}
		}, /* InComplexForm c-tor end */
		{ /* InComplexForm c-tor begin */
			1,
			/* InLGVec */
			{
			}
		} /* InComplexForm c-tor end */
	};

	return buildComplexes(cDef);
}

SysComp::InCFVec * gen_complexforms_k()
{
	const ComplexDef cDef = {
		{ /* InComplexForm c-tor begin */
			0,
			/* InLGVec */
			{
			}
		}, /* InComplexForm c-tor end */
		{ /* InComplexForm c-tor begin */
			1,
			/* InLGVec */
			{
			}
		} /* InComplexForm c-tor end */
	};

	return buildComplexes(cDef);
}

int main(int , char ** )
{
	SysComp::InConstituent formic_acid{
		SysComp::ConstituentType::NUCLEUS,
		createFixedString("Formic acid"),
		-1,
		0,
		mkRealVec( { 3.752 } ),
		mkRealVec( { 56.6, 0.0 } ),
		gen_complexforms_formic_acid(),
		0.0
	};

	SysComp::InConstituent na{
		SysComp::ConstituentType::NUCLEUS,
		createFixedString("Na"),
		0,
		1,
		mkRealVec( { 13.7 } ),
		mkRealVec( { 0.0, 51.9 } ),
		gen_complexforms_na(),
		0.0
	};

	SysComp::InConstituent k{
		SysComp::ConstituentType::NUCLEUS,
		createFixedString("K"),
		0,
		1,
		mkRealVec( { 13.0 } ),
		mkRealVec( { 0.0, 76.2 } ),
		gen_complexforms_k(),
		0.0
	};

	CMapping cBGE = {
		{ "Formic acid", 17.0 },
		{ "Na", 8.0 }
	};

	CMapping cSample = {
		{ "Formic acid", 5.0 },
		{ "Na", 5.0 },
		{ "K", 2.0 }
	};

	const auto r = calculate(
		{
			formic_acid,
			na
		},
		{
			formic_acid,
			na,
			k
		},
		cBGE, cSample,
		false, false, false, false);

	checkBGE(r, 3.7203482361, 0.091448211939, 0.0081903933441, 10.211400274);

	checkEigenzone(r.eigenzones, -2.4976645929e-07, -3.6792715438e-07, 4.2389183872, 0.051285698972);

	checkEigenzone(r.eigenzones, 15.320782056, -9.7671951499, 4.0164006602, 0.10766747473);

	checkEigenzone(r.eigenzones, 76.19999996, -14.172314094, 3.7474959277, 0.099974144085);

	return EXIT_SUCCESS;
}

