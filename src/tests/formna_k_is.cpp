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

	auto r = calculate(
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
		true, true, false, false);

	checkBGE(r, 3.6842351542, 0.087109803668, 0.008226777682, 10.298050561);

	checkEigenzone(1, r.eigenzones, -2.3647765812e-07, -3.4854650466e-07, 1.3407500778, 4.2087522655, 0.049043151237);

	checkEigenzone(2, r.eigenzones, 16.976779795, -10.110721518, 1.2938794741, 3.972135427, 0.10077837249);

	checkEigenzone(3, r.eigenzones, 72.0217234, -13.633815392, 1.8504232006, 3.713550978, 0.09532464026);

	LEMNG::releaseResults(r);
	SysComp::releaseInConstituent(formic_acid);
	SysComp::releaseInConstituent(na);
	SysComp::releaseInConstituent(k);

	return EXIT_SUCCESS;
}

