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

SysComp::InCFVec * gen_complexforms_li()
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

SysComp::InCFVec * gen_complexforms_s_BGE()
{
	const ComplexDef cDef = {
		{ /* InComplexForm c-tor begin */
			0,
			/* InLGVec */
			{
			}
		} /* InComplexForm c-tor end */
	};

	return buildComplexes(cDef);
}

SysComp::InCFVec * gen_complexforms_s()
{
	const ComplexDef cDef = {
		{ /* InComplexForm c-tor begin */
			0,
			/* InLGVec */
			{
				{ /* InLigandGroup c-tor begin */
					/* InLFVec */
					{
						{ /* InLigandForm c-tor begin */
							"X",
							-1,
							2,
							{ -3.778151250383644, -3.477121254719662 },
							{ 10.0, 5.0 }
						} /* InLigandForm c-tor end */
					}
				} /* InLigandGroup c-tor end */
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

	SysComp::InConstituent li{
		SysComp::ConstituentType::NUCLEUS,
		createFixedString("Li"),
		0,
		1,
		mkRealVec( { 13.8 } ),
		mkRealVec( { 0.0, 40.1 } ),
		gen_complexforms_li(),
		0.0
	};

	SysComp::InConstituent x{
		SysComp::ConstituentType::LIGAND,
		createFixedString("X"),
		-1,
		-1,
		mkRealVec( {  } ),
		mkRealVec( { 20.0 } ),
		nullptr,
		0.0
	};

	/* This one was added manually - automated filtering
	 * of ligand analytes would be too much of a hassle */
	SysComp::InConstituent s_BGE{
		SysComp::ConstituentType::NUCLEUS,
		createFixedString("S"),
		0,
		0,
		mkRealVec( {  } ),
		mkRealVec( { 0.0 } ),
		gen_complexforms_s_BGE(),
		0.0
	};

	SysComp::InConstituent s{
		SysComp::ConstituentType::NUCLEUS,
		createFixedString("S"),
		0,
		0,
		mkRealVec( {  } ),
		mkRealVec( { 0.0 } ),
		gen_complexforms_s(),
		0.0
	};

	CMapping cBGE = {
		{ "Formic acid", 10.0 },
		{ "Li", 5.0 },
		{ "S", 10.0 }
	};

	CMapping cSample = {
		{ "Formic acid", 5.0 },
		{ "Li", 2.5 },
		{ "X", 0.2 },
		{ "S", 10.0 }
	};

	auto r = calculate(
		{
			formic_acid,
			li,
			s_BGE
		},
		{
			formic_acid,
			li,
			x,
			s
		},
		cBGE, cSample,
		true, true, false, false);

	checkBGE(r, 3.7525300653, 0.051399578409, 0.0051905732156, 6.1869113609);

	checkEigenzone(r.eigenzones, -1.6124697662e-15, -1.0844943009e-16, 0.51385140826, 3.7525300653, 0.051399578407);

	checkEigenzone(r.eigenzones, -4.3684723305e-07, -5.0886101286e-07, 1.2680782754, 4.050382941, 0.033254934945);

	checkEigenzone(r.eigenzones, -10.128331975, -0.74266661947, 0.26022288244, 3.7512294009, 0.048914937445);

	checkEigenzone(r.eigenzones, 21.206795057, 5.5081075415, 1.091388496, 3.6042029176, 0.048954771308);

	LEMNG::releaseResults(r);
	SysComp::releaseInConstituent(formic_acid);
	SysComp::releaseInConstituent(li);
	SysComp::releaseInConstituent(x);
	SysComp::releaseInConstituent(s_BGE);
	SysComp::releaseInConstituent(s);

	return EXIT_SUCCESS;
}

