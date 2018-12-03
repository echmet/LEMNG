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

SysComp::InCFVec * gen_complexforms_x()
{
	const ComplexDef cDef = {
		{ /* InComplexForm c-tor begin */
			-1,
			/* InLGVec */
			{
				{ /* InLigandGroup c-tor begin */
					/* InLFVec */
					{
						{ /* InLigandForm c-tor begin */
							"S",
							0,
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
		SysComp::ConstituentType::NUCLEUS,
		createFixedString("X"),
		-1,
		-1,
		mkRealVec( {  } ),
		mkRealVec( { 20.0 } ),
		gen_complexforms_x(),
		0.0
	};

	SysComp::InConstituent s{
		SysComp::ConstituentType::LIGAND,
		createFixedString("S"),
		0,
		0,
		mkRealVec( {  } ),
		mkRealVec( { 0.0 } ),
		nullptr,
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

	const auto r = calculate(
		{
			formic_acid,
			li,
			s
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

	checkEigenzone(r.eigenzones, -3.8691823131e-15, -2.5786171514e-16, 0.51385140826, 3.7525300656, 0.051399578379);

	checkEigenzone(r.eigenzones, -4.3961669061e-07, -4.4507499949e-07, 1.2680782754, 3.9803824687, 0.035786599538);

	checkEigenzone(r.eigenzones, -5.1680968568, -1.4356099061, 0.1327816924, 3.7713378643, 0.046563498835);

	checkEigenzone(r.eigenzones, 21.20679506, 5.4496673708, 1.091388496, 3.6056789558, 0.048977534536);

	return EXIT_SUCCESS;
}

