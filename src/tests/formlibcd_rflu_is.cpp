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

SysComp::InCFVec * gen_complexforms_r__flu()
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
							"b-CD",
							0,
							1,
							{ -3.699837725867246 },
							{ 9.24 }
						} /* InLigandForm c-tor end */
					}
				} /* InLigandGroup c-tor end */
			}
		}, /* InComplexForm c-tor end */
		{ /* InComplexForm c-tor begin */
			0,
			/* InLGVec */
			{
				{ /* InLigandGroup c-tor begin */
					/* InLFVec */
					{
						{ /* InLigandForm c-tor begin */
							"b-CD",
							0,
							1,
							{ -3.980003371583746 },
							{ 0.0 }
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

	SysComp::InConstituent b__cd{
		SysComp::ConstituentType::LIGAND,
		createFixedString("b-CD"),
		0,
		0,
		mkRealVec( {  } ),
		mkRealVec( { 0.0 } ),
		nullptr,
		0.0
	};

	SysComp::InConstituent r__flu{
		SysComp::ConstituentType::NUCLEUS,
		createFixedString("R-FLU"),
		-1,
		0,
		mkRealVec( { 4.1 } ),
		mkRealVec( { 20.4, 0.0 } ),
		gen_complexforms_r__flu(),
		0.0
	};

	CMapping cBGE = {
		{ "Formic acid", 10.0 },
		{ "Li", 5.0 },
		{ "b-CD", 10.0 }
	};

	CMapping cSample = {
		{ "Formic acid", 5.0 },
		{ "Li", 2.5 },
		{ "b-CD", 10.0 },
		{ "R-FLU", 0.1 }
	};

	const auto r = calculate(
		{
			formic_acid,
			li,
			b__cd
		},
		{
			formic_acid,
			li,
			b__cd,
			r__flu
		},
		cBGE, cSample,
		true, true, false, false);

	checkBGE(r, 3.7525300653, 0.051399578409, 0.0051905732156, 6.1869113609);

	checkEigenzone(r.eigenzones, -6.715753817e-16, -1.2987600378e-17, 0.51385140826, 3.7525300653, 0.051399578406);

	checkEigenzone(r.eigenzones, -4.3684704003e-07, -5.4153304259e-07, 1.2680782754, 4.0937519747, 0.031966230541);

	checkEigenzone(r.eigenzones, -1.9218729088, -0.15945839744, 0.24022475218, 3.762025126, 0.050124563128);

	checkEigenzone(r.eigenzones, 21.206795057, 5.1960316583, 1.091388496, 3.6121052418, 0.049077189674);

	return EXIT_SUCCESS;
}

