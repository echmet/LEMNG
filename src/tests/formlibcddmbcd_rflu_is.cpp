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
				}, /* InLigandGroup c-tor end */
				{ /* InLigandGroup c-tor begin */
					/* InLFVec */
					{
						{ /* InLigandForm c-tor begin */
							"DM-b-CD",
							0,
							1,
							{ -3.847572659142112 },
							{ 7.47 }
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
				}, /* InLigandGroup c-tor end */
				{ /* InLigandGroup c-tor begin */
					/* InLFVec */
					{
						{ /* InLigandForm c-tor begin */
							"DM-b-CD",
							0,
							1,
							{ -4.0 },
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

	SysComp::InConstituent dm__b__cd{
		SysComp::ConstituentType::LIGAND,
		createFixedString("DM-b-CD"),
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
		mkRealVec( { 4.19 } ),
		mkRealVec( { 20.4, 0.0 } ),
		gen_complexforms_r__flu(),
		0.0
	};

	CMapping cBGE = {
		{ "Formic acid", 10.0 },
		{ "Li", 5.0 },
		{ "b-CD", 5.0 },
		{ "DM-b-CD", 5.0 }
	};

	CMapping cSample = {
		{ "Formic acid", 5.0 },
		{ "Li", 2.5 },
		{ "b-CD", 5.0 },
		{ "DM-b-CD", 5.0 },
		{ "R-FLU", 0.1 }
	};

	auto r = calculate(
		{
			formic_acid,
			li,
			b__cd,
			dm__b__cd
		},
		{
			formic_acid,
			li,
			b__cd,
			dm__b__cd,
			r__flu
		},
		cBGE, cSample,
		true, true, false, false);

	checkBGE(r, 3.7525300653, 0.051399578407, 0.0051905732156, 6.1869113609);

	checkEigenzone(1, r.eigenzones, -5.6304324029e-16, -1.0382233868e-17, 0.51385140826, 3.7525300653, 0.051399578406);

	checkEigenzone(2, r.eigenzones, -1.5779375901e-14, -1.309953715e-18, 0.51385140826, 3.7525300653, 0.051399578407);

	checkEigenzone(3, r.eigenzones, -4.3673515708e-07, -5.3682658999e-07, 1.2680782754, 4.0871952346, 0.032149385576);

	checkEigenzone(4, r.eigenzones, -1.6425123101, -0.15500431692, 0.21638148291, 3.7634763458, 0.049966974598);

	checkEigenzone(5, r.eigenzones, 21.206795057, 5.1917446861, 1.091388496, 3.6122141438, 0.04907888601);

	LEMNG::releaseResults(r);
	SysComp::releaseInConstituent(formic_acid);
	SysComp::releaseInConstituent(li);
	SysComp::releaseInConstituent(b__cd);
	SysComp::releaseInConstituent(dm__b__cd);
	SysComp::releaseInConstituent(r__flu);

	return EXIT_SUCCESS;
}

