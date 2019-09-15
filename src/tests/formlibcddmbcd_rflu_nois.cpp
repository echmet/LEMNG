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

	const auto r = calculate(
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
		false, false, false, false);

	checkBGE(r, 3.7807882258, 0.053349365977, 0.0051656577556, 6.1315850881);

	checkEigenzone(r.eigenzones, -6.2261000774e-16, -1.1431151474e-17, 0.51385140826, 3.7807882258, 0.053349365975);

	checkEigenzone(r.eigenzones, -1.5653481447e-14, -1.3445010186e-18, 0.51385140826, 3.7807882258, 0.053349365976);

	checkEigenzone(r.eigenzones, -4.5657340426e-07, -5.5788050222e-07, 1.3420079091, 4.1069236186, 0.033423472052);

	checkEigenzone(r.eigenzones, -1.6362978402, -0.16185703553, 0.21701802171, 3.7917921482, 0.051754096246);

	checkEigenzone(r.eigenzones, 19.791803183, 5.3115949695, 1.1569824945, 3.6378241546, 0.050359139765);

	return EXIT_SUCCESS;
}

