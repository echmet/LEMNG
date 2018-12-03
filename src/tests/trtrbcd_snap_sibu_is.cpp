#include <cstdlib>
#include "barsarkagang_tests.h"


using namespace ECHMET;
using namespace ECHMET::Barsarkagang;


SysComp::InCFVec * gen_complexforms_tricine__ani()
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

SysComp::InCFVec * gen_complexforms_tris()
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

SysComp::InCFVec * gen_complexforms_s__nap()
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
							{ -2.704150516839799 },
							{ 9.12 }
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
							{ -3.436162647040756 },
							{ 0.0 }
						} /* InLigandForm c-tor end */
					}
				} /* InLigandGroup c-tor end */
			}
		} /* InComplexForm c-tor end */
	};

	return buildComplexes(cDef);
}

SysComp::InCFVec * gen_complexforms_s__ibu()
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
							{ -3.818225893613955 },
							{ 9.48 }
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
							{ -4.264345507050092 },
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
	SysComp::InConstituent tricine__ani{
		SysComp::ConstituentType::NUCLEUS,
		createFixedString("Tricine-ANI"),
		-1,
		0,
		mkRealVec( { 8.15 } ),
		mkRealVec( { 30.0, 0.0 } ),
		gen_complexforms_tricine__ani(),
		0.0
	};

	SysComp::InConstituent tris{
		SysComp::ConstituentType::NUCLEUS,
		createFixedString("TRIS"),
		0,
		1,
		mkRealVec( { 8.076 } ),
		mkRealVec( { 0.0, 29.5 } ),
		gen_complexforms_tris(),
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

	SysComp::InConstituent s__nap{
		SysComp::ConstituentType::NUCLEUS,
		createFixedString("S-NAP"),
		-1,
		0,
		mkRealVec( { 4.33 } ),
		mkRealVec( { 21.0, 0.0 } ),
		gen_complexforms_s__nap(),
		0.0
	};

	SysComp::InConstituent s__ibu{
		SysComp::ConstituentType::NUCLEUS,
		createFixedString("S-IBU"),
		-1,
		0,
		mkRealVec( { 4.45 } ),
		mkRealVec( { 19.6, 0.0 } ),
		gen_complexforms_s__ibu(),
		0.0
	};

	CMapping cBGE = {
		{ "Tricine-ANI", 20.0 },
		{ "TRIS", 20.0 },
		{ "b-CD", 15.0 }
	};

	CMapping cSample = {
		{ "Tricine-ANI", 15.0 },
		{ "TRIS", 15.0 },
		{ "b-CD", 13.0 },
		{ "S-NAP", 0.12 },
		{ "S-IBU", 0.12 }
	};

	const auto r = calculate(
		{
			tricine__ani,
			tris,
			b__cd
		},
		{
			tricine__ani,
			tris,
			b__cd,
			s__nap,
			s__ibu
		},
		cBGE, cSample,
		true, true, false, false);

	checkBGE(r, 8.1129381247, 0.051406531056, 0.010074614856, 23.02791337);

	checkEigenzone(r.eigenzones, -1.2506693969e-14, -2.9906850417e-15, 8.1129381247, 0.051406531055);

	checkEigenzone(r.eigenzones, 2.9282055721e-05, 2.2928095754e-05, 7.9891111057, 0.039712786546);

	checkEigenzone(r.eigenzones, -0.014871781076, -0.0053224091981, 8.2060634036, 0.050093783482);

	checkEigenzone(r.eigenzones, -9.5459923245, -0.12490381773, 8.1103325802, 0.051082495614);

	checkEigenzone(r.eigenzones, -10.139624181, -0.120653929, 8.1103344615, 0.051114301945);

	return EXIT_SUCCESS;
}
