#include <cstdlib>
#include "barsarkagang_tests.h"


using namespace ECHMET;
using namespace ECHMET::Barsarkagang;


SysComp::InCFVec * gen_complexforms_benzoic_acid()
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
							{ -1.348304863048161 },
							{ 9.9 }
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
							{ -2.518513939877887 },
							{ 0 }
						} /* InLigandForm c-tor end */
					}
				} /* InLigandGroup c-tor end */
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

int main(int , char ** )
{
	SysComp::InConstituent benzoic_acid{
		SysComp::ConstituentType::NUCLEUS,
		createFixedString("Benzoic acid"),
		-1,
		0,
		mkRealVec( { 4.203 } ),
		mkRealVec( { 33.6, 0.0 } ),
		gen_complexforms_benzoic_acid(),
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
		0.003
	};

	CMapping cBGE = {
		{ "Benzoic acid", 35.0 },
		{ "Li", 10.0 },
		{ "b-CD", 10.0 }
	};

	CMapping cSample = {
		{ "Benzoic acid", 30.0 },
		{ "Li", 8.0 },
		{ "b-CD", 10.0 }
	};

	const auto r = calculate(
		{
			benzoic_acid,
			li,
			b__cd
		},
		{
			benzoic_acid,
			li,
			b__cd
		},
		cBGE, cSample,
		false, false, false, false);

	checkBGE(r, 3.9706700188, 0.074457266563, 0.010106986745, 12.614407061);

	checkEigenzone(r.eigenzones, -8.7070560173e-08, -2.054598111e-08, 4.0001277163, 0.069854385856);

	checkEigenzone(r.eigenzones, -0.60094175857, -0.094933755149, 3.9905125767, 0.069862350483);

	checkEigenzone(r.eigenzones, 6.372230561, 1.4133711229, 3.9097796852, 0.069633985664);

	return EXIT_SUCCESS;
}

