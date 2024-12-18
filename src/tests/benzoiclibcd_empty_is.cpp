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

	auto r = calculate(
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
		true, true, false, false);

	checkBGE(r, 3.9286783875, 0.068453027608, 0.01013026897, 12.67373577);

	checkEigenzone(1, r.eigenzones, -8.3000447542e-08, -4.2578811524e-05, 0.53089249474, 3.9609903249, 0.064158451661);

	checkEigenzone(2, r.eigenzones, -0.57626859052, -0.086945550121, 0.54118030447, 3.9494636186, 0.064681343272);

	checkEigenzone(3, r.eigenzones, 7.186583962, 1.5900667866, 0.93040426982, 3.8681000409, 0.064339288371);

	LEMNG::releaseResults(r);
	SysComp::releaseInConstituent(benzoic_acid);
	SysComp::releaseInConstituent(li);
	SysComp::releaseInConstituent(b__cd);

	return EXIT_SUCCESS;
}

