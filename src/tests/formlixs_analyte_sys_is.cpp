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
		{ "X", 2.0 },
		{ "S", 10.0 }
	};

	CMapping cSample = {
		{ "Formic acid", 5.0 },
		{ "Li", 2.5 },
		{ "X", 1e-12 },
		{ "S", 1e-12 }
	};

	auto r = calculate(
		{
			formic_acid,
			li,
			x,
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

	checkBGE(r, 3.4308189976, 0.050030135937, 0.0054002657348, 6.0888426984);

	checkEigenzone(r.eigenzones, -0.012316464932, -0.071355066634, 0.48439926692, 3.3886626501, 0.058292552093);

	checkEigenzone(r.eigenzones, 1.6600164232, 1.6372366765, 0.75444110326, 11.120187788, 0.033617883459);

	checkEigenzone(r.eigenzones, -10.510789143, 2.9183155611, 0.68044613415, 3.4355298748, 0.056380660082);

	checkEigenzone(r.eigenzones, 35.059371768, 3.2697570635, 1.0878447758, 3.3705858237, 0.049991341944);

	LEMNG::releaseResults(r);
	SysComp::releaseInConstituent(formic_acid);
	SysComp::releaseInConstituent(li);
	SysComp::releaseInConstituent(x);
	SysComp::releaseInConstituent(s);

	return EXIT_SUCCESS;
}

