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
		false, false, false, false);

	checkBGE(r, 3.4544554959, 0.050948456835, 0.0053511919113, 5.9391385953);

	checkEigenzone(1, r.eigenzones, -0.012266406199, -0.070470519961, 0.48236767171, 3.4132556214, 0.059908842156);

	checkEigenzone(2, r.eigenzones, 1.7080460734, 1.7339668748, 0.79204805167, 11.151567718, 0.034717917459);

	checkEigenzone(3, r.eigenzones, -11.033657892, 3.2374034446, 0.72045431567, 3.4646248412, 0.058084284903);

	checkEigenzone(4, r.eigenzones, 33.42827756, 4.2339445103, 1.1643726698, 3.3831603772, 0.050576459137);

	LEMNG::releaseResults(r);
	SysComp::releaseInConstituent(formic_acid);
	SysComp::releaseInConstituent(li);
	SysComp::releaseInConstituent(x);
	SysComp::releaseInConstituent(s);

	return EXIT_SUCCESS;
}

