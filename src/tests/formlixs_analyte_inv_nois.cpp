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
		false, false, false, false);

	checkBGE(r, 3.7807882258, 0.053349365978, 0.0051656577556, 6.1315850881);

	checkEigenzone(1, r.eigenzones, -1.6124697662e-15, -1.2080938932e-16, 0.51385140826, 3.7807882258, 0.053349365977);

	checkEigenzone(2, r.eigenzones, -4.5669268989e-07, -5.2707166877e-07, 1.3420079092, 4.0689839245, 0.034619983611);

	checkEigenzone(3, r.eigenzones, -10.163934426, -0.87428467248, 0.26113760092, 3.7777349796, 0.050532309836);

	checkEigenzone(4, r.eigenzones, 19.791803184, 5.5740910363, 1.1569824945, 3.6311800827, 0.050225836356);

	LEMNG::releaseResults(r);
	SysComp::releaseInConstituent(formic_acid);
	SysComp::releaseInConstituent(li);
	SysComp::releaseInConstituent(x);
	SysComp::releaseInConstituent(s_BGE);
	SysComp::releaseInConstituent(s);

	return EXIT_SUCCESS;
}

