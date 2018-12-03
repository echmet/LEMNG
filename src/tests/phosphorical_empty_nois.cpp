#include <cstdlib>
#include "barsarkagang_tests.h"


using namespace ECHMET;
using namespace ECHMET::Barsarkagang;


SysComp::InCFVec * gen_complexforms_phosphoric_acid()
{
	const ComplexDef cDef = {
		{ /* InComplexForm c-tor begin */
			-3,
			/* InLGVec */
			{
			}
		}, /* InComplexForm c-tor end */
		{ /* InComplexForm c-tor begin */
			-2,
			/* InLGVec */
			{
			}
		}, /* InComplexForm c-tor end */
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

SysComp::InCFVec * gen_complexforms_al()
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
		}, /* InComplexForm c-tor end */
		{ /* InComplexForm c-tor begin */
			2,
			/* InLGVec */
			{
			}
		}, /* InComplexForm c-tor end */
		{ /* InComplexForm c-tor begin */
			3,
			/* InLGVec */
			{
			}
		} /* InComplexForm c-tor end */
	};

	return buildComplexes(cDef);
}

int main(int , char ** )
{
	SysComp::InConstituent phosphoric_acid{
		SysComp::ConstituentType::NUCLEUS,
		createFixedString("Phosphoric acid"),
		-3,
		0,
		mkRealVec( { 12.67, 7.21, 2.16 } ),
		mkRealVec( { 71.5, 61.4, 34.6, 0.0 } ),
		gen_complexforms_phosphoric_acid(),
		0.0
	};

	SysComp::InConstituent al{
		SysComp::ConstituentType::NUCLEUS,
		createFixedString("Al"),
		0,
		3,
		mkRealVec( { 7.0, 6.0, 4.9 } ),
		mkRealVec( { 0.0, 21.1, 42.1, 63.2 } ),
		gen_complexforms_al(),
		0.0
	};

	CMapping cBGE = {
		{ "Phosphoric acid", 35.0 },
		{ "Al", 10.0 }
	};

	CMapping cSample = {
		{ "Phosphoric acid", 30.0 },
		{ "Al", 8.0 }
	};

	const auto r = calculate(
		{
			phosphoric_acid,
			al
		},
		{
			phosphoric_acid,
			al
		},
		cBGE, cSample,
		false, false, false, false);

	checkBGE(r, 3.0261380593, 0.31738941404, 0.06054709832, 10.972386383);

	checkEigenzone(r.eigenzones, 0.0088429041208, 0.0027173902486, 0.81583614629, 3.1025402684, 0.27055462645);

	checkEigenzone(r.eigenzones, 69.03696619, 3.6031523874, 1.102010988, 2.8522966522, 0.32016146344);

	return EXIT_SUCCESS;
}

