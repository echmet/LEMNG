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

	auto r = calculate(
		{
			phosphoric_acid,
			al
		},
		{
			phosphoric_acid,
			al
		},
		cBGE, cSample,
		true, true, false, false);

	checkBGE(r, 2.9887541279, 0.22812330017, 0.061116910773, 10.790315473);

	checkEigenzone(r.eigenzones, -0.0064462519304, -0.0020480216542, 0.5633025422, 3.0695578841, 0.19782992674);

	checkEigenzone(r.eigenzones, 75.226960638, -3.6656655989, 0.81931011765, 2.8166063213, 0.24076472939);

	LEMNG::releaseResults(r);
	SysComp::releaseInConstituent(phosphoric_acid);
	SysComp::releaseInConstituent(al);

	return EXIT_SUCCESS;
}

