#include <cstdlib>
#include "barsarkagang_tests.h"


using namespace ECHMET;
using namespace ECHMET::Barsarkagang;


int main(int , char ** )
{
	SysComp::InConstituent sebacic_acid{
		SysComp::ConstituentType::LIGAND,
		createFixedString("Sebacic acid"),
		-2,
		0,
		mkRealVec( { 5.38, 4.53 } ),
		mkRealVec( { 44.9, 20.7, 0.0 } ),
		nullptr,
		0.0
	};

	SysComp::InConstituent imidazole{
		SysComp::ConstituentType::LIGAND,
		createFixedString("Imidazole"),
		0,
		1,
		mkRealVec( { 7.15 } ),
		mkRealVec( { 0.0, 52 } ),
		nullptr,
		0.0
	};

	CMapping cBGE = {
		{ "Sebacic acid", 0.21 },
		{ "Imidazole", 0.323 }
	};

	CMapping cSample = {
		{ "Sebacic acid", 0.1 },
		{ "Imidazole", 0.15 }
	};

	const auto r = calculate(
		{
			sebacic_acid,
			imidazole
		},
		{
			sebacic_acid,
			imidazole
		},
		cBGE, cSample,
		true, true, false, true);

	(void)r;

	return EXIT_SUCCESS;
}

