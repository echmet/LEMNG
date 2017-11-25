#include "json_input_processor.h"
#include <algorithm>
#include <cstring>
#include <stdexcept>
#include <iostream>

namespace ECHMET {
namespace LEMNG {

void JsonInputProcessor::cleanupInLigandForm(SysComp::InLigandForm &lf)
{
	if (lf.ligandName)
		lf.ligandName->destroy();
	if (lf.pBs)
		lf.pBs->destroy();
	if (lf.mobilities)
		lf.mobilities->destroy();
}

void JsonInputProcessor::cleanupInLigands(SysComp::InLFVec *lfVec)
{
	for (size_t idx = 0; idx < lfVec->size(); idx++) {
		SysComp::InLigandForm &lf = (*lfVec)[idx];

		cleanupInLigandForm(lf);
	}

	lfVec->destroy();
}

void JsonInputProcessor::cleanupInLigandGroups(SysComp::InLGVec *lgVec)
{
	for (size_t idx = 0; idx < lgVec->size(); lgVec++) {
		if (lgVec->at(idx).ligands != NULL)
			cleanupInLigands((*lgVec)[idx].ligands);
	}

	lgVec->destroy();
}

void JsonInputProcessor::cleanupInComplexForms(SysComp::InCFVec *cfVec)
{
	for (size_t idx = 0; idx < cfVec->size(); idx++) {
		if (cfVec->at(idx).ligandGroups != NULL)
			cleanupInLigandGroups((*cfVec)[idx].ligandGroups);
	}

	cfVec->destroy();
}

void JsonInputProcessor::cleanupInConstituent(SysComp::InConstituent &c)
{
	if (c.mobilities)
		c.mobilities->destroy();

	if (c.pKas)
		c.pKas->destroy();

	if (c.name)
		c.name->destroy();

	if (c.complexForms)
		c.complexForms->destroy();
}

void JsonInputProcessor::cleanupInConstituentVector(SysComp::InConstituentVec *inCtuentVec)
{
	for (size_t idx = 0; idx < inCtuentVec->size(); idx++)
		cleanupInConstituent((*inCtuentVec)[idx]);

	inCtuentVec->destroy();
}

void JsonInputProcessor::makeSysCompComplexForms(SysComp::InConstituent &scCtuent, const constituent_t *ctuent, const std::vector<std::string> &listOfAnalytes)
{
	for (size_t idx = 0; idx < ctuent->complexFormsCount; idx++) {
		const complex_form_t *cForm = &ctuent->complexForms[idx];
		SysComp::InComplexForm scCF;

		scCF.ligandGroups = SysComp::createInLGVec(cForm->count);
		if (scCF.ligandGroups == NULL)
			throw std::runtime_error("Cannot create InLGVec");

		try {
			makeSysCompLigandGroups(scCF, cForm, listOfAnalytes);
		} catch (std::runtime_error &up) {
			cleanupInLigandGroups(scCF.ligandGroups);
			throw up;
		}

		scCF.nucleusCharge = cForm->nucleusCharge;

		if (scCtuent.complexForms->push_back(scCF) != OK) {
			cleanupInLigandGroups(scCF.ligandGroups);
			throw std::runtime_error("Cannot push back complex form");
		}
	}
}
bool isAnalyte(const std::vector<std::string> &listOfAnalytes, const std::string &name)
{
	for (std::vector<std::string>::const_iterator cit = listOfAnalytes.begin(); cit != listOfAnalytes.end(); cit++)
		if (*cit == name) return true;
	return false;
}

void zeroInitializeINC(SysComp::InConstituent &c)
{
	c.complexForms = NULL;
	c.mobilities = NULL;
	c.name = NULL;
	c.pKas = NULL;
	c.chargeLow = 0;
	c.chargeHigh = 0;
}

void JsonInputProcessor::makeSysCompInputInternal(SysComp::InConstituentVec *inCtuentVec, const constituent_array_t *input, const std::vector<std::string> &listOfAnalytes)
{
	for (size_t ctuentIdx = 0; ctuentIdx < input->count; ctuentIdx++) {
		const constituent_t *ctuent = &input->constituents[ctuentIdx];
		const int numpKas = ctuent->chargeHigh - ctuent->chargeLow;
		const int numMobilities = numpKas + 1;
		SysComp::InConstituent scCtuent;

		if (isAnalyte(listOfAnalytes, ctuent->name))
			continue;

		zeroInitializeINC(scCtuent);

		scCtuent.ctype = (ctuent->ctype == ::LIGAND) ? SysComp::LIGAND : SysComp::NUCLEUS;
		scCtuent.chargeLow = ctuent->chargeLow;
		scCtuent.chargeHigh = ctuent->chargeHigh;

		if (numpKas < 0) {
			cleanupInConstituentVector(inCtuentVec);
			throw std::runtime_error("Invalid charges");
		}

		scCtuent.pKas = createRealVec(numpKas);
		if (scCtuent.pKas == NULL) {
			cleanupInConstituentVector(inCtuentVec);
			throw std::runtime_error("Cannot create pKa vector");
		}

		scCtuent.mobilities = createRealVec(numMobilities);
		if (scCtuent.mobilities == NULL) {
			cleanupInConstituent(scCtuent);
			cleanupInConstituentVector(inCtuentVec);
			throw std::runtime_error("Cannot create mobilities vector");
		}

		scCtuent.complexForms = SysComp::createInCFVec(ctuent->complexFormsCount);
		if (scCtuent.complexForms == NULL) {
			cleanupInConstituent(scCtuent);
			cleanupInConstituentVector(inCtuentVec);
			throw std::runtime_error("Cannot create complexForms vector");
		}

		scCtuent.name = createFixedString(ctuent->name);
		if (scCtuent.name == NULL) {
			cleanupInConstituent(scCtuent);
			cleanupInConstituentVector(inCtuentVec);
			throw std::runtime_error("Cannot create constituent name");
		}

		try {
			if (scCtuent.ctype == SysComp::NUCLEUS)
				makeSysCompComplexForms(scCtuent, ctuent, listOfAnalytes);
		} catch (std::runtime_error &up) {
			cleanupInConstituent(scCtuent);
			cleanupInConstituentVector(inCtuentVec);
			throw up;
		}

		for (int pidx = 0; pidx < numpKas; pidx++) {
			if (scCtuent.pKas->push_back(ctuent->pKas[pidx]) != OK) {
				cleanupInConstituent(scCtuent);
				cleanupInConstituentVector(inCtuentVec);
				throw std::runtime_error("Cannot push back pKa");
			}
		}

		for (int midx = 0; midx < numMobilities; midx++) {
			if (scCtuent.mobilities->push_back(ctuent->mobilities[midx]) != OK) {
				cleanupInConstituent(scCtuent);
				cleanupInConstituentVector(inCtuentVec);
				throw std::runtime_error("Cannot push back constituent mobility");
			}
		}

		if (inCtuentVec->push_back(scCtuent) != OK) {
			cleanupInConstituent(scCtuent);
			cleanupInConstituentVector(inCtuentVec);
			throw std::runtime_error("Cannot push back constituent");
		}

	}
}

void JsonInputProcessor::makeSysCompInput(SysComp::InConstituentVec *&inCtuentVecBGE, SysComp::InConstituentVec *&inCtuentVecFull, ConcentrationMap &inBGEConcentrations, ConcentrationMap &inFullConcentrations, const constituent_array_t *input)
{
	std::vector<std::string> listOfAnalytes;

	if (input->count < 1)
		throw std::runtime_error("Input array does not contain any constituents");

	inCtuentVecBGE = SysComp::createInConstituentVec(input->count);
	if (inCtuentVecBGE == NULL)
		throw std::runtime_error("Cannot create SysComp::InConstituentVec for BGE");

	inCtuentVecFull = SysComp::createInConstituentVec(input->count);
	if (inCtuentVecFull == NULL) {
		inCtuentVecBGE->destroy();
		throw std::runtime_error("Cannot create SysComp::InConstituentVec for full system");
	}

	for (size_t idx = 0; idx < input->count; idx++) {
		const constituent_t *ctuent = &input->constituents[idx];

		try {
			if (ctuent->crole == ::ANALYTE) {
				listOfAnalytes.push_back(ctuent->name);
				inFullConcentrations[ctuent->name] = ctuent->concentrationSample;
			} else {
				inBGEConcentrations[ctuent->name] = ctuent->concentrationBGE;
				inFullConcentrations[ctuent->name] = ctuent->concentrationSample;
			}
		} catch (std::bad_alloc &) {
			throw std::runtime_error("Cannot emplace concentration mapping");
		}
	}

	try {
		makeSysCompInputInternal(inCtuentVecBGE, input, listOfAnalytes);
	} catch (...) {
		inCtuentVecBGE->destroy();
		inCtuentVecFull->destroy();
		throw;
	}
	try {
		makeSysCompInputInternal(inCtuentVecFull, input, std::vector<std::string>());
	} catch (...) {
		releaseInputData(inCtuentVecBGE);
		inCtuentVecBGE->destroy();
		inCtuentVecFull->destroy();
		throw;
	}
}

void zeroInitialize(SysComp::InLigandForm &scLF) {
	scLF.ligandName = NULL;
	scLF.pBs = NULL;
	scLF.mobilities = NULL;
	scLF.maxCount = 0;
}

void JsonInputProcessor::makeSysCompLigands(SysComp::InLigandGroup &scLG, const ligand_group_t *lGroup, const std::vector<std::string> &listOfAnalytes)
{
	for (size_t idx = 0; idx < lGroup->count; idx++) {
		const ligand_form_t *lForm = &lGroup->ligandForms[idx];
		SysComp::InLigandForm scLF;

		if (std::find(listOfAnalytes.begin(), listOfAnalytes.end(), std::string(lForm->name)) != listOfAnalytes.end())
			continue;

		zeroInitialize(scLF);

		scLF.pBs = createRealVec(lForm->maxCount);
		if (scLF.pBs == NULL)
			throw std::runtime_error("Cannot create pBs vector");

		for (int pbidx = 0; pbidx < lForm->maxCount; pbidx++) {
			if (scLF.pBs->push_back(lForm->pBs[pbidx]) != OK) {
				cleanupInLigandForm(scLF);
				throw std::runtime_error("Cannot push back pB");
			}
		}

		scLF.mobilities = createRealVec(lForm->maxCount);
		if (scLF.mobilities == NULL) {
			cleanupInLigandForm(scLF);
			throw std::runtime_error("Cannot create mobilities vector");
		}

		for (int midx = 0; midx < lForm->maxCount; midx++) {
			if (scLF.mobilities->push_back(lForm->mobilities[midx]) != OK) {
				cleanupInLigandForm(scLF);
				throw std::runtime_error("Cannot push back ligand form mobility");
			}
		}

		scLF.ligandName = createFixedString(lForm->name);
		if (scLF.ligandName == NULL) {
			cleanupInLigandForm(scLF);
			throw std::runtime_error("Cannot create ligand name");
		}

		scLF.charge = lForm->charge;
		scLF.maxCount = lForm->maxCount;

		if (scLG.ligands->push_back(scLF) != OK) {
			cleanupInLigandForm(scLF);
			throw std::runtime_error("Cannot push back ligand form");
		}
	}
}

void JsonInputProcessor::makeSysCompLigandGroups(SysComp::InComplexForm &scCF, const complex_form_t *cForm, const std::vector<std::string> &listOfAnalytes)
{
	for (size_t idx = 0; idx < cForm->count; idx++) {
		const ligand_group_t *lGroup = &cForm->ligandGroups[idx];
		SysComp::InLigandGroup scLG;

		scLG.ligands = SysComp::createInLFVec(lGroup->count);
		if (scLG.ligands == NULL) {
			cleanupInLigandGroups(scCF.ligandGroups);
			throw std::runtime_error("Cannot create InLFVec");
		}

		try {
			makeSysCompLigands(scLG, lGroup, listOfAnalytes);
		} catch (std::runtime_error &up) {
			cleanupInLigandGroups(scCF.ligandGroups);
			throw up;
		}

		if (scCF.ligandGroups->push_back(scLG) != OK) {
			cleanupInLigands(scLG.ligands);
			throw std::runtime_error("Cannot push back ligand group");
		}
	}
}
JsonInputProcessor::InputDescription JsonInputProcessor::process(const constituent_array_t *input)
{
	SysComp::InConstituentVec *inCtuentVecBGE;
	SysComp::InConstituentVec *inCtuentVecFull;
	ConcentrationMap inBGEConcentrations;
	ConcentrationMap inFullConcentrations;

	makeSysCompInput(inCtuentVecBGE, inCtuentVecFull, inBGEConcentrations, inFullConcentrations, input);

	return InputDescription(inCtuentVecBGE, inCtuentVecFull, inBGEConcentrations, inFullConcentrations);
}

} // namespace LEMNG
} // namespace ECHMET
