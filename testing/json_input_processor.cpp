#include "json_input_processor.h"
#include <algorithm>
#include <cstring>
#include <stdexcept>
#include <iostream>

namespace ECHMET {

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

void JsonInputProcessor::makeSysCompComplexForms(SysComp::InConstituent &scCtuent, const constituent_t *ctuent)
{
	for (size_t idx = 0; idx < ctuent->complexFormsCount; idx++) {
		const complex_form_t *cForm = &ctuent->complexForms[idx];
		SysComp::InComplexForm scCF;

		scCF.ligandGroups = SysComp::createInLGVec(cForm->count);
		if (scCF.ligandGroups == NULL)
			throw std::runtime_error("Cannot create InLGVec");

		try {
			makeSysCompLigandGroups(scCF, cForm);
		} catch (std::runtime_error &up) {
			cleanupInLigandGroups(scCF.ligandGroups);
			throw up;
		}

		scCF.nucleusCharge = cForm->nucleusCharge;

		if (scCtuent.complexForms->push_back(scCF) != RetCode::OK) {
			cleanupInLigandGroups(scCF.ligandGroups);
			throw std::runtime_error("Cannot push back complex form");
		}
	}
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

void JsonInputProcessor::makeSysCompInputInternal(SysComp::InConstituentVec *inCtuentVec, const constituent_array_t *input)
{
	for (size_t ctuentIdx = 0; ctuentIdx < input->count; ctuentIdx++) {
		const constituent_t *ctuent = &input->constituents[ctuentIdx];
		const int numpKas = ctuent->chargeHigh - ctuent->chargeLow;
		const int numMobilities = numpKas + 1;
		SysComp::InConstituent scCtuent;

		zeroInitializeINC(scCtuent);

		scCtuent.ctype = (ctuent->ctype == ::LIGAND) ? SysComp::ConstituentType::LIGAND : SysComp::ConstituentType::NUCLEUS;
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
			if (scCtuent.ctype == SysComp::ConstituentType::NUCLEUS)
				makeSysCompComplexForms(scCtuent, ctuent);
		} catch (std::runtime_error &up) {
			cleanupInConstituent(scCtuent);
			cleanupInConstituentVector(inCtuentVec);
			throw up;
		}

		for (int pidx = 0; pidx < numpKas; pidx++) {
			if (scCtuent.pKas->push_back(ctuent->pKas[pidx]) != RetCode::OK) {
				cleanupInConstituent(scCtuent);
				cleanupInConstituentVector(inCtuentVec);
				throw std::runtime_error("Cannot push back pKa");
			}
		}

		for (int midx = 0; midx < numMobilities; midx++) {
			if (scCtuent.mobilities->push_back(ctuent->mobilities[midx]) != RetCode::OK) {
				cleanupInConstituent(scCtuent);
				cleanupInConstituentVector(inCtuentVec);
				throw std::runtime_error("Cannot push back constituent mobility");
			}
		}

		if (inCtuentVec->push_back(scCtuent) != RetCode::OK) {
			cleanupInConstituent(scCtuent);
			cleanupInConstituentVector(inCtuentVec);
			throw std::runtime_error("Cannot push back constituent");
		}

	}
}

void JsonInputProcessor::makeSysCompInput(SysComp::InConstituentVec *&inCtuentVecBGE, ConcentrationMap &inBGEConcentrations, const constituent_array_t *input)
{
	//if (input->count < 1)
	//	throw std::runtime_error("Input array does not contain any constituents");

	inCtuentVecBGE = SysComp::createInConstituentVec(input->count);
	if (inCtuentVecBGE == NULL)
		throw std::runtime_error("Cannot create SysComp::InConstituentVec for BGE");

	for (size_t idx = 0; idx < input->count; idx++) {
		const constituent_t *ctuent = &input->constituents[idx];
		inBGEConcentrations[ctuent->name] = ctuent->concentrationBGE;
	}

	try {
		makeSysCompInputInternal(inCtuentVecBGE, input);
	} catch (...) {
		inCtuentVecBGE->destroy();
		throw;
	}
}

void zeroInitialize(SysComp::InLigandForm &scLF)
{
	scLF.ligandName = NULL;
	scLF.pBs = NULL;
	scLF.mobilities = NULL;
	scLF.maxCount = 0;
}

void JsonInputProcessor::makeSysCompLigands(SysComp::InLigandGroup &scLG, const ligand_group_t *lGroup)
{
	for (size_t idx = 0; idx < lGroup->count; idx++) {
		const ligand_form_t *lForm = &lGroup->ligandForms[idx];
		SysComp::InLigandForm scLF;

		zeroInitialize(scLF);

		scLF.pBs = createRealVec(lForm->maxCount);
		if (scLF.pBs == NULL)
			throw std::runtime_error("Cannot create pBs vector");

		for (int pbidx = 0; pbidx < lForm->maxCount; pbidx++) {
			if (scLF.pBs->push_back(lForm->pBs[pbidx]) != RetCode::OK) {
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
			if (scLF.mobilities->push_back(lForm->mobilities[midx]) != RetCode::OK) {
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

		if (scLG.ligands->push_back(scLF) != RetCode::OK) {
			cleanupInLigandForm(scLF);
			throw std::runtime_error("Cannot push back ligand form");
		}
	}
}

void JsonInputProcessor::makeSysCompLigandGroups(SysComp::InComplexForm &scCF, const complex_form_t *cForm)
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
			makeSysCompLigands(scLG, lGroup);
		} catch (std::runtime_error &up) {
			cleanupInLigandGroups(scCF.ligandGroups);
			throw up;
		}

		if (scCF.ligandGroups->push_back(scLG) != RetCode::OK) {
			cleanupInLigands(scLG.ligands);
			throw std::runtime_error("Cannot push back ligand group");
		}
	}
}
JsonInputProcessor::InputDescription JsonInputProcessor::process(const constituent_array_t *input)
{
	SysComp::InConstituentVec *inCtuentVecBGE;
	ConcentrationMap inBGEConcentrations;

	makeSysCompInput(inCtuentVecBGE, inBGEConcentrations, input);

	return InputDescription(inCtuentVecBGE, inBGEConcentrations);
}

} // namespace ECHMET
