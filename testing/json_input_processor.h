#ifndef JSON_INPUT_PROCESSOR_H
#define JSON_INPUT_PROCESSOR_H

#include "jsonloader/constituents_json_ldr.h"
#include <echmetsyscomp.h>
#include <map>
#include <string>
#include <vector>

namespace ECHMET {

class JsonInputProcessor {
public:
	typedef std::map<std::string, double> ConcentrationMap;
	class InputDescription {
	public:
		InputDescription() :
			BGEComposition(NULL),
			BGEConcentrations() {}
		InputDescription(SysComp::InConstituentVec *_BGEComposition,
				 const ConcentrationMap  &_BGEConcentrations) :
		BGEComposition(_BGEComposition),
		BGEConcentrations(_BGEConcentrations)
		{
		}

		SysComp::InConstituentVec *BGEComposition;
		ConcentrationMap BGEConcentrations;
	};

	InputDescription process(const constituent_array_t *input);

private:
	static void cleanupInComplexForms(SysComp::InCFVec *cfVec);
	static void cleanupInLigandForm(SysComp::InLigandForm &lf);
	static void cleanupInLigands(SysComp::InLFVec *lfVec);
	static void cleanupInLigandGroups(SysComp::InLGVec *lgVec);
	static void cleanupInConstituent(SysComp::InConstituent &c);
	static void cleanupInConstituentVector(SysComp::InConstituentVec *inCtuentVec);
	static void makeSysCompComplexForms(SysComp::InConstituent &scCtuent, const constituent_t *ctuent);
	static void makeSysCompInput(SysComp::InConstituentVec *&inCtuentVecBGE, ConcentrationMap &inBGEConcentrations, const constituent_array_t *input);
	static void makeSysCompInputInternal(SysComp::InConstituentVec *inCtuentVec, const constituent_array_t *input);
	static void makeSysCompLigands(SysComp::InLigandGroup &scLG, const ligand_group_t *lGroup);
	static void makeSysCompLigandGroups(SysComp::InComplexForm &scCF, const complex_form_t *cForm);
};

} // namespace ECHMET

#endif // JSON_INPUT_PROCESSOR_H

