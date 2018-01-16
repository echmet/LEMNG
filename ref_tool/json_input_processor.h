#ifndef JSON_INPUT_PROCESSOR_H
#define JSON_INPUT_PROCESSOR_H

#include "jsonloader/constituents_json_ldr.h"
#include <echmetsyscomp.h>
#include <map>
#include <string>
#include <vector>

namespace ECHMET {
namespace LEMNG {

class JsonInputProcessor {
public:
	typedef std::map<std::string, double> ConcentrationMap;
	class InputDescription {
	public:
		InputDescription() :
			BGEComposition(NULL),
			SampleComposition(NULL) {}
		InputDescription(SysComp::InConstituentVec *_BGEComposition,
				 SysComp::InConstituentVec *_SampleComposition,
				 const ConcentrationMap  &_BGEConcentrations,
				 const ConcentrationMap &_SampleConcentrations) :
		BGEComposition(_BGEComposition),
		SampleComposition(_SampleComposition),
		BGEConcentrations(_BGEConcentrations),
		SampleConcentrations(_SampleConcentrations) {}

		SysComp::InConstituentVec *BGEComposition;
		SysComp::InConstituentVec *SampleComposition;
		ConcentrationMap BGEConcentrations;
		ConcentrationMap SampleConcentrations;
	};

	InputDescription process(const constituent_array_t *input);

private:
	static void cleanupInComplexForms(SysComp::InCFVec *cfVec);
	static void cleanupInLigandForm(SysComp::InLigandForm &lf);
	static void cleanupInLigands(SysComp::InLFVec *lfVec);
	static void cleanupInLigandGroups(SysComp::InLGVec *lgVec);
	static void cleanupInConstituent(SysComp::InConstituent &c);
	static void cleanupInConstituentVector(SysComp::InConstituentVec *inCtuentVec);
	static void makeSysCompComplexForms(SysComp::InConstituent &scCtuent, const constituent_t *ctuent, const std::vector<std::string> &listOfAnalytes);
	static void makeSysCompInput(SysComp::InConstituentVec *&inCtuentVecBGE, SysComp::InConstituentVec *&inCtuentFull, ConcentrationMap &inBGEConcentrations, ConcentrationMap &inFullConcentrations, const constituent_array_t *input);
	static void makeSysCompInputInternal(SysComp::InConstituentVec *inCtuentVec, const constituent_array_t *input, const std::vector<std::string> &listOfAnalytes);
	static void makeSysCompLigands(SysComp::InLigandGroup &scLG, const ligand_group_t *lGroup, const std::vector<std::string> &listOfAnalytes);
	static void makeSysCompLigandGroups(SysComp::InComplexForm &scCF, const complex_form_t *cForm, const std::vector<std::string> &listOfAnalytes);
};

} // namespace LEMNG
} // namespace ECHMET

#endif // JSON_INPUT_PROCESSOR_H

