#ifndef LEMNG_BASE_TYPES_H
#define LEMNG_BASE_TYPES_H

#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef ECHMET_IMPORT_INTERNAL
#define ECHMET_IMPORT_INTERNAL
#endif // ECHMET_IMPORT_INTERNAL
#include <echmetsyscomp.h>

namespace ECHMET {
namespace LEMNG {

void chemicalSystemDeleter(SysComp::ChemicalSystem *p);
void calculatedPropertiesDeleter(SysComp::CalculatedProperties *p);
void echmetRealVecDeleter(RealVec *p);

typedef std::map<std::string, bool> IsAnalyteMap;
typedef std::unique_ptr<SysComp::ChemicalSystem, decltype(&chemicalSystemDeleter)> ChemicalSystemPtr;
typedef std::unique_ptr<SysComp::CalculatedProperties, decltype(&calculatedPropertiesDeleter)> CalculatedPropertiesPtr;
typedef std::unique_ptr<RealVec, decltype(&echmetRealVecDeleter)> RealVecPtr;

enum class ConstituentRole {
	ANALYTE,
	BACKGROUND
};

class SysCompException : std::exception {
public:
	SysCompException(const std::string &message, const ::ECHMET::RetCode errorCode);
	::ECHMET::RetCode errorCode() const noexcept;
	const char * what() const noexcept;

private:
	const std::string m_message;
	const ::ECHMET::RetCode m_errorCode;
};

RealVecPtr makeAnalyticalConcentrationsVec(const ChemicalSystemPtr &chemSystem);
RealVecPtr makeAnalyticalConcentrationsVec(const SysComp::ChemicalSystem *chemSystemRaw);
CalculatedPropertiesPtr makeCalculatedProperties(const SysComp::ChemicalSystem *chemSystemRaw);
IsAnalyteMap makeIsAnalyteMap(const SysComp::InConstituentVec *BGEVec, const SysComp::InConstituentVec *FullVec);

} // namespace LEMNG
} // namespace ECHMET

#endif // LEMNG_BASE_TYPES_H
