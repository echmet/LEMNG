#include "base_types.h"
#include <algorithm>
#include <cstring>

namespace ECHMET {
namespace LEMNG {

SysCompException::SysCompException(const std::string &message, const ::ECHMET::RetCode errorCode) :
	m_message{message + ": " + ::ECHMET::errorToString(errorCode)},
	m_errorCode{errorCode}
{}

::ECHMET::RetCode SysCompException::errorCode() const noexcept
{
	return m_errorCode;
}

const char * SysCompException::what() const noexcept
{
	return m_message.c_str();
}

void chemicalSystemDeleter(SysComp::ChemicalSystem *p)
{
	SysComp::releaseChemicalSystem(*p);
	delete p;
}

void calculatedPropertiesDeleter(SysComp::CalculatedProperties *p)
{
	SysComp::releaseCalculatedProperties(*p);
	delete p;
}

void echmetRealVecDeleter(RealVec *p)
{
	if (p != nullptr)
		p->destroy();
}

RealVecPtr makeAnalyticalConcentrationsVec(const ChemicalSystemPtr &chemSystem)
{
	return makeAnalyticalConcentrationsVec(chemSystem.get());
}

RealVecPtr makeAnalyticalConcentrationsVec(const SysComp::ChemicalSystem *chemSystemRaw)
{
	RealVec *vec;

	::ECHMET::RetCode tRet = SysComp::makeAnalyticalConcentrationsVec(vec, *chemSystemRaw);
	if (tRet != ::ECHMET::RetCode::OK)
		throw std::bad_alloc{};

	return std::unique_ptr<RealVec, decltype(&echmetRealVecDeleter)>{vec, echmetRealVecDeleter};
}

CalculatedPropertiesPtr makeCalculatedProperties(const SysComp::ChemicalSystem *chemSystem)
{
	SysComp::CalculatedProperties *calcProps = new SysComp::CalculatedProperties{};

	::ECHMET::RetCode tRet = SysComp::initializeCalculatedProperties(*calcProps, *chemSystem);
	if (tRet != ::ECHMET::RetCode::OK) {
		delete calcProps;
		throw SysCompException{"Cannot initialize CalculatedProperties", tRet};
	}

	return std::unique_ptr<SysComp::CalculatedProperties, decltype(&calculatedPropertiesDeleter)>(calcProps, calculatedPropertiesDeleter);
}

IsAnalyteMap makeIsAnalyteMap(const SysComp::InConstituentVec *BGEVec, const SysComp::InConstituentVec *FullVec)
{
	auto isAnalyte = [BGEVec](const FixedString *name) {
		for (size_t idx = 0; idx < BGEVec->size(); idx++) {
			const SysComp::InConstituent &c = BGEVec->at(idx);
			if (*(c.name) == (*name))
				return false;
		}

		return true;
	};

	IsAnalyteMap iaMap{};

	for (size_t idx = 0; idx < FullVec->size(); idx++) {
		const SysComp::InConstituent &scC = FullVec->at(idx);
		const FixedString *name = scC.name;

		iaMap.emplace(std::string(name->c_str()), isAnalyte(name));
	}

	return iaMap;
}

} // namespace LEMNG
} // namespace ECHMET
