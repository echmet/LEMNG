#ifndef ECHMET_ECHMET_LEMNG_P_H
#define ECHMET_ECHMET_LEMNG_P_H

#include <lemng.h>
#include "base_types.h"
#include "calculator_types.h"

namespace ECHMET {
namespace LEMNG {

class CZESystemImpl : public CZESystem {
public:
	explicit CZESystemImpl();
	explicit CZESystemImpl(CZESystemImpl &&other) noexcept;
	explicit CZESystemImpl(const SysComp::ChemicalSystem &chemicalSystemBGE, const SysComp::CalculatedProperties &calcPropsBGE, const SysComp::ChemicalSystem &chemicalSystemFull, const SysComp::CalculatedProperties &calcPropsFull, const IsAnalyteMap &iaMap);
	explicit CZESystemImpl(const SysComp::ChemicalSystem &chemicalSystemBGE, const SysComp::CalculatedProperties &calcPropsBGE, const SysComp::ChemicalSystem &chemicalSystemFull, const SysComp::CalculatedProperties& calcPropsFull, IsAnalyteMap &&iaMap);
	virtual ~CZESystemImpl() noexcept override;
	virtual RetCode ECHMET_CC evaluate(const InAnalyticalConcentrationsMap *acBGE, const InAnalyticalConcentrationsMap *acFull,
					   const NonidealityCorrections corrections, Results &results) noexcept override;
	virtual const char * ECHMET_CC lastErrorString() const noexcept override;
	virtual RetCode ECHMET_CC makeAnalyticalConcentrationsMaps(InAnalyticalConcentrationsMap *&acMapBGE, InAnalyticalConcentrationsMap *&acMapFull) const noexcept override;

	static CZESystemImpl * make(const SysComp::InConstituentVec *inCtuentVecBGE, const SysComp::InConstituentVec *inCtuentVecSample);

private:
	bool isAnalyte(const std::string &name);
	void setupInternal(const SysComp::ChemicalSystem &chemicalSystemBGE, const SysComp::CalculatedProperties &calcPropsBGE, const SysComp::ChemicalSystem &chemicalSystemFull, const SysComp::CalculatedProperties &calcPropsFull);

	ChemicalSystemPtr m_chemicalSystemBGE;
	ChemicalSystemPtr m_chemicalSystemFull;
	CalculatedPropertiesPtr m_calcPropsBGE;
	CalculatedPropertiesPtr m_calcPropsFull;
	Calculator::CalculatorSystemPack m_systemPack;
	Calculator::CalculatorSystemPack m_systemPackUncharged;

	IsAnalyteMap m_isAnalyteMap;

	std::string m_lastErrorString;
};

} // namespace LEMNG
} // namesoace ECHMET

#endif // ECHMET_ECHMET_LEMNG_P_H
