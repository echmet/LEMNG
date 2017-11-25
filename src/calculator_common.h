#ifndef ECHMET_LEMNG_CALCULATOR_COMMON_H
#define ECHMET_LEMNG_CALCULATOR_COMMON_H

#include "lemng_p.h"
#include "calculator_types.h"

namespace ECHMET {
namespace LEMNG {
namespace Calculator {

template<typename T>
constexpr T cxsgn(const T &v)
{
	return (v > 0) - (v < 0);
}

class CalculationException : public std::exception {
public:
    explicit CalculationException(const std::string &message, const RetCode errorCode);
    explicit CalculationException(std::string &&message, const RetCode errorCode) noexcept;

    const char * what() const noexcept;
    RetCode errorCode() const noexcept;

private:
    const RetCode m_errorCode;
    const std::string m_message;
};

SolutionProperties calculateSolutionProperties(const SysComp::ChemicalSystem *chemSystem, const RealVecPtr &concentrations, SysComp::CalculatedProperties *calcProps, const bool correctForIonicStrength, const bool calcBufferCapacity = false);
SolutionProperties calculateSolutionProperties(const ChemicalSystemPtr &chemSystem, const RealVecPtr &concentrations, CalculatedPropertiesPtr &calcProps, const bool correctForIonicStrength, const bool calcBufferCapacity = false);

template <typename T>
bool isComplex(const T &I);

RealVecPtr makeAnalyticalConcentrationsForDerivator(const CalculatorSystemPack &systemPack);
CalculatorSystemPack makeSystemPack(const ChemicalSystemPtr &chemSystem, const CalculatedPropertiesPtr &calcProps,
				    const std::function<bool (const std::string &)> &isAnalyte);
void prepareModelData(CalculatorSystemPack &systemPack, DeltaPackVec &deltaPacks, const RealVecPtr &analConcsBGELike, const RealVecPtr &analConcsSample, const bool correctForIonicStrength);
void solveChemicalSystem(const SysComp::ChemicalSystem *chemSystem, const RealVecPtr &concentrations, SysComp::CalculatedProperties *calcProps, const bool correctForIonicStrength);
void solveChemicalSystem(const ChemicalSystemPtr &chemSystem, const RealVecPtr &concentrations, CalculatedPropertiesPtr &calcProps, const bool correctForIonicStrength);
std::vector<const SysComp::Constituent *> sysCompToLEMNGOrdering(const ChemicalSystemPtr &chemSystem, const std::function<bool (const std::string &)> &isAnalyte);

#ifdef ECHMET_LEMNG_SENSITIVE_NUMDERS	/*!< Use much finer delta and lower analytes concentrations to calculate numerical derivatives. This comes with some additional memory and performance overhead */
static const ECHMETReal DELTA_H = 1.0e-33;
static const double ANALYTE_CONCENTRATION_NUMDERS = 1.0e-13;
#else
static const ECHMETReal DELTA_H = 1.0e-17;
static const double ANALYTE_CONCENTRATION_NUMDERS = 1.0e-13;
#endif // ECHMET_LEMNG_SENSITIVE_NUMDERS

static const double ANALYTE_CONCENTRATION = 1.0e-13;

} // namespace Calculator
} // namespace LEMNG
} // namespace ECHMET

#endif // ECHMET_LEMNG_CALCULATOR_COMMON_H
