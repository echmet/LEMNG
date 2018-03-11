#ifndef ECHMET_LEMNG_CALCULATOR_TYPES_H
#define ECHMET_LEMNG_CALCULATOR_TYPES_H

#include <complex>
#include <memory>
#include <vector>
#include <Eigen/Dense>

namespace ECHMET {

/* Forward-declare SysComp data types */
namespace SysComp {
	class CalculatedProperties;
	class Constituent;
	class ChemicalSystem;
	class IonicForm;
} // namespace SysComp

namespace LEMNG {
namespace Calculator {

typedef std::vector<double> ERVector;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> EMMatrix;
typedef Eigen::Matrix<double, 1, Eigen::Dynamic>  EMVector;
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>  EMMatrixC;
typedef Eigen::Matrix<std::complex<double>, 1, Eigen::Dynamic>  EMVectorC;

typedef std::vector<const SysComp::Constituent *> InternalConstituentVec;

/*!
 * Ionic form representation used internally by the \p Calculator
 */
class CalculatorIonicForm {
public:
	CalculatorIonicForm(const std::string &name, const int32_t charge,
			    const SysComp::IonicForm *internalIonicForm,
			    const size_t internalIonicFormConcentrationIdx,
			    const size_t globalIonicFormConcentrationIdx,
			    const std::vector<size_t> &containedConstituents, const bool isAnaylte);
	CalculatorIonicForm(std::string &&name, const int32_t charge,
			    const SysComp::IonicForm *internalIonicForm,
			    const size_t internalIonicFormConcentrationIdx,
			    const size_t globalIonicFormConcentrationIdx,
			    std::vector<size_t> &&containedConstituents, const bool isAnalyte) noexcept;

	const std::string name;					/*!< Name of the ionic form, useful only for debugging purposes */
	const int32_t charge;					/*!< Total electric charge of the ionic form */
	const SysComp::IonicForm *internalIonicForm;		/*!< Pointer to the underlying \p SysComp::IonicForm */
	const size_t internalIonicFormConcentrationIdx;		/*!< Corresponding index in vector of ionic concentrations used by CoreLibs */
	const size_t globalIonicFormConcentrationIdx;		/*!< Corresponding index in vector of ionic concentrations used by Calculator */
	const std::vector<size_t> containedConstituents;
	const bool isAnalyte;					/*!< Ionic form is a form of an analyte */

	double concentration;					/*!< Concentration of the ionic form. This cannot be set by the c-tor because the
								     internal system representation is reusable. */
	double mobility;					/*!< Actual ionic mobility of the ionic form. This cannot be set by the c-tor because
								     the internal system representation is reusable. */

	CalculatorIonicForm & operator=(const CalculatorIonicForm &other);
	CalculatorIonicForm & operator=(CalculatorIonicForm &&other) noexcept;
};
typedef std::vector<CalculatorIonicForm *> CalculatorIonicFormVec;

/*!
 * Constituent representation used internally by the \p Calculator
 */
class CalculatorConstituent {
public:
	CalculatorConstituent(const std::string &name, const CalculatorIonicFormVec &ifVec, const SysComp::Constituent *internalConstituent, const bool isAnalyte);
	CalculatorConstituent(std::string &&name, CalculatorIonicFormVec &&ifVec, const SysComp::Constituent *internalConstituent, const bool isAnalyte) noexcept;
	CalculatorConstituent(const CalculatorConstituent &other);
	CalculatorConstituent(CalculatorConstituent &&other) noexcept;

	CalculatorConstituent & operator=(const CalculatorConstituent &other);
	CalculatorConstituent & operator=(CalculatorConstituent &&other) noexcept;

	const std::string name;					/*!< Name of the constituent, used only for debugging purposes */
	const CalculatorIonicFormVec ionicForms;		/*!< Ionic forms that contain a given constituent */
	const SysComp::Constituent * const internalConstituent;	/*!< Pointer to the underlying \p SysComp::Constituent object */
	const bool isAnalyte;					/*!< Set to true if the constituent is an analyte in a given system */

	double concentrationBGE;				/*<! Actual analytical concentration of the constituent in the background electrolyte.
								     This cannot be set by the c-tor because the internal system representation is reusable. */
	double concentrationSample;				/*<! Actual analytical concentration of the constituent in the background electrolyte.
								     This cannot be set by the c-tor because the internal system representation is reusable. */
};
typedef std::vector<CalculatorConstituent> CalculatorConstituentVec;

/*!
 * Representation of the entire chemical system used internally by the calculator
 */
class CalculatorSystemPack {
public:
	CalculatorSystemPack();
	CalculatorSystemPack(const CalculatorConstituentVec &ccVec, const CalculatorIonicFormVec &ifVec, const SysComp::ChemicalSystem *chemSystemRaw, SysComp::CalculatedProperties *calcProps);
	CalculatorSystemPack(CalculatorConstituentVec &&ccVec, CalculatorIonicFormVec &&ifVec, const SysComp::ChemicalSystem *chemSystemRaw, SysComp::CalculatedProperties *calcProps) noexcept;
	CalculatorSystemPack(const CalculatorSystemPack &other) = delete;
	CalculatorSystemPack(CalculatorSystemPack &&other) noexcept;
	~CalculatorSystemPack();

	CalculatorSystemPack & operator=(const CalculatorSystemPack &other) = delete;
	CalculatorSystemPack & operator=(CalculatorSystemPack &&other) noexcept;

	CalculatorConstituentVec constituents;		/*!< Vector of all constituents in the system. */
	const InternalConstituentVec internalConstituents;
	const CalculatorIonicFormVec ionicForms;	/*!< Vector of all ionic forms in the system. */
	const SysComp::ChemicalSystem *chemSystemRaw;	/*!< Raw pointer to \p SysComp::ChemicalSystem object corresponding to
							     ECHMETCoreLibs representation of the system. */
	SysComp::CalculatedProperties *calcPropsRaw;	/*!< Raw pointer to \p SysComp::CalculatedProperties object corresponding to
							     ECHMETCoreLibs representation of the system.  */
	double conductivity;				/*!< Electric conductivity of the system. This cannot be set by the c-tor
							     because the internal system representation is reusable. */

private:
	InternalConstituentVec makeInternalConstituentVec();
	bool m_moved;					/*!< Used internally by move c-tor. */
};

class DeltaPack {
public:
	explicit DeltaPack();
	DeltaPack(EMVector &&concentrationDeltas, const double conductivityDelta, const SysComp::Constituent *perturbedConstituent) noexcept;

	DeltaPack & operator=(const DeltaPack &other);
	const EMVector concentrationDeltas;			/*!< Vector of first-derivative deltas of ionic form concentrations sorted in the same order
								     as the ionic concentrations in \p CalculatorSystemPack. */
	const double conductivityDelta;				/*!< Delta of the overall system conductivity. */
	const SysComp::Constituent *perturbedConstituent;	/*!< Constituent whose concentration was perturbed to calculate the deltas. */
};
typedef std::vector<DeltaPack> DeltaPackVec;

class QLQRPack {
public:
	QLQRPack(const EMMatrixC &QL, const EMMatrixC &QR);
	QLQRPack(const QLQRPack &other);
	QLQRPack(QLQRPack &&other) noexcept;

	const EMMatrixC & QL() const;
	const EMMatrixC & QR() const;

private:
	std::unique_ptr<EMMatrixC> m_QL;
	std::unique_ptr<EMMatrixC> m_QR;
};

class SolutionProperties {
public:
	SolutionProperties();
	SolutionProperties(const double bufferCapacity, const double conductivity, const double ionicStrength,
			   std::vector<double> &&analyticalConcentration, std::vector<double> &&ionicConcentrations, std::vector<double> &&effectiveMobilities) noexcept;
	SolutionProperties(const SolutionProperties &other);
	SolutionProperties(SolutionProperties &&other) noexcept;

	SolutionProperties & operator=(const SolutionProperties &other);
	SolutionProperties & operator=(SolutionProperties &&other) noexcept;

	const double bufferCapacity;
	const double conductivity;
	const double ionicStrength;
	const std::vector<double> analyticalConcentrations; /* SysComp ordering */
	const std::vector<double> ionicConcentrations;	    /* SysComp ordering */
	const std::vector<double> effectiveMobilities;	    /* SysComp ordering */
};

} // namespace Calculator
} // namespace LEMNG
} // namespace ECHMET

#endif // ECHMET_LEMNG_CALCULATOR_TYPES_H
