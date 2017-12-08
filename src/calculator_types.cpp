#include "calculator_types.h"
#include "tracing/lemng_tracer_impl.h"

namespace ECHMET {
namespace LEMNG {
namespace Calculator {

CalculatorIonicForm::CalculatorIonicForm(const std::string &name, const int32_t charge,
	const SysComp::IonicForm *internalIonicForm,
	const size_t internalIonicFormConcentrationIdx,
	const size_t globalIonicFormConcentrationIdx,
	const std::vector<size_t> &containedConstituents, const bool isAnalyte) :
	name{ name },
	charge{ charge },
	internalIonicForm{ internalIonicForm },
	internalIonicFormConcentrationIdx{ internalIonicFormConcentrationIdx },
	globalIonicFormConcentrationIdx{ globalIonicFormConcentrationIdx },
	containedConstituents(containedConstituents),
	isAnalyte{ isAnalyte },
	concentration{ -1 },
	mobility{ -1 }
{
}

CalculatorIonicForm::CalculatorIonicForm(std::string &&name, const int32_t charge,
	const SysComp::IonicForm *internalIonicForm,
	const size_t internalIonicFormConcentrationIdx,
	const size_t globalIonicFormConcentrationIdx,
	std::vector<size_t> &&containedConstituents, const bool isAnalyte) noexcept :
	name{ name },
	charge{ charge },
	internalIonicForm{ internalIonicForm },
	internalIonicFormConcentrationIdx{ internalIonicFormConcentrationIdx },
	globalIonicFormConcentrationIdx{ globalIonicFormConcentrationIdx },
	containedConstituents(containedConstituents),
	isAnalyte{ isAnalyte },
	concentration{ -1 },
	mobility{ -1 }
{
}

CalculatorIonicForm & CalculatorIonicForm::operator=(const CalculatorIonicForm &other)
{
	const_cast<std::string&>(name) = other.name;
	const_cast<int32_t&>(charge) = other.charge;
	internalIonicForm = other.internalIonicForm;
	const_cast<size_t&>(internalIonicFormConcentrationIdx) = other.internalIonicFormConcentrationIdx;
	const_cast<size_t&>(globalIonicFormConcentrationIdx) = other.globalIonicFormConcentrationIdx;
	const_cast<std::vector<size_t>&>(containedConstituents) = other.containedConstituents;
	const_cast<bool&>(isAnalyte) = other.isAnalyte;
	mobility = other.mobility;
	concentration = other.concentration;

	return *this;
}

CalculatorIonicForm & CalculatorIonicForm::operator=(CalculatorIonicForm &&other) noexcept
{
	const_cast<std::string&>(name) = std::move(other.name);
	const_cast<int32_t&>(charge) = other.charge;
	internalIonicForm = other.internalIonicForm;
	const_cast<size_t&>(internalIonicFormConcentrationIdx) = other.internalIonicFormConcentrationIdx;
	const_cast<size_t&>(globalIonicFormConcentrationIdx) = other.globalIonicFormConcentrationIdx;
	const_cast<std::vector<size_t>&>(containedConstituents) = std::move(other.containedConstituents);
	const_cast<bool&>(isAnalyte) = other.isAnalyte;
	mobility = other.mobility;
	concentration = other.concentration;

	return *this;
}

CalculatorConstituent::CalculatorConstituent(const std::string &name, const CalculatorIonicFormVec &ifVec, const SysComp::Constituent *internalConstituent, const bool isAnalyte) :
	name{ name },
	ionicForms(ifVec),
	internalConstituent{ internalConstituent },
	isAnalyte{ isAnalyte },
	concentrationBGE{ -1 },
	concentrationSample{ -1 }
{
}

CalculatorConstituent::CalculatorConstituent(std::string &&name, CalculatorIonicFormVec &&ifVec, const SysComp::Constituent *internalConstituent, const bool isAnalyte) noexcept :
	name{ name },
	ionicForms(ifVec),
	internalConstituent{ internalConstituent },
	isAnalyte{ isAnalyte },
	concentrationBGE{ -1 },
	concentrationSample{ -1 }
{
}

CalculatorConstituent::CalculatorConstituent(const CalculatorConstituent &other) :
	name{ other.name },
	ionicForms(other.ionicForms),
	internalConstituent{ other.internalConstituent },
	isAnalyte{ other.isAnalyte },
	concentrationBGE{ other.concentrationBGE },
	concentrationSample{ other.concentrationSample }
{
}

CalculatorConstituent::CalculatorConstituent(CalculatorConstituent &&other) noexcept :
	name{ std::move(other.name) },
	ionicForms(std::move(other.ionicForms)),
	internalConstituent{ other.internalConstituent },
	isAnalyte{ other.isAnalyte },
	concentrationBGE{ other.concentrationBGE },
	concentrationSample{ other.concentrationSample }
{
}

CalculatorSystemPack::~CalculatorSystemPack()
{
	if (!m_moved) {
		for (auto && iF : ionicForms)
			delete iF;
	}
}

CalculatorConstituent & CalculatorConstituent::operator=(const CalculatorConstituent &other)
{
	const_cast<std::string&>(name) = other.name;
	const_cast<CalculatorIonicFormVec&>(ionicForms) = other.ionicForms;
	const_cast<const SysComp::Constituent *&>(internalConstituent) = other.internalConstituent;
	const_cast<bool&>(isAnalyte) = other.isAnalyte;
	concentrationBGE = other.concentrationBGE;
	concentrationSample = other.concentrationSample;

	return *this;
}

CalculatorConstituent & CalculatorConstituent::operator=(CalculatorConstituent &&other) noexcept
{
	const_cast<std::string&>(name) = std::move(other.name);
	const_cast<CalculatorIonicFormVec&>(ionicForms) = std::move(other.ionicForms);
	const_cast<const SysComp::Constituent *&>(internalConstituent) = other.internalConstituent;
	const_cast<bool&>(isAnalyte) = other.isAnalyte;
	concentrationBGE = other.concentrationBGE;
	concentrationSample = other.concentrationSample;

	return *this;
}

CalculatorSystemPack::CalculatorSystemPack() :
	constituents{},
	ionicForms{},
	conductivity{ -1 },
	m_moved{ false }
{
}

CalculatorSystemPack::CalculatorSystemPack(const CalculatorConstituentVec &ccVec, const CalculatorIonicFormVec &ifVec, const SysComp::ChemicalSystem *chemSystemRaw, SysComp::CalculatedProperties *calcPropsRaw) :
	constituents(ccVec),
	ionicForms(ifVec),
	chemSystemRaw{ chemSystemRaw },
	calcPropsRaw{ calcPropsRaw },
	conductivity{ -1 },
	m_moved{ false }
{
}

CalculatorSystemPack::CalculatorSystemPack(CalculatorConstituentVec &&ccVec, CalculatorIonicFormVec &&ifVec, const SysComp::ChemicalSystem *chemSystemRaw, SysComp::CalculatedProperties *calcPropsRaw) noexcept:
	constituents(ccVec),
	ionicForms(ifVec),
	chemSystemRaw{ chemSystemRaw },
	calcPropsRaw{ calcPropsRaw },
	conductivity{ -1 },
	m_moved{ false }
{
}

CalculatorSystemPack::CalculatorSystemPack(CalculatorSystemPack &&other) noexcept :
constituents(std::move(other.constituents)),
	ionicForms(std::move(other.ionicForms)),
	chemSystemRaw{ other.chemSystemRaw },
	calcPropsRaw{ other.calcPropsRaw },
	conductivity{ other.conductivity },
	m_moved{ false }
{
	other.m_moved = true;
}

CalculatorSystemPack & CalculatorSystemPack::operator=(CalculatorSystemPack &&other) noexcept
{
	constituents = std::move(other.constituents);
	const_cast<CalculatorIonicFormVec&>(ionicForms) = std::move(other.ionicForms);
	chemSystemRaw = other.chemSystemRaw;
	calcPropsRaw = other.calcPropsRaw;
	conductivity = other.conductivity;

	other.m_moved = true;

	return *this;
}

DeltaPack::DeltaPack() :
	concentrationDeltas{EMVector{0}},
	conductivityDelta{0},
	perturbedConstituent{nullptr}
{
}

DeltaPack::DeltaPack(EMVector &&concentrationDeltas, const double conductivityDelta, const SysComp::Constituent *perturbedConstituent) noexcept :
	concentrationDeltas{concentrationDeltas},
	conductivityDelta{conductivityDelta},
	perturbedConstituent{perturbedConstituent}
{
}

DeltaPack & DeltaPack::operator=(const DeltaPack &other)
{
	const_cast<EMVector&>(concentrationDeltas) = other.concentrationDeltas;
	const_cast<double&>(conductivityDelta) = other.conductivityDelta;
	perturbedConstituent = other.perturbedConstituent;

	return *this;
}

QLQRPack::QLQRPack(const EMMatrixC &QL, const EMMatrixC &QR)
{
	m_QL = std::unique_ptr<EMMatrixC>(new EMMatrixC{QL});
	m_QR = std::unique_ptr<EMMatrixC>(new EMMatrixC{QR});
}

QLQRPack::QLQRPack(const QLQRPack &other)
{
	m_QL = std::unique_ptr<EMMatrixC>(new EMMatrixC{*other.m_QL});
	m_QR = std::unique_ptr<EMMatrixC>(new EMMatrixC{*other.m_QR});
}

QLQRPack::QLQRPack(QLQRPack &&other) noexcept :
	m_QL{std::move(other.m_QL)},
	m_QR{std::move(other.m_QR)}
{
}

const EMMatrixC & QLQRPack::QL() const
{
	return *m_QL;
}

const EMMatrixC & QLQRPack::QR() const
{
	return *m_QR;
}

SolutionProperties::SolutionProperties() :
	bufferCapacity{-1},
	conductivity{-1},
	ionicStrength{-1},
	analyticalConcentrations{},
	ionicConcentrations{}
{
}

SolutionProperties::SolutionProperties(const double bufferCapacity, const double conductivity, const double ionicStrength,
				       std::vector<double> &&analyticalConcentrations, std::vector<double> &&ionicConcentrations) noexcept :
	bufferCapacity{bufferCapacity},
	conductivity{conductivity},
	ionicStrength{ionicStrength},
	analyticalConcentrations(analyticalConcentrations),
	ionicConcentrations(ionicConcentrations)
{
}

SolutionProperties::SolutionProperties(const SolutionProperties &other) :
	bufferCapacity{other.bufferCapacity},
	conductivity{other.conductivity},
	ionicStrength{other.ionicStrength},
	analyticalConcentrations(other.analyticalConcentrations),
	ionicConcentrations(other.ionicConcentrations)
{
	_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::CALC_OBJECT_CONSTRUCTION, const std::string&>(ECHMET_S("SolutionProperties copy c-tor"));
}

SolutionProperties::SolutionProperties(SolutionProperties &&other) noexcept :
	bufferCapacity{other.bufferCapacity},
	conductivity{other.conductivity},
	ionicStrength{other.ionicStrength},
	analyticalConcentrations(std::move(other.analyticalConcentrations)),
	ionicConcentrations(std::move(other.ionicConcentrations))
{
	_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::CALC_OBJECT_CONSTRUCTION, const std::string&>(ECHMET_S("SolutionProperties move c-tor"));
}

SolutionProperties & SolutionProperties::operator=(const SolutionProperties &other)
{
	const_cast<double&>(bufferCapacity) = other.bufferCapacity;
	const_cast<double&>(conductivity) = other.conductivity;
	const_cast<double&>(ionicStrength) = other.ionicStrength;
	const_cast<std::vector<double>&>(analyticalConcentrations) = other.analyticalConcentrations;
	const_cast<std::vector<double>&>(ionicConcentrations) = other.ionicConcentrations;

	_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::CALC_OBJECT_CONSTRUCTION, const std::string&>(ECHMET_S("SolutionProperties copy assignment"));

	return *this;
}

SolutionProperties & SolutionProperties::operator=(SolutionProperties &&other) noexcept
{
	const_cast<double&>(bufferCapacity) = other.bufferCapacity;
	const_cast<double&>(conductivity) = other.conductivity;
	const_cast<double&>(ionicStrength) = other.ionicStrength;
	const_cast<std::vector<double>&>(analyticalConcentrations) = std::move(other.analyticalConcentrations);
	const_cast<std::vector<double>&>(ionicConcentrations) = std::move(other.ionicConcentrations);

	_ECHMET_TRACE<LEMNGTracing, LEMNGTracing::CALC_OBJECT_CONSTRUCTION, const std::string&>(ECHMET_S("SolutionProperties move assignment"));

	return *this;
}

} // namespace Calculator
} // namespace LEMNG

#ifndef ECHMET_TRACER_DISABLE_TRACING

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, CALC_OBJECT_CONSTRUCTION, "Calculator object construction/assignment")
ECHMET_MAKE_LOGGER(LEMNGTracing, CALC_OBJECT_CONSTRUCTION, const std::string &msg)
{
	return msg;
}

#endif // ECHMET_TRACER_DISABLE_TRACING

} // namespace ECHMET
