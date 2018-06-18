#include <lemng.h>
#include <cmath>
#include <future>
#include <tuple>
#include <vector>

#define USE_ECHMET_CONTAINERS
#include <containers/echmetvec_p.h>

#include "tracing/lemng_tracer_impl.h"
#include <sstream>

namespace ECHMET {
namespace LEMNG {

static const int POINTS_PER_SEC = 40;
static const double TIME_STEP = 1.0 / POINTS_PER_SEC;

typedef std::vector<std::future<std::vector<double>>> HVLPlotFutureVec;

class EigenzonePlotParams {
public:
	EigenzonePlotParams() :
		vZero{0},
		vEMD{0},
		zoneSignal{0},
		diffCoeff{0},
		visible{false}
	{
	}

	EigenzonePlotParams(const double vZero, const double vEMD, const double zoneSignal, const double diffCoeff) :
		vZero{vZero},
		vEMD{vEMD},
		zoneSignal{zoneSignal},
		diffCoeff{diffCoeff},
		visible{true}
	{
	}

	const double vZero;
	const double vEMD;
	const double zoneSignal;
	const double diffCoeff;
	const bool visible;
};

static
double signalResponse(const RSolutionProperties &solProps, const EFGResponseType respType, const char *constituentName)
{
	switch (respType) {
	case EFGResponseType::RESP_CONDUCTIVITY:
		return solProps.conductivity;
	case EFGResponseType::RESP_CONCENTRATION:
	{
		RConstituent ctuent;
		if (solProps.composition->at(ctuent, constituentName) != ::ECHMET::RetCode::OK)
			return 0;

		return ctuent.concentration;
	}
	case EFGResponseType::RESP_PH:
		return solProps.pH;
	}

	throw std::runtime_error{"Unexpected EFG reponse type"};
}

static
int calcSlice(const int NCpus, const int points) noexcept
{
	const double part = NCpus;

	return static_cast<int>(std::ceil(static_cast<double>(points) / part));
}


static
double guessPlotToTime(const double longestZoneMaximumTime) noexcept
{
	if (longestZoneMaximumTime < 60)
		return 60;
	else if (longestZoneMaximumTime > 3600)
		return 3600;

	return longestZoneMaximumTime * 1.1;
}

static
double inputPlotTimeToTime(const double inputTime, const double longestZoneTime) noexcept
{
	if (inputTime > 0.0)
		return inputTime;
	return guessPlotToTime(longestZoneTime);
}

static
void makeEigenzonePlotParams(const REigenzoneVec *eigenzones,
			     const EFGResponseType respType,
			     const char *constituentName,
			     const double E, const double effectiveLength, const double EOFVelocity,
			     std::vector<EigenzonePlotParams> &ezPlotParams,
			     double &longestZoneMaximumTime)
{
	for (size_t idx = 0; idx < eigenzones->size(); idx++) {
		const REigenzone &ez = eigenzones->at(idx);

		if (!ez.valid) {
			ezPlotParams.emplace_back();
			continue;
		}

		const double ezMob = ez.mobility * 1.0e-9;
		const double diffCoeff = ez.a2t * 1.0e-9;
		const double vZero = ezMob * E;
		const double vEMD = ez.uEMD * E * 1.0e-9;
		const double zoneMaximumTime = effectiveLength / (vZero + EOFVelocity);
		double zoneSignal = signalResponse(ez.solutionProperties, respType, constituentName);

		if (zoneMaximumTime < 0.0) {
			/* Invisible zone */
			ezPlotParams.emplace_back();
			continue;
		}

		if (zoneMaximumTime > longestZoneMaximumTime)
			longestZoneMaximumTime = zoneMaximumTime;

		ezPlotParams.emplace_back(vZero, vEMD, zoneSignal, diffCoeff);
	}
}

static
void makePlotBaseline(const double baselineSignal, const int NCpus, const int points, VecImpl<EFGPair, false>::STLVec &stlEfg)
{
	const auto worker = [&](const int from, const int to) noexcept {
		if (from > to) return;

		for (int _idx = from; _idx < to; _idx++) {
			stlEfg[_idx].time = _idx * TIME_STEP;
			stlEfg[_idx].value = baselineSignal;
		}
	};

	const int slice = calcSlice(NCpus, points);

	std::vector<std::thread> workers{};
	workers.reserve(NCpus);

	for (int idx = 0; idx < NCpus - 1; idx++) {
		const int from = slice * idx;
		const int to = from + slice;

		workers.emplace_back(worker, from, to);
	}
	workers.emplace_back(worker, slice * (NCpus - 1), points);

	for (auto &w : workers)
		w.join();
}

static
double calculateHVLR(const double t, const double x, const double d, const double vZero, const double vEMD, const double L) noexcept
{
	static const double ZERO = 0;
	static const double ONE_HALF = 0.5;
	static const double ONE = 1;
	static const double TWO = 2;
	static const double FOUR = 4;
	static const double ERFC_FLAT_THRESHOLD = 25.0;
	static const double LN_PI = std::log(M_PI);

	static const auto lnErfc = [](const double v) {
		const double vSq = std::pow(v, 2);

		const double A = -vSq;
		const double B = -ONE_HALF * LN_PI;
		const double C = -log(v);
		const double D = -ONE / (2.0 * vSq);
		const double E = 5.0 / (8.0 * std::pow(v, 4));

		return A + B + C + D + E;
	};

	static const auto EME = [](const double E, const double b) {
		if (b > ERFC_FLAT_THRESHOLD)
			return E + lnErfc(b);
		else
			return E + log(std::erfc(b));
	};

	static const auto accum = [](double tPlus, double tMinus, const auto &op) {
		if (tPlus > tMinus)
			std::swap(tPlus, tMinus);

		return tMinus + log(op(ONE, std::exp(tPlus - tMinus)));
	};

	const double den = std::sqrt(FOUR * d * t);
	const double LHalf = L / TWO;
	const double xMvZt = x - vZero * t;
	const double vEMDt = vEMD * t;
	const double twoD = TWO * d;
	const double posMinus = (xMvZt - vEMDt - LHalf);
	const double posPlus = (xMvZt - vEMDt + LHalf);

	double aMinus = posMinus / den;
	double aPlus = posPlus / den;

	if (std::abs(aPlus - aMinus) < 1.0e-13)
		return ZERO;

	if (aPlus < ZERO && aMinus < ZERO) {
		/* Invert the problem to "positive" section */
		const double _t = -aPlus;
		aPlus = -aMinus;
		aMinus = _t;
	}

	const double lnRV = [&]() {
		if (aPlus > ERFC_FLAT_THRESHOLD && aMinus > ERFC_FLAT_THRESHOLD) {
			/* We are in an area where the error function rises so slowly
			 * that the standard double precision cannot represent the delta precisely enough. */

			const double lnErfc_aPlus = EME(ZERO, aPlus);
			const double lnErfc_aMinus = EME(ZERO, aMinus);

			return accum(lnErfc_aPlus, lnErfc_aMinus, std::minus<double>{});
		} else
			return log(std::erfc(aMinus) - std::erfc(aPlus));
	}();

	const double EMinus = vEMD / twoD * (xMvZt - 0.5 * vEMDt - LHalf);
	const double EPlus = vEMD / twoD *  (xMvZt - 0.5 * vEMDt + LHalf);
	const double bMinus = -(xMvZt - LHalf) / den;
	const double bPlus = (xMvZt + LHalf) / den;

	const double EMinusErfc = EME(EMinus, bMinus);
	const double EPlusErfc = EME(EPlus, bPlus);

	const double lnQ = accum(EMinusErfc, EPlusErfc, std::plus<double>{});

	const double F = lnQ - lnRV;

	if (F < log(std::numeric_limits<double>::max()) - 2)
		return ONE / (ONE + std::exp(F));

	return ZERO;
}

static
REigenzoneEnvelope calcZoneEnvelope(const EigenzonePlotParams &params, const double E, const double vEOF, const double effectiveLength, const double zoneLength, const double tLast)
{
	if (!params.visible)
		return { -1.0, -1.0, 0.0 };

	const double THRESHOLD = 0.05;
	double zoneTime = effectiveLength / (params.vZero + vEOF);
	double beginsAt;
	double endsAt;

	/* Find maximum value */
	double t = effectiveLength / (params.vZero + vEOF);
	double yMax = calculateHVLR(t, effectiveLength - vEOF * t, params.diffCoeff, params.vZero, params.vEMD, zoneLength);
	double yX = yMax;

	if (t > tLast)
		return { -1.0, -1.0, 0.0 };

	const double uEMDabs = std::abs(params.vEMD / E);

	if (uEMDabs > 1.0e-13) {
		if (params.vEMD > 0.0) {
			do  {
				t -= TIME_STEP;
				yX = calculateHVLR(t, effectiveLength - vEOF * t, params.diffCoeff, params.vZero, params.vEMD, zoneLength);
				if (yX > yMax)
					yMax = yX;
				else
					break;
			} while (t > 0.0);
		} else if (params.vEMD < -0.0) {
			do {
				t += TIME_STEP;
				yX = calculateHVLR(t, effectiveLength - vEOF * t, params.diffCoeff, params.vZero, params.vEMD, zoneLength);
				if (yX > yMax)
					yMax = yX;
				else
					break;
			} while (t < tLast);
		}

		yX = yMax;
		zoneTime = t;
	}

	/* Find the envelope of the zone */
	/* Left lobe */
	t = zoneTime - TIME_STEP;
	while (yX / yMax > THRESHOLD && t > 0.0) {
		yX = calculateHVLR(t, effectiveLength - vEOF * t, params.diffCoeff, params.vZero, params.vEMD, zoneLength);
		t -= TIME_STEP;
	}
	beginsAt = t;

	/* Right lobe */
	t = zoneTime + TIME_STEP;
	yX = yMax;
	while (yX / yMax > THRESHOLD && t < tLast) {
		yX = calculateHVLR(t, effectiveLength - vEOF * t, params.diffCoeff, params.vZero, params.vEMD, zoneLength);
		t += TIME_STEP;
	}
	endsAt = t;

	return REigenzoneEnvelope{ beginsAt, endsAt, yMax };
}

static
void makePlot(const std::vector<EigenzonePlotParams> &ezPlotParams, const double effectiveLength, const double bslSignal,
	      const double plotToTime, const double zoneLength, const double vEOF,
	      VecImpl<EFGPair, false>::STLVec &stlEfg)
{
	const int NCpus = [](){
		const int n = std::thread::hardware_concurrency();
		if (n < 1)
			return 1;
		return n;
	}();

	const int points = plotToTime * POINTS_PER_SEC;

	const auto worker = [bslSignal, effectiveLength, zoneLength, vEOF, &stlEfg](const EigenzonePlotParams &params, const double amplitude, const int from, const int to) noexcept {
		if (from > to) return;

		for (int _idx = from; _idx < to; _idx++) {
			const double t = _idx * TIME_STEP;
			const double actualEffectiveLength = effectiveLength - vEOF * t;
			const double HVLRy = calculateHVLR(t, actualEffectiveLength, params.diffCoeff, params.vZero, params.vEMD, zoneLength);

			auto &p = stlEfg[_idx];
			p.value += HVLRy * amplitude;
		}
	};


	stlEfg.resize(points);

	makePlotBaseline(bslSignal, NCpus, points, stlEfg);

	const int slice = calcSlice(NCpus, points);
	for (const auto &params : ezPlotParams) {
		if (!params.visible)
			continue;

		std::vector<std::thread> workers{};
		workers.reserve(NCpus);

		const double amplitude = (params.zoneSignal - bslSignal);

		for (int idx = 0; idx < NCpus - 1; idx++) {
			const int from = (idx * slice) + ((idx == 0) ? 1 : 0);
			const int to = from + slice - ((idx == 0) ? 1 : 0);

			workers.emplace_back(worker, params, amplitude, from, to);
		}
		workers.emplace_back(worker, params, amplitude, slice * (NCpus - 1), points);

		for (auto &w : workers)
			w.join();
	}
}

RetCode ECHMET_CC findEigenzoneEnvelopes(REigenzoneEnvelopeVec *&envelopes, const Results &results,
					 const double drivingVoltage, const double totalLength, const double effectiveLength,
					 const double EOFMobility, const double injectionZoneLength, const double plotToTime) noexcept
{
	if (totalLength <= 0)
		return RetCode::E_INVALID_CAPILLARY;
	if (totalLength < effectiveLength || effectiveLength <= 0)
		return RetCode::E_INVALID_DETECTOR_POSITION;
	if (injectionZoneLength <= 0.0)
		return RetCode::E_INVALID_ARGUMENT;

	const double E = drivingVoltage / totalLength;		/* Electric field intensity */
	const double EOFVelocity = EOFMobility * E * 1.0e-9;

	envelopes = createECHMETVec<REigenzoneEnvelope, false>(results.eigenzones->size());
	if (envelopes == nullptr)
		return RetCode::E_NO_MEMORY;

	try {
		double longestZoneTime = 0.0;
		std::vector<EigenzonePlotParams> ezPlotParams{};

		makeEigenzonePlotParams(results.eigenzones,
					EFGResponseType::RESP_CONDUCTIVITY,	/* The response type does not matter here */
					nullptr,
					E, effectiveLength, EOFVelocity,
					ezPlotParams,
					longestZoneTime);

		const double _plotToTime = inputPlotTimeToTime(plotToTime, longestZoneTime);

		for (const auto &params : ezPlotParams) {
			const auto envelope = calcZoneEnvelope(params, E, EOFVelocity, effectiveLength, injectionZoneLength, _plotToTime);
			ECHMET_TRACE(LEMNGTracing, EFGPLOT_ZONE_ENVELOPE, params.vZero, envelope.beginsAt, envelope.endsAt);

			if (envelopes->push_back(envelope) != ::ECHMET::RetCode::OK) {
				envelopes->destroy();

				return RetCode::E_NO_MEMORY;
			}
		}
	} catch (std::bad_alloc &) {
		envelopes->destroy();

		return RetCode::E_NO_MEMORY;
	}

	return RetCode::OK;
}

RetCode ECHMET_CC plotElectrophoregram(EFGPairVec *&electrophoregram,
				       const Results &results,
				       const double drivingVoltage, const double totalLength, const double effectiveLength,
				       const double EOFMobility,
				       const double injectionZoneLength,
				       const EFGResponseType respType,
				       const char *constituentName,
				       const double plotToTime) noexcept
{
	ECHMET_TRACE(LEMNGTracing, EFGPLOT_INPUT_PARAMS, drivingVoltage, totalLength, effectiveLength, EOFMobility, injectionZoneLength, respType, constituentName, plotToTime);

	if (totalLength <= 0)
		return RetCode::E_INVALID_CAPILLARY;
	if (totalLength < effectiveLength || effectiveLength <= 0)
		return RetCode::E_INVALID_DETECTOR_POSITION;
	if (respType == EFGResponseType::RESP_CONCENTRATION && constituentName == nullptr)
		return RetCode::E_INVALID_ARGUMENT;
	if (injectionZoneLength <= 0.0)
		return RetCode::E_INVALID_ARGUMENT;

	const double E = drivingVoltage / totalLength;	/* Electric field intensity */
	const double EOFVelocity = EOFMobility * E * 1.0e-9;

	VecImpl<EFGPair, false> *electrophoregramImpl = createECHMETVec<EFGPair, false>(0);
	if (electrophoregramImpl == nullptr)
		return RetCode::E_NO_MEMORY;

	try {
		double longestZoneTime = 0.0;
		const double bslSignal = signalResponse(results.BGEProperties, respType, constituentName);
		std::vector<EigenzonePlotParams> ezPlotParams{};

		makeEigenzonePlotParams(results.eigenzones,
					respType,
					constituentName,
					E, effectiveLength, EOFVelocity,
					ezPlotParams,
					longestZoneTime);

		const double _plotToTime = inputPlotTimeToTime(plotToTime, longestZoneTime);

		makePlot(ezPlotParams, effectiveLength, bslSignal, _plotToTime, injectionZoneLength, EOFVelocity, electrophoregramImpl->STL());
	} catch (std::bad_alloc &) {
		return RetCode::E_NO_MEMORY;
	} catch (std::runtime_error &) {
		return RetCode::E_INTERNAL_ERROR;
	}

	electrophoregram = electrophoregramImpl;

	return RetCode::OK;
}

} // namespace LEMNG

#ifndef ECHMET_TRACER_DISABLE_TRACING

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, EFGPLOT_INPUT_PARAMS, "Input parameters for EFG plotter")
ECHMET_MAKE_LOGGER(LEMNGTracing, EFGPLOT_INPUT_PARAMS, const double voltage, const double totalLength, const double effectiveLength,
						       const double EOFMobility, const double injectionZoneLength,
						       const ECHMET::LEMNG::EFGResponseType respType, const char *constituentName,
						       const double plotToTime)
{
	std::ostringstream ss{};

	auto addParam = [&](const std::string &param, const double v) {
		ss << param << " = " << v << "\n";
	};

	ss << "EFG plotter input parameters\n";

	addParam("Voltage", voltage);
	addParam("Total length", totalLength);
	addParam("Effective length", effectiveLength);
	addParam("EOF mobility", EOFMobility);
	addParam("Injection zone length", injectionZoneLength);

	ss << "Response type = " << [](ECHMET::LEMNG::EFGResponseType rt) {
		switch (rt) {
		case ECHMET::LEMNG::EFGResponseType::RESP_CONDUCTIVITY:
			return "Conductivity";
		case ECHMET::LEMNG::EFGResponseType::RESP_PH:
			return "pH";
		case ECHMET::LEMNG::EFGResponseType::RESP_CONCENTRATION:
			return "Concentration";
		}
		return "";
	}(respType) << "\n";

	ss << "Constituent name = ";
	if (constituentName == nullptr)
		ss << "<nullptr>";
	else
		ss << constituentName;
	ss << "\n";

	addParam("Plot to time", plotToTime);

	return ss.str();
}

ECHMET_MAKE_TRACEPOINT(LEMNGTracing, EFGPLOT_ZONE_ENVELOPE, "Eigenzone envelopes")
ECHMET_MAKE_LOGGER(LEMNGTracing, EFGPLOT_ZONE_ENVELOPE, const double mobility, const double begins, const double ends)
{
	std::ostringstream ss{};

	ss << "Eigenzone envelope: u = " << mobility << ", beginsAt = " << begins << ", endsAt = " << ends << "\n";

	return ss.str();
}

#endif // ECHMET_TRACER_DISABLE_TRACING

} // namespace ECHMET
