#ifndef ECHMET_LEMNG_LEMNG_H
#define ECHMET_LEMNG_LEMNG_H

#define ECHMET_IMPORT_INTERNAL
#include <echmetelems.h>
#include <echmetsyscomp.h>
#undef ECHMET_IMPORT_INTERNAL
#include <echmetmodule.h>

namespace ECHMET {

/*!
 * Package that implements linear model of electromigration.
 */
namespace LEMNG {


/*!
 * Possible return codes from public API calls.
 */
ECHMET_ST_ENUM(RetCode) {
	OK = 0,
	E_NO_MEMORY = 0x1,			/*!< Insufficient memory to complete operation. */
	E_INVALID_ARGUMENT = 0x2,		/*!< Argument passed to a function was invalid. */
	E_INVALID_CAPILLARY = 0x3,		/*!< Invalid capillary length. */
	E_INVALID_DETECTOR_POSITION = 0x4,	/*!< Position of the detector is invalid. */
	E_NOT_IMPLEMENTED = 0x5,		/*!< Requested function or operation is not implemented */
	E_CANNOT_SOLVE_BGE = 0x6,		/*!< Cannot calculate equilibrium composition of background electrolyte or its properties */
	E_DATA_TOO_LARGE = 0x7,			/*!< System is too large to process */
	E_INVALID_CONSTITUENT = 0x8,		/*!< System contains a constituent with invalid properties */
	E_INVALID_COMPLEXATION = 0x9,		/*!< System contains nonsensical complexation relation */
	E_DUPLICIT_CONSTITUENTS = 0x10,		/*!< System contains multiple constituents with the same name */
	E_UNKW_CORELIBS_ERROR = 0x11,		/*!< Unrecognized ECHMETCoreLibs error */
	E_CHEM_SYSTEM_UNSOLVABLE = 0x12,	/*!< Chemical system equilibrium cannot be solved */
	E_INTERNAL_ERROR = 0x13,		/*!< Unspecified internal error */
	E_COMPLEX_EIGENMOBILITIES = 0x14,	/*!< System contains complex eigenmobilites and therefore exhibits oscillating behavior
						     see Hruška, V; Jaroš, M; Gaš, B, ELECTROPHORESIS 2006 Volume: 27  Issue: 3  Pages: 513-518  Special Issue: SI (DOI: 10.1002/elps.200500731) */
	E_CONCENTRATION_TOO_LOW = 0x15,		/*!< Concentration of a constituent is too low to ensure that
						     the numerical sovler will be able to solve the system */
	E_PARTIAL_EIGENZONES = 0x16		/*!< Some eigenzones in the system could not have been fully resolved */
	ENUM_FORCE_INT32_SIZE(LEMNGRetCode)
};

/*!
 * Possible types of eigenzone
 */
ECHMET_ST_ENUM(EigenzoneType) {
	INVALID = 0x0,	/*!< Should not be returned. */
	ANALYTE = 0x1,	/*!< Zone corresponds to an analyte. */
	SYSTEM = 0x2	/*!< Zone does not correspond to any analyte - a system zone. */
	ENUM_FORCE_INT32_SIZE(LEMNGEigenzoneType)
};

ECHMET_ST_ENUM(EFGResponseType) {
	RESP_CONDUCTIVITY = 0x0,	/*!< Plot conductivity response */
	RESP_CONCENTRATION = 0x1,	/*!< Plot concentration of a given constituent */
	RESP_PH = 0x2			/*!< Plot pH response */
	ENUM_FORCE_INT32_SIZE(EFGType)
};

/*!
 * Description of a tracepoint.
 */
class TracepointInfo {
public:
	int32_t id;			/*!< Internal ID of the tracepoint. Used to
					     set the tracepoint state. */
	FixedString *description;	/*!< Human-readable description of the tracepoint. */
};
IS_POD(TracepointInfo)
typedef Vec<TracepointInfo> TracepointInfoVec;

/*!
 * Description of a dissociation ratio for a given ionic form
 */
class RDissociationRatio {
public:
	FixedString *name;	/*!< Name of the ionic compound */
	double fraction;	/*!< Molar fraction of the ionic form */
};
IS_POD(RDissociationRatio)
typedef Vec<RDissociationRatio> RDissociationRatioVec;

/*!
 * Description of molar fractions of all ionic forms of a dissociated component
 */
class RDissociatedConstituent {
public:
	FixedString *name;		/*!< Name of the constituent */
	double effectiveMobility;	/*!< Effective mobility of the constituent */
	RDissociationRatioVec *ratios;	/*!< Dissociation descriptors of all ionic form of the constituent */
};
IS_POD(RDissociatedConstituent)
typedef Vec<RDissociatedConstituent> RDissociatedConstituentVec;

/*!
 * Description of an ion contained in a chemical compound.
 */
class RIon {
public:
	FixedString *name;	/*!< Name of the chemical element constituting the ion. */
	int32_t charge;		/*!< Charge of the ion. */
	int32_t count;		/*!< Number of the given ion present in a compound. */
};
IS_POD(RIon)
typedef Vec<RIon> RIonVec;

/*!
 * Description of a chemical compound (a form).
 */
class RForm {
public:
	int32_t totalCharge;	/*!< Total electric charge of the compound. */
	double concentration;	/*!< Equilibrum concentration in <tt>mmol/dm<sup>3</sup></tt> of the compound. */
	RIonVec *ions;		/*!< Individual ions that make up the compound. */
};
IS_POD(RForm)
typedef SKMap<RForm> RFormMap;

/*!
 * Description of a constituent.
 */
class RConstituent {
public:
	FixedString *name;		/*!< Name of the chemical element constituting the constituent. */
	double concentration;		/*!< Analytical (total) concentration of the constituent in the system <tt>mmol/dm<sup>3</sup></tt>. */
	double effectiveMobility;	/*!< Effective mobility of the constituent. */
	RFormMap *forms;		/*!< All forms that contain the constituent present in the system. */
};
IS_POD(RConstituent)
typedef SKMap<RConstituent> RConstituentMap;

/*!
 * Description of solution properties.
 * This can describe either propeties of the plain background electrolyte
 * of local properties of a solution in an eigenzone.
 */
class RSolutionProperties {
public:
	double pH;			/*!< pH of the solution. This is calculated either from
					     concentration of <tt>H<sub>3</sub>P<sup>+</sup> ions
					     or its activity, respectively, depending on whether the
					     correction for ionic strength was requested. */
	double conductivity;		/*!< Conductivity of the solution in <tt>S/m</tt>. */
	double bufferCapacity;		/*!< Buffering capacity of the solution - currently unimplemented. */
	double ionicStrength;		/*!< Ionic strength of the solution. */
	RConstituentMap *composition;	/*!< Chemical composition of the solution. */
};
IS_POD(RSolutionProperties)

/*!
 * Description of eigenzone.
 */
class REigenzone {
public:
	EigenzoneType ztype;			/*!< Denotes whether a zone belongs to an analyte or not (a system zone). */
	double mobility;			/*!< Electroforetic mobility of the zone */
	double a2t;				/*!< Time-independent diffusive parameter of the zone. */
	double uEMD;				/*!< Measure of electromigration of dispersion of the zone in mobility units. */
	RSolutionProperties solutionProperties;	/*!< Properties of the solution comprising the zone */
	bool tainted;				/*!< Set to true if the concentrations of the constituents
						     that make up the zones had to be clamped to valid values. */
	bool valid;				/*!< Set to false if the eigenzone could not have been fully resolved
						     by the solver. */
};
IS_POD(REigenzone)
typedef Vec<REigenzone> REigenzoneVec;

/*!
 * Eigenzone envelope.
 *
 * Contains time in seconds where an eigenzone begins and ends
 * for a given set of sytem parameters.
 * Beginning and end of the zone is calculated as a point on the time
 * axis there the zone has less than 5 % of its maximum height.
 */
class REigenzoneEnvelope {
public:
	double beginsAt;			/*!< Beginning of the zone */
	double endsAt;				/*!< End of the zone */
	double HVLRMax;				/*!< Peak value of the HVL-R function */
	double tMax;				/*!< Time of maximum signal value */
};
IS_POD(REigenzoneEnvelope)
typedef Vec<REigenzoneEnvelope> REigenzoneEnvelopeVec;

/*!
 * Description of a fully resolved system
 */
class Results {
public:
	RSolutionProperties BGEProperties;			/*!< Properties of the plain background electrolyte. */
	REigenzoneVec *eigenzones;				/*!< Description of all eigenzones present in the system. */
	RDissociatedConstituentVec *analytesDissociation;	/*!< Description of dissociation degrees of all analytes. */
	bool isBGEValid;					/*!< Set to true if the BGE composition was successfully solved. */
};
IS_POD(Results)
typedef MutSKMap<double> InAnalyticalConcentrationsMap;

/*!
 * Time-value data pair.
 * Vector of these composes the expected detector trace.
 */
class EFGPair {
public:
	double time;
	double value;
};
IS_POD(EFGPair)
typedef Vec<EFGPair> EFGPairVec;

/*!
 * Object representing the CZE system to be solved.
 */
class CZESystem {
public:
	/*!
	 * Solves the system.
	 * If the operation does not complete successfully, it it possible to call
	 * \p lastErrorString() to get more detailed information about the reason of failure.
	 *
	 * @param[in] acBGE Analytical concentrations of constituents in plain background electrolyte.
	 * @param[in] acFull Analytical concentrations of constituents in the sample zone.
	 * @param[in] correctForIonicStrength If <tt>true</tt>, correction for ionic strength will be employed in the calculations.
	 * @param[out] results Results will be stored here if the calculation succeeds.
	 *
	 * @retval RetCode::OK Success.
	 * @retval RetCode::E_NO_MEMORY Insufficient memory to perform the calculations.
	 * @retval RetCode::E_INVALID_ARGUMENT Invalid argument was passed to the function.
	 * @retval RetCode::E_INTERNAL_ERROR Unspecifed failure of internal logic occured or the function
	 *                                   the function was passed invalid data.
	 * @retval RetCode::E_CANNOT_SOLVE_BGE Properties of the background electrolyte could not have been solved.
	 * @retval RetCode::E_CHEM_SYSTEM_UNSOLVABLE Equilibrium composition of some of the eigenzones
	 *					     could not have been calculated.
	 * @retval RetCode::E_UNKW_CORELIBS_ERROR Unspecified error returned from the CoreLibs.
	 * @retval RetCode::E_COMPLEX_EIGENMOBILITIES Linear solver detected complex values of eigenmobilities.
	 *                                            Such a system is likely to exhibit oscillating behavior
	 *                                            and cannot be solved by linear theory of electromigration.
	 */
	virtual RetCode ECHMET_CC evaluate(const InAnalyticalConcentrationsMap *acBGE, const InAnalyticalConcentrationsMap *acFull,
					   const NonidealityCorrections corrections, Results &results) ECHMET_NOEXCEPT = 0;

	/*!
	 * Returns human-readable description of the last error that occured
	 * during an attempt to solve the system.
	 *
	 * @retval Human-readable error description.
	 */
	virtual const char * ECHMET_CC lastErrorString() const ECHMET_NOEXCEPT = 0;

	/*
	 * Initializes data structures that can contain the input analytical concentrations.
	 *
	 * @param[out] acMapBGE Concentration map of plain background electrolyte composition.
	 * @param[out] acMapFull Concentration map of sample zone composition.
	 *
	 * @retval RetCode::E_NO_MEMORY Insufficient memory to initialize the maps.
	 */
	virtual RetCode ECHMET_CC makeAnalyticalConcentrationsMaps(InAnalyticalConcentrationsMap *&acMapBGE, InAnalyticalConcentrationsMap *&acMapFull) const ECHMET_NOEXCEPT = 0;

protected:
	virtual ~CZESystem() ECHMET_NOEXCEPT = 0;
};

extern "C" {

ECHMET_API const char * ECHMET_CC LEMNGerrorToString(const RetCode tRet) ECHMET_NOEXCEPT;

/*!
 * Initializes a CZESystem that can solve a system with a given composition.
 *
 * @param[in] BGE Vector of constituents composing the background electrolyte in <tt>SysComp::InConstituentVec</tt> format.
 * @param[in] sample Vector of constituents composing the sample zone in <tt>SysComp::InConstituentVec</tt> format.
 *                   Note that the sample zone shall contain all constituents of the BGE and - optionally - additional
 *                   constituents as analytes.
 * @param[out] czeSystem The CZESystem to be initialized.
 *
 * @retval RetCode::OK Success.
 * @retval RetCode::E_NO_MEMORY
 * @retval Anything that can be returned by <tt>SysComp::makeComposition()</tt> mapped
 *         to LEMNG return codes.
 */
ECHMET_API RetCode ECHMET_CC makeCZESystem(SysComp::InConstituentVec *BGE, SysComp::InConstituentVec *sample,
					   CZESystem *&czeSystem) ECHMET_NOEXCEPT;

/*!
 * Returns the minimum analytical concentrations of a constituent
 * that is considered safe for use by the numerical solver.
 *
 * @return Minimum safe concentration
 */
ECHMET_API double ECHMET_CC minimumSafeConcentration() ECHMET_NOEXCEPT;

/*!
 * Finds envelopes of eigenzones.
 *
 * @param[out] envelopes Vector of found envelopes. Has the same order as <tt>eigenzones</tt> vector is <tt>Results</tt>.
 * @param[in] results Results to generate the electrophoregram for.
 * @param[in] drivingVoltage Voltage applied to the system in <tt>V</tt>.
 * @param[in] totalLength Total length of the capillary in <tt>m</tt>.
 * @param[in] effectiveLength Distance between the inlet and the detector in <tt>m</tt>.
 * @param[in] EOFMobility Mobility of the electroosmotic flow in <tt>m.m/V/s . 1e-9</tt>.
 * @param[in] injectionZoneLength Length of the injection zone in <tt>m</tt>.
 * @param[in] plotToTime End the plotted electrophoregram at a specified in <tt>sec</tt>. Default value
 *                       ends the plot after the last visible eigenzone.
 */
ECHMET_API RetCode ECHMET_CC findEigenzoneEnvelopes(REigenzoneEnvelopeVec *&envelopes, const Results &results,
						    const double drivingVoltage, const double totalLength, const double effectiveLength,
						    const double EOFMobility, const double injectionZoneLength,
						    const double plotToTime) ECHMET_NOEXCEPT;

/*!
 * Plots expected electrophoregrams for given results.
 *
 * @param[out] electrophoregram Generated electrophoregram.
 * @param[in] results Results to generate the electrophoregram for.
 * @param[in] drivingVoltage Voltage applied to the system in <tt>V</tt>.
 * @param[in] totalLength Total length of the capillary in <tt>m</tt>.
 * @param[in] effectiveLength Distance between the inlet and the detector in <tt>m</tt>.
 * @param[in] EOFMobility Mobility of the electroosmotic flow in <tt>m.m/V/s . 1e-9</tt>.
 * @param[in] injectionZoneLength Length of the injection zone in <tt>m</tt>.
 * @param[in] respType Type of the response to plot.
 * @param[in] constituentName Name of the constituent whose concentration response is to ne plotted.
 *                            This parameter is ignored unless \p respType is \p RESP_CONCENTRATION.
 * @param[in] plotToTime End the plotted electrophoregram at a specified in <tt>sec</tt>. Default value
 *                       ends the plot after the last visible eigenzone.
 *
 * @retval RetCode::OK Success.
 * @retval RetCode::E_NO_MEMORY Insufficient memory to generate electrophoregram.
 * @retval RetCode::E_INTERNAL_ERROR Internal error has occured.
 * @retval RetCode::E_INVALID_ARGUMENT Concentration response was requested but no constituent name was given or nonsensical value of \p injectionZoneLength.
 * @retval RetCode::E_INVALID_CAPILLARY Nonsensical value of \p totalLength.
 * @retval RetCode::E_INVALID_DETECTOR_POSITION \p effectiveLength is greater than \p totalLength.
 */
ECHMET_API RetCode ECHMET_CC plotElectrophoregram(EFGPairVec *&electrophoregram,
						  const Results &results,
						  const double drivingVoltage, const double totalLength, const double effectiveLength,
						  const double EOFMobility,
						  const double injectionZoneLength,
						  const EFGResponseType respType,
						  const char *constituentName = ECHMET_NULLPTR,
						  const double plotToTime = -1) ECHMET_NOEXCEPT;
/*!
 * Frees resources claimed by CZESystem object.
 *
 * @param[in] czeSystem CZESystem to be released.
 */
ECHMET_API void ECHMET_CC releaseCZESystem(const CZESystem *czeSystem) ECHMET_NOEXCEPT;

/*!
 * Frees resources claimed by Results object.
 *
 * @param[in] results Results object to be released.
 */
ECHMET_API void ECHMET_CC releaseResults(Results &results) ECHMET_NOEXCEPT;

/*!
 * Sets all tracepoints to the given state.
 *
 * @param[in] state If \p true all tracepoints will be enabled and vice versa.
 */
ECHMET_API void ECHMET_CC toggleAllTracepoints(const bool state) ECHMET_NOEXCEPT;

/*!
 * Set state of one tracepoint.
 *
 * @param[in] TPID Internal ID of the tracepoint to set.
 * @param[in] state If \p true the tracepoint will be enabled and vice versa.
 */
ECHMET_API void ECHMET_CC toggleTracepoint(const int32_t TPID, const bool state) ECHMET_NOEXCEPT;

/*!
 * Returns the complete trace.
 *
 * @param[in] dontClear If \p true the trace log will not be cleared.
 *
 * @return String containing the whole trace.
 */
ECHMET_API FixedString * ECHMET_CC trace(const bool dontClear = false) ECHMET_NOEXCEPT;

/*!
 * Returns information about available tracepoints.
 *
 * @retval Pointer to a vector of all available tracepoints. May be \p NULL if
 *         the operation fails or if no tracepoins are available.
 */
ECHMET_API TracepointInfoVec * ECHMET_CC tracepointInfo() ECHMET_NOEXCEPT;

/*!
 * Returns the state of a given tracepoint.
 *
 * @param[in] TPID Internal ID of the tracepoint whose state is requested.
 *
 * @retval \p true if the tracepoint is enabled and vice versa.
 */
ECHMET_API bool ECHMET_CC tracepointState(const int32_t TPID) ECHMET_NOEXCEPT;


} // extern "C"

} // namespace LEMNG
} // namespace ECHMET

#endif // ECHMET_LEMNG_LEMNG_H
