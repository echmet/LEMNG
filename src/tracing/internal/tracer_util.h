#ifndef _ECHMET_TRACER_UTIL_H
#define _ECHMET_TRACER_UTIL_H

#include "tracer_types.h"
#include <map>
#include <tuple>
#include <vector>

namespace ECHMET {

/*!
 * Returns human-readable information about a tracepoint
 *
 * @tparam TracerClass Tracer class
 * @tparam TPID ID of the requested tracepoint
 * @return Tracepoint object
 */
template <typename TracerClass, TracerClass TPID>
Tracepoint<TracerClass> TRACEPOINT_INFO();

/*!
 * Returns ID of the first tracepoint for a given tracer
 *
 * @tparam TracerClass Tracer class
 * @return ID of the first tracepoint
 */
template <typename TracerClass>
static constexpr TracerClass FIRST_TRACEPOINT_ID();

/*!
 * Returns ID of the next tracepoint
 *
 * @tparam TracerClass Tracer class
 * @tparam TPID ID of the tracepoint to get the next tracepoint for
 * @return ID of the next tracepoint
 */
template <typename TracerClass, TracerClass TPID>
static constexpr TracerClass NEXT_TRACEPOINT_ID()
{
	typedef typename std::underlying_type<TracerClass>::type UType;
	return static_cast<TracerClass>(static_cast<UType>(TPID) + 1);
}

template <typename TracerClass>
static constexpr TracerClass LAST_TRACEPOINT_ID();

template <typename TracerClass>
struct TUTYPE {
	typedef typename std::underlying_type<TracerClass>::type type;
};

template <typename TracerClass, typename UType = typename std::underlying_type<TracerClass>::type>
static constexpr UType TUTYPE_CAST(const TracerClass &TPID)
{
	return static_cast<UType>(TPID);
}

template <typename RTPID, typename TracerClass>
static bool IS_TPID_VALID(const RTPID &tpid)
{
	static_assert(std::is_same<RTPID, typename TUTYPE<TracerClass>::type>::value, "Incompatible raw tracepoint IDs");

	const typename TUTYPE<TracerClass>::type first = TUTYPE_CAST(FIRST_TRACEPOINT_ID<TracerClass>());
	const typename TUTYPE<TracerClass>::type last = TUTYPE_CAST(LAST_TRACEPOINT_ID<TracerClass>());

	return (tpid >= first && tpid < last);
}

/*!
 * Builds a vector with information about all tracepoints available for
 * a given tracer. This function is recursive, each call appends a new
 * entry to the vector.
 *
 * @tparam TracerClass Tracer class
 * @tparam TPID ID of the currently being added tracepoint
 * @return Vector of Tracepoint<TracerClass> objects
 */
template <typename TracerClass, TracerClass TPID>
static void TRACEPOINT_INFO_BUILD(std::vector<std::tuple<TPIDInt, std::string>> &tracepointInfoVec)
{
#ifndef ECHMET_TRACER_DISABLE_TRACING
	const auto tpinfo = TRACEPOINT_INFO<TracerClass, TPID>();
	tracepointInfoVec.emplace_back(static_cast<TPIDInt>(tpinfo.ID), tpinfo.description);
	TRACEPOINT_INFO_BUILD<TracerClass, NEXT_TRACEPOINT_ID<TracerClass, TPID>()>(tracepointInfoVec);
#else
	(void)tracepointInfoVec;
	return;
#endif // TRACER_DISABLE_TRACING
}

template <typename TracerClass, TracerClass TPID>
static void TOGGLE_ALL_TRACEPOINTS(const bool state, std::map<TracerClass, bool> &enabledTracepoints)
{
#ifndef ECHMET_TRACER_DISABLE_TRACING
	enabledTracepoints[TPID] = state;
	TOGGLE_ALL_TRACEPOINTS<TracerClass, NEXT_TRACEPOINT_ID<TracerClass, TPID>()>(state, enabledTracepoints);
#else
	(void)enabledTracepoints;
	return;
#endif // ECHMET_TRACER_DISABLE_TRACING
}

/*!
 * Logging function.
 *
 * @tparam Tracer class
 * @tparam TPID ID of the tracepoint
 * @tparam Args Argument template of the logging function
 * @param[in] Args... Arguments of the logging function
 *
 * @return String to be logged
 */
template <typename TracerClass, TracerClass TPID, typename... Args>
static std::string TRACEPOINT_LOGGER(Args...);

} // namespace ECHMET

/*!
 * \def MAKE_TRACEPOINT(TracerClass, TPID, description)
 * Defines functions necessary to query information about tracepoints for the given \TracerClass
 *
 * @param TracerClass Tracer class
 * @param TPID ID of the tracepoint
 * @param description Human-readable description of the tracepoint
 */
#define ECHMET_MAKE_TRACEPOINT(TracerClass, TPID, description) \
	template <> \
	::ECHMET::Tracepoint<::TracerClass> ECHMET::TRACEPOINT_INFO<::TracerClass, ::TracerClass::TPID>() { return ::ECHMET::Tracepoint<::TracerClass>{::TracerClass::TPID, description}; }

/*!
 * \def MAKE_LOGGER(TracerClass, TPID, Args...)
 * Defines a logging function for a given tracer and its tracepoint
 *
 * @param TracerClass Tracer class
 * @param TPID ID of the tracepoint whose logging function is being declared
 * @Args... Argument template of the logging function
 */
#define ECHMET_MAKE_LOGGER(TracerClass, TPID, ...) \
	template <> \
	std::string ECHMET::TRACEPOINT_LOGGER<::TracerClass, ::TracerClass::TPID>(__VA_ARGS__)

#ifndef _ECHMET_TRACER_IMPL_SECTION
	#ifndef ECHMET_TRACER_DISABLE_TRACING
	/*!
	 * \def MAKE_TRAECPOINT_IDS(TracerClass, first, last)
	 * Defines functions necessary to build a list of tracepoints for the given \TracerClass
	 *
	 * @param TracerClass Tracer class
	 * @param first ID of the first tracepoint
	 * @param last ID of the last tracepoint
	 */
	#define ECHMET_MAKE_TRACEPOINT_IDS(TracerClass, first, last) \
		template <> \
		constexpr TracerClass ECHMET::FIRST_TRACEPOINT_ID<::TracerClass>() { return ::TracerClass::first; } \
		template <> \
		constexpr TracerClass ECHMET::LAST_TRACEPOINT_ID<::TracerClass>() { return ::TracerClass::last; } \
		template <> \
		void ECHMET::TRACEPOINT_INFO_BUILD<::TracerClass, ::TracerClass::last>(std::vector<std::tuple<TPIDInt, std::string>> &tracepointInfoVec) \
		{ \
			(void)tracepointInfoVec; \
			return; /* No-op for the last dummy tracepoint */ \
		} \
		template <> \
		void ECHMET::TOGGLE_ALL_TRACEPOINTS<::TracerClass, ::TracerClass::last>(const bool state, std::map<TracerClass, bool> &enabledTracepoints) \
		{ \
			enabledTracepoints[TracerClass::last] = state; \
		} \
		template <> \
		bool ECHMET::IS_TPID_VALID<::TracerClass, ::TracerClass>(const ::TracerClass &) { return true; }
	#else
	#define ECHMET_MAKE_TRACEPOINT_IDS(TracerClass, first, last)

	template <>
	constexpr ECHMET::__DUMMY_TRACER_CLASS ECHMET::FIRST_TRACEPOINT_ID<ECHMET::__DUMMY_TRACER_CLASS>() { return ECHMET::__DUMMY_TRACER_CLASS::NONE; }
	template <>
	void ECHMET::TRACEPOINT_INFO_BUILD<ECHMET::__DUMMY_TRACER_CLASS, ECHMET::__DUMMY_TRACER_CLASS::NONE>(std::vector<std::tuple<TPIDInt, std::string>> &)
	{
		return;
	}
	template <>
	void ECHMET::TOGGLE_ALL_TRACEPOINTS<ECHMET::__DUMMY_TRACER_CLASS, ECHMET::__DUMMY_TRACER_CLASS::NONE>(const bool, std::map<ECHMET::__DUMMY_TRACER_CLASS, bool> &)
	{
		return;
	}
	#endif // TRACER_DISABLE_TRACING
#else

#define ECHMET_MAKE_TRACEPOINT_IDS(TracerClass, first, last)

#endif // _ECHMET_TRACER_IMPL_SECTION

#endif // _ECHMET_TRACER_UTIL_H
