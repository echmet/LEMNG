#ifndef _ECHMET_TRACER_BASE_H
#define _ECHMET_TRACER_BASE_H

#include "tracer_types.h"
#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace ECHMET {

template <typename TracepointIDs>
class Tracer
{
	static_assert(std::is_enum<TracepointIDs>::value, "TracepointIDs is not an enum");
	static_assert(sizeof(typename std::underlying_type<TracepointIDs>::type) <= sizeof(TPIDInt), "Cannot represent all tracepoints as int32_t");

public:
	void disableAllTracepoints()
	{
#ifndef ECHMET_TRACER_DISABLE_TRACING
		TOGGLE_ALL_TRACEPOINTS<TracepointIDs, FIRST_TRACEPOINT_ID<TracepointIDs>()>(false, m_enabledTracepoints);
#endif // ECHMET_TRACER_DISABLE_TRACING
	}

	template <typename RTPID>
	void disableTracepoint(const RTPID &tpid)
	{
#ifndef ECHMET_TRACER_DISABLE_TRACING
		if (!IS_TPID_VALID<RTPID, TracepointIDs>(tpid))
			return;

		const TracepointIDs _tpid  = static_cast<TracepointIDs>(tpid);
		m_enabledTracepoints[_tpid] = false;
#else
		(void)tpid;
		return;
#endif // ECHMET_TRACER_DISABLE_TRACING
	}

	void enableAllTracepoints()
	{
#ifndef ECHMET_TRACER_DISABLE_TRACING
		TOGGLE_ALL_TRACEPOINTS<TracepointIDs, FIRST_TRACEPOINT_ID<TracepointIDs>()>(true, m_enabledTracepoints);
#endif // ECHMET_TRACER_DISABLE_TRACING
	}

	template <typename RTPID>
	void enableTracepoint(const RTPID &tpid)
	{
#ifndef ECHMET_TRACER_DISABLE_TRACING
		if (!IS_TPID_VALID<RTPID, TracepointIDs>(tpid))
			return;

		const TracepointIDs _tpid  = static_cast<TracepointIDs>(tpid);
		m_enabledTracepoints[_tpid] = true;
#else
		(void)tpid;
		return;
#endif // TRACER_DISABLE_TRACING
	}

	static Tracer<TracepointIDs> & instance()
	{
		if (s_instance == nullptr)
			s_instance = std::unique_ptr<Tracer<TracepointIDs>>(new Tracer<TracepointIDs>{});

		return *s_instance.get();
	}

	template <typename RTPID>
	bool isTracepointEnabled(const RTPID &tpid) const
	{
#ifndef ECHMET_TRACER_DISABLE_TRACING
		if (!IS_TPID_VALID<RTPID, TracepointIDs>(tpid))
			return false;

		const TracepointIDs _tpid = static_cast<TracepointIDs>(tpid);
		auto item = m_enabledTracepoints.find(_tpid);
		if (item == m_enabledTracepoints.end())
			return false;

		return item->second;
#else
		(void)tpid;
		return false;
#endif // ECHMET_TRACER_DISABLE_TRACING
	}

	void log(const std::string &text)
	{
		std::lock_guard<std::mutex> lk(m_logLock);

		m_log.append(text + "\n");
	}

	std::string logged(const bool dontFlush = false)
	{
		std::lock_guard<std::mutex> lk(m_logLock);

		if (dontFlush)
			return m_log;
		else {
			const std::string log{m_log};
			m_log.clear();

			return log;
		}
	}

	std::vector<std::tuple<TPIDInt, std::string>> tracepoints() const
	{
		std::vector<std::tuple<TPIDInt, std::string>> tpVec{};

		TRACEPOINT_INFO_BUILD<TracepointIDs, FIRST_TRACEPOINT_ID<TracepointIDs>()>(tpVec);
		return tpVec;
	}

private:
	std::map<TracepointIDs, bool> m_enabledTracepoints;
	std::string m_log;
	std::mutex m_logLock;

	static std::unique_ptr<Tracer<TracepointIDs>> s_instance;
};

template <typename TracepointIDs>
Tracer<TracepointIDs> & TRACER_INSTANCE();

#ifndef ECHMET_TRACER_DISABLE_TRACING
template <typename TracepointIDs, TracepointIDs TPID, typename ...Args>
void _ECHMET_TRACE(Args ...args)
{
	auto &tracer = TRACER_INSTANCE<TracepointIDs>();
	if (tracer.isTracepointEnabled(TPID))
		tracer.log(TRACEPOINT_LOGGER<TracepointIDs, TPID, Args...>(args...));
		//tracer.log(TRACEPOINT_LOGGER<TracepointIDs, TPID, Args...>(std::forward<Args>(args)...));
	/* Do nothing */
}
#else
template <typename TracepointIDs, TracepointIDs TPID, typename ...Args>
void _ECHMET_TRACE(Args && ...) {} /* Do nothing */
#endif // ECHMET_TRACER_DISABLE_TRACING

} // namespace ECHMET

/*!
 * \def MAKE_TRACER(TracerClass)
 * Defines an instance of the given \TracerClass
 *
 * @param TracerClass Tracer class
 */
#define ECHMET_MAKE_TRACER(TracerClass) \
	namespace ECHMET { \
		template <> \
		std::unique_ptr<Tracer<::TracerClass>> Tracer<::TracerClass>::s_instance{nullptr}; \
		template <> \
		Tracer<::TracerClass> & TRACER_INSTANCE<::TracerClass>() \
		{ \
			return Tracer<::TracerClass>::instance(); \
		} \
	} // namespace ECHMET

#ifndef ECHMET_TRACER_DISABLE_TRACING
/*!
 * \def GET_TRACER_INSTANCE(TracerClass)
 * Returns an instance of the given \TracerClass
 *
 * @param TracerClass Tracer class
 */
#define ECHMET_GET_TRACER_INSTANCE(TracerClass) \
	TRACER_INSTANCE<::TracerClass>()
#else
#define ECHMET_GET_TRACER_INSTANCE(TracerClass) \
	TRACER_INSTANCE<::ECHMET::__DUMMY_TRACER_CLASS>()
#endif // ECHMET_TRACER_DISABLE_TRACING

#endif // _ECHMET_TRACER_BASE_H
