#ifndef _ECHMET_TRACER_USER_H
#define _ECHMET_TRACER_USER_H

#include "tracer_base.h"

#ifndef ECHMET_TRACER_DISABLE_TRACING

/*!
 * \def ECHMET_S(str)
 * Converts \str to C string
 */
#define ECHMET_S(str) \
	std::string{str}.c_str()

/*!
 * \def ECHMET_TRACE(TracerClass, TPID, Args...)
 * Logs data for the given \TracerClass and tracepoint
 *
 * @param TracerClass Tracer class
 * @param ID of the tracepoint to log
 * @param Args... Arguments for the logging function
 */
#define ECHMET_TRACE(TracerClass, TPID, ...) \
	_ECHMET_TRACE<TracerClass, TracerClass::TPID>(__VA_ARGS__)

/*!
 * \def ECHMET_TRACER_LOG(TracerClass)
 * Returns the complete log from a given \TracerClass
 *
 * @param TracerClass Tracer class
 */
#define ECHMET_TRACER_LOG(TracerClass) \
	TRACER_INSTANCE<TracerClass>().logged()

#else

#define ECHMET_S(str)
#define ECHMET_TRACE(TraceClass, TPID, ...)
#define ECHMET_TRACER_LOG(TracerClass) std::string{}

#endif // ECHMET_TRACER_DISABLE_TRACING

#endif // _ECHMET_TRACER_USER_H
