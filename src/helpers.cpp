#include "helpers.h"

namespace ECHMET {
namespace LEMNG {

RetCode coreLibsErrorToNativeError(const ::ECHMET::RetCode errorCode)
{
	switch (errorCode) {
	case ::ECHMET::RetCode::E_NO_MEMORY:
		return RetCode::E_NO_MEMORY;
	case ::ECHMET::RetCode::E_INVALID_ARGUMENT:
		return RetCode::E_INVALID_ARGUMENT;
	case ::ECHMET::RetCode::E_DATA_TOO_LARGE:
		return RetCode::E_DATA_TOO_LARGE;
	case ::ECHMET::RetCode::E_DUPLICIT_CONSTITUENTS:
		return RetCode::E_DUPLICIT_CONSTITUENTS;
	case ::ECHMET::RetCode::E_NOT_IMPLEMENTED:
		return RetCode::E_NOT_IMPLEMENTED;
	case ::ECHMET::RetCode::E_INVALID_COMPLEXATION:
		return RetCode::E_INVALID_COMPLEXATION;
	case ::ECHMET::RetCode::E_INVALID_CONSTITUENT:
		return RetCode::E_INVALID_CONSTITUENT;
	case ::ECHMET::RetCode::E_NRS_FAILURE:
	case ::ECHMET::RetCode::E_NRS_NO_CONVERGENCE:
	case ::ECHMET::RetCode::E_NRS_STUCK:
	case ::ECHMET::RetCode::E_NRS_NO_SOLUTION:
		return RetCode::E_CHEM_SYSTEM_UNSOLVABLE;
	default:
		return RetCode::E_UNKW_CORELIBS_ERROR;
	}
}

} // namespace LEMNG
} // namespace ECHMET
