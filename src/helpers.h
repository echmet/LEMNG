#ifndef ECHMET_LEMNG_HELPERS_H
#define ECHMET_LEMNG_HELPERS_H

#ifndef ECHMET_IMPORT_INTERNAL
#define ECHMET_IMPORT_INTERNAL
#endif // ECHMET_IMPORT_INTERNAL
#include <lemng.h>

namespace ECHMET {
namespace LEMNG {

RetCode coreLibsErrorToNativeError(const ::ECHMET::RetCode errorCode);

} // namespace LEMNG
} // namespace ECHMET

#endif // ECHMET_LEMNG_HELPERS_H
