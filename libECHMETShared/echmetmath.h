#ifndef ECHMET_MATH_EXPOSED_H
#define ECHMET_MATH_EXPOSED_H

#include <echmetelems.h>


namespace ECHMET {

/*!
 * External math functions that operate on ECHMETSharedLibs'
 * native Real number.
 */
extern "C" {

ECHMET_API ECHMETReal ECHMET_CC Eabs(const ECHMETReal &arg);
ECHMET_API bool ECHMET_CC Eisinf(const ECHMETReal &arg);
ECHMET_API bool ECHMET_CC Eisnan(const ECHMETReal &arg);
ECHMET_API ECHMETReal ECHMET_CC Elog10(const ECHMETReal &arg);
ECHMET_API ECHMETReal ECHMET_CC Epow(const ECHMETReal &base, const ECHMETReal &exponent);
ECHMET_API ECHMETReal ECHMET_CC Esqrt(const ECHMETReal &arg);

}

} // namespace ECHMET

#endif // ECHMET_MATH_EXPOSED_H

