#include <echmetmath.h>
#include "echmetmath_internal.h"

namespace ECHMET {

ECHMETReal ECHMET_CC Eabs(const ECHMETReal &arg)
{
	return VMath::abs<ECHMETReal>(arg);
}

bool ECHMET_CC Eisinf(const ECHMETReal &arg)
{
	return VMath::isinf<ECHMETReal>(arg);
}

bool ECHMET_CC Eisnan(const ECHMETReal &arg)
{
	return VMath::isnan<ECHMETReal>(arg);
}

ECHMETReal ECHMET_CC Elog10(const ECHMETReal &arg)
{
	return VMath::log10<ECHMETReal>(arg);
}

ECHMETReal ECHMET_CC Epow(const ECHMETReal &base, const ECHMETReal &exponent)
{
	return VMath::pow<ECHMETReal, ECHMETReal>(base, exponent);
}

ECHMETReal ECHMET_CC Esqrt(const ECHMETReal &arg)
{
	return VMath::sqrt<ECHMETReal>(arg);
}

} // namespace ECHMET
