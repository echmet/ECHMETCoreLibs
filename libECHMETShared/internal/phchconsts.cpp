#include <echmetphchconsts.h>
#include "phchconsts_calcs.hpp"

namespace ECHMET {
namespace PhChConsts {

ECHMETReal ECHMET_CC pActivityCoefficient(const ECHMETReal &is, const int charge) ECHMET_NOEXCEPT
{
	const ECHMETReal sqrtIs = VMath::sqrt(is);
	const int chSq = charge * charge;

	return pActivityCoefficientInternal<ECHMETReal>(is, sqrtIs, chSq);
}

ECHMETReal ECHMET_CC activityCoefficient(const ECHMETReal &is, const int charge) ECHMET_NOEXCEPT
{
	const ECHMETReal sqrtIs = VMath::sqrt(is);
	const int chSq = charge * charge;

	return activityCoefficientInternal<ECHMETReal>(is, sqrtIs, chSq);
}

} // namespace PhChConsts
} // namespace ECHMET

