#ifndef PHCHCONSTS_CALCS_HPP
#define PHCHCONSTS_CALCS_HPP

#include "echmetmath_internal.h"

namespace ECHMET {

/*!
 * Returns p-Scaled activity coefficent for a given ionic strength in <tt>mol/dm3</tt> and charge.
 * This is an internal templated implementation.
 *
 * @param[in] is Ionic strength in <tt>mol/dm3</tt>.
 * @param[in] sqrtIs Square root of the ionic strength in <tt>mol/dm3</tt>
 * @param[in] chSq Squared electric charge
 */
template <typename XReal>
XReal pActivityCoefficientInternal(const XReal &is, const XReal &sqrtIs, const int chSq) noexcept
{
	return -(-((0.50925 * chSq * sqrtIs) / (1.0 + (1.5 * sqrtIs))) + (0.1 * chSq * is));
}

/*!
 * Returns activity coefficient for a given ionic strength in <tt>mol/dm3</tt> and charge.
 * This is an internal templated implementation.
 *
 * @param[in] is Ionic strength in <tt>mol/dm3</tt>.
 * @param[in] sqrtIs Square root of the ionic strength in <tt>mol/dm3</tt>
 * @param[in] chSq Squared electric charge
 */
template <typename XReal>
XReal activityCoefficientInternal(const XReal &is, const XReal &sqrtIs, const int chSq) noexcept
{
	return ECHMET::VMath::pow<XReal>(10.0, -pActivityCoefficientInternal(is, sqrtIs, chSq));
}

} // ECHMET

#endif // PHCHCONSTS_CALCS_HPP
