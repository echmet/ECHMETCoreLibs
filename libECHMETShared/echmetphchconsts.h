#ifndef ECHMET_PHCHCONSTS_H
#define ECHMET_PHCHCONSTS_H

#include <echmetelems.h>

namespace ECHMET {

/*!
 * Physical-chemical constants used by the libraries
 */
namespace PhChConsts {

	const ECHMETReal bk = 1.38064852e-23;		/*!< Boltzmann constant <tt>(m^2 \. kg \. s^(-2) K^(-1))</tt> */
	const ECHMETReal e = 1.60217662e-19;		/*!< Absolute value of electron charge <tt>(C)</tt> */
	const ECHMETReal eAbsVac = 8.854187817e-12;	/*!< Absolute permittivity of vacuum <tt>(F / m)</tt> */
	const ECHMETReal dkWat = 78.3043;		/*!< Relative permittivity of water at 298.15 K */
	const ECHMETReal F = 96485.3392;		/*!< Faraday constant <tt>(C / mol)</tt> */
	const ECHMETReal KW_298 = 1.0e-14;		/*!< Ionic product of water at 298\.15 K */
	const ECHMETReal NA = 6.022140857e23;		/*!< Avogadro's constant <tt>mol^(-1)</tt> */
	const ECHMETReal mobilityH3O = 362.5;		/*!< Limit ionic mobility of the hydroxonium ion * 1e9 */
	const ECHMETReal mobilityOH = 205.5;		/*!< Limit ionic mobility of the hydroxyl ion * 1e9 */
	const ECHMETReal th = 8.9633e-4;		/*!< Water viscosity coefficient <tt>(Pa \. s)</tt>*/
	const ECHMETReal Tlab = 298.15;			/*!< Laboratory temperature of 298\.15 K */

	const ECHMETReal RCs[] = { 0.2929, -0.3536, 0.0884, -0.0442, 0.0276, -0.0193 };    /*!< Precomputed values of Rn coefficients. Used in Onsager-Fuoss ionic mobility correction */
	const size_t RCsCount = 6;

extern "C" {
	/*!
	 * Returns p-Scaled activity coefficent for a given ionic strength in <tt>mol/dm3</tt> and charge.
	 *
	 * @param[in] is Ionic strength in <tt>mol/dm3</tt>.
	 * @param[in] charge Electric charge
	 *
	 * @return p-Scaled activity coefficient
	 */
	ECHMET_API ECHMETReal ECHMET_CC pActivityCoefficient(const ECHMETReal &is, const int charge) noexcept;

	/*!
	 * Returns activity coefficient for a given ionic strength in <tt>mol/dm3</tt> and charge.
	 *
	 * @param[in] is Ionic strength in <tt>mol/dm3</tt>.
	 * @param[in] charge Electric charge
	 *
	 * @return Activity coefficient
	 */
	ECHMET_API ECHMETReal ECHMET_CC activityCoefficient(const ECHMETReal &is, const int charge) noexcept;
}

} // namespace PhChConsts
} // namespace ECHMET

#endif // ECHMET_PHCHCONSTS_H
