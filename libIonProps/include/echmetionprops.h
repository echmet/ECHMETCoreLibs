#ifndef ECHMET_IONPROPS_IONPROPS_H
#define ECHMET_IONPROPS_IONPROPS_H

#include <echmetelems.h>

#define ECHMET_IMPORT_INTERNAL
#include <echmetsyscomp.h>
#undef ECHMET_IMPORT_INTERNAL

#include <echmetmodule.h>

namespace ECHMET {
/*!
 * Calculates certain electrophoretic properties such as
 * pH, conductivitym effective mobilities and Onsager-Fuoss-corrected
 * ionic mobilities of a solution of ions with known ionic composition.
 */
namespace IonProps {

class ComputationContext {
public:
	/*!
	 * Frees resources claimed by the object.
	 */
	virtual void ECHMET_CC destroy() const ECHMET_NOEXCEPT = 0;

protected:
	virtual ~ComputationContext() ECHMET_NOEXCEPT = 0;
};

extern "C" {

/*!
 * Calculates electric conductivity of a solution.
 *
 * @param[in] ctx Computation context for the given chemical system.
 *
 * @return Conductivity of the solution in <tt>S/m</tt>. Value of the \p conductivity field in the
 *	   corresponding \p ChemicalSystem struct is set as well.
 */
ECHMET_API ECHMETReal ECHMET_CC calculateConductivity(const ComputationContext *ctx) ECHMET_NOEXCEPT;

/*!
 * Calculates effective electrophoretic mobilities of all compounds in the system.
 * Calculated effective mobilities are stored in the \p effectiveMobilities vector
 * of the corresponding \p ChemicalSystem struct.
 *
 * @param[in] ctx Computation context for the given chemical system.
 *
 * @retval RetCode::OK Success.
 * @retval RetCode::E_NO_MEMORY Insufficient memory to complete the calculation.
 * @retval RetCode::E_BAD_INPUT Invalid chemical system.
 */
ECHMET_API RetCode ECHMET_CC calculateEffectiveMobilities(ComputationContext *ctx) ECHMET_NOEXCEPT;

/*!
 * Calculates pH of the system.
 *
 * @param[in] ctx Computation context for the given chemical system.
 * @param[in] isCorrection Correct for ionic strength.
 *
 * @return pH value.
 */
ECHMET_API ECHMETReal ECHMET_CC calculatepH(const ComputationContext *ctx, const bool isCorrection) ECHMET_NOEXCEPT;

/*!
 * Returns pH value corresponding to the given concentration of H<sub>3</sub>O<sup>+</sup> ions.
 *
 * @param[in] cH Concentration of H<sub>3</sub>O<sup>+</sup> ions.
 * @param[in] ionicStrength Ionic strength of the system in <tt>mol/dm<sup>3</sup></tt>. Supply 0 to calculate pH
 *            from H<sub>3</sub>O<sup>+</sup> concentration instead of activity.

 * @return pH value.
 */
ECHMET_API ECHMETReal ECHMET_CC calculatepH_direct(const ECHMETReal &cH, const ECHMETReal &ionicStrength) ECHMET_NOEXCEPT;

/*!
 * Corrects limit electrophoretic mobilities for ionic atmosphere effect
 * using Onsager-Fuoss law.
 * <b>Note: If this correction is required, this function shall be called prior to
 * calling \p calculateConductivity() or \p calculateEffectiveMobilities(). Reversing
 * the calling order will result in the conductivity or effective mobilities being
 * calculated with ionic mobilities that do not include the Onsager-Fuoss correction.
 *
 * @param[in] ctx Computation context for the given chemical system.
 *
 * @retval RetCode::OK Success
 * @retval RetCode::E_NO_MEMORY Insufficient memory to complete operation
 * @retval RetCode::E_DATA_TOO_LARGE System is too large to solve
 */
ECHMET_API RetCode ECHMET_CC correctMobilities(ComputationContext *ctx) ECHMET_NOEXCEPT;

/*!
 * Makes computation context for the given chemical system and analytical concentrations.
 * <b>Note: The \p ComputationContext object holds references to all objects passed as
 * parameters to this function. Using the \p ComputationContext after any of these objects
 * has been destroyed will result in undefined behavior.</b>
 *
 * @param[in] chemSystem The chemical system.
 * @param[in] analyticalConcentration The analytical concentrations.
 * @param[in] calcProps Calculated properties solved by \p CAES.
 *
 * @return A pointer to \p ComputationContext object, \p NULL if the context cannot be created.
 */
ECHMET_API ComputationContext * ECHMET_CC makeComputationContext(const SysComp::ChemicalSystem &chemSystem, const RealVec *analyticalConcentrations, SysComp::CalculatedProperties &calcProps) ECHMET_NOEXCEPT;

}

} // namespace IonProps
} // namespace CAES

#endif // ECHMET_IONPROPS_IONPROPS_H
