#ifndef ECHMET_CAES_CAES_EXTENDED_H
#define ECHMET_CAES_CAES_EXTENDED_H

#define ECHMET_IMPORT_INTERNAL
#include <echmetcaes.h>
#include <echmetsyscomp.h>
#undef ECHMET_IMPORT_INTERNAL

#include <echmetmodule.h>

namespace ECHMET {
namespace CAES {

/*!
 * Public context for dissociation degrees derivatives calculations
 */
class DDSContext {
public:
	virtual void ECHMET_CC destroy() const ECHMET_NOEXCEPT = 0;
	/*!
	 * Returns a dissociation degree derivative for a given ionic form.
	 *
	 * @param[out] value Reference to variable that will contain the result.
	 * @param[in] name Name of the ionic form whose dissociation degree derivative is requested.
	 *
	 * @retval RetCode::OK Success.
	 * @retval RetCode::E_NOT_FOUND The given ionic form is not present in the system.
	 */
	virtual RetCode ECHMET_CC findDissocDegreeDerivative(ECHMETReal &value, const FixedString *name) const ECHMET_NOEXCEPT = 0;

protected:
	virtual ~DDSContext() ECHMET_NOEXCEPT = 0;
};

extern "C" {

/*!
 * Calculates buffer capacity of a solution. Equilibrium of the system shall be solver
 * prior to calling this function.
 *
 * @param[out] bufferCapacity Computed buffer capacity.
 * @param[in] chemSystem The chemical system to solve.
 * @param[in] calcProps Calculated properties of the system.
 * @param[in] analyticalConcentrations Vector of analytical concentrations of constituents.
 * @param[in] isCorrection Correct for ionic strength.
 *
 * @retval RetCode::OK Success
 * @retval RetCode::E_NO_MEMORY Insufficient memory to perform the calculation.
 * @retval RetCode::E_BUFFER_CAPACITY_UNSOLVABLE Buffer capacity cannot be calculated for the given system.
 */
ECHMET_API RetCode ECHMET_CC calculateBufferCapacity(ECHMETReal &bufferCapacity, const NonidealityCorrections corrections, const SysComp::ChemicalSystem &chemSystem, const SysComp::CalculatedProperties &calcProps,
						     const RealVec *analyticalConcentrations) ECHMET_NOEXCEPT;

/*!
 * Returns first derivatives of all ionic concentrations when the system is perturbed by a small
 * change of concentration of one constituent.
 *
 * @param[out] derivatives Output vector of the derivatives. Derivatives are sorted in the same order as ionic concentrations.
 * @param[out] conductivityDerivative First derivative of conductivity.
 * @param[in] H Magnitude of the concentration perturbance in <tt>mmol/dm3</tt>.
 * @param[in] isCorrection Correct for ionic strength.
 * @param[in] chemSystem The chemical system to solve.
 * @param[in] analyticalConcentration Vector of analytical concentrations of all constituents.
 * @param[in] perturbedConstituent Constituent whose concentrations is to be perturbed.
 *
 * @retval RetCode::OK Success.
 * @retval RetCode::E_INVALID_ARGUMENT Invalid parameter was passed to the function.
 * @retval RetCode::E_NO_MEMORY Insufficient memory to perform the calculation.
 * @retval . Any other \p RetCode value returned by \p executor.
 */
ECHMET_API RetCode ECHMET_CC calculateFirstConcentrationDerivatives(RealVec *&derivatives, ECHMETReal &conductivityDerivative, const ECHMETReal &H, const NonidealityCorrections corrections, const SysComp::ChemicalSystem &chemSystem, const RealVec *analyticalConcentrations, const SysComp::Constituent *perturbedConstituent, const ECHMETReal &ionicStrength) ECHMET_NOEXCEPT;

/*!
 * Returns mixed derivatives of all ionic concentrations when the system is perturbed by a small
 * change of concentration of two constituents. If the two constituents are the same, second derivative
 * is returned instead.
 *
 * @param[out] derivatives Output vector of the derivatives. Derivatives are sorted in the same order as ionic concentrations.
 * @param[in] H Magnitude of the concentration perturbance in <tt>mmol/dm3</tt>.
 * @param[in] isCorrection Correct for ionic strength.
 * @param[in] chemSystem The chemical system to solve.
 * @param[in] analyticalConcentration Vector of analytical concentrations of all constituents.
 * @param[in] perturbedConstituentJ First constituent whose concentrations is to be perturbed.
 * @param[in] perturbedConstituentJ Second constituent whose concentrations is to be perturbed.
 *
 * @retval RetCode::OK Success.
 * @retval RetCode::E_INVALID_ARGUMENT Invalid parameter was passed to the function.
 * @retval RetCode::E_NO_MEMORY Insufficient memory to perform the calculation.
 * @retval . Any other \p RetCode value returned by \p executor.
 */
ECHMET_API RetCode ECHMET_CC calculateCrossConcentrationDerivatives(RealVec *&derivatives, const ECHMETReal &H, const NonidealityCorrections corrections, const SysComp::ChemicalSystem &chemSystem, const RealVec *analyticalConcentrations, const SysComp::Constituent *perturbedConstituentJ, const SysComp::Constituent *perturbedConstituentK, const ECHMETReal &ionicStrength) ECHMET_NOEXCEPT;

/*!
 * Returns first derivatives of all ionic concentrations when the system is perturbed by a small
 * change of concentration of one constituent. This function accepts already initialized
 * \p Solver and RealVec to store the derivatives. Passing mismatching Solver, vector of
 * derivatives and \p ChemicalSystem may result in undefined behavior.
 *
 * @param[in,out] derivatives Output vector of the derivatives. Derivatives are sorted in the same order as ionic concentrations.
 * @param[in,out] conductivityDerivative First derivative of conductivity.
 * @param[in] Associated solver.
 * @param[in] H Magnitude of the concentration perturbance in <tt>mmol/dm3</tt>.
 * @param[in] chemSystem The chemical system to solve.
 * @param[in] analyticalConcentration Vector of analytical concentrations of all constituents.
 * @param[in] perturbedConstituent Constituent whose concentrations is to be perturbed.
 *
 * @retval RetCode::OK Success.
 * @retval RetCode::E_INVALID_ARGUMENT Invalid parameter was passed to the function.
 * @retval RetCode::E_NO_MEMORY Insufficient memory to perform the calculation.
 * @retval . Any other \p RetCode value returned by \p executor.
 */
ECHMET_API RetCode ECHMET_CC calculateFirstConcentrationDerivatives_prepared(RealVec *derivatives, ECHMETReal &conductivityDerivative, Solver *solver, const ECHMETReal &H, const NonidealityCorrections corrections, const SysComp::ChemicalSystem &chemSystem, const RealVec *analyticalConcentrations, const SysComp::Constituent *perturbedConstituent, const ECHMETReal &ionicStrength) ECHMET_NOEXCEPT;

/*!
 * Returns mixed derivatives of all ionic concentrations when the system is perturbed by a small
 * change of concentration of two constituents. If the two constituents are the same, second derivative
 * is returned instead. This function accepts already initialized \p Solver and RealVec to store
 * the derivatives. Passing mismatching Solver, vector of derivatives and \p ChemicalSystem may
 * result in undefined behavior.
 *
 * @param[out] derivatives Output vector of the derivatives. Derivatives are sorted in the same order as ionic concentrations.
 * @param[in] H Magnitude of the concentration perturbance in <tt>mmol/dm3</tt>.
 * @param[in] isCorrection Correct for ionic strength.
 * @param[in] chemSystem The chemical system to solve.
 * @param[in] analyticalConcentration Vector of analytical concentrations of all constituents.
 * @param[in] perturbedConstituentJ First constituent whose concentrations is to be perturbed.
 * @param[in] perturbedConstituentJ Second constituent whose concentrations is to be perturbed.
 *
 * @retval RetCode::OK Success.
 * @retval RetCode::E_INVALID_ARGUMENT Invalid parameter was passed to the function.
 * @retval RetCode::E_NO_MEMORY Insufficient memory to perform the calculation.
 * @retval . Any other \p RetCode value returned by \p executor.
 */
ECHMET_API RetCode ECHMET_CC calculateCrossConcentrationDerivatives_prepared(RealVec *derivatives, Solver *solver, const ECHMETReal &H, const SysComp::ChemicalSystem &chemSystem, const RealVec *analyticalConcentrations, const SysComp::Constituent *perturbedConstituentJ, const SysComp::Constituent *perturbedConstituentK, const ECHMETReal &ionicStrength) ECHMET_NOEXCEPT;

/*!
 * Prepares solver and output vector of results for derivators.
 *
 * @param[out] derivatives Output vector of the derivatives with the appropriate size.
 * @param[out] solver Associated solver.
 * @param[in] chemSystem The chemical system to solve.
 * @param[in] isCorrection Correct for ionic strength.
 *
 * @retval RetCode::OK Success.
 * @retval RetCode::E_NO_MEMORY Insufficient memory to prepare the context.
 */
ECHMET_API RetCode ECHMET_CC prepareDerivatorContext(RealVec *&derivatives, Solver *&solver, const SysComp::ChemicalSystem &chemSystem, const NonidealityCorrections corrections) ECHMET_NOEXCEPT;

} // extern "C"

} // namespace CAES
} // namespace ECHMET

#endif // ECHMET_CAES_CAES_EXTENDED_H
