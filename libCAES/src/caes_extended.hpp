#ifndef ECHMET_CAES_CAES_EXTENDED_HPP
#define ECHMET_CAES_CAES_EXTENDED_HPP

#include <internal/echmetmath_internal.h>
#include <algorithm>
#include <cassert>
#include <map>
#include <memory>

#define ECHMET_IMPORT_INTERNAL
#include <containers/echmetvec_p.h>

namespace ECHMET {
namespace CAES {

using MPRealVec = Vec<mpfr::mpreal>;

inline
void releaseMPRealVec(MPRealVec *vec) { vec->destroy(); }

using MPRealVecWrap = std::unique_ptr<MPRealVec, decltype(&releaseMPRealVec)>;

/*!
 * Wrapper function for the derivators.
 * The wrapper performs the necessary initialization and cleanup.
 *
 * @param[out] derivatives Vector to contain the results of derivations.
 * @param[in] H The H by which to perturb the system.
 * @param[in] isCorrection Correct the system composition for ionic strength.
 * @param[in] chemSystem The chemical system to solve.
 * @param[in] analyticalConcentrations Analytical concentrations of constituents.
 * @param[in] inIonicStrength Ionic strength of unperturbed system
 * @param[in] executor Function computing the derivation.
 * @param[in] params List of additional parameters for the \p executor.
 *
 * @retval RetCode::OK Success.
 * @retval RetCode::E_INVALID_ARGUMENT Invalid parameter was passed to the function.
 * @retval RetCode::E_NO_MEMORY Insufficient memory to perform the calculation.
 * @retval . Any other \p RetCode value returned by \p executor.
 */
template <InstructionSet ISet, typename... EParams>
static
RetCode derivatorSkin(RealVec *derivatives, const ECHMETReal &H, Solver *solver, const SysComp::ChemicalSystem &chemSystem,
		      const RealVec *analyticalConcentrations, const ECHMETReal &inIonicStrength,
		      std::function<RetCode (RealVec *, const ECHMETReal &, SolverImpl<mpfr::mpreal, ISet, true> *, const SysComp::ChemicalSystem &, const RealVec *, const MPRealVecWrap &, EParams...)> &executor, EParams... params)
{
	RetCode tRet;
	const long currentMpfrPrec = mpfr::mpreal::get_default_prec();

	if (H <= ECHMETReal{0})
		return RetCode::E_INVALID_ARGUMENT;

	if (derivatives != nullptr) {
		if (derivatives->size() != chemSystem.ionicForms->size())
			return RetCode::E_INVALID_ARGUMENT;
	}

	auto solverImpl = dynamic_cast<SolverImpl<mpfr::mpreal, ISet, true> *>(solver);
	if (solverImpl == nullptr)
		return RetCode::E_INVALID_ARGUMENT;

	mpfr::mpreal::set_default_prec(mpfr::digits2bits(200));

	MPRealVecWrap estimatedConcentrations(createECHMETVec<mpfr::mpreal, false>(solverImpl->contextInternal()->concentrationCount), releaseMPRealVec);
	if (estimatedConcentrations == nullptr)
		return RetCode::E_NO_MEMORY;

	estimatedConcentrations->resize(solverImpl->contextInternal()->concentrationCount);

	mpfr::mpreal ionicStrength = mpfr::mpreal(inIonicStrength);
	tRet = solverImpl->estimateDistributionSafeInternal(analyticalConcentrations, estimatedConcentrations->data(), ionicStrength);
	if (tRet != RetCode::OK) {
		mpfr::mpreal::set_default_prec(currentMpfrPrec);
		return tRet;
	}

	tRet = executor(derivatives, H, solverImpl, chemSystem, analyticalConcentrations, estimatedConcentrations, params...);

	mpfr::mpreal::set_default_prec(currentMpfrPrec);

	return tRet;
}

} // namespace CAES
} // namespace ECHMET

#endif // ECHMET_CAES_CAES_EXTENDED_HPP
