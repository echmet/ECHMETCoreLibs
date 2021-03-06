#ifndef ECHMET_IONPROPS_IONPROPS_HPP
#define ECHMET_IONPROPS_IONPROPS_HPP

/*!
 * Generator of static dispatchers based on the types
 * of \p IPReal and \p ECHMETReal.
 * If \p IPReal and \p ECHMETReal are of the same type,
 * the \p PRIMARY branch is called, otherwise the \p ALTERNATE
 * branch is called.
 */
#define DEF_DISPATCHER(NAME, RETTYPE, PRIMARY, ALTERNATE, ...) \
	template <bool, typename IPReal> \
	struct NAME##Dispatcher { \
		static RETTYPE dispatch(__VA_ARGS__); \
	}; \
	template <typename IPReal> \
	struct NAME##Dispatcher<true, IPReal> { \
		static RETTYPE dispatch(__VA_ARGS__) { PRIMARY; } \
	}; \
	template <typename IPReal> \
	struct NAME##Dispatcher<false, IPReal> { \
		static RETTYPE dispatch(__VA_ARGS__) { ALTERNATE; } \
	}

/*!
 * Calls the target function of the given dispatcher.
 */
#define DISPATCH(DISPATCHER, ...) return DISPATCHER##Dispatcher<std::is_same<IPReal, ECHMETReal>::value, IPReal>::dispatch(__VA_ARGS__)

#include "ionprops_workers.hpp"

namespace ECHMET {
namespace IonProps {

DEF_DISPATCHER(calculateEffectiveMobilities, RetCode,
	       return calculateEffectiveMobilitiesWorker(ctxImpl->chemSystem, analyticalConcentrations, calcProps),
	       return calculateEffectiveMobilitiesWorker(ctxImpl->ionicConcentrations, ctxImpl->chemSystem, analyticalConcentrations, calcProps),
	       const ComputationContextImpl<IPReal> *ctxImpl, const RealVec *analyticalConcentrations, SysComp::CalculatedProperties &calcProps);

/*!
 * Internal \p calculateEffectiveMobilities dispatching function.
 */
template <typename IPReal>
static
RetCode calculateEffectiveMobilitiesInternal(const ComputationContext *ctx, const RealVec *analyticalConcentrations, SysComp::CalculatedProperties &calcProps) noexcept
{
	const ComputationContextImpl<IPReal> *ctxImpl = dynamic_cast<const ComputationContextImpl<IPReal> *>(ctx);
	if (ctxImpl == nullptr)
		return RetCode::E_INVALID_ARGUMENT;

	DISPATCH(calculateEffectiveMobilities, ctxImpl, analyticalConcentrations, calcProps);
}


DEF_DISPATCHER(calculatepH, IPReal,
	       (void)ctxImpl; return calculatepHWorker(corrections, calcProps),
	       return calculatepHWorker<IPReal>(ctxImpl->ionicConcentrations, calcProps, corrections),
	       const ComputationContextImpl<IPReal> *ctxImpl, const NonidealityCorrections corrections, SysComp::CalculatedProperties &calcProps);

/*!
 * Internal \p calculatepH dispatching function.
 */
template <typename IPReal>
static
IPReal calculatepHInternal(const ComputationContext *ctx, const NonidealityCorrections corrections, SysComp::CalculatedProperties &calcProps) noexcept
{
	const ComputationContextImpl<IPReal> *ctxImpl = dynamic_cast<const ComputationContextImpl<IPReal> *>(ctx);
	if (ctxImpl == nullptr)
		return IPReal(0);

	DISPATCH(calculatepH, ctxImpl, corrections, calcProps);
}

DEF_DISPATCHER(calculatepH_direct, IPReal,
	       return calculatepH_directWorker(cH, ionicStrength),
	       return calculatepH_directWorker<IPReal>(cH, ionicStrength),
	       const IPReal &cH, const IPReal &ionicStrength);

/*!
 * Internal \p calculatepH_direct dispatching function.
 */
template <typename IPReal>
static
IPReal calculatepH_directInternal(const IPReal &cH, const IPReal &ionicStrength) noexcept
{
	DISPATCH(calculatepH_direct, cH, ionicStrength);
}

DEF_DISPATCHER(correctMobilities, RetCode,
	       return correctMobilitiesWorker(ctxImpl->onsagerFuossPack(), ctxImpl->ions, ctxImpl->chemSystem, calcProps, analyticalConcentrations, corrections),
	       return correctMobilitiesWorker(ctxImpl->ionicConcentrations, ctxImpl->onsagerFuossPack(), ctxImpl->ions, ctxImpl->chemSystem, calcProps, analyticalConcentrations, corrections),
	       ComputationContextImpl<IPReal> *ctxImpl, const NonidealityCorrections corrections, const RealVec *analyticalConcentrations, SysComp::CalculatedProperties &calcProps);

/*!
 * Internal \p correctMobilities dispatching function.
 */
template <typename IPReal>
static
RetCode correctMobilitiesInternal(ComputationContext *ctx, const NonidealityCorrections corrections, const RealVec *analyticalConcentrations, SysComp::CalculatedProperties &calcProps) noexcept
{
	ComputationContextImpl<IPReal> *ctxImpl = dynamic_cast<ComputationContextImpl<IPReal> *>(ctx);
	if (ctxImpl == nullptr)
		return RetCode::E_INVALID_ARGUMENT;

	DISPATCH(correctMobilities, ctxImpl, corrections, analyticalConcentrations, calcProps);
}

DEF_DISPATCHER(calculateConductivity, IPReal,
	       return calculateConductivityWorker(ctxImpl->chemSystem, calcProps),
	       return calculateConductivityWorker<IPReal>(ctxImpl->ionicConcentrations, ctxImpl->chemSystem, calcProps),
	       const ComputationContextImpl<IPReal> *ctxImpl, SysComp::CalculatedProperties &calcProps);

/*!
 * Internal \p calculateConductivity dispatching function.
 */
template <typename IPReal>
static
IPReal calculateConductivityInternal(const ComputationContext *ctx, SysComp::CalculatedProperties &calcProps) noexcept
{
	const ComputationContextImpl<IPReal> *ctxImpl = dynamic_cast<const ComputationContextImpl<IPReal> *>(ctx);
	if (ctxImpl == nullptr)
		return IPReal(0);

	DISPATCH(calculateConductivity, ctxImpl, calcProps);
}

} // namespace IonProps
} // namespace ECHMET

#endif // ECHMET_IONPROPS_IONPROPS_HPP
