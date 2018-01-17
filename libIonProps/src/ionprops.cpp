#include "ionprops_p.h"

namespace ECHMET {
namespace IonProps {

ECHMETReal ECHMET_CC calculateConductivity(const ComputationContext *ctx) noexcept
{
	return calculateConductivityInternal<ECHMETReal>(ctx);
}

RetCode ECHMET_CC calculateEffectiveMobilities(ComputationContext *ctx) noexcept
{
	return calculateEffectiveMobilitiesInternal<ECHMETReal>(ctx);
}

ECHMETReal ECHMET_CC calculatepH(const ComputationContext *ctx, const NonidealityCorrections corrections) noexcept
{
	return calculatepHInternal<ECHMETReal>(ctx, corrections);
}

ECHMETReal ECHMET_CC calculatepH_direct(const ECHMETReal &cH, const ECHMETReal &ionicStrength) noexcept
{
	return calculatepH_directInternal<ECHMETReal>(cH, ionicStrength);
}

RetCode ECHMET_CC correctMobilities(ComputationContext *ctx, const NonidealityCorrections corrections) noexcept
{
	return correctMobilitiesInternal<ECHMETReal>(ctx, corrections);
}

ComputationContext * ECHMET_CC makeComputationContext(const SysComp::ChemicalSystem &chemSystem, const RealVec *analyticalConcentrations, SysComp::CalculatedProperties &calcProps) noexcept
{
	return new (std::nothrow) ComputationContextImpl<ECHMETReal>(chemSystem, analyticalConcentrations, calcProps);
}

ComputationContext::~ComputationContext() noexcept {}

} // namespace IonProps

} // namespace ECHMET
